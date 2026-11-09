"""Checkpointed MultiSuSiE benchmark and development/evaluation pipeline.

Reuse frozen synthetic cohorts and phenotypes, never the previous estimator's
selection or results. Each work item is an arm/locus; only accepted work is reused.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .fedfm.fed_fine_mapping import SiteIndex, reference_variants, site_aggregates
from .fedfm.fine_mapping import plan_columns
from .fedfm.utils import load_config
from .multisusie import MultiSuSiEOptions, prepare_inputs, run_multisusie
from .multisusie_compare import centralized_blocks, compare_fits, score

KEY = ["locus_id", "architecture_id", "replicate"]


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def jobs_for(config):
    cfg = load_config(Path(config))
    loci = pd.read_csv(cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t")
    manifest = pd.read_csv(cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv", sep="\t")
    jobs = manifest[
        manifest.locus_id.isin(loci.locus_id) & (manifest.replicate < cfg.architecture.replicates)
    ]
    if jobs.empty or jobs.duplicated(KEY).any() or loci.locus_id.duplicated().any():
        raise ValueError("Empty or duplicate input grid")
    return cfg, loci, jobs.sort_values(KEY)


def freeze_config(root, phase, name, source, raw=None):
    cfg = load_config(source)
    raw = yaml.safe_load(source.read_text()) if raw is None else raw
    raw["paths"] = {k: str(cfg.resolved_path(k)) for k in raw["paths"]}
    raw["paths"]["reports_dir"] = str(root / phase / name)
    raw["paths"]["logs_dir"] = str(root / "logs")
    path = root / "configs" / phase / f"{name}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(raw, sort_keys=False))
    return path


def make_tasks(root, phase, arms):
    tasks = []
    for arm in arms:
        _, loci, jobs = jobs_for(arm["config"])
        for lid in sorted(jobs.locus_id.unique()):
            tasks.append(
                {
                    **arm,
                    "phase": phase,
                    "locus_id": lid,
                    "n_instances": int((jobs.locus_id == lid).sum()),
                }
            )
    write(root / f"{phase}_tasks.json", tasks)
    return tasks


def prepare(root, source, checkout):
    root.mkdir(parents=True, exist_ok=False)
    (root / "logs").mkdir()
    shutil.copytree(
        checkout / "src", root / "code/src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc")
    )
    shutil.copy2(checkout / "scripts/fine-mapping/multisusie_pipeline.pbs", root / "job.pbs")
    benchmark = source / "fine-mapping-full-rerun-20260914"
    composition = source / "fine-mapping-final-20260928"
    design = read(composition / "design.json")
    design.update(
        method="MultiSuSiE",
        options=asdict(MultiSuSiEOptions()),
        source_benchmark=str(benchmark),
        source_composition=str(composition),
        execution="local synthetic site stages; aggregate-only coordinator inference",
    )
    write(root / "design.json", design)
    arms = []
    for name in [
        "federation",
        "anl",
        "covenant",
        "mbzuai",
        "federation_50k",
        "federation_50k_seed2",
        "federation_50k_seed3",
    ]:
        source_config = benchmark / (
            "pipeline_config.yaml" if name == "federation" else f"arms/{name}.yaml"
        )
        cfg = freeze_config(root, "benchmark", name, source_config)
        arms.append(dict(name=name, config=str(cfg), central_check=name == "federation"))
    make_tasks(root, "benchmark", arms)
    unique = {t["name"]: t for t in read(composition / "development_tasks.json")}
    arms = []
    for name, task in unique.items():
        cfg = freeze_config(root, "development", name, Path(task["config"]))
        arms.append(
            dict(
                name=name,
                config=str(cfg),
                candidate=task["candidate"],
                sampling_seed=task["sampling_seed"],
                central_check=False,
            )
        )
    make_tasks(root, "development", arms)
    # Source and metadata fingerprints; BED files use stat identity (large files).
    hashes = {str(p): sha(p) for p in (root / "code").rglob("*.py")}
    inputs = set()
    for config in (root / "configs").rglob("*.yaml"):
        cfg, _, _ = jobs_for(config)
        inputs.update(
            [
                config,
                cfg.resolved_path("loci_dir") / "selected_loci.tsv",
                cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv",
            ]
        )
        for site in cfg.sites:
            prefix = cfg.site_dir(site) / f"{site}_chr{cfg.chromosome}"
            inputs.update(prefix.with_suffix(s) for s in [".bed", ".bim", ".fam"])
            inputs.add(cfg.site_dir(site) / f"{site}_manifest.tsv")
    stats = {}
    for p in inputs:
        if p.suffix == ".bed":
            s = p.stat()
            stats[str(p)] = [s.st_size, s.st_mtime_ns]
        else:
            hashes[str(p)] = sha(p)
    from MultiSuSiE import susiepy_ss

    hashes[susiepy_ss.__file__] = sha(susiepy_ss.__file__)
    for p in [root / "design.json", root / "benchmark_tasks.json", root / "development_tasks.json"]:
        hashes[str(p)] = sha(p)
    write(root / "snapshot.json", dict(sha256=hashes, bed_stat_identity=stats))


def verify(root):
    snap = read(root / "snapshot.json")
    for p, expected in snap["sha256"].items():
        if sha(p) != expected:
            raise ValueError(f"Frozen input changed: {p}")
    for p, expected in snap["bed_stat_identity"].items():
        s = Path(p).stat()
        if [s.st_size, s.st_mtime_ns] != expected:
            raise ValueError(f"BED identity changed: {p}")


def metrics(path, truth, geno, pheno, options):
    record = read(path / "fit.json")
    variants = pd.read_csv(path / "variants.tsv", sep="\t")
    cs = pd.read_csv(path / "credible_sets.tsv", sep="\t")
    pips = dict(zip(variants.snp_id, variants.pip, strict=True))
    sets = [set(g.snp_id) for _, g in cs.groupby("component")]
    result = score(pips, sets, truth, record["converged"])
    if not record["converged"]:
        raise RuntimeError(f"Nonconverged fit: {path}")
    captured_sizes = [len(s) for s in sets if s & truth]
    # Purity for display only: maximum, over ancestry, of minimum absolute LD.
    inputs = prepare_inputs(geno, pheno, options)
    idx = {v: i for i, v in enumerate(inputs.variants.snp_id)}
    purities = []
    for members in sets:
        ids = [idx[v] for v in sorted(members)]
        values = []
        for gram, eligible in zip(inputs.xtx, inputs.eligible, strict=True):
            if not eligible[ids].all():
                continue
            sub = gram[np.ix_(ids, ids)]
            den = np.sqrt(np.outer(np.diag(sub), np.diag(sub)))
            values.append(float(np.min(np.abs(sub / den))))
        purities.append(max(values, default=0.0))
    result.update(
        converged=True,
        fit_status=record["fit_status"],
        runtime_s=record["runtime_s"],
        pip_definition=record["pip_definition"],
        cs_sizes_json=json.dumps(list(map(len, sets))),
        best_cs_size=min(captured_sizes) if captured_sizes else None,
        causal_pip_max=max((pips.get(v, 0.0) for v in truth), default=0.0),
        causal_pip_mean=float(np.mean([pips.get(v, 0.0) for v in truth])),
        top_pip=max(pips.values()),
        mean_cs_purity=float(np.mean(purities)) if purities else None,
    )
    return result


def work_item(root, task):
    root = Path(root)
    out = root / task["phase"] / task["name"] / task["locus_id"]
    out.mkdir(parents=True, exist_ok=True)
    cfg, loci, jobs = jobs_for(task["config"])
    jobs = jobs[jobs.locus_id == task["locus_id"]]
    expected = set(jobs[KEY].itertuples(index=False, name=None))
    if (out / "accepted.json").exists():
        acceptance = read(out / "accepted.json")
        if (
            acceptance["task"] != task
            or acceptance["results_sha256"] != sha(out / "results.tsv")
            or acceptance["parity_sha256"] != sha(out / "parity.json")
        ):
            raise ValueError(f"Checkpoint changed: {out}")
        for p, digest in read(out / "phenotypes.sha256.json").items():
            if sha(p) != digest:
                raise ValueError(f"Phenotype changed: {p}")
        frame = pd.read_csv(out / "results.tsv", sep="\t")
        if set(frame[KEY].itertuples(index=False, name=None)) != expected or len(frame) != len(
            expected
        ):
            raise ValueError(f"Invalid checkpoint: {out}")
        return str(out)
    locus = loci.set_index("locus_id", drop=False).loc[task["locus_id"]]
    instances = [f"{r.locus_id}_{r.architecture_id}_rep{r.replicate}" for r in jobs.itertuples()]
    # Hash the actual phenotypes used by each checkpoint, not merely the manifest.
    provenance = {
        str(p): sha(p)
        for site in cfg.sites
        for inst in instances
        for p in [cfg.resolved_path("ground_truth_dir") / "phenotypes" / site / f"{inst}.pheno"]
    }
    if (out / "phenotypes.sha256.json").exists() and read(
        out / "phenotypes.sha256.json"
    ) != provenance:
        raise ValueError(f"Phenotypes changed during restart: {out}")
    write(out / "phenotypes.sha256.json", provenance)
    columns = plan_columns(
        {s: dict(c.composition) for s, c in cfg.sites.items()},
        cfg.superpopulations,
        cfg.fine_mapping.min_gwas_n,
    )
    pops = sorted(c.name for c in columns)
    expected_n = {c.name: c.n for c in columns}
    ref = reference_variants(cfg)
    geno, pheno = [], []
    for site in cfg.sites:
        g, p = site_aggregates(
            cfg, site, locus, pops, instances, index=SiteIndex.load(cfg, site, ref)
        )
        geno.extend(g)
        pheno.extend(p)
    if {pop: sum(g.n for g in geno if g.pop == pop) for pop in pops} != expected_n:
        raise ValueError("Observed sample sizes differ from design")
    if task["central_check"]:
        cg, cp = centralized_blocks(cfg, locus, pops, instances)
    results, parities = [], []
    options = MultiSuSiEOptions()
    for row, inst in zip(jobs.itertuples(), instances, strict=True):
        ps = [p for p in pheno if p.instance == inst]
        target = out / inst
        for directory, g, p in [(target / "multisusie", geno, ps)] + (
            [(target / "centralized", cg, [p for p in cp if p.instance == inst])]
            if task["central_check"]
            else []
        ):
            if directory.exists():
                # Only a complete converged artifact is a reusable fit checkpoint.
                ok = (directory / "fit.json").exists() and read(directory / "fit.json").get(
                    "converged"
                )
                ok = ok and all(
                    (directory / f).exists()
                    for f in ["variants.tsv", "credible_sets.tsv", "posterior.npz"]
                )
                if not ok:
                    archive = directory.with_name(directory.name + f".failed-{time.time_ns()}")
                    directory.rename(archive)
                else:
                    continue
            run_multisusie(g, p, directory, options)
        if task["central_check"]:
            parity = compare_fits(target / "multisusie", target / "centralized")
            parities.append(dict(instance=inst, **parity))
            if not parity["pass"]:
                write(out / "parity.json", parities)
                raise RuntimeError(f"Centralized parity failed: {inst}")
        result = dict(
            instance=inst,
            locus_id=row.locus_id,
            architecture_id=row.architecture_id,
            replicate=row.replicate,
            stratum=locus.stratum,
            causal_mode=getattr(row, "causal_mode", "shared"),
            arm=task["name"],
            analyzed_n=sum(expected_n.values()),
            **metrics(
                target / "multisusie", set(str(row.causal_snp_ids).split(",")), geno, ps, options
            ),
        )
        for k in ["candidate", "sampling_seed"]:
            if k in task:
                result[k] = task[k]
        results.append(result)
        pd.DataFrame(results).to_csv(out / "results.tsv", sep="\t", index=False)
    write(out / "parity.json", parities)
    write(
        out / "accepted.json",
        dict(
            pass_=True,
            n_instances=len(results),
            task=task,
            results_sha256=sha(out / "results.tsv"),
            parity_sha256=sha(out / "parity.json"),
        ),
    )
    return str(out)


def collect(root, phase):
    frames = []
    for task in read(root / f"{phase}_tasks.json"):
        out = root / phase / task["name"] / task["locus_id"]
        acceptance = read(out / "accepted.json")
        if (
            not acceptance["pass_"]
            or acceptance["task"] != task
            or acceptance["results_sha256"] != sha(out / "results.tsv")
            or acceptance["parity_sha256"] != sha(out / "parity.json")
        ):
            raise ValueError(f"Unaccepted task {out}")
        frame = pd.read_csv(out / "results.tsv", sep="\t")
        _, _, jobs = jobs_for(task["config"])
        expected = set(
            jobs.loc[jobs.locus_id == task["locus_id"], KEY].itertuples(index=False, name=None)
        )
        if (
            len(frame) != len(expected)
            or set(frame[KEY].itertuples(index=False, name=None)) != expected
            or not frame.converged.all()
        ):
            raise ValueError(f"Incomplete result grid: {out}")
        frames.append(frame)
    frame = pd.concat(frames, ignore_index=True)
    frame.to_csv(root / f"{phase}_results.tsv", sep="\t", index=False)
    write(root / f"{phase}.accepted.json", dict(pass_=True, n_fits=len(frame), n_tasks=len(frames)))
    return frame


def run(root, phase, workers, shard=0, shards=1):
    verify(root)
    tasks = read(root / f"{phase}_tasks.json")[shard::shards]
    failures = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(work_item, root, task): task for task in tasks}
        for f in as_completed(futures):
            task = futures[f]
            try:
                print("Completed", f.result(), flush=True)
            except Exception as exc:
                failures.append(dict(task=task, error=repr(exc)))
                print("FAILED", task, repr(exc), flush=True)
    write(root / f"{phase}.failures.{shard}.json", failures)
    if failures:
        raise RuntimeError(f"{len(failures)} failed tasks; successful checkpoints preserved")
    if shards == 1:
        collect(root, phase)


def select(root):
    frame = collect(root, "development")
    design = read(root / "design.json")
    ranks = []
    for candidate, g in frame.groupby("candidate"):
        total_cs = g.n_credible_sets.sum()
        discoveries = g.n_true_pip95.sum() + g.n_false_pip95.sum()
        coverage = float(g.n_cs_containing_causal.sum() / total_cs) if total_cs else 0.0
        fdr = float(g.n_false_pip95.sum() / discoveries) if discoveries else 0.0
        choice = next(c for c in design["candidates"] if c["name"] == candidate)
        ranks.append(
            dict(
                candidate=candidate,
                power=float(g.any_causal_captured.mean()),
                cs_coverage=coverage,
                fdr_pip95=fdr,
                solo_fallback=bool(choice["solo"]),
                eligible=coverage >= design["calibration_screen"]["min_cs_coverage"]
                and fdr <= design["calibration_screen"]["max_fdr_pip95"],
            )
        )
    ranks = pd.DataFrame(ranks).sort_values(
        ["eligible", "power", "solo_fallback", "candidate"], ascending=[False, False, False, True]
    )
    ranks.to_csv(root / "development_ranking.tsv", sep="\t", index=False)
    if not ranks.eligible.any():
        raise ValueError("No candidate passes the prespecified calibration screen")
    winner = next(c for c in design["candidates"] if c["name"] == ranks.iloc[0].candidate)
    selected = dict(
        **winner,
        evaluation_used_for_selection=False,
        development_results_sha256=sha(root / "development_results.tsv"),
    )
    if (root / "selection.json").exists() and read(root / "selection.json") != selected:
        raise ValueError("Selection is already frozen and differs")
    write(root / "selection.json", selected)
    source = Path(design["source_composition"])
    eval_raw = yaml.safe_load((source / "evaluation/pipeline_config.yaml").read_text())
    arms = []
    dev = {
        t["name"]: t
        for t in read(root / "development_tasks.json")
        if t["candidate"] == winner["name"]
    }
    for name, task in dev.items():
        raw = yaml.safe_load(Path(task["config"]).read_text())
        raw["architecture"] = eval_raw["architecture"]
        for key in ["loci_dir", "ground_truth_dir"]:
            raw["paths"][key] = eval_raw["paths"][key]
        raw["master_seed"] = eval_raw["master_seed"]
        cfg = root / "configs/evaluation" / f"{name}.yaml"
        cfg.parent.mkdir(parents=True, exist_ok=True)
        raw["paths"]["reports_dir"] = str(root / "evaluation" / name)
        cfg.write_text(yaml.safe_dump(raw, sort_keys=False))
        arms.append(
            dict(
                name=name,
                config=str(cfg),
                candidate=winner["name"],
                sampling_seed=task["sampling_seed"],
                central_check=False,
            )
        )
    for name in ["federation", "anl", "covenant", "mbzuai"]:
        cfg = freeze_config(root, "evaluation", name, source / "evaluation/arms" / f"{name}.yaml")
        arms.append(
            dict(
                name=name,
                config=str(cfg),
                candidate=name,
                sampling_seed=design["sampling_seeds"][0],
                central_check=False,
            )
        )
    make_tasks(root, "evaluation", arms)
    write(
        root / "selection.lock.json",
        {
            str(p): sha(p)
            for p in [
                root / "selection.json",
                root / "evaluation_tasks.json",
                *list((root / "configs/evaluation").glob("*.yaml")),
            ]
        },
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["prepare", "run", "collect", "select"])
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--source", type=Path)
    p.add_argument("--checkout", type=Path, default=Path.cwd())
    p.add_argument(
        "--phase", choices=["benchmark", "development", "evaluation"], default="benchmark"
    )
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--shards", type=int, default=1)
    a = p.parse_args()
    root = a.root.resolve()
    if a.stage == "prepare":
        prepare(root, a.source.resolve(), a.checkout.resolve())
    elif a.stage == "run":
        if a.phase == "evaluation":
            for name, digest in read(root / "selection.lock.json").items():
                if sha(name) != digest:
                    raise ValueError(f"Frozen evaluation selection changed: {name}")
        run(root, a.phase, a.workers, a.shard, a.shards)
    elif a.stage == "collect":
        collect(root, a.phase)
    else:
        select(root)


if __name__ == "__main__":
    main()
