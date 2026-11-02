"""Full paired comparison with identical post-hoc CS purity rules.

Native PIPs are never recomputed from filtered sets. Same frozen cohorts, variants,
phenotypes, L, tolerance and CS coverage; estimator-specific priors and native
interfaces remain explicit. Archived MultiSuSiE runtimes are not speed benchmarks.
"""

from __future__ import annotations

import argparse
import multiprocessing
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from . import multisusie_pipeline as p
from .fedfm.fed_fine_mapping import build_columns, pool_pheno, pooled_sumstats
from .fedfm.fine_mapping import run_susiex
from .inference import parse_outputs, pip_by_snp, read_table
from .multisusie_compare import score

POLICIES = {"unfiltered": 0.0, "common_purity_0.5": 0.5}
METHODS = ["multisusie", "susiex"]
_CONTEXT = None


def common_purity(members, index, grams, eligibility):
    ids = [index[s] for s in sorted(members)]
    values = []
    for gram, mask in zip(grams, eligibility, strict=True):
        if mask[ids].all():
            sub = gram[np.ix_(ids, ids)]
            den = np.sqrt(np.outer(np.diag(sub), np.diag(sub)))
            values.append(float(np.min(np.abs(sub / den))))
    return max(values, default=0.0)


def scored_rows(base, method, pips, sets, truth, purities):
    if not all(np.isfinite(v) and 0 <= v <= 1 for v in pips.values()):
        raise ValueError("Invalid PIPs")
    rows = []
    for policy, cutoff in POLICIES.items():
        kept = [s for s, purity in zip(sets, purities, strict=True) if purity >= cutoff]
        row = dict(
            base, method=method, policy=policy, converged=True, **score(pips, kept, truth, True)
        )
        row.update(
            one=1,
            discoveries95=row["n_true_pip95"] + row["n_false_pip95"],
            high_confidence_yield=int(any(pips.get(s, 0) > 0.95 for s in truth)),
        )
        rows.append(row)
    return rows


def verify_record(path):
    value = p.read(path)
    for file, digest in value["sha256"].items():
        if p.sha(file) != digest:
            raise ValueError(f"Checkpoint artifact changed: {file}")
    return value


def fit_instance(row):
    c = _CONTEXT
    instance = f"{row['locus_id']}_{row['architecture_id']}_rep{row['replicate']}"
    dest = c["out"] / instance
    checkpoint = dest / "comparison.json"
    if checkpoint.exists():
        return verify_record(checkpoint)["rows"]
    dest.mkdir(parents=True, exist_ok=True)
    truth = set(str(row["causal_snp_ids"]).split(","))
    source = c["source"] / instance / "multisusie"
    record = p.read(source / "fit.json")
    if (
        not record["converged"]
        or record["options"] != dict(c["options"], max_iter=record["options"]["max_iter"])
        or record["options"]["max_iter"] not in [1000, 10000]
    ):
        raise ValueError("Unmatched or unconverged source fit")
    if (
        record["populations"] != c["inputs"].populations
        or record["sample_sizes"] != c["inputs"].sample_sizes
    ):
        raise ValueError("Population/sample mismatch")
    variants = pd.read_csv(source / "variants.tsv", sep="\t")
    if variants.snp_id.tolist() != c["inputs"].variants.snp_id.tolist():
        raise ValueError("Source variant order differs")
    for pop, mask in zip(c["inputs"].populations, c["inputs"].eligible, strict=True):
        if not np.array_equal(variants[f"eligible_{pop}"].to_numpy(), mask):
            raise ValueError("Ancestry eligibility differs")
    pheno = [v for v in c["pheno"] if v.instance == instance]
    ss, ld, ns = [], [], []
    for col in c["columns"]:
        m = pool_pheno([v for v in pheno if v.pop == col.pop], col.variants.snp_id.to_numpy())
        frame = pooled_sumstats(col.variants, col.g_diag, col.u, m.c, m.q, m.w, col.n)
        path = dest / f"{col.pop}.sumstats"
        # Preserve float64 text precision; native SuSiEx LD remains float32.
        frame.to_csv(path, sep="\t", index=False, float_format="%.17g")
        ss.append(path)
        ld.append(col.ld_prefix)
        ns.append(col.n)
    target = dest / "susiex"
    parsed = None
    if (target / "complete.json").exists():
        saved = verify_record(target / "complete.json")
        parsed = parse_outputs(target, "cs", sorted(truth))
        elapsed = saved["runtime_s"]
    else:
        if target.exists():
            target.rename(target.with_name(f"susiex.interrupted-{time.time_ns()}"))
        start = time.perf_counter()
        for budget in [1000, 10000]:
            run_susiex(
                ss,
                ld,
                ns,
                int(c["locus"].chrom),
                int(c["locus"].start_bp),
                int(c["locus"].end_bp),
                target,
                "cs",
                str(c["binary"]),
                "plink",
                0.95,
                1.0,
                options=dict(min_purity=0, n_signals=10, max_iter=budget, tol=1e-6),
            )
            parsed = parse_outputs(target, "cs", sorted(truth))
            if parsed["converged"]:
                break
            if parsed["fit_status"] != "nonconverged":
                raise RuntimeError(f"Invalid SuSiEx output: {target}: {parsed['fit_status']}")
            if budget == 1000:
                target.rename(target.with_name(f"susiex.iteration-cap-{time.time_ns()}"))
        elapsed = time.perf_counter() - start
        if not parsed["converged"]:
            raise RuntimeError(f"Nonconverged SuSiEx at 10000 iterations: {target}")
        p.write(
            target / "complete.json",
            dict(
                runtime_s=elapsed,
                sha256={str(f): p.sha(f) for f in target.glob("cs.*") if f.is_file()},
            ),
        )
    if not parsed["converged"]:
        raise RuntimeError("Invalid completed SuSiEx checkpoint")
    snps = read_table(target / "cs.snp")
    if set(snps.SNP.astype(str)) != set(variants.snp_id):
        raise ValueError("SuSiEx variant eligibility differs from MultiSuSiE")
    base = {k: row[k] for k in p.KEY}
    base.update(
        phase=c["task"]["phase"],
        arm=c["task"]["name"],
        instance=instance,
        causal_mode=row.get("causal_mode", "shared"),
        stratum=str(c["locus"].stratum),
    )
    data = []
    for method in METHODS:
        if method == "multisusie":
            frame = pd.read_csv(source / "credible_sets.tsv", sep="\t")
            sets = [set(g.snp_id) for _, g in frame.groupby("component")]
            pips = dict(zip(variants.snp_id, variants.pip, strict=True))
        else:
            frame = read_table(target / "cs.cs")
            sets = (
                [set(g.SNP.astype(str)) for _, g in frame.groupby("CS_ID")]
                if "CS_ID" in frame
                else []
            )
            pips = pip_by_snp(target / "cs.snp")
        purities = [
            common_purity(s, c["index"], c["inputs"].xtx, c["inputs"].eligible) for s in sets
        ]
        result = scored_rows(base, method, pips, sets, truth, purities)
        for value in result:
            value.update(
                runtime_s=record["runtime_s"] if method == "multisusie" else elapsed,
                runtime_comparable=False,
                pip_definition="native_active_components"
                if method == "multisusie"
                else "native_components_with_significance_and_purity_filters_disabled",
            )
        data.extend(result)
    files = [
        source / "fit.json",
        source / "variants.tsv",
        source / "credible_sets.tsv",
        target / "complete.json",
        *target.glob("cs.*"),
    ]
    p.write(checkpoint, dict(rows=data, sha256={str(f): p.sha(f) for f in files if f.is_file()}))
    return data


def run_task(root, task, workers=8, smoke=False):
    global _CONTEXT
    source_root = Path(p.read(root / "design.json")["source"])
    source = source_root / task["phase"] / task["name"] / task["locus_id"]
    out = root / task["phase"] / task["name"] / task["locus_id"]
    if (out / "accepted.json").exists():
        verify_record(out / "accepted.json")
        return
    accepted = p.read(source / "accepted.json")
    if (
        accepted["task"] != task
        or not accepted["pass_"]
        or accepted["results_sha256"] != p.sha(source / "results.tsv")
    ):
        raise ValueError("Source acceptance mismatch")
    for name, digest in p.read(source / "phenotypes.sha256.json").items():
        if p.sha(name) != digest:
            raise ValueError("Source phenotype changed")
    cfg, loci, jobs = p.jobs_for(task["config"])
    jobs = jobs[jobs.locus_id == task["locus_id"]]
    instances = [f"{r.locus_id}_{r.architecture_id}_rep{r.replicate}" for r in jobs.itertuples()]
    locus = loci.set_index("locus_id", drop=False).loc[task["locus_id"]]
    pops = sorted(
        c.name
        for c in p.plan_columns(
            {s: dict(c.composition) for s, c in cfg.sites.items()},
            cfg.superpopulations,
            cfg.fine_mapping.min_gwas_n,
        )
    )
    geno, pheno = [], []
    ref = p.reference_variants(cfg)
    for site in cfg.sites:
        g, h = p.site_aggregates(
            cfg, site, locus, pops, instances, index=p.SiteIndex.load(cfg, site, ref)
        )
        geno.extend(g)
        pheno.extend(h)
    inputs = p.prepare_inputs(
        geno, [v for v in pheno if v.instance == instances[0]], p.MultiSuSiEOptions()
    )
    columns = build_columns(geno, pops, out / "ld", 0.005)
    from dataclasses import asdict

    _CONTEXT = dict(
        source=source,
        out=out,
        task=task,
        locus=locus,
        pheno=pheno,
        inputs=inputs,
        columns=columns,
        binary=root / "code/bin/SuSiEx",
        index={v: i for i, v in enumerate(inputs.variants.snp_id)},
        options=asdict(p.MultiSuSiEOptions()),
    )
    selected = jobs
    if smoke:
        selected = jobs[
            (jobs.replicate == 0)
            & jobs.architecture_id.isin(
                ["ncsl1_h2-0.0005_rg1", "ncsl2_h2-0.001_rg1", "ncsl3_h2-0.005_rg1"]
            )
        ]
        if len(selected) != 3:
            raise ValueError("Incomplete prespecified smoke panel")
    rows = []
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=multiprocessing.get_context("fork")
    ) as pool:
        futures = [pool.submit(fit_instance, row) for row in selected.to_dict("records")]
        for f in as_completed(futures):
            result = f.result()
            rows.extend(result)
            print("Completed", result[0]["instance"], flush=True)
    _CONTEXT = None
    frame = pd.DataFrame(rows)
    validate_grid(frame, selected)
    name = "smoke" if smoke else "results"
    frame.to_csv(out / f"{name}.tsv", sep="\t", index=False)
    p.write(
        out / f"{'smoke' if smoke else 'accepted'}.json",
        dict(
            pass_=True,
            n_instances=len(selected),
            n_rows=len(frame),
            task=task,
            sha256={str(out / f"{name}.tsv"): p.sha(out / f"{name}.tsv")},
        ),
    )


def validate_grid(frame, jobs):
    keys = [*p.KEY, "method", "policy"]
    expected = {
        (r.locus_id, r.architecture_id, r.replicate, m, f)
        for r in jobs.itertuples()
        for m in METHODS
        for f in POLICIES
    }
    if (
        frame.duplicated(keys).any()
        or set(frame[keys].itertuples(index=False, name=None)) != expected
        or not frame.converged.all()
    ):
        raise ValueError("Incomplete or failed paired grid")


def prepare(root, source, checkout, binary):
    p.verify(source)
    root.mkdir(parents=True, exist_ok=False)
    shutil.copytree(
        checkout / "src", root / "code/src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc")
    )
    (root / "code/bin").mkdir()
    shutil.copy2(binary, root / "code/bin/SuSiEx")
    (root / "logs").mkdir()
    tasks = []
    for phase in ["benchmark", "evaluation"]:
        if not p.read(source / f"{phase}.accepted.json")["pass_"]:
            raise ValueError("Unaccepted source phase")
        tasks.extend(p.read(source / f"{phase}_tasks.json"))
    p.write(root / "tasks.json", tasks)
    p.write(
        root / "design.json",
        dict(
            source=str(source),
            n_paired_instances=sum(t["n_instances"] for t in tasks),
            policies=POLICIES,
            primary=(
                "Paired unfiltered and identical post-hoc purity comparisons; calibration "
                "and recall jointly reported"
            ),
            purity=(
                "Maximum across eligible ancestries of minimum absolute within-CS "
                "correlation; same float64 centered Gram matrices for both methods"
            ),
            pip_policy="Native PIPs unchanged by post-hoc CS filtering",
            susiex=dict(
                min_purity=0,
                pval_thresh=1,
                n_signals=10,
                tol=1e-6,
                coverage=0.95,
                max_iter=[1000, 10000],
            ),
            precision=(
                "SuSiEx summary statistics 17 significant digits; native float32 LD. "
                "MultiSuSiE native float64 raw moments."
            ),
            runtime="Archived MultiSuSiE timing; do not claim speed superiority",
            inference_scope=(
                "Fixed configured estimators and frozen synthetic design; no superiority assumption"
            ),
            stratification=["architecture_id", "causal_mode"],
            bootstrap=dict(
                unit="locus",
                draws=2000,
                seed=20261029,
                contrast="MultiSuSiE minus SuSiEx, jointly resampled",
            ),
            selection=(
                "Previously frozen MultiSuSiE-selected cohort held fixed for both methods; "
                "no new selection"
            ),
        ),
    )
    files = [
        *list((root / "code").rglob("*.py")),
        root / "code/bin/SuSiEx",
        root / "tasks.json",
        root / "design.json",
        source / "snapshot.json",
        source / "selection.lock.json",
    ]
    files += [Path(t["config"]) for t in tasks]
    files += [source / t["phase"] / t["name"] / t["locus_id"] / "accepted.json" for t in tasks]
    p.write(root / "snapshot.json", {"sha256": {str(f): p.sha(f) for f in files}})


def verify(root):
    for name, digest in p.read(root / "snapshot.json")["sha256"].items():
        if p.sha(name) != digest:
            raise ValueError(f"Frozen comparison input changed: {name}")
    p.verify(Path(p.read(root / "design.json")["source"]))


def main():
    a = argparse.ArgumentParser()
    a.add_argument("stage", choices=["prepare", "run", "smoke", "collect"])
    a.add_argument("--root", type=Path, required=True)
    a.add_argument("--source", type=Path)
    a.add_argument("--checkout", type=Path, default=Path.cwd())
    a.add_argument("--susiex", type=Path)
    a.add_argument("--index", type=int, default=0)
    a.add_argument("--workers", type=int, default=8)
    a.add_argument("--indices", type=Path)
    args = a.parse_args()
    root = args.root.resolve()
    if args.stage == "prepare":
        prepare(root, args.source.resolve(), args.checkout.resolve(), args.susiex.resolve())
        return
    verify(root)
    if args.stage == "collect":
        collect(root)
        return
    index = p.read(args.indices)[args.index] if args.indices else args.index
    run_task(root, p.read(root / "tasks.json")[index], args.workers, args.stage == "smoke")


def collect(root):
    frames = []
    for task in p.read(root / "tasks.json"):
        out = root / task["phase"] / task["name"] / task["locus_id"]
        accepted = verify_record(out / "accepted.json")
        if accepted["task"] != task or not accepted["pass_"]:
            raise ValueError("Task acceptance mismatch")
        f = pd.read_csv(out / "results.tsv", sep="\t")
        _, _, jobs = p.jobs_for(task["config"])
        validate_grid(f, jobs[jobs.locus_id == task["locus_id"]])
        frames.append(f)
    frame = pd.concat(frames, ignore_index=True)
    frame.to_csv(root / "comparison.tsv", sep="\t", index=False)
    summarize(frame, root)
    summarize(frame, root, strata=("architecture_id", "causal_mode"))
    p.write(
        root / "accepted.json",
        dict(
            pass_=True,
            n_paired_instances=len(frame) // 4,
            n_rows=len(frame),
            sha256={
                str(f): p.sha(f)
                for f in [
                    root / "comparison.tsv",
                    root / "summary.tsv",
                    root / "paired_differences.tsv",
                    root / "summary_stratified.tsv",
                    root / "paired_differences_stratified.tsv",
                ]
            },
        ),
    )


METRICS = {
    "power": ("any_causal_captured", "one"),
    "causal_recall": ("n_causal_captured", "n_causal"),
    "cs_coverage": ("n_cs_containing_causal", "n_credible_sets"),
    "mean_cs_size": ("total_cs_snps", "n_credible_sets"),
    "fdr_pip95": ("n_false_pip95", "discoveries95"),
    "high_confidence_yield": ("high_confidence_yield", "one"),
}


def summarize(frame, root, strata=()):
    summaries = []
    contrasts = []
    groups = ["phase", "arm", "policy", *strata]
    for values, g in frame.groupby(groups):
        identity = dict(zip(groups, values, strict=True))
        loci = sorted(g.locus_id.unique())
        rng = np.random.default_rng(20261029)
        draws = rng.integers(len(loci), size=(2000, len(loci)))
        for metric, (num, den) in METRICS.items():
            estimates = {}
            boots = {}
            for method in METHODS:
                v = (
                    g[g.method == method]
                    .groupby("locus_id")[[num, den]]
                    .sum()
                    .reindex(loci)
                    .to_numpy(dtype=float)
                )
                total = v.sum(axis=0)
                estimate = total[0] / total[1] if total[1] > 0 else np.nan
                b = v[draws].sum(axis=1)
                with np.errstate(divide="ignore", invalid="ignore"):
                    boot = b[:, 0] / b[:, 1]
                finite = np.isfinite(boot)
                ci = np.quantile(boot[finite], [0.025, 0.975]) if finite.all() else [np.nan, np.nan]
                summaries.append(
                    dict(
                        **identity,
                        method=method,
                        metric=metric,
                        estimate=estimate,
                        ci_low=ci[0],
                        ci_high=ci[1],
                    )
                )
                estimates[method] = estimate
                boots[method] = boot
            delta = boots["multisusie"] - boots["susiex"]
            ci = (
                np.quantile(delta, [0.025, 0.975]) if np.isfinite(delta).all() else [np.nan, np.nan]
            )
            contrasts.append(
                dict(
                    **identity,
                    metric=metric,
                    difference=estimates["multisusie"] - estimates["susiex"],
                    ci_low=ci[0],
                    ci_high=ci[1],
                )
            )
    suffix = "_stratified" if strata else ""
    pd.DataFrame(summaries).to_csv(root / f"summary{suffix}.tsv", sep="\t", index=False)
    pd.DataFrame(contrasts).to_csv(root / f"paired_differences{suffix}.tsv", sep="\t", index=False)
    (root / "RESULTS.md").write_text(
        "# Matched-filter MultiSuSiE–SuSiEx comparison\n"
        "\n"
        "Complete paired grids are in comparison.tsv. summary.tsv reports each "
        "method; paired_differences.tsv reports MultiSuSiE minus SuSiEx with paired "
        "locus-bootstrap intervals.\n"
        "\n"
        "Compare recall jointly with coverage, set size and FDR. No winner is "
        "inferred from power alone. Filters and resampling were frozen before the "
        "run. Native priors, PIP definitions and input precision remain method-"
        "specific; see design.json. Archived timings do not support speed claims.\n"
    )


if __name__ == "__main__":
    main()
