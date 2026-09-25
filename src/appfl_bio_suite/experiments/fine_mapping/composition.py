"""Select a fixed-N cohort on development regions, then evaluate one frozen choice.

Existing phenotypes, inference settings and solo baselines are reused. Candidate
selection never reads evaluation outcomes. Overlapping genomic windows stay in the
same split. Three sampling seeds are averaged; they are not extra independent loci.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .arms import KEY, SITE_ORDER, Arm, build_configs, run_arm
from .inference import validate_result_grid
from .reporting import clustered_ratio, paired_arm_differences

SEEDS = (20260921, 20260922, 20260923)
SMART_ARM = "federation_smart_50k"


def candidates() -> list[dict]:
    return [
        {
            "name": f"afr{afr // 1000}_eur{(50000 - afr) // 1000}_{policy}",
            "ancestry_counts": {p: n for p, n in {"AFR": afr, "EUR": 50000 - afr}.items() if n},
            "site_priority": list(priority),
        }
        for afr in (50000, 45000, 40000, 35000, 25000)
        for policy, priority in (
            ("proportional", ()),
            ("covenant_first", ("covenant", "anl", "mbzuai")),
        )
    ]


def split_regions(loci: pd.DataFrame, seed: int = 20260921) -> pd.DataFrame:
    """Stratify connected, overlapping-window components; split before reading results."""
    ordered = loci.sort_values(["chrom", "start_bp", "end_bp", "locus_id"]).copy()
    region, previous_chrom, stop = -1, None, -1
    regions = []
    for row in ordered.itertuples():
        if row.chrom != previous_chrom or row.start_bp > stop:
            region += 1
            stop = row.end_bp
        else:
            stop = max(stop, row.end_bp)
        previous_chrom = row.chrom
        regions.append(region)
    ordered["region"] = regions
    groups = ordered.groupby("region").stratum.agg(lambda x: sorted(x.mode())[0])
    rng = np.random.default_rng(seed)
    dev = set()
    for stratum in sorted(groups.unique()):
        ids = groups.index[groups.eq(stratum)].to_numpy()
        if len(ids) < 2:
            raise ValueError(f"Need at least two independent regions in stratum {stratum}")
        dev.update(rng.permutation(ids)[: max(1, len(ids) // 3)].tolist())
    ordered["split"] = np.where(ordered.region.isin(dev), "development", "evaluation")
    return ordered.sort_values("locus_id").reset_index(drop=True)


def _json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def prepare(source: Path, root: Path, development_reps: int = 2) -> None:
    source, root = source.resolve(), root.resolve()
    if root.exists():
        raise FileExistsError(f"Use a new experiment directory: {root}")
    cfg = yaml.safe_load((source / "pipeline_config.yaml").read_text())
    if not 1 <= development_reps <= cfg["architecture"]["replicates"]:
        raise ValueError("Development replicates must be within the source replicate range")
    split = split_regions(pd.read_csv(source / "loci/selected_loci.tsv", sep="\t"))
    truth_dir = Path(cfg["paths"]["ground_truth_dir"])
    truth = pd.read_csv(truth_dir / "causal_manifest.tsv", sep="\t")
    truth = truth.loc[truth.causal_mode.eq("shared")]
    if truth.empty:
        raise ValueError("This experiment requires shared-causal simulation instances")
    root.mkdir(parents=True)
    (root / "logs").mkdir()
    split.to_csv(root / "split.tsv", sep="\t", index=False)
    _json(
        root / "design.json",
        {
            "source": str(source),
            "seeds": list(SEEDS),
            "candidates": candidates(),
            "selection_metric": "mean_unconditional_power_across_development_instances_and_seeds",
            "tie_break": "candidate_name_ascending",
            "development_replicates": development_reps,
            "evaluation_replicates": cfg["architecture"]["replicates"],
            "causal_mode": "shared",
            "split_seed": 20260921,
            "caution": "Existing simulated panel; evaluation loci excluded from selection.",
        },
    )
    for phase in ("development", "evaluation"):
        phase_root = root / phase
        (phase_root / "loci").mkdir(parents=True)
        (phase_root / "ground_truth").mkdir()
        phase_loci = split.loc[split.split.eq(phase)]
        phase_loci.drop(columns=["region", "split"]).to_csv(
            phase_root / "loci/selected_loci.tsv", sep="\t", index=False
        )
        phase_truth = truth.loc[truth.locus_id.isin(phase_loci.locus_id)]
        phase_truth.to_csv(phase_root / "ground_truth/causal_manifest.tsv", sep="\t", index=False)
        for path in truth_dir.iterdir():
            if path.name != "causal_manifest.tsv":
                (phase_root / "ground_truth" / path.name).symlink_to(path.resolve())
        local = yaml.safe_load((source / "pipeline_config.yaml").read_text())
        local["architecture"].setdefault("ancestry_divergent_causal", {})["enabled"] = False
        for key in ("loci_dir", "ground_truth_dir", "reports_dir", "logs_dir"):
            local["paths"][key] = str(phase_root / key.removesuffix("_dir"))
        local["architecture"]["replicates"] = (
            development_reps if phase == "development" else cfg["architecture"]["replicates"]
        )
        (phase_root / "pipeline_config.yaml").write_text(yaml.safe_dump(local, sort_keys=False))
    _build_tasks(root, "development", candidates(), shards=2)


def _build_tasks(root: Path, phase: str, choices: list[dict], shards: int) -> None:
    arms, metadata = [], {}
    for candidate in choices:
        for seed in SEEDS:
            name = f"{candidate['name']}_seed{seed}"
            arms.append(
                Arm(
                    name,
                    SITE_ORDER,
                    seed=seed,
                    ancestry_counts=tuple(candidate["ancestry_counts"].items()),
                    site_priority=tuple(candidate["site_priority"]),
                )
            )
            metadata[name] = {"candidate": candidate["name"], "sampling_seed": seed}
    configs = build_configs(
        root / phase / "pipeline_config.yaml", root / phase / "arms", tuple(arms)
    )
    tasks = [
        {"name": name, "config": str(path), **metadata[name], "shard": shard, "shards": shards}
        for name, path in configs.items()
        for shard in range(shards)
    ]
    _json(root / f"{phase}_tasks.json", tasks)


def _expected(config: Path, shard: int | None = None, shards: int = 1) -> set[tuple]:
    cfg = yaml.safe_load(config.read_text())
    loci = pd.read_csv(Path(cfg["paths"]["loci_dir"]) / "selected_loci.tsv", sep="\t")
    if shard is not None:
        loci = loci.iloc[shard::shards]
    truth = pd.read_csv(Path(cfg["paths"]["ground_truth_dir"]) / "causal_manifest.tsv", sep="\t")
    truth = truth.loc[
        truth.locus_id.isin(loci.locus_id) & truth.replicate.lt(cfg["architecture"]["replicates"])
    ]
    return set(truth[KEY].itertuples(index=False, name=None))


def run_task(root: Path, phase: str, index: int, workers: int) -> None:
    tasks = json.loads((root / f"{phase}_tasks.json").read_text())
    task = tasks[index]
    config = Path(task["config"])
    cfg = yaml.safe_load(config.read_text())
    result = (
        Path(cfg["paths"]["reports_dir"])
        / "fine_mapping"
        / (f"fm_results.part{task['shard']:03d}of{task['shards']:03d}.tsv")
    )
    expected = _expected(config, task["shard"], task["shards"])
    if result.exists():
        validate_result_grid(pd.read_csv(result, sep="\t"), expected)
        return
    run_arm(
        task["name"], config, n_workers=workers, shard_index=task["shard"], n_shards=task["shards"]
    )
    validate_result_grid(pd.read_csv(result, sep="\t"), expected)


def collect(root: Path, phase: str) -> pd.DataFrame:
    tasks = json.loads((root / f"{phase}_tasks.json").read_text())
    frames = []
    for task in tasks:
        cfg_path = Path(task["config"])
        cfg = yaml.safe_load(cfg_path.read_text())
        result = (
            Path(cfg["paths"]["reports_dir"])
            / "fine_mapping"
            / (f"fm_results.part{task['shard']:03d}of{task['shards']:03d}.tsv")
        )
        frame = pd.read_csv(result, sep="\t")
        validate_result_grid(frame, _expected(cfg_path, task["shard"], task["shards"]))
        frame["candidate"] = task["candidate"]
        frame["sampling_seed"] = task["sampling_seed"]
        frames.append(frame)
    merged = pd.concat(frames, ignore_index=True)
    for _, frame in merged.groupby(["candidate", "sampling_seed"]):
        validate_result_grid(frame, _expected(root / phase / "pipeline_config.yaml"))
    merged.to_csv(root / f"{phase}_results.tsv", sep="\t", index=False)
    return merged


def rank_candidates(frame: pd.DataFrame) -> pd.DataFrame:
    """Predeclared objective: unconditional power; never optimize on test outcomes."""
    table = (
        frame.groupby("candidate")
        .agg(
            power=("any_causal_captured", "mean"),
            high_confidence_yield=("causal_pip_max", lambda x: (x > 0.95).mean()),
            n_instances=("candidate", "size"),
            n_loci=("locus_id", "nunique"),
        )
        .reset_index()
    )
    return table.sort_values(["power", "candidate"], ascending=[False, True]).reset_index(drop=True)


def select(root: Path) -> None:
    if (root / "selection.json").exists():
        raise FileExistsError("Selection already frozen; use its existing evaluation tasks")
    development = collect(root, "development")
    ranks = rank_candidates(development)
    ranks.to_csv(root / "development_ranking.tsv", sep="\t", index=False)
    design = json.loads((root / "design.json").read_text())
    winner = next(c for c in design["candidates"] if c["name"] == ranks.iloc[0].candidate)
    _build_tasks(root, "evaluation", [winner], shards=6)
    _json(
        root / "selection.json",
        {
            **winner,
            "development_power": float(ranks.iloc[0].power),
            "development_results_sha256": hashlib.sha256(
                (root / "development_results.tsv").read_bytes()
            ).hexdigest(),
            "evaluation_used_for_selection": False,
        },
    )


def finish(root: Path) -> None:
    selection = json.loads((root / "selection.json").read_text())
    design = json.loads((root / "design.json").read_text())
    selected = collect(root, "evaluation")
    if set(selected.candidate) != {selection["name"]}:
        raise ValueError("Evaluation candidate differs from the frozen selection")
    selected["arm"] = SMART_ARM
    selected["smart_composition"] = ", ".join(
        f"{n:,} {pop}" for pop, n in selection["ancestry_counts"].items()
    ) + (
        "; Covenant-first recruitment"
        if selection["site_priority"]
        else "; proportional recruitment within ancestry"
    )
    baselines = pd.read_csv(Path(design["source"]) / "arms/fm_results_by_arm.tsv", sep="\t")
    expected = _expected(root / "evaluation/pipeline_config.yaml")
    index = pd.MultiIndex.from_tuples(sorted(expected), names=KEY)
    frames = [selected]
    for name in (*SITE_ORDER, "federation"):
        original = baselines.loc[baselines.arm.eq(name)].copy()
        frame = original.set_index(KEY).loc[index].reset_index()
        validate_result_grid(frame, expected)
        for seed in SEEDS:
            frames.append(frame.assign(sampling_seed=seed))
    headline = pd.concat(frames, ignore_index=True)
    headline["causal_mode"] = "shared"
    headline["evaluation_split"] = "evaluation"
    headline.to_csv(root / "fm_results_by_arm.tsv", sep="\t", index=False)
    paired_arm_differences(headline).to_csv(root / "paired_vs_full.tsv", sep="\t", index=False)
    paired_arm_differences(headline, reference=SMART_ARM).to_csv(
        root / "paired_vs_smart.tsv", sep="\t", index=False
    )
    rows = []
    for arm, group in headline.groupby("arm"):
        group = group.copy()
        group["opportunities"] = 1
        group["hits"] = group.any_causal_captured.astype(int)
        group["confident"] = group.causal_pip_max.gt(0.95).astype(int)
        group["discoveries95"] = group.n_true_pip95 + group.n_false_pip95
        for metric, numerator, denominator in (
            ("power", "hits", "opportunities"),
            ("high_confidence_yield", "confident", "opportunities"),
            ("cs_coverage", "n_cs_containing_causal", "n_credible_sets"),
            ("fdr_pip95", "n_false_pip95", "discoveries95"),
        ):
            point, lo, hi = clustered_ratio(group, numerator, denominator)
            rows.append(dict(arm=arm, metric=metric, estimate=point, ci_low=lo, ci_high=hi))
    pd.DataFrame(rows).to_csv(root / "evaluation_summary.tsv", sep="\t", index=False)
    # Record the realized cohort used for annotations rather than assuming three sites.
    task = json.loads((root / "evaluation_tasks.json").read_text())[0]
    cfg = yaml.safe_load(Path(task["config"]).read_text())
    _json(
        root / "evaluation.accepted.json",
        {
            "selected": selection,
            "n_evaluation_loci": selected.locus_id.nunique(),
            "n_instances_per_sampling_seed": len(expected),
            "n_sampling_seeds": len(SEEDS),
            "cohort": cfg["sites"],
            "analyzed_n": sum(s["n"] for s in cfg["sites"].values()),
        },
    )

    render_headline(root)


def render_headline(root: Path) -> None:
    from .figures.federation import fed6_what_federation_buys
    from .figures.render import arm_cohort_sizes

    source = Path(json.loads((root / "design.json").read_text())["source"])
    headline = pd.read_csv(root / "fm_results_by_arm.tsv", sep="\t", low_memory=False)
    cohort = arm_cohort_sizes(source, headline)
    cfg = yaml.safe_load((source / "pipeline_config.yaml").read_text())
    selection = json.loads((root / "selection.json").read_text())
    cohort["federation"] = (
        sum(s["n"] for s in cfg["sites"].values()),
        len(cfg["superpopulations"]),
    )
    cohort[SMART_ARM] = (
        sum(selection["ancestry_counts"].values()),
        len(selection["ancestry_counts"]),
    )
    path = fed6_what_federation_buys(
        headline, root / "figures/federation/shared/fed6_what_federation_buys.png", cohort
    )
    manifest = root / "plots.json"
    record = json.loads(manifest.read_text()) if manifest.exists() else {"figures": {}}
    record["fed6_ready"] = True
    group = record["figures"].setdefault("federation/shared", [])
    relative = str(path.relative_to(root))
    if relative not in group:
        group.append(relative)
    _json(manifest, record)


def render(root: Path, include_eda: bool = True, out_dir: Path | None = None) -> None:
    """Regenerate every supported plot from saved inputs without changing the source run."""
    from appfl_bio_suite.core import plot_style

    from .figures import eda, federation, figstyle, paper_plots
    from .figures.render import arm_cohort_sizes

    source = Path(json.loads((root / "design.json").read_text())["source"])
    output = out_dir if out_dir is not None else root / "figures"
    output.mkdir(exist_ok=True)
    figstyle.apply_style()
    written = {}
    if include_eda:
        written["eda"] = eda.write_figures(source, output / "eda", strict=True)
    frames = {}
    for name, directory, filename in (
        ("centralized", "fine_mapping", "fm_results.tsv"),
        ("federated", "fed_fine_mapping", "fed_fm_results.tsv"),
    ):
        reports = source / "reports" / directory
        frame = pd.read_csv(reports / filename, sep="\t")
        frame = frame.loc[frame.causal_mode.eq("shared")].copy()
        frames[name] = frame
        for mode, subset in frame.groupby("causal_mode"):
            written[f"{name}/{mode}"] = paper_plots.write_figures(
                subset,
                output / name / mode,
                data_root=source,
                detail_dir=reports / "detail",
                strict=True,
            )
    headline = (
        pd.read_csv(root / "fm_results_by_arm.tsv", sep="\t", low_memory=False)
        if (root / "evaluation.accepted.json").exists()
        else None
    )
    cohort = arm_cohort_sizes(source, headline) if headline is not None else {}
    cfg = yaml.safe_load((source / "pipeline_config.yaml").read_text())
    cohort["federation"] = (
        sum(s["n"] for s in cfg["sites"].values()),
        len(cfg["superpopulations"]),
    )
    if headline is not None:
        selection = json.loads((root / "selection.json").read_text())
        cohort[SMART_ARM] = (
            sum(selection["ancestry_counts"].values()),
            len(selection["ancestry_counts"]),
        )
    for mode, central in frames["centralized"].groupby("causal_mode"):
        written[f"federation/{mode}"] = federation.write_figures(
            central,
            output / "federation" / mode,
            federated=frames["federated"].query("causal_mode == @mode"),
            data_root=source,
            by_arm=headline if mode == "shared" else None,
            cohort=cohort,
            strict=True,
        )
    central = frames["centralized"]
    written["diagnostics"] = [
        paper_plots.pap3_convergence(
            paper_plots.parse_architecture(central), output / "fit_status.png"
        )
    ]
    detail = pd.read_csv(source / "reports/fine_mapping/credible_sets.tsv", sep="\t")
    detail = detail.merge(
        central.loc[central.causal_mode.eq("shared"), KEY], on=KEY, validate="many_to_one"
    )
    written["diagnostics"].append(
        paper_plots.pap4_power_coverage(
            paper_plots.parse_architecture(central.loc[central.causal_mode.eq("shared")]),
            output / "coverage_shared.png",
            detail,
        )
    )
    from .fedfm.utils import read_bim

    first = next(iter(cfg["sites"]))
    bim = Path(cfg["paths"]["processed_dir"]) / first / f"{first}_chr{cfg['chromosome']}.bim"
    written["harmonization"] = [
        federation.fed4_harmonization(
            dict.fromkeys(cfg["sites"], 0),
            len(read_bim(bim)),
            output / "federation/harmonization.png",
            (0, len(central)),
        )
    ]
    _json(
        root / "plots.json",
        {
            "font_requested": "Helvetica",
            "font_rendered": plot_style.apply_style(),
            "palette_sha256": hashlib.sha256(plot_style.PALETTE_PATH.read_bytes()).hexdigest(),
            "fed6_ready": headline is not None,
            "causal_mode": "shared",
            "n_source_instances": len(central),
            "figures": {k: [str(p.relative_to(root)) for p in v] for k, v in written.items()},
        },
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage", choices=["prepare", "development", "select", "evaluation", "finish", "plots"]
    )
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--index", type=int, default=int(os.environ.get("PBS_ARRAY_INDEX", "0")))
    parser.add_argument("--development-reps", type=int, default=2)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.stage == "prepare":
        if args.source is None:
            parser.error("prepare requires --source")
        prepare(args.source, args.root, args.development_reps)
    elif args.stage in {"development", "evaluation"}:
        run_task(args.root, args.stage, args.index, args.workers)
    elif args.stage == "select":
        select(args.root)
    elif args.stage == "finish":
        finish(args.root)
    else:
        render(args.root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
