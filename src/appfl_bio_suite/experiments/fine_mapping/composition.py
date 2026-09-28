"""Covenant-centered fixed-N search with fresh development and validation phenotypes.

The candidate grid, simulation seeds, selection rule and region split are frozen before
execution. Validation phenotypes are generated only after the development choice is
frozen. Every validation arm uses those same new phenotypes and inference settings.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .arms import KEY, SITE_ORDER, Arm, build_configs, run_arm
from .fedfm.utils import load_config
from .inference import validate_result_grid
from .reporting import clustered_ratio, paired_arm_differences

SAMPLING_SEEDS = (20260928, 20260929, 20260930)
PHENOTYPE_SEEDS = {"development": 2026092801, "evaluation": 2026092802}
MIN_COVERAGE = 0.93
MAX_FDR = 0.05
SMART_ARM = "federation_smart_50k"


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


def read_json(path: Path):
    return json.loads(path.read_text())


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ancestry_totals(quotas: dict) -> dict[str, int]:
    totals = {}
    for counts in quotas.values():
        for pop, n in counts.items():
            totals[pop] = totals.get(pop, 0) + n
    return totals


def candidates(cfg: dict) -> list[dict]:
    """Twenty fixed, explicit site/ancestry allocations, including all solo controls."""
    base = {"covenant": {"AFR": 47500, "EUR": 1500, "CSA": 1000}}
    if cfg["sites"]["covenant"]["composition"] != base["covenant"]:
        raise ValueError("This search requires the documented Covenant composition")
    choices = []

    def add(name, quotas, solo=None):
        quotas = {s: {p: n for p, n in c.items() if n} for s, c in quotas.items()}
        totals = ancestry_totals(quotas)
        if sum(totals.values()) != 50000:
            raise ValueError(f"{name}: every candidate must contain exactly 50,000 participants")
        if any(n < cfg["fine_mapping"]["min_gwas_n"] for n in totals.values()):
            raise ValueError(f"{name}: ancestry below the analysis threshold")
        for site, counts in quotas.items():
            if any(
                n < 0 or n > cfg["sites"][site]["composition"].get(p, 0) for p, n in counts.items()
            ):
                raise ValueError(f"{name}: unavailable site/ancestry quota")
        choices.append(dict(name=name, site_counts=quotas, ancestry_counts=totals, solo=solo))

    def swap(name, donor_pop, n, recipient_site, recipient_pop):
        quotas = deepcopy(base)
        quotas["covenant"][donor_pop] -= n
        dest = quotas.setdefault(recipient_site, {})
        dest[recipient_pop] = dest.get(recipient_pop, 0) + n
        add(name, quotas)

    for site in SITE_ORDER:
        add(f"solo_{site}", {site: cfg["sites"][site]["composition"]}, solo=site)
    for n in (500, 1000, 2500):
        swap(f"afr_to_anl_eur_{n}", "AFR", n, "anl", "EUR")
    swap("csa_to_anl_eur_1000", "CSA", 1000, "anl", "EUR")
    swap("csa_to_anl_afr_1000", "CSA", 1000, "anl", "AFR")
    swap("eur_to_anl_afr_500", "EUR", 500, "anl", "AFR")
    for site in ("anl", "mbzuai"):
        swap(f"afr_to_{site}_csa_1000", "AFR", 1000, site, "CSA")
        for n in (1000, 2500):
            swap(f"afr_to_{site}_afr_{n}", "AFR", n, site, "AFR")
    swap("afr_to_anl_amr_1000", "AFR", 1000, "anl", "AMR")
    swap("afr_to_mbzuai_mid_1000", "AFR", 1000, "mbzuai", "MID")
    swap("csa_to_mbzuai_csa_1000", "CSA", 1000, "mbzuai", "CSA")
    add(
        "three_sites_eur500_mid1000",
        {
            "covenant": {"AFR": 46000, "EUR": 1500, "CSA": 1000},
            "anl": {"EUR": 500},
            "mbzuai": {"MID": 1000},
        },
    )
    add(
        "previous_winner_afr45_eur5",
        {
            "covenant": {"AFR": 45000, "EUR": 1500},
            "anl": {"EUR": 3500},
        },
    )
    return choices


def phase_config(root: Path, phase: str) -> Path:
    return root / phase / "pipeline_config.yaml"


def prepare(source: Path, root: Path, development_reps: int = 5, evaluation_reps: int = 10) -> None:
    from .arms import _load
    from .rerun import check_cohort

    source, root = source.resolve(), root.resolve()
    if root.exists():
        raise FileExistsError(f"Use a new experiment directory: {root}")
    if min(development_reps, evaluation_reps) < 1:
        raise ValueError("Both phases require positive replicate counts")
    cfg = _load(source / "pipeline_config.yaml")
    choices = candidates(cfg)
    if cfg["master_seed"] in PHENOTYPE_SEEDS.values():
        raise ValueError("Fresh phenotype seeds must differ from the source simulation")
    check_cohort(load_config(source / "pipeline_config.yaml"))
    loci = pd.read_csv(Path(cfg["paths"]["loci_dir"]) / "selected_loci.tsv", sep="\t")
    split = split_regions(loci)
    root.mkdir(parents=True)
    (root / "logs").mkdir()
    split.to_csv(root / "split.tsv", sep="\t", index=False)
    _json(
        root / "design.json",
        dict(
            protocol="refined-composition-v2",
            source=str(source),
            causal_mode="shared",
            sampling_seeds=list(SAMPLING_SEEDS),
            phenotype_seeds=PHENOTYPE_SEEDS,
            source_master_seed=cfg["master_seed"],
            candidates=choices,
            development_replicates=development_reps,
            evaluation_replicates=evaluation_reps,
            simulation_shards=4,
            development_shards=2,
            evaluation_shards=6,
            selection_metric="unconditional_power",
            calibration_screen=dict(min_cs_coverage=MIN_COVERAGE, max_fdr_pip95=MAX_FDR),
            tie_break="prefer_solo_then_candidate_name",
            uncertainty="paired_connected_region_bootstrap",
            validation="new causal variants, effects and noise, shared by every evaluation arm",
            interpretation=(
                "Previously inspected fixed genotype/locus panel; fresh simulation draws, "
                "not independent genome-wide replication."
            ),
        ),
    )
    for phase, reps in (("development", development_reps), ("evaluation", evaluation_reps)):
        folder = root / phase
        local = deepcopy(cfg)
        local["master_seed"] = PHENOTYPE_SEEDS[phase]
        local["architecture"]["replicates"] = reps
        local["architecture"].setdefault("ancestry_divergent_causal", {})["enabled"] = False
        for key in ("loci_dir", "ground_truth_dir", "reports_dir", "logs_dir"):
            path = folder / key.removesuffix("_dir")
            path.mkdir(parents=True)
            local["paths"][key] = str(path)
        (folder / "processed").symlink_to(cfg["paths"]["processed_dir"], target_is_directory=True)
        split.loc[split.split.eq(phase)].drop(columns=["region", "split"]).to_csv(
            folder / "loci/selected_loci.tsv", sep="\t", index=False
        )
        phase_config(root, phase).write_text(yaml.safe_dump(local, sort_keys=False))
    build_tasks(root, "development", choices)


def build_tasks(root: Path, phase: str, choices: list[dict]) -> None:
    design = read_json(root / "design.json")
    seeds = design["sampling_seeds"]
    arms, metadata = [], {}
    for choice in choices:
        for seed in seeds[:1] if choice["solo"] else seeds:
            name = f"{choice['name']}_seed{seed}"
            if choice["solo"]:
                arm = Arm(name, (choice["solo"],), seed=seed)
            else:
                arm = Arm(
                    name,
                    tuple(choice["site_counts"]),
                    seed=seed,
                    site_counts=tuple(
                        (s, tuple(c.items())) for s, c in choice["site_counts"].items()
                    ),
                )
            arms.append(arm)
            metadata[name] = dict(candidate=choice["name"], sampling_seed=seed, role="selected")
    if phase == "evaluation":
        for name in (*SITE_ORDER, "federation"):
            arms.append(Arm(name, SITE_ORDER if name == "federation" else (name,)))
            metadata[name] = dict(candidate=name, sampling_seed=seeds[0], role="baseline")
    configs = build_configs(phase_config(root, phase), root / phase / "arms", tuple(arms))
    shards = design[f"{phase}_shards"]
    tasks = [
        dict(name=name, config=str(path), **metadata[name], shard=i, shards=shards)
        for name, path in configs.items()
        for i in range(shards)
    ]
    _json(root / f"{phase}_tasks.json", tasks)
    files = [root / f"{phase}_tasks.json", *configs.values()]
    files += list((root / phase / "arms").glob("*_processed/*/*_manifest.tsv"))
    _json(root / f"{phase}_plan.json", {str(p.relative_to(root)): sha(p) for p in files})


def verify_plan(root: Path, phase: str) -> None:
    for name, expected in read_json(root / f"{phase}_plan.json").items():
        if sha(root / name) != expected:
            raise ValueError(f"Frozen {phase} allocation changed: {name}")


def verify_acceptance(root: Path, phase: str) -> None:
    accepted = read_json(root / phase / "phenotypes.accepted.json")
    for name, expected in accepted["files"].items():
        if sha(root / phase / name) != expected:
            raise ValueError(f"Accepted {phase} inputs changed: {name}")


def simulate(root: Path, phase: str, index: int) -> None:
    from .fedfm.phenotype_sim import run_phenotype_sim

    if phase == "evaluation":
        verify_selection(root)
    cfg = load_config(phase_config(root, phase))
    if cfg.architecture.ancestry_divergent_causal.enabled:
        raise ValueError("Only shared-causal simulations are allowed")
    if cfg.master_seed != read_json(root / "design.json")["phenotype_seeds"][phase]:
        raise ValueError("Phenotype seed differs from the frozen design")
    run_phenotype_sim(cfg, index, read_json(root / "design.json")["simulation_shards"])


def check_phenotypes(root: Path, phase: str, workers: int) -> None:
    from .fedfm.phenotype_sim import build_architecture_grid, merge_phenotype_manifests

    config = phase_config(root, phase)
    cfg = load_config(config)
    truth = cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv"
    if not truth.exists():
        merge_phenotype_manifests(cfg, read_json(root / "design.json")["simulation_shards"])
    manifest = pd.read_csv(truth, sep="\t")
    loci = pd.read_csv(cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t")
    architectures = build_architecture_grid(
        cfg.architecture.ncsl,
        cfg.architecture.h2,
        cfg.architecture.rg,
        cfg.architecture.factorial_mode,
    )
    expected = {
        (locus, a.architecture_id, r)
        for locus in loci.locus_id
        for a in architectures
        for r in range(cfg.architecture.replicates)
    }
    if (
        manifest.duplicated(KEY).any()
        or not manifest.causal_mode.eq("shared").all()
        or set(manifest[KEY].itertuples(index=False, name=None)) != expected
    ):
        raise ValueError("Fresh simulation must contain the complete shared-causal grid")
    subprocess.run(
        [
            sys.executable,
            "-m",
            "appfl_bio_suite.experiments.fine_mapping.fedfm.validation",
            "--config",
            str(config),
            "--n-workers",
            str(workers),
            "--skip",
            "ld",
        ],
        check=True,
    )
    summary_path = cfg.resolved_path("reports_dir") / "validation/summary.json"
    summary = read_json(summary_path)
    for key in ("rederivation", "structural", "h2_calibration", "effect_size_correlation"):
        if not summary[key]["pass"] or summary["sampled"] != "all":
            raise ValueError(f"Fresh phenotype validation failed: {key}")
    files = [config, truth, truth.with_name("seeds.tsv"), summary_path]
    _json(
        root / phase / "phenotypes.accepted.json",
        dict(
            pass_checks=True,
            n_instances=len(manifest),
            master_seed=cfg.master_seed,
            files={str(p.relative_to(root / phase)): sha(p) for p in files},
        ),
    )


def run_task(root: Path, phase: str, index: int, workers: int) -> None:
    verify_plan(root, phase)
    verify_acceptance(root, phase)
    if phase == "evaluation":
        verify_selection(root)
    tasks = read_json(root / f"{phase}_tasks.json")
    # A solo fallback needs one selected run rather than three identical draws.
    # PBS reserves the maximum array length before the winner is known.
    if index >= len(tasks):
        if phase == "evaluation":
            return
        raise IndexError(index)
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


def rank_candidates(frame: pd.DataFrame, choices: list[dict], screen: dict) -> pd.DataFrame:
    ranks = (
        frame.groupby("candidate")
        .agg(
            power=("any_causal_captured", "mean"),
            high_confidence_yield=("causal_pip_max", lambda x: (x > 0.95).mean()),
            n_instances=("candidate", "size"),
            n_loci=("locus_id", "nunique"),
        )
        .reset_index()
    )
    grouped = frame.groupby("candidate")
    totals = grouped[
        ["n_cs_containing_causal", "n_credible_sets", "n_true_pip95", "n_false_pip95"]
    ].sum()
    totals["cs_coverage"] = totals.n_cs_containing_causal / totals.n_credible_sets.replace(
        0, np.nan
    )
    discoveries = totals.n_true_pip95 + totals.n_false_pip95
    totals["fdr_pip95"] = totals.n_false_pip95 / discoveries.replace(0, np.nan)
    ranks = ranks.merge(totals[["cs_coverage", "fdr_pip95"]], on="candidate", validate="one_to_one")
    ranks["eligible"] = ranks.cs_coverage.ge(screen["min_cs_coverage"]) & ranks.fdr_pip95.fillna(
        0
    ).le(screen["max_fdr_pip95"])
    ranks["solo_fallback"] = ranks.candidate.isin([c["name"] for c in choices if c["solo"]])
    return ranks.sort_values(
        ["eligible", "power", "solo_fallback", "candidate"], ascending=[False, False, False, True]
    ).reset_index(drop=True)


def verify_selection(root: Path) -> None:
    lock = read_json(root / "selection.lock.json")
    for name, expected in lock.items():
        if sha(root / name) != expected:
            raise ValueError(f"Frozen selection changed: {name}")


def select(root: Path) -> None:
    if (root / "selection.json").exists():
        raise FileExistsError("Selection already frozen")
    verify_plan(root, "development")
    verify_acceptance(root, "development")
    design = read_json(root / "design.json")
    development = collect(root, "development")
    ranks = rank_candidates(development, design["candidates"], design["calibration_screen"])
    ranks.to_csv(root / "development_ranking.tsv", sep="\t", index=False)
    if not ranks.eligible.any():
        raise ValueError("No candidate passes the predeclared development calibration screen")
    winner = next(c for c in design["candidates"] if c["name"] == ranks.iloc[0].candidate)
    build_tasks(root, "evaluation", [winner])
    _json(
        root / "selection.json",
        {
            **winner,
            "development_power": float(ranks.iloc[0].power),
            "development_results_sha256": sha(root / "development_results.tsv"),
            "evaluation_used_for_selection": False,
        },
    )
    names = ["selection.json", "development_results.tsv", "evaluation_plan.json"]
    _json(root / "selection.lock.json", {name: sha(root / name) for name in names})


def metric_summary(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    for keys, group in frame.groupby(group_cols):
        group = group.copy()
        group["opportunities"] = 1
        group["hits"] = group.any_causal_captured.astype(int)
        group["confident"] = group.causal_pip_max.gt(0.95).astype(int)
        group["discoveries95"] = group.n_true_pip95 + group.n_false_pip95
        if not isinstance(keys, tuple):
            keys = (keys,)
        for metric, num, den in (
            ("power", "hits", "opportunities"),
            ("high_confidence_yield", "confident", "opportunities"),
            ("cs_coverage", "n_cs_containing_causal", "n_credible_sets"),
            ("fdr_pip95", "n_false_pip95", "discoveries95"),
        ):
            point, lo, hi = clustered_ratio(group, num, den)
            rows.append(
                dict(zip(group_cols, keys, strict=True))
                | dict(
                    metric=metric,
                    estimate=point,
                    ci_low=lo,
                    ci_high=hi,
                    n_loci=group.locus_id.nunique(),
                    n_regions=group.region.nunique(),
                    n_instances=group[KEY].drop_duplicates().shape[0],
                    nonconverged_fraction=float(group.fit_status.eq("nonconverged").mean()),
                )
            )
    return pd.DataFrame(rows)


def assemble_headline(root: Path, results: pd.DataFrame) -> pd.DataFrame:
    design, selection = read_json(root / "design.json"), read_json(root / "selection.json")
    expected = _expected(phase_config(root, "evaluation"))
    if set(results.candidate) != {selection["name"], *SITE_ORDER, "federation"}:
        raise ValueError("Evaluation must include the frozen winner and all four fresh baselines")
    frames = []
    for candidate, group in results.groupby("candidate"):
        selected = candidate == selection["name"]
        seeds = sorted(group.sampling_seed.unique())
        wanted = (
            design["sampling_seeds"]
            if selected and not selection["solo"]
            else design["sampling_seeds"][:1]
        )
        if seeds != sorted(wanted):
            raise ValueError(f"Wrong sampling seeds for {candidate}")
        for _, draw in group.groupby("sampling_seed"):
            validate_result_grid(draw, expected)
        arm = SMART_ARM if selected else candidate
        if len(wanted) == 1:
            frames.extend(group.assign(arm=arm, sampling_seed=s) for s in design["sampling_seeds"])
        else:
            frames.append(group.assign(arm=arm))
    headline = pd.concat(frames, ignore_index=True)
    headline["causal_mode"] = "shared"
    headline["evaluation_split"] = "evaluation"
    label = ", ".join(f"{n:,} {p}" for p, n in selection["ancestry_counts"].items())
    headline.loc[headline.arm.eq(SMART_ARM), "smart_composition"] = label + (
        f"; solo fallback: {selection['solo']}" if selection["solo"] else "; explicit site quotas"
    )
    split = pd.read_csv(root / "split.tsv", sep="\t")
    return headline.merge(split[["locus_id", "region"]], on="locus_id", validate="many_to_one")


def finish(root: Path) -> None:
    from .figures import federation, figstyle
    from .figures.paper_plots import parse_architecture

    verify_selection(root)
    verify_plan(root, "evaluation")
    verify_acceptance(root, "evaluation")
    results = collect(root, "evaluation")
    headline = assemble_headline(root, results)
    headline.to_csv(root / "fm_results_by_arm.tsv", sep="\t", index=False)
    summary = metric_summary(headline, ["arm"])
    summary.to_csv(root / "evaluation_summary.tsv", sep="\t", index=False)
    metric_summary(parse_architecture(headline), ["arm", "h2", "rg", "ncsl"]).to_csv(
        root / "evaluation_by_architecture.tsv", sep="\t", index=False
    )
    comparisons = paired_arm_differences(headline, reference=SMART_ARM)
    comparisons.to_csv(root / "paired_vs_smart.tsv", sep="\t", index=False)
    paired_arm_differences(headline).to_csv(root / "paired_vs_full.tsv", sep="\t", index=False)
    solo = comparisons.loc[comparisons.arm.isin(SITE_ORDER)]
    full = comparisons.loc[comparisons.arm.eq("federation")].iloc[0]
    selected = summary.loc[summary.arm.eq(SMART_ARM)].set_index("metric").estimate
    screen = read_json(root / "design.json")["calibration_screen"]
    calibrated = selected.cs_coverage >= screen["min_cs_coverage"] and (
        pd.isna(selected.fdr_pip95) or selected.fdr_pip95 <= screen["max_fdr_pip95"]
    )
    selection = read_json(root / "selection.json")
    decision = dict(
        selected=selection,
        desired_order_observed=bool(
            solo.power_difference.gt(0).all() and full.power_difference < 0
        ),
        desired_order_supported_by_95pct_intervals=bool(
            solo.ci_low.gt(0).all() and full.ci_high < 0
        ),
        selected_passes_calibration_screen=bool(calibrated),
        n_evaluation_instances=headline[KEY].drop_duplicates().shape[0],
        n_evaluation_loci=headline.locus_id.nunique(),
        n_evaluation_regions=headline.region.nunique(),
        inference=(
            "same centralized inference driver for all arms; "
            "full federation is the pooled comparator"
        ),
    )
    # Publication records completion independently of whether the desired order holds.
    _json(root / "evaluation.accepted.json", decision)
    cfg = load_config(phase_config(root, "evaluation"))
    cohort = {s: (spec.n, len(spec.composition)) for s, spec in cfg.sites.items()}
    cohort["federation"] = (sum(s.n for s in cfg.sites.values()), len(cfg.superpopulations))
    cohort[SMART_ARM] = (
        sum(selection["ancestry_counts"].values()),
        len(selection["ancestry_counts"]),
    )
    figstyle.apply_style()
    federation.fed6_what_federation_buys(
        headline, root / "figures/federation/shared/fed6_what_federation_buys.png", cohort
    )
    truth = pd.read_csv(cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv", sep="\t")
    counts, n_variants, affected = federation.audit_harmonization(cfg, truth)
    federation.fed4_harmonization(
        counts, n_variants, root / "figures/federation/harmonization.png", affected
    )
    _json(
        root / "plots.json",
        dict(
            causal_mode="shared",
            fresh_evaluation=True,
            fed6_ready=True,
            figures=[str(p.relative_to(root)) for p in sorted((root / "figures").rglob("*.png"))],
        ),
    )
    lines = [
        "# Refined composition evaluation",
        "",
        f"Selected: `{selection['name']}`",
        "",
        "| Arm | Power |",
        "| --- | ---: |",
    ]
    for row in summary.loc[summary.metric.eq("power")].itertuples():
        lines.append(f"| {row.arm} | {row.estimate:.2%} |")
    lines += [
        "",
        f"Desired ordering observed: {decision['desired_order_observed']}.",
        "Supported by paired 95% intervals: "
        f"{decision['desired_order_supported_by_95pct_intervals']}.",
        f"Selected arm passes calibration screen: {calibrated}.",
        "",
        "Intervals resample connected genomic regions. See evaluation_summary.tsv and "
        "evaluation_by_architecture.tsv for coverage, false discoveries and nonconvergence.",
    ]
    (root / "RESULTS.md").write_text("\n".join(lines) + "\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=[
            "development-sim",
            "development-check",
            "development",
            "select",
            "evaluation-sim",
            "evaluation-check",
            "evaluation",
            "finish",
        ],
    )
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--index", type=int, default=int(os.environ.get("PBS_ARRAY_INDEX", "0")))
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.stage.endswith("-sim"):
        simulate(args.root, args.stage.removesuffix("-sim"), args.index)
    elif args.stage.endswith("-check"):
        check_phenotypes(args.root, args.stage.removesuffix("-check"), args.workers)
    elif args.stage in ("development", "evaluation"):
        run_task(args.root, args.stage, args.index, args.workers)
    elif args.stage == "select":
        select(args.root)
    else:
        finish(args.root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
