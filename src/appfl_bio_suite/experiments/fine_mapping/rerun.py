"""Dependency-safe batch stages for the scientifically corrected fixed-panel rerun."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .fedfm.utils import load_config


def check_cohort(cfg) -> None:
    from .fedfm.utils import read_bim, read_fam, read_site_manifest

    reference = read_bim(cfg.resolved_path("hapnest_dir") / f"chr{cfg.chromosome}.bim").set_index(
        "snp_id"
    )
    enrolled = set()
    for site, spec in cfg.sites.items():
        prefix = cfg.site_dir(site) / f"{site}_chr{cfg.chromosome}"
        bim = read_bim(prefix.with_suffix(".bim")).set_index("snp_id")
        ref = reference.loc[bim.index]
        if not bim[["a1", "a2", "bp"]].equals(ref[["a1", "a2", "bp"]]):
            raise ValueError(f"{site}: phenotype-generation allele coding differs from reference")
        man = read_site_manifest(cfg.site_dir(site) / f"{site}_manifest.tsv")
        fam = read_fam(prefix.with_suffix(".fam"))
        ids = set(zip(man.FID, man.IID, strict=True))
        fam_ids = set(zip(fam.FID, fam.IID, strict=True))
        if len(ids) != spec.n or ids != fam_ids or ids & enrolled:
            raise ValueError(f"{site}: cohort size, identity or disjointness violation")
        if man.superpopulation.value_counts().to_dict() != dict(spec.composition):
            raise ValueError(f"{site}: ancestry composition mismatch")
        enrolled.update(ids)


def rescore(cfg, index: int, shards: int) -> None:
    from .fedfm.locus_selection import _load_site_handles, score_window

    directory = cfg.resolved_path("loci_dir")
    original = pd.read_csv(directory / "frozen_selected_loci.tsv", sep="\t")
    handles = _load_site_handles(cfg)
    rows = []
    for _, locus in original.iloc[index::shards].iterrows():
        result = score_window(handles, locus, cfg, directory / "ld_matrices")
        if result["status"] != "ok":
            raise ValueError(f"Selected window cannot be re-scored: {locus.locus_id}")
        result["locus_id"] = locus.locus_id
        result["legacy_stratum"] = locus.stratum
        rows.append(result)
    pd.DataFrame(rows).to_csv(
        directory / f"rescore.part{index:03d}of{shards:03d}.tsv", sep="\t", index=False
    )


def merge_rescore(cfg, shards: int) -> None:
    directory = cfg.resolved_path("loci_dir")
    original = pd.read_csv(directory / "frozen_selected_loci.tsv", sep="\t")
    frames = [
        pd.read_csv(directory / f"rescore.part{i:03d}of{shards:03d}.tsv", sep="\t")
        for i in range(shards)
    ]
    frame = pd.concat(frames, ignore_index=True)
    if frame.locus_id.duplicated().any() or set(frame.locus_id) != set(original.locus_id):
        raise ValueError("Incomplete or duplicated fixed-panel re-scoring")
    # This is a frozen, historically selected 79-locus panel. New strata describe
    # corrected divergence within that panel, not a new random genome-wide sample.
    ordered = frame.sort_values(["divergence_score", "locus_id"]).index
    labels = np.array_split(ordered.to_numpy(), 3)
    for label, indices in zip(["low", "medium", "high"], labels, strict=True):
        frame.loc[indices, "stratum"] = label
    frame["selection_design"] = "frozen_legacy_panel_harmonized_rescore"
    frame.sort_values("locus_id").to_csv(directory / "selected_loci.tsv", sep="\t", index=False)
    (directory / "rescore.accepted.json").write_text(
        json.dumps(
            {
                "n_loci": len(frame),
                "pass": True,
                "design": "frozen panel; no genome-wide representativeness claim",
            },
            indent=2,
        )
    )


def bundle(cfg, root: Path) -> None:
    from .simulation import load_scenario, register_bundles
    from .simulation.bundler import bundle_sites

    check_cohort(cfg)
    validation = json.loads(
        (cfg.resolved_path("reports_dir") / "validation/summary.json").read_text()
    )
    checks = [v["pass"] for v in validation.values() if isinstance(v, dict) and "pass" in v]
    if not checks or not all(checks) or validation.get("sampled") != "all":
        raise ValueError("Bundles require successful full-cohort phenotype validation")
    if not (cfg.resolved_path("loci_dir") / "rescore.accepted.json").exists():
        raise ValueError("Bundles require completed harmonized locus re-scoring")
    destination = root / "bundles"
    paths = bundle_sites(
        cfg, list(cfg.sites), cfg.resolved_path("hapnest_dir") / f"chr{cfg.chromosome}", destination
    )
    (destination / "ground_truth").mkdir(exist_ok=True)
    shutil.copy2(
        cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv",
        destination / "ground_truth/causal_manifest.tsv",
    )
    scenario = load_scenario("three-site-hapnest")
    registry = register_bundles(destination, list(cfg.sites), scenario, "drs.alcf.anl.gov", paths)
    # A replacement configuration is published locally only after all content IDs
    # exist. Remote endpoint IDs and paths remain explicit deployment facts.
    template = root / "federation.pending.yaml"
    config = yaml.safe_load(template.read_text())
    config["ga4gh"]["drs"]["registry"] = str(destination / "drs_registry.json")
    exp = config["experiments"]["fine-mapping"]
    exp["causal_manifest"] = str(destination / "ground_truth/causal_manifest.tsv")
    for site in exp["sites"]:
        name = site["client_id"]
        obj = registry.by_name(name)
        site["drs_uri"] = obj.self_uri
        site["data_use_profile"] = str(paths[name] / "DATA_USE.json")
        if name == "anl":
            site["data_dir"] = str(paths[name])
        site["output_dir"] = (
            str(root / "deployment_output" / name) if name == "anl" else site["output_dir"]
        )
    ready = yaml.safe_dump(config, sort_keys=False)
    (root / "federation.ready.yaml").write_text(ready)
    deployment = json.loads((root / "deployment_target.json").read_text())
    target = Path(deployment["path"])
    if hashlib.sha256(target.read_bytes()).hexdigest() != deployment["expected_sha256"]:
        raise RuntimeError(
            "Deployment configuration changed during rerun; "
            "ready config retained without overwriting it"
        )
    temporary = target.with_suffix(".ready.tmp")
    temporary.write_text(ready)
    temporary.replace(target)
    (root / "bundles.accepted.json").write_text(
        json.dumps({"pass": True, "bundles": {s: str(p) for s, p in paths.items()}}, indent=2)
    )


def analyze(root: Path, index: int, workers: int) -> None:
    from .arms import run_arm
    from .fedfm.fed_fine_mapping import run_fed_fine_mapping
    from .fedfm.fine_mapping import run_fine_mapping

    task = json.loads((root / "analysis_tasks.json").read_text())[index]
    cfg = load_config(task["config"])
    arguments = dict(shard_index=task["shard"], n_shards=task["shards"], n_workers=workers)
    if task["kind"] == "federated":
        run_fed_fine_mapping(cfg, maf=cfg.fine_mapping.maf, **arguments)
    elif task["kind"] == "centralized":
        run_fine_mapping(cfg, maf=cfg.fine_mapping.maf, **arguments)
    else:
        run_arm(task["name"], Path(task["config"]), **arguments)


def finish(root: Path) -> None:
    from .arms import merge_arms
    from .fedfm.fed_fine_mapping import merge_fed_fm_results
    from .fedfm.fine_mapping import merge_fm_results
    from .inference import validate_result_grid
    from .reporting import (
        harvest_archives,
        paired_arm_differences,
        parity_check,
        scientific_summary,
    )

    tasks = json.loads((root / "analysis_tasks.json").read_text())
    runs = {t["name"]: t for t in tasks}
    frames = {}
    arm_paths = {}
    for name, task in runs.items():
        cfg = load_config(task["config"])
        fed = task["kind"] == "federated"
        (merge_fed_fm_results if fed else merge_fm_results)(cfg, task["shards"])
        directory = cfg.resolved_path("reports_dir") / (
            "fed_fine_mapping" if fed else "fine_mapping"
        )
        path = directory / ("fed_fm_results.tsv" if fed else "fm_results.tsv")
        frame = pd.read_csv(path, sep="\t")
        truth = pd.read_csv(cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv", sep="\t")
        expected = set(
            truth.loc[
                truth.replicate < cfg.architecture.replicates,
                ["locus_id", "architecture_id", "replicate"],
            ].itertuples(index=False, name=None)
        )
        validate_result_grid(frame, expected)
        frame = frame.merge(
            truth[["locus_id", "architecture_id", "replicate", "causal_mode"]],
            on=["locus_id", "architecture_id", "replicate"],
            validate="one_to_one",
        )
        frame.to_csv(path, sep="\t", index=False)
        scientific_summary(frame).to_csv(
            directory / "scientific_summary.tsv", sep="\t", index=False
        )
        harvest_archives(directory, frame, directory / "credible_sets.tsv")
        frames[name] = frame
        if not fed:
            arm_paths[name] = path
    parity_check(
        root / "reports/fine_mapping", root / "reports/fed_fine_mapping", root / "parity.json"
    )
    merged = merge_arms(arm_paths, root / "arms/fm_results_by_arm.tsv")
    paired_arm_differences(merged).to_csv(
        root / "arms/paired_differences.tsv", sep="\t", index=False
    )
    from .figures.paper_plots import pap3_convergence, pap4_power_coverage, parse_architecture

    figures = root / "figures"
    figures.mkdir(exist_ok=True)
    central = frames["federation"]
    detail = pd.read_csv(root / "reports/fine_mapping/credible_sets.tsv", sep="\t")
    shared_keys = central.loc[
        central.causal_mode == "shared", ["locus_id", "architecture_id", "replicate"]
    ]
    detail = detail.merge(
        shared_keys, on=["locus_id", "architecture_id", "replicate"], validate="many_to_one"
    )
    pap3_convergence(parse_architecture(central), figures / "fit_status.png")
    pap4_power_coverage(
        parse_architecture(central[central.causal_mode == "shared"]),
        figures / "coverage_shared.png",
        detail,
    )
    (root / "analysis.accepted.json").write_text(
        json.dumps(
            {
                "pass": True,
                "n_runs": len(runs),
                "interpretation": (
                    "fixed-panel composition/participation benchmark; inspect coverage "
                    "and nonconvergence before scientific sign-off"
                ),
            },
            indent=2,
        )
    )


def plots(root: Path) -> None:
    """Rebuild EDA, per-mode inference plots, exemplars and paired-arm figures."""
    import gzip
    import logging

    from .figures import eda, federation, paper_plots
    from .figures.detail import harvest_instance
    from .figures.render import arm_cohort_sizes

    if not (root / "analysis.accepted.json").exists():
        raise ValueError("Plot publication requires completed analysis acceptance")
    cfg = load_config(root / "pipeline_config.yaml")
    output = root / "figures"
    output.mkdir(exist_ok=True)
    log = logging.getLogger(__name__)
    written = {"eda": eda.write_figures(root, output / "eda", logger=log, strict=True)}
    frames = {}
    for name, subdirectory, filename in [
        ("centralized", "fine_mapping", "fm_results.tsv"),
        ("federated", "fed_fine_mapping", "fed_fm_results.tsv"),
    ]:
        directory = root / "reports" / subdirectory
        frame = pd.read_csv(directory / filename, sep="\t")
        frames[name] = frame
        detail = directory / "detail"
        detail.mkdir(exist_ok=True)
        paths = [detail / "fm_credible_sets.tsv", detail / "fm_cs_members.tsv"]
        for p in paths:
            p.unlink(missing_ok=True)
        for row in frame.itertuples(index=False):
            inst = f"{row.locus_id}_{row.architecture_id}_rep{row.replicate}"
            source = directory / "artifacts" / inst
            truth = json.loads((source / "record.json").read_text())["truth"]
            cs, members, _ = harvest_instance(
                source, list(cfg.superpopulations), truth, include_variants=False
            )
            for path, table in zip(paths, [cs, members], strict=True):
                if not table.empty:
                    table.to_csv(path, sep="\t", index=False, mode="a", header=not path.exists())
        # Valid zero-discovery runs still produce readable empty tables.
        for path, columns in zip(
            paths,
            [
                ["instance", "cs_id", "contains_causal", "cs_length"],
                ["instance", "cs_id", "snp", "cs_pip"],
            ],
            strict=True,
        ):
            if not path.exists():
                pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False)
        # Select deterministic illustrative fits in every available mode/stratum.
        # They are not an independent validation sample or an estimate of power.
        exemplars = (
            frame[frame.fit_status.eq("converged_cs")]
            .sort_values(["locus_id", "architecture_id", "replicate"])
            .groupby(["causal_mode", "stratum"], sort=True)
            .head(1)
        )
        for row in exemplars.itertuples(index=False):
            inst = f"{row.locus_id}_{row.architecture_id}_rep{row.replicate}"
            destination = detail / "work" / inst
            destination.mkdir(parents=True, exist_ok=True)
            for source in (directory / "artifacts" / inst).glob("*.gz"):
                with gzip.open(source, "rb") as src, (destination / source.stem).open("wb") as dst:
                    shutil.copyfileobj(src, dst)
        for mode, subset in frame.groupby("causal_mode"):
            key = f"{name}/{mode}"
            written[key] = paper_plots.write_figures(
                subset,
                output / name / mode,
                data_root=root,
                detail_dir=detail,
                logger=log,
                strict=True,
            )
    by_arm = pd.read_csv(root / "arms/fm_results_by_arm.tsv", sep="\t")
    cohort = arm_cohort_sizes(root, by_arm)
    # The full federation uses the base config instead of a redundant arm config.
    cohort["federation"] = (sum(s.n for s in cfg.sites.values()), len(cfg.superpopulations))
    for mode, central in frames["centralized"].groupby("causal_mode"):
        fed = frames["federated"].query("causal_mode == @mode")
        arms = by_arm.query("causal_mode == @mode")
        written[f"federation/{mode}"] = federation.write_figures(
            central,
            output / "federation" / mode,
            federated=fed,
            data_root=root,
            by_arm=arms,
            cohort=cohort,
            logger=log,
            strict=True,
        )
    # No new flips are expected after the reference-coding gate. Show that outcome
    # explicitly rather than borrowing recoding counts from an old run's log.
    from .fedfm.utils import read_bim

    first = next(iter(cfg.sites))
    n_variants = len(read_bim(cfg.site_dir(first) / f"{first}_chr{cfg.chromosome}.bim"))
    harmonization = federation.fed4_harmonization(
        {site: 0 for site in cfg.sites},
        n_variants,
        output / "federation/harmonization.png",
        (0, len(frames["centralized"])),
    )
    written["harmonization"] = [harmonization]
    record = {
        "pass": True,
        "figures": {
            group: [str(p.relative_to(root)) for p in paths] for group, paths in written.items()
        },
        "n_figures": sum(map(len, written.values())),
        "scope": "EDA; centralized/federated shared and divergent fits; all participation arms",
        "omissions": [
            "Genome-wide candidate-window plots are unavailable for a frozen selected-locus panel.",
            "Population/effect and exemplar panels require at least one retained credible set.",
        ],
    }
    (root / "plots.accepted.json").write_text(json.dumps(record, indent=2))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=[
            "rescore",
            "merge-rescore",
            "bundle",
            "analysis",
            "finish",
            "plots",
            "cohort-check",
        ],
    )
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--index", type=int, default=int(os.environ.get("PBS_ARRAY_INDEX", 0)))
    parser.add_argument("--shards", type=int, default=4)
    parser.add_argument("--workers", type=int, default=32)
    args = parser.parse_args(argv)
    import logging

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    cfg = load_config(args.root / "pipeline_config.yaml")
    if args.stage == "rescore":
        rescore(cfg, args.index, args.shards)
    elif args.stage == "merge-rescore":
        merge_rescore(cfg, args.shards)
    elif args.stage == "bundle":
        bundle(cfg, args.root)
    elif args.stage == "analysis":
        analyze(args.root, args.index, args.workers)
    elif args.stage == "finish":
        finish(args.root)
    elif args.stage == "plots":
        plots(args.root)
    else:
        check_cohort(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
