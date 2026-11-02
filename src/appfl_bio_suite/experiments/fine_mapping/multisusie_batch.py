"""Freeze and execute a bounded comparison against the completed September run."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import shutil
from pathlib import Path

import pandas as pd
import yaml

from .fedfm.utils import load_config
from .multisusie import MultiSuSiEOptions
from .multisusie_compare import run_comparison

ARCHITECTURES = [
    "ncsl1_h2-0.0005_rg1",
    "ncsl2_h2-0.001_rg0.5",
    "ncsl2_h2-0.001_rg1",
    "ncsl3_h2-0.005_rg1",
]
METHODS = {"multisusie", "susiex", "susiex_unfiltered"}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def select_panel(loci, manifest):
    """One median-divergence locus per stratum; fixed architectures and replicate 0."""
    selected = []
    for stratum in ("low", "medium", "high"):
        group = loci[loci.stratum == stratum].sort_values(["divergence_score", "locus_id"])
        if group.empty:
            raise ValueError(f"No loci in stratum {stratum}")
        selected.append(group.iloc[len(group) // 2].locus_id)
    panel = pd.DataFrame(
        [
            {"locus_id": locus, "architecture_id": arch, "replicate": 0}
            for locus in selected
            for arch in ARCHITECTURES
        ]
    )
    merged = panel.merge(manifest, how="left", validate="one_to_one", indicator=True)
    if not merged._merge.eq("both").all():
        raise ValueError("Prespecified panel is absent from causal manifest")
    return panel


def prepare(base: Path, root: Path, checkout: Path, susiex: Path):
    """Snapshot source/configuration and fingerprint every selected phenotype input."""
    root = root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    config_paths = {
        "federation": base / "pipeline_config.yaml",
        "federation_50k": base / "arms/federation_50k.yaml",
        "covenant": base / "arms/covenant.yaml",
    }
    full = load_config(config_paths["federation"])
    loci_path = full.resolved_path("loci_dir") / "selected_loci.tsv"
    manifest_path = full.resolved_path("ground_truth_dir") / "causal_manifest.tsv"
    panel = select_panel(pd.read_csv(loci_path, sep="\t"), pd.read_csv(manifest_path, sep="\t"))
    panel.to_csv(root / "panel.tsv", sep="\t", index=False)
    shutil.copytree(
        checkout / "src", root / "code/src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc")
    )
    (root / "code/bin").mkdir()
    shutil.copy2(susiex, root / "code/bin/SuSiEx")
    shutil.copy2(checkout / "scripts/fine-mapping/multisusie_comparison.pbs", root / "job.pbs")
    (root / "configs").mkdir()
    (root / "logs").mkdir()
    hashes = {str(loci_path): digest(loci_path), str(manifest_path): digest(manifest_path)}
    bed_identity = {}
    tasks = []
    for name, path in config_paths.items():
        cfg = load_config(path)
        raw = yaml.safe_load(path.read_text())
        raw["paths"] = {key: str(cfg.resolved_path(key)) for key in raw["paths"]}
        # All new artifacts stay under this run, including any future reused utilities.
        raw["paths"]["reports_dir"] = str(root / name)
        raw["paths"]["logs_dir"] = str(root / "logs")
        frozen = root / "configs" / f"{name}.yaml"
        frozen.write_text(yaml.safe_dump(raw, sort_keys=False))
        tasks.append({"arm": name, "config": str(frozen), "central_check": name == "federation"})
        prefixes = [cfg.resolved_path("hapnest_dir") / f"chr{cfg.chromosome}"]
        prefixes += [cfg.site_dir(site) / f"{site}_chr{cfg.chromosome}" for site in cfg.sites]
        for prefix in prefixes:
            for suffix in (".bed", ".bim", ".fam"):
                source = prefix.with_suffix(suffix)
                stat = source.stat()
                if suffix == ".bed":
                    bed_identity[str(source)] = [stat.st_size, stat.st_mtime_ns]
                else:
                    hashes[str(source)] = digest(source)
        for site in cfg.sites:
            source = cfg.site_dir(site) / f"{site}_manifest.tsv"
            hashes[str(source)] = digest(source)
            for row in panel.itertuples():
                inst = f"{row.locus_id}_{row.architecture_id}_rep{row.replicate}"
                source = (
                    cfg.resolved_path("ground_truth_dir") / "phenotypes" / site / f"{inst}.pheno"
                )
                hashes[str(source)] = digest(source)
    (root / "tasks.json").write_text(json.dumps(tasks, indent=2))
    for path in root.rglob("*"):
        if path.is_file():
            hashes[str(path)] = digest(path)
    packages = {
        name: importlib.metadata.version(name)
        for name in ("MultiSuSiE", "numpy", "scipy", "pandas", "numba", "llvmlite")
    }
    from MultiSuSiE import susiepy_ss

    hashes[str(Path(susiepy_ss.__file__).resolve())] = digest(susiepy_ss.__file__)
    (root / "snapshot.json").write_text(
        json.dumps(
            {
                "sha256": hashes,
                "bed_stat_identity": bed_identity,
                "packages": packages,
                "design": "three median-divergence loci; four fixed architectures; replicate zero",
                "n_instances_per_arm": len(panel),
                "n_arms": len(tasks),
                "options": asdict_options(),
            },
            indent=2,
        )
    )
    return root


def asdict_options():
    from dataclasses import asdict

    return asdict(MultiSuSiEOptions())


def verify(root):
    snapshot = json.loads((root / "snapshot.json").read_text())
    for path, expected in snapshot["sha256"].items():
        if digest(path) != expected:
            raise RuntimeError(f"Frozen input changed: {path}")
    for path, expected in snapshot["bed_stat_identity"].items():
        stat = Path(path).stat()
        if [stat.st_size, stat.st_mtime_ns] != expected:
            raise RuntimeError(f"Genotype file identity changed: {path}")
    for name, expected in snapshot["packages"].items():
        if importlib.metadata.version(name) != expected:
            raise RuntimeError(f"Package changed: {name}")


def finish(root):
    tasks = json.loads((root / "tasks.json").read_text())
    panel = pd.read_csv(root / "panel.tsv", sep="\t")
    keys = ["locus_id", "architecture_id", "replicate"]
    expected = {
        (row.locus_id, row.architecture_id, row.replicate, method)
        for row in panel.itertuples()
        for method in METHODS
    }
    frames = []
    for task in tasks:
        directory = root / task["arm"]
        accepted = json.loads((directory / "accepted.json").read_text())
        if not accepted["pass"] or accepted["centralized_parity_checked"] != task["central_check"]:
            raise RuntimeError("Arm acceptance gate failed")
        frame = pd.read_csv(directory / "results.tsv", sep="\t")
        if (
            frame.duplicated(keys + ["method"]).any()
            or set(frame[keys + ["method"]].itertuples(index=False, name=None)) != expected
            or not frame.fit_status.str.startswith("converged_").all()
        ):
            raise RuntimeError("Comparison grid is incomplete, duplicated or contains failed fits")
        if task["central_check"]:
            parity = json.loads((directory / "parity.json").read_text())
            expected_instances = {
                f"{r.locus_id}_{r.architecture_id}_rep{r.replicate}" for r in panel.itertuples()
            }
            if (
                len(parity) != len(panel)
                or {p["instance"] for p in parity} != expected_instances
                or not all(p["pass"] for p in parity)
            ):
                raise RuntimeError("Incomplete centralized parity evidence")
        frames.append(frame.assign(arm=task["arm"]))
    frame = pd.concat(frames, ignore_index=True)
    frame.to_csv(root / "comparison.tsv", sep="\t", index=False)
    summaries = []
    for (arm, method), group in frame.groupby(["arm", "method"]):
        sums = group.select_dtypes("number").sum()
        summaries.append(
            {
                "arm": arm,
                "method": method,
                "n_instances": len(group),
                "causal_recall": sums.n_causal_captured / sums.n_causal,
                "cs_coverage": sums.n_cs_containing_causal / sums.n_credible_sets
                if sums.n_credible_sets
                else None,
                "mean_cs_size": sums.total_cs_snps / sums.n_credible_sets
                if sums.n_credible_sets
                else None,
                "true_pip95": sums.n_true_pip95,
                "false_pip95": sums.n_false_pip95,
                "runtime_s": sums.runtime_s,
            }
        )
    pd.DataFrame(summaries).to_csv(root / "summary.tsv", sep="\t", index=False)
    (root / "accepted.json").write_text(
        json.dumps(
            {
                "pass": True,
                "n_method_fits": len(frame),
                "scope": "bounded synthetic pilot; not population-level calibration or superiority",
            },
            indent=2,
        )
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "run", "finish"])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--base", type=Path)
    parser.add_argument("--checkout", type=Path)
    parser.add_argument("--susiex", type=Path)
    parser.add_argument("--index", type=int, default=0)
    args = parser.parse_args(argv)
    if args.stage == "prepare":
        if not all([args.base, args.checkout, args.susiex]):
            parser.error("prepare requires --base, --checkout and --susiex")
        prepare(args.base.resolve(), args.root, args.checkout.resolve(), args.susiex.resolve())
    else:
        root = args.root.resolve()
        verify(root)
        if args.stage == "finish":
            finish(root)
        else:
            task = json.loads((root / "tasks.json").read_text())[args.index]
            run_comparison(
                Path(task["config"]),
                root / "panel.tsv",
                root / task["arm"],
                root / "code/bin/SuSiEx",
                MultiSuSiEOptions(),
                task["central_check"],
            )


if __name__ == "__main__":
    main()
