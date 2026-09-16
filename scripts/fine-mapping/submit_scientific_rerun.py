#!/usr/bin/env python3
"""Freeze a corrected fixed-panel run and optionally submit its PBS dependency graph."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import io
import json
import re
import shutil
import subprocess
import sys
import tarfile
from datetime import UTC, datetime
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(root: Path, base_path: Path, python: str) -> None:
    from appfl_bio_suite.experiments.fine_mapping.arms import ARMS, SITE_ORDER, Arm, build_configs
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import load_config
    from appfl_bio_suite.experiments.fine_mapping.rerun import check_cohort

    if (root / "snapshot.json").exists():
        raise FileExistsError(f"Already prepared: {root}; use --submit-only to submit it")
    root.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(base_path.read_text())
    for key in ["loci_dir", "reports_dir", "logs_dir", "ground_truth_dir"]:
        destination = root / key.removesuffix("_dir")
        destination.mkdir(exist_ok=True)
        if key == "loci_dir":
            source = Path(cfg["paths"][key]) / "selected_loci.tsv"
            shutil.copy2(source, destination / "frozen_selected_loci.tsv")
            shutil.copy2(source, destination / "selected_loci.tsv")
        cfg["paths"][key] = str(destination)
    (root / "processed").symlink_to(Path(cfg["paths"]["processed_dir"]), target_is_directory=True)
    cfg["fine_mapping"].update(
        maf=0.005, keep_ambiguous=True, n_signals=10, max_iter=1000, tol=1e-6
    )
    cfg["architecture"]["replicates"] = 10
    base = root / "pipeline_config.yaml"
    base.write_text(yaml.safe_dump(cfg, sort_keys=False))
    check_cohort(load_config(base))
    arms = tuple(a for a in ARMS if a.name != "federation") + (
        Arm("federation_50k_seed2", SITE_ORDER, fraction=1 / 3, seed=20260602),
        Arm("federation_50k_seed3", SITE_ORDER, fraction=1 / 3, seed=20260603),
    )
    paths = build_configs(base, root / "arms", arms=arms, replicates=10)
    runs = [("federation", "centralized", base), ("federated", "federated", base)]
    runs.extend((name, "arm", path) for name, path in paths.items())
    tasks = [
        dict(name=name, kind=kind, config=str(path), shard=i, shards=6)
        for name, kind, path in runs
        for i in range(6)
    ]
    (root / "analysis_tasks.json").write_text(json.dumps(tasks, indent=2))

    code = root / "code"
    code.mkdir()
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    source_paths = [
        "src",
        "scripts/fine-mapping/run_stage.py",
        "scripts/fine-mapping/submit_scientific_rerun.py",
        "scripts/fine-mapping/scientific_rerun.pbs",
    ]
    archive = subprocess.check_output(["git", "archive", revision, *source_paths], cwd=REPO)
    with tarfile.open(fileobj=io.BytesIO(archive)) as files:
        files.extractall(code, filter="data")
    shutil.copytree(REPO / "vendor/bin", code / "vendor/bin", symlinks=False)
    deployment = REPO / "local/federation.yaml"
    before = deployment.read_bytes()
    (root / "federation.before.yaml").write_bytes(before)
    pending = yaml.safe_load(before)
    exp = pending["experiments"]["fine-mapping"]
    exp.update(
        enabled=True,
        keep_ambiguous=True,
        n_signals=10,
        max_iter=1000,
        tol=1e-6,
        pops=cfg["superpopulations"],
        susiex_binary=str(code / "vendor/bin/SuSiEx"),
    )
    (root / "federation.pending.yaml").write_text(yaml.safe_dump(pending, sort_keys=False))
    # Prevent an accidental new run from using the known-invalid legacy bundles.
    disabled = yaml.safe_load(before)
    disabled["experiments"]["fine-mapping"]["enabled"] = False
    temporary = deployment.with_suffix(".disabled.tmp")
    temporary.write_text(yaml.safe_dump(disabled, sort_keys=False))
    temporary.replace(deployment)
    (root / "deployment_target.json").write_text(
        json.dumps(dict(path=str(deployment), expected_sha256=sha(deployment)), indent=2)
    )
    packages = ["numpy", "pandas", "scipy", "joblib", "pydantic", "matplotlib", "PyYAML"]
    (root / "environment.json").write_text(
        json.dumps(
            {
                "python": python,
                "python_version": sys.version,
                "packages": {p: importlib.metadata.version(p) for p in packages},
                "git_head": revision,
                "design": (
                    "frozen 79-locus panel, regenerated effects/phenotypes, "
                    "10 replicates, three matched-N seeds"
                ),
            },
            indent=2,
        )
    )
    subprocess.run(
        ["git", "diff", "--binary"],
        cwd=REPO,
        stdout=(root / "working_tree.patch").open("w"),
        check=True,
    )
    frozen = list(p for p in code.rglob("*") if p.is_file())
    frozen += [
        base,
        root / "analysis_tasks.json",
        root / "federation.pending.yaml",
        root / "deployment_target.json",
        root / "environment.json",
    ]
    frozen += list((root / "arms").glob("*.yaml"))
    frozen += list((root / "arms").glob("*_processed/*/*_manifest.tsv"))
    frozen += [root / "loci/frozen_selected_loci.tsv"]
    (root / "snapshot.json").write_text(
        json.dumps({str(p.relative_to(root)): sha(p) for p in sorted(frozen)}, indent=2)
    )
    print(f"Prepared {len(runs)} runs / {len(tasks)} analysis tasks at {root}", flush=True)


def submit(root: Path, account: str, python: str) -> None:
    manifest = root / "jobs.json"
    jobs = json.loads(manifest.read_text()) if manifest.exists() else {}
    script = root / "code/scripts/fine-mapping/scientific_rerun.pbs"
    manifest.write_text(json.dumps(jobs, indent=2))

    def add(stage, dependency=None, array=None, workers=8):
        if stage in jobs:
            print(f"{stage}: retaining {jobs[stage]['id']}", flush=True)
            return jobs[stage]["id"]
        args = [
            "qsub",
            "-r",
            "y",
            "-W",
            "umask=0022",
            "-A",
            account,
            "-q",
            "single-node",
            "-N",
            "fm-" + stage[:11],
            "-l",
            "select=1:ngpus=8:ncpus=256:mem=960gb",
            "-l",
            "walltime=04:00:00",
            "-l",
            "place=exclhost",
            "-l",
            "filesystems=home:grand",
            "-j",
            "oe",
            "-o",
            str(root / "logs"),
            "-v",
            f"FM_ROOT={root},FM_STAGE={stage},FM_PYTHON={python},FM_WORKERS={workers}",
        ]
        if dependency:
            args += ["-W", "depend=" + dependency]
        if array:
            args += ["-J", array]
        args.append(str(script))
        result = subprocess.run(args, text=True, capture_output=True)
        if result.returncode:
            print(result.stderr, file=sys.stderr, flush=True)
            (root / "submission_error.json").write_text(
                json.dumps(dict(command=args, stdout=result.stdout, stderr=result.stderr), indent=2)
            )
            result.check_returncode()
        jid = result.stdout.strip()
        jobs[stage] = {"id": jid, "command": args, "stderr": result.stderr}
        manifest.write_text(json.dumps(jobs, indent=2))
        print(f"{stage}: {jid}", flush=True)
        return jid

    a = add("rescore", array="0-3%2")
    b = add("merge-rescore", dependency=f"afterok:{a}")
    c = add("phenotypes", dependency=f"afterok:{b}", array="0-3%2")
    d = add("merge-phenotypes", dependency=f"afterok:{c}")
    e = add("validation", dependency=f"afterok:{d}", array="0-3%2")
    f = add("merge-validation", dependency=f"afterok:{e}")
    g = add("bundle", dependency=f"afterok:{f}")
    n_tasks = len(json.loads((root / "analysis_tasks.json").read_text()))
    h = add("analysis", dependency=f"afterok:{g}", array=f"0-{n_tasks - 1}%3", workers=32)
    i = add("finish", dependency=f"afterok:{h}")
    add("plots", dependency=f"afterok:{i}")


def archive_terminal_attempt(root: Path) -> Path:
    """Preserve the old graph and scheduler evidence before an explicit retry.

    Reject retries while any recorded job is still queued, running or held.
    Scientific source, data, completed parts and validated LD caches stay intact.
    """
    manifest = root / "jobs.json"
    jobs = json.loads(manifest.read_text())
    if not jobs:
        raise ValueError("No submitted jobs to retry; use --submit-only")
    history = {}
    for stage, job in jobs.items():
        output = subprocess.check_output(["qstat", "-xf", job["id"]], text=True)
        state = re.search(r"^\s*job_state\s*=\s*(\w+)", output, re.MULTILINE)
        if not state or state.group(1) != "F":
            raise RuntimeError(f"Cannot retry while {stage} ({job['id']}) is not terminal")
        history[stage] = output
    destination = root / "attempts" / datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    destination.mkdir(parents=True)
    (destination / "scheduler.json").write_text(json.dumps(history, indent=2))
    manifest.replace(destination / "jobs.json")
    print(f"Archived terminal job graph at {destination}", flush=True)
    return destination


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--base-config", type=Path)
    p.add_argument("--account", help="PBS project allocation (required when submitting)")
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--submit", action="store_true")
    p.add_argument("--submit-only", action="store_true")
    p.add_argument(
        "--retry-failed",
        action="store_true",
        help="archive a terminal job graph and resubmit using existing data/caches",
    )
    args = p.parse_args()
    if (args.submit or args.submit_only or args.retry_failed) and not args.account:
        p.error("--account is required when submitting PBS jobs")
    root = args.root.resolve()
    if args.retry_failed:
        archive_terminal_attempt(root)
    if not (args.submit_only or args.retry_failed):
        if not args.base_config:
            p.error("--base-config is required for preparation")
        prepare(root, args.base_config.resolve(), args.python)
    if args.submit or args.submit_only or args.retry_failed:
        submit(root, args.account, args.python)


if __name__ == "__main__":
    main()
