#!/usr/bin/env python3
"""Freeze and submit the fixed-N composition search with fresh validation draws."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))


def freeze(root: Path) -> None:
    if (root / "snapshot.json").exists():
        return
    shutil.copytree(REPO / "src", root / "code/src", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(REPO / "vendor/bin", root / "code/vendor/bin", symlinks=False)
    scripts = root / "code/scripts/fine-mapping"
    scripts.mkdir(parents=True)
    for name in ("composition.pbs", "submit_composition.py"):
        shutil.copy2(REPO / "scripts/fine-mapping" / name, scripts / name)
    files = [p for p in (root / "code").rglob("*") if p.is_file()]
    files += [
        root / name
        for name in ("design.json", "split.tsv", "development_tasks.json", "development_plan.json")
    ]
    files += list(root.glob("*/pipeline_config.yaml"))
    files += list(root.glob("*/loci/selected_loci.tsv"))
    files += list(root.glob("development/arms/*.yaml"))
    files += list(root.glob("development/arms/*_processed/*/*_manifest.tsv"))
    snapshot = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (root / "snapshot.json").write_text(json.dumps(snapshot, indent=2) + "\n")
    (root / "environment.json").write_text(
        json.dumps(
            dict(
                python=sys.executable,
                git_head=subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
                ).strip(),
                packages={
                    name: importlib.metadata.version(name)
                    for name in (
                        "numpy",
                        "pandas",
                        "scipy",
                        "joblib",
                        "pydantic",
                        "matplotlib",
                        "PyYAML",
                    )
                },
            ),
            indent=2,
        )
        + "\n"
    )
    (root / "working_tree.patch").write_bytes(
        subprocess.check_output(["git", "diff", "--binary", "HEAD"], cwd=REPO)
    )


def submit(root: Path, account: str, python: str, concurrency: int = 3) -> dict:
    if concurrency < 1:
        raise ValueError("Concurrency must be positive")
    path = root / "jobs.json"
    jobs = json.loads(path.read_text()) if path.exists() else {}
    design = json.loads((root / "design.json").read_text())
    count = len(json.loads((root / "development_tasks.json").read_text()))

    def add(stage, dependency=None, count=None, workers=32):
        if stage in jobs:
            return jobs[stage]["id"]
        walltime = "04:00:00"
        if stage.endswith("-check"):
            walltime = "02:00:00"
        elif stage == "select":
            walltime = "00:10:00"
        elif stage == "finish":
            walltime = "00:20:00"
        command = [
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
            "fm-comp-" + stage,
            "-l",
            "select=1:ngpus=8:ncpus=256:mem=960gb",
            "-l",
            "place=exclhost",
            "-l",
            "walltime=" + walltime,
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
            command += ["-W", "depend=afterok:" + dependency]
        if count:
            command += ["-J", f"0-{count - 1}%{concurrency}"]
        command += [str(root / "code/scripts/fine-mapping/composition.pbs")]
        response = subprocess.run(command, capture_output=True, text=True, check=True)
        jobs[stage] = dict(id=response.stdout.strip(), command=command, stderr=response.stderr)
        path.write_text(json.dumps(jobs, indent=2) + "\n")
        print(stage + ": " + jobs[stage]["id"], flush=True)
        return jobs[stage]["id"]

    sim = add("development-sim", count=design["simulation_shards"], workers=1)
    check = add("development-check", sim)
    development = add("development", check, count=count)
    selection = add("select", development, workers=1)
    sim = add("evaluation-sim", selection, count=design["simulation_shards"], workers=1)
    check = add("evaluation-check", sim)
    # Four fresh baselines plus at most three selected-cohort sampling draws.
    maximum = (4 + len(design["sampling_seeds"])) * design["evaluation_shards"]
    evaluation = add("evaluation", check, count=maximum)
    add("finish", evaluation, workers=1)
    return jobs


def main() -> None:
    from appfl_bio_suite.experiments.fine_mapping.composition import prepare

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--account", help="PBS project allocation, required with --submit")
    parser.add_argument("--development-reps", type=int, default=5)
    parser.add_argument("--evaluation-reps", type=int, default=10)
    parser.add_argument("--concurrency", type=int, default=3)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    if args.concurrency < 1 or (args.submit and not args.account):
        parser.error("Use positive --concurrency and supply --account when submitting")
    root = args.root.resolve()
    if not (root / "design.json").exists():
        if args.source is None:
            parser.error("A new run requires --source")
        prepare(args.source, root, args.development_reps, args.evaluation_reps)
    design = json.loads((root / "design.json").read_text())
    if design.get("protocol") != "refined-composition-v2":
        parser.error("Historical runs require their original frozen launcher; use a new root")
    freeze(root)
    if args.submit:
        submit(root, args.account, sys.executable, args.concurrency)
    else:
        print(f"Frozen experiment: {root}; pass --submit --account <allocation> to launch")


if __name__ == "__main__":
    main()
