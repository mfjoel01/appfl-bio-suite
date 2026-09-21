#!/usr/bin/env python3
"""Freeze and submit the fixed-N composition search; never alter the source experiment."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
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
    files += [root / "design.json", root / "split.tsv", root / "development_tasks.json"]
    files += list(root.glob("*/pipeline_config.yaml"))
    files += list(root.glob("*/loci/selected_loci.tsv"))
    files += list(root.glob("*/ground_truth/causal_manifest.tsv"))
    files += list(root.glob("development/arms/*.yaml"))
    files += list(root.glob("development/arms/*_processed/*/*_manifest.tsv"))
    (root / "snapshot.json").write_text(
        json.dumps(
            {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            indent=2,
        )
    )
    packages = ["numpy", "pandas", "scipy", "joblib", "pydantic", "matplotlib", "PyYAML"]
    (root / "environment.json").write_text(
        json.dumps(
            {
                "python": sys.executable,
                "git_head": subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
                ).strip(),
                "packages": {name: importlib.metadata.version(name) for name in packages},
            },
            indent=2,
        )
    )
    (root / "working_tree.patch").write_bytes(
        subprocess.check_output(["git", "diff", "--binary"], cwd=REPO)
    )


def submit(root: Path, account: str, python: str, concurrency: int = 3) -> dict:
    path = root / "jobs.json"
    jobs = json.loads(path.read_text()) if path.exists() else {}

    def add(stage, dependency=None, array=None, workers=32):
        if stage in jobs:
            return jobs[stage]["id"]
        variables = f"FM_ROOT={root},FM_STAGE={stage},FM_PYTHON={python},FM_WORKERS={workers}"
        for name in ("BIOSIM_FONT_PATH", "BIOSIM_REQUIRE_HELVETICA"):
            if os.environ.get(name):
                variables += f",{name}={os.environ[name]}"
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
            f"fm-smart-{stage[:7]}",
            "-l",
            "select=1:ngpus=8:ncpus=256:mem=960gb",
            "-l",
            "place=exclhost",
            "-l",
            "walltime=04:00:00",
            "-l",
            "filesystems=home:grand",
            "-j",
            "oe",
            "-o",
            str(root / "logs"),
            "-v",
            variables,
        ]
        if dependency:
            command += ["-W", "depend=afterok:" + dependency]
        if array:
            command += ["-J", array]
        command += [str(root / "code/scripts/fine-mapping/composition.pbs")]
        try:
            response = subprocess.run(command, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            (root / "submission_error.json").write_text(
                json.dumps(
                    {"command": command, "stdout": exc.stdout, "stderr": exc.stderr}, indent=2
                )
            )
            print(exc.stderr, file=sys.stderr)
            raise
        jobs[stage] = {"id": response.stdout.strip(), "command": command, "stderr": response.stderr}
        path.write_text(json.dumps(jobs, indent=2))
        print(f"{stage}: {jobs[stage]['id']}", flush=True)
        return jobs[stage]["id"]

    count = len(json.loads((root / "development_tasks.json").read_text()))
    development = add("development", array=f"0-{count - 1}%{concurrency}")
    selection = add("select", dependency=development)
    evaluation = add("evaluation", dependency=selection, array=f"0-17%{concurrency}")
    plots = add("plots", workers=4)
    add("finish", dependency=f"{evaluation}:{plots}", workers=4)
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--account", default="GeomicVar")
    parser.add_argument("--concurrency", type=int, default=3)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    if not (root / "design.json").exists():
        if args.source is None:
            parser.error("A new experiment requires --source")
        from appfl_bio_suite.experiments.fine_mapping.composition import prepare

        prepare(args.source, root)
    if args.concurrency < 1:
        parser.error("--concurrency must be positive")
    freeze(root)
    if args.submit:
        submit(root, args.account, sys.executable, args.concurrency)
    else:
        print(f"Frozen experiment: {root}; pass --submit to launch")


if __name__ == "__main__":
    main()
