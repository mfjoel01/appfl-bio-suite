"""Submit the bounded recovery chain, preserving every previous scheduler record."""

import argparse
import json
import os
import subprocess
from pathlib import Path


def submit(root, recovery, stage, after=None):
    record = recovery / f"submission.{stage}.json"
    if record.exists():
        return json.loads(record.read_text())["stdout"].strip()
    checkout = root.parents[2]
    python = checkout / "local/venv-multisusie/bin/python"
    cmd = [
        "qsub",
        "-A",
        os.environ["FM_ACCOUNT"],
        "-q",
        "single-gpu",
        "-N",
        f"msusie-r-{stage[:8]}",
        "-r",
        "y",
        "-m",
        "n",
        "-l",
        "select=1:ngpus=1:ncpus=32:mem=110gb",
        "-l",
        "walltime=12:00:00",
        "-l",
        "filesystems=home:grand",
        "-j",
        "oe",
        "-o",
        str(recovery / "logs"),
    ]
    if after:
        cmd += ["-W", f"depend=afterok:{after}"]
    if stage in ["select", "finish"]:
        cmd += [
            "-v",
            f"FM_ROOT={root},FM_PYTHON={python},FM_STAGE={stage},FM_CHECKOUT={checkout}",
            str(recovery / "finalize.pbs"),
        ]
    else:
        count = (
            64
            if stage == "evaluation"
            else len(json.loads((recovery / f"{stage}.batches.json").read_text()))
        )
        workers = 4 if stage in ["benchmark", "evaluation"] else 12
        cmd += [
            "-J",
            f"0-{count - 1}",
            "-v",
            f"FM_ROOT={root},FM_PYTHON={python},FM_RECOVERY={recovery},FM_PHASE={stage},FM_WORKERS={workers}",
            str(recovery / "job.pbs"),
        ]
    result = subprocess.run(cmd, text=True, capture_output=True)
    record.write_text(
        json.dumps(
            dict(
                command=cmd,
                returncode=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
            ),
            indent=2,
        )
    )
    result.check_returncode()
    job = result.stdout.strip()
    print(stage, job, flush=True)
    return job


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--recovery", type=Path, required=True)
    a = p.parse_args()
    root, rec = a.root.resolve(), a.recovery.resolve()
    diagnostic = json.loads((rec / "submission.diagnostic.json").read_text())["stdout"].strip()
    benchmark = submit(root, rec, "benchmark", diagnostic)
    development = submit(root, rec, "development", diagnostic)
    selected = submit(root, rec, "select", development)
    evaluated = submit(root, rec, "evaluation", selected)
    submit(root, rec, "finish", f"{benchmark}:{evaluated}")


if __name__ == "__main__":
    main()
