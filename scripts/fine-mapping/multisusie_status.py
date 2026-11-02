"""Summarize durable progress and scheduler state without reading large posteriors."""

import argparse
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--root", type=Path, required=True)
a = p.parse_args()
root = a.root.resolve()
status = {"captured_at_utc": datetime.now(UTC).isoformat()}
for phase in ["benchmark", "development", "evaluation"]:
    taskfile = root / f"{phase}_tasks.json"
    if not taskfile.exists():
        status[phase] = {"state": "waiting for development selection"}
        continue
    tasks = json.loads(taskfile.read_text())
    completed = 0
    accepted = 0
    failures = []
    for task in tasks:
        out = root / phase / task["name"] / task["locus_id"]
        if (out / "results.tsv").exists():
            try:
                completed += len(pd.read_csv(out / "results.tsv", sep="\t"))
            except pd.errors.EmptyDataError:
                pass  # writer may be between truncate and flush
        accepted += (out / "accepted.json").exists()
    for path in [
        *root.glob(f"{phase}.failures.*.json"),
        *root.glob(f"recovery-*/**/{phase}.failures.*.json"),
    ]:
        failures += json.loads(path.read_text())
    historical_count = len(failures)
    unresolved = {}
    for failure in failures:
        task = failure["task"]
        if not (root / phase / task["name"] / task["locus_id"] / "accepted.json").exists():
            unresolved[(task["name"], task["locus_id"])] = failure
    failures = list(unresolved.values())
    status[phase] = dict(
        expected_fits=sum(t["n_instances"] for t in tasks),
        completed_rows=completed,
        accepted_loci=accepted,
        total_loci=len(tasks),
        failures=failures,
        historical_failure_records=historical_count,
    )
r = subprocess.run(
    ["qstat", "-t", "-u", __import__("getpass").getuser()], text=True, capture_output=True
)
status["scheduler"] = r.stdout
status["scheduler_returncode"] = r.returncode
status["recovery_managers"] = {
    str(path.parent.relative_to(root)): json.loads(path.read_text())
    for path in sorted(root.glob("recovery-*/manager.status.json"))
}
status["poster_accepted"] = (root / "poster.accepted.json").exists()
completion = root / "completion.json"
status["completion"] = json.loads(completion.read_text()) if completion.exists() else None
(root / "status.latest.json").write_text(json.dumps(status, indent=2))
print(json.dumps(status, indent=2))
