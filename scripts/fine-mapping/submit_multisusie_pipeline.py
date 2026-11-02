"""Submit frozen full-grid phases. Record every scheduler response immediately."""

import argparse
import json
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--root", type=Path, required=True)
p.add_argument("--python", type=Path, required=True)
p.add_argument("--phase", choices=["benchmark", "development", "evaluation"], required=True)
p.add_argument("--shards", type=int, default=8)
p.add_argument("--workers", type=int, default=8)
p.add_argument("--after")
p.add_argument("--account", required=True)
a = p.parse_args()
root = a.root.resolve()
record = root / f"submission.{a.phase}.json"
if record.exists():
    raise FileExistsError(f"Submission already recorded: {record}")
cmd = [
    "qsub",
    "-A",
    a.account,
    "-q",
    "single-gpu",
    "-N",
    f"msusie-{a.phase[:8]}",
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
    str(root / "logs"),
    "-J",
    f"0-{a.shards - 1}",
    "-v",
    f"FM_ROOT={root},FM_PYTHON={a.python.resolve()},FM_PHASE={a.phase},FM_WORKERS={a.workers},FM_SHARDS={a.shards}",
    str(root / "job.pbs"),
]
# Keep the virtual-environment entry point: resolving its symlink loses its packages.
cmd[cmd.index("-v") + 1] = cmd[cmd.index("-v") + 1].replace(
    str(a.python.resolve()), str(a.python.absolute())
)
if a.after:
    cmd[1:1] = ["-W", f"depend=afterok:{a.after}"]
r = subprocess.run(cmd, text=True, capture_output=True)
record.write_text(
    json.dumps(
        dict(command=cmd, returncode=r.returncode, stdout=r.stdout, stderr=r.stderr), indent=2
    )
)
print(r.stdout, r.stderr)
r.check_returncode()
