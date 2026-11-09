"""Watch the authorized recovery chain and retry interrupted batches conservatively.

No scientific failure is waived. A fit still nonconverged at 10,000 iterations,
a genuine selection/export failure, or three unsuccessful recovery rounds requires
inspection. Scheduler dependency cancellations are replaced once prerequisites pass.
"""

import argparse
import fcntl
import json
import re
import shutil
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path


def read(path):
    return json.loads(path.read_text())


def write(path, data):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=2))
    tmp.replace(path)


def scheduler(job):
    r = subprocess.run(["qstat", "-xf", "-F", "json", job], capture_output=True, text=True)
    r.check_returncode()
    jobs = json.loads(r.stdout)["Jobs"]
    if len(jobs) != 1:
        raise RuntimeError(f"Ambiguous scheduler response for {job}")
    return next(iter(jobs.values()))


def phase_finished(primary, supplementary):
    """Never retry a phase while any of its disjoint batches can still write."""
    return primary["job_state"] == "F" and all(job["job_state"] == "F" for job in supplementary)


def resolved_prerequisite_race(log, root):
    """Only retry a missing task acceptance file once that exact file exists."""
    last = log.read_text().rstrip().splitlines()[-1]
    match = re.fullmatch(
        r"FileNotFoundError: \[Errno 2\] No such file or directory: '([^']+)'", last
    )
    if not match:
        return False
    path = Path(match[1])
    try:
        relative = path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return (
        len(relative.parts) == 4
        and relative.parts[0] in {"benchmark", "development", "evaluation"}
        and relative.name == "accepted.json"
        and path.is_file()
    )


def pending(root, phase):
    path = root / f"{phase}_tasks.json"
    if not path.exists():
        return None
    return [
        t
        for t in read(path)
        if not (root / phase / t["name"] / t["locus_id"] / "accepted.json").exists()
    ]


def cap_failure(root, phase):
    files = [
        *root.glob(f"{phase}.failures.*.json"),
        *root.glob(f"recovery-*/**/{phase}.failures.*.json"),
    ]
    for file in files:
        for failure in read(file):
            if "Nonconverged fit: " not in failure["error"]:
                continue
            task = failure["task"]
            directory = root / phase / task["name"] / task["locus_id"]
            if (directory / "accepted.json").exists():
                continue
            for record in directory.glob("*/multisusie/fit.json"):
                f = read(record)
                if not f.get("converged") and f["options"]["max_iter"] >= 10000:
                    return str(record)
    return None


def submit(root, recovery, stage, round_no, tasks=None):
    old = read(recovery / f"submission.{stage}.json")
    cmd = old["command"][:]
    # Completed/aborted dependencies are replaced by explicit acceptance checks
    # in this manager, immediately before submission.
    if "-W" in cmd:
        i = cmd.index("-W")
        del cmd[i : i + 2]
    target = recovery / f"retry-{stage}-{round_no}"
    target.mkdir(exist_ok=False)
    (target / "logs").mkdir()
    if tasks is not None:
        for name in ["recover_multisusie.py", "job.pbs"]:
            shutil.copy2(recovery / name, target / name)
        size = max(1, (12 if stage == "benchmark" else 32) // (2**round_no))
        batches = [tasks[i : i + size] for i in range(0, len(tasks), size)]
        write(target / f"{stage}.batches.json", batches)
        # Evaluation uses the fixed 64-way stride in the recovery overlay.
        count = 64 if stage == "evaluation" else len(batches)
        cmd[cmd.index("-J") + 1] = f"0-{count - 1}"
        env = cmd[cmd.index("-v") + 1]
        env = env.replace(f"FM_RECOVERY={recovery}", f"FM_RECOVERY={target}")
        parts = dict(x.split("=", 1) for x in env.split(","))
        parts["FM_WORKERS"] = str(max(1, int(parts["FM_WORKERS"]) // (2**round_no)))
        cmd[cmd.index("-v") + 1] = ",".join(f"{k}={v}" for k, v in parts.items())
        cmd[-1] = str(target / "job.pbs")
    cmd[cmd.index("-o") + 1] = str(target / "logs")
    r = subprocess.run(cmd, capture_output=True, text=True)
    record = dict(command=cmd, returncode=r.returncode, stdout=r.stdout, stderr=r.stderr)
    write(target / "submission.json", record)
    r.check_returncode()
    return r.stdout.strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--recovery", type=Path, required=True)
    a = parser.parse_args()
    root = a.root.resolve()
    rec = a.recovery.resolve()
    lock = (rec / "manager.lock").open("w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    path = rec / "manager.status.json"
    state = (
        read(path)
        if path.exists()
        else dict(
            jobs={
                s: read(rec / f"submission.{s}.json")["stdout"].strip()
                for s in [
                    "diagnostic",
                    "benchmark",
                    "development",
                    "select",
                    "evaluation",
                    "finish",
                ]
            },
            rounds={s: 0 for s in ["benchmark", "development", "select", "evaluation", "finish"]},
        )
    )
    try:
        while True:
            state["updated_at_utc"] = datetime.now(UTC).isoformat()
            state["state"] = "monitoring"
            states = {s: scheduler(job) for s, job in state["jobs"].items()}
            state["scheduler"] = {
                s: dict(state=v.get("job_state"), exit_status=v.get("Exit_status"))
                for s, v in states.items()
            }
            supplementary = {
                phase: [scheduler(job) for job in jobs]
                for phase, jobs in state.get("supplementary_jobs", {}).items()
            }
            state["supplementary_scheduler"] = {
                phase: [
                    dict(state=v.get("job_state"), exit_status=v.get("Exit_status")) for v in jobs
                ]
                for phase, jobs in supplementary.items()
            }
            if not (rec / "diagnostics.accepted.json").exists():
                if states["diagnostic"]["job_state"] == "F":
                    raise RuntimeError(
                        "Diagnostic job ended without acceptance; inspect its fit traces"
                    )
                write(path, state)
                time.sleep(60)
                continue
            remaining = {p: pending(root, p) for p in ["benchmark", "development", "evaluation"]}
            state["remaining_loci"] = {
                p: len(t) if t is not None else None for p, t in remaining.items()
            }
            for phase, tasks in remaining.items():
                if not tasks or not phase_finished(states[phase], supplementary.get(phase, [])):
                    continue
                if phase == "evaluation" and not (root / "selection.lock.json").exists():
                    continue
                failure = cap_failure(root, phase)
                if failure:
                    raise RuntimeError(f"Fit remains nonconverged at extended ceiling: {failure}")
                if state["rounds"][phase] >= 3:
                    raise RuntimeError(f"{phase} still incomplete after three recovery rounds")
                n = state["rounds"][phase] + 1
                state["jobs"][phase] = submit(root, rec, phase, n, tasks)
                state["rounds"][phase] = n
                write(path, state)
            for stage, ready, artifact in [
                ("select", remaining["development"] == [], root / "selection.lock.json"),
                (
                    "finish",
                    remaining["benchmark"] == [] and remaining["evaluation"] == [],
                    root / "poster.accepted.json",
                ),
            ]:
                if not ready or artifact.exists() or states[stage]["job_state"] != "F":
                    continue
                job = state["jobs"][stage]
                log = root / "logs" / f"{stage}.{job}.log"
                if (
                    log.exists()
                    and log.stat().st_size
                    and not resolved_prerequisite_race(log, root)
                ):
                    raise RuntimeError(f"{stage} executed but did not pass; inspect {log}")
                if state["rounds"][stage] >= 3:
                    raise RuntimeError(f"{stage} repeatedly failed to execute")
                n = state["rounds"][stage] + 1
                state["jobs"][stage] = submit(root, rec, stage, n)
                state["rounds"][stage] = n
                write(path, state)
            if (root / "poster.accepted.json").exists():
                state["state"] = "artifacts_generated_visual_review_pending"
                write(path, state)
                return
            write(path, state)
            time.sleep(60)
    except Exception as exc:
        state.update(
            state="needs_inspection", error=repr(exc), updated_at_utc=datetime.now(UTC).isoformat()
        )
        write(path, state)
        raise


if __name__ == "__main__":
    main()
