"""Monitor paired-grid acceptance and resume interrupted PBS tasks only."""

import argparse
import fcntl
import hashlib
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


def read(p):
    return json.loads(p.read_text())


def write(p, v):
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(v, indent=2))
    tmp.replace(p)


def submit(root, stage, indices=None, dependency=None, attempt=0):
    design = read(root / "scheduler.json")
    cmd = [
        "qsub",
        "-A",
        design["account"],
        "-q",
        "single-gpu",
        "-N",
        f"msusie-match-{stage}",
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
    ]
    env = dict(
        FM_ROOT=str(root),
        FM_PYTHON=design["python"],
        FM_STAGE=stage,
        FM_WORKERS=str(design["workers"]),
    )
    if indices is not None:
        path = root / f"indices.{attempt}.json"
        if path.exists():
            raise RuntimeError(f"Preserve existing submission plan: {path}")
        write(path, indices)
        env["FM_INDICES"] = str(path)
        if len(indices) > 1:
            cmd += ["-J", f"0-{len(indices) - 1}%{design['concurrency']}"]
    if dependency:
        cmd += ["-W", f"depend=afterok:{dependency}"]
    cmd += ["-v", ",".join(f"{k}={v}" for k, v in env.items()), str(root / "job.pbs")]
    result = subprocess.run(cmd, capture_output=True, text=True)
    record = dict(
        command=cmd, stdout=result.stdout, stderr=result.stderr, returncode=result.returncode
    )
    if indices is not None:
        record["indices_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    write(root / f"submission.{stage}.{attempt}.json", record)
    result.check_returncode()
    return result.stdout.strip()


def scheduler(job):
    result = subprocess.run(
        ["qstat", "-xf", "-F", "json", job], capture_output=True, text=True, timeout=30
    )
    result.check_returncode()
    jobs = json.loads(result.stdout)["Jobs"]
    if len(jobs) != 1:
        raise RuntimeError("Ambiguous PBS response")
    return next(iter(jobs.values()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    lock = (root / "monitor.lock").open("w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    state = read(root / "monitor.status.json")
    tasks = read(root / "tasks.json")
    errors = 0
    try:
        while True:
            state["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
            if (root / "accepted.json").exists():
                accepted = read(root / "accepted.json")
                for path, digest in accepted["sha256"].items():
                    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
                        raise RuntimeError("Final result checksum mismatch")
                if not accepted["pass_"]:
                    raise RuntimeError("Final comparison not accepted")
                state["state"] = "complete"
                write(root / "monitor.status.json", state)
                return
            try:
                jobs = {k: scheduler(v) for k, v in state["jobs"].items()}
                errors = 0
            except (subprocess.SubprocessError, ValueError, KeyError) as exc:
                errors += 1
                state["scheduler_error"] = repr(exc)
                write(root / "monitor.status.json", state)
                if errors >= 5:
                    raise
                time.sleep(60)
                continue
            state["scheduler"] = {
                k: dict(state=v["job_state"], exit_status=v.get("Exit_status"))
                for k, v in jobs.items()
            }
            smoke = (
                root / tasks[0]["phase"] / tasks[0]["name"] / tasks[0]["locus_id"] / "smoke.json"
            )
            if not smoke.exists():
                if jobs["smoke"]["job_state"] == "F":
                    raise RuntimeError("Preflight did not pass; inspect smoke log")
                write(root / "monitor.status.json", state)
                time.sleep(60)
                continue
            pending = [
                i
                for i, t in enumerate(tasks)
                if not (root / t["phase"] / t["name"] / t["locus_id"] / "accepted.json").exists()
            ]
            state["remaining_tasks"] = len(pending)
            state["state"] = "running"
            if jobs["run"]["job_state"] == "F" and pending:
                logid = state["jobs"]["run"].split("[")[0]
                for log in (root / "logs").glob(f"run.{logid}*"):
                    if "Nonconverged SuSiEx at 10000" in log.read_text():
                        raise RuntimeError(f"Persistent nonconvergence: {log}")
                if state["round"] >= 5:
                    raise RuntimeError("Incomplete after five checkpoint-preserving retries")
                state["round"] += 1
                state["jobs"]["run"] = submit(root, "run", pending, attempt=state["round"])
            if not pending and jobs["collect"]["job_state"] == "F":
                logfile = root / "logs" / f"collect.{state['jobs']['collect']}.log"
                if logfile.exists():
                    text = logfile.read_text()
                    if text.strip() and not (
                        "FileNotFoundError" in text and "/accepted.json" in text
                    ):
                        raise RuntimeError(f"Collector failed: {logfile}")
                state["collect_round"] += 1
                if state["collect_round"] > 3:
                    raise RuntimeError("Collector repeatedly failed")
                state["jobs"]["collect"] = submit(root, "collect", attempt=state["collect_round"])
            write(root / "monitor.status.json", state)
            time.sleep(60)
    except Exception as exc:
        state.update(
            state="needs_inspection",
            error=repr(exc),
            updated_at_utc=datetime.now(timezone.utc).isoformat(),
        )
        write(root / "monitor.status.json", state)
        raise


if __name__ == "__main__":
    main()
