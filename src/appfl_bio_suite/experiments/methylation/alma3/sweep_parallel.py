"""Run up to three independent CPU sweep seeds; serialize all progress commits."""

import argparse
import importlib.metadata
import json
import os
import platform
import shlex
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

from .sweep import atomic_json, configuration
from .sweep_report import report_sweep


def commit_snapshot(destination, points):
    handoff = Path("docs/experiments/methylation/ALMA3-HANDOFF.md")
    with handoff.open("a") as stream:
        stream.write(
            f"\nParallel sweep checkpoint: {len(points)} completed seed/level points. "
            f"Durable scores: `{destination}`.\n"
        )
    subprocess.run(["git", "add", str(destination), str(handoff)], check=True)
    staged = subprocess.check_output(
        ["git", "diff", "--cached", "--name-only"], text=True
    ).splitlines()
    if any(
        not (Path(p).resolve().is_relative_to(destination.resolve()) or Path(p) == handoff)
        for p in staged
    ):
        raise RuntimeError("unrelated staged files: scores saved but commit stopped")
    if staged:
        subprocess.run(
            ["git", "commit", "-m", f"Checkpoint ALMA3 sweep: {len(points)} seed-level points"],
            check=True,
        )


def run_parallel(spec, out, release, worker_command, publish, *, resume=False):
    configuration(spec)
    if not 1 <= len(spec["seeds"]) <= 3:
        raise ValueError("this runner is capped at three CPU seed workers")
    out, release = Path(out).resolve(), Path(release).resolve()
    publish = Path(publish)
    if out.is_relative_to(release) or publish.resolve().is_relative_to(release):
        raise ValueError("outputs must be outside immutable release")
    if out.exists():
        if not resume or json.loads((out / "specification.json").read_text()) != spec:
            raise ValueError("output exists or resume specification differs")
    else:
        out.mkdir(parents=True)
        atomic_json(out / "specification.json", spec)
    publish.mkdir(parents=True, exist_ok=True)
    atomic_json(publish / "specification.json", spec)
    children, logs = [], []
    try:
        for seed in spec["seeds"]:
            child_spec = {**spec, "seeds": [seed]}
            config_path = out / f"specification-seed-{seed}.json"
            atomic_json(config_path, child_spec)
            child_out = out / f"worker-seed-{seed}"
            child_publish = publish / f"seed-{seed}"
            command = [
                sys.executable,
                "-m",
                "appfl_bio_suite.experiments.methylation.alma3.sweep",
                "--config",
                str(config_path),
                "--out",
                str(child_out),
                "--release",
                str(release),
                "--worker-command",
                worker_command,
                "--publish",
                str(child_publish),
            ]
            if child_out.exists():
                command.append("--resume")
            log = (out / f"seed-{seed}.log").open("a" if resume else "w")
            logs.append(log)
            child = subprocess.Popen(
                command,
                stdout=log,
                stderr=subprocess.STDOUT,
                env={**os.environ, "NUMPY_MADVISE_HUGEPAGE": "0", "OPENBLAS_NUM_THREADS": "1"},
            )
            children.append((seed, child))
        environment = dict(
            parent_pid=os.getpid(),
            host=platform.node(),
            python=platform.python_version(),
            seed_pids={str(s): p.pid for s, p in children},
            implementation_commit=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            packages={
                name: importlib.metadata.version(name)
                for name in ("appfl", "torch", "numpy", "scipy", "matplotlib")
            },
            worker_command=shlex.split(worker_command),
            release=str(release),
        )
        atomic_json(out / ("resume-execution.json" if resume else "execution.json"), environment)
        previous = None
        while True:
            points = sorted(
                str(p.relative_to(publish))
                for p in publish.glob("seed-*/seed-*-convergence-*.json")
            )
            if points != previous:
                status = dict(
                    completed_points=points,
                    expected_points=len(spec["seeds"]) * len(spec["levels"]),
                    output=str(out),
                    stage="validation-only; no test access",
                )
                atomic_json(out / "progress.json", status)
                atomic_json(publish / "status.json", status)
                if points:
                    commit_snapshot(publish, points)
                    print(
                        f"Committed {len(points)}/{status['expected_points']} seed-level points",
                        flush=True,
                    )
                previous = points
            if all(p.poll() is not None for _, p in children):
                final_points = sorted(
                    str(p.relative_to(publish))
                    for p in publish.glob("seed-*/seed-*-convergence-*.json")
                )
                if final_points != points:
                    continue
                break
            time.sleep(5)
        failures = {str(s): p.returncode for s, p in children if p.returncode != 0}
        if failures:
            atomic_json(out / "FAILED.json", failures)
            raise RuntimeError(f"seed workers failed; see per-seed logs: {failures}")
        records, contracts, identity = [], {}, None
        for seed, _ in children:
            child_out = out / f"worker-seed-{seed}"
            observed = json.loads((child_out / "identity.json").read_text())
            if identity is not None and observed != identity:
                raise ValueError("seed worker identities differ")
            identity = observed
            contracts[str(seed)] = json.loads((child_out / "contract.json").read_text())
            for record in json.loads((child_out / "results.json").read_text()):
                record["state"] = str(Path(f"worker-seed-{seed}") / record["state"])
                records.append(record)
        atomic_json(out / "results.json", records)
        atomic_json(out / "contracts.json", contracts)
        atomic_json(out / "identity.json", identity)
        for name in ("LICENSE.ALMA3", "ATTRIBUTION.txt"):
            shutil.copyfile(out / f"worker-seed-{spec['seeds'][0]}" / name, out / name)
        report_sweep(out)
        atomic_json(out / "COMPLETE.json", dict(records=len(records), stage="validation-only"))
        print(f"Completed {len(records)} evaluations; figures in {out}", flush=True)
    except BaseException:
        for _, process in children:
            if process.poll() is None:
                process.send_signal(signal.SIGINT)
        for _, process in children:
            if process.poll() is None:
                try:
                    process.wait(timeout=45)
                except subprocess.TimeoutExpired:
                    process.terminate()
        raise
    finally:
        for log in logs:
            log.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--release", required=True, type=Path)
    parser.add_argument("--worker-command", required=True)
    parser.add_argument("--publish", required=True, type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    run_parallel(
        json.loads(args.config.read_text()),
        args.out,
        args.release,
        args.worker_command,
        args.publish,
        resume=args.resume,
    )


if __name__ == "__main__":
    main()
