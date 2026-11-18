"""Resumable paired five-site label-ratio sweep; real-member validation only."""

import argparse
import hashlib
import itertools
import json
import os
import shlex
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path

import numpy as np

from .config import Config
from .secure import provision, wire_bytes
from .study import Runtime, aggregate, digest, evaluate, request_base
from .worker import read_tensors, save_tensors


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def interpolate_counts(extreme, target, convergence):
    """Nearest integer interpolation preserving row sizes AND pooled class totals."""
    extreme, target = np.asarray(extreme), np.asarray(target)
    if (
        extreme.ndim != 2
        or target.shape != (extreme.shape[1],)
        or extreme.dtype.kind not in "iu"
        or target.dtype.kind not in "iu"
        or (extreme < 1).any()
        or (target < 1).any()
        or not np.all(extreme.sum(1) == target.sum())
        or not np.array_equal(extreme.sum(0), len(extreme) * target)
        or not 0 <= convergence <= 1
    ):
        raise ValueError("invalid fixed-total sweep design")
    desired = (1 - convergence) * extreme + convergence * target
    floor = np.floor(desired + 1e-10).astype(int)
    row_need = target.sum() - floor.sum(1)
    col_need = tuple((extreme.sum(0) - floor.sum(0)).tolist())

    @lru_cache(None)
    def solve(row, remaining):
        if row == len(extreme):
            return (0.0, ()) if not any(remaining) else (-np.inf, ())
        best = (-np.inf, ())
        for columns in itertools.combinations(range(len(target)), int(row_need[row])):
            if any(remaining[c] <= 0 for c in columns):
                continue
            next_remaining = tuple(v - int(c in columns) for c, v in enumerate(remaining))
            score, rest = solve(row + 1, next_remaining)
            score += sum(desired[row, c] - floor[row, c] for c in columns)
            if score > best[0]:
                best = score, (columns, *rest)
        return best

    score, additions = solve(0, col_need)
    if not np.isfinite(score):
        raise ValueError("cannot round sweep design")
    for row, columns in enumerate(additions):
        floor[row, list(columns)] += 1
    return floor


def configuration(spec):
    classes = len(spec["class_indices"])
    sites = len(spec["extreme_counts"])
    interpolate_counts(spec["extreme_counts"], spec["target_counts"], 0)
    if sites != 5 or spec["patients_per_site"] != sum(spec["target_counts"]):
        raise ValueError("sweep requires five equal-sized sites")
    if len(set(spec["levels"])) != len(spec["levels"]) or not spec["levels"]:
        raise ValueError("convergence levels must be nonempty and unique")
    for level in spec["levels"]:
        interpolate_counts(spec["extreme_counts"], spec["target_counts"], level)
    config = Config(
        seeds=spec["seeds"],
        depths=("heads",),
        noise_multipliers=(0.0,),
        site_sizes=(spec["patients_per_site"],) * sites,
        validation_size=classes * spec["validation_per_class"],
        test_size=classes * spec["validation_per_class"],
        rounds=spec["rounds"],
        local_steps=spec["local_steps"],
        learning_rate=spec["learning_rate"],
        coverage=(spec["coverage"],),
        batch_shift=spec["batch_shift"],
        signature_rotation_degrees=spec["signature_rotation_degrees"],
        class_indices=spec["class_indices"],
        minimum_per_class=1,
        site_class_proportions=tuple((1 / classes,) * classes for _ in range(sites)),
    )
    config.validate()
    if type(spec["gradient_batch_size"]) is not int or spec["gradient_batch_size"] < 1:
        raise ValueError("invalid gradient batch size")
    return config


def compact_metrics(metrics, support, count):
    rows = []
    for m in metrics:
        observed = np.array(m["confusion_matrix"]).sum(1)
        expected = np.zeros(len(observed), dtype=int)
        expected[support] = count
        if not np.array_equal(observed, expected):
            raise ValueError("evaluation cohort is not the fixed balanced design")
        rows.append(
            {
                k: m[k]
                for k in (
                    "site",
                    "coverage",
                    "unshifted",
                    "split",
                    "accuracy",
                    "balanced_accuracy",
                    "support_restricted_balanced_accuracy",
                    "n_test",
                )
            }
        )
        rows[-1]["per_class_recall"] = {str(c): m["per_class_recall"][c] for c in support}
        rows[-1]["class_counts"] = observed[support].tolist()
    return rows


def publish_progress(out, destination, seed, level, records, commit):
    destination.mkdir(parents=True, exist_ok=True)
    scores = [
        {k: r[k] for k in ("seed", "convergence", "arm", "class_counts", "metrics")}
        for r in records
    ]
    atomic_json(destination / f"seed-{seed}-convergence-{level:g}.json", scores)
    shutil.copyfile(out / "specification.json", destination / "specification.json")
    completed = sorted(p.name for p in destination.glob("seed-*.json"))
    atomic_json(
        destination / "status.json",
        {
            "completed_points": completed,
            "output": str(out),
            "stage": "validation-only; no test access",
        },
    )
    if commit:
        handoff = Path("docs/experiments/methylation/ALMA3-HANDOFF.md")
        with handoff.open("a") as stream:
            stream.write(
                f"\nSweep checkpoint: seed {seed}, convergence {level:g} complete. "
                f"Scores in `{destination}`; resume command below remains valid.\n"
            )
        subprocess.run(["git", "add", str(destination), str(handoff)], check=True)
        staged = subprocess.check_output(
            ["git", "diff", "--cached", "--name-only"], text=True
        ).splitlines()
        if any(
            not (Path(p).resolve().is_relative_to(destination.resolve()) or Path(p) == handoff)
            for p in staged
        ):
            raise RuntimeError("unrelated staged files: checkpoint saved, git commit stopped")
        if staged:
            subprocess.run(
                ["git", "commit", "-m", f"Record ALMA3 sweep seed {seed} convergence {level:g}"],
                check=True,
            )


def run_sweep(spec, out, release, command, *, resume=False, publish=None, commit=False):
    config = configuration(spec)
    out, release = Path(out).resolve(), Path(release).resolve()
    if out.is_relative_to(release):
        raise ValueError("output must be outside release")
    code_files = [
        Path(__file__).with_name(name)
        for name in (
            "sweep.py",
            "cohort.py",
            "worker.py",
            "privacy.py",
            "model.py",
            "secure.py",
            "study.py",
            "calibration.py",
            "config.py",
        )
    ]
    contract = {
        "specification": spec,
        "release": str(release),
        "worker_command": command,
        "code_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in code_files},
        "evaluation": "fixed class-balanced validation, paired across convergence levels",
        "training": "class-wise prefixes of fixed independent patient reservoirs",
        "stage": "validation-only",
        "calibration": "not fitted; classification metrics only",
    }
    if out.exists():
        if not resume or json.loads((out / "contract.json").read_text()) != contract:
            raise ValueError("output exists or resume contract differs")
    else:
        out.mkdir(parents=True)
        atomic_json(out / "contract.json", contract)
        atomic_json(out / "specification.json", spec)
    settings = dict(
        config=config.json(),
        fixture=False,
        release=str(release),
        member=spec["member"],
        device="cpu",
        allow_test=False,
    )
    runtime = Runtime(settings, out, command)
    roster = [f"S{i}" for i in range(5)]
    sizes = dict(zip(roster, config.site_sizes, strict=True))
    all_records = []
    try:
        identity = runtime.call(op="identity")
        license_text = identity.pop("license")
        if (out / "identity.json").exists() and json.loads(
            (out / "identity.json").read_text()
        ) != identity:
            raise ValueError("runtime identity changed")
        atomic_json(out / "identity.json", identity)
        (out / "LICENSE.ALMA3").write_text(license_text)
        (out / "ATTRIBUTION.txt").write_text(
            "Analysis performed using ALMA3.\nBased on ALMA3 Version 3.1.2.\n"
        )
        for seed in spec["seeds"]:
            initial = out / f"initial-{seed}.npz"
            if not initial.exists():
                temporary = out / f".initial-{seed}.tmp.npz"
                runtime.call(
                    op="initial", seed=seed, depth="heads", scratch=False, output=str(temporary)
                )
                os.replace(temporary, initial)
            for level in spec["levels"]:
                counts = interpolate_counts(
                    spec["extreme_counts"], spec["target_counts"], level
                ).tolist()
                point_records = []
                for arm in ["zero-shot", "pooled", "fedavg", *[f"local-S{i}" for i in range(5)]]:
                    folder = out / f"seed-{seed}" / f"convergence-{level:g}" / arm
                    folder.mkdir(parents=True, exist_ok=True)
                    record_file = folder / "record.json"
                    if record_file.exists():
                        record = json.loads(record_file.read_text())
                        if digest(out / record["state"]) != record["state_sha256"]:
                            raise ValueError("completed checkpoint changed")
                        point_records.append(record)
                        all_records.append(record)
                        continue
                    shutil.copyfile(out / "LICENSE.ALMA3", folder / "LICENSE.ALMA3")
                    shutil.copyfile(out / "ATTRIBUTION.txt", folder / "ATTRIBUTION.txt")
                    progress_file = folder / "training-progress.json"
                    progress = (
                        json.loads(progress_file.read_text())
                        if progress_file.exists()
                        else {
                            "rounds": 0,
                            "state": str(initial.relative_to(out)),
                            "state_sha256": digest(initial),
                            "training_seconds": 0.0,
                            "wire_bytes": 0,
                        }
                    )
                    state = out / progress["state"]
                    if digest(state) != progress["state_sha256"]:
                        raise ValueError("partial checkpoint changed")
                    if arm != "zero-shot":
                        for round_index in range(progress["rounds"], config.rounds):
                            base = request_base(seed, "heads", False, state)
                            output = folder / f"round-{round_index}.npz"
                            if arm == "fedavg":
                                protocol, keys = provision(
                                    roster, f"sweep:{seed}:{level}:{round_index}"
                                )
                                reference = read_tensors(state)
                                progress["wire_bytes"] += wire_bytes(reference) * 5
                                updates = {}
                                for site, cid in enumerate(roster):
                                    update = folder / f"update-{cid}.npz"
                                    result = runtime.call(
                                        **base,
                                        op="train",
                                        sites=[site],
                                        class_counts=counts,
                                        gradient_batch_size=spec["gradient_batch_size"],
                                        sigma=0.0,
                                        weight=0.2,
                                        output=str(update),
                                        cid=cid,
                                        secure=protocol,
                                        keys=keys[cid],
                                    )
                                    progress["training_seconds"] += result["seconds"]
                                    updates[cid] = read_tensors(update)
                                    progress["wire_bytes"] += wire_bytes(updates[cid])
                                    update.unlink()
                                save_tensors(
                                    output,
                                    aggregate(
                                        updates,
                                        sizes,
                                        secure=True,
                                        context=protocol,
                                        reference=reference,
                                    ),
                                )
                            else:
                                sites = (
                                    list(range(5)) if arm == "pooled" else [int(arm.split("S")[1])]
                                )
                                result = runtime.call(
                                    **base,
                                    op="train",
                                    sites=sites,
                                    class_counts=counts,
                                    gradient_batch_size=spec["gradient_batch_size"],
                                    sigma=0.0,
                                    output=str(output),
                                )
                                progress["training_seconds"] += result["seconds"]
                            state = output
                            progress.update(
                                rounds=round_index + 1,
                                state=str(state.relative_to(out)),
                                state_sha256=digest(state),
                            )
                            atomic_json(progress_file, progress)
                    base = request_base(seed, "heads", False, state)
                    metrics = evaluate(
                        runtime,
                        base,
                        config,
                        {str(config.coverage[0]): [1.01] * identity["levels"]},
                        "validation",
                    )
                    record = dict(
                        seed=seed,
                        convergence=level,
                        arm=arm,
                        class_counts=counts,
                        state=str(state.relative_to(out)),
                        state_sha256=digest(state),
                        training_seconds=progress["training_seconds"],
                        wire_bytes=progress["wire_bytes"],
                        metrics=compact_metrics(
                            metrics, spec["class_indices"], spec["validation_per_class"]
                        ),
                    )
                    atomic_json(record_file, record)
                    point_records.append(record)
                    all_records.append(record)
                    print(f"Completed seed={seed} convergence={level:g} arm={arm}", flush=True)
                atomic_json(out / "results.json", all_records)
                if publish:
                    publish_progress(out, Path(publish), seed, level, point_records, commit)
        atomic_json(
            out / "COMPLETE.json",
            {"records": len(all_records), "contract_sha256": digest(out / "contract.json")},
        )
        from .sweep_report import report_sweep

        report_sweep(out)
    finally:
        runtime.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--release", required=True, type=Path)
    parser.add_argument("--worker-command", required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--publish", type=Path)
    parser.add_argument("--commit-progress", action="store_true")
    args = parser.parse_args()
    if args.commit_progress and not args.publish:
        parser.error("--commit-progress requires --publish")
    run_sweep(
        json.loads(args.config.read_text()),
        args.out,
        args.release,
        shlex.split(args.worker_command),
        resume=args.resume,
        publish=args.publish,
        commit=args.commit_progress,
    )


if __name__ == "__main__":
    main()
