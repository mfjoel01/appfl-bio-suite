"""APPFL coordinator: model updates and calibration aggregates cross the boundary."""

import hashlib
import importlib.metadata
import json
import math
import os
import platform
import subprocess
from pathlib import Path

import numpy as np
import torch

from ..pipeline import write_json
from . import ATTRIBUTION
from .calibration import select_threshold
from .config import Config
from .privacy import accounting, epsilon
from .secure import pool, provision, wire_bytes
from .worker import Worker, read_tensors, save_tensors


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


class Runtime:
    def __init__(self, settings, out, command=None):
        self.process = None
        self.worker = None
        if command:
            settings_path = out / "worker-settings.json"
            write_json(settings_path, settings)
            self.log = (out / "worker.log").open("w")
            self.process = subprocess.Popen(
                [
                    *command,
                    "-m",
                    "appfl_bio_suite.experiments.methylation.alma3.worker",
                    "--settings",
                    str(settings_path),
                ],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=self.log,
                env={**os.environ, "NUMPY_MADVISE_HUGEPAGE": "0"},
                text=True,
                bufsize=1,
            )
        elif settings["fixture"]:
            self.worker = Worker(settings)
        else:
            raise ValueError("real ALMA3 requires an explicit isolated --worker-command")

    def call(self, **request):
        if self.worker is not None:
            return self.worker.handle(request)
        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("ALMA3 worker exited; see worker.log")
        response = json.loads(line)
        if not response["ok"]:
            raise RuntimeError(response["error"])
        return response["result"]

    def close(self):
        if self.process:
            self.process.stdin.close()
            try:
                self.process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                self.process.terminate()
                self.process.wait(timeout=10)
            self.log.close()


def aggregate(updates, sizes, *, secure, context, reference=None):
    """Plain sample-weighted APPFL FedAvg, or APPFL applied to a masked pooled mean.

    Secure clients send preweighted parameters. The coordinator never reconstructs
    an individual update, and APPFL receives one virtual client with total weight.
    """
    from appfl.algorithm.aggregator import FedAvgAggregator
    from omegaconf import OmegaConf

    if secure:
        summed = pool(updates, context)
        updates = {
            "secure-pooled": {
                k: v.to(dtype=reference[k].dtype if reference is not None else torch.float32)
                for k, v in summed.items()
            }
        }
        sizes = {"secure-pooled": sum(sizes.values())}
    aggregator = FedAvgAggregator(
        model=None, aggregator_configs=OmegaConf.create({"client_weights_mode": "sample_size"})
    )
    for cid, size in sizes.items():
        aggregator.set_client_sample_size(cid, size)
    return aggregator.aggregate(updates)


def request_base(seed, depth, scratch, state):
    return {"seed": seed, "depth": depth, "scratch": scratch, "state": str(state)}


def calibrate(runtime, base, config, sites, sigma, levels, folder):
    thresholds = {}
    communication = 0
    for coverage in config.coverage:
        prior = []
        for level in range(levels):
            roster = [f"S{s}" for s in sites]
            secure = config.secure and len(roster) > 1
            protocol, keys = (
                provision(roster, f"calibration:{coverage}:{level}") if secure else (None, {})
            )
            updates = {}
            for site, cid in zip(sites, roster, strict=True):
                output = folder / f"calibration-{cid}.npz"
                runtime.call(
                    **base,
                    op="calibrate",
                    split="calibration",
                    site=site,
                    coverage=coverage,
                    level=level,
                    prior=prior,
                    sigma=sigma,
                    output=str(output),
                    cid=cid,
                    secure=protocol,
                    keys=keys.get(cid),
                )
                updates[cid] = read_tensors(output)
                communication += wire_bytes(updates[cid])
                output.unlink()
            counts = (
                pool(updates, protocol)["histogram"]
                if secure
                else sum(v["histogram"] for v in updates.values())
            )
            threshold = select_threshold(
                counts.numpy(),
                floor=config.precision_floor,
                minimum=config.min_calibration,
                noise_std=math.sqrt(2 * len(sites)) * sigma,
            )
            prior.append(threshold)
        thresholds[str(coverage)] = prior
    return thresholds, communication


def evaluate(runtime, base, config, thresholds, split):
    rows = []
    for coverage in config.coverage:
        for site in range(len(config.site_sizes)):
            for unshifted in (False, True):
                result = runtime.call(
                    **base,
                    op="evaluate",
                    split=split,
                    site=site,
                    coverage=coverage,
                    unshifted=unshifted,
                    thresholds=thresholds[str(coverage)],
                )
                rows.append(
                    {
                        "coverage": coverage,
                        "site": site,
                        "unshifted": unshifted,
                        "split": split,
                        **result,
                    }
                )
    return rows


def guard_output(out, release):
    out = Path(out).resolve()
    if release and (out == Path(release).resolve() or out.is_relative_to(Path(release).resolve())):
        raise ValueError("outputs must be outside the immutable ALMA3 release")
    if out.exists():
        raise FileExistsError(f"output already exists: {out}")
    return out


def run(out, config, *, fixture=False, release=None, member=None, command=None, device="cpu"):
    config.validate()
    if not fixture and (not release or not member or not command):
        raise ValueError("real run requires release, neural member and isolated worker command")
    out = guard_output(out, release)
    out.mkdir(parents=True)
    settings = {
        "config": config.json(),
        "fixture": fixture,
        "release": str(release) if release else None,
        "member": member,
        "device": device,
        "allow_test": False,
    }
    runtime = Runtime(settings, out, command)
    records = []
    try:
        identity = runtime.call(op="identity")
        identity["smallest_site"] = int(np.argmin(config.site_sizes))
        license_text = identity.pop("license", None)
        if license_text:
            (out / "LICENSE.ALMA3").write_text(license_text)
        write_json(out / "config.json", config.json())
        write_json(out / "identity.json", identity)
        if not fixture:
            (out / "ATTRIBUTION.txt").write_text(
                ATTRIBUTION + "\nCheckpoints are based on ALMA3 (Version 3.1.2).\n"
                "ALMA Genomics Inc. (2026). https://github.com/almagx/alma3\n"
            )
        roster = [f"S{i}" for i in range(len(config.site_sizes))]
        sizes = dict(zip(roster, config.site_sizes, strict=True))
        for seed in config.seeds:
            for depth in config.depths:
                initials = {}
                for scratch in (False, True):
                    if scratch and config.arms and "scratch-fedavg" not in config.arms:
                        continue
                    path = out / f"initial-{seed}-{depth}-{scratch}.npz"
                    runtime.call(
                        op="initial",
                        seed=seed,
                        depth="full" if scratch else depth,
                        scratch=scratch,
                        output=str(path),
                    )
                    initials[scratch] = path
                for sigma in config.noise_multipliers:
                    arms = ["zero-shot", "pooled", "fedavg", "fedprox", "scratch-fedavg"] + [
                        f"local-S{i}" for i in range(len(roster))
                    ]
                    if not fixture:
                        arms.append("template-zero-shot")
                        if config.release_ensemble:
                            arms.append("ensemble-zero-shot")
                    if config.arms:
                        if set(config.arms) - set(arms):
                            raise ValueError("requested arm unavailable for this runtime")
                        arms = list(config.arms)
                    for arm in arms:
                        scratch = arm == "scratch-fedavg"
                        actual_depth = "full" if scratch else depth
                        folder = out / f"seed-{seed}" / depth / f"noise-{sigma}" / arm
                        folder.mkdir(parents=True)
                        if license_text:
                            (folder / "LICENSE.ALMA3").write_text(license_text)
                            (folder / "ATTRIBUTION.txt").write_text(
                                (out / "ATTRIBUTION.txt").read_text()
                            )
                        state = initials[scratch]
                        communication = 0
                        compute_seconds = 0.0
                        local_site = int(arm.split("S")[1]) if arm.startswith("local-") else None
                        if not arm.endswith("zero-shot"):
                            for round_index in range(config.rounds):
                                base = request_base(seed, actual_depth, scratch, state)
                                next_state = folder / f"round-{round_index}.npz"
                                if arm in {"fedavg", "fedprox", "scratch-fedavg"}:
                                    protocol, keys = (
                                        provision(roster, f"{folder}:{round_index}")
                                        if (config.secure)
                                        else (None, {})
                                    )
                                    updates = {}
                                    # Count trainable downlink tensors for each site.
                                    reference = read_tensors(state)
                                    communication += wire_bytes(reference) * len(roster)
                                    for site, cid in enumerate(roster):
                                        update_file = folder / f"update-{cid}.npz"
                                        result = runtime.call(
                                            **base,
                                            op="train",
                                            sites=[site],
                                            sigma=sigma,
                                            mu=config.proximal_mu if arm == "fedprox" else 0.0,
                                            weight=sizes[cid] / sum(sizes.values())
                                            if config.secure
                                            else 1.0,
                                            output=str(update_file),
                                            cid=cid,
                                            secure=protocol,
                                            keys=keys.get(cid),
                                        )
                                        compute_seconds += result["seconds"]
                                        updates[cid] = read_tensors(update_file)
                                        communication += wire_bytes(updates[cid])
                                        update_file.unlink()
                                    values = aggregate(
                                        updates,
                                        sizes,
                                        secure=config.secure,
                                        context=protocol,
                                        reference=reference,
                                    )
                                    save_tensors(next_state, values)
                                else:
                                    sites = (
                                        list(range(len(roster)))
                                        if local_site is None
                                        else [local_site]
                                    )
                                    result = runtime.call(
                                        **base,
                                        op="train",
                                        sites=sites,
                                        sigma=sigma,
                                        output=str(next_state),
                                    )
                                    compute_seconds += result["seconds"]
                                state = next_state
                        base = request_base(seed, actual_depth, scratch, state)
                        baseline = {
                            "template-zero-shot": "template",
                            "ensemble-zero-shot": "ensemble",
                        }.get(arm)
                        if baseline:
                            base["baseline"] = baseline
                        sites = list(range(len(roster))) if local_site is None else [local_site]
                        thresholds, calibration_bytes = calibrate(
                            runtime, base, config, sites, sigma, identity["levels"], folder
                        )
                        policy = {
                            "kind": "research-precision-floor-v1",
                            "thresholds": thresholds,
                            "floor": config.precision_floor,
                            "calibration_split": "calibration",
                            "original_release_policy": False,
                        }
                        write_json(folder / "calibration.json", policy)
                        record = {
                            "seed": seed,
                            "depth": depth,
                            "sigma": sigma,
                            "arm": arm,
                            "scratch": scratch,
                            "actual_depth": actual_depth,
                            "baseline": baseline,
                            "state": str(state.relative_to(out)),
                            "calibration": str((folder / "calibration.json").relative_to(out)),
                            "trainable_parameters": sum(
                                v.numel() for v in read_tensors(state).values()
                            ),
                            "training_wire_bytes": communication,
                            "calibration_wire_bytes": calibration_bytes,
                            "training_seconds": compute_seconds,
                            "privacy": accounting(
                                sigma,
                                0
                                if arm.endswith("zero-shot")
                                else (config.rounds * config.local_steps),
                                len(config.coverage) * identity["levels"],
                                config.delta,
                            ),
                            "metrics": evaluate(runtime, base, config, thresholds, "validation"),
                        }
                        records.append(record)
                        write_json(out / "results.json", records)
                        print(
                            f"Completed seed={seed} depth={depth} noise={sigma} arm={arm}",
                            flush=True,
                        )
        private = [r["privacy"] for r in records if r["privacy"]["private"]]
        # Conservative sequential composition over every published private arm/seed/depth.
        rho = (
            max(sum(p["training_rho"] for p in private), sum(p["calibration_rho"] for p in private))
            if private
            else 0
        )
        try:
            commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
        except subprocess.CalledProcessError:
            commit = None
        write_json(
            out / "manifest.json",
            {
                "experiment": "methylation/alma3",
                "stage": "validation-only",
                "settings": settings,
                "suite_commit": commit,
                "appfl": importlib.metadata.version("appfl"),
                "coordinator_torch": torch.__version__,
                "coordinator_numpy": np.__version__,
                "python": platform.python_version(),
                "worker_command": command,
                "identity": identity,
                "private_subset_composed_epsilon": epsilon(rho, config.delta),
                "private_subset_composed_rho": rho,
                "delta": config.delta,
                "whole_study_private": False,
                "limitations": "Synthetic public evaluation; nonprivate reference arms. "
                "Single-process "
                "protocol simulation, no host isolation. Template-derived data bias. "
                "Single neural member in FP32, not certified ensemble inference. "
                "Scratch arms train all parameters; depth denotes the paired pretrained arm. "
                "Communication excludes release provisioning and transport envelopes.",
                "attribution": "Tiny stand-in; no ALMA3 weights" if fixture else ATTRIBUTION,
            },
        )
        from .report import report

        report(out, records, identity)
        return {"output": str(out), "models": len(records), "stage": "validation-only"}
    finally:
        runtime.close()


def freeze(directory):
    directory = Path(directory).resolve()
    manifest = json.loads((directory / "manifest.json").read_text())
    if manifest["stage"] != "validation-only":
        raise ValueError("only completed validation studies can be frozen")
    records = json.loads((directory / "results.json").read_text())
    files = {"config.json", "identity.json", "manifest.json", "results.json"}
    files.update(r[k] for r in records for k in ("state", "calibration"))
    files.update(p.name for p in directory.glob("LICENSE*"))
    if (directory / "ATTRIBUTION.txt").exists():
        files.add("ATTRIBUTION.txt")
    hashes = {name: digest(directory / name) for name in sorted(files)}
    # Exclusive create: no silently replacing an already frozen specification.
    with (directory / "frozen.json").open("x") as stream:
        json.dump({"kind": "alma3-frozen-study-v1", "hashes": hashes}, stream, indent=2)
    return directory / "frozen.json"


def test_frozen(frozen, out, *, command=None):
    frozen = Path(frozen).resolve()
    directory = frozen.parent
    contract = json.loads(frozen.read_text())
    if contract.get("kind") != "alma3-frozen-study-v1":
        raise ValueError("invalid frozen study")
    for name, expected in contract["hashes"].items():
        path = (directory / name).resolve()
        if not path.is_relative_to(directory) or digest(path) != expected:
            raise ValueError(f"frozen artifact changed: {name}")
    manifest = json.loads((directory / "manifest.json").read_text())
    settings = {**manifest["settings"], "allow_test": True}
    config = Config(**settings["config"])
    out = guard_output(out, settings["release"])
    # Reserve before any test generation. An interrupted test remains marked as inspected.
    fd = os.open(directory / "TEST_STARTED", os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    with os.fdopen(fd, "w") as stream:
        stream.write(digest(frozen) + "\n")
    out.mkdir(parents=True)
    runtime = Runtime(settings, out, command or manifest["worker_command"])
    try:
        observed_identity = runtime.call(op="identity")
        observed_identity.pop("license", None)
        expected_identity = dict(manifest["identity"])
        expected_identity.pop("smallest_site", None)
        if observed_identity != expected_identity:
            raise ValueError("test runtime or release differs from the frozen study")
        records = json.loads((directory / "results.json").read_text())
        for record in records:
            policy = json.loads((directory / record["calibration"]).read_text())
            base = request_base(
                record["seed"],
                record["actual_depth"],
                record["scratch"],
                directory / record["state"],
            )
            if record.get("baseline"):
                base["baseline"] = record["baseline"]
            record["metrics"] = evaluate(runtime, base, config, policy["thresholds"], "test")
        write_json(out / "results.json", records)
        write_json(
            out / "manifest.json",
            {**manifest, "stage": "test-once", "frozen_sha256": digest(frozen)},
        )
        from .report import report

        report(out, records, manifest["identity"])
        return {"output": str(out), "stage": "test-once"}
    finally:
        runtime.close()
