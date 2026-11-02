"""Single-process E1–E4 experiment using APPFL's sample-weighted FedAvg aggregator."""

from __future__ import annotations

import copy
import importlib.metadata
import json
import numbers
import platform
import subprocess
from dataclasses import asdict, dataclass, fields
from pathlib import Path

import numpy as np
import torch
import yaml

from appfl_bio_suite.core.simulation import file_checksum

from .dataset import MethylationData, finite_number, synthetic_data, validate_synthetic_options
from .evaluation import fit_temperature, metrics
from .model import MethylationMLP
from .preprocessing import prepare, simulate_calls


@dataclass
class RunConfig:
    seed: int = 42
    samples: int = 240
    cpgs: int = 512
    classes: int = 3
    signal_fraction: float = 1.0
    class_weights: tuple[float, ...] | None = None
    profile_noise: float = 0.0
    missing_fraction: float = 0.0
    min_site_class_patients: int = 0
    panel_size: int = 512
    min_class: int = 20
    alpha: float = 0.5
    small_fraction: float = 0.1
    batch_shift: float = 0.03
    rounds: int = 3
    local_epochs: int = 1
    batch_size: int = 32
    hidden: int = 64
    bottleneck: int = 32
    dropout: float = 0.5
    learning_rate: float = 0.001
    weight_decay: float = 0.0001
    min_coverage: float = 0.015
    caller_intercept: float = -0.155
    caller_slope: float = 1.517
    coverage_counts: tuple[int, ...] = (16, 64, 128, 256, 512)
    device: str = "cpu"
    threads: int = 1
    # Score every model on validation patients after each round. Never touches test.
    track_validation: bool = False
    # A development config can declare itself validation-only, as --validation-only does.
    validation_only: bool = False

    def validate(self):
        if not _integer(self.seed) or not 0 <= self.seed < 2**32 - 100000:
            raise ValueError("seed must be a nonnegative integer below 2**32 - 100000")
        for field in (
            "alpha",
            "small_fraction",
            "batch_shift",
            "dropout",
            "learning_rate",
            "weight_decay",
            "min_coverage",
            "caller_intercept",
            "caller_slope",
            "signal_fraction",
        ):
            if not finite_number(getattr(self, field)):
                raise ValueError(f"{field} must be a finite number")
        for field in (
            "samples",
            "cpgs",
            "classes",
            "panel_size",
            "min_class",
            "rounds",
            "local_epochs",
            "batch_size",
            "hidden",
            "bottleneck",
            "threads",
        ):
            value = getattr(self, field)
            if not _integer(value) or value < 1:
                raise ValueError(f"{field} must be a positive integer")
        if not _integer(self.min_site_class_patients) or self.min_site_class_patients < 0:
            raise ValueError("min_site_class_patients must be a nonnegative integer")
        validate_synthetic_options(
            self.classes, self.class_weights, self.profile_noise, self.missing_fraction
        )
        if not 0 < self.signal_fraction <= 1:
            raise ValueError("signal_fraction must be in (0, 1]")
        if self.classes < 2 or not 0 <= self.dropout < 1:
            raise ValueError("require at least two classes and dropout in [0, 1)")
        if not 0 < self.min_coverage <= 1 or not 0 < self.small_fraction < 0.5:
            raise ValueError("invalid coverage or small-site fraction")
        if self.alpha <= 0 or self.learning_rate <= 0 or self.weight_decay < 0:
            raise ValueError("invalid alpha or optimizer configuration")
        if not 0 <= self.batch_shift <= 1:
            raise ValueError("batch_shift must lie in [0, 1]")
        if (
            not isinstance(self.coverage_counts, list | tuple)
            or not self.coverage_counts
            or any(not _integer(k) or k < 1 for k in self.coverage_counts)
        ):
            raise ValueError("coverage_counts must be a list of positive integers")
        if self.device != "cpu" and self.device != "cuda":
            raise ValueError("device must be cpu or cuda")
        for field in ("track_validation", "validation_only"):
            if not isinstance(getattr(self, field), bool):
                raise ValueError(f"{field} must be true or false")


def _integer(value) -> bool:
    return isinstance(value, numbers.Integral) and not isinstance(value, bool)


def load_config(name: str = "ci-tiny") -> RunConfig:
    """Load a packaged scenario by name, or a YAML file of RunConfig fields by path.

    Every failure is a ValueError naming the problem, which the CLI reports cleanly.
    """
    path = Path(name)
    if not path.is_file():
        packaged = Path(__file__).parent / "configs"
        path = packaged / f"{name}.yaml"
        if not path.is_file():
            names = ", ".join(sorted(p.stem for p in packaged.glob("*.yaml")))
            raise ValueError(f"no config file or packaged scenario {name!r}; scenarios: {names}")
    values = yaml.safe_load(path.read_text()) or {}
    if not isinstance(values, dict):
        raise ValueError(f"{path} must be a mapping of RunConfig fields")
    unknown = sorted(set(values) - {field.name for field in fields(RunConfig)})
    if unknown:
        raise ValueError(f"unknown setting(s) in {path}: {', '.join(unknown)}")
    config = RunConfig(**values)
    config.validate()
    return config


def write_json(path: Path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _state(model):
    return {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}


def train_epoch(model, X, y, config, *, weights, seed):
    """Reset AdamW each communication round for all comparison arms."""
    model.train()
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay
    )
    loss_fn = torch.nn.CrossEntropyLoss(weight=torch.tensor(weights, device=config.device))
    total_loss, total = 0.0, 0
    for _ in range(config.local_epochs):
        order = rng.permutation(len(y))
        for start in range(0, len(order), config.batch_size):
            idx = order[start : start + config.batch_size]
            calls = simulate_calls(
                X[idx],
                rng,
                min_fraction=config.min_coverage,
                intercept=config.caller_intercept,
                slope=config.caller_slope,
            )
            inputs = torch.tensor(calls, device=config.device)
            targets = torch.tensor(y[idx], device=config.device, dtype=torch.long)
            optimizer.zero_grad()
            loss = loss_fn(model(inputs), targets)
            if not torch.isfinite(loss):
                raise ValueError("training loss became nonfinite")
            loss.backward()
            optimizer.step()
            total_loss += float(loss.detach()) * len(idx)
            total += len(idx)
    return total_loss / total


def predict(model, calls, config):
    model.eval()
    outputs = []
    with torch.no_grad():
        for start in range(0, len(calls), config.batch_size):
            inputs = torch.tensor(calls[start : start + config.batch_size], device=config.device)
            outputs.append(model(inputs).cpu().numpy())
    return np.concatenate(outputs)


def run_experiment(
    out: Path, config: RunConfig, input_path: Path | None = None, *, validation_only: bool = False
) -> dict:
    config.validate()
    validation_only = validation_only or config.validation_only
    is_synthetic = input_path is None
    if input_path and (input_path.parent / "manifest.json").is_file():
        source_manifest = json.loads((input_path.parent / "manifest.json").read_text())
        is_synthetic = (
            source_manifest.get("experiment") == "methylation"
            and source_manifest.get("stage") == "synthetic-cohort"
            and source_manifest.get("dataset_sha256") == file_checksum(input_path)
        )
    if out.exists():
        raise FileExistsError(f"output directory already exists: {out}")
    data = (
        MethylationData.load(input_path)
        if input_path
        else synthetic_data(
            samples=config.samples,
            cpgs=config.cpgs,
            classes=config.classes,
            seed=config.seed,
            signal_fraction=config.signal_fraction,
            class_weights=config.class_weights,
            profile_noise=config.profile_noise,
            missing_fraction=config.missing_fraction,
        )
    )
    prepared = prepare(
        data,
        panel_size=config.panel_size,
        seed=config.seed,
        min_class=config.min_class,
        alpha=config.alpha,
        small_fraction=config.small_fraction,
        min_site_class_patients=config.min_site_class_patients,
    )
    if max(config.coverage_counts) > len(prepared.panel):
        raise ValueError("coverage count exceeds selected panel size; adjust coverage_counts")
    if config.device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA requested but unavailable")
    # A true APPFL component, without network/remote-agent semantics in this serial PoC.
    from appfl.algorithm.aggregator import FedAvgAggregator
    from omegaconf import OmegaConf

    out.mkdir(parents=True)
    write_json(out / "status.json", {"stage": "prepared", "complete": False})
    write_json(out / "config.json", asdict(config))
    data.save(out / "cohort.npz")
    record = {
        "classes": prepared.classes.tolist(),
        "label_mapping": prepared.label_mapping,
        "cpg_ids": data.cpg_ids[prepared.panel].tolist(),
        "panel_indices": prepared.panel.tolist(),
        "splits": {key: getattr(prepared, key).tolist() for key in ("train", "validation", "test")},
        "patient_ids": data.patient_ids.tolist(),
        "sample_ids": data.sample_ids.tolist(),
        "sites": {f"S{i + 1}": idx.tolist() for i, idx in enumerate(prepared.sites)},
        "site_class_counts": {
            f"S{i + 1}": dict(
                zip(
                    prepared.classes.tolist(),
                    np.bincount(prepared.y[idx], minlength=len(prepared.classes)).tolist(),
                    strict=True,
                )
            )
            for i, idx in enumerate(prepared.sites)
        },
    }
    # Batch shifts are a synthetic training-site perturbation. Held-out data is unchanged.
    shifts = np.random.default_rng(config.seed + 1).normal(0, config.batch_shift, size=4)
    if not is_synthetic:
        shifts[:] = 0
    for idx, shift in zip(prepared.sites, shifts, strict=True):
        prepared.X[idx] = np.clip(prepared.X[idx] + shift, 0, 1)
    record["synthetic_site_shifts"] = shifts.tolist()
    write_json(out / "preprocessing.json", record)
    torch.set_num_threads(config.threads)
    torch.use_deterministic_algorithms(True)
    torch.manual_seed(config.seed)
    initial = MethylationMLP(
        len(prepared.panel), len(prepared.classes), config.hidden, config.bottleneck, config.dropout
    ).to(config.device)
    global_model = copy.deepcopy(initial)
    aggregator = FedAvgAggregator(
        model=copy.deepcopy(global_model).cpu(),
        aggregator_configs=OmegaConf.create({"client_weights_mode": "sample_size"}),
    )
    for i, idx in enumerate(prepared.sites):
        aggregator.set_client_sample_size(f"S{i + 1}", len(idx))
    models = {"centralized": copy.deepcopy(initial), "federated": global_model}
    models.update({f"local_S{i + 1}": copy.deepcopy(initial) for i in range(4)})
    n_classes = len(prepared.classes)

    def weights(indices):
        counts = np.bincount(prepared.y[indices], minlength=n_classes)
        result = np.zeros(n_classes, dtype=np.float32)
        present = counts > 0
        result[present] = len(indices) / (present.sum() * counts[present])
        return result

    pooled_weights = weights(prepared.train)

    def scoring_calls(indices, k):
        """Fixed calls at k observed CpGs, identical for every arm and every round."""
        return simulate_calls(
            prepared.X[indices],
            np.random.default_rng(config.seed + 3000 + k),
            observed=k,
            intercept=config.caller_intercept,
            slope=config.caller_slope,
        )

    # Per-round scoring reuses a validation-only run's final masks, so its last round
    # reproduces that evaluation exactly. Scoring draws no training randomness.
    tracked = {}
    if config.track_validation:
        tracked = {k: scoring_calls(prepared.validation, k) for k in config.coverage_counts}
    history, validation_history = [], []
    for round_id in range(config.rounds):
        losses = {}
        for name, idx in [("centralized", prepared.train)] + [
            (f"local_S{i + 1}", idx) for i, idx in enumerate(prepared.sites)
        ]:
            # Local baselines use only their own class counts.
            losses[name] = train_epoch(
                models[name],
                prepared.X[idx],
                prepared.y[idx],
                config,
                weights=weights(idx),
                seed=config.seed + 1000 + round_id * 10,
            )
        updates = {}
        for i, idx in enumerate(prepared.sites):
            client = copy.deepcopy(global_model)
            losses[f"client_S{i + 1}"] = train_epoch(
                client,
                prepared.X[idx],
                prepared.y[idx],
                config,
                weights=pooled_weights,
                seed=config.seed + 1000 + round_id * 10 + i,
            )
            updates[f"S{i + 1}"] = _state(client)
            del client
        global_model.load_state_dict(aggregator.aggregate(updates))
        history.append({"round": round_id + 1, "losses": losses})
        torch.save(
            {name: _state(model) for name, model in models.items()}, out / "checkpoint.pt.tmp"
        )
        (out / "checkpoint.pt.tmp").replace(out / "checkpoint.pt")
        write_json(out / "history.json", history)
        if tracked:
            scores = []
            for k, calls in tracked.items():
                for name, model in models.items():
                    result = metrics(predict(model, calls, config), prepared.y[prepared.validation])
                    scores.append(
                        {
                            "model": name,
                            "requested_cpgs": k,
                            "balanced_accuracy": result["balanced_accuracy"],
                            "confusion_matrix": result["confusion_matrix"],
                        }
                    )
            validation_history.append({"round": round_id + 1, "scores": scores})
            write_json(out / "validation_history.json", validation_history)
        write_json(
            out / "status.json", {"stage": "training", "round": round_id + 1, "complete": False}
        )
        print(f"Round {round_id + 1}/{config.rounds} complete", flush=True)
    # Validation masks and test masks are fixed and shared across all arms.
    validation_calls = np.concatenate(
        [
            simulate_calls(
                prepared.X[prepared.validation],
                np.random.default_rng(config.seed + 2000 + k),
                observed=k,
                intercept=config.caller_intercept,
                slope=config.caller_slope,
            )
            for k in config.coverage_counts
        ]
    )
    validation_y = np.tile(prepared.y[prepared.validation], len(config.coverage_counts))
    temperatures = {
        name: (
            1.0
            if validation_only
            else fit_temperature(predict(model, validation_calls, config), validation_y)
        )
        for name, model in models.items()
    }
    results = []
    evaluation_indices = prepared.validation if validation_only else prepared.test
    evaluation_split = "validation" if validation_only else "test"
    for k in config.coverage_counts:
        calls = scoring_calls(evaluation_indices, k)
        for name, model in models.items():
            logits = predict(model, calls, config)
            results.append(
                {
                    "model": name,
                    "evaluation_split": evaluation_split,
                    "requested_cpgs": k,
                    "mean_observed_cpgs": float(np.count_nonzero(calls, axis=1).mean()),
                    "temperature": temperatures[name],
                    **metrics(
                        logits, prepared.y[evaluation_indices], temperature=temperatures[name]
                    ),
                }
            )
    if validation_only:
        for row in results:
            row["n_validation"] = row.pop("n_test")
    write_json(out / "metrics.json", results)
    from .reporting import report

    record["input_kind"] = "synthetic" if is_synthetic else "supplied-beta-cohort"
    if validation_only:
        (out / "REPORT.md").write_text(
            "# Validation-only development run\n\n"
            "Test patients were not evaluated. Temperatures are fixed at 1.\n"
        )
    else:
        report(out, results, record)
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path(__file__).parent,
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=Path(__file__).parent,
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    write_json(
        out / "manifest.json",
        {
            "schema_version": 1,
            "experiment": "methylation",
            "suite_commit": commit,
            "working_tree_dirty": bool(dirty),
            "python": platform.python_version(),
            "packages": {
                name: importlib.metadata.version(name)
                for name in ("appfl", "torch", "numpy", "scikit-learn", "matplotlib")
            },
            "aggregator_source": __import__("inspect").getfile(FedAvgAggregator),
            "aggregator_sha256": file_checksum(__import__("inspect").getfile(FedAvgAggregator)),
            "source_sha256": {
                p.name: file_checksum(p) for p in sorted(Path(__file__).parent.glob("*.py"))
            },
            "input_kind": "synthetic" if is_synthetic else "supplied-beta-cohort",
            "evaluation_split": evaluation_split,
            "input_sha256": file_checksum(input_path) if input_path else None,
            "outputs": {
                p.name: file_checksum(p)
                for p in sorted(out.iterdir())
                if p.is_file() and p.name != "status.json"
            },
            "limitations": [
                "Single-node simulation; no clinical validation",
                "Uniform probe sampling; no lognormal read-depth simulation",
                "Verified caller coefficients; Bernoulli calls omit empirical read-depth fitting",
                "Array-level splits are grouped by patient; metrics are per sample",
            ],
        },
    )
    write_json(out / "status.json", {"stage": "complete", "complete": True})
    return {"output": str(out), "models": len(models), "evaluations": len(results)}
