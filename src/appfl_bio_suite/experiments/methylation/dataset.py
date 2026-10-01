"""Portable methylation input contract and small synthetic development fixture."""

from __future__ import annotations

import numbers
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class MethylationData:
    """Rows are samples; columns follow unique, ordered CpG IDs.

    ``site_id`` is absent until training rows have been split into sites. Loaders
    must preserve patient identity externally until patient-level splits are done.
    """

    X: np.ndarray
    y: np.ndarray
    cpg_ids: np.ndarray
    site_id: np.ndarray | None = None
    patient_ids: np.ndarray | None = None
    sample_ids: np.ndarray | None = None

    def validate(self) -> None:
        if self.X.ndim != 2 or min(self.X.shape) == 0:
            raise ValueError("X must be a nonempty samples x CpGs matrix")
        if self.X.dtype.kind != "f":
            raise ValueError("X must contain floating-point beta values")
        if np.isinf(self.X).any() or ((self.X < 0) | (self.X > 1)).any():
            raise ValueError("beta values must be in [0, 1] or NaN")
        for name, values, size in (
            ("y", self.y, self.X.shape[0]),
            ("cpg_ids", self.cpg_ids, self.X.shape[1]),
            ("site_id", self.site_id, self.X.shape[0]),
            ("patient_ids", self.patient_ids, self.X.shape[0]),
            ("sample_ids", self.sample_ids, self.X.shape[0]),
        ):
            if values is None and name in ("site_id", "patient_ids", "sample_ids"):
                continue
            if values.ndim != 1 or len(values) != size:
                raise ValueError(f"{name} has incompatible dimensions")
            if values.dtype.kind != "U" or np.any(np.char.strip(values) == ""):
                raise ValueError(f"{name} must contain nonempty Unicode strings")
        if self.sample_ids is not None and len(np.unique(self.sample_ids)) != len(self.sample_ids):
            raise ValueError("sample_ids must be unique")
        if len(np.unique(self.cpg_ids)) != len(self.cpg_ids):
            raise ValueError("cpg_ids must be unique")

    def save(self, path: str | Path) -> None:
        """Write without object arrays or pickle; refuse to overwrite a dataset."""
        self.validate()
        arrays = {"X": self.X, "y": self.y, "cpg_ids": self.cpg_ids}
        for name in ("site_id", "patient_ids", "sample_ids"):
            if getattr(self, name) is not None:
                arrays[name] = getattr(self, name)
        with Path(path).open("xb") as stream:
            np.savez_compressed(stream, **arrays)

    @classmethod
    def load(cls, path: str | Path) -> MethylationData:
        with np.load(path, allow_pickle=False) as arrays:
            data = cls(**{name: arrays[name] for name in arrays.files})
        data.validate()
        return data


def synthetic_data(
    *,
    samples: int = 120,
    cpgs: int = 256,
    classes: int = 3,
    seed: int = 42,
    signal_fraction: float = 1.0,
    class_weights: tuple[float, ...] | None = None,
    profile_noise: float = 0.0,
    missing_fraction: float = 0.0,
) -> MethylationData:
    """Toy classes with disjoint methylated panels; optional imbalance and noise.

    Defaults preserve the original balanced fixture. Optional profile_noise assigns
    some patients another class's profile while keeping their labels unchanged.
    This is declared synthetic discordance, not a biological model.

    This fixture intentionally precedes site assignment, batch effects and
    sparsification. Those belong after the held-out patient split.
    """
    if classes < 2 or samples < classes or cpgs < classes:
        raise ValueError("require classes >= 2, samples >= classes, cpgs >= classes")
    if not np.isfinite(signal_fraction) or not 0 < signal_fraction <= 1:
        raise ValueError("signal_fraction must be in (0, 1]")
    if round(cpgs * signal_fraction) < classes:
        raise ValueError("signal_fraction must provide at least one CpG per class")
    validate_synthetic_options(classes, class_weights, profile_noise, missing_fraction)
    rng = np.random.default_rng(seed)
    labels = np.arange(samples) % classes
    if class_weights is not None:
        proportions = np.asarray(class_weights) / np.sum(class_weights)
        counts = np.floor(samples * proportions).astype(int)
        remainder = samples * proportions - counts
        counts[np.argsort(-remainder, kind="stable")[: samples - counts.sum()]] += 1
        labels = np.repeat(np.arange(classes), counts)
    rng.shuffle(labels)
    X = rng.beta(2, 5, size=(samples, cpgs)).astype(np.float32)
    informative = np.arange(cpgs)
    if signal_fraction < 1:
        informative = rng.choice(cpgs, size=round(cpgs * signal_fraction), replace=False)
    panels = np.array_split(informative, classes)
    # Irreducible profile/label discordance is a declared synthetic noise process.
    # It is applied before splitting, identically for every model, never to predictions.
    profile_labels = labels.copy()
    if profile_noise:
        noisy = rng.random(samples) < profile_noise
        profile_labels[noisy] = (
            profile_labels[noisy] + rng.integers(1, classes, noisy.sum())
        ) % classes
    for label, panel in enumerate(panels):
        rows = np.flatnonzero(profile_labels == label)
        X[np.ix_(rows, panel)] = rng.beta(5, 2, size=(len(rows), len(panel)))
    if missing_fraction:
        X[rng.random(X.shape) < missing_fraction] = np.nan
    data = MethylationData(
        X=X,
        patient_ids=np.array([f"patient_{i:06d}" for i in range(samples)]),
        sample_ids=np.array([f"sample_{i:06d}" for i in range(samples)]),
        y=np.array([f"synthetic_{label}" for label in labels]),
        cpg_ids=np.array([f"synthetic_cpg_{i:06d}" for i in range(cpgs)]),
    )
    data.validate()
    return data


def finite_number(value) -> bool:
    """A real, finite scalar. Strings, None and booleans are rejected rather than raising."""
    return (
        isinstance(value, numbers.Real)
        and not isinstance(value, bool | np.bool_)
        and bool(np.isfinite(value))
    )


def validate_synthetic_options(classes, class_weights, profile_noise, missing_fraction):
    """Validate optional uneven-cohort controls without generating an array."""
    if class_weights is not None:
        try:
            values = np.asarray(class_weights, dtype=float)
        except (TypeError, ValueError):
            values = np.array([np.nan])
        if values.shape != (classes,) or not np.isfinite(values).all() or (values <= 0).any():
            raise ValueError("class_weights must give one finite positive weight per class")
    for name, value in (("profile_noise", profile_noise), ("missing_fraction", missing_fraction)):
        if not finite_number(value) or not 0 <= value < 1:
            raise ValueError(f"{name} must lie in [0, 1)")
