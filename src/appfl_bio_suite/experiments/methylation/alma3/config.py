"""Explicit study configuration, frozen before test evaluation."""

import math
from dataclasses import asdict, dataclass


@dataclass
class Config:
    seeds: tuple[int, ...] = (42, 43, 44)
    depths: tuple[str, ...] = ("heads", "trunk")
    noise_multipliers: tuple[float, ...] = (0.0, 2.0, 8.0)
    site_sizes: tuple[int, ...] = (48, 48, 48, 16)
    validation_size: int = 32
    test_size: int = 32
    rounds: int = 3
    local_steps: int = 2
    learning_rate: float = 0.1
    proximal_mu: float = 0.1
    clip_norm: float = 1.0
    delta: float = 1e-5
    coverage: tuple[float, ...] = (0.25, 0.75)
    batch_shift: float = 0.6
    signature_rotation_degrees: float = 0.0
    subtype_drift: float = 0.15
    label_alpha: float = 0.5
    bins: int = 20
    precision_floor: float = 0.9
    min_calibration: int = 20
    secure: bool = True
    threads: int = 1
    release_ensemble: bool = False
    class_indices: tuple[int, ...] = ()
    minimum_per_class: int = 0
    site_class_proportions: tuple[tuple[float, ...], ...] = ()
    arms: tuple[str, ...] = ()

    def validate(self):
        if type(self.minimum_per_class) is not int or self.minimum_per_class < 0:
            raise ValueError("minimum_per_class must be a nonnegative integer")
        if len(set(self.class_indices)) != len(self.class_indices) or any(
            type(c) is not int or c < 0 for c in self.class_indices
        ):
            raise ValueError("class_indices must be unique nonnegative integers")
        if self.minimum_per_class and not self.class_indices:
            raise ValueError("minimum_per_class requires explicit class_indices")
        if self.site_class_proportions:
            if not self.class_indices or self.minimum_per_class < 1:
                raise ValueError(
                    "explicit site proportions require classes and minimum_per_class >= 1"
                )
            if len(self.site_class_proportions) != len(self.site_sizes):
                raise ValueError("require one class-proportion row per site")
            for row in self.site_class_proportions:
                if (
                    len(row) != len(self.class_indices)
                    or any(not math.isfinite(v) or v <= 0 for v in row)
                    or not math.isclose(sum(row), 1.0, abs_tol=1e-8)
                ):
                    raise ValueError("class-proportion rows must be positive and sum to one")
        required = self.minimum_per_class * len(self.class_indices)
        if any(n < required for n in (*self.site_sizes, self.validation_size, self.test_size)):
            raise ValueError("each split must accommodate minimum_per_class")
        allowed = {
            "zero-shot",
            "pooled",
            "fedavg",
            "fedprox",
            "scratch-fedavg",
            "template-zero-shot",
            "ensemble-zero-shot",
        }
        allowed.update(f"local-S{i}" for i in range(len(self.site_sizes)))
        if len(set(self.arms)) != len(self.arms) or set(self.arms) - allowed:
            raise ValueError("invalid or duplicate study arms")
        if "ensemble-zero-shot" in self.arms and not self.release_ensemble:
            raise ValueError("ensemble arm requires release_ensemble")
        for name in (
            "validation_size",
            "test_size",
            "rounds",
            "local_steps",
            "bins",
            "min_calibration",
            "threads",
        ):
            if type(getattr(self, name)) is not int or getattr(self, name) < 1:
                raise ValueError(f"{name} must be a positive integer")
        if len(self.site_sizes) < 2 or any(type(n) is not int or n < 2 for n in self.site_sizes):
            raise ValueError("require at least two sites with at least two patients")
        if not self.seeds or len(set(self.seeds)) != len(self.seeds):
            raise ValueError("seeds must be nonempty and unique")
        if any(type(s) is not int or not 0 <= s < 2**31 for s in self.seeds):
            raise ValueError("invalid seed")
        if not self.depths or len(set(self.depths)) != len(self.depths):
            raise ValueError("depths must be nonempty and unique")
        if set(self.depths) - {"heads", "trunk", "full"}:
            raise ValueError("depth must be heads, trunk or full")
        for name in (
            "learning_rate",
            "proximal_mu",
            "clip_norm",
            "delta",
            "batch_shift",
            "signature_rotation_degrees",
            "subtype_drift",
            "label_alpha",
            "precision_floor",
        ):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f"nonfinite {name}")
        if not 0 <= self.signature_rotation_degrees <= 180:
            raise ValueError("signature_rotation_degrees must be in [0,180]")
        if self.signature_rotation_degrees and len(self.class_indices) < 3:
            raise ValueError("signature rotation requires at least three explicit classes")
        if min(self.learning_rate, self.clip_norm, self.label_alpha) <= 0:
            raise ValueError("learning rate, clipping and alpha must be positive")
        if self.proximal_mu < 0 or self.batch_shift < 0 or not 0 <= self.subtype_drift <= 1:
            raise ValueError("invalid shift or proximal strength")
        if not 0 < self.delta < 1 or not 0 < self.precision_floor < 1:
            raise ValueError("invalid delta or precision floor")
        if not self.coverage or any(not 0 < c <= 1 for c in self.coverage):
            raise ValueError("coverage fractions must be in (0,1]")
        if len(set(self.coverage)) != len(self.coverage):
            raise ValueError("duplicate coverage")
        if (
            not self.noise_multipliers
            or any(not math.isfinite(n) or n < 0 for n in self.noise_multipliers)
            or len(set(self.noise_multipliers)) != len(self.noise_multipliers)
        ):
            raise ValueError("noise multipliers must be finite, nonnegative and unique")
        if type(self.secure) is not bool or type(self.release_ensemble) is not bool:
            raise ValueError("secure must be Boolean")

    def json(self):
        return asdict(self)
