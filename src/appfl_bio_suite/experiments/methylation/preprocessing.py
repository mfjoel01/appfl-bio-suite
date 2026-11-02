"""Patient-grouped holdouts, training-only panel selection, and simulated hospitals."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .dataset import MethylationData


@dataclass
class PreparedData:
    X: np.ndarray
    y: np.ndarray
    classes: np.ndarray
    train: np.ndarray
    validation: np.ndarray
    test: np.ndarray
    sites: list[np.ndarray]
    panel: np.ndarray
    label_mapping: dict[str, str]


def prepare(
    data: MethylationData,
    *,
    panel_size: int,
    seed: int,
    min_class: int = 20,
    alpha: float = 0.5,
    small_fraction: float = 0.1,
    min_site_class_patients: int = 0,
) -> PreparedData:
    data.validate()
    if data.patient_ids is None or data.sample_ids is None:
        raise ValueError("training requires explicit patient_ids and unique sample_ids")
    if data.site_id is not None:
        raise ValueError("this simulation assigns sites; provide a cohort without site_id")
    if panel_size < 1 or min_class < 1 or alpha <= 0 or not 0 < small_fraction < 0.5:
        raise ValueError("invalid panel size, minimum class size, alpha or small-site fraction")
    rng = np.random.default_rng(seed)
    patients, inverse = np.unique(data.patient_ids, return_inverse=True)
    labels = []
    for patient in range(len(patients)):
        values = np.unique(data.y[inverse == patient])
        if len(values) != 1:
            raise ValueError("each patient must have a single subtype label")
        labels.append(values[0])
    labels = np.array(labels)
    groups = {key: [] for key in ("train", "validation", "test")}
    for label in np.unique(labels):
        ids = rng.permutation(np.flatnonzero(labels == label))
        if len(ids) < 3:
            raise ValueError(f"subtype {label!r} needs at least three patients for holdouts")
        n_test = max(1, round(0.2 * len(ids)))
        n_val = max(1, round(0.1 * (len(ids) - n_test)))
        groups["test"].extend(ids[:n_test])
        groups["validation"].extend(ids[n_test : n_test + n_val])
        groups["train"].extend(ids[n_test + n_val :])
    train_groups = np.array(groups["train"])
    counts = {label: int(np.sum(labels[train_groups] == label)) for label in np.unique(labels)}
    mapping = {label: label if count >= min_class else "other" for label, count in counts.items()}
    mapped = np.array([mapping[label] for label in data.y])
    classes, y = np.unique(mapped, return_inverse=True)
    if len(classes) < 2:
        raise ValueError("training class threshold leaves fewer than two classes")
    splits = {key: np.flatnonzero(np.isin(inverse, ids)) for key, ids in groups.items()}
    train_X = data.X[splits["train"]]
    count = np.isfinite(train_X).sum(axis=0)
    sums = np.nansum(train_X, axis=0, dtype=np.float64)
    squares = np.nansum(train_X.astype(np.float64) ** 2, axis=0)
    variance = squares / np.maximum(count, 1) - (sums / np.maximum(count, 1)) ** 2
    eligible = np.flatnonzero((count >= 2) & (variance > 0))
    if not len(eligible):
        raise ValueError("no variable CpGs with at least two training observations")
    panel = eligible[np.argsort(-variance[eligible], kind="stable")[:panel_size]]
    # Capacities constrain patient counts; Dirichlet preferences provide label skew.
    n = len(train_groups)
    if n < 4:
        raise ValueError("four sites require at least four training patients")
    small = max(1, min(n - 3, round(n * small_fraction)))
    capacities = np.array([(n - small) // 3] * 3 + [small])
    capacities[: (n - small) % 3] += 1
    preferences = rng.dirichlet(np.full(4, alpha), size=len(classes))
    assignments = [[] for _ in range(4)]
    if min_site_class_patients < 0:
        raise ValueError("min_site_class_patients must be nonnegative")
    reserved = set()
    if min_site_class_patients:
        if np.any(capacities < min_site_class_patients * len(classes)):
            raise ValueError("site capacity is too small for the requested class minimum")
        for label in classes:
            candidates = rng.permutation([p for p in train_groups if mapping[labels[p]] == label])
            required = 4 * min_site_class_patients
            if len(candidates) < required:
                raise ValueError(f"not enough training patients of class {label!r} for all sites")
            for site in range(4):
                selected = candidates[
                    site * min_site_class_patients : (site + 1) * min_site_class_patients
                ]
                assignments[site].extend(selected)
                reserved.update(selected.tolist())
                capacities[site] -= len(selected)
    remaining = [patient for patient in train_groups if patient not in reserved]
    for patient in rng.permutation(remaining):
        label = np.flatnonzero(classes == mapping[labels[patient]])[0]
        weights = preferences[label] * capacities
        site = int(rng.choice(4, p=weights / weights.sum()))
        assignments[site].append(patient)
        capacities[site] -= 1
    sites = [np.flatnonzero(np.isin(inverse, ids)) for ids in assignments]
    return PreparedData(
        data.X[:, panel].copy(),
        y,
        classes,
        splits["train"],
        splits["validation"],
        splits["test"],
        sites,
        panel,
        mapping,
    )


def simulate_calls(
    X: np.ndarray,
    rng: np.random.Generator,
    *,
    observed: int | None = None,
    min_fraction: float = 0.015,
    intercept: float = -0.155,
    slope: float = 1.517,
) -> np.ndarray:
    """Sample biased binary calls; NaNs remain uncovered.

    Uniform selection of observed probes is a PoC approximation, not a read-depth
    model. Training draws coverage log-uniformly to emphasize sparse input.
    """
    if not 0 < min_fraction <= 1 or (observed is not None and observed < 1):
        raise ValueError("coverage must be positive")
    calls = np.zeros(X.shape, dtype=np.float32)
    for i, row in enumerate(X):
        available = np.flatnonzero(np.isfinite(row))
        target = observed
        if target is None:
            target = max(1, round(X.shape[1] * np.exp(rng.uniform(np.log(min_fraction), 0))))
        selected = rng.choice(available, size=min(target, len(available)), replace=False)
        beta = np.clip(row[selected], 1e-6, 1 - 1e-6)
        probability = 1 / (1 + np.exp(-(intercept + slope * np.log(beta / (1 - beta)))))
        calls[i, selected] = np.where(rng.random(len(selected)) < probability, 1, -1)
    return calls
