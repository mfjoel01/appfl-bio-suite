"""Validation-only calibration and multiclass evaluation with explicit undefined values."""

from __future__ import annotations

import numpy as np
import torch
from sklearn.metrics import confusion_matrix, f1_score, recall_score


def fit_temperature(logits: np.ndarray, labels: np.ndarray) -> float:
    # Bounded grid: deterministic, robust for tiny validation cohorts, includes T=1.
    values = torch.as_tensor(logits, dtype=torch.float64)
    targets = torch.as_tensor(labels, dtype=torch.long)
    candidates = np.unique(np.append(np.geomspace(0.05, 20, 161), 1.0))
    losses = [float(torch.nn.functional.cross_entropy(values / t, targets)) for t in candidates]
    return float(candidates[np.argmin(losses)])


def metrics(
    logits: np.ndarray,
    y: np.ndarray,
    *,
    temperature: float = 1,
    threshold: float = 0.9,
    bins: int = 15,
) -> dict:
    probabilities = torch.softmax(torch.as_tensor(logits) / temperature, dim=1).numpy()
    predictions = probabilities.argmax(axis=1)
    confidence = probabilities.max(axis=1)
    correct = predictions == y
    confident = confidence >= threshold
    labels = np.arange(logits.shape[1])
    recalls = recall_score(y, predictions, labels=labels, average=None, zero_division=0)
    present = np.bincount(y, minlength=len(labels)) > 0
    bucket = np.minimum((confidence * bins).astype(int), bins - 1)
    ece = 0.0
    for i in range(bins):
        mask = bucket == i
        if mask.any():
            ece += mask.mean() * abs(correct[mask].mean() - confidence[mask].mean())
    return {
        "accuracy": float(correct.mean()),
        "balanced_accuracy": float(recalls[present].mean()),
        "macro_f1": float(
            f1_score(y, predictions, labels=labels, average="macro", zero_division=0)
        ),
        "per_class_recall": [
            float(r) if p else None for r, p in zip(recalls, present, strict=True)
        ],
        "ece": float(ece),
        "confident_call_rate": float(confident.mean()),
        "confident_accuracy": float(correct[confident].mean()) if confident.any() else None,
        "confusion_matrix": confusion_matrix(y, predictions, labels=labels).tolist(),
        "n_test": len(y),
    }
