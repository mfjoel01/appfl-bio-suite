"""Aggregate-only, sequential hierarchy calibration on a separate validation fold."""

import math

import numpy as np
import torch

from .privacy import private_generator


def path_scores(logp, paths, taxonomy=None):
    if taxonomy is not None:
        from alma3.terminal import selected_path_scores

        _, chosen, scores, _ = selected_path_scores(logp.exp(), taxonomy)
        return chosen.cpu().numpy(), scores.cpu().numpy()
    probability = logp.exp().cpu().numpy()
    chosen = probability.argmax(1)
    return paths[chosen], probability.max(1)[:, None]


def accepted(chosen, scores, thresholds):
    live = np.ones(len(chosen), dtype=bool)
    decisions = np.zeros(chosen.shape, dtype=bool)
    for level, threshold in enumerate(thresholds):
        live = live & (chosen[:, level] >= 0) & (scores[:, level] >= threshold)
        decisions[:, level] = live
    return decisions


def histogram(chosen, scores, truth, *, level, prior, bins, sigma):
    eligible = chosen[:, level] >= 0
    if level:
        eligible &= accepted(chosen, scores, prior)[:, level - 1]
    bucket = np.minimum((scores[:, level].clip(0, 1) * bins).astype(int), bins - 1)
    correct = chosen[:, level] == truth[:, level]
    counts = np.zeros((2, bins), dtype=np.float64)
    np.add.at(counts, (correct[eligible].astype(int), bucket[eligible]), 1)
    if sigma:
        # Replace-one changes at most two disjoint histogram cells.
        counts += (
            torch.randn(counts.shape, generator=private_generator(), dtype=torch.float64).numpy()
            * math.sqrt(2)
            * sigma
        )
    return counts


def select_threshold(counts, *, floor, minimum, noise_std=0.0):
    """Wilson lower bound, with simultaneous Gaussian count-error margins.

    Abstention is 1.01 (finite JSON), including when all bins are unsupported.
    This is validation evidence, not a clinical population guarantee.
    """
    counts = np.asarray(counts)
    if counts.ndim != 2 or counts.shape[0] != 2 or not np.isfinite(counts).all():
        raise ValueError("invalid calibration histogram")
    bins = counts.shape[1]
    cumulative = np.cumsum(counts[:, ::-1], axis=1)[:, ::-1]
    # Union bound over two count tails at each candidate threshold (99% noise bound).
    z_noise = math.sqrt(2 * math.log(4 * bins / 0.01))
    for i in range(bins):
        margin = z_noise * noise_std * math.sqrt(bins - i)
        correct = max(0.0, cumulative[1, i] - margin)
        incorrect = max(0.0, cumulative[0, i] + margin)
        lower_total = max(0.0, cumulative[:, i].sum() - 2 * margin)
        total = correct + incorrect
        if lower_total < minimum or total == 0:
            continue
        p, z = correct / total, 1.959963984540054
        lower = (
            p
            + z * z / (2 * total)
            - z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total))
        ) / (1 + z * z / total)
        if lower >= floor:
            return i / bins
    return 1.01
