"""Independent synthetic patients, fixed site membership and disjoint splits.

Profiles never leave the worker. Every patient has one measurement. Validation
and test use separate RNG streams; validation-only never generates test inputs.
"""

import numpy as np
import torch


def class_counts(n, proportions, minimum):
    """Round explicit target ratios while retaining every class in each split."""
    target = n * np.asarray(proportions, dtype=np.float64)
    if n < minimum * len(target):
        raise ValueError("cohort too small for requested class support")
    counts = np.maximum(np.floor(target).astype(np.int64), minimum)
    while counts.sum() < n:
        counts[np.argmax(target - counts)] += 1
    while counts.sum() > n:
        eligible = np.where(counts > minimum, counts - target, -np.inf)
        counts[np.argmax(eligible)] -= 1
    return counts


def rotate_signature(logits, profiles, support, degrees):
    """Invertible, label-blind rotation of two template-contrast directions.

    A deliberate synthetic stressor, not a validated biological batch model.
    The transform is fixed from source templates and never fitted to patients.
    """
    means = np.clip(profiles[support], 1e-4, 1 - 1e-4)
    reference = np.log(means / (1 - means))
    center = reference.mean(0)
    contrast = reference - center
    eigenvalues, eigenvectors = np.linalg.eigh(contrast @ contrast.T)
    if eigenvalues[-2] <= 1e-6:
        raise ValueError("signature rotation requires two independent template contrasts")
    basis = (eigenvectors[:, -2:].T @ contrast) / np.sqrt(eigenvalues[-2:])[:, None]
    # Resolve arbitrary eigenvector signs for reproducible rotation orientation.
    for direction in basis:
        if direction[np.argmax(np.abs(direction))] < 0:
            direction *= -1
    angle = np.deg2rad(degrees)
    rotation = np.asarray(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]], dtype=np.float32
    )
    scores = (logits - center) @ basis.T
    return logits + (scores @ rotation.T - scores) @ basis


def nanopore_measurement(beta, observed, rng, site_fraction):
    """Synthetic read-count observation; parameters are illustrative, not fitted.

    Neutral beta=0.5 background is a proxy, not a normal-cell reference.
    Site ordering denotes technical profiles, never countries or ancestry.
    """
    purity = rng.uniform(0.55, 0.95, (len(beta), 1))
    latent = purity * beta + (1 - purity) * 0.5
    false_positive = 0.02 + 0.06 * site_fraction
    false_negative = 0.08 + 0.06 * site_fraction
    probability = false_positive + (1 - false_positive - false_negative) * latent
    depth = rng.poisson(2 + 6 * (1 - site_fraction), beta.shape)
    observed = observed & (depth > 0)
    calls = rng.binomial(depth, probability)
    measured = np.divide(calls, depth, out=np.full(beta.shape, 0.5), where=depth > 0)
    return measured.astype(np.float32), observed


def cohort(
    profiles,
    chrom,
    positions,
    config,
    seed,
    site,
    split,
    coverage,
    *,
    unshifted=False,
    replicate=None,
):
    codes = {"train": 11, "calibration": 22, "validation": 33, "test": 44}
    if split not in codes:
        raise ValueError("unknown cohort split")
    n = (
        config.site_sizes[site]
        if split == "train"
        else (config.test_size if split == "test" else config.validation_size)
    )
    classes, width = profiles.shape
    support = np.asarray(config.class_indices or tuple(range(classes)), dtype=np.int64)
    if np.any(support >= classes):
        raise ValueError("class index exceeds released taxonomy")
    site_rng = np.random.default_rng(np.random.SeedSequence([seed, site, 100]))
    weights = site_rng.dirichlet(np.full(len(support), config.label_alpha))
    # Site mechanism fixed across splits and coverages. Each site retains support.
    weights = 0.85 * weights + 0.15 / len(support)
    shift = site_rng.normal(0, config.batch_shift, width).astype(np.float32)
    drift = site_rng.normal(0, 1, (classes, min(width, 256))).astype(np.float32)
    patient_seed = [seed, site, codes[split]]
    if replicate is not None:
        patient_seed.append(replicate)
    rng = np.random.default_rng(np.random.SeedSequence(patient_seed))
    if config.site_class_proportions:
        counts = class_counts(n, config.site_class_proportions[site], config.minimum_per_class)
        y = np.repeat(support, counts)
        rng.shuffle(y)
    else:
        guaranteed = np.repeat(support, config.minimum_per_class)
        y = np.concatenate((guaranteed, rng.choice(support, n - len(guaranteed), p=weights)))
        if config.minimum_per_class:
            rng.shuffle(y)
    means = np.clip(profiles[y], 1e-4, 1 - 1e-4)
    logits = np.log(means / (1 - means))
    if not unshifted:
        if config.signature_rotation_degrees:
            # Related but distinct site transforms, fixed across all patient splits.
            angle = config.signature_rotation_degrees * (
                0.95 + 0.1 * site / (len(config.site_sizes) - 1)
            )
            logits = rotate_signature(logits, profiles, support, angle)
        logits += shift
        logits[:, : drift.shape[1]] += config.subtype_drift * drift[y]
    beta = (1 / (1 + np.exp(-logits))).astype(np.float32)
    # Beta-distributed within-class variation, no labels altered post-generation.
    beta = rng.beta(np.clip(beta * 40, 0.01, None), np.clip((1 - beta) * 40, 0.01, None)).astype(
        np.float32
    )
    observed = np.zeros((n, width), dtype=bool)
    # Different site coverage; unshifted data uses the requested coverage itself.
    rate = coverage if unshifted else coverage * (1 - 0.4 * site / len(config.site_sizes))
    count = max(1 if width < 5000 else 5000, int(width * rate))
    for i in range(n):
        observed[i, rng.choice(width, min(width, count), replace=False)] = True
    if config.nanopore_like and not unshifted:
        beta, observed = nanopore_measurement(
            beta, observed, rng, site / max(1, len(config.site_sizes) - 1)
        )
    inputs = (
        torch.from_numpy(beta),
        torch.from_numpy(observed),
        torch.zeros_like(torch.from_numpy(beta)),
        chrom.expand(n, -1),
        positions.expand(n, -1),
    )
    return inputs, torch.from_numpy(y).long()
