"""Sweep invariants: fixed pool totals, exact endpoints, equivalent optimization."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest
import torch

from appfl_bio_suite.experiments.methylation.alma3.cohort import cohort
from appfl_bio_suite.experiments.methylation.alma3.config import Config
from appfl_bio_suite.experiments.methylation.alma3.model import select_depth, source
from appfl_bio_suite.experiments.methylation.alma3.privacy import train
from appfl_bio_suite.experiments.methylation.alma3.sweep import configuration, interpolate_counts

SPEC = (
    Path(__file__).parents[1]
    / "src/appfl_bio_suite/experiments/methylation/configs/alma3-skew-sweep.json"
)


def test_interpolation_preserves_site_and_pooled_counts():
    spec = json.loads(SPEC.read_text())
    configuration(spec)
    extreme, target = np.array(spec["extreme_counts"]), np.array(spec["target_counts"])
    for level in np.linspace(0, 1, 101):
        counts = interpolate_counts(extreme, target, level)
        assert (counts >= 1).all()
        np.testing.assert_array_equal(counts.sum(0), extreme.sum(0))
        np.testing.assert_array_equal(counts.sum(1), np.full(5, 100))
        assert np.max(np.abs(counts - ((1 - level) * extreme + level * target))) <= 1 + 1e-8
    np.testing.assert_array_equal(interpolate_counts(extreme, target, 0), extreme)
    np.testing.assert_array_equal(interpolate_counts(extreme, target, 1), np.tile(target, (5, 1)))
    with pytest.raises(ValueError):
        interpolate_counts(extreme, [25, 25, 25, 25], 0.5)


def test_batched_per_patient_clipping_matches_serial():
    torch.set_num_threads(1)
    model, profiles, chrom, positions, _, _ = source(fixture=True)
    model.heads.register_parameter("unused", torch.nn.Parameter(torch.tensor(0.0)))
    select_depth(model, "heads")
    other = copy.deepcopy(model)
    inputs, labels = cohort(
        profiles, chrom, positions, Config(site_sizes=(13, 13)), 41, 0, "train", 0.75
    )
    kwargs = dict(steps=3, lr=0.2, clip=0.05, sigma=0, mu=0.1)
    train(model, inputs, labels, **kwargs)
    train(other, inputs, labels, **kwargs, microbatch_size=4)
    for a, b in zip(model.parameters(), other.parameters(), strict=True):
        torch.testing.assert_close(a, b, atol=3e-6, rtol=1e-6)


def test_patient_reservoir_replicates_are_distinct_and_reproducible():
    _, profiles, chrom, positions, _, _ = source(fixture=True)
    config = Config(
        site_sizes=(12, 12),
        class_indices=(0, 1, 2),
        minimum_per_class=1,
        site_class_proportions=((1 / 3,) * 3,) * 2,
    )
    a, y = cohort(profiles, chrom, positions, config, 41, 0, "train", 0.75, replicate=1)
    b, z = cohort(profiles, chrom, positions, config, 41, 0, "train", 0.75, replicate=1)
    c, _ = cohort(profiles, chrom, positions, config, 41, 0, "train", 0.75, replicate=2)
    assert torch.equal(a[0], b[0]) and torch.equal(y, z)
    assert not torch.equal(a[0], c[0])
    np.testing.assert_array_equal(np.bincount(y), [4, 4, 4])
