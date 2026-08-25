"""Unit tests for src.utils — config loading, seed determinism, .bed reader."""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
import pytest
import yaml

from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import (
    derive_seed,
    load_config,
    make_rng,
    read_bed_variants,
)


def test_derive_seed_deterministic() -> None:
    s1 = derive_seed(42, "site_a", 7)
    s2 = derive_seed(42, "site_a", 7)
    # Generate the same sequence of numbers from both.
    a = np.random.default_rng(s1).integers(0, 1 << 32, size=5)
    b = np.random.default_rng(s2).integers(0, 1 << 32, size=5)
    np.testing.assert_array_equal(a, b)


def test_derive_seed_different_components_differ() -> None:
    s1 = derive_seed(42, "site_a")
    s2 = derive_seed(42, "site_b")
    a = np.random.default_rng(s1).integers(0, 1 << 32, size=5)
    b = np.random.default_rng(s2).integers(0, 1 << 32, size=5)
    assert not np.array_equal(a, b)


def test_make_rng_independent_streams() -> None:
    r1 = make_rng(42, "a")
    r2 = make_rng(42, "b")
    x1 = r1.standard_normal(100)
    x2 = r2.standard_normal(100)
    assert not np.allclose(x1, x2)


def _write_tiny_bed(path_prefix: Path, dosages: np.ndarray) -> None:
    """Write a tiny PLINK 1 .bed at ``path_prefix.bed`` encoding ``dosages``.

    dosages: (n_samples, n_variants) with values in {0, 1, 2} (no NaN here for the test).
    """
    n_samples, n_variants = dosages.shape
    # Reverse map: dosage 2 → code 00, 1 → 10, 0 → 11.
    code_map = {2: 0b00, 1: 0b10, 0: 0b11}
    bed_path = path_prefix.with_suffix(".bed")
    with bed_path.open("wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        bytes_per_variant = (n_samples + 3) // 4
        for v in range(n_variants):
            col = dosages[:, v]
            codes = np.array([code_map[int(d)] for d in col], dtype=np.uint8)
            # Pad to full bytes
            pad = bytes_per_variant * 4 - n_samples
            if pad:
                codes = np.concatenate([codes, np.zeros(pad, dtype=np.uint8)])
            block = np.zeros(bytes_per_variant, dtype=np.uint8)
            for k in range(4):
                block |= (codes[k::4] & 0x03) << (k * 2)
            fh.write(block.tobytes())


def test_bed_reader_roundtrip(tmp_path: Path) -> None:
    n_samples, n_variants = 7, 4
    rng = np.random.default_rng(0)
    dosages = rng.integers(0, 3, size=(n_samples, n_variants)).astype(np.float32)
    prefix = tmp_path / "toy"
    _write_tiny_bed(prefix, dosages)
    out = read_bed_variants(prefix, list(range(n_variants)), n_samples)
    np.testing.assert_array_equal(out, dosages)


def test_load_config_validates(tmp_path: Path) -> None:
    cfg_path = tmp_path / "config" / "cfg.yaml"
    cfg_path.parent.mkdir()
    body = {
        "paths": {
            "hapnest_dir": "data/raw",
            "population_manifest": "data/raw/p.tsv",
            "processed_dir": "data/proc",
            "loci_dir": "data/loci",
            "ground_truth_dir": "data/gt",
            "reports_dir": "reports",
            "logs_dir": "logs",
        },
        "master_seed": 1,
        "chromosome": 1,
        "tools": {"plink": "plink", "plink2": "plink2", "king": "king"},
        "sites": {
            "anl": {"n": 10, "composition": {"EUR": 10}, "dominant": "EUR"},
        },
        "superpopulations": ["EUR", "AFR"],
        "locus_selection": {
            "window_size_bp": 1000, "step_size_bp": 500, "n_loci": 6,
            "strata": {"low": 2, "medium": 2, "high": 2},
            "maf_filter": 0.01, "ld_tag_snp_target": 50, "ld_r2_prune_threshold": 0.99,
            "ld_prune_window_variants": 10, "ld_prune_step_variants": 5, "n_workers": 1,
        },
        "architecture": {
            "ncsl": [1], "h2": [0.001], "rg": [1.0], "factorial_mode": "minimal",
            "replicates": 1,
            "ancestry_specific_causal": {
                "enabled": False, "min_per_stratum": 0,
                "common_maf_threshold": 0.05, "rare_maf_threshold": 0.01,
            },
        },
        "phenotype": {"model": "linear_additive", "h2_tolerance_relative": 0.2},
        "qc": {
            "kinship_threshold": 0.05, "pca_n_components": 5,
            "pca_reference": "1kg", "maf_tolerance_sd": 3.0,
            "ld_decay_max_dist_kb": 500, "ld_decay_n_loci_sample": 1,
        },
    }
    cfg_path.write_text(yaml.safe_dump(body))
    cfg = load_config(cfg_path)
    assert cfg.master_seed == 1
    assert cfg.sites["anl"].n == 10
    # Strata sum mismatch → ValueError.
    body["locus_selection"]["strata"]["high"] = 3
    cfg_path.write_text(yaml.safe_dump(body))
    with pytest.raises(Exception):
        load_config(cfg_path)
