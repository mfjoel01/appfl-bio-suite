"""Unit tests for src.locus_selection (math + stratification only)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.fedfm.locus_selection import (
    compute_maf,
    drop_overlapping_selections,
    frobenius_diff,
    generate_candidate_windows,
    light_ld_prune,
    pairwise_divergence,
    r2_matrix,
    standardize_columns,
    stratify_and_sample,
)


def test_generate_candidate_windows_basic() -> None:
    df = generate_candidate_windows(
        chrom=1, min_bp=1, max_bp=3_000_000, window_size_bp=1_000_000, step_size_bp=500_000
    )
    assert (df["end_bp"] - df["start_bp"] + 1 == 1_000_000).all()
    # Step of 500k across [1, 3M-1M+1] = [1, 2_000_001] in steps of 500k → ~5 windows
    assert len(df) >= 3
    assert df["start_bp"].is_monotonic_increasing


def test_compute_maf_handles_missing() -> None:
    dosage = np.array(
        [
            [0.0, 2.0, 1.0],
            [0.0, 2.0, np.nan],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 1.0],
        ]
    )
    maf = compute_maf(dosage)
    # col 0: mean=0.75 → AF=0.375 → MAF=0.375
    # col 1: mean=1.25 → AF=0.625 → MAF=0.375
    np.testing.assert_allclose(maf[:2], [0.375, 0.375], atol=1e-6)


def test_standardize_columns_handles_monomorphic() -> None:
    X = np.array([[1.0, 5.0], [1.0, 7.0], [1.0, 9.0]])
    Xs, valid = standardize_columns(X)
    assert valid.tolist() == [False, True]
    # Monomorphic column should be all zeros
    assert np.allclose(Xs[:, 0], 0.0)
    # Non-monomorphic column should be unit-variance
    assert np.isclose(Xs[:, 1].std(), 1.0, atol=1e-3)


def test_r2_matrix_diagonal_is_one() -> None:
    rng = np.random.default_rng(0)
    X = rng.normal(size=(200, 5))
    R = r2_matrix(X)
    assert R.shape == (5, 5)
    np.testing.assert_allclose(np.diag(R), 1.0, atol=1e-4)
    # symmetric
    np.testing.assert_allclose(R, R.T, atol=1e-5)


def test_r2_matrix_perfect_correlation() -> None:
    rng = np.random.default_rng(1)
    a = rng.normal(size=100)
    X = np.column_stack([a, a, -a])
    R = r2_matrix(X)
    # All pairwise r² should be 1
    np.testing.assert_allclose(R, np.ones((3, 3)), atol=1e-4)


def test_frobenius_diff_shapes() -> None:
    A = np.eye(3)
    B = np.eye(3)
    assert frobenius_diff(A, B) == 0.0
    B = np.zeros((3, 3))
    assert np.isclose(frobenius_diff(A, B), np.sqrt(3.0))
    with pytest.raises(ValueError):
        frobenius_diff(np.eye(3), np.eye(4))


def test_pairwise_divergence_count() -> None:
    R = {f"s{i}": np.eye(4) for i in range(3)}
    mean, pairs = pairwise_divergence(R)
    assert len(pairs) == 3
    assert mean == 0.0
    # Two of three pairs nonzero
    R["s1"] = np.zeros((4, 4))
    mean2, pairs2 = pairwise_divergence(R)
    assert mean2 > 0


def test_light_ld_prune_drops_collinear() -> None:
    rng = np.random.default_rng(42)
    a = rng.normal(size=200)
    # Build X with col 1 == col 0 (perfectly correlated). Prune should drop one.
    X = np.column_stack([a, a + 1e-6 * rng.normal(size=200), rng.normal(size=200)])
    kept = light_ld_prune(X, window_variants=10, r2_threshold=0.99)
    assert 0 in kept
    assert 1 not in kept  # collinear with 0
    assert 2 in kept


def test_stratify_and_sample_bin_counts() -> None:
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"divergence_score": np.linspace(0, 1, 90), "window_id": [f"w{i}" for i in range(90)]})
    out = stratify_and_sample(df, {"low": 10, "medium": 10, "high": 10}, rng)
    counts = out["stratum"].value_counts().to_dict()
    assert counts == {"low": 10, "medium": 10, "high": 10}


def test_stratify_and_sample_insufficient_raises() -> None:
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"divergence_score": np.linspace(0, 1, 9), "window_id": [f"w{i}" for i in range(9)]})
    with pytest.raises(ValueError, match="windows but"):
        stratify_and_sample(df, {"low": 10, "medium": 10, "high": 10}, rng)


def test_drop_overlapping_selections() -> None:
    # Two overlapping high-stratum windows; the one with higher score wins.
    df = pd.DataFrame(
        {
            "chrom": ["1", "1", "1"],
            "start_bp": [1_000_000, 1_200_000, 5_000_000],
            "end_bp": [2_000_000, 2_200_000, 6_000_000],
            "divergence_score": [0.5, 0.9, 0.8],
            "stratum": ["high", "high", "high"],
        }
    )
    out = drop_overlapping_selections(df)
    assert len(out) == 2
    starts = sorted(out["start_bp"].tolist())
    assert starts == [1_200_000, 5_000_000]  # higher-score overlap kept


def test_drop_overlapping_cross_stratum() -> None:
    # A high-stratum window overlaps a medium-stratum window. The per-stratum
    # pass keeps both (they are in different groups); the global cross-stratum
    # pass must drop the medium one, because the high window is more extreme
    # within its own stratum. This is the case the old per-stratum-only code
    # silently missed (and which had to be patched post-hoc on the real data).
    df = pd.DataFrame(
        {
            "chrom": ["1", "1", "1", "1"],
            #          H            M1(overlaps H)  M2            M_anchor
            "start_bp": [10_000_000, 10_500_000, 20_000_000, 30_000_000],
            "end_bp":   [11_000_000, 11_500_000, 21_000_000, 31_000_000],
            "divergence_score": [100.0, 9.0, 5.1, 5.0],
            "stratum": ["high", "medium", "medium", "medium"],
        }
    )
    out = drop_overlapping_selections(df)
    starts = set(out["start_bp"].tolist())
    assert 10_500_000 not in starts          # overlapping medium window dropped
    assert 10_000_000 in starts              # more-extreme high window survives
    assert {20_000_000, 30_000_000} <= starts  # non-overlapping windows kept
    assert len(out) == 3
    # No two surviving windows overlap.
    rows = out.sort_values("start_bp").to_dict("records")
    for a, b in zip(rows, rows[1:]):
        assert a["end_bp"] < b["start_bp"]
