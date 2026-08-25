"""Unit tests for src.sampling."""

from __future__ import annotations

import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.fedfm.sampling import SITE_ORDER, assign_individuals


def _mock_manifest(per_pop: int = 100, pops=("EUR", "AFR", "AMR", "EAS", "CSA", "MID")) -> pd.DataFrame:
    rows = []
    for pop in pops:
        for i in range(per_pop):
            rows.append({"FID": f"{pop}_{i:04d}", "IID": f"{pop}_{i:04d}", "superpopulation": pop})
    return pd.DataFrame(rows)


def _toy_sites() -> dict[str, dict[str, int]]:
    return {
        "anl": {"EUR": 50, "AFR": 10, "AMR": 10, "EAS": 5, "CSA": 5},
        "covenant": {"AFR": 80, "EUR": 3, "CSA": 2},
        "mbzuai": {"MID": 50, "CSA": 30, "AFR": 5, "EUR": 5},
    }


def test_no_cross_site_overlap() -> None:
    manifest = _mock_manifest(per_pop=150)
    sites = _toy_sites()
    assignments = assign_individuals(manifest, sites, master_seed=42)

    seen: set[tuple[str, str]] = set()
    for site_name, df in assignments.items():
        pairs = set(zip(df["FID"], df["IID"]))
        overlap = pairs & seen
        assert not overlap, f"Site {site_name} overlaps prior sites at {overlap}"
        seen |= pairs


def test_exact_counts() -> None:
    manifest = _mock_manifest(per_pop=150)
    sites = _toy_sites()
    assignments = assign_individuals(manifest, sites, master_seed=42)

    for site_name, composition in sites.items():
        df = assignments[site_name]
        assert len(df) == sum(composition.values())
        counts = df["superpopulation"].value_counts().to_dict()
        for pop, expected in composition.items():
            assert counts.get(pop, 0) == expected, (
                f"Site {site_name} pop {pop}: expected {expected}, got {counts.get(pop, 0)}"
            )


def test_determinism_same_seed() -> None:
    manifest = _mock_manifest(per_pop=150)
    sites = _toy_sites()
    a = assign_individuals(manifest, sites, master_seed=42)
    b = assign_individuals(manifest, sites, master_seed=42)
    for site in SITE_ORDER:
        pd.testing.assert_frame_equal(
            a[site].reset_index(drop=True), b[site].reset_index(drop=True)
        )


def test_different_seeds_differ() -> None:
    manifest = _mock_manifest(per_pop=150)
    sites = _toy_sites()
    a = assign_individuals(manifest, sites, master_seed=1)
    b = assign_individuals(manifest, sites, master_seed=2)
    # At least one site should differ.
    differs = False
    for site in SITE_ORDER:
        if not a[site].equals(b[site]):
            differs = True
            break
    assert differs, "Expected different seeds to produce different assignments"


def test_insufficient_pool_raises() -> None:
    manifest = _mock_manifest(per_pop=10)  # too few
    sites = _toy_sites()
    with pytest.raises(ValueError, match="Insufficient individuals"):
        assign_individuals(manifest, sites, master_seed=42)


def test_site_ordering_deterministic() -> None:
    # Verify that assignments dict iteration follows SITE_ORDER.
    manifest = _mock_manifest(per_pop=150)
    sites = _toy_sites()
    assignments = assign_individuals(manifest, sites, master_seed=42)
    assert list(assignments.keys()) == list(SITE_ORDER)
