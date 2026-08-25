"""Unit tests for src.phenotype_sim.

The critical test (per the spec) is the β-by-superpopulation lookup: every
individual's genetic value must use the β indexed by their *own* superpopulation,
not a site-level average.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.fedfm.phenotype_sim import (
    Architecture,
    _compute_genetic_value,
    _find_ancestry_specific_candidates,
    assemble_divergent_beta,
    build_architecture_grid,
    draw_effect_sizes,
    select_ancestry_divergent_causal_sets,
    select_causal_variants,
    simulate_phenotype_for_site,
)


def test_architecture_grid_minimal() -> None:
    grid = build_architecture_grid(
        ncsl=[1, 2, 3], h2=[0.0005, 0.001, 0.005], rg=[0.5, 0.7, 1.0],
        mode="minimal",
    )
    assert len(grid) == 9
    # All at rg=1.0
    assert {a.rg for a in grid} == {1.0}


def test_architecture_grid_extended() -> None:
    grid = build_architecture_grid(
        ncsl=[1, 2, 3], h2=[0.0005, 0.001, 0.005], rg=[0.5, 0.7, 1.0],
        mode="extended",
    )
    # minimal (9) + ncsl×rg at h²=0.001 (9, of which 3 are dupes of rg=1.0 + h²=0.001)
    assert len(grid) == 15


def test_architecture_grid_full() -> None:
    grid = build_architecture_grid(
        ncsl=[1, 2, 3], h2=[0.0005, 0.001, 0.005], rg=[0.5, 0.7, 1.0],
        mode="full",
    )
    assert len(grid) == 27


def test_beta_by_superpopulation_lookup() -> None:
    """The critical correctness test.

    Build a toy cohort of 10 individuals across 3 populations. Assign known
    β values: pop 0 → β=10, pop 1 → β=20, pop 2 → β=30. Single causal variant,
    dosage=1 for everyone. The genetic value for each individual must equal
    their population's β, not the mean (=20) and not any site-level average.
    """
    n_indiv = 9
    # Pop assignment: 3 from each pop.
    pop_index = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2], dtype=np.int64)
    # Dosage matrix: (n_samples, ncsl=1), all 1.0
    dosage = np.ones((n_indiv, 1), dtype=np.float32)
    # β: (ncsl=1, n_pops=3) — pop 0 → 10, pop 1 → 20, pop 2 → 30
    beta = np.array([[10.0, 20.0, 30.0]], dtype=np.float32)

    g = _compute_genetic_value(dosage, beta, pop_index)
    expected = np.array([10, 10, 10, 20, 20, 20, 30, 30, 30], dtype=np.float32)
    np.testing.assert_allclose(g, expected)


def test_beta_lookup_dosage_zero_individuals_get_zero() -> None:
    pop_index = np.array([0, 1, 2], dtype=np.int64)
    dosage = np.zeros((3, 1), dtype=np.float32)
    beta = np.array([[5.0, 7.0, 11.0]], dtype=np.float32)
    g = _compute_genetic_value(dosage, beta, pop_index)
    np.testing.assert_allclose(g, np.zeros(3))


def test_draw_effect_sizes_correlation_structure() -> None:
    """For rg=0.7, the empirical between-population correlation should approach 0.7
    for a large draw."""
    rng = np.random.default_rng(0)
    n = 50_000
    beta = draw_effect_sizes(n, 3, 0.7, rng)
    # Each row is a multivariate draw across 3 populations.
    # Empirical correlation across rows between cols 0 and 1.
    c = np.corrcoef(beta[:, 0], beta[:, 1])[0, 1]
    assert abs(c - 0.7) < 0.02, f"Expected ~0.7, got {c}"


def test_draw_effect_sizes_rg1() -> None:
    """rg=1.0 → identical β across populations (rank-1 Σ)."""
    rng = np.random.default_rng(0)
    beta = draw_effect_sizes(100, 4, 1.0, rng)
    # All cols equal, within jitter.
    for k in range(1, 4):
        np.testing.assert_allclose(beta[:, 0], beta[:, k], atol=1e-3)


def test_select_causal_basic() -> None:
    variants = pd.DataFrame({"snp_id": [f"rs{i}" for i in range(100)]})
    rng = np.random.default_rng(0)
    idx, asc = select_causal_variants(variants, ncsl=3, rng=rng)
    assert len(idx) == 3
    assert len(set(idx.tolist())) == 3
    assert not asc


def test_select_causal_ancestry_specific_forced() -> None:
    variants = pd.DataFrame({"snp_id": [f"rs{i}" for i in range(20)]})
    rng = np.random.default_rng(0)
    mask = np.zeros(20, dtype=bool)
    mask[5] = True   # only valid candidate
    idx, asc = select_causal_variants(variants, ncsl=2, rng=rng, ancestry_specific_mask=mask)
    assert asc
    assert 5 in idx


def test_ancestry_divergent_sets_structure() -> None:
    """Divergent sets: n_shared shared across all pops, n_private per pop; the
    union has the expected size and per-pop sets genuinely differ."""
    variants = pd.DataFrame({"snp_id": [f"rs{i}" for i in range(50)]})
    rng = np.random.default_rng(1)
    ncsl, n_priv, n_pops = 3, 1, 6
    union, per_pop, used = select_ancestry_divergent_causal_sets(
        variants, ncsl, n_priv, n_pops, rng
    )
    assert used
    n_shared = ncsl - n_priv
    assert len(union) == n_shared + n_pops * n_priv        # union size
    assert len(per_pop) == n_pops
    for s in per_pop:
        assert len(s) == ncsl                              # each pop has ncsl
    # Shared variants: in every pop. Private: in exactly one.
    shared = set(per_pop[0])
    for s in per_pop[1:]:
        shared &= set(s)
    assert len(shared) == n_shared
    # At least two pops have distinct causal sets.
    assert len({frozenset(s.tolist()) for s in per_pop}) > 1


def test_ancestry_divergent_falls_back_when_ncsl_too_small() -> None:
    """ncsl <= n_private leaves no shared variant → fall back to a shared set."""
    variants = pd.DataFrame({"snp_id": [f"rs{i}" for i in range(50)]})
    rng = np.random.default_rng(2)
    union, per_pop, used = select_ancestry_divergent_causal_sets(
        variants, ncsl=1, n_private_per_pop=1, n_pops=6, rng=rng
    )
    assert not used
    assert len(union) == 1
    # every pop shares the identical set in the fallback
    assert all(np.array_equal(per_pop[0], s) for s in per_pop)


def test_assemble_divergent_beta_sparsity() -> None:
    """Shared β rows nonzero in all pops; private rows nonzero in exactly one."""
    variants = pd.DataFrame({"snp_id": [f"rs{i}" for i in range(50)]})
    rng = np.random.default_rng(3)
    ncsl, n_priv, n_pops = 3, 1, 6
    union, per_pop, used = select_ancestry_divergent_causal_sets(
        variants, ncsl, n_priv, n_pops, rng
    )
    beta = assemble_divergent_beta(union, per_pop, n_pops, rg=0.7, rng=rng)
    assert beta.shape == (len(union), n_pops)
    nz = (beta != 0.0).sum(axis=1)
    n_shared = ncsl - n_priv
    assert (nz == n_pops).sum() == n_shared     # shared rows
    assert (nz == 1).sum() == n_pops * n_priv   # private rows
    assert np.all((nz == n_pops) | (nz == 1))   # no other pattern


def test_divergent_beta_isolates_populations() -> None:
    """A pop-private causal variant must not affect other pops' genetic values."""
    variants = pd.DataFrame({"snp_id": [f"rs{i}" for i in range(50)]})
    rng = np.random.default_rng(4)
    ncsl, n_priv, n_pops = 2, 1, 3
    union, per_pop, used = select_ancestry_divergent_causal_sets(
        variants, ncsl, n_priv, n_pops, rng
    )
    beta = assemble_divergent_beta(union, per_pop, n_pops, rg=1.0, rng=rng)
    # One individual per pop, dosage 1 on every union variant.
    dosage = np.ones((n_pops, len(union)), dtype=np.float64)
    pop_index = np.arange(n_pops, dtype=np.int64)
    g = _compute_genetic_value(dosage, beta, pop_index)
    # Each individual's g equals the sum of β over *their* causal set only,
    # which equals the sum of their β column (zeros elsewhere).
    np.testing.assert_allclose(g, beta.sum(axis=0))


def test_select_causal_too_few_raises() -> None:
    variants = pd.DataFrame({"snp_id": [f"rs{i}" for i in range(2)]})
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="ncsl="):
        select_causal_variants(variants, ncsl=5, rng=rng)


def test_find_ancestry_specific_candidates() -> None:
    variants = pd.DataFrame({"snp_id": [f"rs{i}" for i in range(5)]})
    sites = ["anl", "covenant", "mbzuai"]
    dominant = {"anl": "EUR", "covenant": "AFR", "mbzuai": "MID"}
    maf = {
        "anl":     np.array([0.10, 0.001, 0.20, 0.30, 0.005]),
        "covenant":np.array([0.001, 0.20, 0.001, 0.30, 0.40]),
        "mbzuai":  np.array([0.50, 0.40, 0.001, 0.30, 0.50]),
    }
    # variant 0: common in anl/mbzuai, rare in covenant → yes
    # variant 1: common in covenant/mbzuai, rare in anl → yes
    # variant 2: common in anl, rare in covenant/mbzuai → yes
    # variant 3: common in all three (none rare) → no
    # variant 4: common in covenant/mbzuai, rare in anl → yes
    mask = _find_ancestry_specific_candidates(
        variants, sites, dominant, maf,
        common_threshold=0.05, rare_threshold=0.01,
    )
    expected = np.array([True, True, True, False, True])
    np.testing.assert_array_equal(mask, expected)


class _ToySiteState:
    """Mock SiteState for end-to-end phenotype test (avoids .bed I/O)."""

    def __init__(self, n_samples: int, pop_index: np.ndarray, dosage: np.ndarray, name="toy"):
        self.name = name
        self.n_samples = n_samples
        self.pop_index = pop_index
        self._dosage = dosage  # (n_samples, n_variants)
        # Stub fam & bed_prefix so simulate_phenotype_for_site can use it.
        self.fam = pd.DataFrame({"FID": [f"f{i}" for i in range(n_samples)],
                                 "IID": [f"i{i}" for i in range(n_samples)]})


def test_simulate_phenotype_empirical_h2_close_to_target(monkeypatch):
    """End-to-end: simulated phenotypes should hit the target h² within tolerance.

    Uses a stub for read_bed_variants so we never touch the disk.
    """
    rng_data = np.random.default_rng(42)
    n_samples = 5000
    n_pops = 3
    ncsl = 2
    # Synthetic dosage: bernoulli-like centered around varying AFs per population.
    pop_index = rng_data.integers(0, n_pops, size=n_samples)
    af = np.array([0.3, 0.5])  # per-variant
    dosage = rng_data.binomial(2, af[None, :], size=(n_samples, ncsl)).astype(np.float32)

    state = _ToySiteState(n_samples=n_samples, pop_index=pop_index, dosage=dosage)
    state.bed_prefix = "/tmp/toy"  # never read

    # Monkeypatch the reader.
    import appfl_bio_suite.experiments.fine_mapping.fedfm.phenotype_sim as mod
    monkeypatch.setattr(mod, "read_bed_variants", lambda bp, vidx, n: dosage[:, :len(vidx)])

    beta = draw_effect_sizes(ncsl, n_pops, rg=0.7, rng=np.random.default_rng(0))
    noise_rng = np.random.default_rng(7)
    y, emp_h2, _ = simulate_phenotype_for_site(
        state, np.array([0, 1]), beta, target_h2=0.05, rng=noise_rng
    )
    assert y.shape == (n_samples,)
    assert abs(emp_h2 - 0.05) / 0.05 < 0.15  # within 15%


def test_determinism_same_seed_same_output(monkeypatch):
    """Same inputs → same outputs at the simulate-for-site level."""
    rng_data = np.random.default_rng(0)
    n_samples = 1000
    pop_index = rng_data.integers(0, 3, size=n_samples)
    dosage = rng_data.binomial(2, 0.3, size=(n_samples, 1)).astype(np.float32)
    state = _ToySiteState(n_samples=n_samples, pop_index=pop_index, dosage=dosage)
    state.bed_prefix = "/tmp/toy"
    import appfl_bio_suite.experiments.fine_mapping.fedfm.phenotype_sim as mod
    monkeypatch.setattr(mod, "read_bed_variants", lambda bp, vidx, n: dosage[:, :len(vidx)])
    beta = draw_effect_sizes(1, 3, rg=1.0, rng=np.random.default_rng(0))
    y1, _, _ = simulate_phenotype_for_site(state, np.array([0]), beta, 0.05, np.random.default_rng(99))
    y2, _, _ = simulate_phenotype_for_site(state, np.array([0]), beta, 0.05, np.random.default_rng(99))
    np.testing.assert_allclose(y1, y2)
