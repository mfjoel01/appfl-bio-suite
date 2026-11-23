"""Federated covariate adjustment (FWL) reproduces regressing covariates out centrally."""

import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.covariates import CovariateGeno, fwl_adjust_locus
from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import (
    GenoAggregate,
    PhenoAggregate,
    pool_geno,
)
from appfl_bio_suite.experiments.fine_mapping.multisusie import MultiSuSiEOptions, prepare_inputs

SITES = ("s1", "s2")
POPS = ("P1", "P2")
M, K = 6, 3


def _cohort(seed=7):
    rng = np.random.default_rng(seed)
    rows = {}
    start = 0
    for site in SITES:
        for pop in POPS:
            n = int(rng.integers(20, 30))
            rows[(site, pop)] = np.arange(start, start + n)
            start += n
    n_total = start
    X = rng.integers(0, 3, size=(n_total, M)).astype(np.float64)
    P2 = np.concatenate([rows[(s, "P2")] for s in SITES])
    X[P2, 5] = 0.0  # monomorphic in P2: must stay ineligible there after adjustment
    Z = rng.normal(size=(n_total, K))
    Z[:, 0] += X[:, 1]  # a covariate that tracks a genotype, as genotype PCs do
    y = 0.4 * X[:, 2] + Z @ np.array([0.5, -0.3, 0.2]) + rng.normal(size=n_total)
    return rows, X, Z, y


def _blocks(rows, X, Z, y, inst="L1_arch_rep0"):
    variants = pd.DataFrame(
        {
            "chrom": ["1"] * M,
            "snp_id": [f"1:{i + 100}:A:G" for i in range(M)],
            "bp": np.arange(100, 100 + M),
            "a1": ["G"] * M,
            "a2": ["A"] * M,
        }
    )
    geno, pheno, cov_geno, cov_pheno = [], [], {}, {}
    for (site, pop), r in rows.items():
        Xb, Zb, yb = X[r], Z[r], y[r]
        geno.append(
            GenoAggregate(
                site=site,
                pop=pop,
                locus_id="L1",
                variants=variants,
                G=Xb.T @ Xb,
                u=Xb.sum(0),
                n=len(r),
            )
        )
        pheno.append(
            PhenoAggregate(
                site=site,
                pop=pop,
                instance=inst,
                snp_ids=variants.snp_id.to_numpy(),
                c=Xb.T @ yb,
                q=float(yb.sum()),
                w=float(yb @ yb),
                n=len(r),
            )
        )
        cov_geno[(site, pop)] = CovariateGeno(zx=Zb.T @ Xb, zz=Zb.T @ Zb, zu=Zb.sum(0))
        cov_pheno[(site, pop, inst)] = Zb.T @ yb
    return geno, {inst: pheno}, cov_geno, cov_pheno


def test_fwl_moments_equal_moments_of_centrally_residualized_data():
    rows, X, Z, y = _cohort()
    geno, pheno, cov_geno, cov_pheno = _blocks(rows, X, Z, y)
    fwl_adjust_locus(geno, pheno, cov_geno, cov_pheno)

    Z1 = np.column_stack([np.ones(len(y)), Z])
    Q, _ = np.linalg.qr(Z1)
    Xp = X - Q @ (Q.T @ X)
    yp = y - Q @ (Q.T @ y)
    by_key = {(p.site, p.pop): p for p in next(iter(pheno.values()))}
    for g in geno:
        r = rows[(g.site, g.pop)]
        np.testing.assert_allclose(g.G, Xp[r].T @ Xp[r], atol=1e-9)
        np.testing.assert_allclose(g.u, Xp[r].sum(0), atol=1e-9)
        np.testing.assert_array_equal(g.u_raw, X[r].sum(0))
        p = by_key[(g.site, g.pop)]
        np.testing.assert_allclose(p.c, Xp[r].T @ yp[r], atol=1e-9)
        assert p.q == pytest.approx(float(yp[r].sum()), abs=1e-9)
        assert p.w == pytest.approx(float(yp[r] @ yp[r]), abs=1e-9)


def test_single_site_run_is_that_site_adjusting_alone():
    rows, X, Z, y = _cohort()
    geno, pheno, cov_geno, cov_pheno = _blocks(rows, X, Z, y)
    keep = [g for g in geno if g.site == "s1"]
    inst = next(iter(pheno))
    fwl_adjust_locus(keep, {inst: [p for p in pheno[inst] if p.site == "s1"]}, cov_geno, cov_pheno)
    r_all = np.concatenate([rows[("s1", p)] for p in POPS])
    Z1 = np.column_stack([np.ones(len(r_all)), Z[r_all]])
    B = np.linalg.lstsq(Z1, X[r_all], rcond=None)[0]
    Xp = X[r_all] - Z1 @ B
    offset = 0
    for g in keep:
        n = len(rows[("s1", g.pop)])
        np.testing.assert_allclose(
            g.G, Xp[offset : offset + n].T @ Xp[offset : offset + n], atol=1e-9
        )
        offset += n


def test_allele_frequencies_and_eligibility_use_raw_dosages():
    rows, X, Z, y = _cohort()
    geno, pheno, cov_geno, cov_pheno = _blocks(rows, X, Z, y)
    fwl_adjust_locus(geno, pheno, cov_geno, cov_pheno)
    pooled = pool_geno([g for g in geno if g.pop == "P2"])
    assert pooled.allele_totals[5] == 0.0  # raw: monomorphic
    assert pooled.u[5] != 0.0  # adjusted dosage moments are not allele counts
    inputs = prepare_inputs(geno, next(iter(pheno.values())), MultiSuSiEOptions(maf=0.01))
    p2 = inputs.populations.index("P2")
    assert not inputs.eligible[p2][5] and inputs.eligible[inputs.populations.index("P1")][5]


def test_mismatched_variant_lists_are_refused():
    rows, X, Z, y = _cohort()
    geno, pheno, cov_geno, cov_pheno = _blocks(rows, X, Z, y)
    geno[0].snp_ids = geno[0].snp_ids[::-1].copy()
    with pytest.raises(ValueError, match="one variant list"):
        fwl_adjust_locus(geno, pheno, cov_geno, cov_pheno)


def test_federation_config_sets_both_sides_of_covariate_adjustment():
    from appfl_bio_suite.core.config import Experiment, load_federation
    from appfl_bio_suite.core.experiments import repo_root
    from appfl_bio_suite.core.launch import build_client_configs, build_server_config

    fed = load_federation(repo_root() / "federation.yaml.example")
    exp = fed.experiment("fine-mapping")
    for client in build_client_configs(fed, "fine-mapping"):
        assert "covariates" not in client["train_configs"]
    server = build_server_config(fed, "fine-mapping")["server_configs"]["aggregator_kwargs"]
    assert "covariate_adjustment" not in server

    exp.covariate_adjustment = "fwl"
    for client in build_client_configs(fed, "fine-mapping"):
        assert client["train_configs"]["covariates"] == "covariates.tsv"
    server = build_server_config(fed, "fine-mapping")["server_configs"]["aggregator_kwargs"]
    assert server["covariate_adjustment"] == "fwl"
    assert "susiex_finite_sample_adjust" not in server
    exp.susiex_finite_sample_adjust = True
    server = build_server_config(fed, "fine-mapping")["server_configs"]["aggregator_kwargs"]
    assert server["susiex_finite_sample_adjust"] is True
    with pytest.raises(ValueError):
        Experiment.model_validate({**exp.model_dump(), "covariate_adjustment": "residualize"})


def test_finite_sample_adjust_gives_susiex_the_marginal_correlation():
    from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import (
        finite_sample_adjust,
    )

    rng = np.random.default_rng(11)
    n = 300
    x = rng.integers(0, 3, size=(n, 4)).astype(np.float64)
    y = 1.5 * x[:, 0] + rng.normal(size=n)  # a large effect, as an eQTL lead can be
    xc, yc = x - x.mean(0), y - y.mean()
    sxx, sxy, syy = (xc * xc).sum(0), xc.T @ yc, yc @ yc
    beta = sxy / sxx
    se = np.sqrt((syy - beta * sxy) / (n - 2) / sxx)
    ss = pd.DataFrame({"beta": beta, "se": se, "stat": beta / se, "p": 0.5})
    adjusted = finite_sample_adjust(ss, n)
    r = sxy / np.sqrt(sxx * syy)
    # SuSiEx standardizes as (BETA / SE) / sqrt(n); that must now be the correlation.
    np.testing.assert_allclose(adjusted["beta"] / adjusted["se"] / np.sqrt(n), r, rtol=1e-12)
    assert abs(beta[0] / se[0] / np.sqrt(n)) > abs(r[0]) + 0.1  # unadjusted overstates it
    np.testing.assert_array_equal(adjusted["se"], ss["se"])
    np.testing.assert_array_equal(adjusted["p"], ss["p"])


def test_inverse_normal_transform_and_expression_pcs():
    from appfl_bio_suite.experiments.fine_mapping.real_data.mage import (
        expression_pcs,
        inverse_normal_transform,
    )

    values = np.array([[5.0, 1.0, 3.0, 3.0, 9.0]])
    out = inverse_normal_transform(values, axis=1)
    assert out[0, 2] == out[0, 3]  # ties share a rank
    assert np.argsort(out[0]).tolist()[0] == 1 and out[0].argmax() == 4
    untied = inverse_normal_transform(np.array([[4.0, 1.0, 3.0, 2.0, 9.0]]), axis=1)
    np.testing.assert_allclose(untied.mean(), 0.0, atol=1e-12)  # symmetric quantiles
    rng = np.random.default_rng(2)
    scores = expression_pcs(rng.normal(size=(40, 12)), 3)
    assert scores.shape == (40, 3)
    np.testing.assert_allclose(
        scores.T @ scores - np.diag(np.diag(scores.T @ scores)), 0, atol=1e-8
    )
