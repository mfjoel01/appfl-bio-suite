"""Real backend parity, ordered LD alignment, and posterior interpretation."""

import json
import logging
import shutil
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping import multisusie_backend as backend
from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import (
    pooled_ld,
    pooled_sumstats,
    write_ld_panel,
)


def write_inputs(root, order=(0, 1, 2), flip=False):
    root.mkdir(parents=True, exist_ok=True)
    variants = pd.DataFrame(
        {
            "chrom": ["1"] * 3,
            "snp_id": ["v0", "v1", "v2"],
            "bp": [1, 2, 3],
            "a1": ["A"] * 3,
            "a2": ["G"] * 3,
        }
    )
    variants = variants.iloc[list(order)].reset_index(drop=True)
    if flip:
        variants.loc[0, "a1"], variants.loc[0, "a2"] = "G", "A"
    full = np.array([[1, 0.1, 0.2], [0.1, 1, 0.3], [0.2, 0.3, 1]])
    prefix = root / "ld"
    write_ld_panel(prefix, variants, full[np.ix_(order, order)], np.full(len(order), 0.3), 100)
    ss = variants.rename(columns={"chrom": "chr", "snp_id": "snp", "a1": "A1", "a2": "A2"})
    ss = ss[["chr", "snp", "bp", "A1", "A2"]].assign(beta=0.1, se=0.05, stat=2.0, p=0.04)
    path = root / "pop.sumstats"
    ss.to_csv(path, sep="\t", index=False)
    return path, prefix


def test_intersection_reorders_both_ld_axes_and_reports_exclusions(tmp_path):
    a, ar = write_inputs(tmp_path / "a")
    b, br = write_inputs(tmp_path / "b", order=(2, 0))
    ids, frames, matrices, records = backend.aligned_inputs([a, b], [ar, br], [100, 100])
    assert ids == ["v0", "v2"]
    assert frames[1].snp.tolist() == ids
    np.testing.assert_array_equal(matrices[0], matrices[1])
    assert records[0]["excluded_variants"] == ["v1"]


def test_cross_ancestry_allele_mismatch_rejected(tmp_path):
    a, ar = write_inputs(tmp_path / "a")
    b, br = write_inputs(tmp_path / "b", flip=True)
    with pytest.raises(ValueError, match="alleles"):
        backend.aligned_inputs([a, b], [ar, br], [100, 100])


@pytest.mark.parametrize("converged", [True, False])
def test_native_pips_are_not_reconstructed_from_filtered_cs(tmp_path, monkeypatch, converged):
    ss, ld = write_inputs(tmp_path)
    fit = SimpleNamespace(
        converged=converged,
        pip=np.array([0.8, 0.3, 0.1]),
        alpha=np.array([[0.8, 0.1, 0.1], [0.1, 0.8, 0.1]]),
        sets=[[np.array([0]), np.array([1])], np.array([1.0, 1.0]), [0.95, 0.95], [True, False]],
    )
    monkeypatch.setattr(backend, "require_multisusie", lambda: lambda **kwargs: fit)
    monkeypatch.setattr(backend.importlib.metadata, "version", lambda _: "test")
    row = backend.run_multisusie(
        [ss], [ld], [100], [1.0], tmp_path, ["v1", "absent"], 0.95, 1.0, {}
    )
    assert row["pip_definition"] == "multisusie_full_model"
    assert row["n_eligible_causal"] == row["n_excluded_causal"] == 1
    if converged:
        assert row["causal_pip_mean"] == 0.15
        assert row["n_credible_sets"] == 1
        assert row["top_pip"] == 0.8
    else:
        assert row["fit_status"] == "nonconverged"
        assert row["n_credible_sets"] == 0
        assert np.isnan(row["causal_pip_mean"])


def test_multisusie_real_fit_matches_upstream_and_recovers_signal(tmp_path):
    fit_function = pytest.importorskip("MultiSuSiE").multisusie_rss
    rng = np.random.default_rng(24)
    paths, prefixes, sizes, variances, frames, matrices = [], [], [], [], [], []
    variants = pd.DataFrame(
        {
            "chrom": ["1"] * 20,
            "snp_id": [f"v{i}" for i in range(20)],
            "bp": np.arange(1, 21),
            "a1": ["A"] * 20,
            "a2": ["G"] * 20,
        }
    )
    for i, n in enumerate([250, 320]):
        x = rng.binomial(2, 0.3, size=(n, 20)).astype(float)
        y = 2 * x[:, 4] + rng.normal(size=n)
        g, u = x.T @ x, x.sum(axis=0)
        r, _ = pooled_ld(g, u, n)
        frame = pooled_sumstats(variants, np.diag(g), u, x.T @ y, y.sum(), y @ y, n)
        path = tmp_path / f"pop{i}.sumstats"
        frame.to_csv(path, sep="\t", index=False, float_format="%.17g")
        prefix = write_ld_panel(tmp_path / f"ld{i}", variants, r, u / (2 * n), n)
        paths.append(path)
        prefixes.append(prefix)
        sizes.append(n)
        variances.append(float(y.var(ddof=1)))
        frames.append(frame)
        matrices.append(r.astype(np.float32).astype(float))
    options = {"n_signals": 2, "max_iter": 200, "tol": 1e-6}
    row = backend.run_multisusie(
        paths, prefixes, sizes, variances, tmp_path, ["v4"], 0.95, 1e-5, options
    )
    expected = fit_function(
        b_list=[f.beta.to_numpy() for f in frames],
        s_list=[f.se.to_numpy() for f in frames],
        R_list=matrices,
        varY_list=variances,
        population_sizes=sizes,
        rho=0.75,
        L=2,
        max_iter=200,
        tol=1e-6,
        coverage=0.95,
        min_abs_corr=0.5,
        float_type=np.float64,
        low_memory_mode=False,
        variant_ids=variants.snp_id.tolist(),
    )
    actual = pd.read_csv(tmp_path / "cs.snp", sep="\t").PIP.to_numpy()
    np.testing.assert_allclose(actual, expected.pip, atol=1e-10, rtol=1e-8)
    assert row["converged"] and row["any_causal_captured"]
    assert row["causal_pip_max"] > 0.99


@pytest.mark.parametrize("method", ["susiex", "multisusie"])
def test_secure_backend_matches_plaintext(tmp_path, tmp_path_factory, monkeypatch, method):
    # This fixture includes allele flips and two ancestries with unequal site sizes.
    from test_fine_mapping_shipped_parity import POPS
    from test_fine_mapping_shipped_parity import package as make_package

    from appfl_bio_suite.experiments.fine_mapping.aggregator import FineMappingAggregator
    from appfl_bio_suite.experiments.fine_mapping.dataset import get_dataset
    from appfl_bio_suite.experiments.fine_mapping.trainer import SiteFineMappingTrainer

    susiex = shutil.which("SuSiEx")
    if method == "susiex" and not susiex:
        pytest.skip("SuSiEx binary required")
    if method == "multisusie":
        pytest.importorskip("MultiSuSiE")
    root, _, _ = make_package.__wrapped__(tmp_path_factory)
    sites = ["covenant", "mbzuai"]
    config = {"run_id": "test-" + method, "participants": sites}
    payloads = {False: {}, True: {}}
    for site in sites:
        dataset, _ = get_dataset(root / site / "data", site)
        for secure in (False, True):
            settings = {"trainer_output_dirname": str(tmp_path / site), "pops": POPS}
            if secure:
                path = tmp_path / f"{site}.json"
                path.write_text(
                    json.dumps(
                        {
                            **config,
                            "client_id": site,
                            "keys": {peer: "a1" * 32 for peer in sites if peer != site},
                        }
                    )
                )
                monkeypatch.setenv("APPFL_FINEMAP_SECAGG_KEYS", str(path))
                settings["secure_aggregation"] = config
            trainer = SiteFineMappingTrainer(
                train_dataset=dataset,
                train_configs=settings,
                client_id=site,
                logger=logging.getLogger("backend-parity"),
            )
            trainer.train()
            payloads[secure][site] = trainer.get_parameters()
    tables = []
    for secure in (False, True):
        output = tmp_path / str(secure)
        cfg = {
            "backend": method,
            "output_dir": str(output),
            "make_figures": False,
            "n_workers": 1,
            "pval_thresh": 1.0,
            "maf": 0.001,
            "n_signals": 2,
            "susiex_binary": susiex,
            "plink_binary": shutil.which("plink") or "plink",
        }
        if secure:
            cfg["secure_aggregation"] = config
        agg = FineMappingAggregator(
            aggregator_configs=cfg, logger=logging.getLogger("backend-parity")
        )
        agg.aggregate(payloads[secure])
        table = pd.read_csv(output / "data/fed_fm_results.tsv", sep="\t")
        assert len(table) == 2 and table.error.fillna("").eq("").all()
        assert table.converged.all()
        assert list((output / "artifacts").glob("*/cs.snp.gz"))
        if secure:
            assert (output / "data/secure_aggregation.json").is_file()
        tables.append(table)
    from appfl_bio_suite.experiments.fine_mapping.inference import inclusion_probabilities

    for plain_path in (tmp_path / "False/artifacts").glob("*/cs.snp.gz"):
        secure_path = tmp_path / "True/artifacts" / plain_path.parent.name / plain_path.name
        a = pd.read_csv(plain_path, sep="\t").set_index("SNP").sort_index()
        b = pd.read_csv(secure_path, sep="\t").set_index("SNP").sort_index()
        assert a.index.equals(b.index)
        np.testing.assert_allclose(
            inclusion_probabilities(a), inclusion_probabilities(b), atol=1e-5, rtol=1e-4
        )
    np.testing.assert_array_equal(tables[0].n_credible_sets, tables[1].n_credible_sets)
    np.testing.assert_allclose(tables[0].top_pip, tables[1].top_pip, atol=1e-5, rtol=1e-4)


def test_reporting_reads_explicit_multisusie_pips():
    from appfl_bio_suite.experiments.fine_mapping.inference import inclusion_probabilities

    frame = pd.DataFrame(
        {"SNP": ["a", "b"], "PIP": [0.9, 0.2], "PIP_DEFINITION": ["multisusie_full_model"] * 2}
    )
    np.testing.assert_array_equal(inclusion_probabilities(frame), [0.9, 0.2])


def test_preflight_multisusie_does_not_require_susiex(monkeypatch):
    from appfl_bio_suite.core.config import load_federation
    from appfl_bio_suite.core.experiments import repo_root
    from appfl_bio_suite.core.preflight import Level, PreflightReport, _check_binaries

    fed = load_federation(repo_root() / "federation.yaml.example")
    fed.experiment("fine-mapping").fine_mapping_backend = "multisusie"
    monkeypatch.setattr(shutil, "which", lambda _: None)
    monkeypatch.setattr(backend, "require_multisusie", lambda: object())
    report = PreflightReport()
    _check_binaries(report, "fine-mapping", fed)
    assert len(report.checks) == 1
    assert report.checks[0].level == Level.OK

    def unavailable():
        raise ImportError("install MultiSuSiE")

    monkeypatch.setattr(backend, "require_multisusie", unavailable)
    report = PreflightReport()
    _check_binaries(report, "fine-mapping", fed)
    assert report.checks[0].level == Level.FAIL


@pytest.mark.parametrize(
    "options",
    [
        {"multisusie_rho": 1.0},
        {"multisusie_rho": -0.1},
        {"multisusie_rho": float("nan")},
        {"multisusie_min_abs_corr": 1.1},
        {"max_iter": 0},
        {"n_signals": 1.5},
        {"tol": -1.0},
    ],
)
def test_invalid_model_options_rejected(options):
    with pytest.raises(ValueError):
        backend.validate_options(options)


def test_disjoint_ancestry_panels_fail_instead_of_fitting_misaligned_snps(tmp_path):
    a, ar = write_inputs(tmp_path / "a", order=(0,))
    b, br = write_inputs(tmp_path / "b", order=(2,))
    with pytest.raises(ValueError, match="no common"):
        backend.aligned_inputs([a, b], [ar, br], [100, 100])
