"""Both configured backends consume the existing APPFL wire protocol."""

import json
import logging

import numpy as np
import pandas as pd
import pytest
import test_fine_mapping_shipped_parity as site_fixture
from test_fine_mapping_shipped_parity import (
    INSTANCES,
    SITES,
    _shipped,
    _vendored,
)

from appfl_bio_suite.core.config import Experiment, load_federation
from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.launch import build_client_configs, build_server_config
from appfl_bio_suite.core.preflight import PreflightReport, _check_binaries
from appfl_bio_suite.experiments.fine_mapping.aggregator import FineMappingAggregator
from appfl_bio_suite.experiments.fine_mapping.multisusie import run_multisusie


@pytest.fixture(scope="module")
def package(tmp_path_factory):
    return site_fixture.package.__wrapped__(tmp_path_factory)


def aggregator(out, **kwargs):
    return FineMappingAggregator(
        aggregator_configs=dict(
            backend="multisusie", output_dir=str(out), make_figures=False, **kwargs
        ),
        logger=logging.getLogger("multisusie-appfl"),
    )


def test_backend_selection_reaches_server_only():
    fed = load_federation(repo_root() / "federation.yaml.example")
    exp = fed.experiment("fine-mapping")
    before = build_client_configs(fed, "fine-mapping")
    exp.fine_mapping_backend = "multisusie"
    exp.multisusie_rho = 0.5
    cfg = build_server_config(fed, "fine-mapping")["server_configs"]["aggregator_kwargs"]
    assert cfg["backend"] == "multisusie" and cfg["multisusie_rho"] == 0.5
    assert build_client_configs(fed, "fine-mapping") == before
    with pytest.raises(ValueError):
        Experiment.model_validate({**exp.model_dump(), "fine_mapping_backend": "typo"})
    with pytest.raises(ValueError):
        Experiment.model_validate({**exp.model_dump(), "multisusie_rho": 1.0})


def test_multisusie_preflight_does_not_require_susiex(monkeypatch):
    pytest.importorskip("MultiSuSiE")
    fed = load_federation(repo_root() / "federation.yaml.example")
    fed.experiment("fine-mapping").fine_mapping_backend = "multisusie"
    monkeypatch.setattr("shutil.which", lambda _: pytest.fail("Unexpected binary lookup"))
    report = PreflightReport()
    _check_binaries(report, "fine-mapping", fed)
    assert report.ok


@pytest.mark.parametrize("workers,scored", [(1, True), (2, False)])
def test_wire_fit_matches_independent_standalone(package, tmp_path, monkeypatch, workers, scored):
    pytest.importorskip("MultiSuSiE")
    root, cfg, locus = package
    payloads = {site: _shipped(root, site)[0] for site in SITES}
    geno, pheno = [], []
    for site in SITES:
        g, p, _ = _vendored(cfg, site, locus)
        geno.extend(g.values())
        pheno.extend(p.values())
    truth = tmp_path / "truth.tsv"
    pd.DataFrame(
        [
            dict(
                locus_id="L0000",
                architecture_id=i.split("L0000_")[1].rsplit("_rep", 1)[0],
                replicate=0,
                causal_snp_ids="rs0",
            )
            for i in INSTANCES
        ]
    ).to_csv(truth, sep="\t", index=False)
    # Neither binary resolution nor the SuSiEx LD path is used, even with no binaries.
    monkeypatch.setattr(
        "appfl_bio_suite.experiments.fine_mapping.aggregator._resolve_binary",
        lambda *a: pytest.fail("Unexpected SuSiEx dependency"),
    )
    a = aggregator(
        tmp_path / "appfl",
        n_workers=workers,
        causal_manifest=str(truth) if scored else None,
        ga4gh={"drs": {"hostname": "results.example.org"}},
    )
    state = a.aggregate(payloads)
    assert state["fine_mapping_ready"].item() == 1
    results = pd.read_csv(a.data_dir / "fed_fm_results.tsv", sep="\t")
    assert len(results) == len(INSTANCES) and results.converged.all()
    assert (results.backend == "multisusie").all()
    assert (
        results.any_causal_captured.notna().all()
        if scored
        else results.any_causal_captured.isna().all()
    )
    assert (
        json.loads((a.data_dir / "inference_backend.json").read_text())["backend"] == "multisusie"
    )
    provenance = json.loads((a.data_dir / "ga4gh_provenance.json").read_text())
    assert provenance["run"]["backend"] == "multisusie"
    for instance in INSTANCES:
        assert f"instances/{instance}/fit.json" in provenance["outputs"]
        assert f"instances/{instance}/posterior.npz" in provenance["outputs"]
        direct = tmp_path / "direct" / instance
        run_multisusie(
            geno, [p for p in pheno if p.instance == instance], direct, a.multisusie_options
        )
        wire = a.data_dir / "instances" / instance
        assert (
            json.loads((wire / "fit.json").read_text())["moments_sha256"]
            == json.loads((direct / "fit.json").read_text())["moments_sha256"]
        )
        pd.testing.assert_frame_equal(
            pd.read_csv(wire / "variants.tsv", sep="\t"),
            pd.read_csv(direct / "variants.tsv", sep="\t"),
        )
        pd.testing.assert_frame_equal(
            pd.read_csv(wire / "credible_sets.tsv", sep="\t"),
            pd.read_csv(direct / "credible_sets.tsv", sep="\t"),
        )
        with np.load(wire / "posterior.npz") as w, np.load(direct / "posterior.npz") as d:
            for key in w.files:
                np.testing.assert_array_equal(w[key], d[key])


def test_nonconvergence_never_reports_ready(package, tmp_path, monkeypatch):
    pytest.importorskip("MultiSuSiE")
    from appfl_bio_suite.experiments.fine_mapping import multisusie

    root, _, _ = package
    payloads = {site: _shipped(root, site)[0] for site in SITES}

    def nonconverged(*args):
        record = run_multisusie(*args)
        record["converged"] = False
        return record

    monkeypatch.setattr(multisusie, "run_multisusie", nonconverged)
    a = aggregator(tmp_path / "out")
    with pytest.raises(RuntimeError, match="Nonconverged"):
        a.aggregate(payloads)
    assert a.global_state["fine_mapping_ready"].item() == 0
    assert list(a.data_dir.glob("instances/*/fit.json"))
    assert not (a.data_dir / "fed_fm_results.tsv").exists()


def test_mismatched_wire_sample_counts_fail(package, tmp_path):
    pytest.importorskip("MultiSuSiE")
    from appfl_bio_suite.experiments.fine_mapping.trainer import KEY_SCALARS, KEY_SEP

    root, _, _ = package
    payloads = {site: _shipped(root, site)[0] for site in SITES}
    payload = next(iter(payloads.values()))
    key = next(k for k in payload if k.startswith(KEY_SCALARS + KEY_SEP))
    payload[key] = payload[key].clone()
    payload[key][-1] += 1
    a = aggregator(tmp_path / "out")
    with pytest.raises(ValueError, match="cohort or variant order differs"):
        a.aggregate(payloads)
    assert a.global_state["fine_mapping_ready"].item() == 0


def test_missing_optional_backend_fails_preflight_and_construction(monkeypatch, tmp_path):
    import builtins

    original = builtins.__import__

    def no_backend(name, *args, **kwargs):
        if name.startswith("MultiSuSiE"):
            raise ImportError("backend intentionally unavailable")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_backend)
    fed = load_federation(repo_root() / "federation.yaml.example")
    fed.experiment("fine-mapping").fine_mapping_backend = "multisusie"
    report = PreflightReport()
    _check_binaries(report, "fine-mapping", fed)
    assert not report.ok and "multisusie" in report.render().lower()
    with pytest.raises(ImportError, match="coordinator extra"):
        aggregator(tmp_path / "out")


@pytest.mark.parametrize("backend", ["susiex", "multisusie"])
@pytest.mark.parametrize("driver", ["serial", "globus_compute"])
def test_cli_swaps_backend_without_changing_driver(package, tmp_path, monkeypatch, backend, driver):
    from click.testing import CliRunner

    from appfl_bio_suite import cli

    root, _, _ = package
    fed = load_federation(repo_root() / "federation.yaml.example")
    monkeypatch.setattr(cli, "_load", lambda *args, **kwargs: fed)
    captured = {}

    def launch(**kwargs):
        captured.update(kwargs)
        return 0

    monkeypatch.setattr("appfl_bio_suite.core.launch.launch", launch)
    args = [
        "run",
        "fine-mapping",
        "--driver",
        driver,
        "--backend",
        backend,
        "--results-dir",
        str(tmp_path / backend),
    ]
    if driver == "serial":
        args += ["--data-root", str(root)]
    result = CliRunner().invoke(cli.main, args)
    assert result.exit_code == 0, result.output
    assert captured["driver"] == driver
    config = build_server_config(captured["federation"], "fine-mapping")
    settings = config["server_configs"]["aggregator_kwargs"]
    assert settings["backend"] == backend
    assert settings["output_dir"] == str(tmp_path / backend)


def test_both_backends_have_same_result_and_artifact_contract(package, tmp_path):
    pytest.importorskip("MultiSuSiE")
    import shutil

    binary = shutil.which("SuSiEx")
    if not binary:
        pytest.skip("SuSiEx unavailable")
    root, _, _ = package
    payloads = {site: _shipped(root, site)[0] for site in SITES}
    results = {}
    for backend in ["susiex", "multisusie"]:
        output = tmp_path / backend
        a = FineMappingAggregator(
            aggregator_configs=dict(
                backend=backend,
                output_dir=str(output),
                make_figures=False,
                susiex_binary=binary,
                n_workers=1,
            ),
            logger=logging.getLogger("swap-contract"),
        )
        assert a.aggregate(payloads)["fine_mapping_ready"].item() == 1
        results[backend] = pd.read_csv(output / "data/fed_fm_results.tsv", sep="\t")
        for instance in INSTANCES:
            files = output / "data/instances" / instance
            assert {"snp_id", "pip"} <= set(pd.read_csv(files / "variants.tsv", sep="\t"))
            assert list(pd.read_csv(files / "credible_sets.tsv", sep="\t")) == [
                "component",
                "snp_id",
                "component_pip",
                "pip",
                "coverage",
            ]
            assert json.loads((files / "fit.json").read_text())["converged"]
        assert results[backend].any_causal_captured.isna().all()
        assert results[backend].n_false_pip95.isna().all()
        # A second run cannot silently overwrite the first estimator's evidence.
        with pytest.raises(FileExistsError):
            a.aggregate(payloads)
        assert a.global_state["fine_mapping_ready"].item() == 0
    assert list(results["susiex"]) == list(results["multisusie"])
    assert results["susiex"].locus_id.equals(results["multisusie"].locus_id)
