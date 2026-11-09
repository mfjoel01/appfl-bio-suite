"""Real-backend parity and input-contract tests for the MultiSuSiE pilot."""

import importlib.util
import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import (
    GenoAggregate,
    PhenoAggregate,
)
from appfl_bio_suite.experiments.fine_mapping.multisusie import (
    MomentInputs,
    MultiSuSiEOptions,
    fit_multisusie,
    prepare_inputs,
    run_multisusie,
)


def blocks(x, y, pop="EUR", split=True):
    variants = pd.DataFrame(
        {
            "chrom": "1",
            "snp_id": [f"rs{i}" for i in range(x.shape[1])],
            "bp": np.arange(x.shape[1]),
            "a1": "A",
            "a2": "G",
        }
    )
    geno, pheno = [], []
    parts = np.array_split(np.arange(len(x)), 3) if split else [np.arange(len(x))]
    for i, rows in enumerate(parts):
        a, b = x[rows], y[rows]
        site = f"site{i}"
        geno.append(GenoAggregate(site, pop, "L0", variants.copy(), a.T @ a, a.sum(0), len(a)))
        pheno.append(
            PhenoAggregate(
                site,
                pop,
                "L0_trait_rep0",
                variants.snp_id.to_numpy(),
                a.T @ b,
                b.sum(),
                b @ b,
                len(a),
            )
        )
    return geno, pheno


def cohort():
    rng = np.random.default_rng(61)
    data = []
    for pop, n, freq, effect in [("EUR", 400, 0.3, 0.9), ("AFR", 300, 0.4, 0.7)]:
        x = rng.binomial(2, freq, (n, 16)).astype(float)
        y = 3 + effect * x[:, 4] + rng.normal(size=n)
        data.append((pop, x, y))
    return data


def test_null_signal_returns_no_credible_sets_without_backend_label_crash():
    pytest.importorskip("MultiSuSiE")
    inputs = MomentInputs(
        pd.DataFrame({"snp_id": ["a", "b", "c"]}),
        ["EUR"],
        [100],
        [np.eye(3) * 50],
        [np.zeros(3)],
        [100.0],
        [np.ones(3, dtype=bool)],
    )
    fit = fit_multisusie(inputs, MultiSuSiEOptions(n_signals=2))
    assert fit.converged
    assert not np.any(fit.sets[3])
    np.testing.assert_array_equal(fit.pip, np.zeros(3))


def test_federated_fit_matches_independently_centered_centralized_data():
    pytest.importorskip("MultiSuSiE")
    from MultiSuSiE.susiepy_ss import susie_multi_ss

    data = cohort()
    geno, pheno = [], []
    for pop, x, y in data:
        g, p = blocks(x, y, pop)
        geno.extend(g)
        pheno.extend(p)
    options = MultiSuSiEOptions(n_signals=3)
    inputs = prepare_inputs(geno, pheno, options)
    original = [a.copy() for a in inputs.xtx]
    fed = fit_multisusie(inputs, options)
    data.sort(key=lambda row: row[0])
    xs = [x - x.mean(0) for _, x, _ in data]
    ys = [y - y.mean() for _, _, y in data]
    for actual, x in zip(inputs.xtx, xs, strict=True):
        np.testing.assert_allclose(actual, x.T @ x, atol=1e-10)
    central = susie_multi_ss(
        XTX_list=[x.T @ x for x in xs],
        XTY_list=[x.T @ y for x, y in zip(xs, ys, strict=True)],
        YTY_list=[y @ y for y in ys],
        population_sizes=[len(y) for y in ys],
        rho=np.array([[1.0, 0.75], [0.75, 1.0]]),
        L=3,
        standardize=False,
        max_iter=1000,
        tol=1e-6,
        coverage=0.95,
        min_abs_corr=0,
        float_type=np.float64,
        low_memory_mode=True,
    )
    assert fed.converged and central.converged
    np.testing.assert_allclose(fed.pip, central.pip, atol=1e-8, rtol=1e-8)
    assert fed.pip[4] > 0.9
    assert [list(s) for s, keep in zip(fed.sets[0], fed.sets[3], strict=True) if keep] == [
        list(s) for s, keep in zip(central.sets[0], central.sets[3], strict=True) if keep
    ]
    for actual, before in zip(inputs.xtx, original, strict=True):
        np.testing.assert_array_equal(actual, before)


def test_union_censors_only_the_ancestry_where_variant_is_ineligible():
    data = cohort()
    geno, pheno = [], []
    for pop, x, y in data:
        x[:, 0] = 0 if pop == "EUR" else x[:, 0]
        x[:, 1] = 0 if pop == "AFR" else x[:, 1]
        g, p = blocks(x, y, pop)
        geno.extend(g)
        pheno.extend(p)
    inputs = prepare_inputs(geno, pheno, MultiSuSiEOptions())
    assert len(inputs.variants) == 16
    assert inputs.populations == ["AFR", "EUR"]
    assert inputs.eligible[0][0] and not inputs.eligible[1][0]
    assert inputs.eligible[1][1] and not inputs.eligible[0][1]
    assert not inputs.xtx[1][0].any() and inputs.xty[1][0] == 0
    assert not inputs.xtx[0][1].any() and inputs.xty[0][1] == 0
    if importlib.util.find_spec("MultiSuSiE"):
        fit = fit_multisusie(inputs, MultiSuSiEOptions(n_signals=3))
        assert fit.converged and np.isfinite(fit.pip).all()


def test_cohort_mismatch_and_duplicate_blocks_are_rejected():
    pop, x, y = cohort()[0]
    g, p = blocks(x, y, pop)
    with pytest.raises(ValueError, match="Duplicate"):
        prepare_inputs(g + g[:1], p + p[:1], MultiSuSiEOptions())
    with pytest.raises(ValueError, match="membership differs"):
        prepare_inputs(g, p[:-1], MultiSuSiEOptions())
    with pytest.raises(ValueError, match="cohort or variant order"):
        prepare_inputs(g, [replace(p[0], n=p[0].n - 1), *p[1:]], MultiSuSiEOptions())


def test_cross_ancestry_allele_mismatch_is_rejected():
    _, x, y = cohort()[0]
    g, p = blocks(x, y)
    other_g, other_p = blocks(x, y, "AFR")
    for block in other_g:
        block.variants.loc[4, ["a1", "a2"]] = ["G", "A"]
    with pytest.raises(ValueError, match="coding differ"):
        prepare_inputs(g + other_g, p + other_p, MultiSuSiEOptions())


@pytest.mark.parametrize("rho", [1.0, -0.1, float("nan")])
def test_invalid_or_singular_correlation_rejected(rho):
    with pytest.raises(ValueError, match="rho"):
        MultiSuSiEOptions(rho=rho)


def test_nonconvergence_is_retained_and_output_cannot_be_overwritten(tmp_path):
    pytest.importorskip("MultiSuSiE")
    _, x, y = cohort()[0]
    g, p = blocks(x, y)
    out = tmp_path / "fit"
    record = run_multisusie(g, p, out, MultiSuSiEOptions(n_signals=3, max_iter=1))
    assert record["fit_status"] == "nonconverged"
    assert json.loads((out / "fit.json").read_text())["backend_source_sha256"]
    assert (out / "posterior.npz").is_file()
    with pytest.raises(FileExistsError):
        run_multisusie(g, p, out)


def test_pilot_runs_existing_site_stage_without_external_binaries(tmp_path):
    pytest.importorskip("MultiSuSiE")
    from test_fine_mapping_federated import _build_package

    from appfl_bio_suite.experiments.fine_mapping.multisusie_experiment import run_pilot

    config = _build_package(tmp_path / "package")
    out = tmp_path / "pilot"
    rows = run_pilot(config, out, MultiSuSiEOptions(n_signals=3), limit=2)
    assert len(rows) == 2 and rows.converged.all()
    assert rows.any_causal_captured.all()
    assert (out / "results.tsv").is_file()


def test_three_method_comparison_and_source_genotype_parity(tmp_path):
    pytest.importorskip("MultiSuSiE")
    import os
    from pathlib import Path

    from test_fine_mapping_federated import _build_package

    from appfl_bio_suite.experiments.fine_mapping.multisusie_compare import run_comparison

    binary = os.environ.get("SUSIEX_TEST_BINARY")
    if not binary or not Path(binary).is_file():
        pytest.skip("Set SUSIEX_TEST_BINARY for the three-method comparison")
    config = _build_package(tmp_path / "package")
    panel = tmp_path / "panel.tsv"
    pd.DataFrame(
        [{"locus_id": "L0000", "architecture_id": "ncsl1_h2-0.005_rg1", "replicate": 0}]
    ).to_csv(panel, sep="\t", index=False)
    out = tmp_path / "compare"
    rows = run_comparison(config, panel, out, Path(binary), MultiSuSiEOptions(n_signals=3))
    assert set(rows.method) == {"susiex", "susiex_unfiltered", "multisusie"}
    assert rows.fit_status.str.startswith("converged_").all()
    assert rows.any_causal_captured.all()
    assert all(r["pass"] for r in json.loads((out / "parity.json").read_text()))
    execution = json.loads(next(out.glob("*/susiex_unfiltered/cs.execution.json")).read_text())
    command = execution["command"]
    assert command[command.index("--min_purity") + 1] == "0"
    assert command[command.index("--pval_thresh") + 1] == "1"


def test_panel_is_fixed_without_looking_at_inference_outcomes():
    from appfl_bio_suite.experiments.fine_mapping.multisusie_batch import (
        ARCHITECTURES,
        select_panel,
    )

    loci = pd.DataFrame(
        [
            {"locus_id": f"{s}{i}", "stratum": s, "divergence_score": i}
            for s in ["low", "medium", "high"]
            for i in range(3)
        ]
    )
    truth = pd.DataFrame(
        [
            {"locus_id": lid, "architecture_id": arch, "replicate": 0}
            for lid in loci.locus_id
            for arch in ARCHITECTURES
        ]
    )
    panel = select_panel(loci.sample(frac=1, random_state=1), truth)
    assert len(panel) == 12
    assert set(panel.locus_id) == {"low1", "medium1", "high1"}
    assert set(panel.architecture_id) == set(ARCHITECTURES)
    with pytest.raises(ValueError, match="absent"):
        select_panel(loci, truth[truth.locus_id != "low1"])


def test_parity_rejects_different_prior_even_when_pips_match(tmp_path):
    from appfl_bio_suite.experiments.fine_mapping.multisusie_compare import compare_fits

    for name, rho in [("a", 0.5), ("b", 0.75)]:
        directory = tmp_path / name
        directory.mkdir()
        pd.DataFrame({"snp_id": ["x"], "pip": [0.9], "eligible_EUR": [True]}).to_csv(
            directory / "variants.tsv", sep="\t", index=False
        )
        pd.DataFrame({"component": [1], "snp_id": ["x"]}).to_csv(
            directory / "credible_sets.tsv", sep="\t", index=False
        )
        (directory / "fit.json").write_text(
            json.dumps(
                {
                    "converged": True,
                    "populations": ["EUR"],
                    "sample_sizes": [100],
                    "options": {"rho": rho},
                }
            )
        )
    assert not compare_fits(tmp_path / "a", tmp_path / "b")["pass"]


def test_batch_finish_requires_complete_method_grid_and_parity(tmp_path):
    from appfl_bio_suite.experiments.fine_mapping.multisusie_batch import METHODS, finish

    (tmp_path / "tasks.json").write_text(json.dumps([{"arm": "full", "central_check": True}]))
    pd.DataFrame([{"locus_id": "L0", "architecture_id": "a", "replicate": 0}]).to_csv(
        tmp_path / "panel.tsv", sep="\t", index=False
    )
    directory = tmp_path / "full"
    directory.mkdir()
    (directory / "accepted.json").write_text(
        json.dumps(
            {
                "pass": True,
                "centralized_parity_checked": True,
            }
        )
    )
    frame = pd.DataFrame(
        [
            {
                "locus_id": "L0",
                "architecture_id": "a",
                "replicate": 0,
                "method": method,
                "fit_status": "converged_cs",
                "n_causal": 1,
                "n_causal_captured": 1,
                "n_credible_sets": 1,
                "n_cs_containing_causal": 1,
                "total_cs_snps": 2,
                "n_true_pip95": 0,
                "n_false_pip95": 0,
                "runtime_s": 1,
            }
            for method in sorted(METHODS)
        ]
    )
    frame.iloc[:-1].to_csv(directory / "results.tsv", sep="\t", index=False)
    with pytest.raises(RuntimeError, match="grid"):
        finish(tmp_path)
    frame.to_csv(directory / "results.tsv", sep="\t", index=False)
    (directory / "parity.json").write_text("[]")
    with pytest.raises(RuntimeError, match="parity"):
        finish(tmp_path)
    (directory / "parity.json").write_text(json.dumps([{"instance": "L0_a_rep0", "pass": True}]))
    finish(tmp_path)
    assert json.loads((tmp_path / "accepted.json").read_text())["n_method_fits"] == 3
    assert len(pd.read_csv(tmp_path / "summary.tsv", sep="\t")) == 3


def test_full_pipeline_checkpoint_and_merge_reject_tampered_results(tmp_path):
    pytest.importorskip("MultiSuSiE")
    from test_fine_mapping_federated import _build_package

    from appfl_bio_suite.experiments.fine_mapping.multisusie_pipeline import (
        collect,
        make_tasks,
        work_item,
    )

    config = _build_package(tmp_path / "package")
    root = tmp_path / "production"
    root.mkdir()
    tasks = make_tasks(
        root, "benchmark", [dict(name="federation", config=str(config), central_check=True)]
    )
    for task in tasks:
        work_item(root, task)
        work_item(root, task)
    rows = collect(root, "benchmark")
    assert len(rows) == 2 and rows.converged.all()
    assert rows.any_causal_captured.all()
    results = root / "benchmark/federation/L0000/results.tsv"
    results.write_text(results.read_text().replace("converged_cs", "execution_failure"))
    with pytest.raises(ValueError, match="Unaccepted"):
        collect(root, "benchmark")
    with pytest.raises(ValueError, match="Checkpoint changed"):
        work_item(root, tasks[0])


def test_pipeline_selects_on_development_and_builds_heldout_configs(tmp_path, monkeypatch):
    import yaml
    from test_fine_mapping_federated import _build_package

    from appfl_bio_suite.experiments.fine_mapping import multisusie_pipeline as pipeline

    root = tmp_path / "run"
    root.mkdir()
    dev_config = _build_package(tmp_path / "development-package")
    eval_config = _build_package(tmp_path / "evaluation-package")
    # Production prepare freezes absolute paths before relocating configurations.
    for config in [dev_config, eval_config]:
        cfg = pipeline.load_config(config)
        raw = yaml.safe_load(config.read_text())
        raw["paths"] = {k: str(cfg.resolved_path(k)) for k in raw["paths"]}
        config.write_text(yaml.safe_dump(raw))
    source = tmp_path / "source"
    (source / "evaluation/arms").mkdir(parents=True)
    (source / "evaluation/pipeline_config.yaml").write_text(eval_config.read_text())
    for name in ["federation", "anl", "covenant", "mbzuai"]:
        (source / f"evaluation/arms/{name}.yaml").write_text(eval_config.read_text())
    choices = [
        dict(name="solo_covenant", solo="covenant", site_counts={"covenant": {"AFR": 50}}),
        dict(name="mixed", solo=None, site_counts={"covenant": {"AFR": 30}, "anl": {"EUR": 20}}),
    ]
    pipeline.write(
        root / "design.json",
        dict(
            candidates=choices,
            source_composition=str(source),
            sampling_seeds=[1, 2, 3],
            calibration_screen=dict(min_cs_coverage=0.93, max_fdr_pip95=0.05),
        ),
    )
    pipeline.write(
        root / "development_tasks.json",
        [
            dict(
                name=f"mixed_seed{seed}",
                candidate="mixed",
                sampling_seed=seed,
                config=str(dev_config),
            )
            for seed in [1, 2, 3]
        ],
    )
    frame = pd.DataFrame(
        [
            dict(
                candidate=c,
                any_causal_captured=power,
                n_credible_sets=10,
                n_cs_containing_causal=10,
                n_true_pip95=10,
                n_false_pip95=0,
            )
            for c, power in [("solo_covenant", 0.7), ("mixed", 0.8)]
        ]
    )
    phases = []

    def development_only(directory, phase):
        phases.append(phase)
        frame.to_csv(directory / "development_results.tsv", sep="\t", index=False)
        return frame

    monkeypatch.setattr(pipeline, "collect", development_only)
    pipeline.select(root)
    assert phases == ["development"]
    assert pipeline.read(root / "selection.json")["name"] == "mixed"
    assert not pipeline.read(root / "selection.json")["evaluation_used_for_selection"]
    tasks = pipeline.read(root / "evaluation_tasks.json")
    assert len({t["name"] for t in tasks}) == 7
    raw = yaml.safe_load((root / "configs/evaluation/mixed_seed1.yaml").read_text())
    assert (
        raw["paths"]["ground_truth_dir"]
        == yaml.safe_load(eval_config.read_text())["paths"]["ground_truth_dir"]
    )
    assert (
        raw["paths"]["processed_dir"]
        == yaml.safe_load(dev_config.read_text())["paths"]["processed_dir"]
    )
    frame.loc[:, "n_cs_containing_causal"] = 1
    with pytest.raises(ValueError, match="No candidate passes"):
        pipeline.select(root)
