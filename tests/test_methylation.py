"""Data integrity and reproducibility checks for the methylation scaffold."""

import numpy as np
import pytest

from appfl_bio_suite.core.experiments import get_spec, implemented_experiments
from appfl_bio_suite.experiments.methylation.dataset import (
    MethylationData,
    synthetic_data,
)


def test_fixture_reproducible_and_round_trips(tmp_path):
    first = synthetic_data(seed=9)
    second = synthetic_data(seed=9)
    np.testing.assert_array_equal(first.X, second.X)
    np.testing.assert_array_equal(first.y, second.y)
    assert not np.array_equal(first.X, synthetic_data(seed=10).X)
    first.X[0, 0] = np.nan
    first.site_id = np.full(len(first.y), "S1")
    path = tmp_path / "cohort.npz"
    first.save(path)
    loaded = MethylationData.load(path)
    for key in ("X", "y", "site_id", "cpg_ids"):
        np.testing.assert_array_equal(getattr(first, key), getattr(loaded, key))
    with pytest.raises(FileExistsError):
        first.save(path)


@pytest.mark.parametrize("value", [-0.1, 1.1, np.inf, -np.inf])
def test_invalid_betas_rejected(value):
    data = synthetic_data()
    data.X[0, 0] = value
    with pytest.raises(ValueError, match="beta values"):
        data.validate()


def test_duplicate_probes_and_misaligned_labels_rejected():
    data = synthetic_data()
    data.cpg_ids[1] = data.cpg_ids[0]
    with pytest.raises(ValueError, match="unique"):
        data.validate()
    data = synthetic_data()
    data.y = data.y[:-1]
    with pytest.raises(ValueError, match="dimensions"):
        data.validate()


def test_fixture_has_no_site_assignment_and_balanced_classes():
    data = synthetic_data(samples=120, classes=3)
    assert data.site_id is None
    np.testing.assert_array_equal(np.unique(data.y, return_counts=True)[1], [40, 40, 40])


def test_serial_experiment_is_registered():
    assert get_spec("methylation").package_path.is_dir()
    assert "methylation" in implemented_experiments()
    assert get_spec("methylation").serial_runner


def test_patient_holdouts_panel_and_site_integrity():
    from appfl_bio_suite.experiments.methylation.preprocessing import prepare

    data = synthetic_data(samples=120, cpgs=32)
    # A second measurement from one patient must stay with that patient.
    data.X = np.concatenate([data.X, data.X[:1]])
    data.y = np.concatenate([data.y, data.y[:1]])
    data.patient_ids = np.concatenate([data.patient_ids, data.patient_ids[:1]])
    data.sample_ids = np.append(data.sample_ids, "repeat_measurement")
    data.X[:, 0] = np.nan
    prepared = prepare(data, panel_size=16, seed=7)
    ids = [
        set(data.patient_ids[idx]) for idx in (prepared.train, prepared.validation, prepared.test)
    ]
    assert all(not ids[i] & ids[j] for i in range(3) for j in range(i))
    np.testing.assert_array_equal(np.sort(np.concatenate(prepared.sites)), prepared.train)
    assert all(len(site) > 0 for site in prepared.sites)
    assert 0 not in prepared.panel
    assert len(prepared.sites[-1]) < len(prepared.train) * 0.15
    # Held-out feature values cannot change the selected training panel.
    data.X[prepared.test] = 1
    data.X[prepared.validation] = 0
    other = prepare(data, panel_size=16, seed=7)
    np.testing.assert_array_equal(prepared.panel, other.panel)


def test_sparse_calls_preserve_missingness_and_requested_coverage():
    from appfl_bio_suite.experiments.methylation.preprocessing import simulate_calls

    X = np.full((5, 20), 0.5)
    X[:, :3] = np.nan
    calls = simulate_calls(X, np.random.default_rng(4), observed=7)
    assert set(np.unique(calls)) <= {-1, 0, 1}
    assert not calls[:, :3].any()
    np.testing.assert_array_equal(np.count_nonzero(calls, axis=1), [7] * 5)
    np.testing.assert_array_equal(calls, simulate_calls(X, np.random.default_rng(4), observed=7))
    X[:] = np.nan
    assert not simulate_calls(X, np.random.default_rng(4), observed=7).any()


def test_undefined_confident_accuracy_and_calibration():
    from appfl_bio_suite.experiments.methylation.evaluation import (
        fit_temperature,
        metrics,
    )

    logits = np.zeros((6, 3), dtype=np.float32)
    labels = np.array([0, 1, 2, 0, 1, 2])
    result = metrics(logits, labels)
    assert result["confident_call_rate"] == 0
    assert result["confident_accuracy"] is None
    assert result["balanced_accuracy"] == pytest.approx(1 / 3)
    assert result["ece"] == pytest.approx(0, abs=1e-7)
    assert 0.05 <= fit_temperature(logits, labels) <= 20


def test_appfl_aggregates_by_sample_count():
    import torch
    from appfl.algorithm.aggregator import FedAvgAggregator
    from omegaconf import OmegaConf

    model = torch.nn.Linear(1, 1, bias=False)
    aggregator = FedAvgAggregator(model, OmegaConf.create({"client_weights_mode": "sample_size"}))
    aggregator.set_client_sample_size("small", 1)
    aggregator.set_client_sample_size("large", 3)
    result = aggregator.aggregate(
        {"small": {"weight": torch.tensor([[2.0]])}, "large": {"weight": torch.tensor([[6.0]])}}
    )
    assert result["weight"].item() == pytest.approx(5)


def test_serial_pipeline_outputs_and_reproducibility(tmp_path):
    import json

    from appfl_bio_suite.experiments.methylation.pipeline import RunConfig, run_experiment

    config = RunConfig(
        samples=90,
        cpgs=16,
        panel_size=16,
        min_class=10,
        rounds=1,
        hidden=8,
        bottleneck=4,
        coverage_counts=(4, 16),
    )
    first, second = tmp_path / "first", tmp_path / "second"
    run_experiment(first, config)
    run_experiment(second, config)
    assert (first / "metrics.json").read_bytes() == (second / "metrics.json").read_bytes()
    result = json.loads((first / "metrics.json").read_text())
    assert len(result) == 12  # Six models at two coverages.
    assert all(row["n_test"] == 18 for row in result)
    assert json.loads((first / "status.json").read_text())["complete"]
    for name in (
        "checkpoint.pt",
        "site_comparison.png",
        "coverage_curve.png",
        "confusion_matrix.png",
        "manifest.json",
        "REPORT.md",
    ):
        assert (first / name).stat().st_size > 0
    with pytest.raises(FileExistsError):
        run_experiment(first, config)


def test_shared_cli_serial_dry_run_and_remote_rejection(tmp_path):
    from click.testing import CliRunner

    from appfl_bio_suite.cli import main

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "run",
            "methylation",
            "--driver",
            "serial",
            "--dry-run",
            "--out-dir",
            str(tmp_path / "resolved"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "resolved" / "config.json").is_file()
    result = runner.invoke(main, ["run", "methylation"])
    assert result.exit_code != 0
    assert "requires --driver serial" in result.output


def test_simulate_then_run_cli_preserves_synthetic_provenance(tmp_path):
    import json

    from click.testing import CliRunner

    from appfl_bio_suite.cli import main

    config = tmp_path / "tiny.yaml"
    config.write_text(
        "samples: 90\ncpgs: 16\npanel_size: 16\nmin_class: 10\n"
        "rounds: 1\nhidden: 8\nbottleneck: 4\ncoverage_counts: [4, 16]\n"
        "signal_fraction: 0.5\nmin_site_class_patients: 1\n"
    )
    runner = CliRunner()
    cohort = tmp_path / "cohort"
    result = runner.invoke(
        main, ["simulate", "methylation", "--scenario", str(config), "--out", str(cohort)]
    )
    assert result.exit_code == 0, result.output
    result = runner.invoke(
        main, ["simulate", "methylation", "--verify", str(cohort / "manifest.json")]
    )
    assert result.exit_code == 0, result.output
    out = tmp_path / "run"
    result = runner.invoke(
        main,
        [
            "run",
            "methylation",
            "--driver",
            "serial",
            "--config",
            str(config),
            "--data-root",
            str(cohort),
            "--out-dir",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    manifest = json.loads((out / "manifest.json").read_text())
    generated = MethylationData.load(cohort / "cohort.npz")
    expected = synthetic_data(samples=90, cpgs=16, signal_fraction=0.5)
    np.testing.assert_array_equal(generated.X, expected.X)
    assert manifest["input_kind"] == "synthetic"
    assert manifest["input_sha256"]
    assert manifest["aggregator_sha256"]
    assert (
        "Synthetic" in (out / "REPORT.md").read_text()
        or "synthetic" in (out / "REPORT.md").read_text()
    )


def test_partner_bundle_refuses_serial_experiment_before_loading_federation():
    from click.testing import CliRunner

    from appfl_bio_suite.cli import main

    result = CliRunner().invoke(main, ["partner-bundle", "methylation", "--site", "S1"])
    assert result.exit_code != 0
    assert "partner bundles are unavailable" in result.output


def test_class_support_control_preserves_site_sizes_and_holdouts():
    from appfl_bio_suite.experiments.methylation.preprocessing import prepare

    data = synthetic_data(samples=768, cpgs=64, classes=6)
    skewed = prepare(data, panel_size=32, seed=42)
    controlled = prepare(data, panel_size=32, seed=42, min_site_class_patients=2)
    np.testing.assert_array_equal(skewed.test, controlled.test)
    np.testing.assert_array_equal(skewed.panel, controlled.panel)
    assert [len(s) for s in skewed.sites] == [len(s) for s in controlled.sites]
    for site in controlled.sites:
        assert np.bincount(controlled.y[site], minlength=6).min() >= 2
    np.testing.assert_array_equal(np.sort(np.concatenate(controlled.sites)), controlled.train)
    with pytest.raises(ValueError, match="capacity"):
        prepare(data, panel_size=32, seed=42, min_site_class_patients=100)


def test_reduced_signal_leaves_most_features_without_large_class_mean_differences():
    data = synthetic_data(samples=600, cpgs=256, classes=3, signal_fraction=0.05)
    class_means = np.stack([data.X[data.y == c].mean(axis=0) for c in np.unique(data.y)])
    assert np.count_nonzero(np.ptp(class_means, axis=0) > 0.25) <= 20
    dense = synthetic_data(samples=600, cpgs=256, classes=3)
    dense_means = np.stack([dense.X[dense.y == c].mean(axis=0) for c in np.unique(dense.y)])
    assert np.count_nonzero(np.ptp(dense_means, axis=0) > 0.25) >= 230
    with pytest.raises(ValueError, match="one CpG per class"):
        synthetic_data(cpgs=16, signal_fraction=0.001)


def test_uneven_synthetic_controls_are_reproducible():
    options = dict(
        samples=1000,
        cpgs=128,
        classes=3,
        seed=13,
        signal_fraction=0.25,
        class_weights=(5, 3, 2),
        profile_noise=0.05,
        missing_fraction=0.03,
    )
    first = synthetic_data(**options)
    second = synthetic_data(**options)
    np.testing.assert_array_equal(first.X, second.X)
    np.testing.assert_array_equal(np.unique(first.y, return_counts=True)[1], [500, 300, 200])
    assert 0.025 < np.isnan(first.X).mean() < 0.035
    clean = synthetic_data(**{**options, "profile_noise": 0})
    np.testing.assert_array_equal(first.y, clean.y)
    assert not np.array_equal(first.X, clean.X, equal_nan=True)
    assert first.site_id is None


@pytest.mark.parametrize(
    "options",
    [
        {"class_weights": (1, 2)},
        {"class_weights": (1, 0, 2)},
        {"class_weights": (1, np.inf, 2)},
        {"profile_noise": -0.1},
        {"profile_noise": 1},
        {"missing_fraction": np.nan},
        {"missing_fraction": 1},
    ],
)
def test_invalid_uneven_controls(options):
    from appfl_bio_suite.experiments.methylation.pipeline import RunConfig

    with pytest.raises(ValueError):
        synthetic_data(**options)
    with pytest.raises(ValueError):
        RunConfig(**options).validate()


def test_validation_only_never_uses_test_features(tmp_path):
    import json

    from appfl_bio_suite.experiments.methylation.pipeline import RunConfig, run_experiment
    from appfl_bio_suite.experiments.methylation.preprocessing import prepare

    data = synthetic_data(samples=120, cpgs=32)
    config = RunConfig(
        samples=120,
        cpgs=32,
        panel_size=32,
        rounds=1,
        hidden=16,
        bottleneck=8,
        coverage_counts=(16, 32),
        min_class=2,
        track_validation=True,
    )
    p = prepare(data, panel_size=32, seed=config.seed, min_class=2)
    original = tmp_path / "original.npz"
    modified = tmp_path / "modified.npz"
    data.save(original)
    data.X[p.test] = np.nan
    data.save(modified)
    for name, source in [("first", original), ("second", modified)]:
        run_experiment(tmp_path / name, config, source, validation_only=True)
    first = json.loads((tmp_path / "first" / "metrics.json").read_text())
    assert first == json.loads((tmp_path / "second" / "metrics.json").read_text())
    for name in ("validation_history.json", "history.json"):
        assert (tmp_path / "first" / name).read_bytes() == (tmp_path / "second" / name).read_bytes()
    assert all(r["evaluation_split"] == "validation" and "n_test" not in r for r in first)
    assert all(r["n_validation"] == len(p.validation) for r in first)
    assert not (tmp_path / "first" / "confusion_matrix.png").exists()


@pytest.mark.parametrize(
    ("yaml_text", "message"),
    [
        ("rounds: 1\nround: 2\n", "unknown setting(s)"),
        ("alpha: high\n", "alpha must be a finite number"),
        ("profile_noise: lots\n", "profile_noise must lie in [0, 1)"),
        ("coverage_counts: 16\n", "coverage_counts must be a list of positive integers"),
        ("- 1\n- 2\n", "must be a mapping of RunConfig fields"),
    ],
)
def test_malformed_config_is_a_clean_cli_error(tmp_path, yaml_text, message):
    from click.testing import CliRunner

    from appfl_bio_suite.cli import main

    config = tmp_path / "bad.yaml"
    config.write_text(yaml_text)
    runner = CliRunner()
    run = runner.invoke(
        main,
        [
            "run",
            "methylation",
            "--driver",
            "serial",
            "--config",
            str(config),
            "--out-dir",
            str(tmp_path / "run"),
        ],
    )
    assert run.exit_code == 1 and message in run.output and "Traceback" not in run.output
    simulate = runner.invoke(
        main,
        ["simulate", "methylation", "--scenario", str(config), "--out", str(tmp_path)],
    )
    assert simulate.exit_code == 1 and message in simulate.output
    assert not (tmp_path / "run").exists()


def test_unknown_scenario_names_the_packaged_ones():
    from appfl_bio_suite.experiments.methylation.pipeline import load_config

    packaged = "ci-tiny, poc, stress, tcga-laml, uneven, uneven-curves"
    with pytest.raises(ValueError, match=f"scenarios: {packaged}$"):
        load_config("no-such-scenario")


@pytest.fixture(scope="module")
def completed_run(tmp_path_factory):
    """A tiny tested run, a curves run of it, and a curves run that trained differently."""
    from dataclasses import replace

    from appfl_bio_suite.experiments.methylation.pipeline import RunConfig, run_experiment

    root = tmp_path_factory.mktemp("figures")
    config = RunConfig(
        samples=90,
        cpgs=64,
        panel_size=64,
        min_class=10,
        rounds=1,
        hidden=8,
        bottleneck=4,
        coverage_counts=(16, 64),
        min_site_class_patients=1,
        missing_fraction=0.02,
    )
    run_experiment(root / "run", config)
    curves = replace(
        config, rounds=2, coverage_counts=(16, 32, 64), track_validation=True, validation_only=True
    )
    run_experiment(root / "curves", curves)
    run_experiment(root / "other-curves", replace(curves, learning_rate=0.01))
    return root


def test_analysis_figures_rebuild_byte_identically(completed_run):
    from appfl_bio_suite.experiments.methylation.figures import build

    run, curves = completed_run / "run", completed_run / "curves"
    first = build(run, completed_run / "first", curves)
    second = build(run, completed_run / "second", curves)
    assert first["figures"] == [
        "01_who_gains",
        "02_training_rounds",
        "03_sparse_input",
        "04_class_recall",
        "05_site_composition",
        "A1_data_overview",
        "A2_cpg_methylation",
    ]
    assert {name.rsplit(".", 1)[1] for name in first["files"]} == {"png", "svg"}
    assert first["files"] == second["files"]
    assert "<dc:date>" not in (completed_run / "first" / "01_who_gains.svg").read_text()
    # Without a curves run, the two validation figures are simply absent.
    alone = build(run, completed_run / "alone")
    assert "02_training_rounds" not in alone["figures"] and len(alone["figures"]) == 5
    low, high = first["who_gains"]["federated"]["interval"]
    assert low <= first["who_gains"]["federated"]["balanced_accuracy"] <= high


def test_balanced_accuracy_interval_brackets_the_estimate():
    from appfl_bio_suite.experiments.methylation.figures import interval

    assert interval([[5, 0], [0, 5]]) == (1.0, 1.0)
    low, high = interval([[30, 10], [20, 20]])
    assert low < (30 / 40 + 20 / 40) / 2 < high
    assert interval([[30, 10], [20, 20]]) == (low, high)


def test_analysis_figures_refuse_runs_that_do_not_match(completed_run, tmp_path):
    import json
    import shutil

    from appfl_bio_suite.experiments.methylation.figures import build

    run = completed_run / "run"
    with pytest.raises(ValueError, match="trains differently"):
        build(run, tmp_path / "out", completed_run / "other-curves")
    with pytest.raises(ValueError, match="must be validation-only"):
        build(run, tmp_path / "out", run)
    tampered = tmp_path / "tampered"
    shutil.copytree(run, tampered)
    (tampered / "metrics.json").write_text("[]\n")
    with pytest.raises(ValueError, match="does not match its run manifest"):
        build(tampered, tmp_path / "out")
    manifest = json.loads((tampered / "manifest.json").read_text())
    manifest["evaluation_split"] = "validation"
    (tampered / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="validation patients only"):
        build(tampered, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_experiments_listing_installs_what_the_serial_run_needs():
    """A serial experiment has no partner extras; its install line must carry its own."""
    from click.testing import CliRunner

    from appfl_bio_suite.cli import main

    result = CliRunner().invoke(main, ["experiments"])
    assert result.exit_code == 0, result.output
    blocks = result.output.strip().split("\n\n")
    entry = {block.split("\n", 1)[0]: block for block in blocks}["methylation"]
    assert "appfl-bio-suite[methylation]" in entry
    assert "single-machine simulation" in entry


def test_validation_tracking_is_passive_and_ends_at_the_final_scores(tmp_path):
    """Scoring every round must not perturb training, and its last round is the result."""
    import json
    from dataclasses import replace

    from appfl_bio_suite.experiments.methylation.pipeline import RunConfig, run_experiment

    # validation_only set in the config, not the keyword: the config alone must suffice.
    base = RunConfig(
        samples=90,
        cpgs=16,
        panel_size=16,
        min_class=10,
        rounds=2,
        hidden=8,
        bottleneck=4,
        coverage_counts=(4, 16),
        validation_only=True,
    )
    run_experiment(tmp_path / "plain", base)
    run_experiment(tmp_path / "tracked", replace(base, track_validation=True))
    for name in ("history.json", "metrics.json"):
        plain_file, tracked_file = (tmp_path / run / name for run in ("plain", "tracked"))
        assert plain_file.read_bytes() == tracked_file.read_bytes()
    assert not (tmp_path / "plain" / "validation_history.json").exists()
    tracked = json.loads((tmp_path / "tracked" / "validation_history.json").read_text())
    assert [entry["round"] for entry in tracked] == [1, 2]
    final = json.loads((tmp_path / "plain" / "metrics.json").read_text())
    final = {(row["model"], row["requested_cpgs"]): row for row in final}
    assert len(tracked[-1]["scores"]) == len(final) == 12
    for score in tracked[-1]["scores"]:
        row = final[(score["model"], score["requested_cpgs"])]
        assert score["balanced_accuracy"] == row["balanced_accuracy"]
        assert score["confusion_matrix"] == row["confusion_matrix"]
    manifest = json.loads((tmp_path / "plain" / "manifest.json").read_text())
    assert manifest["evaluation_split"] == "validation"


def test_curves_config_trains_exactly_the_frozen_uneven_scenario():
    from dataclasses import asdict

    from appfl_bio_suite.experiments.methylation.pipeline import load_config

    frozen, curves = asdict(load_config("uneven")), asdict(load_config("uneven-curves"))
    differing = {key for key in frozen if frozen[key] != curves[key]}
    assert differing == {"rounds", "coverage_counts", "track_validation", "validation_only"}
    assert curves["validation_only"] and curves["track_validation"]
    assert curves["rounds"] > frozen["rounds"]
    assert max(curves["coverage_counts"]) == max(frozen["coverage_counts"])
