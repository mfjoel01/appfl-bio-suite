"""Recovery must preserve the model and refuse unconverged or changed-prior parity."""

import importlib.util
import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from appfl_bio_suite.experiments.fine_mapping.multisusie import MultiSuSiEOptions


def recovery():
    path = Path(__file__).parents[1] / "scripts/fine-mapping/recover_multisusie.py"
    spec = importlib.util.spec_from_file_location("recovery", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_retry_changes_only_iteration_ceiling_and_preserves_failed_attempt(tmp_path, monkeypatch):
    module = recovery()
    seen = []

    def fit(g, p, dest, options):
        dest.mkdir()
        seen.append(asdict(options))
        record = {"converged": options.max_iter > 1000, "options": asdict(options)}
        (dest / "fit.json").write_text(json.dumps(record))
        return record

    monkeypatch.setattr(module, "run_with_trace", fit)
    record = module.recover_fit([], [], tmp_path / "fit", MultiSuSiEOptions())
    assert record["converged"]
    assert len(seen) == 2
    assert seen[0] | {"max_iter": 10000} == seen[1]
    assert not json.loads((Path(record["recovery"]["previous_attempt"]) / "fit.json").read_text())[
        "converged"
    ]


def test_parity_budget_exception_still_requires_convergence_model_and_pips(tmp_path):
    module = recovery()
    a, b = tmp_path / "a", tmp_path / "b"
    for path, budget in [(a, 1000), (b, 10000)]:
        path.mkdir()
        pd.DataFrame({"snp_id": ["rs1"], "pip": [0.9], "eligible_EUR": [True]}).to_csv(
            path / "variants.tsv", sep="\t", index=False
        )
        pd.DataFrame({"component": [1], "snp_id": ["rs1"]}).to_csv(
            path / "credible_sets.tsv", sep="\t", index=False
        )
        (path / "fit.json").write_text(
            json.dumps(
                dict(
                    converged=True,
                    populations=["EUR"],
                    sample_sizes=[100],
                    options=asdict(MultiSuSiEOptions(max_iter=budget)),
                )
            )
        )
    assert module.recovery_parity(a, b)["pass"]
    record = json.loads((b / "fit.json").read_text())
    record["options"]["rho"] = 0.5
    (b / "fit.json").write_text(json.dumps(record))
    assert not module.recovery_parity(a, b)["pass"]
    record["options"]["rho"] = 0.75
    record["converged"] = False
    (b / "fit.json").write_text(json.dumps(record))
    assert not module.recovery_parity(a, b)["pass"]


def test_manager_retry_uses_smaller_batches_without_failed_dependency(tmp_path, monkeypatch):
    from types import SimpleNamespace

    path = Path(__file__).parents[1] / "scripts/fine-mapping/watch_multisusie_recovery.py"
    spec = importlib.util.spec_from_file_location("manager", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rec = tmp_path / "recovery"
    rec.mkdir()
    (rec / "job.pbs").write_text("worker")
    (rec / "recover_multisusie.py").write_text("source")
    cmd = [
        "qsub",
        "-W",
        "depend=afterok:failed",
        "-J",
        "0-30",
        "-o",
        str(rec / "logs"),
        "-v",
        f"FM_RECOVERY={rec},FM_WORKERS=4",
        str(rec / "job.pbs"),
    ]
    (rec / "submission.benchmark.json").write_text(json.dumps({"command": cmd}))
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(
            returncode=0, stdout="123[].pbs\n", stderr="", check_returncode=lambda: None
        )

    monkeypatch.setattr(module.subprocess, "run", run)
    tasks = [{"name": str(i)} for i in range(13)]
    assert module.submit(tmp_path, rec, "benchmark", 1, tasks) == "123[].pbs"
    actual = calls[0]
    assert "-W" not in actual
    assert actual[actual.index("-J") + 1] == "0-2"
    assert "FM_WORKERS=2" in actual[actual.index("-v") + 1]
    batches = json.loads((rec / "retry-benchmark-1/benchmark.batches.json").read_text())
    assert list(map(len, batches)) == [6, 6, 1]
    assert [t for batch in batches for t in batch] == tasks


def test_failed_calibration_screen_evaluates_baselines_without_forcing_winner(
    tmp_path, monkeypatch
):
    path = Path(__file__).parents[1] / "scripts/fine-mapping/select_multisusie_recovery.py"
    spec = importlib.util.spec_from_file_location("selection_recovery", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    p = module.pipeline
    monkeypatch.setattr(p, "verify", lambda root: None)

    def no_candidate(root):
        (root / "development_results.tsv").write_text("completed development evidence")
        raise ValueError("No candidate passes the prespecified calibration screen")

    monkeypatch.setattr(p, "select", no_candidate)

    def config(root, phase, name, source):
        path = root / "configs" / phase / f"{name}.yaml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}")
        return path

    monkeypatch.setattr(p, "freeze_config", config)

    def tasks(root, phase, arms):
        p.write(root / f"{phase}_tasks.json", arms)

    monkeypatch.setattr(p, "make_tasks", tasks)
    p.write(tmp_path / "design.json", {"source_composition": str(tmp_path), "sampling_seeds": [1]})
    module.select(tmp_path)
    selection = p.read(tmp_path / "selection.json")
    assert selection["name"] is None and selection["eligible"] is False
    assert not selection["evaluation_used_for_selection"]
    assert {t["candidate"] for t in p.read(tmp_path / "evaluation_tasks.json")} == {
        "federation",
        "anl",
        "covenant",
        "mbzuai",
    }
    assert (tmp_path / "selection.lock.json").exists()


def test_diagnostic_reconstructs_original_batched_moments(tmp_path):
    import pytest
    from test_fine_mapping_federated import _build_package

    pytest.importorskip("MultiSuSiE")
    module = recovery()
    pipeline = module.pipeline
    config = _build_package(tmp_path / "package")
    root = tmp_path / "run"
    root.mkdir()
    frozen = pipeline.freeze_config(root, "benchmark", "federation", config)
    tasks = pipeline.make_tasks(
        root, "benchmark", [dict(name="federation", config=str(frozen), central_check=False)]
    )
    pipeline.work_item(root, tasks[0])
    out = root / "diagnostics-test"
    instance = "L0000_ncsl1_h2-0.005_rg1_rep0"
    module.diagnostic(root, out, [("benchmark", "federation", "L0000", instance)])
    new = pipeline.read(out / "diagnostics/benchmark" / instance / "fit.json")
    old = pipeline.read(root / "benchmark/federation/L0000" / instance / "multisusie/fit.json")
    assert new["converged"]
    assert new["moments_sha256"] == old["moments_sha256"]
    assert (out / "diagnostics.accepted.json").exists()


def test_manager_waits_for_primary_and_supplementary_batches():
    path = Path(__file__).parents[1] / "scripts/fine-mapping/watch_multisusie_recovery.py"
    spec = importlib.util.spec_from_file_location("manager", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    done = {"job_state": "F"}
    for active in ["R", "Q", "H", "B"]:
        assert not module.phase_finished(done, [{"job_state": active}])
        assert not module.phase_finished({"job_state": active}, [done])
    assert module.phase_finished(done, [done, done])
    assert module.phase_finished(done, [])


def test_credible_set_difference_is_reported_without_claiming_strict_parity(tmp_path):
    module = recovery()
    a, b = tmp_path / "a", tmp_path / "b"
    for path, snp in [(a, "rs1"), (b, "rs2")]:
        path.mkdir()
        pd.DataFrame({"snp_id": ["rs1", "rs2"], "pip": [0.5, 0.5]}).to_csv(
            path / "variants.tsv", sep="\t", index=False
        )
        pd.DataFrame({"component": [1], "snp_id": [snp]}).to_csv(
            path / "credible_sets.tsv", sep="\t", index=False
        )
        (path / "fit.json").write_text(
            json.dumps(
                dict(
                    converged=True,
                    populations=["EUR"],
                    sample_sizes=[100],
                    options=asdict(MultiSuSiEOptions()),
                )
            )
        )
    result = module.recovery_parity(a, b)
    assert result["pass"]
    assert not result["strict_pass"]
    assert not result["credible_sets_exact"]
    assert result["credible_sets_federated"] != result["credible_sets_centralized"]
    pd.DataFrame({"snp_id": ["rs1", "rs2"], "pip": [0.7, 0.3]}).to_csv(
        b / "variants.tsv", sep="\t", index=False
    )
    assert not module.recovery_parity(a, b)["pass"]


def test_manager_only_retries_resolved_missing_task_acceptance(tmp_path):
    path = Path(__file__).parents[1] / "scripts/fine-mapping/watch_multisusie_recovery.py"
    spec = importlib.util.spec_from_file_location("manager", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    missing = tmp_path / "development" / "arm" / "locus" / "accepted.json"
    log = tmp_path / "select.log"
    log.write_text(f"FileNotFoundError: [Errno 2] No such file or directory: '{missing}'\n")
    assert not module.resolved_prerequisite_race(log, tmp_path)
    missing.parent.mkdir(parents=True)
    missing.write_text("{}")
    assert module.resolved_prerequisite_race(log, tmp_path)
    log.write_text("ValueError: failed calibration\n")
    assert not module.resolved_prerequisite_race(log, tmp_path)
