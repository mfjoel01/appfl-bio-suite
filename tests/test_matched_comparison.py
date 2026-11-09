import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping import matched_comparison as m


def test_common_purity_uses_same_eligible_ancestries():
    index = {"a": 0, "b": 1}
    grams = [np.array([[1.0, 0.2], [0.2, 1.0]]), np.array([[1.0, 0.8], [0.8, 1.0]])]
    assert m.common_purity({"a", "b"}, index, grams, [np.ones(2, bool), np.ones(2, bool)]) == 0.8
    assert (
        m.common_purity({"a", "b"}, index, grams, [np.ones(2, bool), np.array([True, False])])
        == 0.2
    )


def test_filter_does_not_change_native_pips_or_discovery_metrics():
    rows = m.scored_rows({}, "multisusie", {"a": 0.99, "b": 0.1}, [{"a", "b"}], {"a"}, [0.2])
    assert rows[0]["n_credible_sets"] == 1 and rows[1]["n_credible_sets"] == 0
    assert rows[0]["any_causal_captured"] and not rows[1]["any_causal_captured"]
    assert all(r["n_true_pip95"] == 1 and r["high_confidence_yield"] == 1 for r in rows)


def test_paired_bootstrap_identical_methods_has_zero_difference(tmp_path):
    rows = []
    for locus in ["A", "B"]:
        for method in m.METHODS:
            rows += m.scored_rows(
                dict(phase="benchmark", arm="federation", locus_id=locus),
                method,
                {"a": 0.99},
                [{"a"}],
                {"a"},
                [1.0],
            )
    m.summarize(pd.DataFrame(rows), tmp_path)
    d = pd.read_csv(tmp_path / "paired_differences.tsv", sep="\t")
    assert (d[["difference", "ci_low", "ci_high"]] == 0).all().all()


def test_complete_paired_task_and_checkpoint_reuse(tmp_path):
    pytest.importorskip("MultiSuSiE")
    from test_fine_mapping_federated import _build_package

    from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import _resolve_binary

    try:
        binary = _resolve_binary("SuSiEx", Path(__file__).parents[1])
    except FileNotFoundError:
        pytest.skip("SuSiEx executable unavailable")
    cfg = _build_package(tmp_path / "package")
    source = tmp_path / "source"
    source.mkdir()
    cfg = m.p.freeze_config(source, "benchmark", "federation", cfg)
    task = m.p.make_tasks(
        source, "benchmark", [dict(name="federation", config=str(cfg), central_check=False)]
    )[0]
    m.p.work_item(source, task)
    root = tmp_path / "matched"
    (root / "code/bin").mkdir(parents=True)
    shutil.copy2(binary, root / "code/bin/SuSiEx")
    m.p.write(root / "design.json", {"source": str(source)})
    m.run_task(root, task, workers=2)
    out = root / "benchmark/federation/L0000"
    a = m.verify_record(out / "accepted.json")
    assert a["n_rows"] == 8
    f = pd.read_csv(out / "results.tsv", sep="\t")
    assert f.converged.all()
    assert set(f.method) == set(m.METHODS) and set(f.policy) == set(m.POLICIES)
    digest = m.p.sha(out / "results.tsv")
    m.run_task(root, task, workers=2)
    assert m.p.sha(out / "results.tsv") == digest
