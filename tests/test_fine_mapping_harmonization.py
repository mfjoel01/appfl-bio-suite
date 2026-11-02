"""Harmonization figures must distinguish measured zeros from unavailable inputs."""

import warnings
from types import SimpleNamespace

import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.figures import federation as f


@pytest.fixture
def audit_inputs(tmp_path):
    for directory in ("reference", "truth", "alpha", "beta"):
        (tmp_path / directory).mkdir()
    original = "1 rs1 0 10 A G\n1 rs2 0 20 C T\n"
    (tmp_path / "reference/chr1.bim").write_text(original)
    for site in ("alpha", "beta"):
        (tmp_path / site / f"{site}_chr1.bim").write_text(original)
    truth = pd.DataFrame(
        {
            "locus_id": ["l1", "l2"],
            "architecture_id": "a",
            "replicate": 0,
            "causal_snp_ids": ["rs1", "rs2"],
        }
    )
    truth.to_csv(tmp_path / "truth/causal_manifest.tsv", sep="\t", index=False)
    cfg = SimpleNamespace(
        chromosome=1,
        sites={"alpha": {}, "beta": {}},
        site_dir=lambda s: tmp_path / s,
        resolved_path=lambda k: (
            tmp_path / {"hapnest_dir": "reference", "ground_truth_dir": "truth"}[k]
        ),
    )
    return cfg, truth[f.KEY_COLS]


def test_audit_reads_actual_alleles_and_only_selected_instances(audit_inputs, tmp_path):
    cfg, instances = audit_inputs
    assert f.audit_harmonization(cfg, instances) == ({"alpha": 0, "beta": 0}, 2, (0, 2))
    # Reverse row order too: comparison must align by variant ID.
    (tmp_path / "beta/beta_chr1.bim").write_text("1 rs2 0 20 C T\n1 rs1 0 10 G A\n")
    assert f.audit_harmonization(cfg, instances) == ({"alpha": 0, "beta": 1}, 2, (1, 2))
    assert f.audit_harmonization(cfg, instances.iloc[1:])[-1] == (0, 1)


@pytest.mark.parametrize(
    "contents",
    [
        "1 rs1 0 10 A G\n",  # missing site variant
        "1 rs1 0 10 A C\n1 rs2 0 20 C T\n",  # incompatible allele
        "1 rs1 0 11 A G\n1 rs2 0 20 C T\n",  # wrong position
        "1 rs1 0 10 A G\n1 rs1 0 20 C T\n",  # duplicate ID
    ],
)
def test_audit_rejects_incomplete_or_incompatible_inputs(audit_inputs, tmp_path, contents):
    cfg, instances = audit_inputs
    (tmp_path / "beta/beta_chr1.bim").write_text(contents)
    with pytest.raises(ValueError):
        f.audit_harmonization(cfg, instances)


def test_audit_requires_truth_for_every_selected_instance(audit_inputs):
    cfg, instances = audit_inputs
    with pytest.raises(ValueError, match="Missing causal"):
        f.audit_harmonization(cfg, instances.iloc[:1].assign(locus_id="unknown"))


@pytest.mark.parametrize(
    "counts,affected",
    [
        ({"anl": 0, "beta": 0}, (0, 10)),
        ({"anl": 2, "beta": 1}, (3, 10)),
        ({"anl": 2}, (10, 10)),
        ({"anl": 0}, (0, 0)),
        ({"anl": 0}, None),
    ],
)
def test_plot_handles_zero_mixed_and_missing_counts(tmp_path, monkeypatch, counts, affected):
    captured = []
    monkeypatch.setattr(f.fs, "save", lambda fig, path, log: captured.append(fig) or path)
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message="Attempting to set identical low and high ylims")
        f.fed4_harmonization(counts, 10, tmp_path / "harmonization.png", affected)
    fig = captured[0]
    try:
        assert all(ax.get_ylim()[0] >= 0 for ax in fig.axes)
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        if not any(counts.values()):
            assert "All sites match the reference" in texts
        if affected == (0, 10):
            assert "10 of 10 (100.0%)" in texts
            assert "require recoding" not in texts
            assert len(fig.axes[1].patches) == 1
        assert "phenotype stage was run before" not in " ".join(t.get_text() for t in fig.texts)
    finally:
        f.plt.close(fig)
