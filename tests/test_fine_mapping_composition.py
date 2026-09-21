"""Scientific controls for fixed-N allocation, selection and five-arm reporting."""

from pathlib import Path

import matplotlib.colors as colors
import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.arms import build_composed_site_dir
from appfl_bio_suite.experiments.fine_mapping.composition import rank_candidates, split_regions
from appfl_bio_suite.experiments.fine_mapping.reporting import paired_arm_differences


@pytest.fixture
def cohort(tmp_path):
    compositions = {
        "anl": {"AFR": 75, "EUR": 300},
        "covenant": {"AFR": 475, "EUR": 15},
        "mbzuai": {"AFR": 50, "EUR": 50},
    }
    for site, composition in compositions.items():
        path = tmp_path / "processed" / site
        path.mkdir(parents=True)
        pops = [pop for pop, n in composition.items() for _ in range(n)]
        pd.DataFrame(
            {
                "FID": [f"{site}-{i}" for i in range(len(pops))],
                "IID": [f"{site}-{i}" for i in range(len(pops))],
                "superpopulation": pops,
            }
        ).to_csv(path / f"{site}_manifest.tsv", sep="\t", index=False)
        (path / "genotypes.bed").write_bytes(b"original")
    return {
        "sites": compositions,
        "paths": {"processed_dir": str(tmp_path / "processed")},
        "fine_mapping": {"min_gwas_n": 10},
    }


@pytest.mark.parametrize("priority", [(), ("covenant", "anl", "mbzuai")])
def test_exact_composition_and_source_membership(cohort, tmp_path, priority):
    requested = {"AFR": 450, "EUR": 50}
    got = build_composed_site_dir(cohort, tmp_path / "draw", requested, site_priority=priority)
    assert {p: sum(c.get(p, 0) for c in got.values()) for p in requested} == requested
    for site in got:
        source = Path(cohort["paths"]["processed_dir"]) / site
        selected = pd.read_csv(tmp_path / "draw" / site / f"{site}_manifest.tsv", sep="\t")
        original = pd.read_csv(source / f"{site}_manifest.tsv", sep="\t")
        assert not selected.duplicated(["FID", "IID"]).any()
        assert set(selected.IID) <= set(original.IID)
        assert (tmp_path / "draw" / site / "genotypes.bed").resolve() == source / "genotypes.bed"
    if priority:
        assert got["covenant"] == {"AFR": 450, "EUR": 15}
    else:
        assert sum(got["anl"].values()) > 0
        assert sum(got["mbzuai"].values()) > 0


def test_sampling_is_reproducible_but_seeds_change_individuals(cohort, tmp_path):
    for draw, seed in [("a", 1), ("b", 1), ("c", 2)]:
        build_composed_site_dir(cohort, tmp_path / draw, {"AFR": 500}, seed)

    def ids(draw):
        return (tmp_path / draw / "anl/anl_manifest.tsv").read_text()

    assert ids("a") == ids("b")
    assert ids("a") != ids("c")


@pytest.mark.parametrize("counts", [{"AFR": 601}, {"EUR": 5}, {"AFR": 0}, {"AFR": 5.5}])
def test_impossible_or_unanalyzable_allocations_fail(cohort, tmp_path, counts):
    with pytest.raises(ValueError):
        build_composed_site_dir(cohort, tmp_path / "bad", counts)
    assert not (tmp_path / "bad").exists()


def test_overlapping_regions_never_cross_selection_boundary():
    loci = pd.DataFrame(
        {
            "locus_id": list("abcdefgh"),
            "chrom": 1,
            "start_bp": [0, 5, 100, 105, 200, 205, 300, 305],
            "end_bp": [10, 15, 110, 115, 210, 215, 310, 315],
            "stratum": "low",
        }
    )
    split = split_regions(loci)
    assert split.groupby("region").split.nunique().max() == 1
    assert set(split.split) == {"development", "evaluation"}
    pd.testing.assert_frame_equal(split, split_regions(loci.sample(frac=1, random_state=2)))


def test_selection_uses_power_and_deterministic_tiebreak():
    data = pd.DataFrame(
        {
            "candidate": ["b", "a", "c"],
            "any_causal_captured": [1, 1, 0],
            "causal_pip_max": [0.99, 0.1, 1.0],
            "locus_id": "l1",
        }
    )
    assert rank_candidates(data).candidate.tolist() == ["a", "b", "c"]


def arm_rows():
    rows = []
    for arm in ("anl", "covenant", "mbzuai", "federation_smart_50k", "federation", "ld_borrowed"):
        for locus in range(3):
            for seed in (1, 2, 3):
                rows.append(
                    dict(
                        arm=arm,
                        locus_id=f"L{locus}",
                        architecture_id="ncsl1_h2-0.001_rg1",
                        replicate=0,
                        sampling_seed=seed,
                        any_causal_captured=arm == "federation" or seed == 1,
                        n_credible_sets=1,
                        causal_pip_max=0.99 if seed == 1 else 0.2,
                        best_cs_size=3,
                        error="",
                        evaluation_split="evaluation",
                    )
                )
    return pd.DataFrame(rows)


def test_paired_seeds_are_not_cross_joined_or_counted_as_independent_loci():
    d = paired_arm_differences(arm_rows())
    assert (d.n_instances == 9).all()
    assert (d.n_loci == 3).all()
    np.testing.assert_allclose(d.power_difference, 2 / 3)
    with pytest.raises(ValueError, match="Duplicate"):
        paired_arm_differences(pd.concat([arm_rows()] * 2))


def test_fed6_has_exactly_five_bars_and_argonne_colors(tmp_path, monkeypatch):
    from appfl_bio_suite.experiments.fine_mapping.figures import federation as f

    captured = []
    monkeypatch.setattr(f.fs, "save", lambda fig, path, log: captured.append(fig) or path)
    f.fed6_what_federation_buys(arm_rows(), tmp_path / "fed6.png")
    figure = captured[0]
    assert len(figure.axes[0].patches) == len(figure.axes[2].patches) == 5
    assert len(figure.axes[3].get_yticklabels()) == 4
    assert [p.get_facecolor() for p in figure.axes[0].patches] == [
        colors.to_rgba(c) for c in f.arm_colors(f.ARM_ORDER)
    ]
    assert "ld_borrowed" not in " ".join(t.get_text() for t in figure.texts)
    f.plt.close(figure)
    with pytest.raises(ValueError, match="missing"):
        f.fed6_what_federation_buys(
            arm_rows().query("arm != 'federation_smart_50k'"), tmp_path / "bad.png"
        )
    bad = arm_rows()
    bad["evaluation_split"] = "development"
    with pytest.raises(ValueError, match="evaluation loci"):
        f.fed6_what_federation_buys(bad, tmp_path / "bad.png")


def test_composed_config_updates_dominant_ancestry(tmp_path):
    from test_fine_mapping_federated import SITES, _build_package

    from appfl_bio_suite.experiments.fine_mapping.arms import Arm, build_configs
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import load_config

    base = _build_package(tmp_path)
    available = sum(c.get("AFR", 0) for c in SITES.values())
    paths = build_configs(
        base,
        tmp_path / "arms",
        (Arm("afr_only", tuple(SITES), ancestry_counts=(("AFR", available // 2),)),),
    )
    config = load_config(paths["afr_only"])
    assert sum(s.n for s in config.sites.values()) == available // 2
    assert all(s.dominant == "AFR" and set(s.composition) == {"AFR"} for s in config.sites.values())


def test_composed_arm_runs_real_inference_without_resimulating(tmp_path, monkeypatch):
    from test_fine_mapping_federated import REPO_ROOT, SITES, _build_package

    from appfl_bio_suite.experiments.fine_mapping.arms import Arm, build_configs, run_arm
    from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import _resolve_binary

    try:
        for name in ("SuSiEx", "plink", "plink2"):
            _resolve_binary(name, REPO_ROOT)
    except FileNotFoundError as exc:
        pytest.skip(str(exc))
    import os

    monkeypatch.setenv("PATH", f"{REPO_ROOT / 'vendor/bin'}:{os.environ['PATH']}")
    base = _build_package(tmp_path)
    pheno = next((tmp_path / "data/ground_truth/phenotypes").rglob("*.pheno"))
    before = pheno.read_bytes()
    available = sum(c.get("AFR", 0) for c in SITES.values())
    paths = build_configs(
        base,
        tmp_path / "arms",
        (Arm("afr_only", tuple(SITES), ancestry_counts=(("AFR", available // 2),)),),
    )
    result = run_arm("afr_only", paths["afr_only"], n_workers=1)
    frame = pd.read_csv(result, sep="\t")
    assert len(frame) == 2
    assert frame.error.isna().all()
    assert frame.converged.all()
    assert pheno.read_bytes() == before


def test_submit_selection_precedes_evaluation_and_finish_waits_for_plots(tmp_path, monkeypatch):
    import importlib.util
    import json
    from types import SimpleNamespace

    script = Path(__file__).resolve().parents[1] / "scripts/fine-mapping/submit_composition.py"
    spec = importlib.util.spec_from_file_location("composition_submit", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    (tmp_path / "development_tasks.json").write_text(json.dumps([{}, {}]))
    calls = []

    def qsub(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(stdout=f"{len(calls)}.pbs\n", stderr="")

    monkeypatch.setattr(module.subprocess, "run", qsub)
    jobs = module.submit(tmp_path, "account", "/python")
    assert "depend=afterok:1.pbs" in jobs["select"]["command"]
    assert "depend=afterok:2.pbs" in jobs["evaluation"]["command"]
    assert "depend=afterok:3.pbs:4.pbs" in jobs["finish"]["command"]
    module.submit(tmp_path, "account", "/python")
    assert len(calls) == 5
