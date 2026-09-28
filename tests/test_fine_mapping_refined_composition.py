"""Fixed quotas, independent simulations, honest fallbacks and paired validation."""

import importlib.util
import json
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import yaml

from appfl_bio_suite.experiments.fine_mapping import refined_composition as r
from appfl_bio_suite.experiments.fine_mapping.arms import (
    Arm,
    build_composed_site_dir,
    build_configs,
)
from appfl_bio_suite.experiments.fine_mapping.figures.federation import DEFAULT_COMPOSITION
from appfl_bio_suite.experiments.fine_mapping.reporting import (
    clustered_ratio,
    paired_arm_differences,
)


def test_grid_has_exact_n_solo_fallbacks_and_preserves_source():
    cfg = dict(
        sites={s: dict(n=50000, composition=c) for s, c in DEFAULT_COMPOSITION.items()},
        fine_mapping=dict(min_gwas_n=1000),
    )
    before = deepcopy(cfg)
    choices = r.candidates(cfg)
    assert cfg == before
    assert len(choices) == 20
    assert {c["solo"] for c in choices if c["solo"]} == set(r.SITE_ORDER)
    quotas = set()
    for candidate in choices:
        assert sum(candidate["ancestry_counts"].values()) == 50000
        assert min(candidate["ancestry_counts"].values()) >= 1000
        quotas.add(json.dumps(candidate["site_counts"], sort_keys=True))
    assert len(quotas) == len(choices)
    assert any("AMR" in c["ancestry_counts"] for c in choices)
    assert any("MID" in c["ancestry_counts"] and not c["solo"] for c in choices)


@pytest.fixture
def source(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    compositions = {
        "anl": {"AFR": 20, "EUR": 30},
        "covenant": {"AFR": 40, "EUR": 10},
        "mbzuai": {"AFR": 30, "EUR": 20},
    }
    for site, counts in compositions.items():
        folder = root / "processed" / site
        folder.mkdir(parents=True)
        pops = [p for p, n in counts.items() for _ in range(n)]
        pd.DataFrame(
            dict(
                FID=[f"{site}-{i}" for i in range(50)],
                IID=[f"{site}-{i}" for i in range(50)],
                superpopulation=pops,
            )
        ).to_csv(folder / f"{site}_manifest.tsv", sep="\t", index=False)
        (folder / f"{site}_chr1.bed").write_bytes(b"unchanged source")
    (root / "loci").mkdir()
    pd.DataFrame(
        dict(
            locus_id=[f"L{i}" for i in range(6)],
            chrom=1,
            start_bp=np.arange(6) * 100,
            end_bp=np.arange(6) * 100 + 50,
            stratum=["low", "low", "medium", "medium", "high", "high"],
        )
    ).to_csv(root / "loci/selected_loci.tsv", sep="\t", index=False)
    cfg = dict(
        master_seed=7,
        chromosome=1,
        sites={
            s: dict(n=50, composition=c, dominant=max(c, key=c.get))
            for s, c in compositions.items()
        },
        paths={
            k + "_dir": str(root / k)
            for k in ("processed", "loci", "ground_truth", "reports", "logs")
        },
        architecture=dict(replicates=2, ancestry_divergent_causal=dict(enabled=True)),
        fine_mapping=dict(min_gwas_n=10),
    )
    (root / "pipeline_config.yaml").write_text(yaml.safe_dump(cfg))
    return root, cfg


def test_explicit_quotas_are_nested_and_config_matches_realized_cohort(source, tmp_path):
    root, cfg = source
    for name, afr in [("a", 30), ("b", 35)]:
        quotas = {"covenant": {"AFR": afr, "EUR": 10}, "anl": {"EUR": 40 - afr}}
        got = build_composed_site_dir(
            cfg, tmp_path / name, r.ancestry_totals(quotas), seed=9, site_counts=quotas
        )
        assert got == quotas

    def ids(name):
        return set(pd.read_csv(tmp_path / name / "covenant/covenant_manifest.tsv", sep="\t").IID)

    assert ids("a") < ids("b")
    quotas = {"covenant": {"AFR": 35, "EUR": 10}, "anl": {"EUR": 5}}
    arm = Arm(
        "quota",
        ("anl", "covenant"),
        site_counts=tuple((s, tuple(c.items())) for s, c in quotas.items()),
    )
    path = build_configs(root / "pipeline_config.yaml", tmp_path / "arms", (arm,))["quota"]
    config = yaml.safe_load(path.read_text())
    assert {s: c["composition"] for s, c in config["sites"].items()} == quotas
    assert sum(c["n"] for c in config["sites"].values()) == 50
    assert config["paths"]["ground_truth_dir"] == cfg["paths"]["ground_truth_dir"]


@pytest.mark.parametrize(
    "quotas,totals",
    [
        ({"covenant": {"AFR": 41}}, {"AFR": 41}),
        ({"missing": {"AFR": 20}}, {"AFR": 20}),
        ({"covenant": {"AFR": -1}}, {"AFR": 20}),
        ({"covenant": {"AFR": 20}}, {"AFR": 21}),
        ({"covenant": {"EUR": 5}}, {"EUR": 5}),
    ],
)
def test_invalid_quotas_fail_before_writing(source, tmp_path, quotas, totals):
    with pytest.raises(ValueError):
        build_composed_site_dir(source[1], tmp_path / "bad", totals, site_counts=quotas)
    assert not (tmp_path / "bad").exists()


def choices_for_test(cfg):
    choices = [
        dict(
            name="solo_" + s,
            solo=s,
            site_counts={s: c["composition"]},
            ancestry_counts=c["composition"],
        )
        for s, c in cfg["sites"].items()
    ]
    quotas = {"covenant": {"AFR": 35, "EUR": 10}, "anl": {"EUR": 5}}
    return choices + [
        dict(name="mixed", solo=None, site_counts=quotas, ancestry_counts=r.ancestry_totals(quotas))
    ]


@pytest.fixture
def prepared(source, tmp_path, monkeypatch):
    from appfl_bio_suite.experiments.fine_mapping import rerun

    monkeypatch.setattr(r, "candidates", choices_for_test)
    monkeypatch.setattr(r, "load_config", lambda _: None)
    monkeypatch.setattr(rerun, "check_cohort", lambda _: None)
    root = tmp_path / "run"
    before = (source[0] / "pipeline_config.yaml").read_bytes()
    r.prepare(source[0], root)
    assert before == (source[0] / "pipeline_config.yaml").read_bytes()
    return root


def test_prepare_is_shared_only_and_both_phases_have_new_independent_draws(prepared):
    configs = [
        yaml.safe_load(r.phase_config(prepared, phase).read_text())
        for phase in ("development", "evaluation")
    ]
    assert configs[0]["master_seed"] != configs[1]["master_seed"] != 7
    assert [c["architecture"]["replicates"] for c in configs] == [5, 10]
    for cfg in configs:
        assert cfg["architecture"]["ancestry_divergent_causal"]["enabled"] is False
        gt = Path(cfg["paths"]["ground_truth_dir"])
        assert not gt.is_symlink() and not list(gt.iterdir())
    assert not (prepared / "evaluation_tasks.json").exists()
    split = pd.read_csv(prepared / "split.tsv", sep="\t")
    assert split.groupby("region").split.nunique().max() == 1
    assert set(split.split) == {"development", "evaluation"}


def ranking_rows():
    rows = []
    for name, hits, coverage, false in [
        ("solo_covenant", 8, 96, 2),
        ("mixed", 7, 96, 2),
        ("bad_calibration", 9, 80, 20),
    ]:
        for i in range(10):
            rows.append(
                dict(
                    candidate=name,
                    any_causal_captured=i < hits,
                    causal_pip_max=0.99,
                    locus_id=f"L{i}",
                    n_cs_containing_causal=coverage,
                    n_credible_sets=100,
                    n_true_pip95=100 - false,
                    n_false_pip95=false,
                )
            )
    return pd.DataFrame(rows)


def test_solo_fallback_wins_and_calibration_cannot_be_traded_for_power():
    choices = [
        dict(name=n, solo="covenant" if n == "solo_covenant" else None)
        for n in ranking_rows().candidate.unique()
    ]
    ranks = r.rank_candidates(
        ranking_rows(), choices, dict(min_cs_coverage=0.93, max_fdr_pip95=0.05)
    )
    assert ranks.iloc[0].candidate == "solo_covenant"
    assert not ranks.iloc[-1].eligible
    rows = ranking_rows()
    rows.loc[rows.candidate.eq("mixed"), "any_causal_captured"] = rows.loc[
        rows.candidate.eq("solo_covenant"), "any_causal_captured"
    ].to_numpy()
    assert (
        r.rank_candidates(rows, choices, dict(min_cs_coverage=0.93, max_fdr_pip95=0.05))
        .iloc[0]
        .solo_fallback
    )


def test_evaluation_tasks_use_only_fresh_truth_and_rerun_every_baseline(prepared):
    design = r.read_json(prepared / "design.json")
    winner = next(c for c in design["candidates"] if c["name"] == "mixed")
    r.build_tasks(prepared, "evaluation", [winner])
    tasks = r.read_json(prepared / "evaluation_tasks.json")
    assert len(tasks) == 42
    assert {t["candidate"] for t in tasks} == {"mixed", *r.SITE_ORDER, "federation"}
    for task in tasks:
        cfg = yaml.safe_load(Path(task["config"]).read_text())
        assert cfg["master_seed"] == r.PHENOTYPE_SEEDS["evaluation"]
        assert cfg["paths"]["ground_truth_dir"] == str(prepared / "evaluation/ground_truth")
        assert Path(cfg["paths"]["reports_dir"]).is_relative_to(prepared / "evaluation")
    r.verify_plan(prepared, "evaluation")
    Path(tasks[0]["config"]).write_text("changed")
    with pytest.raises(ValueError, match="allocation changed"):
        r.verify_plan(prepared, "evaluation")


def test_no_validation_simulation_before_selection(prepared):
    with pytest.raises(FileNotFoundError):
        r.simulate(prepared, "evaluation", 0)


def test_region_bootstrap_keeps_overlapping_loci_together():
    frame = pd.DataFrame(
        dict(locus_id=["a", "b", "c", "d"], region=[0, 0, 1, 1], hit=[0, 0, 1, 1], one=1)
    )
    result = clustered_ratio(frame, "hit", "one")
    collapsed = (
        frame.groupby("region", as_index=False)[["hit", "one"]]
        .sum()
        .rename(columns={"region": "locus_id"})
    )
    assert result == clustered_ratio(collapsed, "hit", "one")
    paired = pd.concat(
        [
            frame.assign(arm="federation", architecture_id="a", replicate=0, any_causal_captured=1),
            frame.assign(
                arm="smart", architecture_id="a", replicate=0, any_causal_captured=frame.hit
            ),
        ]
    )
    comparison = paired_arm_differences(paired).iloc[0]
    assert comparison.power_difference == 0.5
    assert comparison.uncertainty == "paired_region_cluster_bootstrap"


def test_scheduler_waits_for_selection_and_validation(prepared, monkeypatch):
    path = (
        Path(__file__).resolve().parents[1] / "scripts/fine-mapping/submit_refined_composition.py"
    )
    spec = importlib.util.spec_from_file_location("refined_submit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    calls = []

    def qsub(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(stdout=f"{len(calls)}.pbs\n", stderr="")

    monkeypatch.setattr(module.subprocess, "run", qsub)
    jobs = module.submit(prepared, "test-allocation", "/python")
    assert len(calls) == 8
    assert "depend=afterok:2.pbs" in jobs["development"]["command"]
    assert "depend=afterok:4.pbs" in jobs["evaluation-sim"]["command"]
    assert "depend=afterok:6.pbs" in jobs["evaluation"]["command"]
    assert "depend=afterok:7.pbs" in jobs["finish"]["command"]
    assert "0-41%3" in jobs["evaluation"]["command"]
    module.submit(prepared, "test-allocation", "/python")
    assert len(calls) == 8


def test_fresh_simulation_uses_new_effects_and_preserves_source_phenotypes(tmp_path):
    from test_fine_mapping_federated import _build_package

    from appfl_bio_suite.experiments.fine_mapping.arms import _load
    from appfl_bio_suite.experiments.fine_mapping.fedfm.phenotype_sim import run_phenotype_sim
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import load_config
    from appfl_bio_suite.experiments.fine_mapping.fedfm.validation import _rederive_row

    source = _build_package(tmp_path / "source")
    raw = _load(source)
    old_gt = Path(raw["paths"]["ground_truth_dir"])
    old_hashes = {p: r.sha(p) for p in old_gt.rglob("*.pheno")}
    raw["master_seed"] = r.PHENOTYPE_SEEDS["evaluation"]
    raw["paths"]["ground_truth_dir"] = str(tmp_path / "fresh_truth")
    raw["architecture"]["ancestry_divergent_causal"] = dict(enabled=False)
    config = tmp_path / "fresh.yaml"
    config.write_text(yaml.safe_dump(raw))
    cfg = load_config(config)
    run_phenotype_sim(cfg)
    truth = pd.read_csv(tmp_path / "fresh_truth/causal_manifest.tsv", sep="\t")
    assert len(truth) == 2 and truth.causal_mode.eq("shared").all()
    assert all(r.sha(p) == h for p, h in old_hashes.items())
    new_file = next((tmp_path / "fresh_truth/phenotypes").rglob("*.pheno"))
    assert r.sha(new_file) != r.sha(old_gt / new_file.relative_to(tmp_path / "fresh_truth"))
    for row in truth.to_dict("records"):
        assert all(result["status"] == "ok" for result in _rederive_row(str(config), row))


@pytest.mark.parametrize("winner_name", ["mixed", "solo_covenant"])
def test_headline_pairs_fresh_baselines_and_keeps_solo_fallback_honest(prepared, winner_name):
    design = r.read_json(prepared / "design.json")
    winner = next(c for c in design["candidates"] if c["name"] == winner_name)
    r._json(prepared / "selection.json", winner)
    loci = pd.read_csv(prepared / "evaluation/loci/selected_loci.tsv", sep="\t")
    truth = pd.DataFrame(
        dict(locus_id=loci.locus_id, architecture_id="ncsl1_h2-0.001_rg1", replicate=0)
    )
    truth.to_csv(prepared / "evaluation/ground_truth/causal_manifest.tsv", sep="\t", index=False)
    rows = []
    for candidate in [winner_name, *r.SITE_ORDER, "federation"]:
        seeds = (
            r.SAMPLING_SEEDS
            if candidate == winner_name and winner["solo"] is None
            else r.SAMPLING_SEEDS[:1]
        )
        for seed in seeds:
            for row in truth.to_dict("records"):
                rows.append(
                    row
                    | dict(
                        candidate=candidate,
                        sampling_seed=seed,
                        any_causal_captured=True,
                        causal_pip_max=0.99,
                        n_cs_containing_causal=1,
                        n_credible_sets=1,
                        n_true_pip95=1,
                        n_false_pip95=0,
                        error="",
                        fit_status="converged_cs",
                    )
                )
    raw = pd.DataFrame(rows)
    headline = r.assemble_headline(prepared, raw)
    assert len(headline) == 5 * 3 * len(truth)
    assert headline.region.notna().all()
    assert not headline.duplicated(["arm", *r.KEY, "sampling_seed"]).any()
    assert len(r.metric_summary(headline, ["arm"])) == 20
    assert len(paired_arm_differences(headline)) == 4
    if winner["solo"]:
        label = headline.loc[headline.arm.eq(r.original.SMART_ARM), "smart_composition"].iloc[0]
        assert "solo fallback: covenant" in label
    with pytest.raises(ValueError, match="fresh baselines"):
        r.assemble_headline(prepared, raw.loc[raw.candidate.ne("federation")])


def test_selection_freezes_solo_fallback_before_validation_exists(prepared, monkeypatch):
    r._json(
        prepared / "development/phenotypes.accepted.json",
        dict(files={"pipeline_config.yaml": r.sha(r.phase_config(prepared, "development"))}),
    )
    data = ranking_rows().query("candidate != 'bad_calibration'").copy()

    def collect(root, phase):
        assert phase == "development"
        assert not list((root / "evaluation/ground_truth").iterdir())
        data.to_csv(root / "development_results.tsv", sep="\t", index=False)
        return data

    monkeypatch.setattr(r.original, "collect", collect)
    r.select(prepared)
    selected = r.read_json(prepared / "selection.json")
    assert selected["solo"] == "covenant"
    assert selected["evaluation_used_for_selection"] is False
    assert len(r.read_json(prepared / "evaluation_tasks.json")) == 30
    r.verify_selection(prepared)
    with pytest.raises(FileExistsError, match="already frozen"):
        r.select(prepared)
    (prepared / "selection.json").write_text("{}")
    with pytest.raises(ValueError, match="selection changed"):
        r.verify_selection(prepared)
