"""Regression tests for failures found in the September scientific audit."""

import numpy as np
import pandas as pd

from appfl_bio_suite.experiments.fine_mapping.inference import (
    inclusion_probabilities,
    parse_outputs,
)
from appfl_bio_suite.experiments.fine_mapping.reporting import clustered_ratio, parity_check


def test_overall_pip_can_cross_threshold_when_each_component_does_not():
    df = pd.DataFrame({"PIP(CS1)": [0.4, 0.9], "PIP(CS2)": [0.4, 0.0]})
    np.testing.assert_allclose(inclusion_probabilities(df), [0.64, 0.9])
    assert inclusion_probabilities(df).iloc[0] > 0.5


def test_null_and_fail_are_distinct_and_excluded_truth_is_explicit(tmp_path):
    (tmp_path / "cs.snp").write_text("SNP\neligible\n")
    (tmp_path / "cs.cs").write_text("NULL\n")
    (tmp_path / "cs.summary").write_text("# region\nNULL\n")
    result = parse_outputs(tmp_path, "cs", ["eligible", "excluded"])
    assert result["converged"] and result["fit_status"] == "converged_no_cs"
    assert result["causal_pip_mean"] == 0
    assert result["n_eligible_causal"] == result["n_excluded_causal"] == 1
    (tmp_path / "cs.summary").write_text("# region\nFAIL\n")
    result = parse_outputs(tmp_path, "cs", ["eligible"])
    assert not result["converged"] and result["fit_status"] == "nonconverged"
    assert np.isnan(result["causal_pip_mean"])


def test_phenotype_permutation_fails_even_when_heritability_is_unchanged(tmp_path, monkeypatch):
    from appfl_bio_suite.experiments.fine_mapping.fedfm import validation as v
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import derive_seed

    n = 100
    dosage = np.tile([0.0, 1.0, 2.0, 1.0], n // 4).reshape(-1, 1)
    beta = np.array([[0.7]])
    g = dosage[:, 0] * beta[0, 0]
    target = 0.01
    seed = 123
    row = {
        "locus_id": "L0",
        "architecture_id": "a",
        "replicate": 0,
        "ncsl": 1,
        "rg": 1.0,
        "h2_target": target,
        "causal_snp_ids": "snp",
    }
    rng = np.random.default_rng(derive_seed(seed, "noise", "site", "L0", "a", 0))
    y = g + rng.normal(0, np.sqrt(np.var(g) * (1 - target) / target), n)
    row["empirical_h2_site"] = np.var(g) / np.var(y)
    ids = [str(i) for i in range(n)]
    io = {
        "prefix": "unused",
        "n": n,
        "pop_index": np.zeros(n, dtype=int),
        "dom_idx": 0,
        "fam_iid": ids,
        "fam_keys": list(zip(ids, ids, strict=True)),
        "master_seed": seed,
        "pop_names": ["EUR"],
        "snpid_to_vidx": pd.Series([0], index=["snp"]),
    }
    monkeypatch.setattr(v, "_site_io", lambda _: ({"site": io}, ("site",), str(tmp_path)))
    monkeypatch.setattr(v, "read_bed_variants", lambda *args: dosage)
    (tmp_path / "effect_sizes").mkdir()
    (tmp_path / "phenotypes/site").mkdir(parents=True)
    np.save(tmp_path / "effect_sizes/L0_a_rep0.npy", beta)
    path = tmp_path / "phenotypes/site/L0_a_rep0.pheno"
    frame = pd.DataFrame({"FID": ids, "IID": ids, "y": y})
    frame.to_csv(path, sep="\t", index=False)
    assert v._rederive_row("unused", row)[0]["status"] == "ok"
    frame["y"] = y[::-1]
    frame.to_csv(path, sep="\t", index=False)
    broken = v._rederive_row("unused", row)[0]
    assert broken["h2_relerr"] < 1e-12
    assert broken["status"] == "phenotype_vector_mismatch"


def test_matched_arm_reallocates_excluded_ancestry_and_preserves_exact_n(tmp_path):
    from appfl_bio_suite.experiments.fine_mapping.arms import build_downsampled_site_dir

    comps = {
        "anl": {"EUR": 300, "AFR": 75, "AMR": 75, "EAS": 25, "CSA": 25},
        "covenant": {"AFR": 475, "EUR": 15, "CSA": 10},
        "mbzuai": {"MID": 250, "CSA": 150, "AFR": 50, "EUR": 50},
    }
    for site, comp in comps.items():
        dest = tmp_path / "processed" / site
        dest.mkdir(parents=True)
        pops = [p for p, n in comp.items() for _ in range(n)]
        pd.DataFrame(
            {
                "FID": [f"{site}-{i}" for i in range(500)],
                "IID": [f"{site}-{i}" for i in range(500)],
                "superpopulation": pops,
            }
        ).to_csv(dest / f"{site}_manifest.tsv", sep="\t", index=False)
    cfg = {
        "paths": {"processed_dir": str(tmp_path / "processed")},
        "sites": comps,
        "fine_mapping": {"min_gwas_n": 10},
    }
    got = build_downsampled_site_dir(cfg, tmp_path / "subset", 1 / 3)
    pooled = {}
    for comp in got.values():
        for pop, n in comp.items():
            pooled[pop] = pooled.get(pop, 0) + n
    assert "EAS" not in pooled
    assert sum(pooled.values()) == sum(n for n in pooled.values() if n >= 10) == 500


def test_replicate_duplication_does_not_create_independent_loci():
    df = pd.DataFrame(
        {"locus_id": np.repeat(list("abcd"), 2), "yes": [0, 0, 1, 1, 0, 1, 1, 1], "one": 1}
    )
    a = clustered_ratio(df, "yes", "one")
    b = clustered_ratio(pd.concat([df] * 10), "yes", "one")
    np.testing.assert_array_equal(a, b)


def test_central_and_federated_remain_equivalent_with_missing_genotypes(tmp_path):
    from test_fine_mapping_federated import CHROM, M, _build_package

    from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import run_fed_fine_mapping
    from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import run_fine_mapping
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import load_config

    cfg = load_config(_build_package(tmp_path), repo_root=tmp_path)
    site = next(iter(cfg.sites))
    prefix = cfg.site_dir(site) / f"{site}_chr{CHROM}"
    fam = pd.read_csv(prefix.with_suffix(".fam"), sep=r"\s+", header=None)
    iid = fam.iloc[0, 1]
    pool = cfg.resolved_path("hapnest_dir") / f"chr{CHROM}"
    pfam = pd.read_csv(pool.with_suffix(".fam"), sep=r"\s+", header=None)
    for pref, n, row in [
        (prefix, len(fam), 0),
        (pool, len(pfam), int(np.flatnonzero(pfam[1] == iid)[0])),
    ]:
        offset = 3 + (M - 1) * ((n + 3) // 4) + row // 4
        with pref.with_suffix(".bed").open("r+b") as handle:
            handle.seek(offset)
            old = handle.read(1)[0]
            handle.seek(offset)
            handle.write(bytes([(old & ~(3 << (2 * (row % 4)))) | (1 << (2 * (row % 4)))]))
    run_fine_mapping(cfg, n_workers=1)
    run_fed_fine_mapping(cfg, n_workers=1)
    result = parity_check(
        tmp_path / "reports/fine_mapping",
        tmp_path / "reports/fed_fine_mapping",
        tmp_path / "parity.json",
    )
    assert result["pass"] and result["n_instances"] == 2


def test_regenerated_phenotypes_rederive_and_preserve_every_causal_maf(tmp_path):
    import json

    from test_fine_mapping_federated import _build_package

    from appfl_bio_suite.experiments.fine_mapping.fedfm.phenotype_sim import run_phenotype_sim
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import load_config
    from appfl_bio_suite.experiments.fine_mapping.fedfm.validation import _rederive_row

    path = _build_package(tmp_path)
    cfg = load_config(path, repo_root=tmp_path)
    run_phenotype_sim(cfg)
    gt = cfg.resolved_path("ground_truth_dir")
    manifest = pd.read_csv(gt / "causal_manifest.tsv", sep="\t")
    for row in manifest.to_dict("records"):
        inst = f"{row['locus_id']}_{row['architecture_id']}_rep{row['replicate']}"
        beta = np.load(gt / "effect_sizes" / f"{inst}.npy")
        assert beta.dtype == np.float64
        mafs = json.loads(row["per_site_mafs_json"])
        assert all(len(values) == len(beta) for values in mafs.values())
        assert all(r["phenotype_vector_match"] for r in _rederive_row(str(path), row))


def test_submit_graph_orders_simulation_validation_analysis_and_plots(tmp_path, monkeypatch):
    import importlib.util
    import json
    from pathlib import Path
    from types import SimpleNamespace

    script = Path(__file__).resolve().parents[1] / "scripts/fine-mapping/submit_scientific_rerun.py"
    spec = importlib.util.spec_from_file_location("fm_submit_test", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    (tmp_path / "analysis_tasks.json").write_text(json.dumps([{}, {}]))
    calls = []

    def qsub(args, **kwargs):
        calls.append(args)
        suffix = "[]" if "-J" in args else ""
        return SimpleNamespace(returncode=0, stdout=f"{len(calls)}{suffix}.pbs\n", stderr="")

    monkeypatch.setattr(module.subprocess, "run", qsub)
    module.submit(tmp_path, "test", "/python")
    jobs = json.loads((tmp_path / "jobs.json").read_text())
    stages = [
        "rescore",
        "merge-rescore",
        "phenotypes",
        "merge-phenotypes",
        "validation",
        "merge-validation",
        "bundle",
        "analysis",
        "finish",
        "plots",
    ]
    assert list(jobs) == stages
    assert all("place=exclhost" in job["command"] for job in jobs.values())
    assert all("select=1:ngpus=8:ncpus=256:mem=960gb" in job["command"] for job in jobs.values())
    for before, after in zip(stages, stages[1:], strict=False):
        assert f"depend=afterok:{jobs[before]['id']}" in jobs[after]["command"]
    assert "0-1%3" in jobs["analysis"]["command"]
    module.submit(tmp_path, "test", "/python")
    assert len(calls) == 10  # Resuming submission never duplicates recorded jobs.


def test_credible_set_detail_does_not_mix_shared_and_divergent_instances(tmp_path):
    from appfl_bio_suite.experiments.fine_mapping.figures.paper_plots import (
        load_cs_detail,
        parse_architecture,
    )

    aid = "ncsl2_h2-0.005_rg1"
    pd.DataFrame(
        {"instance": [f"L0000_{aid}_rep0", f"L0001_{aid}_rep0"], "contains_causal": [True, False]}
    ).to_csv(tmp_path / "fm_credible_sets.tsv", sep="\t", index=False)
    shared = pd.DataFrame(
        {
            "locus_id": ["L0000"],
            "architecture_id": [aid],
            "replicate": [0],
            "stratum": ["high"],
            "any_causal_captured": [True],
            "n_credible_sets": [1],
            "error": [""],
        }
    )
    detail = load_cs_detail(tmp_path, parse_architecture(shared))
    assert len(detail) == 1 and detail.contains_causal.all()


def test_sumstats_ld_variant_gate_refuses_mismatch_before_inference(tmp_path):
    import pytest

    from appfl_bio_suite.experiments.fine_mapping.inference import variant_input_manifest

    ss = tmp_path / "EUR.sumstats"
    ss.write_text("chr\tsnp\tbp\tA1\tA2\n1\ts1\t10\tA\tG\n1\ts2\t20\tC\tT\n")
    prefix = tmp_path / "EUR_ld"
    bim = tmp_path / "EUR_ld_ref.bim"
    bim.write_text("1\ts1\t0\t10\tA\tG\n1\ts2\t0\t20\tC\tT\n")
    good = variant_input_manifest([ss], [prefix], [100])[0]
    assert good["sumstats_variants_sha256"] == good["ld_variants_sha256"]
    bim.write_text("1\ts1\t0\t10\tA\tG\n")
    with pytest.raises(ValueError, match="refusing inference"):
        variant_input_manifest([ss], [prefix], [100])
    # The explicitly distinct external-reference sensitivity arm remains auditable.
    borrowed = variant_input_manifest([ss], [prefix], [100], require_matching=False)[0]
    assert borrowed["sumstats_n_variants"] == 2 and borrowed["ld_n_variants"] == 1


def test_sumstats_harmonization_preserves_test_but_reverses_signed_effect(tmp_path):
    import pytest

    from appfl_bio_suite.experiments.fine_mapping.inference import harmonize_sumstats

    bim = tmp_path / "ref.bim"
    bim.write_text("1\ts1\t0\t10\tA\tG\n")
    frame = pd.DataFrame(
        {
            "chr": [1],
            "snp": ["s1"],
            "bp": [10],
            "A1": ["G"],
            "A2": ["A"],
            "beta": [0.4],
            "stat": [2.0],
            "se": [0.2],
            "p": [0.05],
        }
    )
    got = harmonize_sumstats(frame, bim)
    assert got.A1.iloc[0] == "A" and got.A2.iloc[0] == "G"
    assert got.beta.iloc[0] == -0.4 and got.stat.iloc[0] == -2
    assert got.se.iloc[0] == 0.2 and got.p.iloc[0] == 0.05
    with pytest.raises(ValueError, match="incompatible"):
        harmonize_sumstats(frame.assign(A1="T"), bim)


def test_retry_preserves_history_and_refuses_live_job_graph(tmp_path, monkeypatch):
    import importlib.util
    import json
    from pathlib import Path

    import pytest

    script = Path(__file__).resolve().parents[1] / "scripts/fine-mapping/submit_scientific_rerun.py"
    spec = importlib.util.spec_from_file_location("fm_retry_test", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    manifest = tmp_path / "jobs.json"
    original = json.dumps({"rescore": {"id": "123[].pbs"}})
    manifest.write_text(original)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "    job_state = R\n")
    with pytest.raises(RuntimeError, match="not terminal"):
        module.archive_terminal_attempt(tmp_path)
    assert manifest.read_text() == original
    cache = tmp_path / "cached.npz"
    cache.write_bytes(b"previous scientific work")
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "    job_state = F\n")
    archived = module.archive_terminal_attempt(tmp_path)
    assert (archived / "jobs.json").read_text() == original
    assert "job_state = F" in json.loads((archived / "scheduler.json").read_text())["rescore"]
    assert not manifest.exists() and cache.read_bytes() == b"previous scientific work"
