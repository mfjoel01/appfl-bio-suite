"""Unit tests for src.fine_mapping.

The risky, version-sensitive surface is parsing SuSiEx's ``.cs`` / ``.snp`` /
``.summary`` output, so these tests exercise the parser on synthetic files that
mimic real SuSiEx output (perfect single-CS recovery, a multi-CS case, and the
"no credible set" case). No plink/SuSiEx binaries are invoked.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import (
    _LD_SUFFIXES,
    _METRIC_COLS,
    FMColumn,
    _ld_ready,
    _pip_by_snp,
    _resolve_binary,
    parse_susiex,
    plan_columns,
)

# mirrors config/simulation_config.yaml compositions + superpop order
_SUPERPOPS = ["EUR", "AFR", "AMR", "EAS", "CSA", "MID"]
_COMPOSITION = {
    "anl": {"EUR": 30000, "AFR": 7500, "AMR": 7500, "EAS": 2500, "CSA": 2500},
    "covenant": {"AFR": 47500, "EUR": 1500, "CSA": 1000},
    "mbzuai": {"MID": 25000, "CSA": 15000, "AFR": 5000, "EUR": 5000},
}


def _write(p: Path, text: str) -> None:
    p.write_text(text.lstrip("\n"))


def _write_case(d: Path, cs: str, snp: str, summary: str) -> None:
    _write(d / "cs.cs", cs)
    _write(d / "cs.snp", snp)
    _write(d / "cs.summary", summary)


# --- perfect single-CS recovery (mirrors the L0000 smoke result) ------------

CS_ONE = """
CS_ID\tSNP\tBP\tREF_ALLELE\tALT_ALLELE\tREF_FRQ\tBETA\tSE\t-LOG10P\tCS_PIP\tOVRL_PIP
1\tchr1:150276359:A:T\t150276359\tA,NA,A\tT,NA,T\t0.97,NA,0.98\t-0.32,NA,-0.50\t0.019,NA,0.028\t61.6,NA,71.0\t1\t1
"""
SNP_ONE = """
BP\tSNP\tPIP(CS1)\tLogBF(CS1,Pop1)\tLogBF(CS1,Pop2)\tLogBF(CS1,Pop3)
150010529\tchr1:150010529:C:T\t1.08e-130\t-2.37\t-3.2e-10\t-2.33
150276359\tchr1:150276359:A:T\t1\t136.1\t-3.2e-10\t158.4
"""
SUMMARY_ONE = """
# chr1:150010180-151010179
CS_ID\tCS_LENGTH\tCS_PURITY\tMAX_PIP_SNP\tBP\tREF_ALLELE\tALT_ALLELE\tREF_FRQ\tBETA\tSE\t-LOG10P\tMAX_PIP\tPP1\tPP2\tPP3
1\t1\t1\tchr1:150276359:A:T\t150276359\tA,NA,A\tT,NA,T\t0.97,NA,0.98\t-0.32,NA,-0.50\t0.019,NA,0.028\t61.6,NA,71.0\t1\t1\t0.5\t1
"""


def test_parse_perfect_recovery(tmp_path: Path) -> None:
    _write_case(tmp_path, CS_ONE, SNP_ONE, SUMMARY_ONE)
    m = parse_susiex(tmp_path, "cs", ["chr1:150276359:A:T"])
    assert m["n_credible_sets"] == 1
    assert m["best_cs_size"] == 1
    assert m["any_causal_captured"] is True
    assert m["n_causal_captured"] == 1
    assert m["causal_pip_max"] == pytest.approx(1.0)
    assert m["mean_cs_purity"] == pytest.approx(1.0)
    assert m["converged"] is True


# --- two credible sets, causal in the larger one ----------------------------

CS_TWO = """
CS_ID\tSNP\tBP\tCS_PIP\tOVRL_PIP
1\tchr1:100:A:G\t100\t0.9\t0.9
1\tchr1:200:C:T\t200\t0.6\t0.6
2\tchr1:300:G:A\t300\t0.8\t0.8
"""
SNP_TWO = """
BP\tSNP\tPIP(CS1)\tPIP(CS2)
100\tchr1:100:A:G\t0.9\t0.0
200\tchr1:200:C:T\t0.6\t0.0
300\tchr1:300:G:A\t0.0\t0.8
"""
SUMMARY_TWO = """
# region
CS_ID\tCS_LENGTH\tCS_PURITY\tMAX_PIP
1\t2\t0.7\t0.9
2\t1\t0.95\t0.8
"""


def test_parse_two_cs_best_size_and_pip(tmp_path: Path) -> None:
    # truth is the single SNP in CS2 (size 1) -> best_cs_size must be 1, not 2
    _write_case(tmp_path, CS_TWO, SNP_TWO, SUMMARY_TWO)
    m = parse_susiex(tmp_path, "cs", ["chr1:300:G:A"])
    assert m["n_credible_sets"] == 2
    assert m["total_cs_snps"] == 3
    assert m["best_cs_size"] == 1
    assert m["any_causal_captured"] is True
    assert m["causal_pip_max"] == pytest.approx(0.8)
    assert m["top_pip"] == pytest.approx(0.9)
    assert m["min_cs_purity"] == pytest.approx(0.7)


def test_parse_causal_not_captured(tmp_path: Path) -> None:
    _write_case(tmp_path, CS_TWO, SNP_TWO, SUMMARY_TWO)
    m = parse_susiex(tmp_path, "cs", ["chr1:999:T:C"])  # absent variant
    assert m["any_causal_captured"] is False
    assert m["n_causal_captured"] == 0
    import math
    assert math.isnan(m["best_cs_size"])
    assert m["causal_pip_max"] == 0.0
    assert m["n_excluded_causal"] == 1


# --- no credible set (weak signal): header-only .cs, no .snp -----------------

def test_parse_no_credible_set(tmp_path: Path) -> None:
    _write(tmp_path / "cs.cs", "CS_ID\tSNP\tBP\tCS_PIP\tOVRL_PIP\n")
    _write(tmp_path / "cs.summary", "# region\nNULL\n")
    # no cs.snp file at all
    m = parse_susiex(tmp_path, "cs", ["chr1:300:G:A"])
    assert m["n_credible_sets"] == 0
    assert m["total_cs_snps"] == 0
    assert m["any_causal_captured"] is False
    assert m["converged"] is False


def test_pip_by_snp_missing_file(tmp_path: Path) -> None:
    assert _pip_by_snp(tmp_path / "nope.snp") == {}


def test_pip_by_snp_max_over_cs_columns(tmp_path: Path) -> None:
    _write(tmp_path / "cs.snp", SNP_TWO)
    pm = _pip_by_snp(tmp_path / "cs.snp")
    assert pm["chr1:100:A:G"] == pytest.approx(0.9)
    assert pm["chr1:300:G:A"] == pytest.approx(0.8)


def test_metric_cols_match_parse_susiex(tmp_path: Path) -> None:
    # the failure path fills _METRIC_COLS by hand; keep it in sync with the
    # success path so both produce the same columns
    _write_case(tmp_path, CS_ONE, SNP_ONE, SUMMARY_ONE)
    assert set(parse_susiex(tmp_path, "cs", ["chr1:100:A:G"])) == set(_METRIC_COLS)


def test_ld_ready_needs_every_suffix(tmp_path: Path) -> None:
    # SuSiEx only takes its precomputed-LD path when all three files are there;
    # a partial set must re-trigger the build rather than be treated as cached
    prefix = tmp_path / "site_ld"
    assert not _ld_ready(prefix)
    for suffix in _LD_SUFFIXES[:-1]:
        _write(Path(f"{prefix}{suffix}"), "x")
    assert not _ld_ready(prefix)
    _write(Path(f"{prefix}{_LD_SUFFIXES[-1]}"), "x")
    assert _ld_ready(prefix)


# --- SuSiEx column planning ------------------------------------------------

def test_plan_columns_one_per_ancestry_pooled() -> None:
    cols = plan_columns(_COMPOSITION, _SUPERPOPS, min_n=1000)
    # one column per superpop (NOT per site x superpop), in _SUPERPOPS order
    assert [c.name for c in cols] == ["EUR", "AFR", "AMR", "EAS", "CSA", "MID"]
    assert len(cols) == 6


def test_plan_columns_pools_n_across_sites() -> None:
    cols = {c.name: c for c in plan_columns(_COMPOSITION, _SUPERPOPS, min_n=0)}
    # the whole point: an ancestry's n is the sum over every site holding it
    assert cols["AFR"].n == 7500 + 47500 + 5000 == 60000
    assert cols["EUR"].n == 30000 + 1500 + 5000 == 36500
    assert cols["CSA"].n == 2500 + 1000 + 15000 == 18500
    # single-site ancestries pool to themselves
    assert cols["AMR"].n == 7500 and cols["AMR"].sites == ("anl",)
    assert cols["MID"].n == 25000 and cols["MID"].sites == ("mbzuai",)
    assert cols["AFR"].sites == ("anl", "covenant", "mbzuai")
    # every enrolled individual is covered exactly once
    assert sum(c.n for c in cols.values()) == 150000


def test_plan_columns_min_n_gate_drops_small_ancestries() -> None:
    # EAS pools to only 2500 (ANL-only); AMR to 7500
    assert [c.name for c in plan_columns(_COMPOSITION, _SUPERPOPS, min_n=3000)] == [
        "EUR", "AFR", "AMR", "CSA", "MID"
    ]
    assert [c.name for c in plan_columns(_COMPOSITION, _SUPERPOPS, min_n=10000)] == [
        "EUR", "AFR", "CSA", "MID"
    ]


def test_plan_columns_empty_raises() -> None:
    with pytest.raises(ValueError, match="No fine-mapping columns"):
        plan_columns(_COMPOSITION, _SUPERPOPS, min_n=10**9)


def test_fmcolumn_name_is_the_ancestry() -> None:
    col = FMColumn("AFR", 60000, ("anl", "covenant"))
    assert col.name == "AFR" and col.keep_path is None


def test_resolve_binary_prefers_path_then_vendor(tmp_path: Path) -> None:
    # a fake tool only in vendor/bin is found; a truly-absent one raises
    (tmp_path / "vendor" / "bin").mkdir(parents=True)
    tool = tmp_path / "vendor" / "bin" / "FakeTool"
    tool.write_text("#!/bin/sh\n")
    assert _resolve_binary("FakeTool", tmp_path) == str(tool)
    with pytest.raises(FileNotFoundError):
        _resolve_binary("DefinitelyMissingTool123", tmp_path)
