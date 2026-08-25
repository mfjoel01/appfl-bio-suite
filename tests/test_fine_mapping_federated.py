"""Tests for src.fed_fine_mapping — the federated fine-mapping path.

Two layers.

The **unit** layer pins the algebra of ``paper/algo.md`` Part II directly: that
the site aggregates are additive over any partition of individuals (Theorem 1),
that standardizing at the coordinator against pooled moments gives a different
— and correct — answer from standardizing at the sites (§9, the worked example
of A.10), and that the LD panel the coordinator writes by hand has the byte
layout the SuSiEx C++ reader demands.

The **acceptance** layer is the one that matters. Corollary 1 says the
federated fit does not approximate the centralized fit, it *equals* it. So
``test_federated_equals_centralized`` runs both paths — ``src.fine_mapping``
untouched on one side, ``src.fed_fine_mapping`` on the other — over the same
miniature three-site package, and asserts the credible sets and PIPs agree.
Any drift there is a bug in the federated path, not a cost of federating.

The package is built in a temp dir: three ancestry-mixed sites, 60 variants in
one window, a causal variant flanked by two near-perfect LD proxies so a
credible set has to spread its mass rather than collapse onto one variant, and
two phenotype instances (one causal variant, then two). PLINK 1.9, PLINK 2 and
SuSiEx are required for the centralized comparator; the acceptance tests skip
if they are not installed.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import stats

from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import (
    GenoAggregate,
    PhenoAggregate,
    SiteIndex,
    _columns_with_missing,
    _gram,
    build_columns,
    fed_finemap_instance,
    pool_geno,
    pool_pheno,
    pooled_beta_hat,
    pooled_ld,
    pooled_sumstats,
    pooled_tau_sq,
    reference_variants,
    run_fed_fine_mapping,
    site_aggregates,
    write_ld_panel,
)
from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import (
    _LD_SUFFIXES,
    _METRIC_COLS,
    _resolve_binary,
    extract_locus_window,
    finemap_instance,
    materialize_column_keeps,
    plan_columns,
    pooled_ids_path,
    precompute_ld,
    run_susiex,
)
from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import ensure_dir, load_config

REPO_ROOT = Path(__file__).resolve().parents[1]

# --- the miniature package -------------------------------------------------
SITES = {"sitea": {"EUR": 120, "AFR": 60},
         "siteb": {"AFR": 100, "EUR": 40},
         "sitec": {"EUR": 40, "AFR": 40, "CSA": 80}}
SUPERPOPS = ["EUR", "AFR", "CSA"]
CHROM, START, END, M = 1, 1_000_000, 1_030_000, 60
CAUSAL, CAUSAL2 = 31, 12
TWINS = (30, 32)          # near-perfect LD proxies for CAUSAL
# architecture id -> (causal variant indices, target h2)
ARCHS = {"ncsl1_h2-0.005_rg1": ([CAUSAL], 0.20),
         "ncsl2_h2-0.005_rg1": ([CAUSAL, CAUSAL2], 0.34)}
REP = 0
MAF, LEVEL, PVAL_THRESH = 0.01, 0.95, 1e-5


def _snp(j: int) -> str:
    return f"chr{CHROM}:{START + 500 * j}:A:T"


# ---------------------------------------------------------------------------
# Unit layer — the algebra of algo.md Part II
# ---------------------------------------------------------------------------

def _raw_moments(X: np.ndarray, y: np.ndarray) -> dict:
    """The six objects of algo.md §7, computed straight from the data."""
    return {"G": X.T @ X, "c": X.T @ y, "u": X.sum(0),
            "n": len(X), "q": y.sum(), "w": y @ y}


def test_theorem1_aggregates_are_additive_over_sites() -> None:
    # Theorem 1: a matrix product over rows decomposes over ANY partition of
    # those rows, so the pooled sufficient statistics are the sums of the
    # per-site ones. Exactly, not approximately.
    rng = np.random.default_rng(0)
    X = rng.integers(0, 3, size=(90, 7)).astype(float)
    y = rng.normal(size=90)
    whole = _raw_moments(X, y)
    parts = [_raw_moments(X[s], y[s]) for s in (slice(0, 13), slice(13, 51), slice(51, 90))]
    for key in whole:
        assert np.allclose(whole[key], sum(p[key] for p in parts), rtol=0, atol=1e-9)


def test_standardizing_at_the_sites_is_biased() -> None:
    # The worked example of algo.md A.10 / §9. Six people, two SNPs, split 3/3
    # across two sites. SNP j's frequency differs between the sites; SNP k's
    # does not. Pooling the RAW moments and standardizing once gives 0.816;
    # standardizing at each site first gives 0.866 from either site, and hence
    # from any weighted average of them. The gap is the between-site variance
    # that local centring deletes before the coordinator can see it.
    A = np.array([[2., 2.], [2., 1.], [0., 0.]])
    B = np.array([[0., 0.], [0., 1.], [2., 2.]])

    GA, GB = A.T @ A, B.T @ B
    assert np.array_equal(GA, [[8, 6], [6, 5]]) and np.array_equal(GB, [[4, 4], [4, 5]])
    assert np.array_equal(A.sum(0), [4, 3]) and np.array_equal(B.sum(0), [2, 3])

    R, _sd = pooled_ld(GA + GB, A.sum(0) + B.sum(0), n=6)
    assert R[0, 1] == pytest.approx(0.8165, abs=5e-4)

    local = [np.corrcoef(blk, rowvar=False)[0, 1] for blk in (A, B)]
    assert local == pytest.approx([0.866, 0.866], abs=5e-4)
    # and no reweighting of the local answers can recover the pooled one
    for weight in (0.0, 0.5, 1.0):
        blended = weight * local[0] + (1 - weight) * local[1]
        assert abs(blended - R[0, 1]) > 0.04


def test_pooled_ld_equals_the_correlation_of_the_pooled_data() -> None:
    rng = np.random.default_rng(1)
    X = rng.integers(0, 3, size=(200, 6)).astype(float)
    m = _raw_moments(X, np.zeros(200))
    R, sd = pooled_ld(m["G"], m["u"], m["n"])
    np.testing.assert_allclose(R, np.corrcoef(X, rowvar=False), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(sd, X.std(axis=0), rtol=1e-12)


def test_pooled_ld_handles_a_monomorphic_column() -> None:
    # PLINK writes NaN for a zero-variance variant and SuSiEx maps that to 0
    # with a unit diagonal (data.cpp:589-605); the coordinator must do the same
    # rather than emit inf/NaN into the panel.
    X = np.array([[0., 1.], [0., 2.], [0., 1.], [0., 0.]])
    R, sd = pooled_ld(X.T @ X, X.sum(0), n=4)
    assert sd[0] == 0.0
    assert np.array_equal(R, [[1.0, 0.0], [0.0, 1.0]])


def test_pooled_beta_hat_is_eq_92_the_correlation() -> None:
    # algo.md eq (9.2): beta_hat = (c - n xbar ybar) / (n sd_x sd_y), which is
    # exactly corr(dosage, phenotype).
    rng = np.random.default_rng(2)
    X = rng.integers(0, 3, size=(150, 5)).astype(float)
    y = X[:, 2] * 0.4 + rng.normal(size=150)
    m = _raw_moments(X, y)
    got = pooled_beta_hat(np.diag(m["G"]), m["u"], m["c"], m["q"], m["w"], m["n"])
    want = np.array([np.corrcoef(X[:, j], y)[0, 1] for j in range(5)])
    np.testing.assert_allclose(got, want, rtol=1e-12)
    assert pooled_tau_sq(got) == pytest.approx(float(np.max(want ** 2)))


def test_pooled_sumstats_is_the_pooled_ols() -> None:
    # SuSiEx reads BETA/SE, not the correlation, so the closed form the
    # coordinator writes has to be the actual simple-regression fit.
    rng = np.random.default_rng(3)
    X = rng.integers(0, 3, size=(120, 4)).astype(float)
    y = X[:, 1] * 0.5 + rng.normal(size=120)
    m = _raw_moments(X, y)
    variants = pd.DataFrame({"chrom": "1", "snp_id": [_snp(j) for j in range(4)],
                             "bp": [START + j for j in range(4)], "a1": "T", "a2": "A"})
    ss = pooled_sumstats(variants, np.diag(m["G"]), m["u"], m["c"], m["q"], m["w"], m["n"])
    for j in range(4):
        design = np.column_stack([np.ones(120), X[:, j]])
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        resid = y - design @ coef
        se = np.sqrt(resid @ resid / (120 - 2) * np.linalg.inv(design.T @ design)[1, 1])
        assert ss.loc[j, "beta"] == pytest.approx(coef[1], rel=1e-10)
        assert ss.loc[j, "se"] == pytest.approx(se, rel=1e-10)
        assert ss.loc[j, "p"] == pytest.approx(
            2 * stats.t.sf(abs(coef[1] / se), 118), rel=1e-10)


def test_uplink_carries_only_the_six_objects_of_section_7() -> None:
    # A regression guard on scope creep in what crosses the site boundary. The
    # payload is algo.md §7's (G, c, u, n, q, w) and nothing else; the rest is
    # routing labels, the harmonized variant list of Theorem 1 condition (i),
    # and a drop counter. No dosages, no phenotypes, no per-individual anything.
    geno = {f.name for f in GenoAggregate.__dataclass_fields__.values()}
    pheno = {f.name for f in PhenoAggregate.__dataclass_fields__.values()}
    assert geno == {"site", "pop", "locus_id", "variants", "G", "u", "n",
                    "n_incomplete_variants", "snp_ids"}
    assert pheno == {"site", "pop", "instance", "snp_ids", "c", "q", "w", "n"}


def test_gram_is_exact_integers_and_rejects_missing() -> None:
    rng = np.random.default_rng(4)
    X = rng.integers(0, 3, size=(400, 9)).astype(np.float32)
    G = _gram(X)
    assert G.dtype == np.float64 and np.array_equal(G, np.rint(G))
    np.testing.assert_array_equal(G, X.astype(np.float64).T @ X.astype(np.float64))
    X[7, 3] = np.nan
    with pytest.raises(RuntimeError, match="integral"):
        _gram(X)


def test_columns_with_missing_flags_the_right_columns() -> None:
    X = np.zeros((20_000, 4), dtype=np.float32)   # spans several row chunks
    X[15_000, 2] = np.nan
    np.testing.assert_array_equal(_columns_with_missing(X), [False, False, True, False])


def _geno_block(site: str, pop: str, X: np.ndarray, a1: str = "T") -> GenoAggregate:
    variants = pd.DataFrame({"chrom": "1", "snp_id": [_snp(j) for j in range(X.shape[1])],
                             "bp": [START + j for j in range(X.shape[1])],
                             "a1": a1, "a2": "A"})
    return GenoAggregate(site=site, pop=pop, locus_id="L0000", variants=variants,
                         G=X.T @ X, u=X.sum(0), n=len(X))


def test_pool_geno_sums_blocks_and_pool_pheno_follows() -> None:
    rng = np.random.default_rng(5)
    X = rng.integers(0, 3, size=(60, 5)).astype(float)
    y = rng.normal(size=60)
    blocks = [_geno_block("sitea", "EUR", X[:20]), _geno_block("siteb", "EUR", X[20:])]
    pooled = pool_geno(blocks)
    np.testing.assert_allclose(pooled.G, X.T @ X)
    np.testing.assert_allclose(pooled.u, X.sum(0))
    assert pooled.n == 60 and pooled.sites == ("sitea", "siteb")

    pblocks = [PhenoAggregate("sitea", "EUR", "i", pooled.snp_ids, X[:20].T @ y[:20],
                              y[:20].sum(), y[:20] @ y[:20], 20),
               PhenoAggregate("siteb", "EUR", "i", pooled.snp_ids, X[20:].T @ y[20:],
                              y[20:].sum(), y[20:] @ y[20:], 40)]
    pooled_y = pool_pheno(pblocks, pooled.snp_ids)
    np.testing.assert_allclose(pooled_y.c, X.T @ y)
    assert pooled_y.q == pytest.approx(y.sum()) and pooled_y.n == 60


def test_pool_geno_rejects_a_site_that_flipped_a_reference_allele() -> None:
    # Theorem 1 condition (i). A site coding 2-x where the others code x sums in
    # silently and negates every off-diagonal it touches, so it has to be caught
    # at the boundary rather than debugged out of a credible set.
    rng = np.random.default_rng(6)
    X = rng.integers(0, 3, size=(40, 5)).astype(float)
    blocks = [_geno_block("sitea", "EUR", X[:20]), _geno_block("siteb", "EUR", X[20:])]
    blocks[1].variants.loc[2, ["a1", "a2"]] = ["A", "T"]
    with pytest.raises(RuntimeError, match="not harmonized"):
        pool_geno(blocks)


def test_pool_geno_intersects_when_a_site_dropped_a_variant() -> None:
    rng = np.random.default_rng(7)
    X = rng.integers(0, 3, size=(40, 5)).astype(float)
    full = _geno_block("sitea", "EUR", X[:20])
    partial = _geno_block("siteb", "EUR", X[20:])
    keep = [0, 1, 3, 4]                                # site B dropped variant 2
    partial.variants = partial.variants.loc[keep].reset_index(drop=True)
    partial.G = partial.G[np.ix_(keep, keep)]
    partial.u = partial.u[keep]
    partial.__post_init__()
    pooled = pool_geno([full, partial])
    assert list(pooled.snp_ids) == [_snp(j) for j in keep]
    np.testing.assert_allclose(pooled.G, X[:, keep].T @ X[:, keep])


def test_write_ld_panel_matches_the_susiex_reader(tmp_path: Path) -> None:
    # main.cpp:204-222 probes all three suffixes; data.cpp:522 aborts unless the
    # binary is exactly nvar^2 * 4 bytes against the _ref.bim line count, and
    # data.cpp:173-183 aborts unless the .frq repeats each bim line's alleles.
    rng = np.random.default_rng(8)
    m = 7
    R = np.corrcoef(rng.normal(size=(50, m)), rowvar=False)
    variants = pd.DataFrame({"chrom": "1", "snp_id": [_snp(j) for j in range(m)],
                             "bp": [START + 500 * j for j in range(m)],
                             "a1": "T", "a2": "A"})
    freq = rng.uniform(0.05, 0.95, m)
    prefix = write_ld_panel(tmp_path / "EUR_ld", variants, R, freq, n=123)

    for suffix in _LD_SUFFIXES:
        assert Path(f"{prefix}{suffix}").exists()
    assert Path(f"{prefix}.ld.bin").stat().st_size == m * m * 4
    back = np.fromfile(f"{prefix}.ld.bin", dtype=np.float32).reshape(m, m)
    np.testing.assert_allclose(back, R.astype(np.float32), rtol=0, atol=0)

    bim = pd.read_csv(f"{prefix}_ref.bim", sep=r"\s+", header=None)
    frq = pd.read_csv(f"{prefix}_frq.frq", sep=r"\s+")
    assert len(bim) == m and len(frq) == m
    assert list(bim[1]) == list(variants["snp_id"])          # order preserved
    assert list(bim[3]) == list(variants["bp"])
    assert list(frq["SNP"]) == list(bim[1])
    assert list(frq["A1"]) == list(bim[4]) and list(frq["A2"]) == list(bim[5])
    np.testing.assert_allclose(frq["MAF"], freq, rtol=1e-3)
    assert set(frq["NCHROBS"]) == {246}


def test_write_ld_panel_rejects_a_mismatched_matrix(tmp_path: Path) -> None:
    variants = pd.DataFrame({"chrom": "1", "snp_id": [_snp(0), _snp(1)],
                             "bp": [START, START + 1], "a1": "T", "a2": "A"})
    with pytest.raises(ValueError, match=r"expected \(2, 2\)"):
        write_ld_panel(tmp_path / "x", variants, np.eye(3), np.full(2, 0.3), n=10)


# ---------------------------------------------------------------------------
# The miniature three-site package
# ---------------------------------------------------------------------------

def _write_bed(prefix: Path, X: np.ndarray) -> None:
    """Variant-major PLINK 1 .bed: dosage 0/1/2 -> bit codes 11/10/00."""
    code = np.array([0b11, 0b10, 0b00], dtype=np.uint8)
    n, m = X.shape
    pad = (-n) % 4
    with open(f"{prefix}.bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        for j in range(m):
            col = np.concatenate([code[X[:, j]], np.zeros(pad, np.uint8)]).reshape(-1, 4)
            fh.write((col[:, 0] | col[:, 1] << 2 | col[:, 2] << 4
                      | col[:, 3] << 6).astype(np.uint8).tobytes())


def _write_fam(prefix: Path, iids) -> None:
    with open(f"{prefix}.fam", "w") as fh:
        for i in iids:
            fh.write(f"{i}\t{i}\t0\t0\t0\t-9\n")


def _haplo_pop(rng, n: int, freq: np.ndarray, rho: float) -> np.ndarray:
    """Dosages for ``n`` people with the given allele frequencies and AR(1) LD."""
    z = np.empty((2 * n, len(freq)))
    z[:, 0] = rng.normal(size=2 * n)
    for j in range(1, len(freq)):
        z[:, j] = rho * z[:, j - 1] + np.sqrt(1 - rho ** 2) * rng.normal(size=2 * n)
    hap = (z < stats.norm.ppf(freq)).astype(np.int8)
    return (hap[0::2] + hap[1::2]).astype(np.int8)


def _build_package(root: Path, seed: int = 20260820) -> Path:
    rng = np.random.default_rng(seed)
    bp = START + 500 * np.arange(M)
    bim = pd.DataFrame({"chrom": CHROM, "snp_id": [_snp(j) for j in range(M)],
                        "cm": 0, "bp": bp, "a1": "T", "a2": "A"})

    # per-ancestry frequencies; a fifth of the variants are rare enough that the
    # pooled MAF filter has something to do and the two paths have to agree on it
    freqs, rhos = {}, {"EUR": 0.86, "AFR": 0.70, "CSA": 0.80}
    for pop in SUPERPOPS:
        f = rng.uniform(0.06, 0.45, M)
        rare = rng.random(M) < 0.2
        f[rare] = rng.uniform(0.001, 0.03, rare.sum())
        f[[CAUSAL, CAUSAL2]] = rng.uniform(0.2, 0.4, 2)
        freqs[pop] = f

    rows, parts = [], []
    for site, comp in SITES.items():
        for pop, n in comp.items():
            parts.append(_haplo_pop(rng, n, freqs[pop], rhos[pop]))
            rows += [(site, pop)] * n
    parts.append(_haplo_pop(rng, 20, freqs["EUR"], rhos["EUR"]))   # not enrolled
    rows += [(None, "EUR")] * 20

    X = np.vstack(parts)
    for twin in TWINS:                       # near-perfect LD proxies for CAUSAL
        X[:, twin] = np.where(rng.random(len(X)) < 0.03, X[:, twin], X[:, CAUSAL])
    who = pd.DataFrame(rows, columns=["site", "superpopulation"])
    who["IID"] = [f"syn{i:04d}" for i in range(len(who))]
    order = rng.permutation(len(who))        # site filesets are not contiguous
    X, who = X[order], who.iloc[order].reset_index(drop=True)

    def _pheno(idx, h2):
        g = np.zeros(len(X))
        for j in idx:
            z = X[:, j].astype(float)
            g += (z - z.mean()) / z.std()
        g /= g.std()
        return np.sqrt(h2) * g + np.sqrt(1 - h2) * rng.normal(size=len(X))

    phenos = {arch: _pheno(idx, h2) for arch, (idx, h2) in ARCHS.items()}

    hap = ensure_dir(root / "data" / "raw" / "hapnest")
    _write_bed(hap / f"chr{CHROM}", X)
    bim.to_csv(hap / f"chr{CHROM}.bim", sep="\t", header=False, index=False)
    _write_fam(hap / f"chr{CHROM}", who["IID"])
    who.assign(FID=who["IID"])[["FID", "IID", "superpopulation"]].to_csv(
        hap / "population_manifest.tsv", sep="\t", index=False)

    for site in SITES:
        sel = rng.permutation(np.flatnonzero((who["site"] == site).to_numpy()))
        d = ensure_dir(root / "data" / "processed" / site)
        _write_bed(d / f"{site}_chr{CHROM}", X[sel])
        bim.to_csv(d / f"{site}_chr{CHROM}.bim", sep="\t", header=False, index=False)
        _write_fam(d / f"{site}_chr{CHROM}", who["IID"].iloc[sel])
        sub = who.iloc[sel]
        sub.assign(FID=sub["IID"])[["FID", "IID", "superpopulation"]].to_csv(
            d / f"{site}_manifest.tsv", sep="\t", index=False)
        ph = ensure_dir(root / "data" / "ground_truth" / "phenotypes" / site)
        for arch, y in phenos.items():
            pd.DataFrame({"FID": sub["IID"].values, "IID": sub["IID"].values,
                          "y": y[sel]}).to_csv(
                ph / f"L0000_{arch}_rep{REP}.pheno", sep="\t", index=False)

    loci = ensure_dir(root / "data" / "loci")
    pd.DataFrame([{"locus_id": "L0000", "window_id": "w0", "chrom": CHROM,
                   "start_bp": START, "end_bp": END, "divergence_score": 1.0,
                   "stratum": "high", "n_variants_in_window": M, "n_tag_snps": M}]
                 ).to_csv(loci / "selected_loci.tsv", sep="\t", index=False)
    pd.DataFrame([{"locus_id": "L0000", "architecture_id": arch, "replicate": REP,
                   "causal_snp_ids": ",".join(_snp(j) for j in idx)}
                  for arch, (idx, _h2) in ARCHS.items()]).to_csv(
        root / "data" / "ground_truth" / "causal_manifest.tsv", sep="\t", index=False)
    ensure_dir(root / "reports")
    ensure_dir(root / "logs")
    # the driver resolves its binaries against the config's repo root
    vendor = REPO_ROOT / "vendor" / "bin"
    if vendor.exists():
        ensure_dir(root / "vendor")
        (root / "vendor" / "bin").symlink_to(vendor)

    cfg = {
        "paths": {"hapnest_dir": "data/raw/hapnest",
                  "population_manifest": "data/raw/hapnest/population_manifest.tsv",
                  "processed_dir": "data/processed", "loci_dir": "data/loci",
                  "ground_truth_dir": "data/ground_truth", "reports_dir": "reports",
                  "logs_dir": "logs"},
        "master_seed": 1, "chromosome": CHROM,
        "tools": {"plink": "plink", "plink2": "plink2", "king": "king"},
        "sites": {s: {"n": sum(c.values()), "composition": c,
                      "dominant": max(c, key=c.get)} for s, c in SITES.items()},
        "superpopulations": SUPERPOPS,
        "locus_selection": {"window_size_bp": END - START, "step_size_bp": 1000,
                            "n_loci": 3, "strata": {"low": 1, "medium": 1, "high": 1},
                            "maf_filter": 0.01, "ld_tag_snp_target": 10,
                            "ld_r2_prune_threshold": 0.99,
                            "ld_prune_window_variants": 10, "ld_prune_step_variants": 2},
        "architecture": {"ncsl": [1, 2], "h2": [0.005], "rg": [1.0],
                         "factorial_mode": "minimal", "replicates": 1,
                         "ancestry_specific_causal": {
                             "enabled": False, "min_per_stratum": 0,
                             "common_maf_threshold": 0.05, "rare_maf_threshold": 0.01}},
        "phenotype": {"model": "linear_additive", "h2_tolerance_relative": 0.2},
        "qc": {"kinship_threshold": 0.05, "pca_n_components": 2, "pca_reference": "x",
               "maf_tolerance_sd": 3.0, "ld_decay_max_dist_kb": 10,
               "ld_decay_n_loci_sample": 1},
        "fine_mapping": {"min_gwas_n": 10},
    }
    ensure_dir(root / "config")
    path = root / "config" / "simulation_config.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


@pytest.fixture(scope="module")
def package(tmp_path_factory: pytest.TempPathFactory):
    root = tmp_path_factory.mktemp("fedfm")
    cfg_path = _build_package(root)
    return root, load_config(cfg_path, repo_root=root)


@pytest.fixture(scope="module")
def binaries():
    try:
        return {name: _resolve_binary(name, REPO_ROOT)
                for name in ("SuSiEx", "plink", "plink2")}
    except FileNotFoundError as exc:
        pytest.skip(f"fine-mapping binaries unavailable: {exc}")


@pytest.fixture(scope="module")
def locus(package) -> pd.Series:
    root, _cfg = package
    return pd.read_csv(root / "data" / "loci" / "selected_loci.tsv", sep="\t").iloc[0]


def _truth(arch: str) -> list[str]:
    return [_snp(j) for j in ARCHS[arch][0]]


@pytest.fixture(scope="module")
def centralized(package, binaries, locus):
    """Run src.fine_mapping, untouched, on every instance of the package."""
    root, cfg = package
    work = ensure_dir(root / "central")
    refs = ensure_dir(work / "refs")
    columns = plan_columns({s: dict(sc.composition) for s, sc in cfg.sites.items()},
                           cfg.superpopulations, cfg.fine_mapping.min_gwas_n)
    materialize_column_keeps(columns, cfg, work / "keeps")
    enrolled = pooled_ids_path(cfg, work / "keeps")
    window = extract_locus_window(cfg, locus, enrolled, refs, binaries["plink2"])
    ld = precompute_ld(locus, columns, window, refs, binaries["plink"], MAF)
    out = {}
    for arch in ARCHS:
        row = finemap_instance(
            cfg, locus, arch, REP, _truth(arch), columns, window, ld, work,
            binaries["SuSiEx"], binaries["plink"], binaries["plink2"],
            LEVEL, PVAL_THRESH, True)
        out[arch] = (row, work / f"L0000_{arch}_rep{REP}")
    return out, columns, ld


@pytest.fixture(scope="module")
def federated(package, binaries, locus):
    """Run src.fed_fine_mapping on the same instances: sites emit, coordinator fits."""
    root, cfg = package
    work = ensure_dir(root / "fed")
    pops = list(SUPERPOPS)
    instances = [f"L0000_{arch}_rep{REP}" for arch in ARCHS]

    reference = reference_variants(cfg)
    geno, pheno = [], []
    for site in cfg.sites:
        index = SiteIndex.load(cfg, site, reference)
        g, p = site_aggregates(cfg, site, locus, pops, instances, index=index)
        geno += g
        pheno += p
    columns = build_columns(geno, pops, ensure_dir(work / "ld"), MAF)

    by_inst: dict[str, dict[str, list]] = {}
    for block in pheno:
        by_inst.setdefault(block.instance, {}).setdefault(block.pop, []).append(block)

    out = {}
    for arch in ARCHS:
        inst = f"L0000_{arch}_rep{REP}"
        moments = {c.pop: pool_pheno(by_inst[inst][c.pop], c.variants["snp_id"].to_numpy())
                   for c in columns}
        row = fed_finemap_instance(
            locus, arch, REP, _truth(arch), columns, moments, work,
            binaries["SuSiEx"], binaries["plink"], LEVEL, PVAL_THRESH, True)
        out[arch] = (row, work / inst)
    return out, columns, geno, pheno


@pytest.fixture(scope="module")
def hybrid(package, binaries, locus, centralized, federated):
    """SuSiEx on the coordinator's sumstats + PLINK's LD panels.

    Isolates the two halves of the federated reconstruction: same summary
    statistics as the federated run, same LD matrices as the centralized one.
    Whatever this agrees with, it agrees with exactly.
    """
    root, _cfg = package
    _crows, _ccols, cld = centralized
    _frows, fcols, _g, _p = federated
    out = {}
    for arch in ARCHS:
        inst = f"L0000_{arch}_rep{REP}"
        target = ensure_dir(root / "hybrid" / inst)
        run_susiex([root / "fed" / inst / f"{c.pop}.sumstats" for c in fcols],
                   [cld[c.pop] for c in fcols], [c.n for c in fcols],
                   int(locus["chrom"]), int(locus["start_bp"]), int(locus["end_bp"]),
                   target, "cs", binaries["SuSiEx"], binaries["plink"],
                   LEVEL, PVAL_THRESH)
        out[arch] = target
    return out


# ---------------------------------------------------------------------------
# Acceptance — algo.md Corollary 1
# ---------------------------------------------------------------------------
#
# Tolerance, and why it is what it is. The residual between the two paths
# decomposes cleanly, and both halves are pinned separately below:
#
#   * the summary statistics half is EXACT. Feed the coordinator's closed-form
#     sumstats to SuSiEx alongside PLINK's own LD panels and the PIPs come back
#     bit-identical to the centralized run — 0.0 difference, asserted with no
#     tolerance at all in ``test_federated_sumstats_reproduce_the_fit_exactly``.
#   * the LD half is float32 on both sides — PLINK's ``--r square bin4`` writes
#     single precision and so does :func:`write_ld_panel` — and the two agree to
#     within 2 units in the last place, asserted in
#     ``test_federated_ld_panel_matches_the_plink_one``.
#
# So the only thing that can move a PIP is a couple of float32 ULP in R, and
# what it moves is the last of the six significant digits SuSiEx prints into
# its ``.snp`` file (observed worst case: 0.0729835 against 0.0729836). The
# tolerance below is that print resolution. Tightening it past what the file
# format can express would be asserting on noise; loosening it would let a real
# discrepancy through, since anything larger than a last-digit wobble is a bug.
PIP_RTOL, PIP_ATOL = 2e-5, 1e-12


def _credible_sets(out_dir: Path) -> dict[int, frozenset[str]]:
    cs = pd.read_csv(out_dir / "cs.cs", sep="\t", comment="#")
    return {int(cid): frozenset(g["SNP"].astype(str))
            for cid, g in cs.groupby("CS_ID")}


def _pips(out_dir: Path) -> pd.DataFrame:
    snp = pd.read_csv(out_dir / "cs.snp", sep="\t")
    cols = [c for c in snp.columns if c.startswith("PIP(")]
    return snp.set_index(snp["SNP"].astype(str))[cols].sort_index()


@pytest.mark.parametrize("arch", list(ARCHS))
def test_federated_equals_centralized(centralized, federated, arch) -> None:
    """THE acceptance test: algo.md Corollary 1, credible sets and PIPs."""
    (crows, _ccols, _cld), (frows, _fcols, _g, _p) = centralized, federated
    crow, cdir = crows[arch]
    frow, fdir = frows[arch]
    assert crow["error"] == "" and frow["error"] == "", (crow["error"], frow["error"])
    # a fit that found nothing would make the comparison vacuous
    assert crow["n_credible_sets"] >= 1 and crow["converged"]

    cs_c, cs_f = _credible_sets(cdir), _credible_sets(fdir)
    assert cs_f == cs_c, f"credible sets differ: {cs_f} vs {cs_c}"

    pip_c, pip_f = _pips(cdir), _pips(fdir)
    assert list(pip_f.index) == list(pip_c.index)
    assert list(pip_f.columns) == list(pip_c.columns)
    np.testing.assert_allclose(pip_f.to_numpy(), pip_c.to_numpy(),
                               rtol=PIP_RTOL, atol=PIP_ATOL)


@pytest.mark.parametrize("arch", list(ARCHS))
def test_federated_sumstats_reproduce_the_fit_exactly(centralized, hybrid, arch) -> None:
    """The GWAS half of Corollary 1, with no tolerance at all.

    The coordinator never ran a regression on genotypes — it inverted the summed
    moments in closed form. Hand SuSiEx those summary statistics with PLINK's
    own LD panels and every PIP comes back bit-identical to the centralized
    fit. Whatever separates the two full paths is therefore the LD panel and
    nothing else.
    """
    _crow, cdir = centralized[0][arch]
    pip_c, pip_h = _pips(cdir), _pips(hybrid[arch])
    assert list(pip_h.index) == list(pip_c.index)
    np.testing.assert_array_equal(pip_h.to_numpy(), pip_c.to_numpy())
    assert _credible_sets(hybrid[arch]) == _credible_sets(cdir)


@pytest.mark.parametrize("arch", list(ARCHS))
def test_federated_metrics_row_equals_centralized(centralized, federated, arch) -> None:
    """The scored result row — what actually lands in the results table."""
    crow, _ = centralized[0][arch]
    frow, _ = federated[0][arch]
    for key in ("n_causal", "n_credible_sets", "cs_sizes_json", "total_cs_snps",
                "n_causal_captured", "any_causal_captured", "best_cs_size",
                "converged"):
        assert frow[key] == crow[key], key
    for key in ("causal_pip_max", "causal_pip_mean", "top_pip",
                "mean_cs_purity", "min_cs_purity"):
        assert frow[key] == pytest.approx(crow[key], rel=PIP_RTOL, abs=PIP_ATOL), key


def test_federated_ld_panel_matches_the_plink_one(centralized, federated) -> None:
    """The panels the coordinator wrote by hand vs the ones PLINK computed."""
    _crows, _ccols, cld = centralized
    _frows, fcols, _g, _p = federated
    for col in fcols:
        cbim = pd.read_csv(f"{cld[col.pop]}_ref.bim", sep=r"\s+", header=None)
        fbim = pd.read_csv(f"{col.ld_prefix}_ref.bim", sep=r"\s+", header=None)
        # same variants, same order, same alleles: the pooled MAF filter the
        # coordinator applies to u/2n reproduces `plink --maf` exactly
        assert fbim.equals(cbim), col.pop
        m = len(fbim)
        c_ld = np.fromfile(f"{cld[col.pop]}.ld.bin", dtype=np.float32).reshape(m, m)
        f_ld = np.fromfile(f"{col.ld_prefix}.ld.bin", dtype=np.float32).reshape(m, m)
        # float32 storage on both sides (`plink --r square bin4` vs
        # write_ld_panel), so agreement is measured in units in the last place
        ulp = np.abs(f_ld.view(np.int32).astype(np.int64)
                     - c_ld.view(np.int32).astype(np.int64))
        assert ulp.max() <= 2, f"{col.pop}: LD differs by {ulp.max()} ULP"


@pytest.mark.parametrize("arch", list(ARCHS))
def test_federated_sumstats_match_plink2_glm(centralized, federated, arch) -> None:
    """The GWAS the coordinator did in closed form vs the one plink2 ran.

    ``plink2 --glm`` reports the minor allele as A1, so for a variant whose
    minor allele is the .bim A2 its BETA is signed the other way round from the
    coordinator's, which always counts the .bim A1. SuSiEx reconciles that
    itself (``data.hpp:349-360`` negates a reversed effect), so the comparison
    orients both to A1 before checking.
    """
    _crows, _ccols, _cld = centralized
    _frows, fcols, _g, _p = federated
    inst = f"L0000_{arch}_rep{REP}"
    root = Path(_crows[arch][1]).parent.parent
    for col in fcols:
        c = pd.read_csv(root / "central" / inst / f"{col.pop}.sumstats", sep="\t")
        f = pd.read_csv(root / "fed" / inst / f"{col.pop}.sumstats", sep="\t")
        # the federated file covers exactly the LD panel; the centralized one is
        # the whole window and SuSiEx intersects. Compare where both speak.
        m = c.merge(f, on="snp", suffixes=("_c", "_f"))
        assert len(m) == len(f) == len(col.variants)
        sign = np.where(m["A1_c"] == m["A1_f"], 1.0, -1.0)
        np.testing.assert_allclose(sign * m["beta_c"], m["beta_f"], rtol=1e-5, atol=1e-12)
        np.testing.assert_allclose(m["se_c"], m["se_f"], rtol=1e-5, atol=1e-12)
        np.testing.assert_allclose(m["p_c"], m["p_f"], rtol=1e-5, atol=1e-300)


def test_site_aggregates_pool_to_the_planned_cohorts(package, federated) -> None:
    """Each ancestry column is that ancestry's people, wherever they live."""
    _root, cfg = package
    _rows, columns, geno, _pheno = federated
    planned = {c.name: c for c in plan_columns(
        {s: dict(sc.composition) for s, sc in cfg.sites.items()},
        cfg.superpopulations, cfg.fine_mapping.min_gwas_n)}
    assert [c.pop for c in columns] == SUPERPOPS
    for col in columns:
        assert col.n == planned[col.pop].n
        assert col.sites == planned[col.pop].sites
    # 3 sites x their ancestries = 7 blocks here (12 in the production package)
    assert len(geno) == sum(len(c) for c in SITES.values())
    for block in geno:
        assert block.n == SITES[block.site][block.pop]
        assert block.n_incomplete_variants == 0


def test_driver_runs_end_to_end_and_writes_the_result_table(package, binaries) -> None:
    """``run_fed_fine_mapping`` itself: site stage, coordinator, workers, rollups.

    Two workers, so the per-locus columns and pooled moments have to survive the
    trip through loky — and the results table has to come out with the schema
    the centralized one has, since the whole point is that the two are directly
    comparable.
    """
    root, cfg = package
    shutil.rmtree(root / "reports" / "fed_fine_mapping", ignore_errors=True)
    run_fed_fine_mapping(cfg, n_workers=2, level=LEVEL, pval_thresh=PVAL_THRESH, maf=MAF)

    out = root / "reports" / "fed_fine_mapping"
    df = pd.read_csv(out / "fed_fm_results.tsv", sep="\t")
    assert len(df) == len(ARCHS)
    assert set(df["architecture_id"]) == set(ARCHS)
    assert df["error"].isna().all() or (df["error"].fillna("") == "").all()
    assert df["converged"].all() and df["any_causal_captured"].all()
    assert set(_METRIC_COLS).issubset(df.columns)
    assert {f"min_p_{p}" for p in SUPERPOPS}.issubset(df.columns)
    for name in ("by_architecture", "by_stratum_rg"):
        assert (out / f"fed_fm_rollup_{name}.tsv").exists()
    # the coordinator's scratch is cleaned up when not asked to keep it
    assert not list((out / "work").glob("L0000_*"))


def test_a_site_never_reads_another_sites_data(package, locus, monkeypatch) -> None:
    """The site stage touches only its own directory tree.

    Not a proof of privacy — it is a guard that the block a site emits is a
    function of that site's own genotypes and phenotypes and nothing else, which
    is the property the whole federation claim rests on.
    """
    _root, cfg = package
    opened: list[str] = []
    real_open = Path.open

    def spy(self, *args, **kwargs):
        opened.append(str(self))
        return real_open(self, *args, **kwargs)

    monkeypatch.setattr(Path, "open", spy)
    real_read_csv = pd.read_csv

    def spy_csv(path, *args, **kwargs):
        opened.append(str(path))
        return real_read_csv(path, *args, **kwargs)

    monkeypatch.setattr(pd, "read_csv", spy_csv)
    site_aggregates(cfg, "sitea", locus, list(SUPERPOPS),
                    [f"L0000_{next(iter(ARCHS))}_rep{REP}"])
    monkeypatch.undo()
    for other in ("siteb", "sitec"):
        assert not [p for p in opened if other in p], f"read {other}'s data"
    assert not [p for p in opened if "hapnest" in p]


def _copy_package(root: Path, dest: Path):
    shutil.copytree(root / "data", dest / "data")
    shutil.copytree(root / "config", dest / "config")
    ensure_dir(dest / "reports")
    ensure_dir(dest / "logs")
    return load_config(dest / "config" / "simulation_config.yaml", repo_root=dest)


def _flip_bed_variant(bed: Path, j: int, n: int) -> None:
    """Swap which allele a .bed variant counts: hom-A1 <-> hom-A2 (00 <-> 11)."""
    stride = (n + 3) // 4
    with open(bed, "r+b") as fh:
        fh.seek(3 + j * stride)
        block = np.frombuffer(fh.read(stride), dtype=np.uint8)
        codes = np.empty(stride * 4, dtype=np.uint8)
        for k in range(4):
            codes[k::4] = (block >> (2 * k)) & 0x03
        codes = np.where(codes == 0b00, 0b11, np.where(codes == 0b11, 0b00, codes))
        packed = codes.reshape(-1, 4)
        fh.seek(3 + j * stride)
        fh.write((packed[:, 0] | packed[:, 1] << 2 | packed[:, 2] << 4
                  | packed[:, 3] << 6).astype(np.uint8).tobytes())


def test_a_site_with_flipped_reference_alleles_is_harmonized(package, locus, tmp_path):
    """Theorem 1 condition (i), on the defect the real package actually has.

    The production per-site filesets were cut with ``plink --keep --make-bed``
    and no ``--keep-allele-order``, so PLINK 1.9 chose A1 by each site's own
    minor allele and ~6% of chr1 variants ended up coded oppositely at ANL and
    Covenant. Here one site is flipped the same way — .bim swapped and .bed
    recoded — and the block it emits must come out bit-identical to the
    unflipped one, because the site recodes against the reference list before
    forming a single moment.
    """
    root, cfg = package
    cfg2 = _copy_package(root, tmp_path / "flipped")
    flipped = [3, 17, CAUSAL]
    site_dir = tmp_path / "flipped" / "data" / "processed" / "sitea"
    n_site = sum(SITES["sitea"].values())
    for j in flipped:
        _flip_bed_variant(site_dir / f"sitea_chr{CHROM}.bed", j, n_site)
    bim_path = site_dir / f"sitea_chr{CHROM}.bim"
    bim = pd.read_csv(bim_path, sep="	", header=None)
    bim.loc[flipped, [4, 5]] = bim.loc[flipped, [5, 4]].to_numpy()
    bim.to_csv(bim_path, sep="	", header=False, index=False)

    reference = reference_variants(cfg)
    index = SiteIndex.load(cfg2, "sitea", reference)
    assert index.n_flipped == len(flipped)
    assert list(index.bim.loc[flipped, "a1"]) == ["T"] * len(flipped)   # canonical again

    inst = f"L0000_{next(iter(ARCHS))}_rep{REP}"
    args = (list(SUPERPOPS), [inst])
    got, got_pheno = site_aggregates(cfg2, "sitea", locus, *args, index=index)
    want, want_pheno = site_aggregates(
        cfg, "sitea", locus, *args, index=SiteIndex.load(cfg, "sitea", reference))
    assert [b.pop for b in got] == [b.pop for b in want]
    for g, w in zip(got, want):
        np.testing.assert_array_equal(g.G, w.G)
        np.testing.assert_array_equal(g.u, w.u)
        assert g.variants.equals(w.variants)
    for g, w in zip(got_pheno, want_pheno):
        np.testing.assert_allclose(g.c, w.c, rtol=1e-12)

    # ... and without harmonizing, the flip survives the sum and corrupts it
    naive = site_aggregates(cfg2, "sitea", locus, *args,
                            index=SiteIndex.load(cfg2, "sitea"))[0]
    assert not np.array_equal(naive[0].G, want[0].G)


def test_partially_observed_variants_are_dropped_and_counted(package, locus, tmp_path):
    """Theorem 1 condition (ii): a variant a site cannot observe completely.

    ``G`` and ``u`` are only exact integers over a cohort with no missing calls,
    and neither ``plink2 --glm``'s per-variant nor ``plink --r``'s pairwise
    complete-case handling is expressible in one ``(G, u, n)`` triple. So the
    site drops such a variant and says so, and the coordinator intersects.
    """
    root, cfg = package
    copy = tmp_path / "pkg"
    cfg2 = _copy_package(root, copy)

    # blank out one genotype of variant 5 at sitea: bit code 01 = missing
    bed = copy / "data" / "processed" / "sitea" / f"sitea_chr{CHROM}.bed"
    n_site = sum(SITES["sitea"].values())
    offset = 3 + 5 * ((n_site + 3) // 4)
    with open(bed, "r+b") as fh:
        fh.seek(offset)
        byte = fh.read(1)[0]
        fh.seek(offset)
        fh.write(bytes([(byte & 0b11111100) | 0b01]))

    inst = f"L0000_{next(iter(ARCHS))}_rep{REP}"
    geno, _pheno = site_aggregates(cfg2, "sitea", locus, list(SUPERPOPS), [inst])
    hit = [b for b in geno if b.n_incomplete_variants]
    assert len(hit) == 1, "exactly one ancestry block holds that individual"
    assert len(hit[0].snp_ids) == M - 1
    assert _snp(5) not in set(hit[0].snp_ids)
    assert hit[0].n == SITES["sitea"][hit[0].pop]      # the cohort is not reduced

    clean, _ = site_aggregates(cfg, "sitea", locus, list(SUPERPOPS), [inst])
    same = [b for b in clean if b.pop == hit[0].pop][0]
    assert same.n_incomplete_variants == 0 and len(same.snp_ids) == M
    # and the coordinator falls back to what every site can speak about
    others = [b for b in clean if b.pop == hit[0].pop and b.site != "sitea"]
    pooled = pool_geno([hit[0], *others]) if others else pool_geno([hit[0]])
    assert _snp(5) not in set(pooled.snp_ids)
