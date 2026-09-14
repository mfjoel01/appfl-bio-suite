"""Federated SuSiEx cross-ancestry fine-mapping — the site -> coordinator path.

The twin of :mod:`src.fine_mapping`. Same estimator, same SuSiEx binary, same
six ancestry columns — but no individual-level data ever leaves a site. This
module implements the protocol derived in ``paper/algo.md`` Part II (§7 and §9;
Appendix A.8-A.11 is the same thing in prose), in two stages:

**Site stage** (:func:`site_aggregates`). Each of the three sites reads *only*
its own ``data/processed/{site}/{site}_chr{N}.bed`` for the locus window and,
for each superpopulation it holds, emits the raw second-moment aggregates of
algo.md §7 --- unstandardized, on the raw 0/1/2 dosage scale::

    G = X'X   (M x M)      c = X'y   (M)      u = 1'X   (M)
    n                       q = 1'y            w = y'y

Twelve ``(site, ancestry)`` blocks in all: ANL holds five ancestries, Covenant
three, MBZUAI four. ``G`` and ``u`` are exact integers (see
:class:`GenoAggregate`). Nothing else crosses the site boundary — in particular
no genotype row, no phenotype value, and no locally standardized LD matrix.

Only ``c``, ``q`` and ``w`` depend on the phenotype, so the O(M^2) half of the
uplink (:class:`GenoAggregate`) is emitted once per (site, ancestry, locus) and
the O(M) half (:class:`PhenoAggregate`) once per instance. That split is an
economy of transmission, not of information: together the two are exactly the
six objects of §7.

**Coordinator stage** (:func:`build_columns`, :func:`fed_finemap_instance`).
The coordinator sums each ancestry's blocks across sites (Theorem 1: matrix
products decompose over any partition of rows) and only *then* standardizes,
once, against the pooled moments — eq (9.1)-(9.2). This ordering is the one
real trap of the whole design: standardizing at the sites and combining the
results is biased, because each site would centre against its own column means
and delete the between-site frequency variation before the coordinator can see
it (algo.md §9; the worked example in A.10 gives 0.866 where the truth is
0.816 — :func:`~tests.test_fed_fine_mapping` pins that number).

The coordinator holds no genotypes, so unlike
:func:`src.fine_mapping.precompute_ld` it cannot reach for ``plink --r`` to
build SuSiEx's LD reference panel. It writes the three panel files itself, in
the byte layout the SuSiEx C++ reader expects (:func:`write_ld_panel`;
``main.cpp:204-222``, ``data.cpp:126-188,505-610``). SuSiEx is then invoked
exactly as the centralized path invokes it, through
:func:`src.fine_mapping.run_susiex`.

**Exactness.** By algo.md Corollary 1 this is not an approximation of the
centralized fit — it *is* the centralized fit, up to floating-point summation
order, provided the harmonization conditions of Theorem 1 hold (a common
ordered allele-harmonized variant list; a common complete-case cohort per
ancestry; identical phenotype/QC handling). ``tests/test_fed_fine_mapping.py``
runs both paths on one (locus, architecture, replicate) and asserts the
credible sets and PIPs agree. A mismatch there is a bug in *this* module, never
a cost of federating.

**Allele harmonization.** Theorem 1 condition (i) wants one ordered,
allele-harmonized variant list, and in this package it is not free: the
per-site filesets were cut without ``--keep-allele-order``, so PLINK 1.9 set A1
to each site's own minor allele and ~6% of chr1 variants are coded oppositely
at ANL and Covenant. The centralized path never notices — it cuts its window
straight from the HAPNEST fileset — but a flipped site contributes ``2 - x``
where the others contribute ``x``, which sums in silently and negates every
off-diagonal that variant touches (algo.md A.10, trap 2). So each site recodes
against :func:`reference_variants` before forming a single moment, and the
coordinator re-checks (:func:`pool_geno`) rather than trusting it.

**Missing genotypes.** Theorem 1 condition (ii) needs one complete-case cohort
per ancestry, or missing-data handling identical to the comparator's. The
comparator's is per-variant (``plink2 --glm`` regresses on that variant's
complete cases) and pairwise (``plink --r`` correlates on pairwise-complete
observations) — neither is expressible in a single ``(G, u, n)`` triple. So a
site drops any window variant not fully observed in the block, keeping every
retained variant's moments exact over the block's full cohort, and reports the
count. On HAPNEST chr1 this costs about four variants in two thousand (the
missing-call rate is ~2e-7, concentrated in a handful of variants); the
coordinator logs the drop, and where it is nonzero the two paths fine-map
slightly different variant sets and Corollary 1's equality no longer bites.

**Privacy.** Out of scope by ``design.md`` §11.5: the genotypes are synthetic,
so the aggregates ship in the clear — no secure aggregation, no DP. What that
section says a real deployment would need is not implemented here.

Sharding mirrors the centralized stage: ``--n-shards N`` partitions loci by
``index % N``; ``--merge`` reduces the parts into
``reports/fed_fine_mapping/fed_fm_results.tsv``.
"""

from __future__ import annotations

import argparse
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

from .fine_mapping import (
    _METRIC_COLS,
    _SUMSTATS_COLS,
    _resolve_binary,
    parse_susiex,
    plan_columns,
    run_susiex,
)
from .phenotype_sim import build_architecture_grid
from .utils import (
    SimulationConfig,
    ensure_dir,
    get_logger,
    load_config,
    read_bed_variants,
    read_bim,
    read_fam,
    read_site_manifest,
    setup_logging,
)

# ``X'X`` entries are sums of products of dosages in {0,1,2}, so they are
# integers bounded by 4n. float32 holds every integer below 2^24 exactly, and
# every partial sum a BLAS sgemm forms is itself such an integer, so the
# single-precision Gram is exact — and half the memory — while 4n stays under
# that ceiling. Above it we pay for float64.
_EXACT_GRAM_MAX_N = (1 << 24) // 4

# Rows read per pass when accumulating X'y in float64 out of a float32 dosage
# block: big enough to keep BLAS busy, small enough that the upcast temporary
# stays ~100 MB at locus-scale M.
_ROW_CHUNK = 8192


# ---------------------------------------------------------------------------
# Stage 1 — what a site emits
# ---------------------------------------------------------------------------

@dataclass
class GenoAggregate:
    """Phenotype-independent half of one ``(site, ancestry, locus)`` block.

    ``G = X'X`` and ``u = 1'X`` over the raw, **unstandardized** 0/1/2 dosage
    matrix of this site's individuals of this ancestry, restricted to the locus
    window. Both are exact integers (:data:`_EXACT_GRAM_MAX_N`); ``n`` is the
    block's cohort size.

    ``variants`` is the ordered variant list they are indexed by — id, position
    and the A1/A2 coding the dosages count. That is Theorem 1 condition (i)
    made explicit: it is the coordinator's only way to know two sites' matrices
    are aligned, and a site that had flipped a reference allele would show it
    here rather than silently negate the off-diagonals it contributes. It is
    reference metadata (the same ``.bim`` lines every site already holds), not
    anything about individuals.

    Emitted once per locus and reused across every architecture x replicate
    instance there, since none of it depends on the phenotype.
    """

    site: str
    pop: str
    locus_id: str
    variants: pd.DataFrame
    G: np.ndarray
    u: np.ndarray
    n: int
    # Diagnostics, not inputs to the fit: window variants dropped because they
    # were not fully observed in this block (see the module docstring).
    n_incomplete_variants: int = 0
    snp_ids: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.snp_ids = self.variants["snp_id"].to_numpy()


@dataclass
class PhenoAggregate:
    """Phenotype-dependent half of one block, for one instance.

    ``c = X'y``, ``q = 1'y``, ``w = y'y`` over the same rows and the same
    ordered ``snp_ids`` as the matching :class:`GenoAggregate`. Raw scale: ``y``
    is the simulated phenotype as stored, not centred or scaled.
    """

    site: str
    pop: str
    instance: str
    snp_ids: np.ndarray
    c: np.ndarray
    q: float
    w: float
    n: int


def locus_window_variants(bim: pd.DataFrame, locus: pd.Series) -> pd.DataFrame:
    """In-window variants of a ``.bim``, in ``.bed`` order.

    The index is preserved, so it is the variant's row index inside the
    matching ``.bed`` — what :func:`~src.utils.read_bed_variants` wants. Sites
    are disjoint ``--keep`` subsets of the same HAPNEST fileset, so they all
    return the same list in the same order: Theorem 1 condition (i) holds by
    construction here rather than by harmonization effort.
    """
    chrom, start, end = int(locus["chrom"]), int(locus["start_bp"]), int(locus["end_bp"])
    mask = (bim["chrom"].astype(str) == str(chrom)) & bim["bp"].between(start, end)
    return bim.loc[mask]


def _columns_with_missing(X: np.ndarray, chunk: int = _ROW_CHUNK) -> np.ndarray:
    """Boolean mask of columns carrying at least one missing dosage.

    Chunked so the intermediate boolean array stays small next to a
    locus-scale (n x M) block.
    """
    bad = np.zeros(X.shape[1], dtype=bool)
    for start in range(0, X.shape[0], chunk):
        bad |= np.isnan(X[start:start + chunk]).any(axis=0)
    return bad


def _gram(X: np.ndarray) -> np.ndarray:
    """``X'X`` as exact integers in float64."""
    n = X.shape[0]
    work = X if n <= _EXACT_GRAM_MAX_N else X.astype(np.float64)
    G = np.asarray(work.T @ work, dtype=np.float64)
    if not np.array_equal(G, np.rint(G)):
        raise RuntimeError(
            "X'X is not integral — the dosage block is not raw 0/1/2 counts, "
            "or it still carries missing values"
        )
    return G


def _xty(X: np.ndarray, Y: np.ndarray, chunk: int = _ROW_CHUNK) -> np.ndarray:
    """``X'Y`` accumulated in float64 out of a float32 dosage block.

    Done in row chunks rather than by upcasting the whole block: single
    precision would cost ~1e-5 relative on a 50,000-row sum, which is coarser
    than the six significant digits the summary statistics are written at.
    """
    out = np.zeros((X.shape[1], Y.shape[1]), dtype=np.float64)
    for start in range(0, X.shape[0], chunk):
        stop = start + chunk
        out += X[start:stop].astype(np.float64).T @ Y[start:stop]
    return out


def _site_phenotypes(
    cfg: SimulationConfig, site: str, instances: Sequence[str], iids: pd.Series
) -> np.ndarray:
    """Phenotypes for ``instances`` at ``site``, aligned to the site's fam order.

    Returns an ``(n_site, n_instances)`` float64 array. Raises if an instance's
    file does not cover every individual in the fam with a non-missing value:
    a partially phenotyped block would make ``n`` differ between ``G`` and
    ``c``, which no single ``(G, c, u, n, q, w)`` tuple can express (Theorem 1
    condition (ii)).
    """
    root = cfg.resolved_path("ground_truth_dir") / "phenotypes" / site
    Y = np.empty((len(iids), len(instances)), dtype=np.float64)
    for col, inst in enumerate(instances):
        pheno = pd.read_csv(root / f"{inst}.pheno", sep="\t",
                            dtype={"FID": str, "IID": str})
        y = pheno.set_index("IID")["y"].reindex(iids.values).to_numpy(dtype=np.float64)
        if np.isnan(y).any():
            raise RuntimeError(
                f"{site} phenotype {inst}.pheno does not cover "
                f"{int(np.isnan(y).sum())} of its {len(iids)} individuals"
            )
        Y[:, col] = y
    return Y


def reference_variants(cfg: SimulationConfig) -> pd.DataFrame:
    """The study's canonical variant list — Theorem 1 condition (i).

    Indexed by variant id, carrying the position and the A1/A2 coding every
    site's dosages must be expressed against. This is *reference metadata*, the
    same annotation table a public panel or a consortium protocol distributes
    before any data is touched; a real deployment would take it from there
    rather than, as here, from the ``.bim`` of the source fileset the cohorts
    were drawn from. It says nothing about any individual.

    It has to be agreed in advance and it has to be used. In this package it is
    not optional bookkeeping: the per-site filesets were cut with
    ``plink --keep --make-bed`` and *without* ``--keep-allele-order``, so PLINK
    1.9 set A1 to each site's own minor allele — about 6% of chr1 variants are
    coded the opposite way round at ANL and Covenant. Summing those Grams
    without harmonizing first is exactly trap 2 of algo.md A.10: the flipped
    site contributes ``2 - x`` where the others contribute ``x``, which adds up
    silently and corrupts every off-diagonal that variant touches.
    """
    bim = read_bim(cfg.resolved_path("hapnest_dir") / f"chr{cfg.chromosome}.bim")
    if bim["snp_id"].duplicated().any():
        raise RuntimeError("Reference variant list has duplicate variant ids")
    return bim.set_index("snp_id")[["chrom", "bp", "a1", "a2"]]


@dataclass
class SiteIndex:
    """A site's chromosome-wide ``.bim``/``.fam``, read once and harmonized once.

    None of it depends on the locus, and the chr1 ``.bim`` is half a million
    rows, so re-reading it per locus would cost more than every Gram in the run.

    ``flip`` marks the variants whose A1/A2 this site has the other way round
    from :func:`reference_variants`; :func:`site_aggregates` recodes those
    dosages to ``2 - x`` before forming any moment, and ``bim`` is rewritten to
    the canonical alleles so what the site emits is already harmonized.
    """

    site: str
    prefix: Path
    bim: pd.DataFrame
    fam: pd.DataFrame
    manifest: pd.DataFrame
    flip: np.ndarray

    @classmethod
    def load(cls, cfg: SimulationConfig, site: str,
             reference: pd.DataFrame | None = None) -> "SiteIndex":
        prefix = cfg.site_dir(site) / f"{site}_chr{cfg.chromosome}"
        bim = read_bim(prefix.with_suffix(".bim"))
        flip = np.zeros(len(bim), dtype=bool)
        if reference is not None:
            ref = reference.reindex(bim["snp_id"].to_numpy())
            if ref["a1"].isna().any():
                raise RuntimeError(
                    f"{site} holds {int(ref['a1'].isna().sum())} variant(s) absent "
                    f"from the reference variant list")
            a1, a2 = bim["a1"].to_numpy(), bim["a2"].to_numpy()
            r1, r2 = ref["a1"].to_numpy(), ref["a2"].to_numpy()
            same = (a1 == r1) & (a2 == r2)
            flip = (a1 == r2) & (a2 == r1)
            if not (same | flip).all():
                bad = bim.loc[~(same | flip), "snp_id"].to_numpy()
                raise RuntimeError(
                    f"{site} codes {len(bad)} variant(s) with alleles the reference "
                    f"does not list at all, e.g. {bad[0]}: not a strand/order swap")
            bim = bim.assign(a1=r1, a2=r2)
        return cls(site=site, prefix=prefix, bim=bim,
                   fam=read_fam(prefix.with_suffix(".fam")),
                   manifest=read_site_manifest(cfg.site_dir(site) / f"{site}_manifest.tsv"),
                   flip=flip)

    @property
    def n_flipped(self) -> int:
        return int(self.flip.sum())


def site_aggregates(
    cfg: SimulationConfig, site: str, locus: pd.Series,
    pops: Sequence[str], instances: Sequence[str],
    index: SiteIndex | None = None,
) -> tuple[list[GenoAggregate], list[PhenoAggregate]]:
    """THE SITE STAGE: emit this site's blocks for one locus. Reads no other site.

    One pass over the site's own locus window. For each superpopulation the
    site holds (intersected with ``pops``, the ancestry columns the coordinator
    asked for) this returns the §7 aggregates: one :class:`GenoAggregate` and
    one :class:`PhenoAggregate` per instance.

    The only inputs are ``{site}_chr{N}.{bed,bim,fam}``, ``{site}_manifest.tsv``
    (which of this site's individuals are which ancestry) and this site's own
    ``.pheno`` files. Nothing is standardized here — that is the coordinator's
    job and doing it here is the §9 bias.

    ``index`` supplies an already-loaded :class:`SiteIndex`; without one this
    reads the site's ``.bim``/``.fam``/manifest itself.
    """
    idx = index if index is not None else SiteIndex.load(cfg, site)
    if idx.site != site:
        raise ValueError(f"SiteIndex is for {idx.site}, not {site}")
    prefix, bim, fam, manifest = idx.prefix, idx.bim, idx.fam, idx.manifest
    variants = locus_window_variants(bim, locus)
    if variants.empty:
        raise RuntimeError(f"No variants in {locus['locus_id']} window at {site}")

    row_of = pd.Series(np.arange(len(fam)), index=fam["IID"].values)
    Y = _site_phenotypes(cfg, site, instances, fam["IID"])

    X = read_bed_variants(prefix, variants.index.to_numpy(), len(fam))
    # Recode the variants this site has the other way round, so that every
    # site's column j counts the same allele and the Grams are addable.
    flip = idx.flip[variants.index.to_numpy()]
    if flip.any():
        X[:, flip] = 2.0 - X[:, flip]
    variants = variants[["chrom", "snp_id", "bp", "a1", "a2"]].reset_index(drop=True)
    geno: list[GenoAggregate] = []
    pheno: list[PhenoAggregate] = []
    try:
        for pop in pops:
            members = manifest.loc[manifest["superpopulation"] == pop, "IID"]
            if members.empty:
                continue  # this site holds none of this ancestry
            rows = np.sort(row_of.reindex(members.values).to_numpy(dtype=np.int64))
            Xp = X[rows]
            keep = ~_columns_with_missing(Xp)
            if not keep.all():
                Xp = Xp[:, keep]
            n = int(Xp.shape[0])
            block_vars = (variants if keep.all()
                          else variants.loc[keep].reset_index(drop=True))
            geno.append(GenoAggregate(
                site=site, pop=pop, locus_id=str(locus["locus_id"]),
                variants=block_vars, G=_gram(Xp), u=Xp.sum(axis=0, dtype=np.float64),
                n=n, n_incomplete_variants=int((~keep).sum()),
            ))
            Yp = Y[rows]
            C = _xty(Xp, Yp)
            q = Yp.sum(axis=0)
            w = np.einsum("ij,ij->j", Yp, Yp)
            for col, inst in enumerate(instances):
                pheno.append(PhenoAggregate(
                    site=site, pop=pop, instance=inst,
                    snp_ids=block_vars["snp_id"].to_numpy(),
                    c=C[:, col], q=float(q[col]), w=float(w[col]), n=n,
                ))
            del Xp, Yp
    finally:
        del X
    return geno, pheno


# ---------------------------------------------------------------------------
# Stage 2a — the coordinator sums (Theorem 1)
# ---------------------------------------------------------------------------

@dataclass
class PooledGeno:
    """One ancestry's summed genotype moments: ``sum_k (G, u, n)`` over sites."""

    pop: str
    variants: pd.DataFrame
    G: np.ndarray
    u: np.ndarray
    n: int
    sites: tuple[str, ...] = ()
    n_incomplete_variants: int = 0

    @property
    def snp_ids(self) -> np.ndarray:
        return self.variants["snp_id"].to_numpy()


@dataclass
class PooledPheno:
    """One ancestry's summed phenotype moments: ``sum_k (c, q, w, n)``."""

    pop: str
    c: np.ndarray
    q: float
    w: float
    n: int


def _common_variants(blocks: Sequence[GenoAggregate]) -> pd.DataFrame:
    """Ordered variant table every block covers, in the first block's order.

    Identical across sites in the normal case — the fast path checks for that
    first. They differ only when a site dropped a partially observed variant,
    in which case the intersection is what all sites can speak about.
    """
    first = blocks[0].variants
    if all(np.array_equal(first["snp_id"].to_numpy(), b.snp_ids) for b in blocks[1:]):
        return first
    keep = np.ones(len(first), dtype=bool)
    for b in blocks[1:]:
        keep &= np.isin(first["snp_id"].to_numpy(), b.snp_ids)
    return first.loc[keep].reset_index(drop=True)


def _check_harmonized(blocks: Sequence[GenoAggregate], variants: pd.DataFrame) -> None:
    """Theorem 1 condition (i): every site codes these variants the same way.

    A site that had A1 and A2 the other way round would be contributing
    ``2 - x`` where the others contribute ``x``, which survives summation
    silently and negates every off-diagonal that variant touches. Cheap to
    check, impossible to debug downstream.
    """
    ref = variants.set_index("snp_id")
    for b in blocks[1:]:
        got = b.variants.set_index("snp_id").reindex(ref.index)
        bad = (got["a1"] != ref["a1"]) | (got["a2"] != ref["a2"]) | (got["bp"] != ref["bp"])
        if bad.any():
            raise RuntimeError(
                f"{b.site} disagrees with {blocks[0].site} on the coding of "
                f"{int(bad.sum())} variant(s) in the {b.pop} block, e.g. "
                f"{ref.index[bad][0]} — the variant lists are not harmonized"
            )


def _align(snp_ids: np.ndarray, target: np.ndarray) -> np.ndarray | None:
    """Positional index taking ``snp_ids``-indexed data onto ``target`` order."""
    if np.array_equal(snp_ids, target):
        return None
    idx = pd.Index(snp_ids).get_indexer(target)
    if (idx < 0).any():
        raise RuntimeError("A block is missing variants the pooled list requires")
    return idx


def pool_geno(blocks: Sequence[GenoAggregate]) -> PooledGeno:
    """Sum one ancestry's genotype aggregates across sites — algo.md §7/§8.

    ``X_p'X_p = sum_k G_{p,k}`` etc.: a matrix product over rows decomposes
    over any partition of those rows, so this is exact and unconditional. The
    result is still raw and unstandardized — see :func:`pooled_ld`.
    """
    pops = {b.pop for b in blocks}
    if len(pops) != 1:
        raise ValueError(f"pool_geno got mixed ancestries: {sorted(pops)}")
    target = _common_variants(blocks)
    _check_harmonized(blocks, target)
    ids = target["snp_id"].to_numpy()
    M = len(ids)
    G = np.zeros((M, M), dtype=np.float64)
    u = np.zeros(M, dtype=np.float64)
    n = 0
    for b in blocks:
        idx = _align(b.snp_ids, ids)
        G += b.G if idx is None else b.G[np.ix_(idx, idx)]
        u += b.u if idx is None else b.u[idx]
        n += b.n
    return PooledGeno(
        pop=blocks[0].pop, variants=target, G=G, u=u, n=n,
        sites=tuple(b.site for b in blocks),
        n_incomplete_variants=sum(b.n_incomplete_variants for b in blocks),
    )


def pool_pheno(blocks: Sequence[PhenoAggregate], snp_ids: np.ndarray) -> PooledPheno:
    """Sum one ancestry's phenotype aggregates across sites, onto ``snp_ids``."""
    c = np.zeros(len(snp_ids), dtype=np.float64)
    q = w = 0.0
    n = 0
    for b in blocks:
        idx = _align(b.snp_ids, snp_ids)
        c += b.c if idx is None else b.c[idx]
        q += b.q
        w += b.w
        n += b.n
    return PooledPheno(pop=blocks[0].pop, c=c, q=q, w=w, n=n)


# ---------------------------------------------------------------------------
# Stage 2b — the coordinator standardizes ONCE (eq 9.1-9.2)
# ---------------------------------------------------------------------------

def pooled_ld(G: np.ndarray, u: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Pooled LD matrix and column SDs from raw moments — algo.md eq (9.1)-(9.2).

    ``sd_j = sqrt(G_jj/n - xbar_j^2)`` and
    ``R_jk = (G_jk/n - xbar_j xbar_k) / (sd_j sd_k)``, i.e. the ordinary
    correlation over the pooled cohort. The population-vs-sample variance
    convention cancels in the ratio, so this is exactly what ``plink --r``
    computes on the same individuals.

    Standardizing *here*, against pooled moments, rather than at each site
    against its own, is the whole point of shipping raw ``G`` and ``u``
    (algo.md §9, A.10). Zero-variance columns get correlation 0, matching how
    SuSiEx reads PLINK's NaNs (``data.cpp:589-592``); the diagonal is forced to
    1 as SuSiEx forces it (``data.cpp:605``).
    """
    xbar = u / n
    cov = G / n
    cov -= np.outer(xbar, xbar)
    sd = np.sqrt(np.clip(np.diag(cov).copy(), 0.0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        R = cov / np.outer(sd, sd)
    R[~np.isfinite(R)] = 0.0
    np.fill_diagonal(R, 1.0)
    return R, sd


def pooled_beta_hat(
    G_diag: np.ndarray, u: np.ndarray, c: np.ndarray, q: float, w: float, n: int
) -> np.ndarray:
    """Standardized marginal effects from raw moments — algo.md eq (9.2).

    ``beta_hat_j = (c_j - n xbar_j ybar) / (n sd_j sd_y)``, which is the
    Pearson correlation of dosage with phenotype: exactly the quantity SuSiEx's
    model calls ``beta``. SuSiEx does not read it directly — it reconstructs
    ``z/sqrt(n)`` from the ``BETA``/``SE`` columns of the summary-statistics
    file (``data.cpp:311``), which differs from this by the usual
    ``t``-vs-correlation factor of order 1/n. :func:`pooled_sumstats` writes
    those columns; this function is the §9 quantity itself, used for
    ``tau^2`` and for checking the two agree.
    """
    sxx = G_diag - u * u / n
    sxy = c - u * q / n
    syy = w - q * q / n
    with np.errstate(divide="ignore", invalid="ignore"):
        beta_hat = sxy / np.sqrt(sxx * syy)
    return np.nan_to_num(beta_hat, nan=0.0, posinf=0.0, neginf=0.0)


def pooled_tau_sq(beta_hat: np.ndarray) -> float:
    """``tau_p^2 = max_j beta_hat_pj^2`` — algo.md §3.5, computed after pooling.

    Note 1 of §9: because it is derived at the coordinator from already-pooled
    moments it is exact, and the sites never have to agree on a prior variance
    in advance. SuSiEx recomputes its own from the sumstats
    (``data.cpp:482``); this is the same quantity to O(1/n).
    """
    return float(np.max(beta_hat ** 2)) if beta_hat.size else 0.0


def pooled_sumstats(
    variants: pd.DataFrame, G_diag: np.ndarray, u: np.ndarray,
    c: np.ndarray, q: float, w: float, n: int,
) -> pd.DataFrame:
    """The per-ancestry GWAS, in closed form from the pooled moments.

    SuSiEx's file interface wants ``BETA``/``SE``/``P`` per variant, from which
    it forms ``beta_hat = (BETA/SE)/sqrt(n)``. The simple regression of ``y``
    on one dosage column with an intercept has a closed form in exactly the
    moments the sites shipped, so the coordinator reproduces ``plink2 --glm``
    without a genotype in sight::

        beta = Sxy/Sxx,  SSE = Syy - beta Sxy,  se = sqrt(SSE/(n-2)/Sxx)

    with ``t = beta/se`` on ``n-2`` degrees of freedom. Columns and order match
    :data:`src.fine_mapping._SUMSTATS_COLS`, so the same SuSiEx column-index
    arguments serve both paths. Zero-variance columns are dropped, mirroring
    the ``dropna`` in :func:`src.fine_mapping.run_gwas`.
    """
    if n < 3:
        raise ValueError(f"Need n >= 3 for a regression, got {n}")
    sxx = G_diag - u * u / n
    sxy = c - u * q / n
    syy = w - q * q / n
    ok = sxx > 0
    beta = np.divide(sxy, sxx, out=np.zeros_like(sxy), where=ok)
    sse = syy - beta * sxy
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(sse / (n - 2) / sxx)
        stat = beta / se
    pval = 2.0 * stats.t.sf(np.abs(stat), n - 2)
    out = pd.DataFrame({
        "chr": variants["chrom"].astype(str).to_numpy(),
        "snp": variants["snp_id"].to_numpy(),
        "bp": variants["bp"].to_numpy(),
        "A1": variants["a1"].to_numpy(),
        "A2": variants["a2"].to_numpy(),
        "beta": beta, "se": se, "stat": stat, "p": pval,
    })[list(_SUMSTATS_COLS)]
    return out.loc[ok & np.isfinite(se) & (se > 0)].reset_index(drop=True)


def write_sumstats(sumstats: pd.DataFrame, path: Path) -> Path:
    """Write a SuSiEx summary-statistics file.

    Six significant digits: that is what ``plink2 --glm`` writes, so the
    federated file carries the same digits the centralized one does and the two
    paths hand SuSiEx byte-identical numbers rather than two roundings of them.
    """
    sumstats.to_csv(path, sep="\t", index=False, float_format="%.6g")
    return path


# ---------------------------------------------------------------------------
# Stage 2c — the coordinator writes SuSiEx's LD panel by hand
# ---------------------------------------------------------------------------

def write_ld_panel(
    prefix: Path, variants: pd.DataFrame, R: np.ndarray, freq: np.ndarray, n: int
) -> Path:
    """Write the three files SuSiEx probes for at an ``--ld_file`` prefix.

    The centralized path gets these from ``plink --r square bin4`` and
    ``plink --freq`` on a genotype panel it holds. The coordinator holds no
    genotypes, so it writes them directly. Matching the reader exactly
    (``main.cpp:204-222``, ``data.cpp:126-188`` and ``505-610``):

    ``{prefix}.ld.bin``
        ``M x M`` float32, row-major, no header. SuSiEx checks the file is
        exactly ``M^2 * 4`` bytes against the ``_ref.bim`` line count and
        aborts otherwise, so the two must be written from the same variant
        list in the same order.
    ``{prefix}_ref.bim``
        Standard 6-column PLINK ``.bim``. Every line's chromosome must equal
        the ``--chr`` argument or SuSiEx exits.
    ``{prefix}_frq.frq``
        A header line then one line per variant, in ``_ref.bim`` order, whose
        third and fourth tokens must repeat that line's A1 and A2 exactly —
        SuSiEx cross-checks them and aborts on any disagreement. The frequency
        is the A1 frequency, ``u / 2n``, which is what ``plink --freq
        --keep-allele-order`` reports. It feeds only the reported allele
        frequency and per-allele effect columns, not the fit.
    """
    M = len(variants)
    if R.shape != (M, M):
        raise ValueError(f"LD matrix is {R.shape}, expected ({M}, {M})")
    np.asarray(R, dtype=np.float32).tofile(f"{prefix}.ld.bin")

    bim = variants[["chrom", "snp_id", "bp", "a1", "a2"]].to_numpy(dtype=object)
    with open(f"{prefix}_ref.bim", "w") as fh:
        for chrom, snp, bp, a1, a2 in bim:
            fh.write(f"{chrom}\t{snp}\t0\t{int(bp)}\t{a1}\t{a2}\n")
    with open(f"{prefix}_frq.frq", "w") as fh:
        fh.write(f"{'CHR':>4} {'SNP':>15} {'A1':>4} {'A2':>4} {'MAF':>12} {'NCHROBS':>8}\n")
        for (chrom, snp, _bp, a1, a2), f in zip(bim, freq):
            fh.write(f"{chrom:>4} {snp:>15} {a1:>4} {a2:>4} {f:>12.4g} {2 * n:>8}\n")
    return prefix


# ---------------------------------------------------------------------------
# The coordinator's ancestry columns
# ---------------------------------------------------------------------------

@dataclass
class FedColumn:
    """One SuSiEx population column, reconstructed from site aggregates.

    Everything here is derived from summed moments: the LD panel is already on
    disk, and ``u``/``g_diag`` are the O(M) moments a per-instance summary
    statistic still needs. ``n`` is the pooled cohort — the same number the
    centralized path passes as ``--n_gwas``.
    """

    pop: str
    n: int
    variants: pd.DataFrame
    u: np.ndarray
    g_diag: np.ndarray
    ld_prefix: Path
    sites: tuple[str, ...] = ()
    dropped_maf: int = 0
    n_incomplete_variants: int = 0

    @property
    def name(self) -> str:
        return self.pop


def build_columns(
    geno_blocks: Sequence[GenoAggregate], pops: Sequence[str],
    out_dir: Path, maf: float,
) -> list[FedColumn]:
    """Pool, standardize once, and write each ancestry's LD panel.

    This is algo.md §7 then §9 in order: sum the sites' raw moments, *then*
    form ``R_p``. The MAF filter is likewise a pooled decision — ``plink --maf``
    keeps variants whose minor allele frequency is at or above the threshold,
    computed here from the pooled ``u/2n`` — so the panel covers the same
    variants the centralized path's ``plink --extract ... --maf`` panel does.

    Genotype-only, hence built once per locus and shared by every instance
    there, exactly as :func:`src.fine_mapping.precompute_ld` is.
    """
    ensure_dir(out_dir)
    by_pop: dict[str, list[GenoAggregate]] = {}
    for block in geno_blocks:
        by_pop.setdefault(block.pop, []).append(block)

    columns: list[FedColumn] = []
    for pop in pops:
        blocks = by_pop.get(pop)
        if not blocks:
            continue
        pooled = pool_geno(blocks)
        freq = pooled.u / (2 * pooled.n)
        keep = np.minimum(pooled.u, 2 * pooled.n - pooled.u) >= 2 * pooled.n * maf - 1e-10
        if not keep.any():
            raise RuntimeError(f"No variant survives MAF {maf} in the pooled {pop} column")
        idx = np.flatnonzero(keep)
        G = pooled.G[np.ix_(idx, idx)]
        u = pooled.u[idx]
        R, _sd = pooled_ld(G, u, pooled.n)
        variants = pooled.variants.loc[keep].reset_index(drop=True)
        prefix = write_ld_panel(out_dir / f"{pop}_ld", variants, R, freq[idx], pooled.n)
        columns.append(FedColumn(
            pop=pop, n=pooled.n, variants=variants, u=u, g_diag=np.diag(G).copy(),
            ld_prefix=prefix, sites=pooled.sites, dropped_maf=int((~keep).sum()),
            n_incomplete_variants=pooled.n_incomplete_variants,
        ))
        del G, R
    if not columns:
        raise RuntimeError(
            "No ancestry columns were reconstructed from the site aggregates")
    return columns


# ---------------------------------------------------------------------------
# One instance
# ---------------------------------------------------------------------------

def fed_finemap_instance(
    locus: pd.Series, arch_id: str, rep: int, truth_snps: list[str],
    columns: Sequence[FedColumn], pheno: dict[str, PooledPheno],
    work_root: Path, susiex: str, plink: str,
    level: float, pval_thresh: float, keep_work: bool, options: dict | None = None,
) -> dict:
    """Coordinator side of one (locus, architecture, replicate): SuSiEx + eval.

    The summary statistics come out of the pooled moments in closed form; the
    LD panels are already on disk from :func:`build_columns`. The SuSiEx call
    is :func:`src.fine_mapping.run_susiex` unchanged — same binary, same
    arguments, same precomputed-LD path — so the only thing that differs from
    the centralized run is *where the numbers came from*. Result rows carry the
    same schema, so the two tables are directly comparable.
    """
    t0 = time.time()
    chrom, start, end = int(locus["chrom"]), int(locus["start_bp"]), int(locus["end_bp"])
    inst = f"{locus['locus_id']}_{arch_id}_rep{rep}"
    priv = ensure_dir(work_root / inst)
    row = {
        "locus_id": locus["locus_id"], "architecture_id": arch_id, "replicate": rep,
        "stratum": locus.get("stratum", ""),
    }
    min_p: dict[str, float] = {}
    try:
        sumstats, ld_list, n_list = [], [], []
        for col in columns:
            moments = pheno[col.pop]
            if moments.n != col.n:
                raise RuntimeError(
                    f"{col.pop}: phenotype aggregates cover {moments.n} individuals "
                    f"but genotype aggregates cover {col.n}"
                )
            ss = pooled_sumstats(col.variants, col.g_diag, col.u,
                                 moments.c, moments.q, moments.w, col.n)
            min_p[col.pop] = float(ss["p"].min()) if len(ss) else np.nan
            sumstats.append(write_sumstats(ss, priv / f"{col.pop}.sumstats"))
            ld_list.append(col.ld_prefix)
            n_list.append(col.n)
        run_susiex(sumstats, ld_list, n_list, chrom, start, end,
                   priv, "cs", susiex, plink, level, pval_thresh, options=options)
        row.update(parse_susiex(priv, "cs", truth_snps))
        row["error"] = ""
    except Exception as exc:  # keep the shard alive; record the failure
        get_logger().warning("fed fine-map FAILED %s: %s", inst, exc)
        row.update(dict.fromkeys(_METRIC_COLS, np.nan))
        row.update({"any_causal_captured": False, "converged": False,
                    "fit_status": "execution_failure", "error": str(exc)[:200]})
    finally:
        # Sorted by ancestry, to match the centralized baseline's schema exactly. See
        # the same block in fine_mapping.py: column order matters to SuSiEx and differs
        # between the two paths, but it must not reach the results table, which exists
        # to be compared against the centralized one column for column.
        for name in sorted(col.pop for col in columns):
            row[f"min_p_{name}"] = min_p.get(name, np.nan)
        row["runtime_s"] = round(time.time() - t0, 2)
        from ..inference import archive_instance
        archive_instance(priv, work_root.parent / "artifacts" / inst, truth_snps, row)
        if not keep_work:
            shutil.rmtree(priv, ignore_errors=True)
    return row


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def _fed_fm_dir(cfg: SimulationConfig) -> Path:
    return cfg.resolved_path("reports_dir") / "fed_fine_mapping"


def run_fed_fine_mapping(
    cfg: SimulationConfig, shard_index: int = 0, n_shards: int = 1,
    n_workers: int = 8, level: float = 0.95, pval_thresh: float = 1e-5,
    keep_work: bool = False, limit: int | None = None, maf: float = 0.005,
) -> None:
    logger = get_logger()
    susiex = _resolve_binary("SuSiEx", cfg.repo_root)
    # SuSiEx's argument validator demands --plink unconditionally (main.cpp:274)
    # but never runs it once the precomputed LD files are present, which on this
    # path they always are. A coordinator with no genotypes has no reason to
    # hold PLINK, so a missing binary is not an error here.
    try:
        plink = _resolve_binary(cfg.tools.plink, cfg.repo_root)
    except FileNotFoundError:
        plink = cfg.tools.plink

    fm_dir = ensure_dir(_fed_fm_dir(cfg))
    work_root = ensure_dir(fm_dir / "work")

    # Ancestry columns are locus-independent. Planned exactly as the centralized
    # path plans them — but no pooled --keep file is ever written, because the
    # coordinator never assembles a cross-site cohort.
    planned = plan_columns(
        {s: dict(sc.composition) for s, sc in cfg.sites.items()},
        cfg.superpopulations, cfg.fine_mapping.min_gwas_n,
    )
    pops = [c.name for c in planned]
    expected_n = {c.name: c.n for c in planned}
    logger.info("Federated SuSiEx ancestry columns (min_gwas_n=%d): %s",
                cfg.fine_mapping.min_gwas_n,
                ", ".join(f"{c.name}(n={c.n} from {'+'.join(c.sites)})" for c in planned))

    loci = pd.read_csv(
        cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t"
    ).reset_index(drop=True)
    manifest = pd.read_csv(
        cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv", sep="\t")
    truth_map = {
        (r.locus_id, r.architecture_id, int(r.replicate)): str(r.causal_snp_ids).split(",")
        for r in manifest.itertuples(index=False)
    }
    architectures = build_architecture_grid(
        ncsl=cfg.architecture.ncsl, h2=cfg.architecture.h2,
        rg=cfg.architecture.rg, mode=cfg.architecture.factorial_mode,
    )
    n_rep = cfg.architecture.replicates

    # The canonical variant coding is agreed once, then each site opens its own
    # chromosome index once (not once per locus) and harmonizes against it.
    reference = reference_variants(cfg)
    indexes = {site: SiteIndex.load(cfg, site, reference) for site in cfg.sites}
    for site, index in indexes.items():
        if index.n_flipped:
            logger.info("%s: %d of %d variants recoded to the reference allele order",
                        site, index.n_flipped, len(index.bim))

    if n_shards > 1:
        loci = loci.loc[(loci.index % n_shards) == shard_index].copy()
    logger.info("Federated fine-mapping %d loci x %d arch x %d rep "
                "(shard %d/%d, %d workers)",
                len(loci), len(architectures), n_rep, shard_index, n_shards, n_workers)

    all_rows: list[dict] = []
    for _, locus in loci.iterrows():
        jobs = []
        for arch in architectures:
            for rep in range(n_rep):
                truth = truth_map.get((locus["locus_id"], arch.architecture_id, rep))
                if truth is None:
                    continue
                jobs.append((arch.architecture_id, rep, truth))
                if limit is not None and len(all_rows) + len(jobs) >= limit:
                    break
            if limit is not None and len(all_rows) + len(jobs) >= limit:
                break
        if not jobs:
            continue
        instances = [f"{locus['locus_id']}_{aid}_rep{rep}" for (aid, rep, _t) in jobs]

        # --- stage 1: each site, on its own data, emits its blocks -----------
        geno_blocks: list[GenoAggregate] = []
        pheno_blocks: list[PhenoAggregate] = []
        for site in cfg.sites:
            geno, pheno = site_aggregates(cfg, site, locus, pops, instances,
                                          index=indexes[site])
            geno_blocks.extend(geno)
            pheno_blocks.extend(pheno)
        logger.info("  locus %s: %d (site, ancestry) blocks uplinked",
                    locus["locus_id"], len(geno_blocks))

        # --- stage 2: the coordinator sums, standardizes once, fits ----------
        locus_dir = ensure_dir(work_root / f"{locus['locus_id']}_ld")
        columns = build_columns(geno_blocks, pops, locus_dir, maf)
        for col in columns:
            if col.n != expected_n.get(col.pop):
                raise ValueError(f"{col.pop}: pooled n={col.n} differs from planned n={expected_n.get(col.pop)}")
            if col.n_incomplete_variants:
                logger.warning(
                    "  %s: %d window variant(s) dropped as not fully observed; the "
                    "centralized comparator uses the same complete-variant policy", col.pop, col.n_incomplete_variants)
        del geno_blocks

        by_inst: dict[str, dict[str, list[PhenoAggregate]]] = {}
        for block in pheno_blocks:
            by_inst.setdefault(block.instance, {}).setdefault(block.pop, []).append(block)
        del pheno_blocks
        pooled_pheno = {
            inst: {col.pop: pool_pheno(pops_blocks[col.pop],
                                       col.variants["snp_id"].to_numpy())
                   for col in columns}
            for inst, pops_blocks in by_inst.items()
        }

        rows = Parallel(n_jobs=n_workers, backend="loky")(
            delayed(fed_finemap_instance)(
                locus, aid, rep, truth, columns, pooled_pheno[inst], work_root,
                susiex, plink, level, pval_thresh, keep_work, cfg.fine_mapping.inference_options())
            for (aid, rep, truth), inst in zip(jobs, instances)
        )
        all_rows.extend(rows)
        if not keep_work:
            shutil.rmtree(locus_dir, ignore_errors=True)
        logger.info("  locus %s done (%d instances, cumulative %d)",
                    locus["locus_id"], len(jobs), len(all_rows))
        if limit is not None and len(all_rows) >= limit:
            break

    df = pd.DataFrame(all_rows)
    if n_shards > 1:
        tag = f"part{shard_index:03d}of{n_shards:03d}"
        df.to_csv(fm_dir / f"fed_fm_results.{tag}.tsv", sep="\t", index=False)
        logger.info("Shard %d/%d wrote %d result rows (%s)",
                    shard_index, n_shards, len(df), tag)
    else:
        df.to_csv(fm_dir / "fed_fm_results.tsv", sep="\t", index=False)
        logger.info("Wrote %d result rows -> %s", len(df), fm_dir / "fed_fm_results.tsv")
        write_rollup(cfg, df)

    from ..inference import require_successful_fits
    require_successful_fits(df)


def merge_fed_fm_results(cfg: SimulationConfig, n_shards: int) -> None:
    """REDUCE: concat per-shard parts, sort, write fed_fm_results.tsv + rollup."""
    logger = get_logger()
    fm_dir = _fed_fm_dir(cfg)
    parts = []
    for r in range(n_shards):
        p = fm_dir / f"fed_fm_results.part{r:03d}of{n_shards:03d}.tsv"
        if not p.exists():
            raise FileNotFoundError(f"Missing federated fine-mapping result part: {p}")
        parts.append(pd.read_csv(p, sep="\t"))
    df = pd.concat(parts, ignore_index=True).sort_values(
        ["locus_id", "architecture_id", "replicate"]).reset_index(drop=True)
    df.to_csv(fm_dir / "fed_fm_results.tsv", sep="\t", index=False)
    for r in range(n_shards):
        (fm_dir / f"fed_fm_results.part{r:03d}of{n_shards:03d}.tsv").unlink()
    logger.info("Merged %d shards -> %d rows", n_shards, len(df))
    write_rollup(cfg, df)


def write_rollup(cfg: SimulationConfig, df: pd.DataFrame) -> None:
    """Aggregate power/credible-set metrics by architecture and by stratum.

    Same aggregation and same column names as the centralized rollups, written
    beside the federated results so the two tables line up row for row.
    """
    if df.empty:
        return
    fm_dir = _fed_fm_dir(cfg)
    arch = build_architecture_grid(
        ncsl=cfg.architecture.ncsl, h2=cfg.architecture.h2,
        rg=cfg.architecture.rg, mode=cfg.architecture.factorial_mode)
    meta = pd.DataFrame([{"architecture_id": a.architecture_id, "ncsl": a.ncsl,
                          "h2_target": a.h2, "rg": a.rg} for a in arch])
    d = df.merge(meta, on="architecture_id", how="left")

    def _agg(g: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "n_instances": len(g),
            "power_any_causal": g["any_causal_captured"].mean(),
            "mean_n_cs": g["n_credible_sets"].mean(),
            "median_best_cs_size": g["best_cs_size"].median(),
            "mean_causal_pip": g["causal_pip_mean"].mean(),
        })

    for by, name in [(["ncsl", "h2_target", "rg"], "by_architecture"),
                     (["stratum", "rg"], "by_stratum_rg")]:
        cols = [c for c in by if c in d.columns]
        if cols:
            d.groupby(cols).apply(_agg, include_groups=False).reset_index().to_csv(
                fm_dir / f"fed_fm_rollup_{name}.tsv", sep="\t", index=False)
    get_logger().info("Wrote rollups -> %s/fed_fm_rollup_*.tsv", fm_dir)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="FedFM federated SuSiEx fine-mapping (site aggregates -> coordinator)")
    parser.add_argument("--config", default="config/simulation_config.yaml")
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--n-shards", type=int, default=1)
    parser.add_argument("--merge", action="store_true",
                        help="REDUCE: merge per-shard result parts, then exit")
    parser.add_argument("--n-workers", type=int, default=8)
    parser.add_argument("--level", type=float, default=0.95,
                        help="credible-set coverage level (SuSiEx --level)")
    parser.add_argument("--pval-thresh", type=float, default=1e-5,
                        help="SuSiEx marginal p-value filter for credible sets")
    parser.add_argument("--maf", type=float, default=None,
                        help="MAF filter applied at the coordinator to pooled "
                             "frequencies; defaults to the scenario's fine_mapping.maf, "
                             "so that this path and the centralized baseline filter "
                             "identically")
    parser.add_argument("--keep-work", action="store_true",
                        help="keep per-instance sumstats/LD scratch (debug)")
    parser.add_argument("--limit", type=int, default=None,
                        help="stop after N instances (smoke test)")
    args = parser.parse_args(list(argv) if argv is not None else None)
    cfg = load_config(args.config)
    setup_logging(log_dir=cfg.resolved_path("logs_dir"))

    if args.merge:
        merge_fed_fm_results(cfg, args.n_shards)
    else:
        if args.n_shards > 1 and not (0 <= args.shard_index < args.n_shards):
            parser.error("need 0 <= --shard-index < --n-shards")
        run_fed_fine_mapping(cfg, args.shard_index, args.n_shards, args.n_workers,
                             args.level, args.pval_thresh, args.keep_work, args.limit,
                             cfg.fine_mapping.maf if args.maf is None else args.maf)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
