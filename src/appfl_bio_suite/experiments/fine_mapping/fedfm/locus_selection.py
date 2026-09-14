"""LD-divergence scoring and stratified locus selection.

Pipeline:
    1. Generate candidate 1 Mb windows on the target chromosome (sliding 500 kb).
    2. For each window:
         a. From the (shared) variant set, keep variants with MAF > threshold
            in at least one of the three sites.
         b. Thin to ~``ld_tag_snp_target`` tag SNPs via a light in-memory LD prune.
         c. Compute per-site standardized r² matrices among the tag SNPs.
         d. Compute the mean pairwise Frobenius norm of r² differences — the
            "divergence score".
    3. Rank windows by divergence, split into low / medium / high strata, sample
       33 / 33 / 34 loci. Drop overlapping selections, preferring the more extreme
       divergence within the stratum.

Outputs:
    * ``<loci_dir>/candidate_windows.tsv`` — every candidate window scored.
    * ``<loci_dir>/selected_loci.tsv``      — the 100 selected loci with strata.
    * ``<loci_dir>/ld_matrices/<locus_id>_<site>.npz`` — cached per-site r² matrices.
    * ``<loci_dir>/tag_snps/<locus_id>.tsv`` — tag SNP variant IDs used for that locus.

The math (windowing, LD prune, r², Frobenius, stratification) is split out into
pure functions that are unit-testable without any genotype I/O.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .utils import (
    SimulationConfig,
    derive_seed,
    ensure_dir,
    get_logger,
    load_config,
    read_bed_variants,
    read_bim,
    read_fam,
    setup_logging,
)

# Numerical floor for column standard deviation (a monomorphic SNP has sd=0;
# dividing by zero is the silent way to corrupt an r² matrix). Variants below
# this floor are removed from the tag set before standardization.
_MIN_COL_SD = 1e-8


# ---------------------------------------------------------------------------
# Pure-math helpers (unit-testable, no I/O)
# ---------------------------------------------------------------------------


def generate_candidate_windows(
    chrom: str | int,
    min_bp: int,
    max_bp: int,
    window_size_bp: int,
    step_size_bp: int,
) -> pd.DataFrame:
    """Generate sliding-window candidates spanning [min_bp, max_bp].

    The last window is included even if it doesn't fit a full step; this preserves
    coverage at the chromosome end.
    """
    if window_size_bp <= 0 or step_size_bp <= 0:
        raise ValueError("window/step must be positive")
    starts = list(range(min_bp, max(min_bp + 1, max_bp - window_size_bp + 1), step_size_bp))
    if not starts:
        starts = [min_bp]
    rows = []
    for i, start in enumerate(starts):
        end = start + window_size_bp - 1
        rows.append({
            "window_id": f"chr{chrom}_w{i:05d}",
            "chrom": str(chrom),
            "start_bp": int(start),
            "end_bp": int(end),
        })
    return pd.DataFrame(rows)


def compute_maf(dosage: np.ndarray) -> np.ndarray:
    """Per-column minor allele frequency from a dosage matrix.

    NaN entries (missing genotypes) are ignored.
    """
    af = np.nanmean(dosage, axis=0) / 2.0
    return np.minimum(af, 1.0 - af)


def light_ld_prune(
    dosage: np.ndarray,
    window_variants: int,
    r2_threshold: float,
) -> np.ndarray:
    """Greedy LD prune: walk variants left-to-right, drop any with r² above
    threshold against any kept variant within the last ``window_variants``.

    Parameters
    ----------
    dosage : (n_samples, n_variants) float array.
    window_variants : sliding-window size in variant counts.
    r2_threshold : variants with r² ≥ threshold against a recently kept variant
        are dropped.

    Returns
    -------
    np.ndarray of variant indices to keep (sorted ascending).
    """
    n_samples, n_variants = dosage.shape
    if n_variants == 0:
        return np.array([], dtype=np.int64)

    # Standardize once.
    mean = np.nanmean(dosage, axis=0)
    std = np.nanstd(dosage, axis=0)
    valid = std > _MIN_COL_SD
    X = (dosage - mean) / np.where(valid, std, 1.0)
    X = np.where(np.isnan(X), 0.0, X)
    X[:, ~valid] = 0.0  # monomorphic columns contribute 0 correlation

    kept: list[int] = []
    kept_cols: list[np.ndarray] = []
    for j in range(n_variants):
        if not valid[j]:
            continue
        col = X[:, j]
        # Compare only against the last `window_variants` kept variants.
        if kept_cols:
            recent = kept_cols[-window_variants:]
            R = np.column_stack(recent)
            r = (R.T @ col) / n_samples
            if np.any(r * r >= r2_threshold):
                continue
        kept.append(j)
        kept_cols.append(col)
    return np.array(kept, dtype=np.int64)


def standardize_columns(dosage: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (X_std, valid_mask). Monomorphic columns are zeroed out."""
    mean = np.nanmean(dosage, axis=0)
    std = np.nanstd(dosage, axis=0)
    valid = std > _MIN_COL_SD
    X = (dosage - mean) / np.where(valid, std, 1.0)
    X = np.where(np.isnan(X), 0.0, X)
    X[:, ~valid] = 0.0
    return X.astype(np.float32, copy=False), valid


def r2_matrix(dosage: np.ndarray) -> np.ndarray:
    """Compute the variant × variant r² (squared Pearson correlation) matrix.

    Uses the population-variance convention so that ``R[i, i] == 1`` exactly:
    standardize columns by their population std (``nanstd`` ddof=0), then
    correlation = ``(X_std.T @ X_std) / n``.
    """
    n_samples = dosage.shape[0]
    if n_samples < 2:
        raise ValueError(f"Need ≥2 samples for r², got {n_samples}")
    X, _ = standardize_columns(dosage)
    R = (X.T @ X) / n_samples
    return (R * R).astype(np.float32, copy=False)


def frobenius_diff(A: np.ndarray, B: np.ndarray) -> float:
    """Frobenius norm of ``A - B``. Assumes same shape; raises otherwise."""
    if A.shape != B.shape:
        raise ValueError(f"Shape mismatch: {A.shape} vs {B.shape}")
    diff = A - B
    return float(np.sqrt(np.sum(diff * diff)))


def pairwise_divergence(per_site_r2: dict[str, np.ndarray]) -> tuple[float, dict[str, float]]:
    """Mean pairwise Frobenius divergence + per-pair scores.

    Site order is dict iteration order; passing the same dict gives reproducible output.
    """
    site_names = list(per_site_r2.keys())
    pair_scores: dict[str, float] = {}
    for i in range(len(site_names)):
        for j in range(i + 1, len(site_names)):
            a, b = site_names[i], site_names[j]
            pair_scores[f"{a}__{b}"] = frobenius_diff(per_site_r2[a], per_site_r2[b])
    if not pair_scores:
        return 0.0, pair_scores
    return float(np.mean(list(pair_scores.values()))), pair_scores


def stratify_and_sample(
    scored: pd.DataFrame,
    strata: dict[str, int],
    rng: np.random.Generator,
    score_col: str = "divergence_score",
) -> pd.DataFrame:
    """Split ``scored`` into low/medium/high thirds by ``score_col`` and sample.

    Parameters
    ----------
    scored : DataFrame with at least ``score_col``.
    strata : dict with keys ``low``, ``medium``, ``high`` and integer counts.
    rng : numpy Generator used for the within-stratum draw.
    score_col : column name to rank on.

    Returns
    -------
    A copy of the sampled rows with an added ``stratum`` column.
    """
    n = len(scored)
    if n == 0:
        raise ValueError("No windows to stratify")

    sorted_df = scored.sort_values(score_col).reset_index(drop=True)
    # Thirds by row index (ties resolved by sort stability).
    low_cut = n // 3
    high_cut = n - (n // 3)
    sorted_df["stratum"] = "medium"
    sorted_df.loc[: low_cut - 1, "stratum"] = "low"
    sorted_df.loc[high_cut:, "stratum"] = "high"

    selected_parts: list[pd.DataFrame] = []
    for stratum, count in strata.items():
        pool = sorted_df[sorted_df["stratum"] == stratum]
        if len(pool) < count:
            raise ValueError(
                f"Stratum {stratum} has {len(pool)} windows but {count} requested"
            )
        idx = rng.choice(len(pool), size=count, replace=False)
        selected_parts.append(pool.iloc[np.sort(idx)].copy())
    return pd.concat(selected_parts, ignore_index=True)


def _greedy_nonoverlap(df: pd.DataFrame, priority_col: str) -> list[pd.Series]:
    """Greedily admit non-overlapping windows in ascending ``priority_col`` order
    (smaller priority = admitted first). Returns the list of kept rows."""
    kept: list[pd.Series] = []
    for _, row in df.sort_values(priority_col, kind="stable").iterrows():
        overlap = any(
            row["chrom"] == prev["chrom"]
            and not (row["end_bp"] < prev["start_bp"] or row["start_bp"] > prev["end_bp"])
            for prev in kept
        )
        if not overlap:
            kept.append(row)
    return kept


def drop_overlapping_selections(selected: pd.DataFrame, score_col: str = "divergence_score") -> pd.DataFrame:
    """Drop overlapping windows in two passes.

    Pass 1 (within-stratum): for each stratum, keep the window with the more
    extreme within-stratum score (high → highest; low → lowest; medium →
    closest to the stratum median), dropping any window that overlaps an
    already-kept one. This roughly preserves the stratum counts.

    Pass 2 (cross-stratum): windows selected for *different* strata can still
    share a physical interval (the candidate windows overlap 50% by design). A
    final global greedy sweep resolves those, keeping the window that is more
    extreme *within its own stratum* (compared via a per-stratum rank normalised
    to [0, 1), so the comparison is meaningful across strata of different scale).

    Implemented with ``groupby`` iteration rather than ``groupby(...).apply`` so
    it is robust to the pandas≥3.0 change where the grouping column is no longer
    exposed inside the applied group.
    """
    if "stratum" not in selected.columns:
        raise ValueError("selected must have a 'stratum' column")
    if selected.empty:
        return selected.copy()

    work = selected.copy()

    # Per-stratum extremeness priority, normalised to [0, 1): 0 == most extreme.
    pri = pd.Series(np.nan, index=work.index, dtype=float)
    for stratum, sub in work.groupby("stratum"):
        if stratum == "high":
            rank = (-sub[score_col]).rank(method="first")
        elif stratum == "low":
            rank = sub[score_col].rank(method="first")
        else:  # medium: closeness to the stratum median
            median = sub[score_col].median()
            rank = (sub[score_col] - median).abs().rank(method="first")
        pri.loc[sub.index] = (rank - 1.0) / len(sub)
    work["_pri"] = pri

    # Pass 1: within each stratum.
    stage1_rows: list[pd.Series] = []
    for _, sub in work.groupby("stratum"):
        stage1_rows.extend(_greedy_nonoverlap(sub, "_pri"))
    stage1 = pd.DataFrame(stage1_rows)

    # Pass 2: global cross-stratum sweep.
    final_rows = _greedy_nonoverlap(stage1, "_pri")
    final = pd.DataFrame(final_rows).drop(columns="_pri")
    return final.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Window scoring (I/O + math)
# ---------------------------------------------------------------------------


@dataclass
class SiteHandle:
    """All the per-site state the window scorer needs."""

    name: str
    bed_prefix: Path
    bim: pd.DataFrame      # cached .bim DataFrame
    n_samples: int         # number of rows in .fam


def _load_site_handles(cfg: SimulationConfig) -> dict[str, SiteHandle]:
    handles: dict[str, SiteHandle] = {}
    for site_name in cfg.sites:
        site_dir = cfg.site_dir(site_name)
        prefix = site_dir / f"{site_name}_chr{cfg.chromosome}"
        bim = read_bim(str(prefix) + ".bim")
        fam = read_fam(str(prefix) + ".fam")
        handles[site_name] = SiteHandle(
            name=site_name,
            bed_prefix=prefix,
            bim=bim,
            n_samples=len(fam),
        )
    return handles


def _select_window_variants(
    handles: dict[str, SiteHandle],
    window: pd.Series,
    maf_threshold: float,
) -> pd.DataFrame:
    """Return the set of variants in [start_bp, end_bp] passing the
    "MAF > threshold in at least one site" filter, intersected across sites.

    The returned DataFrame has columns ``snp_id``, ``bp``, and one ``vidx_<site>``
    integer column per site giving the variant's row index in that site's .bim.
    """
    chrom = str(window["chrom"])
    start = int(window["start_bp"])
    end = int(window["end_bp"])

    # Intersect SNP IDs across sites (they should be identical since all sites
    # derive from the same HAPNEST source, but guard anyway).
    per_site_in_window: dict[str, pd.DataFrame] = {}
    for site_name, h in handles.items():
        bim = h.bim
        mask = (bim["chrom"].astype(str) == chrom) & (bim["bp"] >= start) & (bim["bp"] <= end)
        sub = bim.loc[mask, ["snp_id", "bp"]].copy()
        sub["vidx"] = bim.index[mask].to_numpy()
        per_site_in_window[site_name] = sub

    site_names = list(handles.keys())
    if not site_names:
        return pd.DataFrame(columns=["snp_id", "bp"])

    merged = per_site_in_window[site_names[0]].rename(columns={"vidx": f"vidx_{site_names[0]}"})
    for site_name in site_names[1:]:
        right = per_site_in_window[site_name][["snp_id", "vidx"]].rename(
            columns={"vidx": f"vidx_{site_name}"}
        )
        merged = merged.merge(right, on="snp_id", how="inner")
    merged.sort_values("bp", inplace=True)
    merged.reset_index(drop=True, inplace=True)

    if merged.empty:
        return merged

    # MAF filter: any-site MAF > threshold. Reads dosages once per site.
    pass_any = np.zeros(len(merged), dtype=bool)
    for site_name, h in handles.items():
        vidx = merged[f"vidx_{site_name}"].to_numpy()
        dosage = read_bed_variants(h.bed_prefix, vidx, h.n_samples)
        maf = compute_maf(dosage)
        pass_any |= maf > maf_threshold
    return merged.loc[pass_any].reset_index(drop=True)


def _pick_tag_snps(
    handles: dict[str, SiteHandle],
    variants: pd.DataFrame,
    target_n: int,
    prune_window: int,
    prune_threshold: float,
    seed_components: tuple,
    master_seed: int,
) -> pd.DataFrame:
    """Down-select variants to ~target_n tag SNPs.

    Strategy:
        1. If already at/below target, return as-is.
        2. Spatial uniform thin to ~2× target.
        3. Light LD prune on the pooled-site dosage (concatenated across sites).
        4. If still over target, uniform-stride down to target.
    """
    if len(variants) <= target_n:
        return variants

    rng = np.random.default_rng(derive_seed(master_seed, "tag_snp", *seed_components))

    # Step 1: spatial uniform thin to ~2× target.
    n = len(variants)
    desired = min(n, max(target_n * 2, target_n + 50))
    pick_idx = np.linspace(0, n - 1, desired).round().astype(np.int64)
    pick_idx = np.unique(pick_idx)
    sub = variants.iloc[pick_idx].reset_index(drop=True)

    # Step 2: light LD prune. Use pooled dosage from all sites (concat rows).
    pooled_blocks: list[np.ndarray] = []
    for site_name, h in handles.items():
        vidx = sub[f"vidx_{site_name}"].to_numpy()
        dosage = read_bed_variants(h.bed_prefix, vidx, h.n_samples)
        alleles = h.bim.iloc[vidx][["a1", "a2"]].to_numpy()
        if not pooled_blocks:
            reference_alleles = alleles
        else:
            same = (alleles == reference_alleles).all(axis=1)
            flip = (alleles[:, ::-1] == reference_alleles).all(axis=1)
            if not (same | flip).all():
                raise ValueError("Incompatible variant alleles during pooled tag pruning")
            dosage[:, flip & ~same] = 2.0 - dosage[:, flip & ~same]
        pooled_blocks.append(dosage)
    pooled = np.concatenate(pooled_blocks, axis=0)
    keep_local = light_ld_prune(pooled, prune_window, prune_threshold)
    sub = sub.iloc[keep_local].reset_index(drop=True)

    # Step 3: if still over target, uniform-stride; ties broken by RNG-permuted order.
    if len(sub) > target_n:
        idx = np.linspace(0, len(sub) - 1, target_n).round().astype(np.int64)
        idx = np.unique(idx)
        sub = sub.iloc[idx].reset_index(drop=True)

    return sub


def score_window(
    handles: dict[str, SiteHandle],
    window: pd.Series,
    cfg: SimulationConfig,
    cache_dir: Path | None,
) -> dict[str, object]:
    """Score a single window. Returns a flat dict of result fields."""
    cfg_ls = cfg.locus_selection

    # Resume from cache: if every site's r² matrix is already on disk for this
    # window, skip the genotype read + LD computation and rebuild the score
    # directly from cached arrays. Safe because the window_id and tag-SNP
    # selection are deterministic in (config, master_seed).
    if cache_dir is not None:
        cached_paths = {
            site: Path(cache_dir) / f"{window['window_id']}_{site}.npz"
            for site in handles
        }
        from ..inference import input_fingerprint
        stamp = Path(cache_dir) / f"{window['window_id']}.inputs.json"
        files = [Path(str(h.bed_prefix) + ext) for h in handles.values()
                 for ext in (".bed", ".bim", ".fam")]
        fingerprint = input_fingerprint(files, {"selection": cfg_ls.model_dump(), "seed": cfg.master_seed,
                                                "protocol": "harmonized-pooling-v2"})
        if (all(p.exists() for p in cached_paths.values()) and stamp.exists()
                and stamp.read_text() == fingerprint):
            loaded = {s: np.load(p, allow_pickle=True) for s, p in cached_paths.items()}
            per_site_r2 = {s: z["r2"] for s, z in loaded.items()}
            n_tags = next(iter(per_site_r2.values())).shape[0]
            # Recover the pre-prune variant count if the cache carries it (newer
            # caches do); fall back to -1 for legacy caches written before this.
            first = next(iter(loaded.values()))
            n_var = int(first["n_variants_in_window"]) if "n_variants_in_window" in first.files else -1
            mean_div, pair_scores = pairwise_divergence(per_site_r2)
            return {
                **window.to_dict(),
                "n_variants_in_window": n_var,
                "n_tag_snps": int(n_tags),
                "divergence_score": float(mean_div),
                "status": "ok",
            }

    variants = _select_window_variants(handles, window, cfg_ls.maf_filter)
    if variants.empty:
        return {
            **window.to_dict(),
            "n_variants_in_window": 0,
            "n_tag_snps": 0,
            "divergence_score": np.nan,
            "status": "no_variants",
        }

    tags = _pick_tag_snps(
        handles,
        variants,
        cfg_ls.ld_tag_snp_target,
        cfg_ls.ld_prune_window_variants,
        cfg_ls.ld_r2_prune_threshold,
        seed_components=(window["window_id"],),
        master_seed=cfg.master_seed,
    )
    if len(tags) < 2:
        return {
            **window.to_dict(),
            "n_variants_in_window": int(len(variants)),
            "n_tag_snps": int(len(tags)),
            "divergence_score": np.nan,
            "status": "too_few_tags",
        }

    per_site_r2: dict[str, np.ndarray] = {}
    for site_name, h in handles.items():
        vidx = tags[f"vidx_{site_name}"].to_numpy()
        dosage = read_bed_variants(h.bed_prefix, vidx, h.n_samples)
        per_site_r2[site_name] = r2_matrix(dosage)

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        for site_name, R in per_site_r2.items():
            np.savez_compressed(
                cache_dir / f"{window['window_id']}_{site_name}.npz",
                r2=R,
                snp_ids=tags["snp_id"].to_numpy(),
                n_variants_in_window=np.int64(len(variants)),
            )

        stamp.write_text(fingerprint)

    mean_div, pair_scores = pairwise_divergence(per_site_r2)
    out = {
        **window.to_dict(),
        "n_variants_in_window": int(len(variants)),
        "n_tag_snps": int(len(tags)),
        "divergence_score": mean_div,
        "status": "ok",
    }
    for pair, score in pair_scores.items():
        out[f"frob_{pair}"] = score
    return out


# ---------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------


def _score_candidate_rows(
    handles: dict,
    candidates: pd.DataFrame,
    cfg: SimulationConfig,
    cache_dir: Path,
) -> list[dict]:
    """Score the given candidate windows with joblib. Each window's per-site r²
    matrices are written to ``cache_dir`` (independent files), so this is safe to
    call from several nodes on disjoint shards of the candidate list at once."""

    def _job(row_dict: dict) -> dict:
        return score_window(handles, pd.Series(row_dict), cfg, cache_dir)

    rows = list(candidates.to_dict(orient="records"))
    n_workers = max(1, cfg.locus_selection.n_workers)
    if n_workers == 1:
        return [_job(r) for r in rows]
    # loky forks workers; each gets its own handles via closure (read-only → safe).
    return Parallel(n_jobs=n_workers, backend="loky", verbose=5)(
        delayed(_job)(r) for r in rows
    )


def _build_candidates(cfg: SimulationConfig, handles: dict) -> pd.DataFrame:
    """Generate the deterministic candidate-window list (shared by the sharded
    map step and the single-node reduce step, so both see identical windows)."""
    logger = get_logger()
    chrom = str(cfg.chromosome)
    mins, maxs = [], []
    for h in handles.values():
        on_chrom = h.bim[h.bim["chrom"].astype(str) == chrom]
        if on_chrom.empty:
            raise RuntimeError(f"Site {h.name} has no variants on chr {chrom}")
        mins.append(int(on_chrom["bp"].min()))
        maxs.append(int(on_chrom["bp"].max()))
    cfg_ls = cfg.locus_selection
    candidates = generate_candidate_windows(
        chrom=chrom,
        min_bp=min(mins),
        max_bp=max(maxs),
        window_size_bp=cfg_ls.window_size_bp,
        step_size_bp=cfg_ls.step_size_bp,
    )
    logger.info("Generated %d candidate windows", len(candidates))
    return candidates


def score_windows_shard(cfg: SimulationConfig, shard_index: int, n_shards: int) -> int:
    """MAP step for multi-node runs: score only the windows assigned to this shard
    (``i % n_shards == shard_index``), populating the shared ``ld_matrices/`` cache.
    A subsequent single-node ``run_locus_selection`` reduce reads the cache (all
    hits) and does the stratify/select. Returns the number of windows scored."""
    logger = get_logger()
    handles = _load_site_handles(cfg)
    if not handles:
        raise RuntimeError("No site handles loaded — did sampling.py run?")
    candidates = _build_candidates(cfg, handles).reset_index(drop=True)
    shard = candidates.iloc[shard_index::n_shards].copy()
    logger.info(
        "Shard %d/%d: scoring %d of %d windows",
        shard_index, n_shards, len(shard), len(candidates),
    )
    cache_dir = cfg.resolved_path("loci_dir") / "ld_matrices"
    ensure_dir(cache_dir)
    ensure_dir(cfg.resolved_path("loci_dir") / "tag_snps")
    _score_candidate_rows(handles, shard, cfg, cache_dir)
    logger.info("Shard %d/%d: done", shard_index, n_shards)
    return len(shard)


def run_locus_selection(cfg: SimulationConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the locus selection stage end-to-end.

    Returns
    -------
    (candidate_windows_df, selected_loci_df)
    """
    logger = get_logger()
    handles = _load_site_handles(cfg)
    if not handles:
        raise RuntimeError("No site handles loaded — did sampling.py run?")

    cfg_ls = cfg.locus_selection
    candidates = _build_candidates(cfg, handles)

    cache_dir = cfg.resolved_path("loci_dir") / "ld_matrices"
    ensure_dir(cache_dir)
    ensure_dir(cfg.resolved_path("loci_dir") / "tag_snps")

    scored = pd.DataFrame(_score_candidate_rows(handles, candidates, cfg, cache_dir))
    loci_dir = cfg.resolved_path("loci_dir")
    candidates_out = loci_dir / "candidate_windows.tsv"
    scored.to_csv(candidates_out, sep="\t", index=False)
    logger.info("Wrote %s (%d rows)", candidates_out, len(scored))

    ok = scored[scored["status"] == "ok"].copy()
    if len(ok) < cfg_ls.n_loci:
        raise RuntimeError(
            f"Only {len(ok)} scorable windows; need {cfg_ls.n_loci}"
        )

    rng = np.random.default_rng(derive_seed(cfg.master_seed, "stratify"))
    strata = {"low": cfg_ls.strata.low, "medium": cfg_ls.strata.medium, "high": cfg_ls.strata.high}
    selected = stratify_and_sample(ok, strata, rng)
    selected = drop_overlapping_selections(selected)

    # Stable ID for downstream code.
    selected["locus_id"] = [f"L{i:04d}" for i in range(len(selected))]
    cols = ["locus_id", "window_id", "chrom", "start_bp", "end_bp",
            "divergence_score", "stratum", "n_variants_in_window", "n_tag_snps"]
    selected = selected[[c for c in cols if c in selected.columns]]

    selected_out = loci_dir / "selected_loci.tsv"
    selected.to_csv(selected_out, sep="\t", index=False)
    logger.info("Wrote %s (%d loci)", selected_out, len(selected))

    return scored, selected


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="FedFM locus selection")
    parser.add_argument("--config", default="config/simulation_config.yaml")
    parser.add_argument("--shard-index", type=int, default=0,
                        help="0-based shard index for multi-node scoring")
    parser.add_argument("--n-shards", type=int, default=1,
                        help="total shards; >1 with --score-only runs the MAP step")
    parser.add_argument("--score-only", action="store_true",
                        help="MAP step: score this shard's windows into the cache, "
                             "then exit (no stratify/select). Run a plain (unsharded) "
                             "invocation afterwards as the reduce step.")
    args = parser.parse_args(list(argv) if argv is not None else None)

    cfg = load_config(args.config)
    setup_logging(log_dir=cfg.resolved_path("logs_dir"))

    if args.score_only:
        if args.n_shards < 1 or not (0 <= args.shard_index < args.n_shards):
            parser.error("--score-only needs 0 <= --shard-index < --n-shards")
        score_windows_shard(cfg, args.shard_index, args.n_shards)
    else:
        run_locus_selection(cfg)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
