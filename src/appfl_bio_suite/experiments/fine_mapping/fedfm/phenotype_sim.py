"""Causal architecture and phenotype simulation.

For each (locus, architecture, replicate) instance:
    1. Pick ``ncsl`` causal variants from the locus window.
       Optionally override one variant to satisfy the ancestry-specific
       MAF criterion for the stratum. When the (orthogonal) ancestry-DIVERGENT
       mode is flagged for the instance, instead pick a *union* causal set whose
       composition differs per superpopulation (shared + per-pop-private).
    2. Draw a per-superpopulation effect-size matrix β ~ MVN(0, Σ_rg).
       Σ has 1 on the diagonal and rg off-diagonal.
       In divergent mode β has structural zeros so each pop only sees its own
       causal set. Scale β so the genetic variance equals h² × Var(y).
    3. For each site, compute the per-individual genetic value using the
       β indexed by **each individual's own superpopulation**, then add
       Gaussian noise sized to hit the target per-locus h².

Outputs:
    * ``<gt>/phenotypes/<site>/locus<L>_arch<A>_rep<R>.pheno``  (PLINK FID/IID/y)
    * ``<gt>/effect_sizes/locus<L>_arch<A>_rep<R>.npy``         (ncsl × n_pops)
    * ``<gt>/causal_manifest.tsv``                               (one row per instance)
    * ``<gt>/seeds.tsv``                                         (one row per instance)
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .utils import (
    SimulationConfig,
    derive_seed,
    ensure_dir,
    get_logger,
    load_config,
    read_bed_variants,
    read_bim,
    read_fam,
    read_site_manifest,
    setup_logging,
)


# ---------------------------------------------------------------------------
# Architecture grid
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Architecture:
    architecture_id: str
    ncsl: int
    h2: float
    rg: float


def build_architecture_grid(
    ncsl: Sequence[int],
    h2: Sequence[float],
    rg: Sequence[float],
    mode: str,
    h2_anchor: float = 0.001,
    rg_anchor: float = 1.0,
) -> list[Architecture]:
    """Build the architecture grid for the chosen factorial mode.

    Modes:
        * ``minimal``:  ncsl × h² at rg=rg_anchor (3 × 3 = 9 combinations).
        * ``extended``: minimal + (ncsl × rg at h²=h2_anchor); deduped.
        * ``full``:     ncsl × h² × rg, the full Cartesian product.
    """
    combos: list[tuple[int, float, float]] = []
    if mode == "minimal":
        for n in ncsl:
            for h in h2:
                combos.append((n, h, rg_anchor))
    elif mode == "extended":
        for n in ncsl:
            for h in h2:
                combos.append((n, h, rg_anchor))
        for n in ncsl:
            for r in rg:
                combos.append((n, h2_anchor, r))
    elif mode == "full":
        for n in ncsl:
            for h in h2:
                for r in rg:
                    combos.append((n, h, r))
    else:
        raise ValueError(f"Unknown factorial mode: {mode}")

    seen = set()
    out: list[Architecture] = []
    for n, h, r in combos:
        key = (n, round(h, 6), round(r, 6))
        if key in seen:
            continue
        seen.add(key)
        aid = f"ncsl{n}_h2-{h:g}_rg{r:g}"
        out.append(Architecture(architecture_id=aid, ncsl=n, h2=h, rg=r))
    return out


# ---------------------------------------------------------------------------
# Per-site state
# ---------------------------------------------------------------------------


@dataclass
class SiteState:
    name: str
    bed_prefix: Path
    bim: pd.DataFrame                          # full chr<N> .bim
    fam: pd.DataFrame                          # FID, IID, ...
    manifest: pd.DataFrame                     # FID, IID, superpopulation
    pop_index: np.ndarray                      # per-sample integer index into superpop list
    dominant: str
    n_samples: int


def _load_site_states(cfg: SimulationConfig) -> dict[str, SiteState]:
    states: dict[str, SiteState] = {}
    pop_to_idx = {p: i for i, p in enumerate(cfg.superpopulations)}
    for site_name, site_cfg in cfg.sites.items():
        site_dir = cfg.site_dir(site_name)
        prefix = site_dir / f"{site_name}_chr{cfg.chromosome}"
        bim = read_bim(str(prefix) + ".bim")
        fam = read_fam(str(prefix) + ".fam")
        manifest = read_site_manifest(site_dir / f"{site_name}_manifest.tsv")

        # Align manifest order to .fam order (PLINK's row order in .bed).
        manifest = manifest.set_index(["FID", "IID"]).loc[
            list(zip(fam["FID"], fam["IID"]))
        ].reset_index()
        pop_index = np.array(
            [pop_to_idx[p] for p in manifest["superpopulation"]], dtype=np.int64
        )
        states[site_name] = SiteState(
            name=site_name,
            bed_prefix=prefix,
            bim=bim,
            fam=fam,
            manifest=manifest,
            pop_index=pop_index,
            dominant=site_cfg.dominant,
            n_samples=len(fam),
        )
    return states


# ---------------------------------------------------------------------------
# Causal selection
# ---------------------------------------------------------------------------


def _variants_in_locus(states: dict[str, SiteState], locus: pd.Series) -> pd.DataFrame:
    """Variants present in ALL sites within [start_bp, end_bp] on the locus chrom.

    Returns merged DataFrame with snp_id, bp, vidx_<site> columns.
    """
    chrom = str(locus["chrom"])
    start = int(locus["start_bp"])
    end = int(locus["end_bp"])

    pieces: list[pd.DataFrame] = []
    for site_name, st in states.items():
        bim = st.bim
        mask = (bim["chrom"].astype(str) == chrom) & (bim["bp"] >= start) & (bim["bp"] <= end)
        sub = bim.loc[mask, ["snp_id", "bp"]].copy()
        sub[f"vidx_{site_name}"] = bim.index[mask].to_numpy()
        pieces.append(sub)

    merged = pieces[0]
    for sub in pieces[1:]:
        merged = merged.merge(sub.drop(columns=["bp"]), on="snp_id", how="inner")
    merged.sort_values("bp", inplace=True)
    merged.reset_index(drop=True, inplace=True)
    return merged


def _per_site_maf(
    states: dict[str, SiteState],
    variants: pd.DataFrame,
) -> pd.DataFrame:
    """Return a copy of ``variants`` with a ``maf_<site>`` column per site."""
    out = variants.copy()
    for site_name, st in states.items():
        vidx = out[f"vidx_{site_name}"].to_numpy()
        dosage = read_bed_variants(st.bed_prefix, vidx, st.n_samples)
        af = np.nanmean(dosage, axis=0) / 2.0
        out[f"maf_{site_name}"] = np.minimum(af, 1.0 - af)
    return out


def _passes_maf_filter(variants: pd.DataFrame, sites: Sequence[str], threshold: float) -> np.ndarray:
    """Any-site MAF > threshold mask."""
    return np.any(
        np.column_stack([variants[f"maf_{s}"].to_numpy() for s in sites]) > threshold,
        axis=1,
    )


def _find_ancestry_specific_candidates(
    variants: pd.DataFrame,
    sites: Sequence[str],
    dominant_pop_per_site: dict[str, str],
    site_dominant_maf: dict[str, np.ndarray],
    common_threshold: float,
    rare_threshold: float,
) -> np.ndarray:
    """Mask of variants that are common in one site's dominant pop and rare in another's.

    Per the spec: common (MAF > common_threshold) in one site's dominant ancestry,
    rare (MAF < rare_threshold) in at least one other site's dominant ancestry.

    Note: we approximate each site's "dominant ancestry MAF" with that site's
    overall MAF, which is reasonable because the dominant superpop is the
    majority of the cohort in all three sites (≥50% by design).
    """
    n_var = len(variants)
    if n_var == 0:
        return np.zeros(0, dtype=bool)

    mafs = np.column_stack([site_dominant_maf[s] for s in sites])  # (n_var, n_sites)
    common_mask = mafs > common_threshold
    rare_mask = mafs < rare_threshold

    # For each variant: ∃ site i common AND ∃ site j ≠ i rare.
    any_common = common_mask.any(axis=1)
    any_rare = rare_mask.any(axis=1)
    return any_common & any_rare


def select_causal_variants(
    variants: pd.DataFrame,
    ncsl: int,
    rng: np.random.Generator,
    ancestry_specific_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, bool]:
    """Pick ``ncsl`` causal variant row indices from ``variants``.

    If ``ancestry_specific_mask`` is provided and non-empty, one of the picks is
    forced to be a variant from that mask (the ancestry-specific override).

    Returns
    -------
    (chosen_local_indices, ancestry_specific_used)
    """
    n_var = len(variants)
    if n_var < ncsl:
        raise ValueError(f"Locus has {n_var} variants but ncsl={ncsl}")

    if ancestry_specific_mask is not None and ancestry_specific_mask.any():
        candidates = np.where(ancestry_specific_mask)[0]
        anchor = int(rng.choice(candidates))
        others_pool = np.setdiff1d(np.arange(n_var), [anchor], assume_unique=False)
        rest = rng.choice(others_pool, size=ncsl - 1, replace=False) if ncsl > 1 else np.array([], dtype=np.int64)
        chosen = np.sort(np.concatenate([[anchor], rest]))
        return chosen, True

    chosen = np.sort(rng.choice(n_var, size=ncsl, replace=False))
    return chosen, False


def select_ancestry_divergent_causal_sets(
    variants: pd.DataFrame,
    ncsl: int,
    n_private_per_pop: int,
    n_pops: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, list[np.ndarray], bool]:
    """Pick a union causal set whose per-superpop composition differs.

    Each superpop's causal set has ``ncsl`` variants: ``n_shared = ncsl -
    n_private_per_pop`` variants common to all pops, plus ``n_private_per_pop``
    variants private to that pop. The union has
    ``n_shared + n_pops * n_private_per_pop`` distinct variants.

    Falls back to a single shared set (returned as ``divergent=False``) when
    ``ncsl <= n_private_per_pop`` (no shared variant would remain) or the locus
    lacks enough variants for the full union.

    Returns
    -------
    (union_local_indices_sorted, per_pop_local_index_lists, divergent_used)
        ``per_pop_local_index_lists[p]`` is the array of local indices (into
        ``variants``) that are causal for population ``p``. In the fallback case
        every pop shares the same set.
    """
    n_var = len(variants)
    n_shared = ncsl - n_private_per_pop
    n_union = n_shared + n_pops * n_private_per_pop

    if n_shared < 1 or n_union > n_var:
        # Fall back to a shared causal set (ordinary mode).
        chosen = np.sort(rng.choice(n_var, size=ncsl, replace=False))
        return chosen, [chosen for _ in range(n_pops)], False

    pool = rng.choice(n_var, size=n_union, replace=False)
    shared = pool[:n_shared]
    privates = pool[n_shared:].reshape(n_pops, n_private_per_pop)
    per_pop = [np.sort(np.concatenate([shared, privates[p]])) for p in range(n_pops)]
    union = np.sort(pool)
    return union, per_pop, True


# ---------------------------------------------------------------------------
# Effect sizes and phenotypes
# ---------------------------------------------------------------------------


def draw_effect_sizes(
    ncsl: int,
    n_pops: int,
    rg: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw (ncsl × n_pops) effect sizes from MVN(0, Σ) per variant.

    Σ has 1 on the diagonal and rg off-diagonal; same Σ across all ncsl variants.

    Returns
    -------
    np.ndarray of shape (ncsl, n_pops). Standardized (per-variant draw is
    MVN(0, Σ); no h² scaling yet — that happens after we know the per-site
    genetic variance).
    """
    Sigma = np.full((n_pops, n_pops), rg, dtype=np.float64)
    np.fill_diagonal(Sigma, 1.0)
    # Ensure SPD even when rg=1 (then Σ is rank-1; add tiny jitter for the cholesky).
    if rg >= 1.0 - 1e-9:
        Sigma = Sigma + 1e-10 * np.eye(n_pops)
    L = np.linalg.cholesky(Sigma)
    Z = rng.standard_normal(size=(ncsl, n_pops))
    return Z @ L.T


def assemble_divergent_beta(
    union_local_idx: np.ndarray,
    per_pop_local_idx: Sequence[np.ndarray],
    n_pops: int,
    rg: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Build a ``(n_union, n_pops)`` β with structural zeros for a divergent set.

    A union row that is causal in *all* pops (a shared variant) gets an
    MVN(0, Σ_rg) draw across the pop columns; a row causal in exactly one pop
    (a private variant) gets a single N(0, 1) in that pop's column and zeros
    elsewhere. Under the per-individual β lookup, individuals of pop ``p`` then
    only receive contributions from their own causal set, because β is zero for
    every variant not in it.
    """
    n_union = len(union_local_idx)
    beta = np.zeros((n_union, n_pops), dtype=np.float64)
    row_of = {int(u): i for i, u in enumerate(union_local_idx)}
    members: list[list[int]] = [[] for _ in range(n_union)]
    for p in range(n_pops):
        for u in per_pop_local_idx[p]:
            members[row_of[int(u)]].append(p)

    shared_rows = [i for i in range(n_union) if len(members[i]) == n_pops]
    private_rows = [i for i in range(n_union) if len(members[i]) == 1]

    if shared_rows:
        shared_beta = draw_effect_sizes(len(shared_rows), n_pops, rg, rng)
        for j, i in enumerate(shared_rows):
            beta[i] = shared_beta[j]
    for i in private_rows:
        beta[i, members[i][0]] = rng.standard_normal()
    return beta


def _compute_genetic_value(
    dosage: np.ndarray,
    beta: np.ndarray,
    pop_index: np.ndarray,
) -> np.ndarray:
    """Compute per-individual genetic value using each individual's own population's β.

    Parameters
    ----------
    dosage : (n_samples, ncsl) genotype matrix for the causal variants.
    beta : (ncsl, n_pops) per-population effect sizes.
    pop_index : (n_samples,) integer population index per individual.

    Returns
    -------
    g : (n_samples,) genetic values.
    """
    # beta[:, pop_index] has shape (ncsl, n_samples) — each individual gets the
    # column of β corresponding to *their* superpopulation. Multiply elementwise
    # with the dosage and sum across causal variants.
    g = (dosage * beta[:, pop_index].T).sum(axis=1)
    # Mean-impute NaN dosages so missing genotypes don't produce NaN phenotypes.
    if np.isnan(g).any():
        # Recompute with NaN→column-mean imputation.
        dosage_imp = dosage.copy()
        col_means = np.nanmean(dosage_imp, axis=0)
        inds = np.where(np.isnan(dosage_imp))
        dosage_imp[inds] = np.take(col_means, inds[1])
        g = (dosage_imp * beta[:, pop_index].T).sum(axis=1)
    return g


def simulate_phenotype_for_site(
    state: SiteState,
    causal_vidx: np.ndarray,
    beta: np.ndarray,
    target_h2: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float, float]:
    """Simulate a phenotype vector for one site / one (locus, arch, rep) instance.

    Returns
    -------
    (y, empirical_h2, genetic_var)
    """
    dosage = read_bed_variants(state.bed_prefix, causal_vidx, state.n_samples)
    g = _compute_genetic_value(dosage, beta, state.pop_index)
    g_var = float(np.var(g))
    if g_var <= 0:
        # Degenerate: all monomorphic causal vars at this site. Phenotype is pure noise.
        y = rng.standard_normal(state.n_samples)
        return y, 0.0, 0.0
    noise_var = g_var * (1.0 - target_h2) / target_h2
    eps = rng.normal(loc=0.0, scale=np.sqrt(noise_var), size=state.n_samples)
    y = g + eps
    emp_h2 = float(g_var / np.var(y))
    return y, emp_h2, g_var


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def _write_phenotype_file(state: SiteState, y: np.ndarray, path: Path) -> None:
    df = state.fam[["FID", "IID"]].copy()
    df["y"] = y
    df.to_csv(path, sep="\t", index=False)


def _assign_ancestry_specific_flags(
    selected_loci: pd.DataFrame,
    architectures: list[Architecture],
    cfg: SimulationConfig,
) -> dict[tuple[str, str], bool]:
    """Pick at least ``min_per_stratum`` (locus, arch) pairs per stratum where
    the ancestry-specific override is required.

    Returns a dict keyed by (locus_id, architecture_id) → bool.
    """
    asc = cfg.architecture.ancestry_specific_causal
    flags: dict[tuple[str, str], bool] = {}
    if not asc.enabled or not architectures:
        return flags
    rng = np.random.default_rng(derive_seed(cfg.master_seed, "ancestry_specific"))
    arch_by_id = {a.architecture_id: a for a in architectures}
    # Prefer architectures with rg<1 (most informative for federation) when picking.
    sorted_arch_ids = sorted(
        arch_by_id.keys(), key=lambda aid: (arch_by_id[aid].rg, arch_by_id[aid].ncsl)
    )
    for stratum, sub in selected_loci.groupby("stratum"):
        loci_ids = sub["locus_id"].tolist()
        picks_needed = asc.min_per_stratum
        if not loci_ids or not sorted_arch_ids:
            continue
        for _ in range(picks_needed):
            loc = loci_ids[rng.integers(len(loci_ids))]
            arch = sorted_arch_ids[rng.integers(len(sorted_arch_ids))]
            flags[(loc, arch)] = True
    return flags


def _assign_ancestry_divergent_flags(
    selected_loci: pd.DataFrame,
    architectures: list[Architecture],
    cfg: SimulationConfig,
) -> dict[tuple[str, str], bool]:
    """Pick ``min_per_stratum`` (locus, arch) pairs per stratum that use the
    ancestry-divergent causal-set mode. Only architectures with
    ``ncsl > n_private_per_pop`` are eligible (else the mode cannot leave a
    shared variant and would silently fall back).

    Like :func:`_assign_ancestry_specific_flags`, this is a GLOBAL seeded draw
    over the full locus set and MUST be computed before sharding, so every shard
    agrees with a single-node run. Returns {(locus_id, architecture_id): True}.
    """
    adc = cfg.architecture.ancestry_divergent_causal
    flags: dict[tuple[str, str], bool] = {}
    if not adc.enabled or not architectures:
        return flags
    eligible = [a for a in architectures if a.ncsl > adc.n_private_per_pop]
    if not eligible:
        get_logger().warning(
            "ancestry_divergent_causal enabled but no architecture has "
            "ncsl > n_private_per_pop=%d; no divergent instances assigned",
            adc.n_private_per_pop,
        )
        return flags
    rng = np.random.default_rng(derive_seed(cfg.master_seed, "ancestry_divergent"))
    # Prefer rg<1 architectures (most informative for federation) when picking.
    sorted_arch_ids = sorted(
        (a.architecture_id for a in eligible),
        key=lambda aid: next(a.rg for a in eligible if a.architecture_id == aid),
    )
    for _stratum, sub in selected_loci.groupby("stratum"):
        loci_ids = sub["locus_id"].tolist()
        if not loci_ids:
            continue
        for _ in range(adc.min_per_stratum):
            loc = loci_ids[rng.integers(len(loci_ids))]
            arch = sorted_arch_ids[rng.integers(len(sorted_arch_ids))]
            flags[(loc, arch)] = True
    return flags


def clean_outputs(cfg: SimulationConfig) -> None:
    """Remove previously-generated phenotype artifacts so a re-run with a
    different locus set leaves no orphans. Everything is regenerated
    deterministically. MUST run exactly once before a sharded multi-node sim
    (else shards would wipe each other's output)."""
    logger = get_logger()
    gt_dir = cfg.resolved_path("ground_truth_dir")
    for sub in ("effect_sizes", "phenotypes"):
        d = gt_dir / sub
        if d.exists():
            shutil.rmtree(d)
    for fname in ("causal_manifest.tsv", "seeds.tsv"):
        f = gt_dir / fname
        if f.exists():
            f.unlink()
    # Also clear any stale per-shard manifest parts from a previous run.
    for part in gt_dir.glob("causal_manifest.part*of*.tsv"):
        part.unlink()
    for part in gt_dir.glob("seeds.part*of*.tsv"):
        part.unlink()
    ensure_dir(gt_dir / "effect_sizes")
    for site in cfg.sites:
        ensure_dir(gt_dir / "phenotypes" / site)
    logger.info("Cleared previous phenotype outputs under %s", gt_dir)


def run_phenotype_sim(
    cfg: SimulationConfig,
    shard_index: int = 0,
    n_shards: int = 1,
) -> None:
    logger = get_logger()
    gt_dir = cfg.resolved_path("ground_truth_dir")

    sharded = n_shards > 1
    if not sharded:
        # Single-node full run owns the clean. Sharded runs rely on a one-time
        # `--clean-only` invocation having run first (see the prod PBS).
        clean_outputs(cfg)
    else:
        ensure_dir(gt_dir / "effect_sizes")
        for site in cfg.sites:
            ensure_dir(gt_dir / "phenotypes" / site)

    states = _load_site_states(cfg)
    site_names = list(states.keys())
    n_pops = len(cfg.superpopulations)

    full_loci = pd.read_csv(
        cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t"
    ).reset_index(drop=True)
    logger.info("Loaded %d selected loci", len(full_loci))

    architectures = build_architecture_grid(
        ncsl=cfg.architecture.ncsl,
        h2=cfg.architecture.h2,
        rg=cfg.architecture.rg,
        mode=cfg.architecture.factorial_mode,
    )
    logger.info(
        "Built %d architectures (mode=%s)", len(architectures), cfg.architecture.factorial_mode
    )

    # Ancestry-specific flags are a GLOBAL decision over the full locus set (a
    # single seeded RNG draws min_per_stratum (locus, arch) pairs per stratum).
    # It MUST be computed on the full set, before sharding, or each shard would
    # see a different locus subset and assign different flags than a single-node
    # run — breaking sharded/single-node equivalence.
    asc_flags = _assign_ancestry_specific_flags(full_loci, architectures, cfg)
    div_flags = _assign_ancestry_divergent_flags(full_loci, architectures, cfg)
    if div_flags:
        logger.info("Ancestry-divergent causal mode active on %d (locus,arch) pairs",
                    len(div_flags))

    if sharded:
        keep = (full_loci.index % n_shards) == shard_index
        selected_loci = full_loci.loc[keep].copy()
        logger.info(
            "Shard %d/%d: simulating %d loci", shard_index, n_shards, len(selected_loci)
        )
    else:
        selected_loci = full_loci
    n_replicates = cfg.architecture.replicates
    total_instances = len(selected_loci) * len(architectures) * n_replicates
    logger.info(
        "Total instances to simulate: %d × %d × %d = %d",
        len(selected_loci), len(architectures), n_replicates, total_instances
    )

    manifest_rows: list[dict] = []
    seed_rows: list[dict] = []
    tol = cfg.phenotype.h2_tolerance_relative

    t_start = time.time()
    instance_idx = 0
    for _, locus in selected_loci.iterrows():
        # Per-locus variant table + per-site MAF (single read across all instances).
        variants = _variants_in_locus(states, locus)
        if variants.empty:
            logger.warning("Locus %s has no shared variants — skipping", locus["locus_id"])
            continue
        variants_with_maf = _per_site_maf(states, variants)
        any_maf_mask = _passes_maf_filter(
            variants_with_maf, site_names, cfg.locus_selection.maf_filter
        )
        usable = variants_with_maf.loc[any_maf_mask].reset_index(drop=True)

        # Ancestry-specific candidate mask (using each site's overall MAF as a proxy
        # for the dominant-pop MAF; see docstring of _find_ancestry_specific_candidates).
        dominant_maf = {s: usable[f"maf_{s}"].to_numpy() for s in site_names}
        dominant_pop = {s: states[s].dominant for s in site_names}
        asc_candidates_mask = _find_ancestry_specific_candidates(
            usable,
            site_names,
            dominant_pop,
            dominant_maf,
            cfg.architecture.ancestry_specific_causal.common_maf_threshold,
            cfg.architecture.ancestry_specific_causal.rare_maf_threshold,
        )

        n_private = cfg.architecture.ancestry_divergent_causal.n_private_per_pop
        for arch in architectures:
            require_asc = asc_flags.get((locus["locus_id"], arch.architecture_id), False)
            require_div = div_flags.get((locus["locus_id"], arch.architecture_id), False)
            for rep in range(n_replicates):
                instance_idx += 1
                seed_seq = derive_seed(
                    cfg.master_seed, "instance", locus["locus_id"], arch.architecture_id, rep
                )
                rng = np.random.default_rng(seed_seq)

                # 1) Pick causal variants and 2) draw β per causal variant × pop.
                #    Divergent mode picks a per-superpop union set + masked β; it
                #    takes precedence over the (orthogonal) asc MAF override.
                per_pop_snp_ids: dict[str, list[str]] | None = None
                div_used = False
                if require_div:
                    local_idx, per_pop_local, div_used = (
                        select_ancestry_divergent_causal_sets(
                            usable, arch.ncsl, n_private, n_pops, rng
                        )
                    )
                    causal = usable.iloc[local_idx]
                    asc_used = False
                    if div_used:
                        beta = assemble_divergent_beta(
                            local_idx, per_pop_local, n_pops, arch.rg, rng
                        )
                        per_pop_snp_ids = {
                            cfg.superpopulations[p]:
                                usable.iloc[per_pop_local[p]]["snp_id"].tolist()
                            for p in range(n_pops)
                        }
                    else:
                        beta = draw_effect_sizes(arch.ncsl, n_pops, arch.rg, rng)
                else:
                    mask = asc_candidates_mask if require_asc else None
                    local_idx, asc_used = select_causal_variants(
                        usable, arch.ncsl, rng, mask
                    )
                    causal = usable.iloc[local_idx]
                    beta = draw_effect_sizes(arch.ncsl, n_pops, arch.rg, rng)

                # 3) For each site: phenotype vector. Sites get independent noise streams.
                empirical_h2_per_site: dict[str, float] = {}
                for site_name, state in states.items():
                    vidx = causal[f"vidx_{site_name}"].to_numpy()
                    noise_rng = np.random.default_rng(
                        derive_seed(
                            cfg.master_seed, "noise", site_name, locus["locus_id"],
                            arch.architecture_id, rep
                        )
                    )
                    y, emp_h2, _ = simulate_phenotype_for_site(
                        state, vidx, beta, arch.h2, noise_rng
                    )
                    empirical_h2_per_site[site_name] = emp_h2

                    out_path = (
                        gt_dir / "phenotypes" / site_name
                        / f"{locus['locus_id']}_{arch.architecture_id}_rep{rep}.pheno"
                    )
                    _write_phenotype_file(state, y, out_path)

                # 4) Save β matrix.
                beta_path = gt_dir / "effect_sizes" / (
                    f"{locus['locus_id']}_{arch.architecture_id}_rep{rep}.npy"
                )
                np.save(beta_path, beta.astype(np.float32))

                # 5) Validate empirical h² is within tolerance (per site). Log warnings
                # rather than hard-fail, since degenerate monomorphic-at-a-site cases
                # are real and recorded with empirical_h2=0.
                for site_name, emp in empirical_h2_per_site.items():
                    if emp == 0.0:
                        continue
                    rel = abs(emp - arch.h2) / arch.h2
                    if rel > tol:
                        logger.warning(
                            "h² out of tolerance: locus=%s arch=%s site=%s target=%.4f emp=%.4f rel_err=%.3f",
                            locus["locus_id"], arch.architecture_id, site_name,
                            arch.h2, emp, rel,
                        )

                # 6) Record manifest + seed rows. ``causal`` spans the union set
                #    in divergent mode (n_union rows), so iterate its actual rows.
                per_site_maf_json = json.dumps(
                    {s: float(causal[f"maf_{s}"].iloc[k])
                     for k in range(len(causal))
                     for s in site_names},
                    sort_keys=True,
                )
                manifest_rows.append({
                    "locus_id": locus["locus_id"],
                    "architecture_id": arch.architecture_id,
                    "replicate": rep,
                    "ncsl": arch.ncsl,
                    "h2_target": arch.h2,
                    "rg": arch.rg,
                    "causal_mode": "divergent" if div_used else "shared",
                    "causal_snp_ids": ",".join(causal["snp_id"].tolist()),
                    "causal_bp": ",".join(str(int(b)) for b in causal["bp"]),
                    "causal_snp_ids_by_pop_json": (
                        json.dumps(per_pop_snp_ids, sort_keys=True)
                        if per_pop_snp_ids is not None else ""
                    ),
                    "per_site_mafs_json": per_site_maf_json,
                    "ancestry_specific_required": require_asc,
                    "ancestry_specific_used": asc_used,
                    "ancestry_divergent_required": require_div,
                    "ancestry_divergent_used": div_used,
                    "empirical_h2_anl": empirical_h2_per_site.get("anl", np.nan),
                    "empirical_h2_covenant": empirical_h2_per_site.get("covenant", np.nan),
                    "empirical_h2_mbzuai": empirical_h2_per_site.get("mbzuai", np.nan),
                })
                seed_rows.append({
                    "locus_id": locus["locus_id"],
                    "architecture_id": arch.architecture_id,
                    "replicate": rep,
                    "seed_entropy": str(list(seed_seq.entropy) if isinstance(seed_seq.entropy, (list, tuple)) else seed_seq.entropy),
                })

                if instance_idx % 100 == 0 or instance_idx == total_instances:
                    elapsed = time.time() - t_start
                    rate = instance_idx / elapsed
                    eta_s = (total_instances - instance_idx) / max(rate, 1e-9)
                    logger.info(
                        "[%d/%d] %.1f inst/s, ETA %.1f min",
                        instance_idx, total_instances, rate, eta_s / 60.0
                    )

    manifest_df = pd.DataFrame(manifest_rows)
    seeds_df = pd.DataFrame(seed_rows)
    if sharded:
        tag = f"part{shard_index:03d}of{n_shards:03d}"
        manifest_df.to_csv(gt_dir / f"causal_manifest.{tag}.tsv", sep="\t", index=False)
        seeds_df.to_csv(gt_dir / f"seeds.{tag}.tsv", sep="\t", index=False)
        logger.info("Shard %d/%d wrote %d manifest rows (%s)",
                    shard_index, n_shards, len(manifest_df), tag)
    else:
        manifest_df.to_csv(gt_dir / "causal_manifest.tsv", sep="\t", index=False)
        seeds_df.to_csv(gt_dir / "seeds.tsv", sep="\t", index=False)
        logger.info("Wrote %d manifest rows", len(manifest_df))


def merge_phenotype_manifests(cfg: SimulationConfig, n_shards: int) -> None:
    """REDUCE step: concatenate the per-shard manifest/seed parts into the final
    causal_manifest.tsv / seeds.tsv (sorted deterministically), then remove the
    parts. Errors if any expected part is missing."""
    logger = get_logger()
    gt_dir = cfg.resolved_path("ground_truth_dir")
    man_parts, seed_parts = [], []
    for r in range(n_shards):
        tag = f"part{r:03d}of{n_shards:03d}"
        mp = gt_dir / f"causal_manifest.{tag}.tsv"
        sp = gt_dir / f"seeds.{tag}.tsv"
        if not mp.exists() or not sp.exists():
            raise FileNotFoundError(f"Missing manifest part for shard {r}: {mp} / {sp}")
        man_parts.append(pd.read_csv(mp, sep="\t"))
        seed_parts.append(pd.read_csv(sp, sep="\t"))
    sort_cols = ["locus_id", "architecture_id", "replicate"]
    manifest = pd.concat(man_parts, ignore_index=True).sort_values(sort_cols).reset_index(drop=True)
    seeds = pd.concat(seed_parts, ignore_index=True).sort_values(sort_cols).reset_index(drop=True)
    manifest.to_csv(gt_dir / "causal_manifest.tsv", sep="\t", index=False)
    seeds.to_csv(gt_dir / "seeds.tsv", sep="\t", index=False)
    for r in range(n_shards):
        tag = f"part{r:03d}of{n_shards:03d}"
        (gt_dir / f"causal_manifest.{tag}.tsv").unlink()
        (gt_dir / f"seeds.{tag}.tsv").unlink()
    logger.info("Merged %d shards -> %d manifest rows", n_shards, len(manifest))


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="FedFM phenotype simulation")
    parser.add_argument("--config", default="config/simulation_config.yaml")
    parser.add_argument("--shard-index", type=int, default=0,
                        help="0-based locus shard for multi-node sim")
    parser.add_argument("--n-shards", type=int, default=1,
                        help="total shards; >1 enables sharded (no-clean) mode")
    parser.add_argument("--clean-only", action="store_true",
                        help="wipe + recreate output dirs once, then exit "
                             "(run before sharded sim)")
    parser.add_argument("--merge", action="store_true",
                        help="REDUCE: merge per-shard manifest parts, then exit")
    args = parser.parse_args(list(argv) if argv is not None else None)
    cfg = load_config(args.config)
    setup_logging(log_dir=cfg.resolved_path("logs_dir"))

    if args.clean_only:
        clean_outputs(cfg)
    elif args.merge:
        merge_phenotype_manifests(cfg, args.n_shards)
    else:
        if args.n_shards > 1 and not (0 <= args.shard_index < args.n_shards):
            parser.error("need 0 <= --shard-index < --n-shards")
        run_phenotype_sim(cfg, args.shard_index, args.n_shards)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
