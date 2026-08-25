"""Step 2: split the pooled cohort into disjoint per-site datasets.

Each site receives a complete, self-contained dataset: its own genotype subset plus the
matching phenotype, covariate and score rows. Sites share the full variant set -- which
meta-analysis requires -- but no samples, so the federation is genuinely partitioned.

=============================================================================
NUMERICAL BEHAVIOUR IS FROZEN. Verified against the pre-migration pipeline.
=============================================================================

The single most dangerous line in this file:

    np.random.seed(split_seed)
    np.random.shuffle(indices)

That is numpy's LEGACY global generator, and it must stay that way. The obvious
modernization --

    rng = np.random.default_rng(split_seed)
    rng.shuffle(indices)

-- is a different algorithm (PCG64 rather than Mersenne Twister). It produces a
completely different permutation, and therefore a completely different assignment of
individuals to sites, with no error and nothing visibly wrong. Every per-site result
would change while every summary statistic still looked plausible.

``tests/test_gwas_simulation_parity.py`` catches exactly this.

WHAT CHANGED FROM THE ORIGINAL
------------------------------
The site sizes. They were five integers inline --

    SITE_SIZES = [18032, 9237, 25028, 40752, 6951]
    assert sum(SITE_SIZES) == 100000, "Site sizes must sum to 100,000"

-- with an assertion that welded the pipeline to one cohort size and one split. They are
now a scenario config, so "five sites, this skew" and "two sites, balanced" are
configuration rather than a source edit. The published split is the default, so the
published behaviour is what you get by not choosing.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ["bundle_sites", "SITE_DATA_FILES", "REQUIRED_BY_LOADER"]

log = logging.getLogger(__name__)

# The genotype stem every site's files carry. Kept from the original cohort's name so
# that already-distributed bundles and the client loader continue to agree.
PLINK_STEM = "EUR.synthetic.100k.ld.maf"

# Everything written into a site's data directory.
SITE_DATA_FILES = (
    f"{PLINK_STEM}.bed",
    f"{PLINK_STEM}.bim",
    f"{PLINK_STEM}.fam",
    "phenotypes_gwas.csv",
    "phenotypes_pgs_eval.csv",
    "covariates.csv",
    "pgs_scores.csv",
    "score_t2d.txt",
    "score_bmi.txt",
)

# The subset the client-side loader actually validates and reads. The rest are simulation
# and central-baseline artifacts, useful to a coordinator and inert at a partner site.
#
# tests/test_gwas_loader_contract.py asserts this list matches what the loader requires.
# Keeping the two in agreement used to be a manual job across two directories, and a
# mismatch surfaced on a partner's cluster rather than in CI.
REQUIRED_BY_LOADER = (
    f"{PLINK_STEM}.bed",
    f"{PLINK_STEM}.bim",
    f"{PLINK_STEM}.fam",
    "phenotypes_gwas.csv",
    "phenotypes_pgs_eval.csv",
    "covariates.csv",
)


def _ext(prefix: Path, extension: str) -> Path:
    """Append an extension to a PLINK prefix. Not with_suffix() -- the stem ends '.maf'."""
    return Path(str(prefix) + extension)


def _subset_rows_by_iid(df: pd.DataFrame, ordered_iids: pd.Series) -> pd.DataFrame:
    """Take one site's rows, in the genotype file's sample order.

    Row order must follow the genotypes, not the original table: the trainer merges these
    on FID/IID but the genotype matrix is positional, and a mismatch would silently pair
    each individual's genotypes with someone else's phenotype.
    """
    return df.set_index("IID", drop=False).loc[ordered_iids].reset_index(drop=True)


def bundle_sites(
    plink_prefix: Path,
    pooled_dir: Path,
    sites_dir: Path,
    site_sizes: list[int],
    split_seed: int = 42,
    site_ids: list[str] | None = None,
    data_sim_scaling: float = 1.0,
    data_sim_scaling_seed: int = 7,
) -> dict[str, Path]:
    """Split the cohort and write each site's dataset. Returns {site_id: data_dir}.

    Args:
        plink_prefix: pooled PLINK1 triple prefix, WITHOUT extension.
        pooled_dir: where step 1 wrote its outputs.
        sites_dir: parent directory for per-site data.
        site_sizes: samples per site, in order.
        split_seed: seed for the sample shuffle. Uses numpy's LEGACY global RNG --
            see this module's docstring before changing anything about it.
    """
    from pandas_plink import read_plink1_bin, write_plink1_bin

    from appfl_bio_suite.experiments.gwas.simulation.phenotypes import apply_data_sim_scaling

    plink_prefix = Path(plink_prefix)
    pooled_dir = Path(pooled_dir)
    sites_dir = Path(sites_dir)

    site_ids = site_ids or [f"Site{i}" for i in range(1, len(site_sizes) + 1)]
    if len(site_ids) != len(site_sizes):
        raise ValueError(
            f"{len(site_ids)} site ids but {len(site_sizes)} sizes -- they must correspond."
        )

    required = [
        _ext(plink_prefix, ".bed"),
        _ext(plink_prefix, ".bim"),
        _ext(plink_prefix, ".fam"),
        pooled_dir / "phenotypes_gwas.csv",
        pooled_dir / "phenotypes_pgs_eval.csv",
        pooled_dir / "covariates.csv",
        pooled_dir / "pgs_scores.csv",
        pooled_dir / "score_t2d.txt",
        pooled_dir / "score_bmi.txt",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "cannot bundle -- step 1 outputs are missing:\n  - " + "\n  - ".join(missing)
            + "\n\nRun the phenotype simulation first."
        )

    log.info("loading pooled genotypes")
    G_full = read_plink1_bin(str(_ext(plink_prefix, ".bed")), verbose=False)
    n_variants_orig = G_full.sizes["variant"]
    variant_df = pd.DataFrame({"SNP": G_full.snp.values})
    variant_df, G_full = apply_data_sim_scaling(
        variant_df, G_full, fraction=data_sim_scaling, seed=data_sim_scaling_seed
    )
    if data_sim_scaling < 1.0:
        log.info(
            "  data_sim_scaling=%s: variants reduced to %s / %s",
            data_sim_scaling,
            f"{len(variant_df):,}",
            f"{n_variants_orig:,}",
        )

    fam_full = pd.DataFrame(
        {"FID": G_full.fid.values.astype(str), "IID": G_full.iid.values.astype(str)}
    )
    total_available = len(fam_full)
    requested = sum(site_sizes)
    if requested > total_available:
        raise ValueError(
            f"site sizes sum to {requested:,} but the cohort has {total_available:,} "
            "samples. Lower the site sizes or generate a larger cohort."
        )

    log.info("pooled genotypes: %s samples x %s variants", *G_full.shape)

    pheno_gwas = pd.read_csv(pooled_dir / "phenotypes_gwas.csv")
    pheno_pgs_eval = pd.read_csv(pooled_dir / "phenotypes_pgs_eval.csv")
    covariates = pd.read_csv(pooled_dir / "covariates.csv")
    pgs_scores = pd.read_csv(pooled_dir / "pgs_scores.csv")

    # ---------------------------------------------------------------------
    # LEGACY GLOBAL RNG. Do not modernize -- see the module docstring.
    # np.random.default_rng(seed).shuffle() is a different algorithm and would
    # silently produce a different assignment of individuals to sites.
    # ---------------------------------------------------------------------
    np.random.seed(split_seed)
    indices = np.arange(len(fam_full))
    np.random.shuffle(indices)

    written: dict[str, Path] = {}
    start = 0
    for site_id, size in zip(site_ids, site_sizes, strict=True):
        log.info("writing %s (%s samples)", site_id, f"{size:,}")
        site_indices = indices[start : start + size]

        site_dir = sites_dir / site_id / "data"
        site_dir.mkdir(parents=True, exist_ok=True)

        G_sub = G_full.isel(sample=site_indices)
        out_prefix = site_dir / PLINK_STEM
        write_plink1_bin(G_sub, str(_ext(out_prefix, ".bed")), verbose=False)

        site_iids = pd.Series(G_sub.iid.values.astype(str), name="IID")

        _subset_rows_by_iid(pheno_gwas, site_iids).to_csv(
            site_dir / "phenotypes_gwas.csv", index=False
        )
        _subset_rows_by_iid(pheno_pgs_eval, site_iids).to_csv(
            site_dir / "phenotypes_pgs_eval.csv", index=False
        )
        _subset_rows_by_iid(covariates, site_iids).to_csv(
            site_dir / "covariates.csv", index=False
        )
        _subset_rows_by_iid(pgs_scores, site_iids).to_csv(
            site_dir / "pgs_scores.csv", index=False
        )

        # Scoring weights are cohort-wide, not per-site.
        shutil.copy(pooled_dir / "score_t2d.txt", site_dir / "score_t2d.txt")
        shutil.copy(pooled_dir / "score_bmi.txt", site_dir / "score_bmi.txt")

        written[site_id] = site_dir
        start += size

    log.info("bundling complete: %s site(s) under %s", len(written), sites_dir)
    return written


def verify_disjoint(site_dirs: dict[str, Path]) -> dict[tuple[str, str], int]:
    """Confirm no individual appears at two sites. Returns any overlaps found.

    Cheap, and worth doing: a partition bug would make the federation train twice on the
    same people while every per-site number still looked reasonable, and the meta-analysis
    would report a sample size it does not have.
    """
    members = {}
    for site_id, data_dir in site_dirs.items():
        fam = pd.read_csv(_ext(data_dir / PLINK_STEM, ".fam"), sep=r"\s+", header=None, dtype=str)
        members[site_id] = set(fam.iloc[:, 1])

    overlaps: dict[tuple[str, str], int] = {}
    ids = sorted(members)
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            shared = members[a] & members[b]
            if shared:
                overlaps[(a, b)] = len(shared)
    return overlaps
