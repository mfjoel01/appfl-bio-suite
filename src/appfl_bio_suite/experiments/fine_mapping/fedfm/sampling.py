"""Cohort subsampling for the FedFM simulation pipeline.

Assigns HAPNEST individuals to the three site cohorts (ANL, Covenant, MBZUAI)
with **no individual reused across sites**. Sampling is stratified by
superpopulation according to per-site target compositions in the YAML config.

Outputs (per site, under ``<processed_dir>/<site>/``):
    * ``<site>_ids.txt``         — PLINK --keep file (FID, IID, no header)
    * ``<site>_manifest.tsv``    — FID, IID, superpopulation (with header)
    * ``<site>_chr<N>.{bed,bim,fam}`` — PLINK binaries extracted via --keep
                                       (written when ``extract_plink=True``)

The cross-site uniqueness invariant is asserted **before** any PLINK extraction.

Module entry point: ``main(config_path)`` — wired up by ``scripts/run_sampling.sh``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .utils import (
    SimulationConfig,
    derive_seed,
    ensure_dir,
    get_logger,
    load_config,
    read_population_manifest,
    run_plink,
    setup_logging,
    write_ids_file,
    write_site_manifest,
)


# Deterministic site ordering for reproducibility. Anything that reads the
# config should use this list — do NOT iterate `cfg.sites` directly.
SITE_ORDER: tuple[str, ...] = ("anl", "covenant", "mbzuai")


def _validate_pool(
    manifest: pd.DataFrame, sites: dict[str, dict[str, int]]
) -> None:
    """Verify the manifest has enough individuals in each superpop to satisfy demand."""
    demand: dict[str, int] = {}
    for site_name in SITE_ORDER:
        if site_name not in sites:
            continue
        for pop, count in sites[site_name].items():
            demand[pop] = demand.get(pop, 0) + count

    available = manifest["superpopulation"].value_counts().to_dict()
    for pop, needed in demand.items():
        have = available.get(pop, 0)
        if have < needed:
            raise ValueError(
                f"Insufficient individuals for superpop {pop}: need {needed}, have {have}"
            )


def assign_individuals(
    manifest: pd.DataFrame,
    sites: dict[str, dict[str, int]],
    master_seed: int,
) -> dict[str, pd.DataFrame]:
    """Assign manifest individuals to sites with no cross-site overlap.

    Parameters
    ----------
    manifest : pd.DataFrame
        Columns FID, IID, superpopulation. One row per HAPNEST individual.
    sites : dict[str, dict[str, int]]
        Mapping site_name → {superpop: count}.
    master_seed : int
        Deterministic seed anchor.

    Returns
    -------
    dict[str, pd.DataFrame]
        site_name → DataFrame (FID, IID, superpopulation).
    """
    logger = get_logger()
    _validate_pool(manifest, sites)

    # Per-superpop pool of (FID, IID) indices that have NOT been assigned yet.
    # Indexed by integer position into a per-pop array so we can sample
    # quickly with replacement-free numpy choice.
    pool: dict[str, np.ndarray] = {}
    for pop, sub in manifest.groupby("superpopulation"):
        # Sort to make ordering deterministic before we shuffle below.
        idx = np.array(sub.index, dtype=np.int64)
        idx.sort()
        pool[pop] = idx

    # Shuffle each pop's pool once with a deterministic, pop-keyed seed so the
    # "remaining" set after one site doesn't have a structural bias for the next.
    for pop in sorted(pool):
        rng = np.random.default_rng(derive_seed(master_seed, "shuffle", pop))
        rng.shuffle(pool[pop])

    assignments: dict[str, pd.DataFrame] = {}
    # Track absolute manifest indices used anywhere — for the cross-site invariant.
    used_indices: set[int] = set()

    for site_name in SITE_ORDER:
        if site_name not in sites:
            continue
        composition = sites[site_name]
        chosen_indices: list[int] = []
        for pop in sorted(composition):  # sort for determinism
            count = composition[pop]
            available_idx = pool[pop]
            if len(available_idx) < count:
                raise RuntimeError(
                    f"Pool exhausted for {pop}: need {count}, "
                    f"have {len(available_idx)} after prior site assignments"
                )
            # Take the first `count` from the (already-shuffled) pool — this is
            # uniform without replacement and reproducible.
            take = available_idx[:count].tolist()
            chosen_indices.extend(take)
            # Shrink the pool in place.
            pool[pop] = available_idx[count:]

        # Build the site DataFrame and check no overlap with prior sites.
        overlap = used_indices.intersection(chosen_indices)
        if overlap:
            raise RuntimeError(
                f"Cross-site overlap detected for site {site_name}: {len(overlap)} ids"
            )
        used_indices.update(chosen_indices)

        site_df = manifest.loc[chosen_indices, ["FID", "IID", "superpopulation"]].copy()
        # Sort for stable downstream output.
        site_df.sort_values(["superpopulation", "FID", "IID"], inplace=True)
        site_df.reset_index(drop=True, inplace=True)

        # Per-pop count assertion.
        counts = site_df["superpopulation"].value_counts().to_dict()
        for pop, expected in composition.items():
            actual = counts.get(pop, 0)
            if actual != expected:
                raise RuntimeError(
                    f"Site {site_name} pop {pop}: expected {expected}, got {actual}"
                )
        if len(site_df) != sum(composition.values()):
            raise RuntimeError(
                f"Site {site_name} total mismatch: got {len(site_df)}, "
                f"expected {sum(composition.values())}"
            )

        assignments[site_name] = site_df
        logger.info(
            "Site %s: assigned %d individuals (%s)",
            site_name,
            len(site_df),
            ", ".join(f"{pop}={counts.get(pop, 0)}" for pop in sorted(composition)),
        )

    # Final global no-overlap check.
    all_pairs = pd.concat(
        [df.assign(site=name) for name, df in assignments.items()], ignore_index=True
    )
    dupes = all_pairs.duplicated(subset=["FID", "IID"], keep=False)
    if dupes.any():
        n_dupes = int(dupes.sum())
        raise RuntimeError(
            f"Cross-site uniqueness violated: {n_dupes} duplicate (FID, IID) rows"
        )

    return assignments


def _extract_plink_for_site(
    cfg: SimulationConfig,
    site_name: str,
    ids_file: Path,
    out_prefix: Path,
) -> None:
    """Call PLINK --keep to materialize the site's chr<chromosome> binary."""
    logger = get_logger()
    hapnest_dir = cfg.resolved_path("hapnest_dir")
    chrom = cfg.chromosome
    # Convention: HAPNEST chromosome files are <hapnest_dir>/chr<N>.{bed,bim,fam}.
    # If the user has different naming, they can adapt this single line.
    bfile = hapnest_dir / f"chr{chrom}"
    if not (bfile.with_suffix(".bed")).exists():
        raise FileNotFoundError(
            f"Expected HAPNEST PLINK binary at {bfile}.bed — adapt sampling.py "
            f"_extract_plink_for_site if your naming differs"
        )
    logger.info(
        "PLINK --keep extract: site=%s bfile=%s out=%s", site_name, bfile, out_prefix
    )
    run_plink(
        [
            "--bfile",
            str(bfile),
            "--keep",
            str(ids_file),
            # WITHOUT THIS FLAG THE SITES DISAGREE ABOUT WHAT A DOSAGE MEANS.
            # PLINK 1.9's --make-bed re-derives A1 as the minor allele *of the sample it
            # is given*. Each site is a different subset of individuals with a different
            # ancestry mix, so a variant whose minor allele differs between cohorts is
            # written A1=G in one site's .bim and A1=A in another's. The .bed dosage is
            # then "count of G" at one site and "count of A" at the other -- x against
            # 2 - x -- with no error and no warning.
            #
            # It corrupts the ground truth, not just the analysis: phenotype_sim reads
            # these per-site filesets, so a causal variant flipped at one site gives that
            # site's individuals beta * (2 - x) while the manifest records beta * x.
            # Every downstream number inherits it.
            "--keep-allele-order",
            "--make-bed",
            "--out",
            str(out_prefix),
        ],
        binary=cfg.tools.plink,
    )


def run_sampling(cfg: SimulationConfig, extract_plink: bool = True) -> dict[str, Path]:
    """Run the full sampling stage.

    Writes ID files + manifests per site, asserts the no-replacement invariant,
    and (optionally) calls PLINK to extract the per-site PLINK binaries for the
    target chromosome.

    Returns
    -------
    dict[str, Path]
        site_name → path to its ID file.
    """
    logger = get_logger()
    manifest_path = cfg.resolved_path("population_manifest")
    logger.info("Loading population manifest: %s", manifest_path)
    manifest = read_population_manifest(manifest_path)
    logger.info(
        "Loaded %d individuals across %d superpopulations",
        len(manifest),
        manifest["superpopulation"].nunique(),
    )

    # Convert pydantic SiteConfig → plain dict so the assigner stays library-agnostic.
    sites = {
        name: dict(cfg.sites[name].composition)
        for name in SITE_ORDER
        if name in cfg.sites
    }

    assignments = assign_individuals(manifest, sites, cfg.master_seed)

    id_files: dict[str, Path] = {}
    for site_name, site_df in assignments.items():
        site_dir = ensure_dir(cfg.site_dir(site_name))
        ids_file = site_dir / f"{site_name}_ids.txt"
        manifest_file = site_dir / f"{site_name}_manifest.tsv"
        write_ids_file(site_df, ids_file)
        write_site_manifest(site_df, manifest_file)
        id_files[site_name] = ids_file
        logger.info("Wrote %s and %s", ids_file, manifest_file)

    if extract_plink:
        for site_name in SITE_ORDER:
            if site_name not in assignments:
                continue
            site_dir = cfg.site_dir(site_name)
            out_prefix = site_dir / f"{site_name}_chr{cfg.chromosome}"
            _extract_plink_for_site(cfg, site_name, id_files[site_name], out_prefix)

    return id_files


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="FedFM cohort subsampling")
    parser.add_argument(
        "--config",
        default="config/simulation_config.yaml",
        help="Path to YAML config",
    )
    parser.add_argument(
        "--no-plink",
        action="store_true",
        help="Skip the PLINK extract step (only write ID files + manifests)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    cfg = load_config(args.config)
    setup_logging(log_dir=cfg.resolved_path("logs_dir"))
    run_sampling(cfg, extract_plink=not args.no_plink)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
