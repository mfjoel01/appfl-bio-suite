"""Step 3: turn the pipeline's per-site outputs into the bundles partners receive.

The vendored pipeline writes its results in the layout the standalone repository reads:
site filesets under ``processed/<site>/``, loci under ``loci/``, phenotypes under
``ground_truth/phenotypes/<site>/``, and the answer key beside them. That layout assumes
one machine can see all of it, which is exactly the assumption a federation removes.

This module rearranges it into one self-contained directory per site, in the layout
``dataset.py`` validates on arrival::

    <out>/<site>/data/
    ├── site_genotypes.{bed,bim,fam}   the site's own cohort, nobody else's
    ├── site_manifest.tsv              its individuals' ancestry labels
    ├── reference_variants.tsv         the agreed canonical allele coding
    ├── selected_loci.tsv              which windows this study fine-maps
    └── phenotypes/<instance>.pheno    one file per locus x architecture x replicate

WHAT IS DELIBERATELY NOT IN A BUNDLE
------------------------------------
``causal_manifest.tsv`` -- the ground-truth causal variants, the whole answer key. It
stays at the coordinator. A site that held it could score its own credible sets, and more
to the point a benchmark whose answers travel with its inputs is not measuring what it
claims to. The aggregator reads it from the coordinator's copy; see its
``causal_manifest`` setting.

Also not in a bundle: any other site's anything. :func:`verify_disjoint` checks that
before the bundles are considered finished, because the premise of the whole federation
is that no individual appears twice -- and a duplicated individual would be counted twice
by the pooling while still being reported once in ``n``.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pandas as pd

__all__ = ["bundle_sites", "verify_disjoint", "write_reference_variants", "BUNDLE_SUBDIR"]

log = logging.getLogger(__name__)

# Site data sits one level down, so `<out>/<site>/` can also hold per-site artifacts that
# are NOT sent -- QC reports, for instance. The loader is handed `<out>/<site>/data`.
BUNDLE_SUBDIR = "data"


def write_reference_variants(pool_prefix: Path, destination: Path) -> Path:
    """Write the canonical variant list every site harmonizes against.

    This is *reference metadata*: variant id, position, and the A1/A2 coding that the
    dosages must count. It is the same annotation table a public panel or a consortium
    protocol distributes before any data is touched, and it says nothing about any
    individual.

    It is not optional bookkeeping. The per-site filesets are cut with
    ``plink --keep --make-bed`` and *without* ``--keep-allele-order``, so PLINK 1.9 sets
    A1 to each site's own minor allele -- on real chr1 data about 6% of variants end up
    coded oppositely at two sites. A site that summed its Gram without harmonizing first
    would contribute ``2 - x`` where the others contribute ``x``, which survives addition
    silently and negates every off-diagonal that variant touches.

    Taken here from the source pool's ``.bim``, which is what upstream's
    ``reference_variants()`` does. A real deployment would take it from the panel.
    """
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import read_bim

    bim = read_bim(Path(str(pool_prefix) + ".bim"))
    if bim["snp_id"].duplicated().any():
        n = int(bim["snp_id"].duplicated().sum())
        raise RuntimeError(
            f"the source pool's .bim has {n} duplicate variant id(s). The reference list "
            "must be keyed by variant id, so duplicates make harmonization ambiguous."
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    bim[["snp_id", "chrom", "bp", "a1", "a2"]].to_csv(destination, sep="\t", index=False)
    return destination


def bundle_sites(cfg, site_ids, pool_prefix: Path, out_dir: Path) -> dict[str, Path]:
    """Assemble one bundle per site. Returns ``{site_id: <its data directory>}``.

    ``cfg`` is the vendored :class:`SimulationConfig`, which knows where the pipeline put
    everything.
    """
    out_dir = Path(out_dir)
    loci_source = cfg.resolved_path("loci_dir") / "selected_loci.tsv"
    if not loci_source.is_file():
        raise FileNotFoundError(
            f"{loci_source} does not exist -- the locus-selection stage did not run or "
            "did not finish."
        )

    # Written once and copied, rather than regenerated per site: it is identical for
    # everyone by definition, and a per-site regeneration is a per-site chance to differ.
    reference = write_reference_variants(pool_prefix, out_dir / "reference_variants.tsv")

    bundles: dict[str, Path] = {}
    for site in site_ids:
        data_dir = out_dir / site / BUNDLE_SUBDIR
        data_dir.mkdir(parents=True, exist_ok=True)

        site_prefix = cfg.site_dir(site) / f"{site}_chr{cfg.chromosome}"
        for extension in (".bed", ".bim", ".fam"):
            source = Path(str(site_prefix) + extension)
            if not source.is_file():
                raise FileNotFoundError(
                    f"{source} does not exist -- the sampling stage did not produce a "
                    f"fileset for site '{site}'."
                )
            # Hard-link where the filesystem allows it. A site fileset is gigabytes and
            # copying it doubles the run's disk for no benefit; the bundle is read-only
            # from here on. Falls back to a copy across devices.
            _link_or_copy(source, data_dir / f"site_genotypes{extension}")

        manifest_source = cfg.site_dir(site) / f"{site}_manifest.tsv"
        if not manifest_source.is_file():
            raise FileNotFoundError(f"{manifest_source} does not exist")
        shutil.copy2(manifest_source, data_dir / "site_manifest.tsv")

        shutil.copy2(loci_source, data_dir / "selected_loci.tsv")
        _link_or_copy(reference, data_dir / "reference_variants.tsv")

        pheno_source = cfg.resolved_path("ground_truth_dir") / "phenotypes" / site
        pheno_dest = data_dir / "phenotypes"
        if not pheno_source.is_dir():
            raise FileNotFoundError(
                f"{pheno_source} does not exist -- the phenotype stage did not produce "
                f"phenotypes for site '{site}'."
            )
        pheno_dest.mkdir(parents=True, exist_ok=True)
        n_pheno = 0
        for source in sorted(pheno_source.glob("*.pheno")):
            _link_or_copy(source, pheno_dest / source.name)
            n_pheno += 1
        if n_pheno == 0:
            raise FileNotFoundError(f"{pheno_source} contains no .pheno files")

        bundles[site] = data_dir
        log.info("bundled %s: %s (%d phenotype instance(s))", site, data_dir, n_pheno)

    overlaps = verify_disjoint(bundles)
    if overlaps:
        raise RuntimeError(
            f"sites share individuals, which breaks the premise of the federation: "
            f"{overlaps}.\n"
            "Every individual must appear at exactly one site, or the pooled moments "
            "count them twice while `n` reports them once -- which biases every "
            "standardized quantity downstream without any check firing."
        )
    return bundles


def _link_or_copy(source: Path, destination: Path) -> None:
    if destination.exists():
        destination.unlink()
    try:
        destination.hardlink_to(source)
    except (OSError, AttributeError):
        shutil.copy2(source, destination)


def verify_disjoint(bundles: dict[str, Path]) -> dict[tuple[str, str], int]:
    """Return ``{(site_a, site_b): n_shared}`` for every pair that shares individuals."""
    members: dict[str, set[str]] = {}
    for site, data_dir in bundles.items():
        manifest = pd.read_csv(
            data_dir / "site_manifest.tsv", sep="\t", dtype={"IID": str}
        )
        members[site] = set(manifest["IID"])

    overlaps: dict[tuple[str, str], int] = {}
    names = sorted(members)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            shared = members[a] & members[b]
            if shared:
                overlaps[(a, b)] = len(shared)
    return overlaps
