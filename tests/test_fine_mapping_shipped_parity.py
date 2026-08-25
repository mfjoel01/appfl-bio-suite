"""The shipped site stage must equal the vendored one, exactly.

WHY THIS TEST IS THE LOAD-BEARING ONE FOR THIS EXPERIMENT
----------------------------------------------------------
The fine-mapping site stage exists twice, deliberately:

* ``fedfm/fed_fine_mapping.py::site_aggregates`` -- vendored verbatim from the standalone
  repository, and the implementation whose numbers the paper's proofs and tests are about.
* ``trainer.py::SiteFineMappingTrainer`` -- a self-contained re-implementation, because
  shipped source may not import from ``appfl_bio_suite`` and so cannot call the first one.

Duplicated statistics normally rot. This one does not get to, because the whole claim of
the experiment rests on it: algo.md Corollary 1 says the federated fit *equals* the
centralized fit rather than approximating it, and that equality is only inherited by the
APPFL path if the APPFL path computes the same aggregates. A trainer that was merely
close would turn an exact result into an approximate one, silently, and the credible sets
would still look plausible.

So the comparison here is at zero tolerance. Not ``allclose``: ``array_equal``. Every
number the site emits -- ``G``, ``u``, ``c``, ``q``, ``w``, ``n`` -- must be bit-identical
between the two implementations, along with the variant list they are indexed by and the
count of variants dropped as incompletely observed.

Bit-identity is achievable rather than aspirational because both sides do the same
operations in the same order: the same float32 Gram while ``4n`` stays under float32's
exact-integer ceiling, and the same ``_ROW_CHUNK`` when accumulating ``X'y`` in float64.
Floating-point addition is not associative, so the chunk size is part of the answer. If
this test starts failing on the ``c`` values alone, look there first.

THE FIXTURE
-----------
A miniature three-site package built in a temp directory: two ancestries, one locus
window, two phenotype instances, and -- importantly -- a site whose reference alleles are
flipped and a variant that is missing in one block, so the harmonization and complete-case
paths are both exercised rather than skipped.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.dataset import get_dataset
from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import (
    SiteIndex,
    reference_variants,
    site_aggregates,
)
from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import SimulationConfig
from appfl_bio_suite.experiments.fine_mapping.simulation.bundler import (
    write_reference_variants,
)
from appfl_bio_suite.experiments.fine_mapping.simulation.cohort import write_plink1
from appfl_bio_suite.experiments.fine_mapping.trainer import (
    KEY_GRAM,
    KEY_MANIFEST,
    KEY_SCALARS,
    KEY_SEP,
    KEY_USUM,
    KEY_XTY,
    SiteFineMappingTrainer,
)

CHROM = 1
START, END = 1_000_000, 1_400_000
M = 40  # variants in the window
POPS = ["EUR", "AFR"]
INSTANCES = ["L0000_ncsl1_h2-0.005_rg1_rep0", "L0000_ncsl2_h2-0.005_rg1_rep0"]

# Two ancestries at each site, unequal, so pooling has something to do. `covenant` is the
# site whose alleles get flipped relative to the reference.
SITES = {
    "anl": {"EUR": 40, "AFR": 15},
    "covenant": {"AFR": 35, "EUR": 10},
    "mbzuai": {"EUR": 20, "AFR": 20},
}
FLIPPED_SITE = "covenant"
FLIPPED_VARIANTS = (3, 11, 12)  # indices into the window
MISSING_CELL = (2, 7)  # (individual index within anl, variant index)


def _snp(j: int) -> str:
    return f"rs{j + 1}"


@pytest.fixture(scope="module")
def package(tmp_path_factory: pytest.TempPathFactory):
    """Build a miniature package in BOTH layouts, from one set of genotypes.

    The vendored code reads ``processed/<site>/<site>_chr1.*`` plus
    ``ground_truth/phenotypes/<site>/``; the shipped loader reads ``<site>/data/`` with
    fixed filenames. Writing both from the same arrays is what makes the comparison a
    comparison of the two implementations rather than of two datasets.
    """
    root = tmp_path_factory.mktemp("fm-parity")
    rng = np.random.default_rng(20260820)

    total = sum(sum(comp.values()) for comp in SITES.values())
    # A few extra variants outside the window, so the window selection has to do work.
    n_variants = M + 10
    positions = START - 5 * 1000 + np.arange(n_variants) * 10_000
    genotypes = rng.integers(0, 3, size=(total, n_variants)).astype(np.int8)

    iids = [f"IND{i:05d}" for i in range(total)]
    labels: list[str] = []
    site_of: list[str] = []
    for site, comp in SITES.items():
        for pop, count in comp.items():
            labels += [pop] * count
            site_of += [site] * count

    # The pool: what the reference variant list is cut from.
    pool_dir = root / "raw" / "hapnest"
    pool_dir.mkdir(parents=True)
    write_plink1(pool_dir / f"chr{CHROM}", genotypes, CHROM, positions, iids)
    with (pool_dir / "population_manifest.tsv").open("w") as handle:
        handle.write("FID\tIID\tsuperpopulation\n")
        for iid, pop in zip(iids, labels, strict=True):
            handle.write(f"{iid}\t{iid}\t{pop}\n")

    reference = write_reference_variants(pool_dir / f"chr{CHROM}", root / "reference_variants.tsv")

    loci = pd.DataFrame(
        [
            {
                "locus_id": "L0000",
                "chrom": CHROM,
                "start_bp": START,
                "end_bp": END,
                "stratum": "high",
            }
        ]
    )
    (root / "loci").mkdir()
    loci.to_csv(root / "loci" / "selected_loci.tsv", sep="\t", index=False)

    # Phenotypes: one column per instance, over every individual.
    phenotypes = {inst: rng.normal(size=total) for inst in INSTANCES}

    site_index = {site: [i for i, s in enumerate(site_of) if s == site] for site in SITES}
    for site, rows in site_index.items():
        site_geno = genotypes[rows].copy()
        alleles = [("A", "G")] * n_variants
        if site == FLIPPED_SITE:
            # PLINK 1.9 without --keep-allele-order sets A1 to the site's own minor
            # allele, so a site can legitimately hold the opposite coding. Both
            # implementations must recode it back before forming a single moment.
            for j in FLIPPED_VARIANTS:
                site_geno[:, j] = 2 - site_geno[:, j]
                alleles[j] = ("G", "A")
        if site == "anl":
            site_geno[MISSING_CELL[0], MISSING_CELL[1]] = -1

        site_iids = [iids[i] for i in rows]

        # -- upstream layout ------------------------------------------------
        processed = root / "processed" / site
        processed.mkdir(parents=True)
        write_plink1(
            processed / f"{site}_chr{CHROM}",
            site_geno,
            CHROM,
            positions,
            site_iids,
            alleles,
        )
        manifest = pd.DataFrame(
            {
                "FID": site_iids,
                "IID": site_iids,
                "superpopulation": [labels[i] for i in rows],
            }
        )
        manifest.to_csv(processed / f"{site}_manifest.tsv", sep="\t", index=False)

        pheno_dir = root / "ground_truth" / "phenotypes" / site
        pheno_dir.mkdir(parents=True)
        for inst, values in phenotypes.items():
            pd.DataFrame(
                {
                    "FID": site_iids,
                    "IID": site_iids,
                    "y": values[rows],
                }
            ).to_csv(pheno_dir / f"{inst}.pheno", sep="\t", index=False)

        # -- bundle layout --------------------------------------------------
        bundle = root / site / "data"
        bundle.mkdir(parents=True)
        write_plink1(bundle / "site_genotypes", site_geno, CHROM, positions, site_iids, alleles)
        manifest.to_csv(bundle / "site_manifest.tsv", sep="\t", index=False)
        loci.to_csv(bundle / "selected_loci.tsv", sep="\t", index=False)
        (bundle / "reference_variants.tsv").write_text(
            reference.read_text(encoding="utf-8"), encoding="utf-8"
        )
        bundle_pheno = bundle / "phenotypes"
        bundle_pheno.mkdir()
        for inst, values in phenotypes.items():
            pd.DataFrame(
                {
                    "FID": site_iids,
                    "IID": site_iids,
                    "y": values[rows],
                }
            ).to_csv(bundle_pheno / f"{inst}.pheno", sep="\t", index=False)

    cfg = SimulationConfig(
        repo_root=root,
        paths={
            "hapnest_dir": str(pool_dir),
            "population_manifest": str(pool_dir / "population_manifest.tsv"),
            "processed_dir": str(root / "processed"),
            "loci_dir": str(root / "loci"),
            "ground_truth_dir": str(root / "ground_truth"),
            "reports_dir": str(root / "reports"),
            "logs_dir": str(root / "logs"),
        },
        master_seed=1,
        chromosome=CHROM,
        tools={"plink": "plink", "plink2": "plink2", "king": "king"},
        sites={
            site: {
                "n": sum(comp.values()),
                "composition": comp,
                "dominant": max(comp, key=comp.get),
            }
            for site, comp in SITES.items()
        },
        superpopulations=POPS,
        locus_selection={
            "window_size_bp": END - START + 1,
            "step_size_bp": 100_000,
            "n_loci": 3,
            "strata": {"low": 1, "medium": 1, "high": 1},
            "maf_filter": 0.01,
            "ld_tag_snp_target": 10,
            "ld_r2_prune_threshold": 0.995,
            "ld_prune_window_variants": 50,
            "ld_prune_step_variants": 5,
            "n_workers": 1,
        },
        architecture={
            "ncsl": [1],
            "h2": [0.005],
            "rg": [1.0],
            "factorial_mode": "minimal",
            "replicates": 1,
            "ancestry_specific_causal": {
                "enabled": False,
                "min_per_stratum": 0,
                "common_maf_threshold": 0.05,
                "rare_maf_threshold": 0.01,
            },
        },
        phenotype={"model": "linear_additive", "h2_tolerance_relative": 0.5},
        qc={
            "kinship_threshold": 0.05,
            "pca_n_components": 2,
            "pca_reference": "1kg_phase3",
            "maf_tolerance_sd": 3.0,
            "ld_decay_max_dist_kb": 1000,
            "ld_decay_n_loci_sample": 1,
        },
        fine_mapping={"min_gwas_n": 0},
    )
    return root, cfg, loci.iloc[0]


def _vendored(cfg, site, locus):
    """What ``site_aggregates`` emits, keyed for comparison."""
    reference = reference_variants(cfg)
    index = SiteIndex.load(cfg, site, reference)
    geno, pheno = site_aggregates(cfg, site, locus, POPS, INSTANCES, index=index)
    return (
        {b.pop: b for b in geno},
        {(b.instance, b.pop): b for b in pheno},
        index.n_flipped,
    )


def _shipped(root, site):
    """What the shipped trainer emits, decoded back into comparable pieces."""
    dataset, _ = get_dataset(root / site / "data", site, instances=INSTANCES)
    trainer = SiteFineMappingTrainer(
        train_dataset=dataset,
        train_configs={"pops": POPS, "trainer_output_dirname": str(root / "out" / site)},
        logger=logging.getLogger(f"parity.{site}"),
        client_id=site,
    )
    trainer.train()
    payload = trainer.get_parameters()

    import json

    manifest = json.loads(payload[KEY_MANIFEST].numpy().tobytes().decode("utf-8"))
    return payload, manifest


@pytest.mark.parametrize("site", sorted(SITES))
def test_shipped_site_stage_equals_the_vendored_one(package, site):
    root, cfg, locus = package
    geno, pheno, n_flipped = _vendored(cfg, site, locus)
    payload, manifest = _shipped(root, site)

    assert manifest["n_flipped"] == n_flipped
    assert {b["pop"] for b in manifest["geno_blocks"]} == set(geno)

    for block in manifest["geno_blocks"]:
        pop = block["pop"]
        want = geno[pop]
        key = f"{KEY_GRAM}{KEY_SEP}L0000{KEY_SEP}{pop}"

        # Zero tolerance, deliberately. See the module docstring.
        np.testing.assert_array_equal(payload[key].numpy(), want.G)
        np.testing.assert_array_equal(
            payload[f"{KEY_USUM}{KEY_SEP}L0000{KEY_SEP}{pop}"].numpy(), want.u
        )
        assert block["n"] == want.n
        assert block["n_incomplete_variants"] == want.n_incomplete_variants
        assert block["variants"]["snp_id"] == list(want.variants["snp_id"])
        assert block["variants"]["a1"] == list(want.variants["a1"])
        assert block["variants"]["a2"] == list(want.variants["a2"])

    for block in manifest["pheno_blocks"]:
        want = pheno[(block["instance"], block["pop"])]
        suffix = f"{KEY_SEP}{block['instance']}{KEY_SEP}{block['pop']}"
        np.testing.assert_array_equal(payload[f"{KEY_XTY}{suffix}"].numpy(), want.c)
        q, w, n = payload[f"{KEY_SCALARS}{suffix}"].numpy()
        assert q == want.q
        assert w == want.w
        assert int(n) == want.n


def test_the_flipped_site_is_actually_flipped(package):
    """Guard the fixture, not the code.

    If the harmonization branch stopped being exercised -- a changed fixture, a changed
    allele convention -- the parity test above would still pass while checking strictly
    less. A fixture that quietly stops testing the hard case is worse than no fixture.
    """
    root, cfg, locus = package
    flipped = {site: _vendored(cfg, site, locus)[2] for site in SITES}
    assert flipped[FLIPPED_SITE] == len(FLIPPED_VARIANTS)
    assert all(n == 0 for site, n in flipped.items() if site != FLIPPED_SITE)


def test_the_incomplete_variant_is_actually_dropped(package):
    """Same, for the complete-case path."""
    root, cfg, locus = package
    geno, _pheno, _flipped = _vendored(cfg, "anl", locus)
    dropped = sum(b.n_incomplete_variants for b in geno.values())
    assert dropped == 1, (
        "the fixture no longer contains a partially observed variant, so neither "
        "implementation's complete-case handling is being compared"
    )


def test_float32_gram_transport_is_lossless(package):
    """The uplink halving must not cost a bit.

    ``uplink_gram_dtype: float32`` exists because the Gram dominates the payload and its
    entries are integers bounded by 4n, which float32 holds exactly below 2^24. That is a
    proof, and this is the check that the proof was implemented rather than assumed.
    """
    root, cfg, locus = package
    dataset, _ = get_dataset(root / "anl" / "data", "anl", instances=INSTANCES)
    trainer = SiteFineMappingTrainer(
        train_dataset=dataset,
        train_configs={
            "pops": POPS,
            "uplink_gram_dtype": "float32",
            "trainer_output_dirname": str(root / "out" / "anl-f32"),
        },
        logger=logging.getLogger("parity.f32"),
        client_id="anl",
    )
    trainer.train()
    payload = trainer.get_parameters()

    geno, _pheno, _flipped = _vendored(cfg, "anl", locus)
    for pop, want in geno.items():
        got = payload[f"{KEY_GRAM}{KEY_SEP}L0000{KEY_SEP}{pop}"]
        assert got.dtype.itemsize == 4, "the Gram was not actually downcast"
        np.testing.assert_array_equal(got.numpy().astype(np.float64), want.G)


def test_the_shipped_loader_rejects_an_instance_no_site_holds(package):
    """Every site must hold the same instance set, or pooling shrinks a cohort silently."""
    root, _cfg, _locus = package
    with pytest.raises(FileNotFoundError, match="no phenotype for"):
        get_dataset(root / "anl" / "data", "anl", instances=["L0000_nosuch_rep0"])


def test_the_shipped_loader_rejects_an_unlabelled_individual(package, tmp_path):
    """An individual with no ancestry label would be dropped from every aggregate."""
    import shutil

    root, _cfg, _locus = package
    broken = tmp_path / "broken"
    shutil.copytree(root / "anl" / "data", broken)
    manifest = pd.read_csv(broken / "site_manifest.tsv", sep="\t")
    manifest.iloc[:-1].to_csv(broken / "site_manifest.tsv", sep="\t", index=False)

    with pytest.raises(ValueError, match="no ancestry label"):
        get_dataset(broken, "anl")
