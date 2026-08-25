"""Fine-mapping site data loader. SHIPPED TO WORKERS -- self-contained by design.

=============================================================================
THIS FILE'S SOURCE IS SENT OVER THE WIRE AND EXECUTED ON A PARTNER'S CLUSTER.
=============================================================================

It may import only the standard library and packages the partner installed with
``appfl-bio-suite[finemapping]``. It must NOT import from ``appfl_bio_suite``, and it must
not import a sibling module by bare name -- the shipped source runs from a temporary
working directory where its siblings do not exist.

Enforced by tests/test_shipped_modules.py. The rule exists because the alternative was
tried on the GWAS experiment: sibling imports meant every partner had to add a
``PYTHONPATH`` entry to their endpoint's ``worker_init``, and a missing or wrong one was
the single largest source of partner-side breakage on that project.

WHAT THIS LOADER DOES NOT DO
----------------------------
It does not read the genotypes. A site's chromosome fileset is tens of gigabytes and the
trainer reads only the few thousand variants inside each locus window, seeking straight
to them in the ``.bed``. Loading the matrix here would exhaust a worker before any
statistics were computed.

What it does instead is validate -- eagerly, at construction, before a scheduler
allocation has been spent -- that every file the trainer will need is present and mutually
consistent: the manifest covers the fam, the phenotypes cover the manifest, and the loci
are inside the chromosome. Each of those has failed on a real bundle, and each is far
cheaper to discover here than three minutes into a run on someone else's cluster.
"""

import csv
from pathlib import Path

# The stem of the PLINK1 triple inside a site's bundle. Fixed rather than derived from
# the client id: a bundle should be inspectable without knowing which site it was cut
# for, and a stem that encodes the site name is a stem that goes stale when a site is
# renamed in federation.yaml.
PLINK_STEM = "site_genotypes"

# Files that must sit DIRECTLY in data_dir. A bundle unpacked one level too deep is the
# most common partner-side mistake, and naming the files in the error is what makes that
# diagnosable without a round trip.
REQUIRED_FILES = (
    f"{PLINK_STEM}.bed",
    f"{PLINK_STEM}.bim",
    f"{PLINK_STEM}.fam",
    "site_manifest.tsv",
    "reference_variants.tsv",
    "selected_loci.tsv",
)

# Phenotypes live one level down because there is one file per instance and a real run
# has hundreds -- 100 loci x 15 architectures x 10 replicates. Flattening them into
# data_dir would make the directory unlistable.
PHENOTYPE_DIRNAME = "phenotypes"


def _read_tsv(path, required_columns):
    """Read a small TSV into a list of dicts, checking its header names.

    Deliberately ``csv`` rather than pandas. These files are metadata -- a few thousand
    rows at most -- and the loader's job is to fail with a clear message before the
    trainer imports anything heavy. A header typo caught here names the column; the same
    typo caught inside pandas surfaces as a KeyError three functions deeper.
    """
    path = Path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in required_columns if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(
                f"{path} is missing column(s) {missing}.\n"
                f"  found: {reader.fieldnames}\n"
                f"  expected at least: {list(required_columns)}"
            )
        return list(reader)


class SiteFineMappingDataset:
    """One site's cohort for fine-mapping: paths, ancestry composition, and loci.

    Validation is strict and happens at construction, for the reason above. Everything it
    reads is small; the genotype matrix is never touched.
    """

    def __init__(self, data_dir, site_id, instances=None):
        self.data_dir = Path(data_dir).resolve()
        self.site_id = str(site_id)

        if not self.data_dir.is_dir():
            raise FileNotFoundError(
                f"{self.site_id}: data_dir does not exist: {self.data_dir}\n"
                "\n"
                "This path is resolved on YOUR cluster, by the worker. It must be an "
                "absolute path readable by the service account the endpoint maps tasks "
                "to -- not a path on the coordinator's machine."
            )

        missing = [name for name in REQUIRED_FILES if not (self.data_dir / name).is_file()]
        if missing:
            present = sorted(p.name for p in self.data_dir.iterdir())
            raise FileNotFoundError(
                f"{self.site_id}: missing required input files in {self.data_dir}\n"
                "  missing: " + ", ".join(missing) + "\n"
                "  present: " + (", ".join(present) if present else "(directory is empty)")
                + "\n\n"
                "Unpack the bundle you were sent so that these files sit directly in "
                "data_dir, not in a nested subdirectory."
            )

        self.plink_prefix = self.data_dir / PLINK_STEM
        self.plink_bed = self.data_dir / f"{PLINK_STEM}.bed"
        self.plink_bim = self.data_dir / f"{PLINK_STEM}.bim"
        self.plink_fam = self.data_dir / f"{PLINK_STEM}.fam"
        self.manifest_path = self.data_dir / "site_manifest.tsv"
        self.reference_path = self.data_dir / "reference_variants.tsv"
        self.loci_path = self.data_dir / "selected_loci.tsv"
        self.phenotype_dir = self.data_dir / PHENOTYPE_DIRNAME

        # -- the fam: one line per individual, in .bed column order ------------
        with self.plink_fam.open("r", encoding="utf-8") as handle:
            fam_ids = [
                line.split()[1] for line in handle if line.strip()
            ]
        self.sample_size = len(fam_ids)
        if self.sample_size == 0:
            raise ValueError(
                f"{self.site_id}: {self.plink_fam} lists no samples. The bundle is "
                "empty or truncated -- check the transfer."
            )

        # -- the manifest: which of those individuals are which ancestry ------
        manifest = _read_tsv(self.manifest_path, ("FID", "IID", "superpopulation"))
        fam_set = set(fam_ids)
        manifest_ids = {row["IID"] for row in manifest}
        unknown = manifest_ids - fam_set
        if unknown:
            raise ValueError(
                f"{self.site_id}: site_manifest.tsv names {len(unknown)} individual(s) "
                f"that are not in {self.plink_fam.name}, e.g. {sorted(unknown)[:3]}.\n"
                "The manifest and the fileset must describe the same cohort. This "
                "usually means two bundles' files were mixed."
            )
        unlabelled = fam_set - manifest_ids
        if unlabelled:
            raise ValueError(
                f"{self.site_id}: {len(unlabelled)} individual(s) in "
                f"{self.plink_fam.name} have no ancestry label in site_manifest.tsv, "
                f"e.g. {sorted(unlabelled)[:3]}.\n"
                "Every individual must be assigned to a superpopulation: the site stage "
                "partitions its cohort by ancestry and an unlabelled individual would be "
                "silently dropped from every aggregate."
            )

        # Ancestries this site actually holds, and how many of each. Reported to the
        # coordinator so it can check its column plan against reality rather than
        # assuming the composition it recorded when the bundle was cut.
        composition = {}
        for row in manifest:
            pop = row["superpopulation"]
            composition[pop] = composition.get(pop, 0) + 1
        self.composition = dict(sorted(composition.items()))

        # -- the loci to fine-map ---------------------------------------------
        self.loci = _read_tsv(
            self.loci_path, ("locus_id", "chrom", "start_bp", "end_bp")
        )
        if not self.loci:
            raise ValueError(f"{self.site_id}: selected_loci.tsv lists no loci")

        # -- the instances (locus x architecture x replicate) ------------------
        # A site holds one .pheno per instance. The coordinator decides which subset a
        # run covers; absent an explicit list, everything present is offered.
        if not self.phenotype_dir.is_dir():
            raise FileNotFoundError(
                f"{self.site_id}: no {PHENOTYPE_DIRNAME}/ directory in {self.data_dir}.\n"
                "It holds one <instance>.pheno per (locus, architecture, replicate)."
            )
        available = sorted(p.stem for p in self.phenotype_dir.glob("*.pheno"))
        if not available:
            raise FileNotFoundError(
                f"{self.site_id}: {self.phenotype_dir} contains no .pheno files."
            )
        self.available_instances = available

        if instances is None:
            self.instances = available
        else:
            requested = [str(i) for i in instances]
            absent = [i for i in requested if i not in set(available)]
            if absent:
                raise FileNotFoundError(
                    f"{self.site_id}: the coordinator asked for {len(absent)} instance(s) "
                    f"this site has no phenotype for, e.g. {absent[:3]}.\n"
                    "Every site must hold the same instance set -- the coordinator pools "
                    "one instance's aggregates across all of them, and a site missing an "
                    "instance would silently shrink that instance's cohort."
                )
            self.instances = requested

    def __len__(self):
        return self.sample_size

    def __repr__(self):
        return (
            f"SiteFineMappingDataset({self.site_id!r}, n={self.sample_size}, "
            f"pops={list(self.composition)}, loci={len(self.loci)}, "
            f"instances={len(self.instances)}, dir={self.data_dir})"
        )


def get_dataset(data_dir, site_id, instances=None, **kwargs):
    """APPFL dataset entry point.

    Returns ``(train_dataset, val_dataset)``. There is no validation split, for the same
    reason the GWAS experiment has none: this is single-round aggregate federated
    learning, not iterative training, so there is no global model to score against a
    held-out set each round. Fine-mapping accuracy is measured against the simulated
    ground truth at the coordinator, which is the only place that knows it.
    """
    return SiteFineMappingDataset(data_dir=data_dir, site_id=site_id, instances=instances), None
