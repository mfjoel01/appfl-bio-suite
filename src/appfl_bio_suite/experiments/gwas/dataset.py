"""GWAS site data loader. SHIPPED TO WORKERS -- self-contained by design.

=============================================================================
THIS FILE'S SOURCE IS SENT OVER THE WIRE AND EXECUTED ON A PARTNER'S CLUSTER.
=============================================================================

It may import only the standard library and packages the partner installed with
``appfl-bio-suite[gwas]``. It must NOT import from ``appfl_bio_suite``.

Enforced by tests/test_shipped_modules.py. The rule exists because the alternative was
tried: the original GWAS modules imported two sibling helpers by bare name, which meant
every partner had to add a ``PYTHONPATH`` entry to their endpoint's ``worker_init``, and
a missing or wrong one was the single largest source of partner-side breakage on that
project. Self-contained modules removed that line from the setup guide entirely.

This loader does not read the genotypes. It validates that the site's files are present
and reports the sample count, then hands the paths to the trainer, which streams the
genotype matrix in chunks. Loading a multi-gigabyte matrix eagerly here would blow up a
worker before training started.
"""

from pathlib import Path

# The six files the trainer actually requires. The bundler writes three more
# (pgs_scores.csv, score_t2d.txt, score_bmi.txt) which are simulation and
# central-baseline artifacts -- useful to a coordinator, inert at a partner site.
#
# tests/test_gwas_loader_contract.py asserts this list agrees with what the bundler
# produces. That agreement used to be maintained by hand across two directories, and a
# mismatch surfaced on a partner's cluster rather than in CI.
REQUIRED_FILES = (
    "EUR.synthetic.100k.ld.maf.bed",
    "EUR.synthetic.100k.ld.maf.bim",
    "EUR.synthetic.100k.ld.maf.fam",
    "phenotypes_gwas.csv",
    "phenotypes_pgs_eval.csv",
    "covariates.csv",
)

PLINK_STEM = "EUR.synthetic.100k.ld.maf"


class SiteGWASDataset:
    """One site's cohort: paths plus a sample count, validated eagerly.

    Validation is deliberately strict and happens at construction. A missing file
    discovered mid-analysis costs a scheduler allocation and produces a confusing
    traceback; discovered here it produces a message naming the file.
    """

    def __init__(self, data_dir, site_id):
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

        self.plink_prefix = self.data_dir / PLINK_STEM
        self.plink_bed = Path(str(self.plink_prefix) + ".bed")
        self.plink_bim = Path(str(self.plink_prefix) + ".bim")
        self.plink_fam = Path(str(self.plink_prefix) + ".fam")
        self.pheno_gwas = self.data_dir / "phenotypes_gwas.csv"
        self.pheno_eval = self.data_dir / "phenotypes_pgs_eval.csv"
        self.covariates = self.data_dir / "covariates.csv"

        missing = [name for name in REQUIRED_FILES if not (self.data_dir / name).is_file()]
        if missing:
            present = sorted(p.name for p in self.data_dir.iterdir() if p.is_file())
            raise FileNotFoundError(
                f"{self.site_id}: missing required input files in {self.data_dir}\n"
                "  missing: " + ", ".join(missing) + "\n"
                "  present: " + (", ".join(present) if present else "(directory is empty)") + "\n\n"
                "Unpack the bundle you were sent so that these files sit directly in "
                "data_dir, not in a nested subdirectory."
            )

        # Counted from the .fam rather than by loading genotypes: one line per sample,
        # and it avoids touching a file that may be several gigabytes.
        with self.plink_fam.open("r", encoding="utf-8") as handle:
            self.sample_size = sum(1 for line in handle if line.strip())

        if self.sample_size == 0:
            raise ValueError(
                f"{self.site_id}: {self.plink_fam} lists no samples. The bundle is "
                "empty or truncated -- check the transfer."
            )

    def __len__(self):
        return self.sample_size

    def __repr__(self):
        return f"SiteGWASDataset({self.site_id!r}, n={self.sample_size}, dir={self.data_dir})"


def get_dataset(data_dir, site_id, **kwargs):
    """APPFL dataset entry point.

    Returns ``(train_dataset, val_dataset)``. There is no validation split: this is
    single-round summary-statistic federated learning, not iterative training, so there
    is no held-out set to score a global model against each round. Polygenic scores are
    evaluated against an independent phenotype replicate, which the trainer handles.
    """
    return SiteGWASDataset(data_dir=data_dir, site_id=site_id), None
