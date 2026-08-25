"""Parameters for the GWAS simulation, and the scenario file that carries them.

Everything here was a module-level constant in the pre-migration scripts. Moving them
into a versioned configuration object is the point of the migration: the site split in
particular was five integers inline in the bundler, with an assertion that they summed to
exactly 100,000, which made "five sites with this skew" the only split the code could
express.

THE DEFAULTS ARE NOT ARBITRARY. They reproduce the simulation design of the paper this
experiment supports, and changing one changes what the experiment measures. They are
defaults so that the published configuration is what you get by not making a decision.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

__all__ = ["PhenotypeParams", "CohortParams", "GwasScenario"]


@dataclass(frozen=True)
class PhenotypeParams:
    """Genetic architecture and covariate effects for the two simulated traits.

    Two traits, chosen because they exercise different analyses: BMI is continuous and
    tested by covariate-adjusted linear regression; T2D is binary and tested by a
    logistic score test.
    """

    # Heritability -- the fraction of trait variance attributable to the polygenic score.
    # Sets the environmental noise variance as (1 - h2) / h2.
    h2_t2d: float = 0.20
    h2_bmi: float = 0.25

    # Target T2D prevalence. Case status is assigned by thresholding the liability at
    # this quantile, which is the standard liability-threshold model.
    t2d_prevalence: float = 0.119

    # Covariate effects on each trait's liability / value.
    age_coef_t2d: float = 0.05
    sex_coef_t2d: float = 0.35
    age_coef_bmi: float = 0.02
    sex_coef_bmi: float = 2.0
    bmi_intercept: float = 25.0

    # Multiplicative scale on the genetic effect. At 1.0 the polygenic score explains its
    # full nominal share; 0.5 halves it, which is what the published design uses -- it
    # keeps the simulated effects detectable without making the GWAS trivially easy.
    effect_scale_t2d: float = 0.50
    effect_scale_bmi: float = 0.50

    # Seeds. Three distinct streams, deliberately:
    #   - covariates are drawn once and shared by both replicates
    #   - the GWAS replicate and the PGS-evaluation replicate get INDEPENDENT noise, so
    #     polygenic scores are evaluated on a fresh draw rather than on the data that
    #     produced them. Reusing one draw would inflate every reported R2 and AUROC.
    seed_covariates: int = 42
    seed_gwas_replicate: int = 42
    seed_pgs_replicate: int = 43

    # Fraction of variants retained at simulation time. 1.0 is a full run; smaller values
    # are a smoke fraction. Compounds with the analysis-time variant scaling.
    data_sim_scaling: float = 1.0
    data_sim_scaling_seed: int = 7

    # Chunk size for the polygenic-score dot product.
    #
    # NUMERICALLY LOAD-BEARING: floating-point addition is not associative, so changing
    # this changes the accumulation order and therefore the low-order bits of every
    # score. Byte-identical output requires this value to stay at 2000.
    pgs_chunk_size: int = 2000

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CohortParams:
    """How the input genotype cohort is obtained.

    Two providers, because the original cohort and a reproducible one are different
    things and the repository has to be honest about which is in use:

    ``plink_prefix``
        Use an existing PLINK1 triple. This is the path for the original 100,000-sample
        cohort -- and that cohort is NOT regenerable from this repository. See
        docs/experiments/gwas/DATA.md.

    ``synthetic_ld``
        Generate a cohort from the documented model in cohort.py. Deterministic from a
        seed, needs no download, and scales from a CI-sized fixture to full size. This is
        what makes the pipeline runnable end to end by someone who does not have the
        original data -- which is most people.
    """

    provider: str = "synthetic_ld"

    # provider == "plink_prefix"
    plink_prefix: str | None = None

    # provider == "synthetic_ld"
    n_samples: int = 10_000
    n_variants: int = 20_000
    seed: int = 20260101
    n_chromosomes: int = 22
    # Variants per LD block. Real genotypes are correlated in blocks; independent
    # variants would make every association test cleanly independent and the simulation
    # unrealistically easy.
    ld_block_size: int = 25
    ld_correlation: float = 0.65
    # Beta distribution shaping the minor-allele frequency spectrum, then clipped.
    maf_beta_a: float = 0.8
    maf_beta_b: float = 2.4
    maf_min: float = 0.01
    maf_max: float = 0.5
    missing_rate: float = 0.001

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class GwasScenario:
    """A complete, named GWAS simulation scenario."""

    name: str
    description: str = ""
    cohort: CohortParams = field(default_factory=CohortParams)
    phenotypes: PhenotypeParams = field(default_factory=PhenotypeParams)
    # Site sizes, in order. Site N gets sites[N-1] samples. Deliberately a list rather
    # than a count plus a split rule: the published design uses a specific unequal split
    # that mimics real cross-site heterogeneity, and that is data, not an algorithm.
    sites: list[int] = field(default_factory=lambda: [18032, 9237, 25028, 40752, 6951])
    # Seed for the sample shuffle that assigns individuals to sites.
    split_seed: int = 42

    @property
    def total_samples(self) -> int:
        return sum(self.sites)

    @property
    def site_count(self) -> int:
        return len(self.sites)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> GwasScenario:
        return cls(
            name=data["name"],
            description=data.get("description", ""),
            cohort=CohortParams(**(data.get("cohort") or {})),
            phenotypes=PhenotypeParams(**(data.get("phenotypes") or {})),
            sites=list(data.get("sites") or []),
            split_seed=int(data.get("split_seed", 42)),
        )

    @classmethod
    def load(cls, path) -> GwasScenario:
        from pathlib import Path

        import yaml

        raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            raise ValueError(f"{path} must contain a YAML mapping")
        return cls.from_dict(raw)

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "cohort": self.cohort.to_dict(),
            "phenotypes": self.phenotypes.to_dict(),
            "sites": list(self.sites),
            "split_seed": self.split_seed,
        }

    def validate(self) -> None:
        if not self.sites:
            raise ValueError(
                f"scenario '{self.name}' declares no sites. Add at least one size to "
                "`sites:` -- the list is the per-site sample counts, in order."
            )
        if any(size <= 0 for size in self.sites):
            raise ValueError(f"scenario '{self.name}' has a non-positive site size: {self.sites}")

        if self.cohort.provider == "plink_prefix":
            if not self.cohort.plink_prefix:
                raise ValueError(
                    f"scenario '{self.name}' uses provider 'plink_prefix' but sets no "
                    "`cohort.plink_prefix`. Point it at your PLINK1 triple, without the "
                    "file extension."
                )
        elif self.cohort.provider == "synthetic_ld":
            if self.total_samples > self.cohort.n_samples:
                raise ValueError(
                    f"scenario '{self.name}' splits {self.total_samples} samples across "
                    f"sites, but the cohort only generates {self.cohort.n_samples}. "
                    "Either raise cohort.n_samples or lower the site sizes."
                )
        else:
            raise ValueError(
                f"scenario '{self.name}' has unknown cohort provider "
                f"'{self.cohort.provider}'. Known: synthetic_ld, plink_prefix."
            )
