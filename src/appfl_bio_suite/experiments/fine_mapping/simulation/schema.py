"""What a fine-mapping simulation scenario declares, and how it becomes a pipeline config.

TWO CONFIG OBJECTS, DELIBERATELY NOT MERGED
-------------------------------------------
The vendored pipeline has its own configuration model -- ``fedfm/utils.py``'s
:class:`SimulationConfig`, a pydantic object loaded from the standalone repository's
``config/simulation_config.yaml``. It is thorough, validated, and already the single
source of truth for every parameter the science depends on.

Rewriting it into the suite's idiom would mean re-deriving forty validated fields and
would break the guarantee the vendoring exists to provide: that the code producing the
numbers is byte-identical to upstream. So a scenario **carries that config verbatim**,
under a ``pipeline:`` key, and this module supplies only the two things the suite owns:

* ``paths:`` -- injected, rooted at whatever ``simulate --out`` was given, so a scenario
  never hardcodes a directory and two runs of the same scenario cannot collide.
* ``cohort:`` -- where the input genotypes come from, which upstream assumed was a
  manually staged HAPNEST download and the suite has to be able to answer for itself.

A scenario is therefore readable side by side with the upstream YAML, and diffing them
answers "did the migration change any parameter?" mechanically.

THE THREE-SITE CONSTRAINT IS REAL, AND CHECKED
----------------------------------------------
``fedfm/sampling.py`` enumerates ``SITE_ORDER = ("anl", "covenant", "mbzuai")`` and
assigns individuals in that fixed order; a site named anything else is silently dropped
from the cohort assignment rather than rejected. That is a genuine limitation of the
vendored pipeline, not of the suite, and it is not worked around here -- working around it
would mean editing code whose value is being unedited.

What is done instead is to make the silent drop loud: :meth:`FineMappingScenario.validate`
refuses a scenario naming any other site, with a message saying where the list lives. A
fourth site is a change to make upstream and re-vendor.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

__all__ = ["CohortParams", "FineMappingScenario", "KNOWN_SITE_IDS"]

# Mirrors fedfm/sampling.py::SITE_ORDER. Duplicated rather than imported so that this
# module stays importable without pulling in numpy/pandas/pydantic, and kept honest by
# tests/test_fine_mapping_configs.py, which asserts the two agree.
KNOWN_SITE_IDS = ("anl", "covenant", "mbzuai")

# Where a coordinator points at their staged HAPNEST copy, so a committed scenario
# never has to carry one machine's path.
_HAPNEST_ENV = "FEDFM_HAPNEST_DIR"

# Subdirectories of the run's --out directory. The pipeline writes to all of them, and
# they are what `--verify` re-checksums.
_PATH_LAYOUT = {
    "hapnest_dir": "raw/hapnest",
    "processed_dir": "processed",
    "loci_dir": "loci",
    "ground_truth_dir": "ground_truth",
    "reports_dir": "reports",
    "logs_dir": "logs",
}


@dataclass(frozen=True)
class CohortParams:
    """How the input genotype pool is obtained.

    Two providers, because a staged 135 GB download and a runnable fixture are different
    things and the repository has to be honest about which produced a given result.

    ``hapnest``
        Use a staged copy of the HAPNEST pre-generated synthetic cohort (EBI BioStudies
        ``S-BSST936``). This is what the published run used. It is a 135 GB download for
        chr1 alone; ``scripts/fine-mapping/download_hapnest.sh`` fetches it.

    ``synthetic_hapnest``  (default)
        Generate a small HAPNEST-*shaped* pool from the documented model in cohort.py:
        the same file layout, the same superpopulation labelling, a fraction of the size.
        Deterministic from a seed and needs no download, which is what makes the whole
        chain runnable from this repository alone.

    **The substitute is not HAPNEST.** It is a different pool from a different generative
    process, so it does not reproduce the published figures and must not be described as
    doing so. It reproduces the *pipeline*, not the *numbers*. The run manifest records
    which provider was used, so no output is ever ambiguous about its own provenance.
    """

    provider: str = "synthetic_hapnest"

    # provider == "hapnest": directory holding chr<N>.{bed,bim,fam} and
    # population_manifest.tsv, as download_hapnest.sh lays them out.
    #
    # Shipped as null, and resolvable from $FEDFM_HAPNEST_DIR instead. A staged pool is a
    # fact about a machine, not about a scenario: hardcoding one coordinator's path into a
    # committed scenario means every other coordinator edits a tracked file to run it, and
    # then carries that edit forever.
    hapnest_dir: str | None = None

    # provider == "synthetic_hapnest"
    n_individuals: int = 3_000
    n_variants: int = 6_000
    seed: int = 20260101
    # Relative sizes of the six superpopulations in the generated pool. Absolute counts
    # come from the per-site compositions; this only has to be large enough to serve them.
    superpopulation_weights: dict[str, float] = field(
        default_factory=lambda: {
            "EUR": 0.30, "AFR": 0.30, "AMR": 0.10, "EAS": 0.08, "CSA": 0.14, "MID": 0.08
        }
    )
    # Variants per LD block, and the within-block correlation. Fine-mapping is entirely
    # about correlation structure, so a pool of independent variants would make credible
    # sets trivially small and the whole experiment meaningless.
    ld_block_size: int = 25
    ld_correlation: float = 0.70
    # How much the block structure and allele frequencies differ BETWEEN superpopulations.
    # This is the knob the locus-selection stage stratifies on: at 0 every ancestry has
    # identical LD and the low/medium/high divergence strata are indistinguishable.
    ld_divergence: float = 0.35
    maf_beta_a: float = 0.8
    maf_beta_b: float = 2.4
    maf_min: float = 0.01
    maf_max: float = 0.5
    missing_rate: float = 0.0005
    # Base-pair spacing between consecutive variants. With the default locus window of
    # 1 Mb, this sets how many variants land in a window.
    bp_spacing: int = 2_000
    bp_start: int = 100_000

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def resolved_hapnest_dir(self) -> Path | None:
        """The staged pool's directory, from the scenario or from the environment."""
        raw = self.hapnest_dir or os.environ.get(_HAPNEST_ENV)
        return Path(raw).expanduser().resolve() if raw else None


@dataclass(frozen=True)
class FineMappingScenario:
    """A complete, named fine-mapping simulation scenario."""

    name: str
    description: str = ""
    cohort: CohortParams = field(default_factory=CohortParams)
    # The upstream SimulationConfig, verbatim, minus `paths`. Kept as a plain dict rather
    # than parsed here: fedfm/utils.py validates it, and validating it twice in two
    # places is how the two definitions drift apart.
    pipeline: dict[str, Any] = field(default_factory=dict)

    @property
    def site_ids(self) -> list[str]:
        return list((self.pipeline.get("sites") or {}).keys())

    @property
    def site_count(self) -> int:
        return len(self.site_ids)

    @property
    def total_samples(self) -> int:
        return sum(int(s["n"]) for s in (self.pipeline.get("sites") or {}).values())

    @property
    def chromosome(self) -> int:
        return int(self.pipeline.get("chromosome", 1))

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> FineMappingScenario:
        return cls(
            name=data["name"],
            description=data.get("description", ""),
            cohort=CohortParams(**(data.get("cohort") or {})),
            pipeline=dict(data.get("pipeline") or {}),
        )

    @classmethod
    def load(cls, path) -> FineMappingScenario:
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
            "pipeline": self.pipeline,
        }

    # -- validation --------------------------------------------------------

    def validate(self) -> None:
        """Refuse a scenario that would fail late, or worse, quietly succeed wrongly."""
        if "sites" not in self.pipeline or not self.pipeline["sites"]:
            raise ValueError(
                f"scenario '{self.name}' declares no sites under `pipeline.sites`. "
                "Each entry is {n, composition, dominant}, as in the upstream "
                "simulation_config.yaml."
            )

        unknown = [s for s in self.site_ids if s not in KNOWN_SITE_IDS]
        if unknown:
            raise ValueError(
                f"scenario '{self.name}' names site(s) {unknown}, which the vendored "
                f"cohort sampler does not know about. It assigns individuals in a fixed "
                f"order over exactly {list(KNOWN_SITE_IDS)} and would silently drop the "
                "others, producing a cohort smaller than the scenario asked for.\n"
                "\n"
                "The list is SITE_ORDER in "
                "src/appfl_bio_suite/experiments/fine_mapping/fedfm/sampling.py. Adding a "
                "site is a change to make in the upstream repository and re-vendor, not "
                "here -- see fedfm/__init__.py."
            )

        for key in ("master_seed", "chromosome", "superpopulations", "locus_selection",
                    "architecture", "phenotype", "qc"):
            if key not in self.pipeline:
                raise ValueError(
                    f"scenario '{self.name}' is missing `pipeline.{key}`. The pipeline "
                    "block is the upstream SimulationConfig verbatim; compare it against "
                    "docs/experiments/fine-mapping/DATA.md."
                )

        if self.cohort.provider == "hapnest":
            # Deliberately NOT checking that the pool is actually staged. validate()
            # answers "is this scenario well-formed", which is a property of the file;
            # whether a 135 GB download is present is a property of the machine, and
            # conflating them means `--list-scenarios` cannot even list the published
            # scenario on a machine that has not downloaded it yet. materialize_pool()
            # raises the provisioning error, where the fix is actionable.
            pass
        elif self.cohort.provider == "synthetic_hapnest":
            if self.total_samples > self.cohort.n_individuals:
                raise ValueError(
                    f"scenario '{self.name}' enrols {self.total_samples:,} individuals "
                    f"across its sites, but the generated pool holds only "
                    f"{self.cohort.n_individuals:,}. Sites are disjoint -- no individual "
                    "is reused -- so the pool must be at least as large. Raise "
                    "cohort.n_individuals or lower the site sizes."
                )
            # Each superpopulation must be able to serve total demand for it, or the
            # sampler raises deep inside the assignment with a less helpful message.
            demand: dict[str, int] = {}
            for site in self.pipeline["sites"].values():
                for pop, count in (site.get("composition") or {}).items():
                    demand[pop] = demand.get(pop, 0) + int(count)
            for pop, needed in sorted(demand.items()):
                weight = self.cohort.superpopulation_weights.get(pop, 0.0)
                available = int(self.cohort.n_individuals * weight)
                if available < needed:
                    raise ValueError(
                        f"scenario '{self.name}' needs {needed:,} {pop} individuals "
                        f"across its sites, but the generated pool would hold about "
                        f"{available:,} ({pop} weight {weight:g} of "
                        f"{self.cohort.n_individuals:,}). Raise cohort.n_individuals, "
                        f"raise the {pop} weight, or shrink the site compositions."
                    )
        else:
            raise ValueError(
                f"scenario '{self.name}' has unknown cohort provider "
                f"'{self.cohort.provider}'. Known: synthetic_hapnest, hapnest."
            )

    # -- becoming a pipeline config ----------------------------------------

    def paths_for(self, out_dir: Path, pool_dir: Path | None = None) -> dict[str, str]:
        """The `paths` block the upstream config wants, rooted at this run's output.

        Absolute, because the pipeline's own ``resolved_path`` would otherwise resolve
        them against a repo root that means nothing here, and because the per-site
        directories end up named in generated configs.

        ``pool_dir`` overrides where the genotype pool lives. A staged HAPNEST copy sits
        outside the run directory and is shared by every run on the machine -- symlinking
        135 GB into each one, or copying it, would both be absurd.
        """
        out_dir = Path(out_dir).resolve()
        paths = {key: str(out_dir / rel) for key, rel in _PATH_LAYOUT.items()}
        pool = Path(pool_dir).resolve() if pool_dir else out_dir / _PATH_LAYOUT["hapnest_dir"]
        paths["hapnest_dir"] = str(pool)
        paths["population_manifest"] = str(pool / "population_manifest.tsv")
        return paths

    def pool_dir(self, out_dir: Path) -> Path:
        """Where this scenario's genotype pool lives, staged or generated."""
        if self.cohort.provider == "hapnest":
            staged = self.cohort.resolved_hapnest_dir()
            if staged is None:
                raise ValueError(
                    f"scenario '{self.name}' uses cohort provider 'hapnest' but no "
                    "staged pool is configured.\n"
                    "Set the directory holding "
                    f"chr{self.chromosome}.{{bed,bim,fam}} and population_manifest.tsv, "
                    "either as `cohort.hapnest_dir` in your own copy of the scenario or "
                    f"as ${_HAPNEST_ENV} in the environment.\n"
                    "scripts/fine-mapping/download_hapnest.sh stages it; see "
                    "docs/experiments/fine-mapping/DATA.md."
                )
            return staged
        return Path(out_dir).resolve() / _PATH_LAYOUT["hapnest_dir"]

    def to_pipeline_config(self, out_dir: Path):
        """Build the vendored :class:`SimulationConfig` for a run rooted at ``out_dir``."""
        from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import SimulationConfig

        return SimulationConfig(
            repo_root=Path(out_dir).resolve(),
            paths=self.paths_for(out_dir, self.pool_dir(out_dir)),
            **self.pipeline,
        )

    def seeds(self) -> dict[str, int]:
        """Every seed that affects the output, for the run manifest."""
        return {
            "master": int(self.pipeline.get("master_seed", 0)),
            "cohort": int(self.cohort.seed),
        }
