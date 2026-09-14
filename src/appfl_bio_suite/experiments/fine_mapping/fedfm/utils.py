"""Shared utilities for the FedFM simulation pipeline.

Provides:
- Typed YAML config loading with validation (pydantic).
- Deterministic seed derivation via numpy SeedSequence.
- A safe PLINK subprocess wrapper that fails loudly.
- Logging setup (console + file).
- Small I/O helpers for PLINK fam/manifest files.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import yaml
from pydantic import BaseModel, Field, PositiveFloat, PositiveInt, model_validator


# ---------------------------------------------------------------------------
# Config models
# ---------------------------------------------------------------------------


class PathsConfig(BaseModel):
    hapnest_dir: str
    population_manifest: str
    processed_dir: str
    loci_dir: str
    ground_truth_dir: str
    reports_dir: str
    logs_dir: str


class ToolsConfig(BaseModel):
    plink: str = "plink"
    plink2: str = "plink2"
    king: str = "king"


class SiteConfig(BaseModel):
    n: PositiveInt
    composition: dict[str, PositiveInt]
    dominant: str

    @model_validator(mode="after")
    def _check_composition(self) -> "SiteConfig":
        total = sum(self.composition.values())
        if total != self.n:
            raise ValueError(
                f"Site composition sums to {total} but n={self.n}"
            )
        if self.dominant not in self.composition:
            raise ValueError(
                f"Dominant superpop {self.dominant} not in composition"
            )
        return self


class StrataConfig(BaseModel):
    low: PositiveInt
    medium: PositiveInt
    high: PositiveInt


class LocusSelectionConfig(BaseModel):
    window_size_bp: PositiveInt
    step_size_bp: PositiveInt
    n_loci: PositiveInt
    strata: StrataConfig
    maf_filter: PositiveFloat
    ld_tag_snp_target: PositiveInt
    ld_r2_prune_threshold: PositiveFloat
    ld_prune_window_variants: PositiveInt
    ld_prune_step_variants: PositiveInt
    n_workers: PositiveInt = 1

    @model_validator(mode="after")
    def _check_strata(self) -> "LocusSelectionConfig":
        s = self.strata.low + self.strata.medium + self.strata.high
        if s != self.n_loci:
            raise ValueError(
                f"Strata sum {s} != n_loci {self.n_loci}"
            )
        return self


class AncestrySpecificCausalConfig(BaseModel):
    enabled: bool
    min_per_stratum: int
    common_maf_threshold: PositiveFloat
    rare_maf_threshold: PositiveFloat


class AncestryDivergentCausalConfig(BaseModel):
    """Causal *set* differs across superpopulations (distinct from
    ``AncestrySpecificCausalConfig``, which only diverges one variant's MAF).

    A divergent instance uses a union causal set of ``n_shared`` variants common
    to all superpops plus ``n_private_per_pop`` variants private to each. Only
    fires when ``ncsl > n_private_per_pop`` (so at least one shared variant
    remains); otherwise the instance falls back to a shared causal set.
    """

    enabled: bool = False
    n_private_per_pop: PositiveInt = 1
    min_per_stratum: int = 1


class ArchitectureConfig(BaseModel):
    ncsl: list[PositiveInt]
    h2: list[PositiveFloat]
    rg: list[float]
    factorial_mode: str
    replicates: PositiveInt
    ancestry_specific_causal: AncestrySpecificCausalConfig
    ancestry_divergent_causal: AncestryDivergentCausalConfig = (
        AncestryDivergentCausalConfig()
    )

    @model_validator(mode="after")
    def _check_mode(self) -> "ArchitectureConfig":
        if self.factorial_mode not in {"minimal", "extended", "full"}:
            raise ValueError(f"Unknown factorial_mode {self.factorial_mode}")
        for r in self.rg:
            if not -1.0 <= r <= 1.0:
                raise ValueError(f"rg {r} outside [-1, 1]")
        return self


class PhenotypeConfig(BaseModel):
    model: str
    h2_tolerance_relative: PositiveFloat


class QCConfig(BaseModel):
    kinship_threshold: PositiveFloat
    pca_n_components: PositiveInt
    pca_reference: str
    maf_tolerance_sd: PositiveFloat
    ld_decay_max_dist_kb: PositiveInt
    ld_decay_n_loci_sample: PositiveInt


class FineMappingConfig(BaseModel):
    """Settings for the centralized SuSiEx fine-mapping baseline.

    One SuSiEx population column per superpopulation, pooling that ancestry's
    individuals across every site (AFR at ANL + Covenant + MBZUAI = one AFR
    column). This is the centralized reference: it is what running SuSiEx on all
    the data held in one place would produce, and it matches the simulator's
    generative model, where effect sizes are indexed by superpopulation rather
    than by site.

    ``min_gwas_n`` drops any ancestry whose pooled cohort is smaller than this,
    keeping SuSiEx from being handed a column too small to inform a GWAS or a
    stable LD panel. At 1000 all six superpopulations are retained (the smallest
    pooled cohort is EAS at 2,500).

    ``maf`` is the frequency filter applied to the POOLED columns. It belongs to
    the scenario rather than to either stage's argv because it is the one setting
    the centralized baseline and the federated path have to agree on: the two
    result tables are compared column for column, and filtering at two different
    thresholds fine-maps two different variant sets, which shows up as a
    credible-set disagreement indistinguishable from a bug in the federated path.
    A small scenario needs a lower value than the published 0.005 -- at ci-tiny's
    cohort size that filter can leave a pooled column with too few variants for
    SuSiEx to do anything with.
    """

    min_gwas_n: int = 0
    maf: float = 0.005
    keep_ambiguous: bool = True
    n_signals: PositiveInt = 10
    max_iter: PositiveInt = 1000
    tol: PositiveFloat = 1e-6

    def inference_options(self) -> dict:
        return {k: getattr(self, k) for k in ("keep_ambiguous", "n_signals", "max_iter", "tol")}


    @model_validator(mode="after")
    def _check(self) -> "FineMappingConfig":
        if self.min_gwas_n < 0:
            raise ValueError("min_gwas_n must be >= 0")
        if not 0.0 <= self.maf < 0.5:
            raise ValueError("maf must be in [0, 0.5)")
        return self


class SimulationConfig(BaseModel):
    paths: PathsConfig
    master_seed: int
    chromosome: PositiveInt
    tools: ToolsConfig
    sites: dict[str, SiteConfig]
    superpopulations: list[str]
    locus_selection: LocusSelectionConfig
    architecture: ArchitectureConfig
    phenotype: PhenotypeConfig
    qc: QCConfig
    # Optional so configs predating the fine-mapping stage still validate.
    fine_mapping: FineMappingConfig = FineMappingConfig()

    # Resolved at load time so all callers see absolute paths.
    repo_root: Path = Field(default_factory=lambda: Path.cwd())

    @model_validator(mode="after")
    def _check_pops(self) -> "SimulationConfig":
        known = set(self.superpopulations)
        for site_name, site in self.sites.items():
            for pop in site.composition:
                if pop not in known:
                    raise ValueError(
                        f"Site {site_name} references unknown superpop {pop}"
                    )
        lower = -1.0 / (len(self.superpopulations) - 1) if len(self.superpopulations) > 1 else -1.0
        if any(r <= lower for r in self.architecture.rg):
            raise ValueError(f"Equicorrelation must be greater than {lower}")
        if any(not 0 < h < 1 for h in self.architecture.h2):
            raise ValueError("Heritability must lie strictly between zero and one")
        return self

    def resolved_path(self, key: str) -> Path:
        """Return an absolute path for a paths.* entry, resolved against repo_root."""
        raw = getattr(self.paths, key)
        p = Path(raw)
        if not p.is_absolute():
            p = self.repo_root / p
        return p

    def site_dir(self, site: str) -> Path:
        return self.resolved_path("processed_dir") / site


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------


def load_config(
    config_path: str | os.PathLike[str], repo_root: str | os.PathLike[str] | None = None
) -> SimulationConfig:
    """Load and validate a simulation YAML config.

    Parameters
    ----------
    config_path : path-like
        Path to YAML config.
    repo_root : path-like, optional
        Base directory for resolving relative paths. Defaults to the parent
        of the config file (which should normally be the repo root if the
        file lives at ``<repo>/config/simulation_config.yaml``).

    Returns
    -------
    SimulationConfig
    """
    config_path = Path(config_path).resolve()
    with config_path.open() as fh:
        raw = yaml.safe_load(fh)

    if repo_root is None:
        # Assume config lives at <repo>/config/<file>.yaml
        repo_root_p = config_path.parent.parent
    else:
        repo_root_p = Path(repo_root).resolve()

    cfg = SimulationConfig(repo_root=repo_root_p, **raw)
    return cfg


# ---------------------------------------------------------------------------
# Deterministic seeding
# ---------------------------------------------------------------------------


def _to_seed_entropy(component: Any) -> int:
    """Convert an arbitrary deterministic key component to a non-negative int.

    Strings are hashed via a stable algorithm (not Python's randomized hash).
    """
    if isinstance(component, (int, np.integer)):
        return int(component) & 0xFFFFFFFF
    if isinstance(component, (float, np.floating)):
        # Convert to a stable integer representation
        return int(np.float64(component).view(np.int64)) & 0xFFFFFFFF
    # Stable string hash: fold bytes
    s = str(component).encode("utf-8")
    h = 1469598103934665603  # FNV-1a 64-bit offset basis
    for b in s:
        h ^= b
        h = (h * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h & 0xFFFFFFFF


def derive_seed(master_seed: int, *components: Any) -> np.random.SeedSequence:
    """Derive a child SeedSequence from a master seed and deterministic components.

    Same (master_seed, components) → same SeedSequence on every invocation,
    regardless of process or Python version.
    """
    entropy = [int(master_seed) & 0xFFFFFFFF] + [_to_seed_entropy(c) for c in components]
    return np.random.SeedSequence(entropy)


def make_rng(master_seed: int, *components: Any) -> np.random.Generator:
    """Construct a deterministic numpy Generator from master seed + components."""
    return np.random.default_rng(derive_seed(master_seed, *components))


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------


_LOGGER_NAME = "fedfm"


def setup_logging(
    log_dir: str | os.PathLike[str] | None = None,
    name: str = _LOGGER_NAME,
    level: int = logging.INFO,
) -> logging.Logger:
    """Configure a structured logger writing to console (+ optional file)."""
    logger = logging.getLogger(name)
    if logger.handlers:
        # Already configured.
        return logger

    logger.setLevel(level)
    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    ch = logging.StreamHandler(sys.stderr)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    if log_dir is not None:
        log_dir_p = Path(log_dir)
        log_dir_p.mkdir(parents=True, exist_ok=True)
        ts = time.strftime("%Y%m%dT%H%M%S")
        fh = logging.FileHandler(log_dir_p / f"{name}_{ts}.log")
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    logger.propagate = False
    return logger


def get_logger(name: str = _LOGGER_NAME) -> logging.Logger:
    return logging.getLogger(name)


# ---------------------------------------------------------------------------
# PLINK wrapper
# ---------------------------------------------------------------------------


class PlinkError(RuntimeError):
    """Raised when PLINK returns a nonzero exit code."""


@dataclass
class PlinkResult:
    cmd: list[str]
    stdout: str
    stderr: str
    returncode: int


def run_plink(
    args: Sequence[str],
    binary: str = "plink",
    check: bool = True,
    cwd: str | os.PathLike[str] | None = None,
) -> PlinkResult:
    """Run a PLINK command and capture stdout/stderr.

    Parameters
    ----------
    args : sequence of str
        Arguments to pass to PLINK (without the binary name).
    binary : str
        PLINK executable. Either ``"plink"`` or ``"plink2"``, or an absolute path.
    check : bool
        If True (default), raise PlinkError on nonzero exit.
    cwd : path-like, optional
        Working directory for the subprocess.

    Returns
    -------
    PlinkResult
    """
    if shutil.which(binary) is None and not Path(binary).exists():
        raise PlinkError(f"PLINK binary not found on PATH: {binary}")

    cmd = [binary, *list(args)]
    logger = get_logger()
    logger.debug("Running: %s", " ".join(cmd))
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    result = PlinkResult(
        cmd=cmd,
        stdout=proc.stdout,
        stderr=proc.stderr,
        returncode=proc.returncode,
    )
    if check and proc.returncode != 0:
        logger.error("PLINK failed (exit %d): %s", proc.returncode, " ".join(cmd))
        logger.error("stderr: %s", proc.stderr)
        raise PlinkError(
            f"PLINK command failed (exit {proc.returncode}): {' '.join(cmd)}\n"
            f"stderr:\n{proc.stderr}"
        )
    return result


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def read_population_manifest(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read the HAPNEST per-individual superpopulation manifest.

    Expected columns: FID, IID, superpopulation.
    """
    df = pd.read_csv(path, sep="\t", dtype={"FID": str, "IID": str, "superpopulation": str})
    expected = {"FID", "IID", "superpopulation"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"Manifest {path} missing columns: {missing}")
    return df


def write_ids_file(df: pd.DataFrame, path: str | os.PathLike[str]) -> None:
    """Write a PLINK --keep style ID file (FID\\tIID, no header)."""
    # Several locus shards publish the same cohort lists. Atomic replacement
    # prevents another shard's PLINK process from reading a truncated keep file.
    import tempfile
    destination = Path(path)
    with tempfile.NamedTemporaryFile(mode="w", dir=destination.parent,
                                     prefix=f".{destination.name}.", delete=False) as handle:
        temporary = Path(handle.name)
        try:
            df[["FID", "IID"]].to_csv(handle, sep="\t", index=False, header=False)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(destination)


def write_site_manifest(df: pd.DataFrame, path: str | os.PathLike[str]) -> None:
    """Write a per-site individual → superpopulation manifest (TSV with header)."""
    df[["FID", "IID", "superpopulation"]].to_csv(path, sep="\t", index=False)


def read_site_manifest(path: str | os.PathLike[str]) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"FID": str, "IID": str, "superpopulation": str})


def ensure_dir(path: str | os.PathLike[str]) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def chunks(seq: Sequence[Any], size: int) -> Iterable[Sequence[Any]]:
    """Yield successive chunks of ``size`` from ``seq``."""
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


# ---------------------------------------------------------------------------
# Minimal PLINK 1 bed/bim/fam reader
# ---------------------------------------------------------------------------

BIM_COLUMNS = ["chrom", "snp_id", "cm", "bp", "a1", "a2"]
FAM_COLUMNS = ["FID", "IID", "father", "mother", "sex", "pheno"]


def read_bim(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read a PLINK .bim file. Standard 6 columns, no header.

    Returns a DataFrame indexed 0..N-1 where the row index is the variant
    index inside the matching .bed file.
    """
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=BIM_COLUMNS,
        dtype={"chrom": str, "snp_id": str, "cm": float, "bp": np.int64, "a1": str, "a2": str},
        engine="python",
    )
    return df


def read_fam(path: str | os.PathLike[str]) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=FAM_COLUMNS,
        dtype={"FID": str, "IID": str, "father": str, "mother": str, "sex": str, "pheno": str},
        engine="python",
    )
    return df


# PLINK 1 .bed bit code → dosage of A1 allele.
# Codes: 00=hom A1 (2), 01=missing (NaN), 10=het (1), 11=hom A2 (0).
_BED_CODE_TO_DOSAGE = np.array([2, np.nan, 1, 0], dtype=np.float32)


def read_bed_variants(
    bed_prefix: str | os.PathLike[str],
    variant_indices: Sequence[int] | np.ndarray,
    n_samples: int,
) -> np.ndarray:
    """Read a subset of variants from a PLINK 1 .bed file as a dosage matrix.

    Parameters
    ----------
    bed_prefix : path-like
        Path to the binary without the ``.bed`` suffix.
    variant_indices : sequence of int
        Zero-based variant indices (rows in the .bim file) to extract.
    n_samples : int
        Number of samples in the .fam file. Required to compute the per-variant
        byte block size.

    Returns
    -------
    np.ndarray of shape (n_samples, len(variant_indices)), dtype float32.
        Dosage of the A1 allele. Missing genotypes are NaN.

    Notes
    -----
    Only variant-major (mode byte 0x01) .bed files are supported. PLINK 1.9+
    always writes variant-major; this is the universal format.
    """
    bed_path = Path(str(bed_prefix) + ".bed")
    if not bed_path.exists():
        raise FileNotFoundError(f"No .bed file at {bed_path}")

    bytes_per_variant = (n_samples + 3) // 4
    variant_indices = np.asarray(variant_indices, dtype=np.int64)
    n_variants = len(variant_indices)

    out = np.empty((n_samples, n_variants), dtype=np.float32)

    with bed_path.open("rb") as fh:
        magic = fh.read(3)
        if len(magic) != 3:
            raise ValueError(f"Truncated .bed header in {bed_path}")
        if magic[0] != 0x6C or magic[1] != 0x1B:
            raise ValueError(f"Bad .bed magic in {bed_path}: {magic!r}")
        if magic[2] != 0x01:
            raise ValueError(
                f"Sample-major .bed not supported ({bed_path}); "
                f"re-encode with `plink --bfile <x> --make-bed --out <x>`"
            )

        for out_col, vidx in enumerate(variant_indices):
            offset = 3 + int(vidx) * bytes_per_variant
            fh.seek(offset)
            block = fh.read(bytes_per_variant)
            if len(block) != bytes_per_variant:
                raise ValueError(
                    f"Truncated variant {vidx} in {bed_path}: "
                    f"expected {bytes_per_variant} bytes, got {len(block)}"
                )
            arr = np.frombuffer(block, dtype=np.uint8)
            # Unpack 4 samples per byte (LSB first).
            codes = np.empty(bytes_per_variant * 4, dtype=np.uint8)
            codes[0::4] = arr & 0x03
            codes[1::4] = (arr >> 2) & 0x03
            codes[2::4] = (arr >> 4) & 0x03
            codes[3::4] = (arr >> 6) & 0x03
            out[:, out_col] = _BED_CODE_TO_DOSAGE[codes[:n_samples]]

    return out

