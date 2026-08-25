"""Step 0: obtain the input genotype cohort.

WHY THIS STEP EXISTS -- read this before using the output for anything published
--------------------------------------------------------------------------------
The GWAS pipeline used to begin at step 1, reading a ~6 GB PLINK1 triple of 100,000
synthetic European-ancestry samples over ~240,000 variants. That file was referenced only
as a path constant. Nothing in the code or the documentation said where it came from.

It was generated in a Google Colab notebook that was never version-controlled and is no
longer recoverable. So the original cohort is **not regenerable**, and a stranger
starting from this repository could not run the GWAS experiment at all -- the pipeline
started one step too late.

This module closes that gap honestly rather than pretending it is closed. It offers two
providers:

``plink_prefix``
    Use an existing PLINK1 triple. This is how the original cohort is used, by anyone who
    has a copy. It is reproducible in the weak sense that the same file gives the same
    answers, and not reproducible in the sense that matters: nobody can regenerate it.

``synthetic_ld``  (default)
    Generate a cohort from the documented model below. Deterministic from a seed, needs
    no download, and scales from a CI-sized fixture to full size. This makes the whole
    chain -- from nothing to per-site data -- runnable from this repository alone.

**The substitute is not the original.** It is a different cohort from a different
generative process, so it does not reproduce the published figures and must not be
described as doing so. It reproduces the *pipeline*, not the *numbers*.
docs/experiments/gwas/DATA.md states this plainly, and the run manifest records which
provider was used so no output is ever ambiguous about its own provenance.

THE GENERATIVE MODEL
--------------------
Deliberately simple and fully described here, because a substitute whose own provenance
is unclear would repeat the problem it exists to solve.

1. **Allele frequencies.** Each variant's minor-allele frequency is drawn from a Beta
   distribution and clipped. Beta(0.8, 2.4) is skewed toward low frequencies, which is
   the qualitative shape of a real site-frequency spectrum -- most variants are rare.

2. **Linkage disequilibrium.** Real genotypes are correlated in blocks; independent
   variants would make every association test cleanly independent and the simulation
   unrealistically easy. Variants are grouped into blocks, and within a block each
   variant's latent liability is a mixture of a shared block factor and its own noise,
   weighted so the within-block correlation is approximately ``ld_correlation``.

3. **Genotypes.** Each latent value is thresholded against the variant's allele frequency
   under a liability-threshold model to give an allele count of 0, 1 or 2.

4. **Missingness.** A small fraction of calls are set missing, so the mean-imputation
   path downstream is genuinely exercised rather than dead code.

This is a coarse model of population structure. It is adequate for exercising and
demonstrating the pipeline, and it is not a substitute for real or well-validated
synthetic genotypes if the goal is a claim about genetics rather than about federation.
"""

from __future__ import annotations

import logging
import struct
from pathlib import Path

import numpy as np

from appfl_bio_suite.experiments.gwas.simulation.schema import CohortParams

__all__ = ["materialize_cohort", "generate_synthetic_cohort", "write_plink1", "PLINK_STEM"]

log = logging.getLogger(__name__)

# Kept from the original cohort so that already-distributed bundles and the client loader
# continue to agree on filenames.
PLINK_STEM = "EUR.synthetic.100k.ld.maf"

# PLINK1 .bed packs two bits per call, in a deliberately non-obvious order:
#   00 -> homozygous A1   01 -> missing   10 -> heterozygous   11 -> homozygous A2
_DOSAGE_TO_BITS = {0: 0b00, 1: 0b10, 2: 0b11}
_MISSING_BITS = 0b01


def _ext(prefix: Path, extension: str) -> Path:
    """Append an extension to a PLINK prefix.

    NOT ``Path.with_suffix()``: the stem ends in ".maf", which with_suffix() treats as a
    suffix to replace, silently writing a differently-named file.
    """
    return Path(str(prefix) + extension)


def write_plink1(
    prefix: Path,
    genotypes: np.ndarray,
    chroms: list[str] | np.ndarray,
    positions: list[int] | np.ndarray,
    sample_ids: list[str] | None = None,
) -> None:
    """Write a PLINK1 .bed/.bim/.fam triple.

    ``genotypes`` is (n_samples, n_variants) of allele counts 0/1/2, with -1 for missing.
    Written SNP-major, which is what PLINK1 and pandas-plink expect.
    """
    n_samples, n_variants = genotypes.shape
    prefix.parent.mkdir(parents=True, exist_ok=True)

    ids = sample_ids or [f"IID{i + 1}" for i in range(n_samples)]

    with open(_ext(prefix, ".fam"), "w", encoding="utf-8") as handle:
        for i, iid in enumerate(ids):
            handle.write(f"FID{i + 1} {iid} 0 0 0 -9\n")

    with open(_ext(prefix, ".bim"), "w", encoding="utf-8") as handle:
        for j in range(n_variants):
            handle.write(f"{chroms[j]}\trs{j + 1}\t0\t{positions[j]}\tA\tG\n")

    with open(_ext(prefix, ".bed"), "wb") as handle:
        handle.write(struct.pack("BBB", 0x6C, 0x1B, 0x01))
        for j in range(n_variants):
            column = genotypes[:, j]
            packed = bytearray()
            for start in range(0, n_samples, 4):
                byte = 0
                for offset in range(4):
                    index = start + offset
                    if index >= n_samples:
                        break
                    value = int(column[index])
                    byte |= (_MISSING_BITS if value < 0 else _DOSAGE_TO_BITS[value]) << (2 * offset)
                packed.append(byte)
            handle.write(bytes(packed))


def generate_synthetic_cohort(params: CohortParams, out_prefix: Path) -> Path:
    """Generate a PLINK1 cohort with LD block structure. Deterministic from the seed."""
    rng = np.random.default_rng(params.seed)
    n = params.n_samples
    m = params.n_variants

    log.info("generating cohort: %s samples x %s variants", f"{n:,}", f"{m:,}")

    # 1. Allele frequencies, skewed toward rare as a real spectrum is.
    freqs = np.clip(
        rng.beta(params.maf_beta_a, params.maf_beta_b, m), params.maf_min, params.maf_max
    )

    # 2/3. LD blocks, then thresholding to allele counts.
    #
    # Within a block, each variant's latent value is
    #     sqrt(rho) * shared_block_factor  +  sqrt(1 - rho) * own_noise
    # so that any two variants in the block correlate at approximately rho. The
    # coefficients are square roots because correlation composes through the variance.
    genotypes = np.empty((n, m), dtype=np.int8)
    rho = float(np.clip(params.ld_correlation, 0.0, 0.999))
    shared_w = np.sqrt(rho)
    own_w = np.sqrt(1.0 - rho)
    block = max(1, params.ld_block_size)

    for start in range(0, m, block):
        stop = min(start + block, m)
        width = stop - start
        block_factor = rng.standard_normal((n, 1))
        own = rng.standard_normal((n, width))
        latent = shared_w * block_factor + own_w * own

        # Two independent haplotype draws per individual, thresholded at the frequency's
        # normal quantile. Summing them gives an allele count in {0, 1, 2} whose
        # distribution is Binomial(2, freq) marginally, while retaining the block
        # correlation.
        for k in range(width):
            j = start + k
            threshold = _norm_ppf(1.0 - freqs[j])
            hap_a = (latent[:, k] > threshold).astype(np.int8)
            second = shared_w * block_factor[:, 0] + own_w * rng.standard_normal(n)
            hap_b = (second > threshold).astype(np.int8)
            genotypes[:, j] = hap_a + hap_b

    # 4. Missingness, so the downstream mean-imputation path is exercised.
    if params.missing_rate > 0:
        n_missing = int(round(n * m * params.missing_rate))
        if n_missing:
            rows = rng.integers(0, n, n_missing)
            cols = rng.integers(0, m, n_missing)
            genotypes[rows, cols] = -1

    # Positions spread across chromosomes, GRCh37-style numeric coding.
    per_chrom = int(np.ceil(m / max(1, params.n_chromosomes)))
    chroms, positions = [], []
    for j in range(m):
        chroms.append(str(min(j // per_chrom + 1, params.n_chromosomes)))
        positions.append(100_000 + (j % per_chrom) * 1_500)

    write_plink1(out_prefix, genotypes, chroms, positions)
    log.info("cohort written: %s", out_prefix)
    return out_prefix


def _norm_ppf(q: np.ndarray | float):
    """Standard-normal inverse CDF.

    Uses scipy when available and falls back to an Acklam rational approximation
    otherwise, so that cohort generation does not force scipy into an environment that
    does not otherwise need it.
    """
    try:
        from scipy.stats import norm

        return norm.ppf(q)
    except ImportError:
        q = np.asarray(q, dtype=np.float64)
        a = [-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
             1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00]
        b = [-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
             6.680131188771972e01, -1.328068155288572e01]
        c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
             -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00]
        d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
             3.754408661907416e00]
        plow, phigh = 0.02425, 1 - 0.02425
        out = np.empty_like(q)
        lo, hi = q < plow, q > phigh
        mid = ~(lo | hi)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.sqrt(-2 * np.log(np.where(lo, q, plow)))
            out = np.where(
                lo,
                (((((c[0]*r+c[1])*r+c[2])*r+c[3])*r+c[4])*r+c[5])
                / ((((d[0]*r+d[1])*r+d[2])*r+d[3])*r+1),
                out,
            )
            r = np.sqrt(-2 * np.log(np.where(hi, 1 - q, plow)))
            out = np.where(
                hi,
                -(((((c[0]*r+c[1])*r+c[2])*r+c[3])*r+c[4])*r+c[5])
                / ((((d[0]*r+d[1])*r+d[2])*r+d[3])*r+1),
                out,
            )
            qm = np.where(mid, q, 0.5) - 0.5
            r = qm * qm
            out = np.where(
                mid,
                (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*qm
                / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1),
                out,
            )
        return out


def materialize_cohort(params: CohortParams, work_dir: Path) -> Path:
    """Resolve the cohort for a run and return its PLINK prefix (no extension).

    Either points at an existing triple or generates one, depending on the provider.
    """
    work_dir = Path(work_dir)

    if params.provider == "plink_prefix":
        prefix = Path(params.plink_prefix).expanduser()
        missing = [
            str(_ext(prefix, e)) for e in (".bed", ".bim", ".fam") if not _ext(prefix, e).is_file()
        ]
        if missing:
            raise FileNotFoundError(
                "cohort provider is 'plink_prefix' but the triple is incomplete:\n  - "
                + "\n  - ".join(missing)
                + "\n\nThe prefix must NOT include a file extension.\n"
                "If you do not have this cohort, switch the scenario to provider "
                "'synthetic_ld', which generates one. See docs/experiments/gwas/DATA.md."
            )
        log.info("using existing cohort: %s", prefix)
        return prefix

    if params.provider == "synthetic_ld":
        prefix = work_dir / "cohort" / PLINK_STEM
        if all(_ext(prefix, e).is_file() for e in (".bed", ".bim", ".fam")):
            log.info("cohort already present, reusing: %s", prefix)
            return prefix
        return generate_synthetic_cohort(params, prefix)

    raise ValueError(
        f"unknown cohort provider '{params.provider}'. Known: synthetic_ld, plink_prefix."
    )
