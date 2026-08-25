"""Step 0: obtain the genotype pool the site cohorts are drawn from.

WHY THIS STEP EXISTS
--------------------
The upstream pipeline begins at step 1, subsampling a staged copy of HAPNEST. Staging
HAPNEST means a 135 GB download for chromosome 1 alone (1.5 TB for the full set), which
is entirely reasonable for the production run and completely unreasonable as the price of
finding out whether your install works.

So this module offers the same two-provider shape the GWAS experiment uses:

``hapnest``
    Use the staged real pool. This is what the published run used, and the only provider
    whose output means anything about genetics.

``synthetic_hapnest``  (default)
    Generate a small pool with HAPNEST's *shape* -- the same file layout, the same
    superpopulation labelling, the same variant-major PLINK1 encoding -- from the model
    documented below. Deterministic from a seed, needs no download, and runs in seconds.

**The substitute is not HAPNEST.** It does not reproduce the published figures and must
not be described as doing so; it reproduces the pipeline, not the numbers. The run
manifest records the provider, so no output is ambiguous about its own provenance.

THE GENERATIVE MODEL
--------------------
Described in full here, because a substitute whose own provenance is unclear repeats the
problem it exists to solve. What matters for *this* experiment is different from what
mattered for the GWAS one: fine-mapping is entirely a statement about correlation
structure, and credible sets are only interesting when ancestries disagree about it. A
pool with one shared LD pattern would make the low/medium/high divergence strata
indistinguishable and the whole locus-selection stage vacuous.

1. **Ancestry labels.** Individuals are assigned to the six HAPNEST superpopulations by
   the configured weights.

2. **Ancestry-differentiated allele frequencies.** Each variant gets an ancestral
   frequency from a Beta spectrum skewed toward rare, then a per-superpopulation
   frequency drawn from the Balding-Nichols model::

       p_pop ~ Beta( p (1-F)/F , (1-p)(1-F)/F )

   ``F`` plays the role of Fst and comes from ``ld_divergence``. This is the standard
   model for population differentiation, and it is what makes allele frequencies -- and
   therefore LD -- genuinely differ between ancestries rather than by simulation noise.

3. **Ancestry-differentiated LD blocks.** Variants are grouped into blocks. Within a
   block each variant's latent liability mixes a shared block factor with its own noise,
   and the mixing weight is itself perturbed per superpopulation. So two ancestries
   differ in *how correlated* a block is, not only in its frequencies -- which is the
   divergence the locus-selection stage scores with a Frobenius norm.

4. **Genotypes.** Two haplotype draws per individual, thresholded against that
   individual's own superpopulation frequency, summed to an allele count in {0, 1, 2}.

5. **Missingness.** A small fraction of calls are set missing, so the complete-case
   handling in the site stage is genuinely exercised rather than dead code.

This is a coarse model of population structure. It is adequate for exercising and
demonstrating the pipeline, and it is not a substitute for HAPNEST if the goal is a claim
about genetics rather than about federation.
"""

from __future__ import annotations

import logging
import struct
from pathlib import Path

import numpy as np

from appfl_bio_suite.experiments.fine_mapping.simulation.schema import CohortParams

__all__ = ["materialize_pool", "generate_synthetic_pool", "write_plink1", "POOL_STEM"]

log = logging.getLogger(__name__)


def POOL_STEM(chromosome: int) -> str:
    """The filename stem the pipeline expects, ``chr<N>``.

    Not configurable: ``fedfm/sampling.py`` builds this name itself and says so in a
    comment. Deriving it here rather than repeating the literal keeps the two in step.
    """
    return f"chr{int(chromosome)}"


# PLINK1 .bed packs two bits per call, in a deliberately non-obvious order:
#   00 -> homozygous A1   01 -> missing   10 -> heterozygous   11 -> homozygous A2
# Indexed by dosage 0/1/2, with missing handled separately.
_DOSAGE_TO_CODE = np.array([0b11, 0b10, 0b00], dtype=np.uint8)
_MISSING_CODE = np.uint8(0b01)


def write_plink1(prefix: Path, genotypes, chrom, positions, sample_ids, alleles=None) -> Path:
    """Write a variant-major PLINK1 .bed/.bim/.fam triple.

    ``genotypes`` is (n_samples, n_variants) of allele counts 0/1/2, with -1 for missing.

    The bit packing is vectorized rather than looped. That is not premature optimization:
    the obvious nested-loop version costs one Python-level operation per (sample, variant)
    pair, which at even a modest fixture size is minutes of pure interpreter overhead in
    the middle of a command whose whole purpose is to be quick enough to run casually.
    """
    genotypes = np.asarray(genotypes)
    n_samples, n_variants = genotypes.shape
    prefix.parent.mkdir(parents=True, exist_ok=True)

    with open(_ext(prefix, ".fam"), "w", encoding="utf-8") as handle:
        for iid in sample_ids:
            handle.write(f"{iid} {iid} 0 0 0 -9\n")

    a1a2 = alleles if alleles is not None else [("A", "G")] * n_variants
    with open(_ext(prefix, ".bim"), "w", encoding="utf-8") as handle:
        for j in range(n_variants):
            a1, a2 = a1a2[j]
            handle.write(f"{chrom}\trs{j + 1}\t0\t{int(positions[j])}\t{a1}\t{a2}\n")

    # Pad the sample axis up to a multiple of 4; PLINK ignores the trailing bits, and
    # padding lets the whole column reshape cleanly into bytes.
    padded = ((n_samples + 3) // 4) * 4
    with open(_ext(prefix, ".bed"), "wb") as handle:
        handle.write(struct.pack("BBB", 0x6C, 0x1B, 0x01))
        for j in range(n_variants):
            column = genotypes[:, j]
            codes = np.empty(padded, dtype=np.uint8)
            codes[:n_samples] = np.where(
                column < 0, _MISSING_CODE, _DOSAGE_TO_CODE[np.clip(column, 0, 2)]
            )
            codes[n_samples:] = _MISSING_CODE
            quads = codes.reshape(-1, 4)
            packed = quads[:, 0] | (quads[:, 1] << 2) | (quads[:, 2] << 4) | (quads[:, 3] << 6)
            handle.write(packed.astype(np.uint8).tobytes())

    return prefix


def _ext(prefix: Path, extension: str) -> Path:
    return Path(str(prefix) + extension)


def _assign_superpopulations(params: CohortParams, rng) -> np.ndarray:
    """Ancestry label per individual, in the order they appear in the .fam.

    Allocated by rounded proportion rather than by multinomial draw, so a scenario asking
    for a specific composition gets it exactly. A multinomial would make the smallest
    superpopulation's count wobble between runs of the same seed-independent size, and
    the site compositions are checked against these counts.
    """
    pops = list(params.superpopulation_weights)
    weights = np.array([params.superpopulation_weights[p] for p in pops], dtype=np.float64)
    weights = weights / weights.sum()

    counts = np.floor(weights * params.n_individuals).astype(int)
    # Hand out the rounding remainder to the largest groups, largest first.
    shortfall = params.n_individuals - counts.sum()
    for i in np.argsort(-weights)[:shortfall]:
        counts[i] += 1

    labels = np.repeat(np.array(pops, dtype=object), counts)
    rng.shuffle(labels)
    return labels


def _balding_nichols(ancestral, fst, rng):
    """Per-population allele frequencies around an ancestral frequency.

    ``Beta(p(1-F)/F, (1-p)(1-F)/F)`` has mean ``p`` and variance ``p(1-p)F``, so ``F`` is
    exactly the between-population differentiation. Clipped away from 0 and 1: a variant
    monomorphic in one ancestry contributes a zero-variance column, which the site stage
    handles but which is not what a divergence knob is meant to produce.
    """
    fst = float(np.clip(fst, 1e-4, 0.5))
    a = ancestral * (1.0 - fst) / fst
    b = (1.0 - ancestral) * (1.0 - fst) / fst
    return np.clip(rng.beta(a, b), 0.005, 0.995)


def generate_synthetic_pool(params: CohortParams, chromosome: int, out_dir: Path) -> Path:
    """Generate a HAPNEST-shaped pool. Deterministic from ``params.seed``."""
    from scipy.stats import norm

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(params.seed)

    n, m = int(params.n_individuals), int(params.n_variants)
    log.info("generating pool: %s individuals x %s variants", f"{n:,}", f"{m:,}")

    labels = _assign_superpopulations(params, rng)
    pops = sorted(set(labels.tolist()))
    pop_rows = {pop: np.flatnonzero(labels == pop) for pop in pops}

    # 2. Frequencies: one ancestral spectrum, then per-ancestry differentiation.
    ancestral = np.clip(
        rng.beta(params.maf_beta_a, params.maf_beta_b, m), params.maf_min, params.maf_max
    )
    freq = {pop: _balding_nichols(ancestral, params.ld_divergence, rng) for pop in pops}

    # 3. Block correlation, perturbed per ancestry so LD *structure* diverges too.
    block = max(1, int(params.ld_block_size))
    n_blocks = (m + block - 1) // block
    base_rho = float(np.clip(params.ld_correlation, 0.0, 0.99))
    spread = float(np.clip(params.ld_divergence, 0.0, 1.0)) * 0.5
    rho = {
        pop: np.clip(base_rho + rng.uniform(-spread, spread, n_blocks), 0.05, 0.98) for pop in pops
    }

    genotypes = np.empty((n, m), dtype=np.int8)
    for pop in pops:
        rows = pop_rows[pop]
        if rows.size == 0:
            continue
        thresholds = norm.ppf(1.0 - freq[pop])
        for b in range(n_blocks):
            start, stop = b * block, min((b + 1) * block, m)
            width = stop - start
            shared_w = np.sqrt(rho[pop][b])
            own_w = np.sqrt(1.0 - rho[pop][b])
            # Two haplotypes per individual: independent own-noise, shared block factor
            # drawn once per haplotype so the block correlation survives the sum.
            counts = np.zeros((rows.size, width), dtype=np.int8)
            for _ in range(2):
                factor = rng.standard_normal((rows.size, 1))
                latent = shared_w * factor + own_w * rng.standard_normal((rows.size, width))
                counts += (latent > thresholds[start:stop]).astype(np.int8)
            genotypes[np.ix_(rows, np.arange(start, stop))] = counts

    # 5. Missingness, so the complete-case path in the site stage is exercised.
    if params.missing_rate > 0:
        n_missing = int(round(n * m * params.missing_rate))
        if n_missing:
            genotypes[rng.integers(0, n, n_missing), rng.integers(0, m, n_missing)] = -1

    sample_ids = [f"SYN{i + 1:07d}" for i in range(n)]
    positions = params.bp_start + np.arange(m, dtype=np.int64) * int(params.bp_spacing)

    prefix = out_dir / POOL_STEM(chromosome)
    write_plink1(prefix, genotypes, chromosome, positions, sample_ids)

    manifest_path = out_dir / "population_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8") as handle:
        handle.write("FID\tIID\tsuperpopulation\n")
        for iid, pop in zip(sample_ids, labels, strict=True):
            handle.write(f"{iid}\t{iid}\t{pop}\n")

    log.info("pool written: %s (+ %s)", prefix, manifest_path.name)
    return prefix


def materialize_pool(params: CohortParams, chromosome: int, hapnest_dir: Path) -> Path:
    """Resolve the genotype pool for a run and return its PLINK prefix.

    Either points the pipeline at a staged HAPNEST copy or generates a substitute in the
    run's own directory.
    """
    hapnest_dir = Path(hapnest_dir)

    if params.provider == "hapnest":
        staged = params.resolved_hapnest_dir()
        if staged is None:
            raise FileNotFoundError(
                "cohort provider is 'hapnest' but no staged pool is configured. Set "
                "`cohort.hapnest_dir` in your copy of the scenario, or $FEDFM_HAPNEST_DIR "
                "in the environment.\nscripts/fine-mapping/download_hapnest.sh stages it; "
                "see docs/experiments/fine-mapping/DATA.md."
            )
        prefix = staged / POOL_STEM(chromosome)
        missing = [
            str(_ext(prefix, e)) for e in (".bed", ".bim", ".fam") if not _ext(prefix, e).is_file()
        ]
        if not (staged / "population_manifest.tsv").is_file():
            missing.append(str(staged / "population_manifest.tsv"))
        if missing:
            raise FileNotFoundError(
                "cohort provider is 'hapnest' but the staged pool is incomplete:\n  - "
                + "\n  - ".join(missing)
                + "\n\nStage it with scripts/fine-mapping/download_hapnest.sh, which "
                "fetches the EBI BioStudies accession S-BSST936 and lays out exactly "
                "these names. If you do not have it, switch the scenario to provider "
                "'synthetic_hapnest', which generates a substitute.\n"
                "See docs/experiments/fine-mapping/DATA.md."
            )
        log.info("using staged HAPNEST pool: %s", prefix)
        return prefix

    if params.provider == "synthetic_hapnest":
        prefix = hapnest_dir / POOL_STEM(chromosome)
        if (
            all(_ext(prefix, e).is_file() for e in (".bed", ".bim", ".fam"))
            and (hapnest_dir / "population_manifest.tsv").is_file()
        ):
            log.info("pool already present, reusing: %s", prefix)
            return prefix
        return generate_synthetic_pool(params, chromosome, hapnest_dir)

    raise ValueError(
        f"unknown cohort provider '{params.provider}'. Known: synthetic_hapnest, hapnest."
    )
