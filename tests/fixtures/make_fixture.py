"""Build a tiny, deterministic PLINK1 cohort plus PGS weight files, for tests.

Generated rather than committed, for two reasons: the repo's .gitignore excludes PLINK
extensions on purpose (real genotype data must never be committable by accident), and a
generated fixture cannot drift from the generator that documents it.

The fixture is deliberately small enough to run the whole simulation pipeline in a couple
of seconds, and deliberately structured like the real cohort -- same file naming, same
PGS column layout, same chromosome coding -- so that code exercised against it is
exercising the real code paths.

This is fixture-building only. It is NOT the substitute cohort a user runs the pipeline
on; that is ``experiments/gwas/simulation/cohort.py``, which is a documented scientific
artifact rather than a test convenience.
"""

from __future__ import annotations

import gzip
import struct
from pathlib import Path

import numpy as np

# Same stem as the real cohort. The legacy scripts hardcode this name, and the harness
# has to run them unmodified.
PLINK_STEM = "EUR.synthetic.100k.ld.maf"

# PLINK1 .bed packs two bits per genotype. The codes are not in a natural order:
#   00 -> homozygous A1      01 -> missing
#   10 -> heterozygous       11 -> homozygous A2
# pandas-plink with ref="a1" then reports the dosage of the SECOND allele, so
# 0/1/2 maps to 00/10/11 and missing maps to 01.
_DOSAGE_TO_BITS = {0: 0b00, 1: 0b10, 2: 0b11}
_MISSING_BITS = 0b01


def _ext(prefix: Path, extension: str) -> Path:
    """Append an extension to a PLINK prefix.

    Deliberately not Path.with_suffix(): the real cohort stem ends in ".maf", which
    with_suffix() treats as a suffix to replace rather than part of the name.
    """
    return Path(str(prefix) + extension)


def write_plink1(
    prefix: Path, genotypes: np.ndarray, chroms: list[str], positions: list[int]
) -> None:
    """Write a PLINK1 .bed/.bim/.fam triple.

    ``genotypes`` is (n_samples, n_variants) with values 0/1/2, or -1 for missing.
    Written SNP-major, which is what PLINK1 and pandas-plink both expect.
    """
    n_samples, n_variants = genotypes.shape
    prefix.parent.mkdir(parents=True, exist_ok=True)

    # NOTE: string concatenation, not Path.with_suffix(). The cohort stem is
    # "EUR.synthetic.100k.ld.maf", whose final dotted component is ".maf" -- so
    # with_suffix(".bed") REPLACES it and silently writes "EUR.synthetic.100k.ld.bed".
    # The legacy pipeline builds these paths the same way for the same reason.
    with open(_ext(prefix, ".fam"), "w", encoding="utf-8") as handle:
        for i in range(n_samples):
            handle.write(f"FID{i + 1} IID{i + 1} 0 0 0 -9\n")

    # .bim -- CHR SNP CM BP A1 A2
    with open(_ext(prefix, ".bim"), "w", encoding="utf-8") as handle:
        for j in range(n_variants):
            handle.write(f"{chroms[j]}\trs{j + 1}\t0\t{positions[j]}\tA\tG\n")

    # .bed -- magic, SNP-major mode byte, then packed genotypes
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
                    bits = _MISSING_BITS if value < 0 else _DOSAGE_TO_BITS[value]
                    byte |= bits << (2 * offset)
                packed.append(byte)
            handle.write(bytes(packed))


def write_pgs_catalog(path: Path, chroms, positions, weights, seed_label: str) -> None:
    """Write a gzipped PGS Catalog scoring file.

    The real files carry a block of '#'-prefixed metadata before the header, and the
    parser skips comments -- so the fixture carries one too. A fixture that omits it
    would not exercise that path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("###PGS CATALOG SCORING FILE - see www.pgscatalog.org\n")
        handle.write(f"#pgs_id={seed_label}\n")
        handle.write("#genome_build=GRCh37\n")
        handle.write("rsID\thm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n")
        for i, (chrom, pos, weight) in enumerate(zip(chroms, positions, weights, strict=True)):
            # Alternate the effect allele so the allele-flip branch in the scoring code
            # is genuinely exercised rather than always taking the same path.
            effect, other = ("G", "A") if i % 2 == 0 else ("A", "G")
            handle.write(f"rs{i + 1}\t{chrom}\t{pos}\t{effect}\t{other}\t{weight:.6f}\n")


def build(
    root: Path,
    n_samples: int = 240,
    n_variants: int = 300,
    seed: int = 20260814,
) -> dict[str, Path]:
    """Create the full fixture under ``root`` and return the paths written.

    ``n_samples`` is a multiple of 4 so the .bed has no partial trailing byte in the
    common case, and large enough that regressions produce visible numeric differences.
    """
    rng = np.random.default_rng(seed)

    # Realistic-ish allele frequencies: a beta draw clipped away from monomorphic, so
    # every variant is actually polymorphic and MAF calculations are meaningful.
    freqs = np.clip(rng.beta(1.4, 1.4, n_variants), 0.05, 0.95)
    genotypes = np.empty((n_samples, n_variants), dtype=np.int8)
    for j, freq in enumerate(freqs):
        genotypes[:, j] = rng.binomial(2, freq, n_samples)

    # A little missingness, because the pipeline mean-imputes and that branch should run.
    n_missing = max(1, (n_samples * n_variants) // 500)
    rows = rng.integers(0, n_samples, n_missing)
    cols = rng.integers(0, n_variants, n_missing)
    genotypes[rows, cols] = -1

    # Spread across a few chromosomes with increasing positions, matching the real
    # cohort's GRCh37 coding.
    chroms, positions = [], []
    for j in range(n_variants):
        chrom = str((j % 5) + 1)
        chroms.append(chrom)
        positions.append(100_000 + j * 1_000)

    prefix = root / "data_sim" / "input" / PLINK_STEM
    write_plink1(prefix, genotypes, chroms, positions)

    # PGS files overlap the genotypes only partially, as real ones do -- the pipeline's
    # overlap step is load-bearing and a fully-overlapping fixture would not test it.
    overlap = int(n_variants * 0.6)
    idx = np.sort(rng.choice(n_variants, overlap, replace=False))
    t2d_path = prefix.parent / "PGS003443_hmPOS_GRCh37.txt.gz"
    bmi_path = prefix.parent / "PGS004994_hmPOS_GRCh37.txt.gz"

    write_pgs_catalog(
        t2d_path,
        [chroms[i] for i in idx],
        [positions[i] for i in idx],
        rng.normal(0, 0.35, overlap),
        "PGS003443",
    )
    write_pgs_catalog(
        bmi_path,
        [chroms[i] for i in idx],
        [positions[i] for i in idx],
        rng.normal(0, 0.30, overlap),
        "PGS004994",
    )

    return {
        "plink_prefix": prefix,
        "bed": _ext(prefix, ".bed"),
        "bim": _ext(prefix, ".bim"),
        "fam": _ext(prefix, ".fam"),
        "pgs_t2d": t2d_path,
        "pgs_bmi": bmi_path,
        "n_samples": n_samples,
        "n_variants": n_variants,
    }


if __name__ == "__main__":
    import sys

    target = Path(sys.argv[1] if len(sys.argv) > 1 else "./fixture")
    written = build(target)
    for key, value in written.items():
        print(f"{key}: {value}")
