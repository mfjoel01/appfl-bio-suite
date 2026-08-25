"""Step 1: simulate phenotypes, covariates, and polygenic scores from genotypes.

Reproduces the simulation design of "Simulation Benchmarking of Differential Privacy and
Federated Learning for Secure Genetic Analysis and Polygenic Risk Scoring" (Joel et al.).

=============================================================================
NUMERICAL BEHAVIOUR IS FROZEN. This code backs published figures.
=============================================================================

``tests/test_gwas_simulation_parity.py`` runs the pre-migration scripts and this port
against the same fixture and compares checksums of all 35 outputs. Any change that alters
a single byte fails that test.

That is deliberately strict, because the ways to break this silently are not obvious:

* **RNG stream order.** ``rng.normal(...)`` then ``rng.integers(...)`` consumes the
  generator in that order. Reordering the two draws, or hoisting one out of a loop,
  changes every value with no error.
* **Population vs. sample standard deviation.** ``ndarray.std()`` is ddof=0.
  ``pandas.Series.std()`` is ddof=1. Substituting one for the other rescales every
  polygenic score by sqrt(n/(n-1)).
* **Chunked accumulation.** The scoring dot product is chunked at 2000 variants.
  Floating-point addition is not associative, so changing the chunk size changes the
  low-order bits of every score.
* **Operation order in the PGS overlap.** The groupby-then-join-then-filter sequence
  determines which duplicate variants survive and in what row order.

What DID change from the original, all of it non-numerical: module-level constants became
an explicit configuration object; the ``GWAS_PROJECT_DIR`` environment variable and its
``sys.path.insert`` with a guessed relative fallback became ordinary package imports; and
the script became a function that takes its inputs and returns its outputs.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from appfl_bio_suite.experiments.gwas.simulation.schema import PhenotypeParams

__all__ = ["simulate_phenotypes", "apply_data_sim_scaling", "PhenotypeOutputs"]

log = logging.getLogger(__name__)

# PGS Catalog uses letter codes for the non-autosomes; the BIM uses numbers.
CHR_NORM = {"X": "23", "Y": "24", "XY": "25", "MT": "26", "M": "26"}

# Written in this order by the original. Order matters only for readability, but the
# names and contents are the loader's and the bundler's contract.
OUTPUT_FILES = (
    "covariates.csv",
    "phenotypes_gwas.csv",
    "phenotypes_pgs_eval.csv",
    "pgs_scores.csv",
    "simulation_summary.csv",
    "variant_counts.csv",
    "score_t2d.txt",
    "score_bmi.txt",
)


class PhenotypeOutputs(dict):
    """Result paths, keyed by filename."""


def _ext(prefix: Path, extension: str) -> Path:
    """Append an extension to a PLINK prefix.

    NOT ``Path.with_suffix()``: the cohort stem ends in ".maf", which with_suffix()
    treats as a suffix to replace, silently producing a different filename.
    """
    return Path(str(prefix) + extension)


def apply_data_sim_scaling(
    variant_df: pd.DataFrame, G=None, *, fraction: float = 1.0, seed: int = 7
):
    """Deterministically subsample variants at simulation time.

    Ported verbatim from the original ``gwas_config.apply_data_sim_scaling``, including
    the ``round``-then-``sort`` behaviour, because the chosen indices must match.

    Uses a different seed from the analysis-time variant scaling so the two can be
    applied independently and compound.
    """
    if fraction >= 1.0:
        return variant_df, G

    n_total = len(variant_df)
    n_keep = max(1, round(n_total * fraction))
    rng = np.random.default_rng(seed)
    chosen_idx = np.sort(rng.choice(n_total, size=n_keep, replace=False))

    filtered_df = variant_df.iloc[chosen_idx].reset_index(drop=True)
    filtered_G = G.isel(variant=chosen_idx) if G is not None else None
    return filtered_df, filtered_G


def _load_pgs(path: Path) -> pd.DataFrame:
    """Parse a gzipped PGS Catalog scoring file.

    Every line starting with '#' is metadata; the first non-'#' line is the header.
    """
    return pd.read_csv(path, sep="\t", comment="#", compression="gzip", low_memory=False)


def _prep_pgs(df: pd.DataFrame, bim: pd.DataFrame, bim_set: set) -> pd.DataFrame:
    """Normalise PGS coordinates, match to the BIM, collapse duplicates, check alleles.

    The order of these steps is part of the frozen behaviour: collapsing duplicates
    before joining the BIM, and filtering on allele match after the join, determines
    which variants survive.
    """
    df = df.copy()

    df["effect_weight"] = pd.to_numeric(df["effect_weight"], errors="coerce")

    df["hm_chr"] = df["hm_chr"].astype(str).str.strip()
    df["hm_pos"] = pd.to_numeric(df["hm_pos"], errors="coerce")
    df = df.dropna(subset=["hm_chr", "hm_pos", "effect_weight"])
    df["hm_pos"] = df["hm_pos"].astype(int)

    df["hm_chr"] = df["hm_chr"].replace(CHR_NORM)
    df["chr_pos"] = df["hm_chr"] + ":" + df["hm_pos"].astype(str)

    df = df[df["chr_pos"].isin(bim_set)].copy()

    # Sum the weights of duplicate (position, effect allele) pairs. A catalog can list
    # the same variant more than once.
    df = df.groupby(["chr_pos", "effect_allele"], as_index=False).agg(
        effect_weight=("effect_weight", "sum")
    )

    bim_sub = bim.set_index("chr_pos")[["SNP", "A1", "A2"]]
    df = df.join(bim_sub, on="chr_pos")

    # Drop variants whose effect allele is neither of the genotyped alleles -- a strand
    # or build mismatch, and scoring them would be meaningless.
    ea = df["effect_allele"].str.upper()
    a1 = df["A1"].str.upper()
    a2 = df["A2"].str.upper()
    df = df[ea.eq(a1) | ea.eq(a2)].copy()

    df = df.dropna(subset=["SNP"])
    return df.reset_index(drop=True)


def _compute_pgs(pgs_df: pd.DataFrame, G, snp_to_pos: dict, chunk_size: int = 2000) -> np.ndarray:
    """Weighted allele-dosage sum, chunked.

    Dosage encoding with ``ref="a1"`` is the count of the SECOND allele (the BIM's A2):
    0 = A1/A1, 1 = A1/A2, 2 = A2/A2, NaN = missing. If the PGS effect allele is A1, the
    dosage is flipped to ``2 - G``.

    Missing genotypes are mean-imputed per variant, within the chunk. That is what the
    original does, and the per-chunk (rather than global) mean is part of the frozen
    behaviour -- with a fixed chunk size the two coincide anyway, since chunks partition
    variants, not samples.

    ``chunk_size`` changes floating-point accumulation order and therefore the result.
    """
    ea = pgs_df["effect_allele"].str.upper().values
    a1v = pgs_df["A1"].str.upper().values
    wts = pgs_df["effect_weight"].values.astype(np.float64)
    snps = pgs_df["SNP"].values

    positions, keep = [], []
    for i, snp in enumerate(snps):
        pos = snp_to_pos.get(snp)
        if pos is not None:
            positions.append(pos)
            keep.append(i)

    ea = ea[keep]
    a1v = a1v[keep]
    wts = wts[keep]
    flip = ea == a1v

    scores = np.zeros(G.sizes["sample"], dtype=np.float64)

    for start in range(0, len(positions), chunk_size):
        sl = slice(start, start + chunk_size)
        idx = positions[sl]
        w = wts[sl]
        flip_sl = flip[sl]

        chunk = G.isel(variant=idx).values.astype(np.float64)

        if flip_sl.any():
            chunk[:, flip_sl] = 2.0 - chunk[:, flip_sl]

        col_means = np.nanmean(chunk, axis=0)
        r, c = np.where(np.isnan(chunk))
        chunk[r, c] = col_means[c]

        scores += chunk @ w

    return scores


def _standardise(v: np.ndarray) -> np.ndarray:
    """Z-score using the POPULATION standard deviation (ddof=0).

    ``ndarray.std()`` is ddof=0; ``pandas.Series.std()`` is ddof=1. Using the latter
    would rescale every score by sqrt(n/(n-1)) with no error and no obvious symptom.
    """
    sd = v.std()
    return (v - v.mean()) / sd if sd > 0 else np.zeros_like(v)


def simulate_phenotypes(
    plink_prefix: Path,
    pgs_t2d_path: Path,
    pgs_bmi_path: Path,
    out_dir: Path,
    params: PhenotypeParams | None = None,
) -> PhenotypeOutputs:
    """Run step 1 and write its outputs to ``out_dir``.

    Args:
        plink_prefix: PLINK1 triple prefix, WITHOUT extension.
        pgs_t2d_path, pgs_bmi_path: gzipped PGS Catalog scoring files.
        out_dir: destination for the CSV and score outputs.
        params: simulation parameters; defaults reproduce the published design.
    """
    from pandas_plink import read_plink1_bin

    params = params or PhenotypeParams()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plink_prefix = Path(plink_prefix)

    # -- 1. FAM and BIM ----------------------------------------------------
    log.info("reading FAM")
    fam_raw = pd.read_csv(_ext(plink_prefix, ".fam"), sep=r"\s+", header=None, dtype=str)
    # Only the first two columns, regardless of trailing phenotype columns.
    fam = fam_raw.iloc[:, :2].copy()
    fam.columns = ["FID", "IID"]
    n = len(fam)
    log.info("  %s samples", f"{n:,}")

    log.info("reading BIM")
    bim = pd.read_csv(
        _ext(plink_prefix, ".bim"),
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 3, 4, 5],
        names=["CHR", "SNP", "BP", "A1", "A2"],
        dtype={"CHR": str, "SNP": str, "BP": str, "A1": str, "A2": str},
    )
    n_bim_full = len(bim)
    bim, _ = apply_data_sim_scaling(
        bim, fraction=params.data_sim_scaling, seed=params.data_sim_scaling_seed
    )
    if params.data_sim_scaling < 1.0:
        log.info(
            "  data_sim_scaling=%s: BIM reduced to %s / %s variants",
            params.data_sim_scaling,
            f"{len(bim):,}",
            f"{n_bim_full:,}",
        )
    bim["chr_pos"] = bim["CHR"].str.strip() + ":" + bim["BP"].str.strip()
    bim_chr_pos_set = set(bim["chr_pos"])
    log.info("  %s variants", f"{len(bim):,}")

    # -- 2. PGS catalog files ----------------------------------------------
    log.info("loading PGS files")
    with ThreadPoolExecutor(max_workers=2) as executor:
        fut_t2d = executor.submit(_load_pgs, pgs_t2d_path)
        fut_bmi = executor.submit(_load_pgs, pgs_bmi_path)
        raw_t2d = fut_t2d.result()
        raw_bmi = fut_bmi.result()
    n_t2d_total = len(raw_t2d)
    n_bmi_total = len(raw_bmi)

    pgs_t2d = _prep_pgs(raw_t2d, bim, bim_chr_pos_set)
    pgs_bmi = _prep_pgs(raw_bmi, bim, bim_chr_pos_set)
    log.info("  T2D overlap: %s   BMI overlap: %s", f"{len(pgs_t2d):,}", f"{len(pgs_bmi):,}")

    if len(pgs_t2d) == 0 or len(pgs_bmi) == 0:
        raise RuntimeError(
            "no PGS variants matched the genotypes.\n"
            "Almost always a chromosome-coding mismatch: the BIM uses numeric codes "
            "(1..22, then 23/24/25/26) and the PGS file's hm_chr must normalise to the "
            "same. Check a few rows of each."
        )

    # -- 3. plink2 --score input files -------------------------------------
    # Written for optional external use. Nothing in this pipeline reads them back;
    # scoring below is done in numpy. No plink2 binary is required.
    score_t2d_path = out_dir / "score_t2d.txt"
    score_bmi_path = out_dir / "score_bmi.txt"
    pgs_t2d[["SNP", "effect_allele", "effect_weight"]].to_csv(
        score_t2d_path, sep="\t", header=False, index=False
    )
    pgs_bmi[["SNP", "effect_allele", "effect_weight"]].to_csv(
        score_bmi_path, sep="\t", header=False, index=False
    )

    # -- 4. polygenic scores -----------------------------------------------
    log.info("loading genotypes (memory-mapped)")
    G = read_plink1_bin(str(_ext(plink_prefix, ".bed")), verbose=False, ref="a1")
    snp_to_pos = {snp: i for i, snp in enumerate(G.snp.values)}

    log.info("scoring T2D")
    g_raw_t2d = _compute_pgs(pgs_t2d, G, snp_to_pos, params.pgs_chunk_size)
    log.info("scoring BMI")
    g_raw_bmi = _compute_pgs(pgs_bmi, G, snp_to_pos, params.pgs_chunk_size)

    g_std_t2d = _standardise(g_raw_t2d)
    g_std_bmi = _standardise(g_raw_bmi)

    # -- 5. covariates -----------------------------------------------------
    # Draw order is part of the frozen behaviour: normal() then integers(), from one
    # generator. Swapping them changes both.
    log.info("simulating covariates")
    rng_cov = np.random.default_rng(params.seed_covariates)
    age_raw = rng_cov.normal(40, 15, n)
    age = np.clip(age_raw, 18, 90).round().astype(int)
    sex = rng_cov.integers(0, 2, n)

    # -- 6. phenotypes, two independent replicates --------------------------
    log.info("simulating phenotypes")
    sd_e_t2d = np.sqrt((1.0 - params.h2_t2d) / params.h2_t2d)
    sd_e_bmi = np.sqrt((1.0 - params.h2_bmi) / params.h2_bmi)
    threshold = 1.0 - params.t2d_prevalence

    def _draw(seed: int):
        rng = np.random.default_rng(seed)
        eps_t2d = rng.normal(0.0, sd_e_t2d, n)
        eps_bmi = rng.normal(0.0, sd_e_bmi, n)

        liability = (
            params.effect_scale_t2d * g_std_t2d
            + params.age_coef_t2d * age
            + params.sex_coef_t2d * sex
            + eps_t2d
        )
        t2d = (liability > np.quantile(liability, threshold)).astype(np.int8)

        bmi = (
            params.bmi_intercept
            + params.effect_scale_bmi * g_std_bmi
            + params.age_coef_bmi * age
            + params.sex_coef_bmi * sex
            + eps_bmi
        )
        return t2d, bmi

    t2d_gwas, bmi_gwas = _draw(params.seed_gwas_replicate)
    t2d_pgs, bmi_pgs = _draw(params.seed_pgs_replicate)

    # -- 7. write ----------------------------------------------------------
    fid = fam["FID"].values
    iid = fam["IID"].values

    pd.DataFrame({"FID": fid, "IID": iid, "age": age, "sex": sex}).to_csv(
        out_dir / "covariates.csv", index=False
    )
    pd.DataFrame({"FID": fid, "IID": iid, "T2D": t2d_gwas, "BMI": bmi_gwas}).to_csv(
        out_dir / "phenotypes_gwas.csv", index=False
    )
    pd.DataFrame({"FID": fid, "IID": iid, "T2D": t2d_pgs, "BMI": bmi_pgs}).to_csv(
        out_dir / "phenotypes_pgs_eval.csv", index=False
    )
    pd.DataFrame(
        {
            "FID": fid,
            "IID": iid,
            "pgs_T2D": g_std_t2d,
            "pgs_BMI": g_std_bmi,
            "pgs_T2D_raw": g_raw_t2d,
            "pgs_BMI_raw": g_raw_bmi,
        }
    ).to_csv(out_dir / "pgs_scores.csv", index=False)

    summary = {
        "n_samples": n,
        "n_snps_bim": len(bim),
        "n_pgs_t2d_total": n_t2d_total,
        "n_pgs_t2d_overlap": len(pgs_t2d),
        "n_pgs_bmi_total": n_bmi_total,
        "n_pgs_bmi_overlap": len(pgs_bmi),
        "h2_t2d": params.h2_t2d,
        "h2_bmi": params.h2_bmi,
        "t2d_prevalence_gwas": float(t2d_gwas.mean()),
        "t2d_prevalence_pgs": float(t2d_pgs.mean()),
        "bmi_mean_gwas": float(bmi_gwas.mean()),
        "bmi_std_gwas": float(bmi_gwas.std()),
        "bmi_mean_pgs": float(bmi_pgs.mean()),
        "bmi_std_pgs": float(bmi_pgs.std()),
        "g_std_t2d_mean": float(g_std_t2d.mean()),
        "g_std_t2d_std": float(g_std_t2d.std()),
        "g_std_bmi_mean": float(g_std_bmi.mean()),
        "g_std_bmi_std": float(g_std_bmi.std()),
    }
    pd.DataFrame([summary]).to_csv(out_dir / "simulation_summary.csv", index=False)

    pd.DataFrame(
        [
            {"Step": "Total SNPs (BIM)", "Count": len(bim)},
            {"Step": "PGS003443 (T2D) — catalog total", "Count": n_t2d_total},
            {"Step": "PGS003443 (T2D) — overlap w/ genotypes", "Count": len(pgs_t2d)},
            {"Step": "PGS004994 (BMI) — catalog total", "Count": n_bmi_total},
            {"Step": "PGS004994 (BMI) — overlap w/ genotypes", "Count": len(pgs_bmi)},
        ]
    ).to_csv(out_dir / "variant_counts.csv", index=False)

    log.info("step 1 complete: %s", out_dir)
    return PhenotypeOutputs({name: out_dir / name for name in OUTPUT_FILES})
