"""Experimental coordinator-only MultiSuSiE fit from the existing raw moments.

Center after pooling within ancestry. Preserve the union of eligible variants across
ancestries, censoring an unavailable variant only in the ancestry that excludes it.
No truth, simulated rg, individual records, or rounded GWAS files enter inference.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import inspect
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .fedfm.fed_fine_mapping import GenoAggregate, PhenoAggregate, pool_geno, pool_pheno


@dataclass(frozen=True)
class MultiSuSiEOptions:
    rho: float = 0.75
    n_signals: int = 10
    max_iter: int = 1000
    tol: float = 1e-6
    coverage: float = 0.95
    min_abs_corr: float = 0.0
    maf: float = 0.005
    min_mac: float = 0.0  # Match the existing pooled-MAF protocol, not upstream's MAC 20.

    def __post_init__(self):
        if not np.isfinite(self.rho) or not 0 <= self.rho < 1:
            raise ValueError("rho must be in [0, 1); rho=1 is singular")
        if self.n_signals < 1 or self.max_iter < 1 or not np.isfinite(self.tol) or self.tol <= 0:
            raise ValueError("Signal count, iteration limit and tolerance must be positive")
        if not 0 < self.coverage < 1 or not 0 <= self.min_abs_corr <= 1:
            raise ValueError("Invalid coverage or purity threshold")
        if not 0 <= self.maf < 0.5 or not np.isfinite(self.min_mac) or self.min_mac < 0:
            raise ValueError("Invalid MAF or MAC threshold")


@dataclass
class MomentInputs:
    variants: pd.DataFrame
    populations: list[str]
    sample_sizes: list[int]
    xtx: list[np.ndarray]
    xty: list[np.ndarray]
    yty: list[float]
    eligible: list[np.ndarray]


def prepare_inputs(
    geno: list[GenoAggregate],
    pheno: list[PhenoAggregate],
    options: MultiSuSiEOptions,
) -> MomentInputs:
    """Validate matched site blocks, then pool, center and censor by ancestry."""
    if not geno or not pheno:
        raise ValueError("Both genotype and phenotype blocks are required")
    if len({g.locus_id for g in geno}) != 1 or len({p.instance for p in pheno}) != 1:
        raise ValueError("Fit exactly one locus and phenotype instance at a time")
    gkeys = [(g.site, g.pop) for g in geno]
    pkeys = [(p.site, p.pop) for p in pheno]
    if len(set(gkeys)) != len(gkeys) or len(set(pkeys)) != len(pkeys):
        raise ValueError("Duplicate site/ancestry blocks")
    if set(gkeys) != set(pkeys):
        raise ValueError("Genotype and phenotype site/ancestry membership differs")
    lookup = dict(zip(pkeys, pheno, strict=True))
    metadata = ["chrom", "snp_id", "bp", "a1", "a2"]
    for g in geno:
        p = lookup[g.site, g.pop]
        m = len(g.snp_ids)
        if g.n != p.n or g.n <= 0 or not np.array_equal(g.snp_ids, p.snp_ids):
            raise ValueError("Genotype and phenotype cohort or variant order differs")
        if g.G.shape != (m, m) or g.u.shape != (m,) or p.c.shape != (m,):
            raise ValueError("Moment dimensions do not match variant metadata")
        if g.variants[metadata].isna().any().any() or g.variants.snp_id.duplicated().any():
            raise ValueError("Missing or duplicate variant metadata")
        if any(not np.isfinite(a).all() for a in (g.G, g.u, p.c, [p.q, p.w])):
            raise ValueError("Nonfinite moments")
        if not np.allclose(g.G, g.G.T, rtol=0, atol=1e-8):
            raise ValueError("Genotype Gram matrix is not symmetric")
        if ((g.u < 0) | (g.u > 2 * g.n)).any():
            raise ValueError("Allele totals are outside dosage bounds")

    # Ancestry order must not depend on the order in which sites report.
    populations = sorted({g.pop for g in geno})
    pooled = []
    frames = []
    for pop in populations:
        g = pool_geno(sorted((g for g in geno if g.pop == pop), key=lambda g: g.site))
        p = pool_pheno(
            sorted((p for p in pheno if p.pop == pop), key=lambda p: p.site),
            g.variants.snp_id.to_numpy(),
        )
        if g.n <= 2:
            raise ValueError("Each pooled ancestry needs at least three samples")
        centered = g.G - np.outer(g.u, g.u) / g.n
        mac = np.minimum(g.u, 2 * g.n - g.u)
        keep = (mac >= max(options.min_mac, 2 * g.n * options.maf) - 1e-10) & (
            np.diag(centered) > 1e-10
        )
        yty = p.w - p.q**2 / p.n
        if not keep.any() or not np.isfinite(yty) or yty <= 0:
            raise ValueError(f"{pop}: no eligible variants or nonpositive phenotype variance")
        variants = g.variants.loc[keep, metadata].copy()
        variants["chrom"] = variants.chrom.astype(str)
        frames.append(variants)
        pooled.append((g.n, centered[np.ix_(keep, keep)], (p.c - g.u * p.q / g.n)[keep], yty))

    union = pd.concat(frames, ignore_index=True).drop_duplicates()
    if union.snp_id.duplicated().any():
        raise ValueError("Variant positions or allele coding differ across ancestries")
    union = union.sort_values(["chrom", "bp", "snp_id"]).reset_index(drop=True)
    index = pd.Index(union.snp_id)
    xtx, xty, eligible = [], [], []
    for frame, (_, gram, cross, _) in zip(frames, pooled, strict=True):
        idx = index.get_indexer(frame.snp_id)
        matrix = np.zeros((len(index), len(index)), dtype=np.float64)
        vector = np.zeros(len(index), dtype=np.float64)
        mask = np.zeros(len(index), dtype=bool)
        matrix[np.ix_(idx, idx)], vector[idx], mask[idx] = gram, cross, True
        xtx.append(matrix)
        xty.append(vector)
        eligible.append(mask)
    return MomentInputs(
        union, populations, [p[0] for p in pooled], xtx, xty, [p[3] for p in pooled], eligible
    )


def fit_multisusie(inputs: MomentInputs, options: MultiSuSiEOptions):
    """Fit without mutating caller-owned moments; import the optional backend lazily."""
    try:
        from MultiSuSiE.susiepy_ss import susie_multi_ss
    except ImportError as exc:
        raise ImportError(
            "Install the coordinator extra appfl-bio-suite[multisusie] and the backend "
            "from requirements-multisusie.txt"
        ) from exc
    k = len(inputs.populations)
    rho = np.full((k, k), options.rho, dtype=np.float64)
    np.fill_diagonal(rho, 1.0)
    return susie_multi_ss(
        XTX_list=[x.copy() for x in inputs.xtx],
        XTY_list=[x.copy() for x in inputs.xty],
        YTY_list=inputs.yty,
        population_sizes=inputs.sample_sizes,
        rho=rho,
        L=options.n_signals,
        standardize=False,
        estimate_residual_variance=True,
        estimate_prior_variance=True,
        estimate_prior_method="early_EM",
        pop_spec_effect_priors=True,
        max_iter=options.max_iter,
        tol=options.tol,
        coverage=options.coverage,
        min_abs_corr=options.min_abs_corr,
        float_type=np.float64,
        low_memory_mode=True,
        # The pinned backend returns a tuple when no CS survives, but attempts
        # list.append if variant_ids is supplied. We map integer CS indices to
        # our validated metadata below; optional backend labels are unnecessary.
        variant_ids=None,
    )


def run_multisusie(geno, pheno, destination: Path, options=None) -> dict:
    """Write native PIPs, component probabilities, CS membership and fit provenance.

    A new directory is required so an exploratory run cannot overwrite a previous fit.
    Nonconverged estimates are retained for diagnosis but are explicitly flagged.
    """
    options = options if options is not None else MultiSuSiEOptions()
    destination.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    record = {"method": "MultiSuSiE", "options": asdict(options), "fit_status": "running"}
    try:
        inputs = prepare_inputs(geno, pheno, options)
        record.update(populations=inputs.populations, sample_sizes=inputs.sample_sizes)
        fingerprint = hashlib.sha256(inputs.variants.to_csv(index=False).encode())
        for pop, n, x, c, w in zip(
            inputs.populations, inputs.sample_sizes, inputs.xtx, inputs.xty, inputs.yty, strict=True
        ):
            fingerprint.update(json.dumps([pop, n, w]).encode())
            fingerprint.update(x.tobytes())
            fingerprint.update(c.tobytes())
        record["moments_sha256"] = fingerprint.hexdigest()
        fit = fit_multisusie(inputs, options)
        from MultiSuSiE import susiepy_ss

        record["backend_version"] = importlib.metadata.version("MultiSuSiE")
        record["backend_source_sha256"] = hashlib.sha256(
            Path(inspect.getfile(susiepy_ss)).read_bytes()
        ).hexdigest()
        record["converged"] = bool(fit.converged)
        record["niter"] = int(fit.niter)
        record["pip_definition"] = "multisusie_native_active_components_before_cs_filtering"
        if not np.isfinite(fit.pip).all() or ((fit.pip < 0) | (fit.pip > 1)).any():
            raise ValueError("MultiSuSiE returned invalid inclusion probabilities")
        variants = inputs.variants.copy()
        variants["pip"] = fit.pip
        for pop, mask in zip(inputs.populations, inputs.eligible, strict=True):
            variants[f"eligible_{pop}"] = mask
        variants.to_csv(destination / "variants.tsv", sep="\t", index=False)
        members = []
        for component, keep in enumerate(fit.sets[3]):
            if keep:
                for idx in fit.sets[0][component]:
                    members.append(
                        {
                            "component": component + 1,
                            "snp_id": variants.iloc[idx].snp_id,
                            "component_pip": float(fit.alpha[component, idx]),
                            "pip": float(fit.pip[idx]),
                            "coverage": float(fit.sets[2][component]),
                        }
                    )
        pd.DataFrame(
            members, columns=["component", "snp_id", "component_pip", "pip", "coverage"]
        ).to_csv(destination / "credible_sets.tsv", sep="\t", index=False)
        np.savez_compressed(
            destination / "posterior.npz",
            alpha=fit.alpha,
            prior_variance=fit.V,
            coef=fit.coef,
            coef_sd=fit.coef_sd,
            residual_variance=fit.sigma2,
        )
        record["n_credible_sets"] = int(np.sum(fit.sets[3]))
        record["fit_status"] = (
            ("converged_cs" if members else "converged_no_cs") if fit.converged else "nonconverged"
        )
    except Exception as exc:
        record.update(fit_status="execution_failure", error=str(exc))
        raise
    finally:
        record["runtime_s"] = time.perf_counter() - started
        (destination / "fit.json").write_text(json.dumps(record, indent=2, allow_nan=False))
    return record
