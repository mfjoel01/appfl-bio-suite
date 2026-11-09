"""Coordinator-only MultiSuSiE RSS adapter over pooled, allele-aligned moments."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .inference import DEFAULT_OPTIONS, summarize_fit, variant_input_manifest

PROTOCOL = "multisusie-common-variants-pooled-maf-v1"


def require_multisusie():
    try:
        from MultiSuSiE import multisusie_rss
    except ImportError as exc:
        raise ImportError(
            "MultiSuSiE backend unavailable. Install the coordinator extra: "
            "pip install '.[finemapping-multisusie]' -c constraints.txt"
        ) from exc
    return multisusie_rss


def validate_options(options):
    rho = float(options.get("multisusie_rho", 0.75))
    purity = float(options.get("multisusie_min_abs_corr", 0.5))
    if not 0 <= rho < 1:
        raise ValueError("multisusie_rho must satisfy 0 <= rho < 1")
    if not 0 <= purity <= 1:
        raise ValueError("multisusie_min_abs_corr must be between 0 and 1")
    for key in ("n_signals", "max_iter"):
        value = options.get(key, DEFAULT_OPTIONS[key])
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError(f"{key} must be a positive integer")
    if (
        not np.isfinite(options.get("tol", DEFAULT_OPTIONS["tol"]))
        or options.get("tol", DEFAULT_OPTIONS["tol"]) <= 0
    ):
        raise ValueError("tol must be finite and positive")


def aligned_inputs(sumstats, ld_prefixes, sample_sizes, keep_ambiguous=True):
    """Intersect ancestry panels, explicitly reorder both LD axes, and reject allele drift.

    Within each ancestry the LD and summary-statistic lists must already match. Across
    ancestries, pooled MAF filtering can leave different lists. Only their intersection
    is fitted; excluded variants are recorded and counted as ineligible for scoring.
    """
    records = variant_input_manifest(sumstats, ld_prefixes, sample_sizes)
    frames = [pd.read_csv(p, sep="\t", dtype={"snp": str, "chr": str}) for p in sumstats]
    common = set(frames[0].snp)
    for frame in frames[1:]:
        common.intersection_update(frame.snp)
    ids = [s for s in frames[0].snp if s in common]
    if not keep_ambiguous:
        ref = frames[0].set_index("snp").loc[ids]
        ids = [
            s
            for s, a, b in zip(ids, ref.A1, ref.A2, strict=True)
            if frozenset((a.upper(), b.upper())) not in (frozenset("AT"), frozenset("CG"))
        ]
    if not ids:
        raise ValueError("MultiSuSiE has no common eligible variants across ancestries")
    reference = frames[0].set_index("snp").loc[ids, ["chr", "bp", "A1", "A2"]]
    fitted = set(ids)
    aligned, matrices = [], []
    for frame, prefix, record in zip(frames, ld_prefixes, records, strict=True):
        index = pd.Index(frame.snp).get_indexer(ids)
        selected = frame.iloc[index].reset_index(drop=True)
        if not selected.set_index("snp")[["chr", "bp", "A1", "A2"]].equals(reference):
            raise ValueError("MultiSuSiE ancestry alleles or positions disagree; harmonize first")
        matrix = np.fromfile(f"{prefix}.ld.bin", dtype=np.float32)
        if matrix.size != len(frame) ** 2:
            raise ValueError("MultiSuSiE LD matrix size does not match its variant list")
        matrix = matrix.reshape(len(frame), len(frame))[np.ix_(index, index)].astype(np.float64)
        if (
            not np.isfinite(matrix).all()
            or not np.allclose(matrix, matrix.T, atol=1e-6)
            or not np.allclose(np.diag(matrix), 1, atol=1e-6)
        ):
            raise ValueError(
                "MultiSuSiE requires finite symmetric correlation matrices with unit diagonal"
            )
        values = selected[["beta", "se", "p"]].to_numpy(dtype=float)
        if (
            not np.isfinite(values).all()
            or (values[:, 1] <= 0).any()
            or (values[:, 2] < 0).any()
            or (values[:, 2] > 1).any()
        ):
            raise ValueError("MultiSuSiE requires finite effects, positive SEs and valid p-values")
        record["excluded_variants"] = [s for s in frame.snp if s not in fitted]
        aligned.append(selected)
        matrices.append(matrix)
    return ids, aligned, matrices, records


def run_multisusie(
    sumstats,
    ld_prefixes,
    sample_sizes,
    var_y,
    out_dir: Path,
    truth_snps,
    level,
    pval_thresh,
    options,
):
    """Run the upstream model and retain full-model PIPs separately from filtered CSs."""
    fit_function = require_multisusie()
    opts = {**DEFAULT_OPTIONS, "multisusie_rho": 0.75, "multisusie_min_abs_corr": 0.5, **options}
    validate_options(opts)
    if not 0 < level < 1 or not 0 < pval_thresh <= 1:
        raise ValueError("credible-set level must be in (0,1), pval_thresh in (0,1]")
    if len(var_y) != len(sample_sizes) or any(n < 3 for n in sample_sizes):
        raise ValueError("MultiSuSiE requires one sample variance per ancestry and n >= 3")
    if not np.isfinite(var_y).all() or (np.asarray(var_y) <= 0).any():
        raise ValueError("MultiSuSiE requires positive finite phenotype sample variances")
    ids, frames, matrices, records = aligned_inputs(
        sumstats, ld_prefixes, sample_sizes, opts["keep_ambiguous"]
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    execution = {
        "backend": "multisusie",
        "version": importlib.metadata.version("MultiSuSiE"),
        "protocol": PROTOCOL,
        "options": opts,
        "coverage": level,
        "pval_thresh": pval_thresh,
        "sample_sizes": sample_sizes,
        "phenotype_sample_variances": var_y,
        "variant_inputs": records,
        "fitted_variants": ids,
        "float_type": "float64",
        "input_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sumstats},
    }
    (out_dir / "cs.execution.json").write_text(json.dumps(execution, indent=2))
    fit = fit_function(
        b_list=[f.beta.to_numpy(dtype=float) for f in frames],
        s_list=[f.se.to_numpy(dtype=float) for f in frames],
        R_list=matrices,
        varY_list=var_y,
        population_sizes=sample_sizes,
        rho=opts["multisusie_rho"],
        L=opts["n_signals"],
        max_iter=opts["max_iter"],
        tol=opts["tol"],
        coverage=level,
        min_abs_corr=opts["multisusie_min_abs_corr"],
        float_type=np.float64,
        low_memory_mode=False,
        variant_ids=ids,
    )
    converged = bool(fit.converged)
    probabilities = np.asarray(fit.pip, dtype=float)
    if (
        probabilities.shape != (len(ids),)
        or not np.isfinite(probabilities).all()
        or (probabilities < 0).any()
        or (probabilities > 1).any()
    ):
        raise ValueError("MultiSuSiE returned invalid posterior inclusion probabilities")
    groups, purities, cs_rows = [], [], []
    marginal_p = np.min([f.p.to_numpy() for f in frames], axis=0)
    if converged:
        for component, included in enumerate(fit.sets[3]):
            if not included:
                continue
            indices = np.asarray(fit.sets[0][component], dtype=int)
            # Same marginal-p CS filter as SuSiEx; this never changes full-model PIPs.
            if not len(indices) or marginal_p[indices].min() >= pval_thresh:
                continue
            groups.append({ids[i] for i in indices})
            purity = float(fit.sets[1][component])
            purities.append(purity)
            cs_rows.extend(
                {"CS_ID": component + 1, "SNP": ids[i], "CS_PURITY": purity} for i in indices
            )
    status = ("converged_cs" if groups else "converged_no_cs") if converged else "nonconverged"
    # Preserve the native posterior even for nonconvergence, but do not score it as a fit.
    np.savez_compressed(
        out_dir / "cs.posterior.npz", alpha=fit.alpha, pip=probabilities, converged=converged
    )
    pd.DataFrame(
        {
            "SNP": ids,
            "PIP": probabilities if converged else np.nan,
            "PIP_DEFINITION": "multisusie_full_model",
        }
    ).to_csv(out_dir / "cs.snp", sep="\t", index=False)
    pd.DataFrame(cs_rows, columns=["CS_ID", "SNP", "CS_PURITY"]).to_csv(
        out_dir / "cs.cs", sep="\t", index=False
    )
    result = summarize_fit(
        truth_snps,
        ids,
        dict(zip(ids, probabilities, strict=True)),
        groups,
        purities,
        status,
        pip_definition="multisusie_full_model",
        protocol=PROTOCOL,
    )
    result["fine_mapping_backend"] = "multisusie"
    (out_dir / "cs.summary.json").write_text(json.dumps(result, indent=2))
    return result
