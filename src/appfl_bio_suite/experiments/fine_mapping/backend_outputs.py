"""Shared APPFL result schema and portable posterior tables for both estimators."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .inference import METRIC_COLS, inclusion_probabilities

TRUTH_METRICS = {
    "n_causal",
    "n_eligible_causal",
    "n_excluded_causal",
    "n_causal_captured",
    "any_causal_captured",
    "n_cs_containing_causal",
    "best_cs_size",
    "causal_pip_max",
    "causal_pip_mean",
    "causal_pip_mean_eligible",
    "n_true_pip50",
    "n_false_pip50",
    "n_true_pip95",
    "n_false_pip95",
}


def result_table(rows, backend, populations, truth_map):
    """Use identical columns; unknown truth is never a negative finding."""
    columns = [
        "locus_id",
        "architecture_id",
        "replicate",
        "stratum",
        *METRIC_COLS,
        *[f"min_p_{pop}" for pop in sorted(populations)],
        "runtime_s",
        "error",
        "backend",
    ]
    result = pd.DataFrame(rows).reindex(columns=columns)
    result["backend"] = backend
    unknown = [
        idx
        for idx, row in result.iterrows()
        if (row.locus_id, row.architecture_id, row.replicate) not in truth_map
    ]
    if unknown:
        # Missing truth must work on pandas versions that reject NaN in bool/int columns.
        for name in TRUTH_METRICS:
            result[name] = result[name].astype(float)
        result.loc[unknown, sorted(TRUTH_METRICS)] = np.nan
    return result


def export_susiex(archive: Path, destination: Path, columns, row, options):
    """Expose native SuSiEx retained-component PIPs in the shared table layout.

    The compressed native files remain in the existing artifacts directory. Failed
    fits get diagnostic metadata and missing PIPs, never apparent posterior evidence.
    """
    destination.mkdir(parents=True, exist_ok=False)

    def read(name):
        path = archive / name
        if not path.is_file():
            return pd.DataFrame()
        try:
            return pd.read_csv(path, sep="\t", comment="#")
        except pd.errors.EmptyDataError:
            return pd.DataFrame()

    snps, cs = read("cs.snp.gz"), read("cs.cs.gz")
    metadata = ["chrom", "snp_id", "bp", "a1", "a2"]
    variants = pd.concat([c.variants[metadata] for c in columns], ignore_index=True)
    variants = variants.drop_duplicates().sort_values(["chrom", "bp", "snp_id"])
    converged = bool(row.get("converged", False))
    pips = (
        dict(zip(snps.SNP.astype(str), inclusion_probabilities(snps), strict=True))
        if converged and "SNP" in snps
        else {}
    )
    variants["pip"] = variants.snp_id.map(pips)
    for col in sorted(columns, key=lambda c: c.pop):
        variants[f"eligible_{col.pop}"] = variants.snp_id.isin(col.variants.snp_id)
    variants.to_csv(destination / "variants.tsv", sep="\t", index=False)
    members = []
    if converged and {"CS_ID", "SNP", "CS_PIP"}.issubset(cs.columns):
        for component, group in cs.groupby("CS_ID"):
            coverage = float(group.CS_PIP.sum())
            for entry in group.itertuples(index=False):
                members.append(
                    dict(
                        component=int(component),
                        snp_id=entry.SNP,
                        component_pip=entry.CS_PIP,
                        pip=pips.get(entry.SNP, np.nan),
                        coverage=coverage,
                    )
                )
    pd.DataFrame(
        members, columns=["component", "snp_id", "component_pip", "pip", "coverage"]
    ).to_csv(destination / "credible_sets.tsv", sep="\t", index=False)
    (destination / "fit.json").write_text(
        json.dumps(
            {
                "method": "SuSiEx",
                "options": options,
                "converged": converged,
                "fit_status": row.get("fit_status"),
                "pip_definition": row.get("pip_definition"),
                "runtime_s": row.get("runtime_s"),
                "n_credible_sets": int(cs.CS_ID.nunique()) if members else 0,
                "native_artifacts": str(archive),
                "error": row.get("error", ""),
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
