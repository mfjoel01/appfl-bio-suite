"""Scientific acceptance checks and locus-clustered outcome summaries."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

KEY = ["locus_id", "architecture_id", "replicate"]


def clustered_ratio(
    frame: pd.DataFrame, numerator: str, denominator: str, seed: int = 20260914, draws: int = 2000
) -> tuple[float, float, float]:
    """Bootstrap whole loci, preserving all replicates/sets within each draw."""
    if "locus_id" not in frame:
        raise ValueError("Locus identifiers are required for scientific uncertainty intervals")
    totals = frame.groupby("locus_id")[[numerator, denominator]].sum()
    n = len(totals)
    den = totals[denominator].sum()
    point = float(totals[numerator].sum() / den) if den else np.nan
    if n < 2 or not den:
        return point, np.nan, np.nan
    rng = np.random.default_rng(seed)
    picks = rng.integers(n, size=(draws, n))
    a = totals[numerator].to_numpy()[picks].sum(axis=1)
    b = totals[denominator].to_numpy()[picks].sum(axis=1)
    ratios = np.divide(a, b, out=np.full(draws, np.nan), where=b != 0)
    lo, hi = np.nanquantile(ratios, [0.025, 0.975])
    return point, min(float(lo), point), max(float(hi), point)


def scientific_summary(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["architecture_id", "stratum"]
    if "causal_mode" in results:
        group_cols.append("causal_mode")
    for keys, group in results.groupby(group_cols, dropna=False):
        group = group.copy()
        group["opportunities"] = 1
        # Keep algorithm failure in unconditional power; probabilities themselves stay NA.
        group["hits"] = group.any_causal_captured.astype(int)
        measures = {
            "power_any": ("hits", "opportunities"),
            "causal_recall": ("n_causal_captured", "n_causal"),
            "cs_coverage": ("n_cs_containing_causal", "n_credible_sets"),
            "causal_eligibility": ("n_eligible_causal", "n_causal"),
        }
        group["eligibility_denominator"] = group.n_causal.where(group.n_eligible_causal.notna(), 0)
        measures["causal_eligibility"] = ("n_eligible_causal", "eligibility_denominator")
        group["discoveries50"] = group.n_true_pip50 + group.n_false_pip50
        measures["fdr_pip50"] = ("n_false_pip50", "discoveries50")
        for metric, (num, den) in measures.items():
            point, lo, hi = clustered_ratio(group, num, den)
            rows.append(
                {
                    **dict(zip(group_cols, keys, strict=True)),
                    "metric": metric,
                    "estimate": point,
                    "ci_low": lo,
                    "ci_high": hi,
                    "n_loci": group.locus_id.nunique(),
                    "n_instances": len(group),
                    "n_nonconverged": group.fit_status.eq("nonconverged").sum(),
                    "uncertainty": "locus_cluster_bootstrap",
                }
            )
    return pd.DataFrame(rows)


def _artifact_table(root: Path, inst: str, suffix: str) -> pd.DataFrame:
    path = root / "artifacts" / inst / f"cs.{suffix}.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    try:
        return pd.read_csv(path, sep="\t", comment="#")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def parity_check(central_root: Path, federated_root: Path, destination: Path) -> dict:
    """Compare every result key, CS membership and emitted component probability.

    Existing print-resolution tolerances are retained; inputs may not silently
    differ. No missing or failed instance can make a comparison vacuously pass.
    """
    c = pd.read_csv(central_root / "fm_results.tsv", sep="\t")
    f = pd.read_csv(federated_root / "fed_fm_results.tsv", sep="\t")
    if c.duplicated(KEY).any() or f.duplicated(KEY).any():
        raise ValueError("Duplicate parity keys")
    ci = c.set_index(KEY)
    fi = f.set_index(KEY)
    if set(ci.index) != set(fi.index) or ci.empty:
        raise ValueError("Parity requires identical nonempty instance grids")
    mismatches = []
    max_delta = 0.0
    for key in ci.index:
        cr, fr = ci.loc[key], fi.loc[key]
        inst = f"{key[0]}_{key[1]}_rep{key[2]}"
        reasons = []
        execution = []
        for directory in [central_root, federated_root]:
            with gzip.open(directory / "artifacts" / inst / "cs.execution.json.gz", "rt") as handle:
                execution.append(json.load(handle))
        for field in ["variant_inputs", "options", "sample_sizes", "protocol"]:
            if field not in execution[0] or field not in execution[1]:
                reasons.append(f"missing_{field}")
            elif execution[0][field] != execution[1][field]:
                reasons.append(f"input_{field}")
        for col in ["fit_status", "n_causal_captured", "n_eligible_causal", "n_credible_sets"]:
            if not (pd.isna(cr[col]) and pd.isna(fr[col])) and cr[col] != fr[col]:
                reasons.append(col)
        if pd.notna(cr.error) or pd.notna(fr.error):
            reasons.append("execution_error")
        for suffix in ["cs", "snp"]:
            a = _artifact_table(central_root, inst, suffix)
            b = _artifact_table(federated_root, inst, suffix)
            if suffix == "cs":

                def sets(df):
                    return (
                        sorted(tuple(sorted(g.SNP.astype(str))) for _, g in df.groupby("CS_ID"))
                        if "CS_ID" in df
                        else []
                    )

                if sets(a) != sets(b):
                    reasons.append("cs_membership")
            else:
                if cr.fit_status == fr.fit_status == "nonconverged" and "FAIL" in a and "FAIL" in b:
                    continue
                if "SNP" not in a or "SNP" not in b:
                    reasons.append("missing_snp_output")
                    continue
                a = a.set_index("SNP").sort_index()
                b = b.set_index("SNP").sort_index()
                if not a.index.equals(b.index):
                    reasons.append("variant_eligibility")
                    continue
                # Compare overall probabilities: component numbering can permute.
                from .inference import inclusion_probabilities

                pa = inclusion_probabilities(a).to_numpy()
                pb = inclusion_probabilities(b).to_numpy()
                if len(pa):
                    max_delta = max(max_delta, float(np.max(np.abs(pa - pb))))
                if not np.allclose(pa, pb, rtol=2e-5, atol=1e-12):
                    reasons.append("pip")
        if reasons:
            mismatches.append({"instance": inst, "reasons": sorted(set(reasons))})
    result = {
        "n_instances": len(ci),
        "n_loci": c.locus_id.nunique(),
        "n_converged_pairs": int(
            (
                ci.fit_status.str.startswith("converged_")
                & fi.reindex(ci.index).fit_status.str.startswith("converged_")
            ).sum()
        ),
        "n_nonconverged_pairs": int(
            (
                ci.fit_status.eq("nonconverged")
                & fi.reindex(ci.index).fit_status.eq("nonconverged")
            ).sum()
        ),
        "n_mismatches": len(mismatches),
        "max_absolute_pip_difference": max_delta,
        "mismatches": mismatches,
        "pass": not mismatches,
        "rtol": 2e-5,
        "atol": 1e-12,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2))
    if mismatches:
        raise RuntimeError(
            f"Production parity failed for {len(mismatches)} instances: {destination}"
        )
    return result


def harvest_archives(root: Path, results: pd.DataFrame, destination: Path) -> pd.DataFrame:
    """Retained credible-set membership, including instance/locus identifiers."""
    rows = []
    for record in results.itertuples(index=False):
        inst = f"{record.locus_id}_{record.architecture_id}_rep{record.replicate}"
        metadata = json.loads((root / "artifacts" / inst / "record.json").read_text())
        truth = set(metadata["truth"])
        table = _artifact_table(root, inst, "cs")
        if "CS_ID" not in table:
            continue
        for cid, group in table.groupby("CS_ID"):
            members = set(group.SNP.astype(str))
            rows.append(
                {
                    "instance": inst,
                    "locus_id": record.locus_id,
                    "architecture_id": record.architecture_id,
                    "replicate": record.replicate,
                    "stratum": record.stratum,
                    "cs_id": cid,
                    "cs_length": len(members),
                    "contains_causal": bool(members & truth),
                    "members_json": json.dumps(sorted(members)),
                }
            )
    output = pd.DataFrame(
        rows,
        columns=[
            "instance",
            "locus_id",
            "architecture_id",
            "replicate",
            "stratum",
            "cs_id",
            "cs_length",
            "contains_causal",
            "members_json",
        ],
    )
    output.to_csv(destination, sep="\t", index=False)
    return output


def paired_arm_differences(frame: pd.DataFrame, reference: str = "federation") -> pd.DataFrame:
    """Paired locus bootstrap of power differences on identical simulation instances."""
    base = frame[frame.arm == reference].set_index(KEY)
    rows = []
    for arm, group in frame.groupby("arm"):
        if arm == reference:
            continue
        other = group.set_index(KEY)
        if set(base.index) != set(other.index):
            raise ValueError(f"Arm {arm} is not paired with {reference}")
        other = other.loc[base.index]
        d = (
            (base.any_causal_captured.astype(int) - other.any_causal_captured.astype(int))
            .rename("difference")
            .reset_index()
        )
        d["one"] = 1
        point, lo, hi = clustered_ratio(d, "difference", "one")
        rows.append(
            {
                "reference": reference,
                "arm": arm,
                "power_difference": point,
                "ci_low": lo,
                "ci_high": hi,
                "n_loci": d.locus_id.nunique(),
                "n_instances": len(d),
                "uncertainty": "paired_locus_cluster_bootstrap",
            }
        )
    return pd.DataFrame(rows)
