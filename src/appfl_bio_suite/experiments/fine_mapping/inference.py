"""Auditable SuSiEx outcomes and retained-component inclusion probabilities."""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

PROTOCOL = "complete-variants-pooled-maf-v3"
DEFAULT_OPTIONS = {"keep_ambiguous": True, "n_signals": 10, "max_iter": 1000, "tol": 1e-6}
METRIC_COLS = (
    "n_causal",
    "n_credible_sets",
    "cs_sizes_json",
    "total_cs_snps",
    "n_causal_captured",
    "any_causal_captured",
    "best_cs_size",
    "causal_pip_max",
    "causal_pip_mean",
    "top_pip",
    "mean_cs_purity",
    "min_cs_purity",
    "converged",
    "fit_status",
    "n_eligible_causal",
    "n_excluded_causal",
    "causal_pip_mean_eligible",
    "n_cs_containing_causal",
    "n_true_pip50",
    "n_false_pip50",
    "n_true_pip95",
    "n_false_pip95",
    "pip_definition",
    "analysis_protocol",
)


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", comment="#")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def inclusion_probabilities(df: pd.DataFrame) -> pd.Series:
    """Overall PIP over retained CS components, matching SuSiEx OVRL_PIP.

    Filtered components are not emitted by SuSiEx; this is explicitly not a
    reconstruction of the full, unfiltered model posterior.
    """
    cols = [c for c in df if c.startswith("PIP(")]
    if not cols:
        return pd.Series(0.0, index=df.index)
    alpha = df[cols].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(alpha.to_numpy()).all() or ((alpha < 0) | (alpha > 1)).any().any():
        raise ValueError("SuSiEx emitted invalid component probabilities")
    return 1.0 - (1.0 - alpha).prod(axis=1)


def pip_by_snp(path: Path) -> dict[str, float]:
    df = read_table(path)
    if "SNP" not in df:
        return {}
    return dict(zip(df.SNP.astype(str), inclusion_probabilities(df), strict=True))


def parse_outputs(out_dir: Path, out_name: str, truth_snps: list[str]) -> dict:
    cs = read_table(out_dir / f"{out_name}.cs")
    summ = read_table(out_dir / f"{out_name}.summary")
    snps = read_table(out_dir / f"{out_name}.snp")
    truth = set(truth_snps)
    has_cs = not cs.empty and {"CS_ID", "SNP"}.issubset(cs)
    # The installed writer emits FAIL on nonconvergence and NULL on a successful
    # fit with all sets filtered. File absence is neither of these outcomes.
    if "FAIL" in cs.columns or "FAIL" in summ.columns:
        status = "nonconverged"
    elif has_cs:
        status = "converged_cs"
    elif "NULL" in summ.columns and "SNP" in snps.columns:
        status = "converged_no_cs"
    else:
        status = "invalid_output"
    converged = status.startswith("converged_")
    eligible = set(snps.SNP.astype(str)) if "SNP" in snps else set()
    pip = pip_by_snp(out_dir / f"{out_name}.snp") if converged else {}
    groups = [set(g.SNP.astype(str)) for _, g in cs.groupby("CS_ID")] if has_cs else []
    sizes = [len(g) for g in groups]
    captured = set().union(*groups) & truth if groups else set()
    successful = [len(g) for g in groups if g & truth]
    pur = (
        pd.to_numeric(summ["CS_PURITY"], errors="raise")
        if "CS_PURITY" in summ
        else pd.Series(dtype=float)
    )
    # Unconditional recovery assigns zero to excluded truths; eligibility-conditioned
    # PIP is a separate metric. A failed fit never receives posterior probabilities.
    causal_pips = [pip.get(s, 0.0) for s in sorted(truth)] if converged else []
    eligible_pips = [pip.get(s, 0.0) for s in sorted(truth & eligible)] if converged else []
    result = {
        "n_causal": len(truth),
        "n_credible_sets": len(groups),
        "cs_sizes_json": json.dumps(sorted(sizes)),
        "total_cs_snps": sum(sizes),
        "n_causal_captured": len(captured),
        "any_causal_captured": bool(captured),
        "best_cs_size": min(successful) if successful else np.nan,
        "causal_pip_max": max(causal_pips) if causal_pips else np.nan,
        "causal_pip_mean": float(np.mean(causal_pips)) if causal_pips else np.nan,
        "top_pip": max(pip.values()) if pip else (0.0 if converged else np.nan),
        "mean_cs_purity": float(pur.mean()),
        "min_cs_purity": float(pur.min()),
        "converged": converged,
        "fit_status": status,
        "n_eligible_causal": len(truth & eligible) if "SNP" in snps else np.nan,
        "n_excluded_causal": len(truth - eligible) if "SNP" in snps else np.nan,
        "causal_pip_mean_eligible": float(np.mean(eligible_pips)) if eligible_pips else np.nan,
        "n_cs_containing_causal": sum(bool(g & truth) for g in groups),
        "pip_definition": "overall_retained_components",
        "analysis_protocol": PROTOCOL,
    }
    for label, threshold in [("50", 0.5), ("95", 0.95)]:
        discovered = {s for s, prob in pip.items() if prob > threshold}
        result[f"n_true_pip{label}"] = len(discovered & truth)
        result[f"n_false_pip{label}"] = len(discovered - truth)
    return result


def archive_instance(priv: Path, destination: Path, truth: list[str], row: dict) -> None:
    """Preserve fit evidence even when genotype/sumstat scratch is discarded."""
    destination.mkdir(parents=True, exist_ok=True)
    for source in [*priv.glob("cs.*"), *priv.glob("*.sumstats")]:
        if source.is_file():
            with (
                source.open("rb") as src,
                gzip.open(destination / (source.name + ".gz"), "wb") as dst,
            ):
                shutil.copyfileobj(src, dst)
    record = {"instance": priv.name, "truth": truth, "result": row}
    (destination / "record.json").write_text(json.dumps(record, indent=2, default=str))


def input_fingerprint(paths: list[Path], options: dict) -> str:
    """Content fingerprint for small metadata and file identity for genotype payloads."""
    digest = hashlib.sha256(json.dumps(options, sort_keys=True).encode())
    for path in paths:
        path = path.resolve()
        stat = path.stat()
        digest.update(str((str(path), stat.st_size, stat.st_mtime_ns)).encode())
        if path.suffix != ".bed":
            digest.update(path.read_bytes())
    return digest.hexdigest()


def validate_result_grid(df: pd.DataFrame, expected: set[tuple]) -> None:
    keys = ["locus_id", "architecture_id", "replicate"]
    if df.empty or df.duplicated(keys).any():
        raise ValueError("Missing result table or duplicate analysis instances")
    actual = set(df[keys].itertuples(index=False, name=None))
    if actual != expected:
        raise ValueError(
            f"Incomplete result grid: missing={len(expected - actual)}, "
            f"extra={len(actual - expected)}"
        )
    if df.error.fillna("").astype(str).ne("").any():
        raise RuntimeError("Execution failures in fine-mapping results; inspect retained artifacts")
    if df.fit_status.eq("invalid_output").any():
        raise RuntimeError("Invalid SuSiEx output; inspect retained artifacts")


def variant_input_manifest(
    sumstats: list[Path],
    ld_prefixes: list[Path],
    sample_sizes: list[int],
    require_matching: bool = True,
) -> list[dict]:
    """Fingerprint both ordered, allele-coded variant lists for each ancestry.

    The borrowed external-LD sensitivity arm explicitly permits different lists;
    its manifests are still recorded and it is excluded from the equivalence claim.
    """
    if not sumstats or not (len(sumstats) == len(ld_prefixes) == len(sample_sizes)):
        raise ValueError("SuSiEx population columns and sample sizes do not align")
    records = []
    for ss, prefix, n in zip(sumstats, ld_prefixes, sample_sizes, strict=True):
        s = pd.read_csv(ss, sep="\t", dtype=str).iloc[:, :5]
        s.columns = ["chr", "snp", "bp", "a1", "a2"]
        r = pd.read_csv(f"{prefix}_ref.bim", sep=r"\s+", header=None, dtype=str)
        r = r.iloc[:, [0, 1, 3, 4, 5]]
        r.columns = s.columns
        for table in [s, r]:
            if table.empty or table.snp.duplicated().any() or table.isna().any().any():
                raise ValueError(f"Invalid or duplicated input variants for {ss.stem}")
            table["bp"] = pd.to_numeric(table.bp, errors="raise").astype(np.int64).astype(str)
        if require_matching and not s.reset_index(drop=True).equals(r.reset_index(drop=True)):
            raise ValueError(
                f"{ss.stem}: sumstats/LD ordered variants or alleles differ "
                f"(sumstats={len(s)}, LD={len(r)}); refusing inference"
            )
        records.append(
            {
                "population": ss.stem,
                "sample_size": int(n),
                "sumstats_n_variants": len(s),
                "ld_n_variants": len(r),
                "sumstats_variants_sha256": hashlib.sha256(
                    s.to_csv(index=False, header=False).encode()
                ).hexdigest(),
                "ld_variants_sha256": hashlib.sha256(
                    r.to_csv(index=False, header=False).encode()
                ).hexdigest(),
                "require_matching": require_matching,
            }
        )
    return records


def require_successful_fits(frame: pd.DataFrame) -> None:
    """Preserve failed rows for diagnosis, but never report a successful batch."""
    if frame.empty:
        return
    if (
        frame.error.fillna("").astype(str).ne("").any()
        or frame.fit_status.eq("invalid_output").any()
    ):
        raise RuntimeError("Fine-mapping batch failed; inspect result rows and archived inputs")


def harmonize_sumstats(frame: pd.DataFrame, reference_bim: Path) -> pd.DataFrame:
    """Align GWAS effect alleles to the genotype reference, including signed statistics."""
    ref = pd.read_csv(
        reference_bim,
        sep=r"\s+",
        header=None,
        names=["chr", "snp", "cm", "bp", "a1", "a2"],
        dtype={"snp": str},
    )
    if frame.snp.duplicated().any() or ref.snp.duplicated().any():
        raise ValueError("Duplicate variants cannot be harmonized")
    ref = ref.set_index("snp").loc[frame.snp]
    out = frame.copy()
    a1, a2 = ref.a1.to_numpy(), ref.a2.to_numpy()
    same = (out.A1.to_numpy() == a1) & (out.A2.to_numpy() == a2)
    swapped = (out.A1.to_numpy() == a2) & (out.A2.to_numpy() == a1)
    if not (same | swapped).all():
        raise ValueError("GWAS alleles incompatible with the genotype reference")
    if not np.array_equal(pd.to_numeric(out.bp), pd.to_numeric(ref.bp)):
        raise ValueError("GWAS positions differ from the genotype reference")
    out.loc[swapped, ["beta", "stat"]] *= -1
    out["A1"], out["A2"] = a1, a2
    return out
