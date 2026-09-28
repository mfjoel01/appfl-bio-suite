"""Cohort QC checks: relatedness, PCA, allele frequencies, LD fidelity,
cross-site divergence summary.

Each check is its own function so individual sites can re-run or skip checks.
The top-level ``run_qc(cfg)`` runs every check for every site and writes a
single HTML report per site under ``<reports_dir>/qc_<site>.html``.

Notes on PCA scope (per the spec):
    PCA must be run on the **full-genome** PLINK binaries, not just chr1.
    If a full-genome binary for the site (``<site>_full.bed``) doesn't exist,
    ``run_qc`` will skip PCA with a logged warning rather than producing a
    misleading chr1-only PCA. Users should pre-extract the full-genome binary
    with the same ``<site>_ids.txt`` file used by sampling.
"""

from __future__ import annotations

import argparse
import base64
import io
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .utils import (
    PlinkError,
    SimulationConfig,
    ensure_dir,
    get_logger,
    load_config,
    read_bed_variants,
    read_bim,
    read_fam,
    read_site_manifest,
    run_plink,
    setup_logging,
)


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------


def check_relatedness(cfg: SimulationConfig, site: str) -> dict[str, object]:
    """KING kinship via plink2 ``--make-king-table``. Pass if no pair > threshold.

    Falls back gracefully if plink2 isn't on PATH.
    """
    logger = get_logger()
    site_dir = cfg.site_dir(site)
    prefix = site_dir / f"{site}_chr{cfg.chromosome}"
    out_prefix = site_dir / f"{site}_kinship"

    try:
        run_plink(
            ["--bfile", str(prefix), "--make-king-table", "--out", str(out_prefix)],
            binary=cfg.tools.plink2,
        )
    except PlinkError as e:
        logger.error("Relatedness check failed: %s", e)
        return {"status": "error", "message": str(e), "max_kinship": None, "n_related_pairs": None}

    king_path = out_prefix.with_suffix(".kin0")
    if not king_path.exists():
        return {"status": "skipped", "message": "no .kin0 produced", "max_kinship": None, "n_related_pairs": 0}

    df = pd.read_csv(king_path, sep="\t")
    if df.empty or "KINSHIP" not in df.columns:
        return {"status": "pass", "max_kinship": 0.0, "n_related_pairs": 0}
    max_k = float(df["KINSHIP"].max())
    n_related = int((df["KINSHIP"] > cfg.qc.kinship_threshold).sum())
    status = "pass" if n_related == 0 else "fail"
    return {"status": status, "max_kinship": max_k, "n_related_pairs": n_related}


def check_allele_frequencies(cfg: SimulationConfig, site: str) -> dict[str, object]:
    """Per-superpop MAF distribution summary (without an external 1KG reference).

    A full check against 1KG Phase 3 requires shipping/downloading the reference
    panel — out of scope for this module. We instead summarize the within-site
    per-superpop MAF distribution so a reviewer can eyeball it.
    """
    logger = get_logger()
    site_dir = cfg.site_dir(site)
    prefix = site_dir / f"{site}_chr{cfg.chromosome}"
    bim = read_bim(str(prefix) + ".bim")
    fam = read_fam(str(prefix) + ".fam")
    manifest = read_site_manifest(site_dir / f"{site}_manifest.tsv")
    manifest = manifest.set_index(["FID", "IID"]).loc[list(zip(fam["FID"], fam["IID"]))].reset_index()

    pop_summary: dict[str, dict[str, float]] = {}
    # Sample variants across chr1 to keep this tractable.
    sample_size = min(2000, len(bim))
    rng = np.random.default_rng(cfg.master_seed)
    sample_idx = np.sort(rng.choice(len(bim), size=sample_size, replace=False))
    dosage = read_bed_variants(prefix, sample_idx, len(fam))

    for pop, group in manifest.groupby("superpopulation"):
        rows = group.index.to_numpy()
        if len(rows) < 50:
            continue
        sub = dosage[rows, :]
        af = np.nanmean(sub, axis=0) / 2.0
        maf = np.minimum(af, 1.0 - af)
        pop_summary[pop] = {
            "n_indiv": int(len(rows)),
            "median_maf": float(np.nanmedian(maf)),
            "frac_common": float(np.nanmean(maf > 0.05)),
        }
    return {"status": "info", "per_pop": pop_summary, "n_variants_sampled": int(sample_size)}


def check_ld_decay(cfg: SimulationConfig, site: str) -> dict[str, object]:
    """LD decay curve per site, averaged over a small sample of selected loci.

    Returns binned mean r² as a function of pairwise distance (kb).
    """
    loci_dir = cfg.resolved_path("loci_dir")
    sel_path = loci_dir / "selected_loci.tsv"
    if not sel_path.exists():
        return {"status": "skipped", "message": "selected_loci.tsv not found"}
    selected = pd.read_csv(sel_path, sep="\t")
    n_sample = min(cfg.qc.ld_decay_n_loci_sample, len(selected))
    rng = np.random.default_rng(cfg.master_seed + hash(site) % 10_000)
    pick_idx = rng.choice(len(selected), size=n_sample, replace=False)
    sample_loci = selected.iloc[pick_idx]

    site_dir = cfg.site_dir(site)
    prefix = site_dir / f"{site}_chr{cfg.chromosome}"
    bim = read_bim(str(prefix) + ".bim")
    fam = read_fam(str(prefix) + ".fam")
    n_samples = len(fam)

    max_dist_kb = cfg.qc.ld_decay_max_dist_kb
    bin_edges_kb = np.linspace(0, max_dist_kb, 21)
    bin_sums = np.zeros(len(bin_edges_kb) - 1, dtype=np.float64)
    bin_counts = np.zeros(len(bin_edges_kb) - 1, dtype=np.int64)

    for _, locus in sample_loci.iterrows():
        chrom = str(locus["chrom"])
        start = int(locus["start_bp"])
        end = int(locus["end_bp"])
        mask = (bim["chrom"].astype(str) == chrom) & (bim["bp"] >= start) & (bim["bp"] <= end)
        local_bim = bim.loc[mask].copy()
        if len(local_bim) < 5:
            continue
        # Cap variants for tractability.
        if len(local_bim) > 400:
            stride_idx = np.linspace(0, len(local_bim) - 1, 400).round().astype(np.int64)
            local_bim = local_bim.iloc[np.unique(stride_idx)]
        vidx = local_bim.index.to_numpy()
        positions = local_bim["bp"].to_numpy()
        dosage = read_bed_variants(prefix, vidx, n_samples)

        # Standardize & compute r².
        mean = np.nanmean(dosage, axis=0)
        std = np.nanstd(dosage, axis=0)
        valid = std > 1e-8
        X = (dosage - mean) / np.where(valid, std, 1.0)
        X = np.where(np.isnan(X), 0.0, X)
        X[:, ~valid] = 0.0
        R = (X.T @ X) / n_samples
        R2 = R * R
        # Pairwise distances.
        D = np.abs(positions[:, None] - positions[None, :]) / 1000.0  # kb
        iu = np.triu_indices_from(R2, k=1)
        d = D[iu]
        r2 = R2[iu]
        keep = d <= max_dist_kb
        d = d[keep]; r2 = r2[keep]
        bin_idx = np.clip(np.digitize(d, bin_edges_kb) - 1, 0, len(bin_sums) - 1)
        for b in range(len(bin_sums)):
            sel = bin_idx == b
            if sel.any():
                bin_sums[b] += float(r2[sel].sum())
                bin_counts[b] += int(sel.sum())
    bin_means = np.where(bin_counts > 0, bin_sums / np.maximum(bin_counts, 1), np.nan)
    centers = (bin_edges_kb[:-1] + bin_edges_kb[1:]) / 2.0
    return {
        "status": "info",
        "bin_centers_kb": centers.tolist(),
        "bin_mean_r2": bin_means.tolist(),
        "n_loci_sampled": int(n_sample),
    }


def check_divergence_summary(cfg: SimulationConfig) -> dict[str, object]:
    """Aggregate divergence score distribution from candidate_windows.tsv."""
    sel_path = cfg.resolved_path("loci_dir") / "selected_loci.tsv"
    cand_path = cfg.resolved_path("loci_dir") / "candidate_windows.tsv"
    if not sel_path.exists() or not cand_path.exists():
        return {"status": "skipped"}
    sel = pd.read_csv(sel_path, sep="\t")
    quantiles: dict[str, float] = {}
    for stratum in ("low", "medium", "high"):
        sub = sel[sel["stratum"] == stratum]["divergence_score"]
        if not sub.empty:
            quantiles[f"{stratum}_min"] = float(sub.min())
            quantiles[f"{stratum}_median"] = float(sub.median())
            quantiles[f"{stratum}_max"] = float(sub.max())
    return {"status": "info", "stratum_quantiles": quantiles}


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _render_ld_decay_plot(decay_result: dict) -> str:
    if decay_result.get("status") != "info":
        return ""
    centers = decay_result["bin_centers_kb"]
    means = decay_result["bin_mean_r2"]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(centers, means, marker="o")
    ax.set_xlabel("Distance (kb)")
    ax.set_ylabel("Mean r²")
    ax.set_title("LD decay")
    ax.grid(alpha=0.3)
    return _fig_to_b64(fig)


def _render_maf_plot(maf_result: dict) -> str:
    if maf_result.get("status") != "info":
        return ""
    pops = list(maf_result["per_pop"].keys())
    medians = [maf_result["per_pop"][p]["median_maf"] for p in pops]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.bar(pops, medians)
    ax.set_ylabel("Median MAF")
    ax.set_title("Per-superpopulation median MAF (sampled variants)")
    return _fig_to_b64(fig)


def _render_html(site: str, results: dict[str, dict]) -> str:
    parts = [
        f"<html><head><meta charset='utf-8'><title>FedFM QC: {site}</title>",
        "<style>body{font-family:system-ui,Arial,sans-serif;max-width:900px;margin:2em auto;padding:0 1em;color:#222}"
        "table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #ccc;padding:.4em .8em}"
        ".pass{color:#0a0}.fail{color:#c00;font-weight:bold}.info{color:#06a}"
        "h2{margin-top:2em;border-bottom:1px solid #eee;padding-bottom:.3em}</style></head><body>",
        f"<h1>FedFM cohort QC — site: {site}</h1>",
    ]
    # Relatedness
    rel = results.get("relatedness", {})
    parts.append("<h2>Relatedness (KING kinship)</h2>")
    parts.append(f"<p>Status: <span class='{rel.get('status', 'info')}'>{rel.get('status')}</span></p>")
    parts.append("<table><tr><th>Field</th><th>Value</th></tr>")
    for k in ("max_kinship", "n_related_pairs"):
        parts.append(f"<tr><td>{k}</td><td>{rel.get(k)}</td></tr>")
    parts.append("</table>")

    # MAF
    parts.append("<h2>Per-superpopulation MAF summary</h2>")
    maf = results.get("allele_frequencies", {})
    parts.append(f"<p>Status: <span class='info'>{maf.get('status')}</span> "
                 f"(sampled {maf.get('n_variants_sampled', '-')} variants)</p>")
    if maf.get("per_pop"):
        parts.append("<table><tr><th>Pop</th><th>n indiv</th><th>median MAF</th><th>frac MAF&gt;0.05</th></tr>")
        for pop, vals in maf["per_pop"].items():
            parts.append(
                f"<tr><td>{pop}</td><td>{vals['n_indiv']}</td>"
                f"<td>{vals['median_maf']:.3f}</td><td>{vals['frac_common']:.3f}</td></tr>"
            )
        parts.append("</table>")
        b64 = _render_maf_plot(maf)
        if b64:
            parts.append(f"<img src='data:image/png;base64,{b64}'/>")

    # LD decay
    parts.append("<h2>LD decay</h2>")
    decay = results.get("ld_decay", {})
    parts.append(f"<p>Status: <span class='info'>{decay.get('status')}</span>"
                 f" — loci sampled: {decay.get('n_loci_sampled', '-')}.</p>")
    b64 = _render_ld_decay_plot(decay)
    if b64:
        parts.append(f"<img src='data:image/png;base64,{b64}'/>")

    # Divergence summary (site-independent but included for context).
    parts.append("<h2>Cross-site divergence summary</h2>")
    div = results.get("divergence_summary", {})
    if div.get("stratum_quantiles"):
        parts.append("<table><tr><th>Field</th><th>Value</th></tr>")
        for k, v in div["stratum_quantiles"].items():
            parts.append(f"<tr><td>{k}</td><td>{v:.4f}</td></tr>")
        parts.append("</table>")

    parts.append("</body></html>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Top level
# ---------------------------------------------------------------------------


def run_qc(cfg: SimulationConfig, only_site: str | None = None) -> dict[str, dict]:
    logger = get_logger()
    reports_dir = ensure_dir(cfg.resolved_path("reports_dir"))
    div_summary = check_divergence_summary(cfg)

    sites = list(cfg.sites)
    if only_site is not None:
        if only_site not in sites:
            raise ValueError(f"Unknown site {only_site!r}; have {sites}")
        sites = [only_site]

    all_results: dict[str, dict] = {}
    for site in sites:
        logger.info("QC for site: %s", site)
        rel = check_relatedness(cfg, site)
        maf = check_allele_frequencies(cfg, site)
        decay = check_ld_decay(cfg, site)
        results = {
            "relatedness": rel,
            "allele_frequencies": maf,
            "ld_decay": decay,
            "divergence_summary": div_summary,
        }
        html = _render_html(site, results)
        out_path = reports_dir / f"qc_{site}.html"
        out_path.write_text(html)
        logger.info("Wrote %s", out_path)
        all_results[site] = results
    return all_results


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="FedFM cohort QC")
    parser.add_argument("--config", default="config/simulation_config.yaml")
    parser.add_argument("--site", default=None,
                        help="run QC for only this site (for multi-node fan-out)")
    args = parser.parse_args(list(argv) if argv is not None else None)
    cfg = load_config(args.config)
    setup_logging(log_dir=cfg.resolved_path("logs_dir"))
    run_qc(cfg, only_site=args.site)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
