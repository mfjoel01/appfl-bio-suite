"""Manhattan, QQ and hits-table rendering. COORDINATOR-SIDE ONLY.

Not shipped to workers, so it may import freely and may depend on matplotlib.

Plotting used to happen at each site, inside the shipped trainer. It was moved here
because the coordinator receives per-variant beta, standard error and MAF from every
site, so it can render identical per-site plots -- and doing it here means matplotlib is
not on any partner's dependency list, and the code exists once rather than being
duplicated between a shipped module and a server-side one that cannot import it.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

# Set before pyplot is imported. Chooses a non-interactive backend, without which
# importing pyplot on a headless node tries to reach a display and fails.
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

__all__ = ["plot_manhattan", "plot_qq", "write_hits_table", "normalize_chr"]

CHROM_MAP = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}

_BLUE = "#1f5a96"
_ORANGE = "#d76f30"
_RED = "#9b1c31"


def normalize_chr(chrom: pd.Series) -> pd.Series:
    return chrom.astype(str).str.strip().str.upper().replace(CHROM_MAP).astype(int)


def plot_manhattan(gwas_df: pd.DataFrame, trait: str, threshold: float, out_path, label: str = ""):
    """Manhattan plot: -log10(P) against genomic position, chromosomes alternating."""
    n_snps = len(gwas_df)
    n_col = "N_META" if "N_META" in gwas_df.columns else ("N" if "N" in gwas_df.columns else None)
    n_samples = int(gwas_df[n_col].iloc[0]) if n_col else 0

    plot_df = gwas_df[["CHR", "BP", "P"]].copy().sort_values(["CHR", "BP"]).reset_index(drop=True)

    # Lay chromosomes end to end with a gap, so position is continuous across the genome.
    chrom_offsets, tick_pos, tick_labels = {}, [], []
    offset = 0
    for chrom, group in plot_df.groupby("CHR", sort=True):
        chrom = int(chrom)
        chrom_offsets[chrom] = offset
        bp = group["BP"].to_numpy(dtype=np.int64)
        tick_pos.append(offset + 0.5 * (bp.min() + bp.max()))
        tick_labels.append(str(chrom))
        offset += int(bp.max()) + 1_000_000

    plot_df["X"] = plot_df["BP"] + plot_df["CHR"].map(chrom_offsets)
    # Clipped at the smallest positive float: a p-value of 0 would be infinite here.
    plot_df["LOGP"] = -np.log10(
        np.clip(plot_df["P"].to_numpy(dtype=np.float64), np.finfo(np.float64).tiny, 1.0)
    )

    fig, ax = plt.subplots(figsize=(14, 6))
    colors = [_BLUE, _ORANGE]
    for idx, (_, group) in enumerate(plot_df.groupby("CHR", sort=True)):
        ax.scatter(group["X"], group["LOGP"], s=4, color=colors[idx % 2], alpha=0.8, linewidths=0)
    ax.axhline(-np.log10(threshold), color=_RED, linestyle="--", linewidth=1.2)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=9)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(P)")
    ax.set_title(
        f"{trait} Manhattan Plot – {label}\nSNPs = {n_snps:,}  |  N = {n_samples:,}", fontsize=11
    )
    ax.grid(axis="y", color="#dddddd", linewidth=0.8, alpha=0.8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_qq(p_values, trait: str, max_points: int, out_path):
    """QQ plot of observed against expected -log10(P).

    Deviation from the diagonal at the tail is signal; deviation along its whole length
    suggests uncontrolled population structure or another systematic problem.
    """
    p_values = np.asarray(p_values, dtype=np.float64)
    p_values = p_values[np.isfinite(p_values)]
    p_values = np.clip(p_values, np.finfo(np.float64).tiny, 1.0)
    p_values = np.sort(p_values)

    # The expected quantile of a p-value is its rank within the FULL set of tests, so the
    # total is captured before any downsampling. Recomputing it over the truncated array
    # would rank the 250,000th-smallest p-value as if it were the largest of 250,000
    # tests, shifting every expected value down and lifting the whole curve off the
    # diagonal -- roughly 0.6 log units on perfectly null data, which reads as
    # genome-wide inflation that is not there.
    n_tests = p_values.size

    # Keeps the most significant points when downsampling -- the tail is the part
    # anyone reads, and rendering millions of points is slow and illegible.
    if n_tests > max_points:
        p_values = p_values[:max_points]

    expected = -np.log10((np.arange(1, p_values.size + 1) - 0.5) / n_tests)
    observed = -np.log10(p_values)
    upper = max(float(expected.max()), float(observed.max()), 1.0)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=8, color=_BLUE, alpha=0.75, linewidths=0)
    ax.plot([0, upper], [0, upper], color=_RED, linestyle="--", linewidth=1.2)
    ax.set_xlim(0, upper * 1.03)
    ax.set_ylim(0, upper * 1.03)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(f"{trait} QQ Plot")
    ax.grid(color="#dddddd", linewidth=0.8, alpha=0.8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def write_hits_table(bmi_df, t2d_df, hit_threshold: float, out_path: Path):
    """Significant variants, or the top 100 per trait when there are none.

    The fallback matters: an empty file is indistinguishable from a failed run.
    """
    tables = []
    for trait_df in (bmi_df, t2d_df):
        hits = trait_df.loc[trait_df["P"] < hit_threshold].copy()
        if hits.empty:
            hits = trait_df.nsmallest(100, "P").copy()
            hits["HIT_SET"] = "top100"
        else:
            hits = hits.sort_values("P").copy()
            hits["HIT_SET"] = f"p<{hit_threshold:g}"
        tables.append(hits)
    pd.concat(tables, ignore_index=True).to_csv(out_path, index=False)
