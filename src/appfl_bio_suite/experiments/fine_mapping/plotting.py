"""Publication figures for the fine-mapping stage. COORDINATOR-SIDE.

Ported from the standalone repository's ``scripts/plot_fine_mapping.py``, with the figure
code unchanged. Two things differ, both about where it runs rather than what it draws:

* :func:`write_figures` renders straight from a results DataFrame, so the aggregator can
  produce the figures in-process at the end of a run instead of leaving someone to
  remember a second command.
* Progress goes to a logger rather than to stdout, because during a federated run stdout
  belongs to APPFL.

Plotting is deliberately coordinator-side, exactly as it is for the GWAS experiment. It
keeps matplotlib off every partner's dependency surface, and it puts the diagnostics
where the person who needs them can actually see them.

Four figures, from one row per (locus, architecture, replicate):

  fig1_power_grid.png     recall (P[any causal captured]) across the star grid,
                          two OFAT arms: power vs h2 (rg=1.0) and power vs rg
                          (h2=0.001), one line per ncsl, Wilson 95% CIs.
  fig2_cs_size.png        resolution: distribution of the best (causal-containing)
                          credible-set size, by h2 arm.
  fig3_causal_pip.png     calibration: distribution of the causal variant's PIP,
                          by h2 arm.
  fig4_stratum_power.png  power by inter-site LD-divergence stratum x rg.

Architecture parameters are parsed out of ``architecture_id`` rather than looked up, so
this reads a results table and nothing else. That is what lets the same code draw the
centralized baseline's figures and the federated run's, and lets the two be compared by
overlaying files rather than by re-deriving anything.

    python -m appfl_bio_suite.experiments.fine_mapping.plotting \
        --results local/output/fine-mapping/data/fed_fm_results.tsv
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Progress goes here rather than to stdout: during a federated run stdout belongs to
# APPFL, and the aggregator passes its own logger in.
log = logging.getLogger(__name__)

# --- Okabe-Ito: the canonical CVD-safe qualitative palette for science ------
# ncsl identity (categorical, fixed order — never cycled):
NCSL_COLORS = {1: "#0072B2", 2: "#E69F00", 3: "#009E73"}  # blue / orange / green
# rg is ordinal (genetic correlation 0.5 -> 1.0): sequential single-hue blue ramp
RG_COLORS = {0.5: "#9ecae1", 0.7: "#4292c6", 1.0: "#08519c"}  # light -> dark
INK, MUTED, GRID = "#222222", "#666666", "#dddddd"

CENTER_H2 = 0.001  # the star-design centre: rg arm is taken at this h2
CENTER_RG = 1.0    # the h2 arm is taken at this rg
STRATUM_ORDER = ["high", "medium", "low"]  # inter-site LD divergence

_ARCH_RE = re.compile(r"ncsl(\d+)_h2-([\d.]+)_rg([\d.]+)")


def _apply_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "legend.frameon": False, "figure.facecolor": "white",
    })


def _parse_arch(df: pd.DataFrame) -> pd.DataFrame:
    """Add ncsl / h2 / rg columns parsed from architecture_id."""
    parsed = df["architecture_id"].astype(str).str.extract(_ARCH_RE)
    df = df.copy()
    df["ncsl"] = pd.to_numeric(parsed[0], errors="coerce").astype("Int64")
    df["h2"] = pd.to_numeric(parsed[1], errors="coerce")
    df["rg"] = pd.to_numeric(parsed[2], errors="coerce")
    bad = df["ncsl"].isna().sum()
    if bad:
        log.warning("%d row(s) had an unparseable architecture_id", bad)
    return df


def _wilson(k: int, n: int, z: float = 1.959963985) -> tuple[float, float, float]:
    """Wilson score interval for a binomial proportion (point, lo, hi)."""
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, centre - half), min(1.0, centre + half)


def _power_ci(g: pd.DataFrame) -> tuple[float, float, float, int]:
    n = len(g)
    k = int(g["any_causal_captured"].fillna(False).astype(bool).sum())
    p, lo, hi = _wilson(k, n)
    return p, lo, hi, n


# --------------------------------------------------------------------------- #
# Figure 1 — power across the star grid (two OFAT arms)
# --------------------------------------------------------------------------- #
def fig_power_grid(df: pd.DataFrame, out: Path) -> None:
    fig, (axh, axr) = plt.subplots(1, 2, figsize=(11, 4.4))

    # (a) power vs h2 at rg = CENTER_RG
    arm_h = df[df["rg"] == CENTER_RG]
    for ncsl in sorted(arm_h["ncsl"].dropna().unique()):
        sub = arm_h[arm_h["ncsl"] == ncsl]
        rows = [(h2, *_power_ci(g)) for h2, g in sub.groupby("h2")]
        rows.sort()
        if not rows:
            continue
        h2s = [r[0] for r in rows]
        p = np.array([r[1] for r in rows])
        lo = np.array([r[2] for r in rows])
        hi = np.array([r[3] for r in rows])
        c = NCSL_COLORS.get(int(ncsl), INK)
        axh.errorbar(h2s, p, yerr=[p - lo, hi - p], marker="o", ms=7, lw=2,
                     capsize=3, color=c, label=f"ncsl={int(ncsl)}")
    axh.set_xscale("log")
    h2_ticks = sorted(arm_h["h2"].dropna().unique())
    if h2_ticks:
        axh.set_xticks(h2_ticks)
        axh.set_xticklabels([f"{v:g}" for v in h2_ticks])
        axh.minorticks_off()
    axh.set_xlabel(f"target per-locus $h^2$  (rg = {CENTER_RG:.1f})")
    axh.set_ylabel("power  =  P(any causal captured)")
    axh.set_title("(a)  Recall vs heritability", loc="left")
    axh.set_ylim(-0.02, 1.02)
    axh.legend(title=None)

    # (b) power vs rg at h2 = CENTER_H2
    arm_r = df[df["h2"] == CENTER_H2]
    for ncsl in sorted(arm_r["ncsl"].dropna().unique()):
        sub = arm_r[arm_r["ncsl"] == ncsl]
        rows = [(rg, *_power_ci(g)) for rg, g in sub.groupby("rg")]
        rows.sort()
        if not rows:
            continue
        rgs = [r[0] for r in rows]
        p = np.array([r[1] for r in rows])
        lo = np.array([r[2] for r in rows])
        hi = np.array([r[3] for r in rows])
        c = NCSL_COLORS.get(int(ncsl), INK)
        axr.errorbar(rgs, p, yerr=[p - lo, hi - p], marker="o", ms=7, lw=2,
                     capsize=3, color=c, label=f"ncsl={int(ncsl)}")
    axr.set_xlabel(
        f"cross-ancestry genetic correlation $r_g$  ($h^2$ = {CENTER_H2:.4g})"
    )
    axr.set_title("(b)  Recall vs cross-ancestry $r_g$", loc="left")
    axr.set_ylim(-0.02, 1.02)
    axr.legend(title=None)

    fig.suptitle("SuSiEx cross-ancestry fine-mapping — recall across the causal grid",
                 fontsize=13, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out)
    plt.close(fig)
    log.info("wrote %s", out)


# --------------------------------------------------------------------------- #
# helper: box/strip distribution of a metric by h2 arm
# --------------------------------------------------------------------------- #
def _box_by_h2(df: pd.DataFrame, value: str, mask: pd.Series,
               ax: plt.Axes, ylabel: str) -> None:
    arm = df[(df["rg"] == CENTER_RG) & mask]
    h2s = sorted(arm["h2"].dropna().unique())
    data, labels, colors = [], [], []
    ramp = ["#c6dbef", "#6baed6", "#2171b5"]  # sequential blue by h2 rank
    for i, h2 in enumerate(h2s):
        v = arm[arm["h2"] == h2][value].dropna().values
        data.append(v)
        labels.append(f"{h2:g}\n(n={len(v)})")
        colors.append(ramp[min(i, len(ramp) - 1)])
    if not any(len(d) for d in data):
        ax.text(0.5, 0.5, "no captured instances", ha="center", va="center",
                transform=ax.transAxes, color=MUTED)
        return
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False,
                    medianprops=dict(color=INK, lw=1.6),
                    whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED),
                    boxprops=dict(edgecolor=MUTED))
    for patch, c in zip(bp["boxes"], colors, strict=False):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
    # jittered points on top for the real distribution
    rng = np.random.default_rng(0)
    for i, v in enumerate(data, start=1):
        if len(v):
            x = i + rng.uniform(-0.14, 0.14, size=len(v))
            ax.scatter(x, v, s=8, color=INK, alpha=0.18, zorder=3, linewidths=0)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlabel(f"target per-locus $h^2$  (rg = {CENTER_RG:.1f})")
    ax.set_ylabel(ylabel)


def fig_cs_size(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    captured = df["any_causal_captured"].fillna(False).astype(bool)
    _box_by_h2(df, "best_cs_size", captured, ax,
               "best credible-set size (variants)\n— smaller = sharper resolution")
    ax.set_title("Fine-mapping resolution among captured causals", loc="left")
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    log.info("wrote %s", out)


def fig_causal_pip(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    has_pip = df["causal_pip_max"].notna()
    _box_by_h2(df, "causal_pip_max", has_pip, ax,
               "causal-variant PIP (max over credible sets)")
    ax.axhline(0.95, color="#D55E00", lw=1.2, ls="--", zorder=1)
    ax.text(0.02, 0.955, "0.95", color="#D55E00", va="bottom",
            transform=ax.get_yaxis_transform(), fontsize=9)
    ax.set_title("Posterior confidence in the true causal variant", loc="left")
    ax.set_ylim(-0.02, 1.02)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    log.info("wrote %s", out)


# --------------------------------------------------------------------------- #
# Figure 4 — power by LD-divergence stratum x rg
# --------------------------------------------------------------------------- #
def fig_stratum_power(df: pd.DataFrame, out: Path) -> None:
    if "stratum" not in df.columns or df["stratum"].isna().all():
        log.info("skipping fig4: the results carry no stratum column")
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    strata = [s for s in STRATUM_ORDER if s in set(df["stratum"].dropna())]
    rgs = sorted(df["rg"].dropna().unique())
    width = 0.8 / max(len(rgs), 1)
    x = np.arange(len(strata))
    for j, rg in enumerate(rgs):
        p, lo, hi = [], [], []
        for s in strata:
            g = df[(df["stratum"] == s) & (df["rg"] == rg)]
            point, low, high, _ = _power_ci(g)
            p.append(point)
            lo.append(point - low if not np.isnan(point) else 0)
            hi.append(high - point if not np.isnan(point) else 0)
        c = RG_COLORS.get(round(float(rg), 1), "#4292c6")
        ax.bar(x + j * width - 0.4 + width / 2, p, width * 0.92,
               color=c, label=f"$r_g$ = {rg:g}", yerr=[lo, hi],
               capsize=2, ecolor=MUTED, error_kw=dict(lw=1))
    ax.set_xticks(x)
    ax.set_xticklabels([s.capitalize() for s in strata])
    ax.set_xlabel("inter-site LD-divergence stratum")
    ax.set_ylabel("power  =  P(any causal captured)")
    ax.set_ylim(0, 1.02)
    ax.set_title("Recall by cross-site LD divergence and genetic correlation", loc="left")
    ax.legend(title=None, ncol=len(rgs))
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    log.info("wrote %s", out)


def write_figures(results, out_dir, logger=None) -> list[Path]:
    """Render all four figures from a results DataFrame. The aggregator's entry point.

    Returns the paths written. A figure that cannot be drawn from these results -- no
    stratum column, no captured causals -- is skipped with a note rather than failing the
    run: the numbers are the deliverable and the figures are a convenience, so a plotting
    problem must never lose a fine-mapping run that took hours.
    """
    active = logger or log
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    frame = _parse_arch(pd.DataFrame(results))
    if frame.empty:
        active.warning("no results to plot")
        return []

    _apply_style()
    written = []
    for name, draw in (
        ("fig1_power_grid.png", fig_power_grid),
        ("fig2_cs_size.png", fig_cs_size),
        ("fig3_causal_pip.png", fig_causal_pip),
        ("fig4_stratum_power.png", fig_stratum_power),
    ):
        target = out_dir / name
        try:
            draw(frame, target)
        except Exception as exc:  # noqa: BLE001 - a figure must not lose a run
            active.warning("could not draw %s: %s", name, exc)
            continue
        if target.exists():
            written.append(target)
    active.info("wrote %d figure(s) -> %s", len(written), out_dir)
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Render the fine-mapping figures from a results table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--results",
        default="local/output/fine-mapping/data/fed_fm_results.tsv",
        help="per-instance results TSV, from a federated or a centralized run",
    )
    ap.add_argument(
        "--out-dir", default=None,
        help="figure output directory (default: <results dir>/../graphs)",
    )
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    results = Path(args.results)
    if not results.exists():
        print(
            f"ERROR: {results} not found -- no fine-mapping results to plot yet.\n"
            "Produce some with:\n"
            "    appfl-bio-suite run fine-mapping --config loopback --driver serial \\\n"
            "        --data-root <the directory `simulate fine-mapping --out` wrote>",
            file=sys.stderr,
        )
        return 2

    out_dir = Path(args.out_dir) if args.out_dir else results.parent.parent / "graphs"
    frame = pd.read_csv(results, sep="\t")
    log.info(
        "loaded %d instance(s), %d architecture(s)",
        len(frame), frame["architecture_id"].nunique(),
    )
    write_figures(frame, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
