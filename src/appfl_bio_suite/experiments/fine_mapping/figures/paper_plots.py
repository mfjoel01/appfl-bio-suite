"""Figures in the form the SuSiEx paper uses, over this experiment's results.

COORDINATOR-SIDE. Reference: Yuan et al., "Fine-mapping across diverse ancestries
drives the discovery of putative causal variants underlying human complex traits and
diseases" (medRxiv 2023.01.07.23284293v4).

WHAT MAPS ACROSS, AND WHAT DOES NOT
-------------------------------------
The paper varies the DISCOVERY CONFIGURATION -- which populations were combined, at
what sample size -- and reports calibration, power and resolution against it. This
experiment holds the configuration fixed (six ancestry columns over 150,000 people)
and varies the LOCUS: how much the three sites' LD disagrees there, and what causal
architecture was simulated on it. So the paper's x-axis becomes ours:

    paper                             here
    which populations were combined   cross-site LD-divergence stratum
    discovery sample size             per-locus h2 (signal available)
    cross-population r_g              effect-draw correlation parameter
    method (SuSiEx / PAINTOR / ...)   analysis path (federated / centralized)

Its metrics carry over unchanged, and are the reason to follow its forms at all:

    coverage    the share of credible sets containing a true causal variant. THE
                calibration number. Not derivable from the summary results table --
                see the note on pap3 -- so it comes from the harvested per-CS table.
    power       the share of instances where a causal variant was captured.
    resolution  credible-set size, and how many causals reach PIP > 0.95.

FIGURES, AND WHICH PAPER FIGURE EACH ANSWERS TO

    pap1_recall_grid          Supp. Figs 10-24  power along the two design arms
    pap2_diversity_panels     Figure 2 (a-d)    power / high-PIP / coverage /
                                                resolution against LD divergence
    pap3_convergence          Figure 3a         did the fit succeed, and how often
    pap4_power_coverage       Figure 3b         the power-calibration trade-off
    pap5_resolution           Figure 4a, 4b     PIP and credible-set distributions
    pap6_high_confidence      Figure 4c, 4d     counts at PIP > 0.95, incl. the
                                                single-credible-set restriction
    pap7_ld_divergence        Supp. Fig 39      the LD-mismatch analogue
    pap8_scalability          Supp. Figs 8-9    runtime, and what drives it
    pap9_coding_quality       Ext. Data Fig 6   quality issues against PIP

Extended Data Figure 7 (VEP functional impact) has no analogue and is not attempted:
the genotypes are synthetic, so no variant carries a functional annotation to bin by.
Figure 1 and Extended Data Figure 1 are method schematics rather than plots of data.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from appfl_bio_suite.experiments.fine_mapping.figures import figstyle as fs
from appfl_bio_suite.experiments.fine_mapping.figures.figstyle import plt

log = logging.getLogger(__name__)

STRATUM_ORDER = ["low", "medium", "high"]
CENTER_RG = 1.0  # the h2 arm is taken at this rg
CENTER_H2 = 0.001  # the rg arm is taken at this h2
_ARCH_RE = re.compile(r"ncsl(\d+)_h2-([\d.]+)_rg([\d.]+)")
# Marker areas for the three heritability levels, by rank rather than by value.
_H2_SIZE = {0.0005: 45, 0.001: 120, 0.005: 300}


# --------------------------------------------------------------------------- #
# preparation
# --------------------------------------------------------------------------- #
def parse_architecture(df: pd.DataFrame) -> pd.DataFrame:
    """Add ncsl / h2 / rg columns parsed out of ``architecture_id``.

    Parsed rather than joined against the manifest so that this module reads a results
    table and nothing else -- which is what lets the same code draw the centralized
    baseline's figures and a federated run's, and lets the two be compared by
    overlaying files rather than by re-deriving anything.
    """
    parsed = df["architecture_id"].astype(str).str.extract(_ARCH_RE)
    df = df.copy()
    df["ncsl"] = pd.to_numeric(parsed[0], errors="coerce").astype("Int64")
    df["h2"] = pd.to_numeric(parsed[1], errors="coerce")
    df["rg"] = pd.to_numeric(parsed[2], errors="coerce")
    bad = int(df["ncsl"].isna().sum())
    if bad:
        log.warning("%d row(s) had an unparseable architecture_id", bad)
    # Derived flags every figure below wants, defined once so "captured" means the
    # same thing in all of them.
    df["captured"] = df.get("any_causal_captured", False).fillna(False).astype(bool)
    df["returned_cs"] = df.get("n_credible_sets", 0).fillna(0) > 0
    df["failed"] = df.get("error", "").fillna("").astype(str).str.len() > 0
    for col, thresh in (("pip95", 0.95), ("pip50", 0.50)):
        has = "causal_pip_max" in df
        df[col] = (df["causal_pip_max"].fillna(-1) > thresh) if has else False
    df["nonconverged"] = df.get("fit_status", pd.Series("unknown", index=df.index)).eq(
        "nonconverged"
    )
    df["converged_no_cs"] = df.get("fit_status", pd.Series("unknown", index=df.index)).eq(
        "converged_no_cs"
    )
    df["single_cs"] = df.get("n_credible_sets", 0) == 1
    return df


def _strata_present(df: pd.DataFrame) -> list[str]:
    have = set(df.get("stratum", pd.Series(dtype=object)).dropna())
    return [s for s in STRATUM_ORDER if s in have]


def _rate(sub: pd.DataFrame, col: str) -> tuple[float, float, float, int]:
    n = len(sub)
    from ..reporting import clustered_ratio

    d = sub.copy()
    d["_hits"] = sub[col].fillna(False).astype(int)
    d["_one"] = 1
    p, lo, hi = clustered_ratio(d, "_hits", "_one")
    return p, lo, hi, n


def _err_bars(rows):
    """``[(point, lo, hi), ...]`` -> the asymmetric yerr matplotlib wants."""
    p = np.array([r[0] for r in rows], dtype=float)
    lo = np.array([r[1] for r in rows], dtype=float)
    hi = np.array([r[2] for r in rows], dtype=float)
    return p, np.vstack([p - lo, hi - p])


# --------------------------------------------------------------------------- #
# pap1 -- recall along the two design arms   (paper Supp. Figs 10-24)
# --------------------------------------------------------------------------- #
def pap1_recall_grid(df: pd.DataFrame, out: Path) -> Path:
    fig, (axh, axr) = plt.subplots(1, 2, figsize=(11.4, 4.6))
    ncsls = sorted(int(v) for v in df["ncsl"].dropna().unique())
    colors = dict(zip(ncsls, fs.ordinal_colors(len(ncsls)), strict=False))

    # Log scale FIRST: set_xscale rebuilds the locator and discards ticks set before
    # it, which is how the 0.0005 and 0.005 labels went missing.
    axh.set_xscale("log")
    axh.minorticks_off()
    for ax, arm_col, fixed_col, fixed_val, xlabel, letter, title in (
        (
            axh,
            "h2",
            "rg",
            CENTER_RG,
            f"target per-locus $h^2$        ($r_g$ = {CENTER_RG:g})",
            "a",
            "Recall against heritability",
        ),
        (
            axr,
            "rg",
            "h2",
            CENTER_H2,
            f"cross-ancestry $r_g$        ($h^2$ = {CENTER_H2:g})",
            "b",
            "Recall against genetic correlation",
        ),
    ):
        arm = df[df[fixed_col] == fixed_val]
        ends = []
        for ncsl in ncsls:
            sub = arm[arm["ncsl"] == ncsl]
            rows = sorted((x, *_rate(g, "captured")) for x, g in sub.groupby(arm_col))
            if not rows:
                continue
            xs = [r[0] for r in rows]
            p, yerr = _err_bars([(r[1], r[2], r[3]) for r in rows])
            ax.errorbar(
                xs, p, yerr=yerr, marker="o", ms=7, lw=2, capsize=3, color=colors[ncsl], zorder=4
            )
            ends.append((float(p[-1]), f"$n_{{csl}}$ = {ncsl}", colors[ncsl], xs[-1]))
        ax.set_xlabel(xlabel)
        ax.set_ylim(0.0, 1.02)
        ticks = sorted(arm[arm_col].dropna().unique())
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{v:g}" for v in ticks])
        # Direct labels rather than a legend box -- three series, and a reader should
        # not have to cross the panel to learn which line is which. Spread, because
        # power saturates and the three ends land within a pixel of each other.
        if ends:
            fs.spread_labels(ax, [(y, txt, c) for y, txt, c, _ in ends], ends[0][3])
        fs.panel_letter(ax, letter, title)
    axh.set_xlim(right=max(sorted(df[df["rg"] == CENTER_RG]["h2"].dropna().unique())) * 2.6)
    axr.set_xlim(right=max(sorted(df[df["h2"] == CENTER_H2]["rg"].dropna().unique())) + 0.22)
    axh.set_ylabel("power   =   P(a causal variant is captured)")

    fs.footnote(
        fig,
        f"{len(df):,} instances. Bars are Wilson 95% intervals. The two "
        "arms share their centre cell, so the rightmost point of a is the "
        "same data as the rightmost point of b for each $n_{csl}$.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap2 -- the paper's Figure 2, against LD divergence
# --------------------------------------------------------------------------- #
def pap2_diversity_panels(
    df: pd.DataFrame, out: Path, cs_detail: pd.DataFrame | None = None
) -> Path:
    """Power, high-confidence yield, calibration and resolution against LD divergence.

    The paper's Figure 2 puts "which populations were combined" on the x-axis, because
    that is the lever it controls. Here the lever is the locus: how far apart the three
    sites' LD is where the signal sits. Colour is r_g, the one factor the paper and this
    design share outright.

    Panel c is the calibration panel and it is the one that needs the harvested per-CS
    table. See :func:`pap3_convergence` for why the summary table cannot supply it.
    """
    arm = df[df["h2"] == CENTER_H2]
    if arm.empty:
        arm = df
    strata = _strata_present(arm)
    rgs = sorted(arm["rg"].dropna().unique())
    colors = dict(zip(rgs, fs.ordinal_colors(len(rgs)), strict=False))
    x = np.arange(len(strata))
    width = 0.8 / max(len(rgs), 1)

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.2))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    def _grouped_bar(ax, value_fn, ylabel, letter, title, pct=False, ref=None):
        for j, rg in enumerate(rgs):
            pts, los, his, ns = [], [], [], []
            for s in strata:
                g = arm[(arm["stratum"] == s) & (arm["rg"] == rg)]
                p, lo, hi, n = value_fn(g)
                pts.append(p)
                los.append(p - lo if np.isfinite(p) else 0)
                his.append(hi - p if np.isfinite(p) else 0)
                ns.append(n)
            ax.bar(
                x + j * width - 0.4 + width / 2,
                pts,
                width * 0.88,
                color=colors[rg],
                edgecolor=fs.SURFACE,
                linewidth=1.4,
                yerr=[los, his],
                capsize=2,
                ecolor=fs.INK2,
                error_kw=dict(lw=1.0),
                zorder=3,
            )
        if ref is not None:
            ax.axhline(ref, color=fs.INK, lw=1.3, ls="--", zorder=5)
        ax.set_xticks(x)
        ax.set_xticklabels([s.capitalize() for s in strata])
        ax.set_ylabel(ylabel)
        if pct:
            ax.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
        fs.panel_letter(ax, letter, title)

    # (a) power -- the share of instances where a causal variant landed in a set
    _grouped_bar(
        ax_a,
        lambda g: _rate(g, "captured"),
        "power   =   P(a causal variant is captured)",
        "a",
        "Power",
        pct=True,
    )
    ax_a.set_ylim(0, 1.18)
    fs.legend_swatches(
        ax_a,
        {f"$r_g$ = {r:g}": colors[r] for r in rgs},
        loc="upper center",
        ncol=len(rgs),
        bbox_to_anchor=(0.5, 1.02),
    )

    # (b) resolution at the top end -- how often the causal variant reaches PIP > 0.95
    _grouped_bar(
        ax_b,
        lambda g: _rate(g, "pip95"),
        "P(causal variant at PIP > 0.95)",
        "b",
        "High-confidence yield",
        pct=True,
    )
    ax_b.set_ylim(0, max(0.35, ax_b.get_ylim()[1]))

    # (c) calibration -- credible sets containing a true causal variant
    if cs_detail is not None and len(cs_detail):
        d = cs_detail.copy()
        for j, rg in enumerate(rgs):
            pts, los, his = [], [], []
            for s in strata:
                g = d[(d["stratum"] == s) & (d["rg"] == rg)]
                p, lo, hi, n = _rate(g, "contains_causal")
                pts.append(p)
                los.append(p - lo if np.isfinite(p) else 0)
                his.append(hi - p if np.isfinite(p) else 0)
            ax_c.bar(
                x + j * width - 0.4 + width / 2,
                pts,
                width * 0.88,
                color=colors[rg],
                edgecolor=fs.SURFACE,
                linewidth=1.4,
                yerr=[los, his],
                capsize=2,
                ecolor=fs.INK2,
                error_kw=dict(lw=1.0),
                zorder=3,
            )
        ax_c.axhline(0.95, color=fs.INK, lw=1.3, ls="--", zorder=5)
        ax_c.text(
            0.005,
            0.95,
            " nominal 95%",
            transform=ax_c.get_yaxis_transform(),
            ha="left",
            va="bottom",
            fontsize=8.5,
            color=fs.INK2,
        )
        ax_c.set_xticks(x)
        ax_c.set_xticklabels([s.capitalize() for s in strata])
        ax_c.set_ylabel("coverage   =   P(credible set holds a causal variant)")
        ax_c.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
        ax_c.set_ylim(0.5, 1.02)
        fs.panel_letter(ax_c, "c", f"Calibration  ({len(d):,} credible sets)")
    else:
        fs.no_data(
            ax_c,
            "coverage needs the harvested per-credible-set table\n"
            "(figures.detail); the summary results table\n"
            "cannot attribute a causal variant to a set",
        )
        fs.panel_letter(ax_c, "c", "Calibration")

    # (d) resolution -- credible-set size, among sets that caught something
    data, positions, cols = [], [], []
    for i, s in enumerate(strata):
        for j, rg in enumerate(rgs):
            g = arm[(arm["stratum"] == s) & (arm["rg"] == rg) & arm["captured"]]
            v = g["best_cs_size"].dropna().values
            positions.append(i + j * width - 0.4 + width / 2)
            data.append(v)
            cols.append(colors[rg])
    fs.boxplot(ax_d, data, positions=positions, colors=cols, widths=width * 0.82, points=False)
    ax_d.set_yscale("log")
    # Most of these boxes have Q1 = median = 1. On a log axis whose floor is also 1
    # they render as solid bars rising from the spine and stop reading as boxes.
    ax_d.set_ylim(bottom=0.82)
    ax_d.set_xticks(x)
    ax_d.set_xticklabels([s.capitalize() for s in strata])
    ax_d.set_xlim(-0.55, len(strata) - 0.45)
    ax_d.set_ylabel("credible-set size (variants)\nsmaller = sharper resolution")
    fs.panel_letter(ax_d, "d", "Resolution")

    for ax in (ax_a, ax_b, ax_c, ax_d):
        ax.set_xlabel("cross-site LD-divergence stratum")
        ax.grid(axis="x", visible=False)

    fig.suptitle(
        "Fine-mapping against cross-site LD divergence",
        x=0.005,
        ha="left",
        fontsize=13,
        fontweight="semibold",
    )
    cs_note = (
        f" Panel c is computed over the {len(cs_detail):,} credible sets of the "
        "harvested sample rather than the full sweep, because the summary "
        "results table cannot attribute a causal variant to a particular set; "
        "its wider intervals are that smaller sample, not worse calibration."
        if cs_detail is not None and len(cs_detail)
        else ""
    )
    fs.footnote(
        fig,
        f"Panels a, b and d: the $r_g$ arm of the design "
        f"($h^2$ = {CENTER_H2:g}), {len(arm):,} instances. Error bars are "
        "Wilson 95% intervals. Panel d is restricted to instances that "
        "captured a causal variant, since a set that missed has no "
        "meaningful size." + cs_note,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap3 -- did the fit succeed?   (paper Figure 3a)
# --------------------------------------------------------------------------- #
def pap3_convergence(df: pd.DataFrame, out: Path) -> Path:
    """Job outcomes per design cell, the paper's robustness panel.

    The paper's categories are "successful / unreasonable / unfinished", where
    "unreasonable" means a PIP sum far from the number of simulated causals -- a
    diagnosis aimed at PAINTOR and MsCAVIAR, whose failures are silent. SuSiEx's
    failure mode is different and visible: it either returns credible sets, or
    converges and returns none, or errors. Those three are the categories here, and
    relabelling them to match the paper's would misreport what happened.

    "Converged, no credible set" is not a failure. At h2 = 0.0005 over a 1 Mb window
    there is often genuinely no resolvable signal, and returning nothing is the correct
    answer; the panel exists to show that this tracks heritability rather than
    appearing at random, which is what it would do if it were a bug.
    """
    order = (
        df.groupby("architecture_id")[["h2", "rg", "ncsl"]]
        .first()
        .sort_values(["h2", "rg", "ncsl"])
        .index.tolist()
    )
    cats = [
        (
            "returned credible set(s)",
            fs.STATUS["good"],
            lambda g: (g["returned_cs"] & ~g["failed"]).sum(),
        ),
        (
            "converged, no credible set",
            fs.STATUS["warning"],
            lambda g: g["converged_no_cs"].sum(),
        ),
        ("nonconverged", fs.STATUS["critical"], lambda g: g["nonconverged"].sum()),
        ("error", fs.STATUS["critical"], lambda g: g["failed"].sum()),
    ]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.6, 5.2), gridspec_kw={"width_ratios": [2.1, 1]}
    )
    y = np.arange(len(order))
    left = np.zeros(len(order))
    for label, color, fn in cats:
        vals = np.array([fn(df[df["architecture_id"] == a]) for a in order], dtype=float)
        tot = np.array([len(df[df["architecture_id"] == a]) for a in order], dtype=float)
        frac = np.divide(vals, tot, out=np.zeros_like(vals), where=tot > 0)
        ax1.barh(
            y,
            frac,
            left=left,
            height=0.66,
            color=color,
            edgecolor=fs.SURFACE,
            linewidth=1.4,
            label=label,
        )
        for yi, (f, offset, v) in enumerate(zip(frac, left, vals, strict=False)):
            if f >= 0.055:
                ax1.text(
                    offset + f / 2,
                    yi,
                    f"{int(v)}",
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    color=fs.SURFACE,
                    fontweight="semibold",
                )
        left += frac
    ax1.set_yticks(y)
    ax1.set_yticklabels(order, fontsize=8)
    ax1.set_xlim(0, 1)
    ax1.xaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    ax1.set_xlabel("share of instances")
    ax1.invert_yaxis()
    ax1.grid(axis="y", visible=False)
    ax1.legend(loc="upper left", bbox_to_anchor=(0, -0.13), ncol=3)
    fs.title_with_letter(ax1, "a", "Outcome by design cell", dx=-0.30)

    # (b) the same outcome against heritability alone -- the panel that shows the
    # empty-result rate is a property of the signal, not of the software.
    h2s = sorted(df["h2"].dropna().unique())
    bottom = np.zeros(len(h2s))
    for _label, color, fn in cats:
        vals = np.array([fn(df[df["h2"] == h]) for h in h2s], dtype=float)
        tot = np.array([len(df[df["h2"] == h]) for h in h2s], dtype=float)
        frac = np.divide(vals, tot, out=np.zeros_like(vals), where=tot > 0)
        ax2.bar(
            range(len(h2s)),
            frac,
            bottom=bottom,
            width=0.6,
            color=color,
            edgecolor=fs.SURFACE,
            linewidth=1.4,
        )
        for xi, (f, b) in enumerate(zip(frac, bottom, strict=False)):
            if f >= 0.06:
                ax2.text(
                    xi,
                    b + f / 2,
                    f"{f * 100:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=8.5,
                    color=fs.SURFACE,
                    fontweight="semibold",
                )
        bottom += frac
    ax2.set_xticks(range(len(h2s)))
    ax2.set_xticklabels([f"{h:g}" for h in h2s])
    ax2.set_xlabel("target per-locus $h^2$")
    ax2.set_ylim(0, 1)
    ax2.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    ax2.grid(axis="x", visible=False)
    fs.panel_letter(ax2, "b", "Outcome against heritability")

    n_err = int(df["failed"].sum())
    fs.footnote(
        fig,
        f"{len(df):,} instances, {n_err} error(s). An empty result is the "
        "correct answer where no signal is resolvable, which is why it "
        "concentrates at the lowest heritability rather than spreading "
        "evenly across the design.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap4 -- power against calibration   (paper Figure 3b)
# --------------------------------------------------------------------------- #
def pap4_power_coverage(df: pd.DataFrame, out: Path, cs_detail: pd.DataFrame | None = None) -> Path:
    """The paper's most compact summary: yield on one axis, calibration on the other.

    A method can buy power by inflating credible sets, and buy calibration by refusing
    to commit; only the joint position says whether it did either. Each point is one
    (design cell x LD-divergence stratum) group.
    """
    if cs_detail is None or cs_detail.empty:
        fig, ax = plt.subplots(figsize=(8.0, 6.0))
        fs.no_data(ax, "Calibration unavailable: per-credible-set membership is required")
        return fs.save(fig, out, log)
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    strata = _strata_present(df)
    colors = dict(zip(strata, fs.ordinal_colors(len(strata), fs.VIOLET_RAMP), strict=False))
    use_exact = cs_detail is not None and len(cs_detail) > 0

    rows = []
    for (arch, s), g in df.groupby(["architecture_id", "stratum"]):
        if s not in colors or len(g) < 20:
            continue
        yield_rate, _, _, n = _rate(g, "pip50")
        if use_exact:
            d = cs_detail[(cs_detail["architecture_id"] == arch) & (cs_detail["stratum"] == s)]
            if len(d) < 10:
                continue
            cov, lo, hi, ncs = _rate(d, "contains_causal")
        else:
            cov, lo, hi, ncs = _rate(g, "captured")
        rows.append((yield_rate, cov, lo, hi, n, ncs, s, arch, g["h2"].iloc[0]))

    if not rows:
        fs.no_data(ax)
        return fs.save(fig, out, log)

    for s in strata:
        pts = [r for r in rows if r[6] == s]
        if not pts:
            continue
        ax.errorbar(
            [p[0] for p in pts],
            [p[1] for p in pts],
            yerr=[[p[1] - p[2] for p in pts], [p[3] - p[1] for p in pts]],
            fmt="none",
            ecolor=colors[s],
            elinewidth=0.9,
            alpha=0.6,
            zorder=2,
        )
        # Size carries the group's h2. Mapped through the RANK of the level rather
        # than its value: h2 spans 0.0005 to 0.005, so any area proportional to the
        # value either collapses the small levels to a dot or makes the large one
        # swamp the panel. Three levels, three legible sizes.
        sizes = [_H2_SIZE.get(round(float(p[8]), 6), 90) for p in pts]
        ax.scatter(
            [p[0] for p in pts],
            [p[1] for p in pts],
            s=sizes,
            color=colors[s],
            alpha=0.85,
            lw=0.8,
            edgecolor=fs.SURFACE,
            zorder=4,
        )

    ax.axhline(0.95, color=fs.INK, lw=1.3, ls="--", zorder=5)
    ax.text(
        0.995,
        0.95,
        "nominal 95% coverage ",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=9,
        color=fs.INK2,
    )
    ax.set_xlabel("yield   =   P(causal variant reaches PIP > 0.5)")
    ax.set_ylabel(
        "coverage" + ("   (per credible set)" if use_exact else "   (per instance, power)")
    )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(min(0.55, min(r[1] for r in rows) - 0.05), 1.02)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    ax.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    ax.set_title("Yield against calibration, per design cell", loc="left")
    leg = fs.legend_swatches(
        ax, colors, title="LD-divergence stratum", loc="lower right", marker="o"
    )
    ax.add_artist(leg)
    from matplotlib.lines import Line2D

    h2s = sorted(df["h2"].dropna().unique())
    ax.legend(
        handles=[
            Line2D(
                [],
                [],
                ls="none",
                marker="o",
                color=fs.MUTED,
                markersize=np.sqrt(_H2_SIZE.get(round(float(h), 6), 90)),
                label=f"$h^2$ = {h:g}",
            )
            for h in h2s
        ],
        title="marker size",
        loc="lower left",
        labelspacing=1.3,
        borderpad=0.8,
        handletextpad=1.4,
    )

    src = (
        "credible sets, harvested"
        if use_exact
        else "instances (the summary table cannot attribute a causal variant to a "
        "specific credible set, so this axis is per-instance power; run "
        "figures.detail for true per-set coverage)"
    )
    fs.footnote(
        fig,
        f"One point per design cell x stratum, {len(rows)} points. "
        f"Coverage computed over {src}. Bars are Wilson 95% intervals.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap5 -- resolution distributions   (paper Figure 4a / 4b)
# --------------------------------------------------------------------------- #
def pap5_resolution(df: pd.DataFrame, out: Path) -> Path:
    """Where the posterior mass ends up, and how wide the sets are.

    Violin rather than box alone because the SHAPE is the finding: the PIP
    distribution is strongly bimodal -- fine-mapping either resolves a locus or it
    does not, and the middle is nearly empty. A box plot reports a median around
    which almost no instance actually sits.
    """
    h2s = sorted(df["h2"].dropna().unique())
    arm = df[df["rg"] == CENTER_RG]
    colors = fs.ordinal_colors(len(h2s))

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.7))
    ax_a, ax_b, ax_c = axes

    # (a) causal-variant PIP
    data = [arm[(arm["h2"] == h)]["causal_pip_max"].dropna().values for h in h2s]
    fs.violin_box(ax_a, data, colors=colors)
    ax_a.axhline(0.95, color=fs.STATUS["critical"], lw=1.1, ls="--", zorder=6)
    ax_a.text(
        0.01,
        0.955,
        " 0.95",
        transform=ax_a.get_yaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=fs.STATUS["critical"],
    )
    ax_a.set_ylim(-0.03, 1.05)
    ax_a.set_ylabel("PIP of the true causal variant")
    fs.panel_letter(ax_a, "a", "Posterior confidence")

    # (b) credible-set size, log scale -- the distribution spans 1 to >100
    data = [arm[(arm["h2"] == h) & arm["captured"]]["best_cs_size"].dropna().values for h in h2s]
    fs.violin_box(ax_b, data, colors=colors)
    ax_b.set_yscale("log")
    ax_b.set_ylabel("credible-set size (variants)")
    fs.panel_letter(ax_b, "b", "Resolution")

    # (c) purity -- the paper's own filter criterion, worth showing rather than
    # assuming: a set below 0.5 purity carries no inferential value.
    data = [arm[arm["h2"] == h]["mean_cs_purity"].dropna().values for h in h2s]
    fs.violin_box(ax_c, data, colors=colors)
    ax_c.axhline(0.5, color=fs.STATUS["critical"], lw=1.1, ls="--", zorder=6)
    ax_c.text(
        0.99,
        0.505,
        "0.5 purity floor ",
        transform=ax_c.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=8.5,
        color=fs.STATUS["critical"],
    )
    ax_c.set_ylabel("mean credible-set purity")
    fs.panel_letter(ax_c, "c", "Purity")

    for ax in axes:
        ax.set_xticks(range(1, len(h2s) + 1))
        ax.set_xticklabels([f"{h:g}" for h in h2s])
        ax.set_xlabel(f"target per-locus $h^2$   ($r_g$ = {CENTER_RG:g})", labelpad=20)
        ax.grid(axis="x", visible=False)
    for ax, ds in (
        (ax_a, [arm[arm["h2"] == h]["causal_pip_max"].dropna() for h in h2s]),
        (ax_b, [arm[(arm["h2"] == h) & arm["captured"]]["best_cs_size"].dropna() for h in h2s]),
        (ax_c, [arm[arm["h2"] == h]["mean_cs_purity"].dropna() for h in h2s]),
    ):
        fs.count_labels(ax, range(1, len(ds) + 1), [len(d) for d in ds])

    fig.suptitle(
        "Where fine-mapping puts its confidence",
        x=0.005,
        ha="left",
        fontsize=12.5,
        fontweight="semibold",
    )
    fs.footnote(
        fig,
        "Red dot is the mean, the dark bar the interquartile range, the "
        "violin the full distribution. The heritability arm only "
        f"($r_g$ = {CENTER_RG:g}).",
        y=-0.055,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap6 -- high-confidence counts   (paper Figure 4c / 4d)
# --------------------------------------------------------------------------- #
def pap6_high_confidence(df: pd.DataFrame, out: Path) -> Path:
    """Counts of confidently fine-mapped variants, with the paper's own restriction.

    Panel b repeats panel a over single-credible-set loci only. The paper introduces
    that restriction (its Figure 4d) so that the comparison is not confounded by
    multiple causal variants in LD or by the algorithm that aligns sets across
    analyses; the same confound applies here whenever ncsl > 1.
    """
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.6))
    ax_a, ax_b, ax_c = axes
    h2s = sorted(df["h2"].dropna().unique())
    arm = df[df["rg"] == CENTER_RG]
    colors = fs.ordinal_colors(len(h2s))

    for ax, sub, letter, title in (
        (ax_a, arm, "a", "All loci"),
        (ax_b, arm[arm["single_cs"]], "b", "Single-credible-set loci only"),
    ):
        x = np.arange(len(h2s))
        for k, (col, _lbl, alpha) in enumerate(
            (("pip50", "PIP > 0.5", 0.45), ("pip95", "PIP > 0.95", 1.0))
        ):
            vals = [int(sub[sub["h2"] == h][col].sum()) for h in h2s]
            ax.bar(
                x + (k - 0.5) * 0.36,
                vals,
                0.33,
                color=[colors[i] for i in range(len(h2s))],
                alpha=alpha,
                edgecolor=fs.SURFACE,
                linewidth=1.4,
                zorder=3,
            )
            for xi, v in zip(x + (k - 0.5) * 0.36, vals, strict=False):
                ax.text(xi, v, f"{v:,}", ha="center", va="bottom", fontsize=8, color=fs.INK2)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{h:g}" for h in h2s])
        ax.set_xlabel(f"target per-locus $h^2$   ($r_g$ = {CENTER_RG:g})")
        ax.set_ylabel("causal variants identified")
        ax.grid(axis="x", visible=False)
        # Bottom-left: the tall bar is at the right of these panels, so a caption
        # anchored top-left is the only corner the counts never reach.
        ax.set_ylim(top=ax.get_ylim()[1] * 1.10)
        ax.text(
            0.02,
            0.02,
            "left bar: PIP > 0.5 (pale)\nright bar: PIP > 0.95 (solid)",
            transform=ax.transAxes,
            va="bottom",
            ha="left",
            fontsize=8.5,
            color=fs.INK2,
        )
        fs.panel_letter(ax, letter, title)

    # (c) the same thing as a rate, so the two panels above can be compared despite
    # resting on different denominators.
    x = np.arange(len(h2s))
    for k, (sub, _lbl, hatch) in enumerate(
        ((arm, "all loci", None), (arm[arm["single_cs"]], "single-CS loci", "///"))
    ):
        rows = [_rate(sub[sub["h2"] == h], "pip95") for h in h2s]
        p, yerr = _err_bars([(r[0], r[1], r[2]) for r in rows])
        ax_c.bar(
            x + (k - 0.5) * 0.36,
            p,
            0.33,
            yerr=yerr,
            capsize=2,
            ecolor=fs.INK2,
            color=[colors[i] for i in range(len(h2s))],
            edgecolor=fs.SURFACE,
            linewidth=1.4,
            hatch=hatch,
            zorder=3,
            error_kw=dict(lw=1.0),
        )
    ax_c.set_xticks(x)
    ax_c.set_xticklabels([f"{h:g}" for h in h2s])
    ax_c.set_xlabel(f"target per-locus $h^2$   ($r_g$ = {CENTER_RG:g})")
    ax_c.set_ylabel("P(causal variant at PIP > 0.95)")
    ax_c.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    ax_c.grid(axis="x", visible=False)
    ax_c.text(
        0.02,
        0.97,
        "left bar: all loci\nright bar: single-CS loci (hatched)",
        transform=ax_c.transAxes,
        va="top",
        ha="left",
        fontsize=8.5,
        color=fs.INK2,
    )
    ax_c.set_ylim(top=ax_c.get_ylim()[1] * 1.28)
    fs.panel_letter(ax_c, "c", "As a rate")

    fig.suptitle(
        "Confidently fine-mapped variants", x=0.005, ha="left", fontsize=12.5, fontweight="semibold"
    )
    fs.footnote(
        fig,
        f"{len(arm):,} instances in the heritability arm; "
        f"{int(arm['single_cs'].sum()):,} of them resolved to exactly one "
        "credible set. Colour is heritability, repeating the x position.",
        y=-0.05,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap7 -- the LD-mismatch analogue   (paper Supp. Fig 39)
# --------------------------------------------------------------------------- #
def locus_ld_strength(
    data_root: Path, loci: pd.DataFrame, cache: Path | None = None
) -> pd.DataFrame | None:
    """Mean off-diagonal r-squared per locus, averaged over sites.

    The covariate the divergence score turns out to be mostly measuring, so it has to be
    computed to say so. Reads the cached per-site tag-SNP r-squared matrices the locus
    selection stage already wrote; cached because it opens 3 files per locus.
    """
    if cache and Path(cache).exists():
        return pd.read_csv(cache, sep="\t")
    ld_dir = Path(data_root) / "loci" / "ld_matrices"
    if not ld_dir.exists() or "window_id" not in loci.columns:
        return None
    rows = []
    for r in loci.itertuples(index=False):
        vals = []
        for p in sorted(ld_dir.glob(f"{r.window_id}_*.npz")):
            try:
                with np.load(p, allow_pickle=True) as z:
                    m = np.asarray(z["r2"], dtype=float)
            except (OSError, ValueError, KeyError):
                continue
            vals.append(float(np.nanmean(m[np.triu_indices_from(m, k=1)])))
        if vals:
            rows.append((r.locus_id, float(np.mean(vals))))
    if not rows:
        return None
    out = pd.DataFrame(rows, columns=["locus_id", "mean_r2"])
    if cache:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(cache, sep="\t", index=False)
    return out


def pap7_ld_divergence(
    df: pd.DataFrame,
    out: Path,
    loci: pd.DataFrame | None = None,
    ld_strength: pd.DataFrame | None = None,
) -> Path:
    """Resolution against the LD-divergence stratum, and what that stratum actually is.

    THE STRATIFICATION IS CONFOUNDED, AND PANEL D IS THE POINT OF THE FIGURE.

    Read left to right, panels a-c look like a finding: power is flat but credible sets
    grow and posterior confidence falls as the cross-site LD-divergence stratum rises --
    the opposite of the premise that LD diversity sharpens fine-mapping.

    It is not a finding. The divergence score is the Frobenius norm of pairwise r-squared
    differences, and that quantity is bounded by the magnitudes it differences: a window
    where every site has r-squared near zero *cannot* score high. So the score cannot
    separate "the sites disagree about LD" from "there is a great deal of LD here", and
    on this package it is overwhelmingly the latter -- panel d measures the correlation.
    What a-c show is therefore the textbook result that more LD means bigger credible
    sets, arriving under a misleading axis label.

    A magnitude-free version of the same question -- one minus the mean correlation
    between the sites' r-squared vectors, which asks whether they RANK variant pairs the
    same way -- is near-independent of LD strength and shows no relationship with
    resolution at all. The honest statement is a null on LD disagreement, not a reversal.

    Neither is it the allele-coding defect: flip rate is flat across the strata, and
    removing every affected instance leaves the pattern intact.

    None of this indicts the method. Coverage holds near nominal in every stratum
    (pap2, panel c): returning wider credible sets where the data genuinely supports less
    resolution is correct behaviour, and a method that stayed sharp there would be
    miscalibrated.
    """
    strata = _strata_present(df)
    if not strata:
        fig, ax = plt.subplots(figsize=(7, 4))
        fs.no_data(ax, "results carry no stratum column")
        return fs.save(fig, out, log)
    rgs = sorted(df["rg"].dropna().unique())
    colors = dict(zip(rgs, fs.ordinal_colors(len(rgs)), strict=False))
    x = np.arange(len(strata))
    width = 0.8 / max(len(rgs), 1)

    fig, axes = plt.subplots(1, 4, figsize=(16.6, 4.8))
    ax_a, ax_b, ax_c, ax_d = axes

    for j, rg in enumerate(rgs):
        rows = [_rate(df[(df["stratum"] == s) & (df["rg"] == rg)], "captured") for s in strata]
        p, yerr = _err_bars([(r[0], r[1], r[2]) for r in rows])
        ax_a.bar(
            x + j * width - 0.4 + width / 2,
            p,
            width * 0.88,
            color=colors[rg],
            edgecolor=fs.SURFACE,
            linewidth=1.4,
            yerr=yerr,
            capsize=2,
            ecolor=fs.INK2,
            error_kw=dict(lw=1.0),
            zorder=3,
        )
    ax_a.set_ylim(0, 1.20)
    ax_a.set_ylabel("power   =   P(a causal variant is captured)")
    ax_a.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    fs.legend_swatches(
        ax_a,
        {f"$r_g$ = {r:g}": colors[r] for r in rgs},
        loc="upper center",
        ncol=len(rgs),
        bbox_to_anchor=(0.5, 1.01),
    )
    fs.panel_letter(ax_a, "a", "Power — flat")

    for ax, value, ylab, letter, title, logy in (
        (ax_b, "best_cs_size", "credible-set size (variants)", "b", "Resolution — falls", True),
        (
            ax_c,
            "causal_pip_max",
            "PIP of the true causal variant",
            "c",
            "Confidence — falls",
            False,
        ),
    ):
        data, positions, cols = [], [], []
        for i, s in enumerate(strata):
            for j, rg in enumerate(rgs):
                g = df[(df["stratum"] == s) & (df["rg"] == rg)]
                if value == "best_cs_size":
                    g = g[g["captured"]]
                positions.append(i + j * width - 0.4 + width / 2)
                data.append(g[value].dropna().values)
                cols.append(colors[rg])
        fs.boxplot(ax, data, positions=positions, colors=cols, widths=width * 0.82, points=False)
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(bottom=0.82)
        ax.set_ylabel(ylab)
        ax.set_xlim(-0.55, len(strata) - 0.45)
        fs.panel_letter(ax, letter, title)

    for ax in (ax_a, ax_b, ax_c):
        ax.set_xticks(x)
        ax.set_xticklabels([s.capitalize() for s in strata])
        ax.set_xlabel("cross-site LD-divergence stratum")
        ax.grid(axis="x", visible=False)

    # (d) what the stratum actually measures.
    merged = None
    if loci is not None and ld_strength is not None:
        merged = loci.merge(ld_strength, on="locus_id", how="inner")
    if merged is not None and len(merged) > 3:
        scol = dict(zip(strata, fs.ordinal_colors(len(strata), fs.VIOLET_RAMP), strict=False))
        for s in strata:
            sub = merged[merged["stratum"] == s]
            ax_d.scatter(
                sub["mean_r2"],
                sub["divergence_score"],
                s=40,
                color=scol[s],
                alpha=0.85,
                lw=0.7,
                edgecolor=fs.SURFACE,
                zorder=4,
                label=s,
            )
        r = float(np.corrcoef(merged["mean_r2"], merged["divergence_score"])[0, 1])
        ax_d.set_xlabel("mean $r^2$ in the window (LD strength)")
        ax_d.set_ylabel("cross-site LD-divergence score")
        ax_d.text(
            0.04,
            0.96,
            f"Pearson r = {r:+.3f}\n({len(merged)} loci)",
            transform=ax_d.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            color=fs.INK,
            fontweight="semibold",
            bbox=dict(boxstyle="round,pad=0.4", facecolor=fs.SURFACE, edgecolor=fs.GRID),
        )
        ax_d.legend(title="stratum", loc="lower right")
        fs.panel_letter(ax_d, "d", "…but the axis is LD strength")
    else:
        fs.no_data(ax_d, "needs the cached per-site r² matrices\n(pass --data-root)")
        fs.panel_letter(ax_d, "d", "…but the axis is LD strength")

    fig.suptitle(
        "Panels a–c are an LD-strength effect wearing a diversity label",
        x=0.005,
        ha="left",
        fontsize=12.5,
        fontweight="semibold",
    )
    fs.footnote(
        fig,
        f"{len(df):,} instances, all design cells. The divergence score is "
        "a Frobenius norm of r² differences and is bounded by the "
        "magnitudes it differences, so it cannot separate LD disagreement "
        "from LD abundance — panel d shows it tracking LD strength almost "
        "exactly. Panels a–c are therefore the expected result that more LD "
        "gives bigger credible sets. A magnitude-free measure of the same "
        "disagreement is near-independent of LD strength and shows no "
        "relationship with resolution, so the honest statement is a null. "
        "Coverage stays near nominal in every stratum (pap2c): wider sets "
        "where the data supports less resolution is correct behaviour.",
        y=-0.05,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap8 -- scalability   (paper Supp. Figs 8-9)
# --------------------------------------------------------------------------- #
def pap8_scalability(df: pd.DataFrame, out: Path, loci: pd.DataFrame | None = None) -> Path:
    """What one instance costs, and what drives the cost.

    The paper's scalability claim is about SuSiEx against PAINTOR and MsCAVIAR, which
    cannot run this problem size at all. With only SuSiEx in play the useful question
    changes: for a federated deployment the number that matters is what the coordinator
    spends per instance and how it scales with the locus, because that is what a
    scheduler allocation is sized against.
    """
    if "runtime_s" not in df.columns:
        fig, ax = plt.subplots(figsize=(7, 4))
        fs.no_data(ax, "results carry no runtime_s column")
        return fs.save(fig, out, log)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5))
    ax_a, ax_b, ax_c = axes
    rt = df["runtime_s"].dropna()

    # (a) the distribution, log x -- it is heavy-tailed and a linear axis shows one bar
    ax_a.hist(
        rt,
        bins=np.logspace(np.log10(max(rt.min(), 0.1)), np.log10(rt.max()), 50),
        color=fs.BLUE_RAMP[1],
        edgecolor=fs.SURFACE,
        linewidth=0.5,
    )
    ax_a.set_xscale("log")
    for q, style in ((0.5, "-"), (0.95, "--")):
        v = float(rt.quantile(q))
        ax_a.axvline(v, color=fs.INK, lw=1.3, ls=style, zorder=5)
        ax_a.text(
            v,
            ax_a.get_ylim()[1] * 0.97,
            f" p{int(q * 100)} = {v:.0f}s",
            rotation=90,
            va="top",
            ha="left",
            fontsize=8.5,
            color=fs.INK,
        )
    ax_a.set_xlabel("wall-clock per instance (s)")
    ax_a.set_ylabel("instances")
    fs.panel_letter(ax_a, "a", "Cost per instance")

    # (b) against the number of causal variants -- more single effects, more iterations
    ncsls = sorted(int(v) for v in df["ncsl"].dropna().unique())
    fs.boxplot(
        ax_b,
        [df[df["ncsl"] == n]["runtime_s"].dropna().values for n in ncsls],
        colors=fs.ordinal_colors(len(ncsls)),
        points=False,
    )
    ax_b.set_yscale("log")
    ax_b.set_xticks(range(1, len(ncsls) + 1))
    ax_b.set_xticklabels([str(n) for n in ncsls])
    ax_b.set_xlabel("causal variants simulated ($n_{csl}$)")
    ax_b.set_ylabel("wall-clock per instance (s)")
    ax_b.grid(axis="x", visible=False)
    fs.panel_letter(ax_b, "b", "Against model size")

    # (c) against the locus's variant count -- the quadratic term in the uplink and
    # in SuSiEx's own linear algebra
    if loci is not None and "n_variants_in_window" in loci.columns:
        merged = df.merge(loci[["locus_id", "n_variants_in_window"]], on="locus_id", how="left")
        g = (
            merged.groupby("locus_id")
            .agg(m=("n_variants_in_window", "first"), t=("runtime_s", "median"))
            .dropna()
        )
        ax_c.scatter(
            g["m"],
            g["t"],
            s=32,
            color=fs.BLUE_RAMP[2],
            alpha=0.75,
            lw=0.6,
            edgecolor=fs.SURFACE,
            zorder=4,
        )
        if len(g) > 3:
            r = float(np.corrcoef(g["m"], g["t"])[0, 1])
            # A quadratic reference anchored at the median locus. Not a fit: the point
            # is whether the observed growth is near M^2, and a fitted curve would hide
            # that by absorbing any exponent into its coefficients.
            mm = np.linspace(g["m"].min(), g["m"].max(), 100)
            anchor = float(g["t"].median()) / float(g["m"].median()) ** 2
            ax_c.plot(mm, anchor * mm**2, ls="--", lw=1.3, color=fs.INK2, zorder=3)
            fs.direct_label(ax_c, mm[-1], anchor * mm[-1] ** 2, " $\\propto M^2$", fs.INK2, dx=2)
            ax_c.text(
                0.035,
                0.96,
                f"Pearson r = {r:+.2f}   ({len(g)} loci)",
                transform=ax_c.transAxes,
                ha="left",
                va="top",
                fontsize=9,
                color=fs.INK,
            )
        ax_c.set_xlabel("variants in the locus window ($M$)")
        ax_c.set_ylabel("median wall-clock per instance (s)")
    else:
        fs.no_data(ax_c, "locus table not supplied")
    fs.panel_letter(ax_c, "c", "Against locus size")

    total = float(rt.sum())
    fs.footnote(
        fig,
        f"{len(rt):,} instances, {total / 3600:.0f} CPU-hours in total. "
        "Single-threaded; the sweep parallelises across instances.",
        y=-0.05,
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap9 -- coding quality against PIP   (paper Extended Data Figure 6)
# --------------------------------------------------------------------------- #
def pap9_coding_quality(df: pd.DataFrame, out: Path, flips: pd.DataFrame | None = None) -> Path:
    """Data-quality contamination binned by the confidence it produced.

    The paper bins variants by how much their PIP FELL between single- and
    cross-population fine-mapping, and shows that the biggest falls are enriched for
    variants with quality problems -- low-complexity regions, Hardy-Weinberg
    violations, multi-allelic sites. Its conclusion is that joint modelling across
    populations surfaces bad data.

    This package has its own quality defect of exactly that shape, and it is documented
    rather than hypothesised: the per-site filesets were cut without
    ``--keep-allele-order``, so about 7% of chr1 variants are coded against opposite
    alleles at different sites. Where that hits a causal variant, the sites' phenotypes
    carry opposite effect directions and the pooled signal is attenuated. This figure
    is that defect measured against the PIP it produced -- the same relationship the
    paper reports, on a contaminant whose identity is known exactly rather than
    inferred.
    """
    if flips is None or not len(flips):
        fig, ax = plt.subplots(figsize=(7.4, 4.6))
        fs.no_data(ax, "allele-coding flags not supplied\n(see load_flip_flags)")
        return fs.save(fig, out, log)

    d = df.merge(flips, on=["locus_id", "architecture_id", "replicate"], how="left")
    d["has_flip"] = d["has_flip"].fillna(False).astype(bool)
    d = d[d["causal_pip_max"].notna()]
    if d.empty:
        fig, ax = plt.subplots(figsize=(7.4, 4.6))
        fs.no_data(ax)
        return fs.save(fig, out, log)

    bins = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 1.0001]
    labels = ["0–0.1", "0.1–0.3", "0.3–0.5", "0.5–0.7", "0.7–0.9", "0.9–0.95", ">0.95"]
    d["bin"] = pd.cut(d["causal_pip_max"], bins=bins, labels=labels, include_lowest=True)

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.0, 4.7), gridspec_kw={"width_ratios": [1.45, 1]}
    )
    # (a) the paper's form: a stacked share per PIP bin
    share, ns = [], []
    for lab in labels:
        g = d[d["bin"] == lab]
        share.append(float(g["has_flip"].mean()) if len(g) else np.nan)
        ns.append(len(g))
    x = np.arange(len(labels))
    ax1.bar(
        x,
        [1 - (s if np.isfinite(s) else 0) for s in share],
        0.66,
        bottom=[s if np.isfinite(s) else 0 for s in share],
        color=fs.STATUS["good"],
        edgecolor=fs.SURFACE,
        linewidth=1.4,
    )
    ax1.bar(
        x,
        [s if np.isfinite(s) else 0 for s in share],
        0.66,
        color=fs.STATUS["critical"],
        edgecolor=fs.SURFACE,
        linewidth=1.4,
    )
    baseline = float(d["has_flip"].mean())
    ax1.axhline(baseline, color=fs.INK, lw=1.3, ls="--", zorder=6)
    ax1.text(
        -0.45,
        baseline,
        f"overall {baseline * 100:.1f}% ",
        ha="right",
        va="bottom",
        fontsize=8.5,
        color=fs.INK,
    )
    for xi, s in enumerate(share):
        if np.isfinite(s):
            ax1.text(
                xi,
                s / 2,
                f"{s * 100:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color=fs.SURFACE,
                fontweight="semibold",
            )

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=30, ha="right")
    fs.count_labels(ax1, x, ns, pad_pt=38)
    ax1.set_ylim(0, 1)
    ax1.set_xlabel("PIP reached by the true causal variant", labelpad=34)
    ax1.set_ylabel("share of instances")
    ax1.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    ax1.grid(axis="x", visible=False)
    # No legend box: the two categories are already named in the panel, and the red
    # segments carry their own percentage. A legend here has nowhere to sit that does
    # not cover either the bars or the baseline annotation.
    ax1.text(
        0.5,
        1.055,
        "red: at least one causal variant coded inconsistently across "
        "sites    ·    green: consistent",
        transform=ax1.transAxes,
        ha="center",
        va="bottom",
        fontsize=8.5,
        color=fs.INK2,
    )
    fs.title_with_letter(ax1, "a", "Coding contamination against confidence", dx=-0.11, dy=1.10)

    # (b) the effect stated directly, as the PIP distribution on each side
    data = [d[~d["has_flip"]]["causal_pip_max"].values, d[d["has_flip"]]["causal_pip_max"].values]
    fs.violin_box(ax2, data, colors=[fs.STATUS["good"], fs.STATUS["critical"]])
    ax2.set_xticks([1, 2])
    ax2.set_xticklabels(
        [f"consistent\n(n={len(data[0]):,})", f"inconsistent\n(n={len(data[1]):,})"]
    )
    ax2.set_ylabel("PIP of the true causal variant")
    ax2.set_ylim(-0.03, 1.05)
    ax2.grid(axis="x", visible=False)
    med = [float(np.median(v)) if len(v) else np.nan for v in data]
    ax2.text(
        0.5,
        0.02,
        f"median {med[0]:.3f} vs {med[1]:.3f}",
        transform=ax2.transAxes,
        ha="center",
        va="bottom",
        fontsize=9,
        color=fs.INK,
    )
    fs.panel_letter(ax2, "b", "Stated directly")

    fs.footnote(
        fig,
        f"{len(d):,} instances with a scored causal PIP. Inconsistent "
        "coding is detected by comparing the three sites' .bim A1 alleles "
        "at each causal variant. The defect attenuates pooled signal, so "
        "direction and size of any bias require a paired corrected rerun.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap10 -- population-specific causal probability   (paper Figure 4 / Discussion)
# --------------------------------------------------------------------------- #
def pap10_population_probability(
    cs_detail: pd.DataFrame, out: Path, cohort_n: dict[str, int] | None = None
) -> Path:
    """Is this signal causal in this population? SuSiEx's per-population answer.

    Once converged, SuSiEx estimates a population-specific causal probability for every
    credible set, and the paper thresholds it at 0.8 to call a signal causal in a given
    population. It is the only per-population statement the method makes, and it is the
    quantity most easily over-read: the paper is explicit that a low probability is
    "akin to a non-significant P-value from GWAS when the sample size is limited" and
    does not mean the variant is not causal there.

    Panel b is the check on exactly that. If the probability were measuring biology,
    it would vary with the ancestries' genetics; if it is measuring power, it will
    track cohort size. Plotting it against each column's n settles which.
    """
    pops = [c[len("prob_causal_") :] for c in cs_detail.columns if c.startswith("prob_causal_")]
    if not pops:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        fs.no_data(ax, "harvested table carries no POST-HOC_PROB columns")
        return fs.save(fig, out, log)
    # Ordered by cohort size when known, so panel a already carries the finding that
    # panel b then tests, rather than hiding it behind an alphabetical axis.
    if cohort_n:
        pops = sorted(pops, key=lambda p: -cohort_n.get(p, 0))

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.4, 4.8), gridspec_kw={"width_ratios": [1.5, 1]}
    )
    data = [cs_detail[f"prob_causal_{p}"].dropna().values for p in pops]
    fs.violin_box(ax1, data, colors=fs.ancestry_colors(pops))
    ax1.axhline(0.8, color=fs.INK, lw=1.3, ls="--", zorder=6)
    ax1.text(
        0.995,
        0.8,
        "0.8 — called causal in this population ",
        transform=ax1.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=8.5,
        color=fs.INK2,
    )
    labels = []
    for p in pops:
        n = f"\nn = {cohort_n[p]:,}" if cohort_n and p in cohort_n else ""
        labels.append(f"{p}{n}")
    ax1.set_xticks(range(1, len(pops) + 1))
    ax1.set_xticklabels(labels)
    for tick, p in zip(ax1.get_xticklabels(), pops, strict=False):
        tick.set_color(fs.ANCESTRY.get(p, fs.INK))
        tick.set_fontweight("semibold")
    ax1.set_ylim(-0.03, 1.05)
    ax1.set_ylabel("population-specific causal probability")
    ax1.grid(axis="x", visible=False)
    fs.count_labels(
        ax1,
        range(1, len(pops) + 1),
        [int((d > 0.8).sum()) for d in data],
        pad_pt=34,
        fmt=">0.8 in {:,}",
    )
    ax1.set_xlabel("SuSiEx population column", labelpad=32)
    fs.panel_letter(ax1, "a", "Per-population causal probability")

    if cohort_n:
        xs = [cohort_n.get(p, np.nan) for p in pops]
        ys = [float((d > 0.8).mean()) for d in data]
        ax2.scatter(
            xs, ys, s=110, color=fs.ancestry_colors(pops), lw=0.9, edgecolor=fs.SURFACE, zorder=5
        )
        for p, x, y in zip(pops, xs, ys, strict=False):
            fs.direct_label(ax2, x, y, f" {p}", fs.ANCESTRY.get(p, fs.INK))
        ok = [(x, y) for x, y in zip(xs, ys, strict=False) if np.isfinite(x)]
        if len(ok) > 2:
            lx = np.log10([x for x, _ in ok])
            r = float(np.corrcoef(lx, [y for _, y in ok])[0, 1])
            ax2.text(
                0.035,
                0.96,
                f"Pearson r = {r:+.2f}\n(against $\\log_{{10}}$ n)",
                transform=ax2.transAxes,
                ha="left",
                va="top",
                fontsize=9.5,
                color=fs.INK,
                bbox=dict(boxstyle="round,pad=0.35", facecolor=fs.SURFACE, edgecolor=fs.GRID),
            )
        ax2.set_xscale("log")
        ax2.set_xlabel("individuals in that ancestry column")
        ax2.set_ylabel("share of credible sets called causal (> 0.8)")
        ax2.set_ylim(-0.03, 1.05)
        ax2.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
        ax2.set_xlim(right=max(x for x in xs if np.isfinite(x)) * 2.4)
        fs.panel_letter(ax2, "b", "It tracks power, not biology")
    else:
        fs.no_data(ax2, "column sizes not supplied")

    fs.footnote(
        fig,
        f"{len(cs_detail):,} credible sets over "
        f"{cs_detail['instance'].nunique():,} instances. The causal "
        "architecture is stated in the figure group. Differences between populations "
        "can reflect sample size, frequencies and genetic effects; "
        "this association with N does not isolate a sample-size effect.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap11 -- the LocusZoom-style panel   (paper Figure 5, Ext. Data Figs 4-5)
# --------------------------------------------------------------------------- #
def pap11_locuszoom(
    work_dir: Path,
    out: Path,
    truth_snps: list[str] | None = None,
    title: str = "",
    pops: list[str] | None = None,
) -> Path | None:
    """One locus, every ancestry's association track, and the joint PIP beneath them.

    The paper's Figure 5 form, and the reason it uses it: the argument for
    cross-ancestry fine-mapping is that no single population's association plot picks
    the causal variant out, and the joint posterior does. That is an argument about
    several plots stacked on a shared genomic axis, and it does not survive being
    summarised into a table.

    Reads one instance's kept working directory: the per-ancestry ``.sumstats`` that
    went into SuSiEx, and the ``.snp`` PIPs that came out.
    """
    work_dir = Path(work_dir)
    snp_file = work_dir / "cs.snp"
    if not snp_file.exists():
        return None
    sumstats = sorted(work_dir.glob("*.sumstats"))
    if pops:
        order = {p: i for i, p in enumerate(pops)}
        sumstats = sorted(sumstats, key=lambda p: order.get(p.stem, 99))
    if not sumstats:
        return None
    truth = set(truth_snps or [])

    snp = pd.read_csv(snp_file, sep="\t")
    pip_cols = [c for c in snp.columns if c.startswith("PIP(")]
    if not pip_cols or "BP" not in snp.columns:
        return None
    from ..inference import inclusion_probabilities

    snp["pip"] = inclusion_probabilities(snp)

    n = len(sumstats)
    fig, axes = plt.subplots(
        n + 1,
        1,
        figsize=(9.6, 1.55 * n + 2.4),
        sharex=True,
        gridspec_kw={"height_ratios": [1] * n + [1.5]},
    )
    lead_bp = None
    for i, ss_path in enumerate(sumstats):
        ax = axes[i]
        pop = ss_path.stem
        ss = pd.read_csv(ss_path, sep="\t")
        if not {"bp", "p"}.issubset(ss.columns):
            fs.no_data(ax, f"{pop}: unreadable")
            continue
        logp = -np.log10(pd.to_numeric(ss["p"], errors="coerce").clip(lower=np.finfo(float).tiny))
        ax.scatter(
            ss["bp"], logp, s=5, color=fs.ANCESTRY.get(pop, fs.MUTED), alpha=0.45, lw=0, zorder=3
        )
        ax.axhline(-np.log10(5e-8), color=fs.MUTED, lw=1.0, ls="--", zorder=2)
        if "snp" in ss.columns and truth:
            hit = ss[ss["snp"].astype(str).isin(truth)]
            if len(hit):
                ax.scatter(
                    hit["bp"],
                    -np.log10(
                        pd.to_numeric(hit["p"], errors="coerce").clip(lower=np.finfo(float).tiny)
                    ),
                    s=64,
                    marker="D",
                    color=fs.STATUS["critical"],
                    zorder=6,
                    lw=0.9,
                    edgecolor=fs.SURFACE,
                )
                lead_bp = float(hit["bp"].iloc[0])
        # The ancestry's name, in its own colour, inside its own track: with six
        # stacked panels a shared legend would be six lookups per read.
        ax.text(
            0.006,
            0.93,
            f"{pop} GWAS",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9.5,
            fontweight="semibold",
            color=fs.ANCESTRY.get(pop, fs.INK),
        )
        ax.set_ylabel("$-\\log_{10}P$", fontsize=8.5)
        ax.set_ylim(bottom=0)
        ax.grid(axis="x", visible=False)
        ax.tick_params(labelsize=8)

    axp = axes[-1]
    axp.vlines(snp["BP"], 0, snp["pip"], color=fs.BLUE_RAMP[2], lw=0.8, zorder=3)
    top = snp.loc[snp["pip"].idxmax()] if len(snp) else None
    if top is not None:
        axp.scatter(
            [top["BP"]],
            [top["pip"]],
            s=90,
            marker="D",
            color=fs.STATUS["critical"],
            zorder=6,
            lw=0.9,
            edgecolor=fs.SURFACE,
        )
        name = str(top.get("SNP", ""))
        mark = "  ← the true causal variant" if name in truth else ""
        axp.annotate(
            f"{name}\nPIP = {top['pip']:.3f}{mark}",
            xy=(top["BP"], top["pip"]),
            xytext=(8, -6),
            textcoords="offset points",
            fontsize=9,
            color=fs.INK,
            ha="left",
            va="top",
            fontweight="semibold",
        )
    axp.set_ylim(0, 1.12)
    axp.set_ylabel("PIP\n(all ancestries jointly)", fontsize=8.5)
    axp.text(
        0.006,
        0.93,
        "SuSiEx",
        transform=axp.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        fontweight="semibold",
        color=fs.INK,
    )
    axp.grid(axis="x", visible=False)
    axp.tick_params(labelsize=8)
    axp.set_xlabel("position on chromosome 1 (Mb)")
    axp.xaxis.set_major_formatter(lambda v, _: f"{v / 1e6:.2f}")
    if lead_bp is not None:
        for ax in axes:
            ax.axvline(lead_bp, color=fs.STATUS["critical"], lw=0.9, ls=":", zorder=1, alpha=0.7)

    fig.suptitle(title or work_dir.name, x=0.005, ha="left", fontsize=12, fontweight="semibold")
    fs.footnote(
        fig,
        "Dashed horizontal line is genome-wide significance (5e-8); the "
        "red diamond is the simulated causal variant, and the dotted "
        "vertical line carries its position through every track. The "
        "bottom panel is one joint posterior over all ancestry columns, "
        "not a per-ancestry result stacked.",
        y=-0.015,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# pap12 -- estimated effect-size concordance   (paper Figure 4f / 4g)
# --------------------------------------------------------------------------- #
def pap12_effect_concordance(
    members: pd.DataFrame, out: Path, reference: str = "EUR", cohort_n: dict[str, int] | None = None
) -> Path:
    """Per-allele effect sizes of fine-mapped variants, one ancestry pair per panel.

    The paper's Figure 4f/4g, and its finding: effect sizes agree closely between the
    two well-powered populations (r = 0.79 for EUR vs EAS) and much less well against
    the small one (r = 0.21 for EUR vs AFR), and downsampling shows the difference is
    sample size rather than biology.

    This design can make that argument cleanly, because here the truth is known: the
    simulator drew the ancestry effects from one MVN with a specified r_g, so any
    departure from the identity line beyond that r_g is estimation error. Each panel
    prints the column's n, so the widening with smaller n can be read directly.
    """
    pops = [c[len("beta_") :] for c in members.columns if c.startswith("beta_")]
    others = [p for p in pops if p != reference]
    if reference not in pops or not others:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        fs.no_data(ax, "harvested members carry no per-population betas")
        return fs.save(fig, out, log)
    if cohort_n:
        others = sorted(others, key=lambda p: -cohort_n.get(p, 0))

    # Only confidently fine-mapped variants: a marginal effect at a variant with PIP
    # 0.02 is an estimate of a tag's effect, and comparing those across ancestries
    # measures LD rather than biology. The paper applies the same restriction.
    d = members[members["cs_pip"] > 0.5] if "cs_pip" in members.columns else members
    if len(d) < 10:
        d = members

    ncol = min(len(others), 5)
    fig, axes = plt.subplots(1, ncol, figsize=(2.9 * ncol + 0.6, 3.5))
    axes = np.atleast_1d(axes)
    for k, (pop, ax) in enumerate(zip(others[:ncol], axes, strict=False)):
        x = pd.to_numeric(d[f"beta_{reference}"], errors="coerce")
        y = pd.to_numeric(d[f"beta_{pop}"], errors="coerce")
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        if len(x) < 3:
            fs.no_data(ax, f"{pop}: too few")
            continue
        lim = float(max(np.abs(x).max(), np.abs(y).max())) * 1.12
        ax.plot([-lim, lim], [-lim, lim], ls="--", lw=1.2, color=fs.INK2, zorder=3)
        ax.axhline(0, color=fs.GRID, lw=0.8, zorder=1)
        ax.axvline(0, color=fs.GRID, lw=0.8, zorder=1)
        ax.scatter(
            x,
            y,
            s=22,
            color=fs.ANCESTRY.get(pop, fs.MUTED),
            alpha=0.55,
            lw=0.5,
            edgecolor=fs.SURFACE,
            zorder=4,
        )
        r = float(np.corrcoef(x, y)[0, 1])
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"$\\beta$  {reference}", fontsize=9)
        ax.set_ylabel(f"$\\beta$  {pop}", fontsize=9, color=fs.ANCESTRY.get(pop, fs.INK))
        sub = f"n = {cohort_n[pop]:,}" if cohort_n and pop in cohort_n else ""
        ax.text(
            0.04,
            0.96,
            f"r = {r:.2f}\n{sub}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            color=fs.INK,
        )
        ax.tick_params(labelsize=8)
        fs.panel_letter(ax, "abcdef"[k], f"{reference} vs {pop}", dx=-0.20, dy=1.04)

    fig.suptitle(
        "Estimated per-allele effect sizes of fine-mapped variants",
        x=0.005,
        ha="left",
        fontsize=12,
        fontweight="semibold",
    )
    fs.footnote(
        fig,
        f"{len(d):,} credible-set members with CS-PIP > 0.5. Effects are "
        "the marginal per-allele estimates SuSiEx reports per population. "
        "Panels are ordered by column size, so the widening left to right "
        "is the loss of precision as the cohort shrinks -- the causal "
        "effects themselves were drawn from one distribution.",
        y=-0.04,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# inputs that are not the results table
# --------------------------------------------------------------------------- #
def load_flip_flags(data_root: Path, cache: Path | None = None) -> pd.DataFrame | None:
    """Per-instance flag: is any of its causal variants coded inconsistently across sites?

    Recomputed from the site ``.bim`` files rather than read from a report, because the
    report of record is prose in ABOUT.md. Cached because it reads three 533k-line
    files, and every figure that wants it wants the same answer.
    """
    data_root = Path(data_root)
    if cache and Path(cache).exists():
        return pd.read_csv(cache, sep="\t")
    sites = (
        sorted(p.name for p in (data_root / "processed").iterdir() if p.is_dir())
        if (data_root / "processed").exists()
        else []
    )
    frames = {}
    for site in sites:
        bim = data_root / "processed" / site / f"{site}_chr1.bim"
        if bim.exists():
            frames[site] = pd.read_csv(
                bim, sep="\t", header=None, usecols=[1, 4], names=["snp", f"a1_{site}"]
            )
    if len(frames) < 2:
        return None
    merged = None
    for f in frames.values():
        merged = f if merged is None else merged.merge(f, on="snp")
    a1_cols = [c for c in merged.columns if c.startswith("a1_")]
    flipped = set(merged.loc[merged[a1_cols].nunique(axis=1) > 1, "snp"])

    man_path = data_root / "ground_truth" / "causal_manifest.tsv"
    if not man_path.exists():
        return None
    man = pd.read_csv(
        man_path, sep="\t", usecols=["locus_id", "architecture_id", "replicate", "causal_snp_ids"]
    )
    man["has_flip"] = (
        man["causal_snp_ids"].astype(str).map(lambda s: any(x in flipped for x in s.split(",")))
    )
    out = man[["locus_id", "architecture_id", "replicate", "has_flip"]]
    log.info(
        "allele-coding: %d of %d variants inconsistent; %d of %d instances affected",
        len(flipped),
        len(merged),
        int(out["has_flip"].sum()),
        len(out),
    )
    if cache:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(cache, sep="\t", index=False)
    return out


# Pooled cohort size per SuSiEx column in the published federation. Used only to
# ORDER and annotate figures, never to compute one; a mismatch would mislabel a panel
# but cannot change a statistic.
COHORT_N = {"AFR": 60000, "EUR": 36500, "MID": 25000, "CSA": 18500, "AMR": 7500, "EAS": 2500}


def load_cs_detail(detail_dir: Path, results: pd.DataFrame) -> pd.DataFrame | None:
    """The harvested per-credible-set table, joined to its instance's design cell."""
    p = Path(detail_dir) / "fm_credible_sets.tsv"
    if not p.exists():
        return None
    cs = pd.read_csv(p, sep="\t")
    if cs.empty:
        return None
    # instance -> (locus, architecture, replicate). Split on the LAST "_rep" so that an
    # architecture id containing underscores -- which they all do -- stays intact.
    keys = cs["instance"].astype(str).str.rsplit("_rep", n=1, expand=True)
    cs["replicate"] = pd.to_numeric(keys[1], errors="coerce").astype("Int64")
    head = keys[0].str.split("_", n=1, expand=True)
    cs["locus_id"] = head[0]
    cs["architecture_id"] = head[1]
    meta = results[
        ["locus_id", "architecture_id", "replicate", "stratum", "ncsl", "h2", "rg"]
    ].drop_duplicates()
    meta["replicate"] = meta["replicate"].astype("Int64")
    return cs.merge(
        meta, on=["locus_id", "architecture_id", "replicate"], how="inner", validate="many_to_one"
    )


# --------------------------------------------------------------------------- #
def write_figures(
    results, out_dir, data_root=None, detail_dir=None, loci=None, logger=None, strict=False
) -> list[Path]:
    """Draw every paper-analogue figure. A figure whose inputs are absent is skipped."""
    active = logger or log
    fs.apply_style()
    out_dir = Path(out_dir)
    df = parse_architecture(pd.DataFrame(results))
    if df.empty:
        active.warning("no results to plot")
        return []

    instance_ids = {
        f"{r.locus_id}_{r.architecture_id}_rep{r.replicate}" for r in df.itertuples(index=False)
    }
    flips = cs_detail = ld_strength = None
    exemplar_root = Path(detail_dir).parent / "exemplars" / "work" if detail_dir else None
    if data_root:
        try:
            flips = load_flip_flags(
                Path(data_root), Path(out_dir).parent / "cache" / "allele_flips.tsv"
            )
        except Exception as exc:  # noqa: BLE001
            active.warning("could not derive allele-coding flags: %s", exc)
        if loci is None:
            lp = Path(data_root) / "loci" / "selected_loci.tsv"
            if lp.exists():
                loci = pd.read_csv(lp, sep="\t")
        if loci is not None:
            try:
                ld_strength = locus_ld_strength(
                    Path(data_root), loci, Path(out_dir).parent / "cache" / "locus_ld_strength.tsv"
                )
            except Exception as exc:  # noqa: BLE001
                active.warning("could not measure locus LD strength: %s", exc)
    members = None
    if detail_dir:
        try:
            cs_detail = load_cs_detail(Path(detail_dir), df)
            if cs_detail is not None:
                active.info(
                    "per-credible-set detail: %d sets over %d instances",
                    len(cs_detail),
                    cs_detail["instance"].nunique(),
                )
        except Exception as exc:  # noqa: BLE001
            active.warning("could not load per-credible-set detail: %s", exc)
        mp = Path(detail_dir) / "fm_cs_members.tsv"
        if mp.exists():
            members = pd.read_csv(mp, sep="\t")
            instance_ids = {
                f"{r.locus_id}_{r.architecture_id}_rep{r.replicate}"
                for r in df.itertuples(index=False)
            }
            members = members[members.instance.isin(instance_ids)]
            active.info("credible-set members: %d rows", len(members))

    jobs = [
        ("pap1_recall_grid.png", lambda o: pap1_recall_grid(df, o)),
        ("pap2_diversity_panels.png", lambda o: pap2_diversity_panels(df, o, cs_detail)),
        ("pap3_convergence.png", lambda o: pap3_convergence(df, o)),
        ("pap4_power_coverage.png", lambda o: pap4_power_coverage(df, o, cs_detail)),
        ("pap5_resolution.png", lambda o: pap5_resolution(df, o)),
        ("pap6_high_confidence.png", lambda o: pap6_high_confidence(df, o)),
        ("pap7_ld_divergence.png", lambda o: pap7_ld_divergence(df, o, loci, ld_strength)),
        ("pap8_scalability.png", lambda o: pap8_scalability(df, o, loci)),
        ("pap9_coding_quality.png", lambda o: pap9_coding_quality(df, o, flips)),
        (
            "pap10_population_probability.png",
            lambda o: (
                pap10_population_probability(cs_detail, o, COHORT_N)
                if cs_detail is not None and len(cs_detail)
                else None
            ),
        ),
        (
            "pap12_effect_concordance.png",
            lambda o: (
                pap12_effect_concordance(members, o, "EUR", COHORT_N)
                if members is not None and len(members)
                else None
            ),
        ),
    ]
    # One LocusZoom-style panel per kept exemplar. Discovered rather than configured:
    # which instances were kept is a property of the detail run, not of this call.
    truth_by_instance = {}
    if data_root:
        mp = Path(data_root) / "ground_truth" / "causal_manifest.tsv"
        if mp.exists():
            man = pd.read_csv(
                mp, sep="\t", usecols=["locus_id", "architecture_id", "replicate", "causal_snp_ids"]
            )
            truth_by_instance = {
                f"{r.locus_id}_{r.architecture_id}_rep{int(r.replicate)}": str(
                    r.causal_snp_ids
                ).split(",")
                for r in man.itertuples(index=False)
            }
    for root in [
        r for r in (exemplar_root, Path(detail_dir) / "work" if detail_dir else None) if r
    ]:
        if not Path(root).exists():
            continue
        for d in sorted(Path(root).iterdir()):
            if not d.is_dir() or d.name in {"keeps"} or d.name.endswith("_refs"):
                continue
            if d.name not in instance_ids:
                continue
            if not (d / "cs.snp").exists():
                continue
            jobs.append(
                (
                    f"pap11_locuszoom_{d.name}.png",
                    lambda o, d=d: pap11_locuszoom(
                        d,
                        o,
                        truth_by_instance.get(d.name),
                        title=f"{d.name}   —   one locus, six ancestries, one posterior",
                        pops=list(COHORT_N),
                    ),
                )
            )
        break  # the first root that has exemplars wins; do not draw both

    written = []
    for name, draw in jobs:
        try:
            p = draw(out_dir / name)
        except Exception as exc:  # noqa: BLE001 - a figure must never lose a run
            if strict:
                raise RuntimeError(f"Could not draw {name}") from exc
            active.warning("could not draw %s: %s", name, exc)
            continue
        if p is None:
            active.info("skipped %s (inputs absent)", name)
        else:
            written.append(p)
    active.info("wrote %d paper-analogue figure(s) -> %s", len(written), out_dir)
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Render the SuSiEx-paper analogue figures.")
    ap.add_argument("--results", required=True, help="per-instance results TSV")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument(
        "--data-root",
        default=None,
        help="simulated package, for the locus table and coding-quality flags",
    )
    ap.add_argument(
        "--detail-dir",
        default=None,
        help="output of figures.detail, for true per-credible-set coverage",
    )
    args = ap.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )
    p = Path(args.results)
    if not p.exists():
        print(f"ERROR: {p} not found", file=sys.stderr)
        return 2
    df = pd.read_csv(p, sep="\t")
    log.info("loaded %d instance(s), %d architecture(s)", len(df), df["architecture_id"].nunique())
    write_figures(df, Path(args.out_dir), args.data_root, args.detail_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
