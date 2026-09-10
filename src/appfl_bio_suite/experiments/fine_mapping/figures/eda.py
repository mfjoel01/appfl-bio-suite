"""Exploratory figures over the simulated package, BEFORE any fine-mapping runs.

COORDINATOR-SIDE. These describe the data the experiment is about to be run on, which
is a different job from describing its results, and they are worth drawing first: every
claim the results figures make is conditional on this package being what it says it is.

Each figure answers a question a reader will otherwise ask in review:

  eda1_site_composition        who is enrolled where, and how unlike each other are
                               the three sites really?
  eda2_locus_selection         how was cross-site LD divergence distributed over the
                               candidate windows, and where did the chosen loci land?
  eda3_window_size             does the divergence score just track window density?
                               (if it did, the stratification would be confounded)
  eda4_ld_divergence_maps      what does "high divergence" actually look like as LD?
  eda5_causal_maf              are the causal variants at frequencies all three sites
                               can actually see?
  eda6_h2_calibration          did the simulator hit the heritability it was asked for?
  eda7_rg_recovery             did it hit the cross-ancestry genetic correlation?
  eda8_design_grid             which cells of the ncsl x h2 x rg star design exist,
                               and how many instances landed in each?

They read the simulated package directly rather than a summary of it, because a summary
is exactly the thing under test here.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from appfl_bio_suite.experiments.fine_mapping.figures import figstyle as fs
from appfl_bio_suite.experiments.fine_mapping.figures.figstyle import plt

log = logging.getLogger(__name__)

# The site order is fixed by declared dominance, not by size: it is the axis the whole
# federation argument runs along (EUR-dominant / AFR-dominant / MID-dominant).
SITE_ORDER = ["anl", "covenant", "mbzuai"]
SITE_DOMINANT = {"anl": "EUR", "covenant": "AFR", "mbzuai": "MID"}
ANCESTRY_ORDER = ["EUR", "AFR", "AMR", "EAS", "CSA", "MID"]
STRATUM_ORDER = ["low", "medium", "high"]


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #
def load_package(root: Path) -> dict:
    """Read every table the EDA figures need. Missing pieces are None, not fatal.

    A partial package still produces the figures it can support. That matters because
    these are also the figures someone draws while a long simulation is still running.
    """
    root = Path(root)
    out: dict = {"root": root}

    def _tsv(rel):
        p = root / rel
        return pd.read_csv(p, sep="\t") if p.exists() else None

    out["candidates"] = _tsv("loci/candidate_windows.tsv")
    out["loci"] = _tsv("loci/selected_loci.tsv")
    out["manifest"] = _tsv("ground_truth/causal_manifest.tsv")

    comp = {}
    for site in SITE_ORDER:
        p = root / "processed" / site / f"{site}_manifest.tsv"
        if p.exists():
            m = pd.read_csv(p, sep="\t", usecols=["superpopulation"])
            comp[site] = m["superpopulation"].value_counts().to_dict()
    out["composition"] = comp or None
    out["ld_dir"] = root / "loci" / "ld_matrices"
    out["effects_dir"] = root / "ground_truth" / "effect_sizes"
    return out


# --------------------------------------------------------------------------- #
# eda1 -- who is enrolled where
# --------------------------------------------------------------------------- #
def eda1_site_composition(pkg: dict, out: Path) -> Path | None:
    comp = pkg.get("composition")
    if not comp:
        return None
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(11.4, 4.4), gridspec_kw={"width_ratios": [1.25, 1]}
    )

    # (a) absolute headcount, stacked by ancestry. Stacked because the question is
    # "how is one site's 50,000 divided", and a stack answers that in one bar.
    sites = [s for s in SITE_ORDER if s in comp]
    bottoms = np.zeros(len(sites))
    for anc in ANCESTRY_ORDER:
        vals = np.array([comp[s].get(anc, 0) for s in sites], dtype=float)
        if not vals.any():
            continue
        ax1.bar(
            range(len(sites)),
            vals,
            bottom=bottoms,
            width=0.62,
            color=fs.ANCESTRY[anc],
            edgecolor=fs.SURFACE,
            linewidth=1.4,
        )
        # Direct labels inside the segment carry identity where the segment is big
        # enough to hold text; the legend covers the rest. This is the relief rule
        # for the three sub-3:1 ancestry slots.
        for i, (v, b) in enumerate(zip(vals, bottoms, strict=False)):
            if v >= 4500:
                ax1.text(
                    i,
                    b + v / 2,
                    f"{anc}\n{int(v):,}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=fs.SURFACE,
                    fontweight="semibold",
                )
        bottoms += vals
    ax1.set_xticks(range(len(sites)))
    ax1.set_xticklabels([f"{s}\n({SITE_DOMINANT.get(s, '?')}-dominant)" for s in sites])
    ax1.set_ylabel("individuals enrolled")
    ax1.yaxis.set_major_formatter(lambda v, _: f"{v / 1000:g}k")
    for i in range(len(sites)):
        ax1.text(
            i,
            bottoms[i] + bottoms.max() * 0.015,
            f"n = {int(bottoms[i]):,}",
            ha="center",
            va="bottom",
            fontsize=9,
            color=fs.INK2,
        )
    ax1.set_ylim(0, bottoms.max() * 1.10)
    fs.panel_letter(ax1, "a", "Site cohorts, by ancestry")
    fs.legend_swatches(
        ax1,
        {a: fs.ANCESTRY[a] for a in ANCESTRY_ORDER if any(comp[s].get(a) for s in sites)},
        loc="upper left",
        bbox_to_anchor=(1.005, 1.0),
        title="ancestry",
    )

    # (b) the same numbers as within-site SHARE. The absolute panel hides the point --
    # every site holds 50,000 -- and the point is that the three shares barely overlap.
    for i, s in enumerate(sites):
        total = sum(comp[s].values())
        left = 0.0
        for anc in ANCESTRY_ORDER:
            v = comp[s].get(anc, 0) / total
            if v <= 0:
                continue
            ax2.barh(
                i,
                v,
                left=left,
                height=0.55,
                color=fs.ANCESTRY[anc],
                edgecolor=fs.SURFACE,
                linewidth=1.4,
            )
            if v >= 0.10:
                ax2.text(
                    left + v / 2,
                    i,
                    f"{v * 100:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=fs.SURFACE,
                    fontweight="semibold",
                )
            left += v
    ax2.set_yticks(range(len(sites)))
    ax2.set_yticklabels(sites)
    ax2.set_xlim(0, 1)
    ax2.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax2.xaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    ax2.set_xlabel("share of the site's cohort")
    ax2.invert_yaxis()
    ax2.grid(axis="y", visible=False)
    fs.panel_letter(ax2, "b", "The same cohorts as shares")

    fs.footnote(
        fig,
        "Every individual appears at exactly one site; disjointness is "
        "asserted before extraction and again after bundling.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# eda2 -- how the loci were chosen
# --------------------------------------------------------------------------- #
def eda2_locus_selection(pkg: dict, out: Path) -> Path | None:
    cand, loci = pkg.get("candidates"), pkg.get("loci")
    if cand is None or loci is None:
        return None
    ok = cand[cand.get("status", "ok") == "ok"] if "status" in cand else cand
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9.2, 6.0), sharex=True, gridspec_kw={"height_ratios": [2.4, 1]}
    )

    # (a) the candidate pool as a histogram, with the selected loci laid over it.
    # One hue: this is a magnitude distribution, not three categories competing.
    ax1.hist(
        ok["divergence_score"],
        bins=44,
        color=fs.BLUE_RAMP[0],
        alpha=0.9,
        edgecolor=fs.SURFACE,
        linewidth=0.6,
        label=f"candidate windows (n={len(ok):,})",
    )
    strata = [s for s in STRATUM_ORDER if s in set(loci["stratum"])]
    colors = dict(zip(strata, fs.ordinal_colors(len(strata), fs.VIOLET_RAMP), strict=False))
    # Headroom first, then the rug: placing the rug at a fraction of the CURRENT ylim
    # and then extending it is what leaves the marks floating mid-panel.
    ymax = ax1.get_ylim()[1] * 1.16
    ax1.set_ylim(0, ymax)
    for s in strata:  # noqa: B007 - s indexes colors[s] below
        sub = loci[loci["stratum"] == s]
        ax1.plot(
            sub["divergence_score"],
            np.full(len(sub), ymax * 0.965),
            "|",
            color=colors[s],
            ms=11,
            mew=1.6,
            zorder=5,
        )
    ax1.set_ylabel("candidate windows")
    ax1.text(
        0.008,
        0.985,
        "selected loci  │││",
        transform=ax1.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color=fs.INK2,
    )
    fs.panel_letter(ax1, "a", "Cross-site LD divergence over the candidate pool")
    fs.legend_swatches(
        ax1, colors, title="selected stratum", loc="upper right", bbox_to_anchor=(1.0, 0.90), ncol=1
    )

    # (b) the selected loci alone, as a strip per stratum -- shows the bands actually
    # separate, which the overlay above cannot because the histogram dominates it.
    for i, s in enumerate(strata):
        sub = loci[loci["stratum"] == s]
        rng = np.random.default_rng(i)
        jitter = np.full(len(sub), i) + rng.uniform(-0.16, 0.16, len(sub))
        ax2.scatter(
            sub["divergence_score"],
            jitter,
            s=26,
            color=colors[s],
            alpha=0.85,
            lw=0.6,
            edgecolor=fs.SURFACE,
        )
        fs.direct_label(
            ax2, sub["divergence_score"].max(), i, f"{s}  (n={len(sub)})", colors[s], dx=10
        )
    ax2.set_yticks(range(len(strata)))
    ax2.set_yticklabels([])
    ax2.set_ylim(-0.6, len(strata) - 0.15)
    ax2.grid(axis="y", visible=False)
    ax2.set_xlabel(
        "cross-site LD divergence score   (Frobenius norm of pairwise $r^2$ differences)"
    )
    fs.panel_letter(ax2, "b", "")

    fs.footnote(
        fig,
        f"{len(loci)} loci retained of {len(ok):,} scored windows; the rest "
        "dropped greedily by the non-overlap constraint. Divergence is the "
        "experiment's independent variable, so loci are chosen to span it "
        "rather than sampled at random.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# eda3 -- is the score confounded with window density?
# --------------------------------------------------------------------------- #
def eda3_window_size(pkg: dict, out: Path) -> Path | None:
    cand, loci = pkg.get("candidates"), pkg.get("loci")
    if cand is None:
        return None
    ok = cand[cand.get("status", "ok") == "ok"] if "status" in cand else cand
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.scatter(
        ok["n_variants_in_window"],
        ok["divergence_score"],
        s=14,
        color=fs.MUTED,
        alpha=0.45,
        lw=0,
        label="candidate window",
    )
    if loci is not None:
        strata = [s for s in STRATUM_ORDER if s in set(loci["stratum"])]
        colors = dict(zip(strata, fs.ordinal_colors(len(strata), fs.VIOLET_RAMP), strict=False))
        for s in strata:
            sub = loci[loci["stratum"] == s]
            ax.scatter(
                sub["n_variants_in_window"],
                sub["divergence_score"],
                s=34,
                color=colors[s],
                lw=0.7,
                edgecolor=fs.SURFACE,
                zorder=4,
                label=s,
            )
    r = np.corrcoef(ok["n_variants_in_window"], ok["divergence_score"])[0, 1]
    ax.set_xlabel("variants in the 1 Mb window")
    ax.set_ylabel("cross-site LD divergence score")
    ax.set_title("Divergence is not a restatement of window density", loc="left")
    # The correlation IS the finding here, so it is stated in the panel rather than
    # left for the reader to estimate from the cloud.
    ax.text(
        0.985,
        0.03,
        f"Pearson r = {r:+.2f}   (n = {len(ok):,})",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
        color=fs.INK,
        bbox=dict(boxstyle="round,pad=0.4", facecolor=fs.SURFACE, edgecolor=fs.GRID),
    )
    ax.legend(loc="upper left", ncol=2)
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# eda4 -- what divergence looks like
# --------------------------------------------------------------------------- #
def eda4_ld_divergence_maps(pkg: dict, out: Path, n_show: int = 120) -> Path | None:
    """LD heatmaps for the same window at all three sites, for a high- and a
    low-divergence locus, plus the pairwise difference that the score integrates.

    Diverging colours for the difference panels only. A signed quantity gets a neutral
    grey at zero so "the sites agree here" reads as nothing rather than as a value.
    """
    loci, ld_dir = pkg.get("loci"), pkg.get("ld_dir")
    if loci is None or not Path(ld_dir).exists():
        return None

    def _pick(stratum):
        sub = loci[loci["stratum"] == stratum]
        if not len(sub):
            return None
        row = sub.loc[
            sub["divergence_score"].idxmax()
            if stratum == "high"
            else sub["divergence_score"].idxmin()
        ]
        return row

    picks = [(s, _pick(s)) for s in ("high", "low")]
    picks = [(s, r) for s, r in picks if r is not None]
    if not picks:
        return None

    import matplotlib.colors as mcolors

    seq = mcolors.LinearSegmentedColormap.from_list("r2", ["#fcfcfb", *fs.BLUE_RAMP])
    div = mcolors.LinearSegmentedColormap.from_list(
        "d", [fs.DIVERGING[0], fs.DIVERGING[1], fs.DIVERGING[2]]
    )

    # Load every panel BEFORE drawing any, so the difference panels can share one
    # symmetric colour limit. Rescaling each row to its own maximum -- the obvious way
    # to write this -- would make the low-divergence difference panel as saturated as
    # the high-divergence one, and the whole figure exists to show that it is not.
    loaded = []
    for stratum, row in picks:
        mats = {}
        for site in SITE_ORDER:
            p = Path(ld_dir) / f"{row['window_id']}_{site}.npz"
            if p.exists():
                with np.load(p, allow_pickle=True) as z:
                    mats[site] = np.asarray(z["r2"], dtype=float)[:n_show, :n_show]
        loaded.append((stratum, row, mats))
    dlim = 0.0
    for _, _, mats in loaded:
        pair = [x for x in SITE_ORDER if x in mats][:2]
        if len(pair) == 2:
            dlim = max(dlim, float(np.abs(mats[pair[0]] - mats[pair[1]]).max()))
    dlim = dlim or 1.0

    fig, axes = plt.subplots(len(picks), 4, figsize=(13.4, 3.35 * len(picks)))
    axes = np.atleast_2d(axes)
    for r_i, (stratum, row, mats) in enumerate(loaded):
        if len(mats) < 2:
            for ax in axes[r_i]:
                fs.no_data(ax)
            continue
        for c_i, site in enumerate(SITE_ORDER):
            ax = axes[r_i][c_i]
            if site not in mats:
                fs.no_data(ax, f"{site}: absent")
                continue
            im = ax.imshow(mats[site], cmap=seq, vmin=0, vmax=1, interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            ax.set_title(site, loc="left", fontsize=10, color=fs.SITE[site])
            if c_i == 0:
                ax.set_ylabel(
                    f"{stratum} divergence\n{row['locus_id']}  "
                    f"(score {row['divergence_score']:.1f})",
                    fontsize=9,
                )
        # the difference panel: the two most unlike sites at this locus
        ax = axes[r_i][3]
        pair = [x for x in SITE_ORDER if x in mats][:2]
        d = mats[pair[0]] - mats[pair[1]]
        imd = ax.imshow(d, cmap=div, vmin=-dlim, vmax=dlim, interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        ax.set_title(f"{pair[0]} − {pair[1]}", loc="left", fontsize=10)
        # The mean absolute difference is the panel's one number; without it a reader
        # has to judge saturation by eye across two rows.
        ax.text(
            0.5,
            -0.045,
            f"mean |$\\Delta r^2$| = {np.abs(d).mean():.3f}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8.5,
            color=fs.INK,
        )
        cb = fig.colorbar(imd, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("$\\Delta r^2$", fontsize=8)
        cb.ax.tick_params(labelsize=7)
        cbs = fig.colorbar(im, ax=axes[r_i][2], fraction=0.046, pad=0.03)
        cbs.set_label("$r^2$", fontsize=8)
        cbs.ax.tick_params(labelsize=7)

    fig.suptitle(
        "Linkage disequilibrium at the same window, seen from three sites",
        x=0.005,
        ha="left",
        fontsize=12,
        fontweight="semibold",
    )
    fs.footnote(
        fig,
        f"First {n_show} tag SNPs of each window, ordered by position. "
        "Both difference panels share one colour scale "
        f"(±{dlim:.2f}), so the low-divergence row is paler because the "
        "sites agree there, not because it was rescaled. The divergence "
        "score integrates that panel over all site pairs and all tag SNPs.",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# eda5 -- are the causal variants visible at every site?
# --------------------------------------------------------------------------- #
def eda5_causal_maf(pkg: dict, out: Path) -> Path | None:
    man = pkg.get("manifest")
    if man is None or "per_site_mafs_json" not in man.columns:
        return None
    rows = []
    for raw in man["per_site_mafs_json"].dropna():
        try:
            d = json.loads(raw)
        except (TypeError, ValueError):
            continue
        for site, maf in d.items():
            if maf is not None and np.isfinite(maf):
                rows.append((site, float(maf)))
    if not rows:
        return None
    maf = pd.DataFrame(rows, columns=["site", "maf"])

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(11.2, 4.4), gridspec_kw={"width_ratios": [1.5, 1]}
    )
    # (a) the spectra, one step-outline per site so three overlapping distributions
    # stay readable -- filled histograms at this overlap would hide two of the three.
    bins = np.linspace(0, 0.5, 51)
    for site in SITE_ORDER:
        v = maf.loc[maf["site"] == site, "maf"]
        if not len(v):
            continue
        ax1.hist(v, bins=bins, histtype="step", lw=2.0, color=fs.SITE[site], label=site)
        fs.direct_label(ax1, float(np.median(v)), 0, "", fs.SITE[site])
    ax1.set_xlabel("minor allele frequency of the causal variant, at that site")
    ax1.set_ylabel("causal variants")
    ax1.legend(title="site")
    fs.panel_letter(ax1, "a", "Causal-variant frequency spectra")

    # (b) the operational question: how often is a causal variant rare enough at one
    # site to be filtered out there? The MAF floor the pipeline uses is 0.005.
    thresholds = [0.005, 0.01, 0.05]
    x = np.arange(len(thresholds))
    width = 0.26
    for i, site in enumerate(SITE_ORDER):
        v = maf.loc[maf["site"] == site, "maf"]
        if not len(v):
            continue
        frac = [float((v < t).mean()) for t in thresholds]
        ax2.bar(
            x + (i - 1) * width,
            frac,
            width * 0.92,
            color=fs.SITE[site],
            edgecolor=fs.SURFACE,
            linewidth=1.4,
            label=site,
        )
        for xi, f in zip(x + (i - 1) * width, frac, strict=False):
            ax2.text(xi, f, f"{f * 100:.1f}%", ha="center", va="bottom", fontsize=8, color=fs.INK2)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"MAF < {t:g}" for t in thresholds])
    ax2.set_ylabel("share of causal variants")
    ax2.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    fs.panel_letter(ax2, "b", "Below the frequency floor")

    fs.footnote(
        fig,
        f"{len(man):,} instances x sites. A causal variant below a site's "
        "MAF floor is dropped from that site's moments, so the ancestry "
        "column it belongs to loses it.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# eda6 -- did the simulator hit the heritability it was asked for?
# --------------------------------------------------------------------------- #
def eda6_h2_calibration(pkg: dict, out: Path) -> Path | None:
    man = pkg.get("manifest")
    cols = [c for c in (man.columns if man is not None else []) if c.startswith("empirical_h2_")]
    if man is None or not cols:
        return None
    long = man.melt(id_vars=["h2_target"], value_vars=cols, var_name="site", value_name="h2_emp")
    long["site"] = long["site"].str.replace("empirical_h2_", "", regex=False)
    long = long[np.isfinite(long["h2_emp"])]
    long["rel_err"] = (long["h2_emp"] - long["h2_target"]) / long["h2_target"]
    targets = sorted(long["h2_target"].unique())
    sites = [s for s in SITE_ORDER if s in set(long["site"])]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.5))

    # (a) realised against target. Log y because the targets span an order of
    # magnitude; a linear axis puts the two smaller ones on top of each other.
    positions, data, colors = [], [], []
    width = 0.26
    for i, t in enumerate(targets):
        for j, site in enumerate(sites):
            v = long[(long["h2_target"] == t) & (long["site"] == site)]["h2_emp"].values
            if len(v):
                positions.append(i + (j - 1) * width)
                data.append(v)
                colors.append(fs.SITE[site])
    # Wider boxes and a hairline edge: at 0.22 units the 1.4pt surface-coloured edge
    # of the house boxplot covers the fill entirely and every box reads grey.
    bp = ax1.boxplot(
        data,
        positions=positions,
        widths=width * 0.80,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color=fs.SURFACE, lw=1.2),
        whiskerprops=dict(color=fs.MUTED, lw=0.9),
        capprops=dict(color=fs.MUTED, lw=0.9),
        boxprops=dict(edgecolor="none"),
    )
    for patch, c in zip(bp["boxes"], colors, strict=False):
        patch.set_facecolor(c)
    for i, t in enumerate(targets):
        ax1.plot([i - 0.46, i + 0.46], [t, t], color=fs.INK, lw=1.4, ls="--", zorder=6)
    ax1.set_yscale("log")
    ax1.set_xticks(range(len(targets)))
    ax1.set_xticklabels([f"{t:g}" for t in targets])
    ax1.set_xlabel("target per-locus $h^2$")
    ax1.set_ylabel("realised per-locus $h^2$")
    ax1.set_ylim(min(targets) * 0.55, max(targets) * 1.9)
    fs.legend_swatches(
        ax1,
        {s: fs.SITE[s] for s in sites},
        title="site",
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        ncol=3,
    )
    ax1.text(
        0.99,
        0.03,
        "dashed = target",
        transform=ax1.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        color=fs.INK2,
    )
    fs.panel_letter(ax1, "a", "Realised against target heritability")

    # (b) the error, on the scale it actually occupies. The validation harness's +/-20%
    # tolerance is an order of magnitude outside this range, so drawing it would flatten
    # the distribution into a line at zero and say nothing; the number goes in the text
    # instead, where it can be read.
    data = [long[long["site"] == s]["rel_err"].values * 100 for s in sites]
    fs.boxplot(ax2, data, colors=[fs.SITE[s] for s in sites], points=True, point_alpha=0.05)
    ax2.axhline(0, color=fs.INK, lw=1.2, ls="--", zorder=5)
    ax2.set_xticks(range(1, len(sites) + 1))
    ax2.set_xticklabels(sites)
    ax2.set_ylabel("relative error   $(h^2_{obs} - h^2_{target}) / h^2_{target}$   (%)")
    lim = float(np.percentile(np.abs(np.concatenate(data)), 99.9)) * 1.25
    ax2.set_ylim(-lim, lim)
    worst = float(np.abs(np.concatenate(data)).max())
    for i, (site, v) in enumerate(zip(sites, data, strict=False), start=1):
        fs.direct_label(
            ax2,
            i,
            np.percentile(np.abs(v), 95),
            f"p95 |err| {np.percentile(np.abs(v), 95):.2f}%",
            fs.SITE[site],
            dx=0,
            dy=9,
            ha="center",
        )
    ax2.text(
        0.99,
        0.02,
        f"worst single instance: {worst:.2f}%\nvalidation tolerance: ±20%",
        transform=ax2.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        color=fs.INK2,
    )
    fs.panel_letter(ax2, "b", "Calibration error, on its own scale")

    fs.footnote(
        fig,
        f"{len(man):,} instances x {len(sites)} sites. Realised $h^2$ is "
        "computed from the simulated genetic values at each site, so it "
        "carries that site's own allele frequencies.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# eda7 -- did it hit the genetic correlation?
# --------------------------------------------------------------------------- #
def eda7_rg_recovery(pkg: dict, out: Path, max_files: int = 4000) -> Path | None:
    """Cross-ancestry effect-size correlation, recovered from the drawn betas.

    THE ESTIMATOR MATTERS AND THE OBVIOUS ONE IS WRONG. Correlating a single instance's
    ancestry columns uses only its ncsl causal variants -- two or three points -- and a
    correlation over three points is almost always +/-1 regardless of the truth. The
    betas are drawn i.i.d. per causal variant from one MVN, so the correlation is a
    property of the DRAW, not of an instance: pooling every causal variant at a given
    target rg into one (N x n_ancestries) matrix and taking that matrix's off-diagonals
    is the estimator with the sample size to say anything. This is also what the
    validation harness does, which is why its numbers and these agree.

    Panels b-c are the SuSiEx paper's Figure 4f/4g form: one ancestry pair per panel
    with a y = x reference. Faceting by pair rather than colouring six ancestries in one
    scatter is not a style preference -- the palette validator reports the six-ancestry
    all-pairs separation as a hard failure, so those hues cannot carry identity here.
    """
    man, edir = pkg.get("manifest"), Path(pkg.get("effects_dir", ""))
    if man is None or not edir.exists():
        return None
    rng = np.random.default_rng(0)
    sub = man.sample(min(max_files, len(man)), random_state=0)
    stacked: dict[float, list[np.ndarray]] = {}
    for row in sub.itertuples(index=False):
        p = edir / f"{row.locus_id}_{row.architecture_id}_rep{int(row.replicate)}.npy"
        if not p.exists():
            continue
        beta = np.load(p)  # (ncsl, n_ancestries)
        if beta.ndim == 2 and beta.shape[1] >= 2:
            stacked.setdefault(float(row.rg), []).append(beta)
    if not stacked:
        return None
    pooled = {rg: np.vstack(v) for rg, v in stacked.items()}

    fig = plt.figure(figsize=(12.4, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.34)
    ax1 = fig.add_subplot(gs[0, 0])
    rgs = sorted(pooled)
    colors = fs.ordinal_colors(len(rgs))
    for i, (rg, c) in enumerate(zip(rgs, colors, strict=False), start=1):
        b = pooled[rg]
        corr = np.corrcoef(b.T)
        off = corr[np.triu_indices_from(corr, k=1)]
        off = off[np.isfinite(off)]
        # Every ancestry pair drawn as its own point: 15 pairs is few enough to show
        # individually, and the spread across pairs is the quantity of interest.
        ax1.scatter(
            np.full(len(off), i) + rng.uniform(-0.13, 0.13, len(off)),
            off,
            s=34,
            color=c,
            alpha=0.85,
            lw=0.6,
            edgecolor=fs.SURFACE,
            zorder=4,
        )
        ax1.plot([i - 0.34, i + 0.34], [rg, rg], color=fs.INK, lw=1.5, ls="--", zorder=6)
        # Below the tick labels, not on top of them: get_xaxis_transform puts y=0 at
        # the axes floor, which is where the tick text already is.
        ax1.text(
            i,
            -0.105,
            f"{len(b):,} $\\beta$\n{len(off)} pairs",
            transform=ax1.transAxes if False else ax1.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=8,
            color=fs.MUTED,
        )
        fs.direct_label(ax1, i + 0.34, float(off.mean()), f"{off.mean():.3f}", c, dx=6)
    ax1.set_xticks(range(1, len(rgs) + 1))
    ax1.set_xticklabels([f"{r:g}" for r in rgs])
    ax1.set_xlim(0.4, len(rgs) + 0.75)
    ax1.set_ylim(0.0, 1.06)
    ax1.set_xlabel("target cross-ancestry $r_g$", labelpad=26)
    ax1.set_ylabel("realised $\\beta$ correlation, per ancestry pair")
    ax1.text(
        0.02,
        0.02,
        "dashed = target",
        transform=ax1.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=fs.INK2,
    )
    fs.panel_letter(ax1, "a", "Recovered genetic correlation")

    pairs = [("EUR", "AFR"), ("EUR", "MID")]
    for k, (pair, ax) in enumerate(
        zip(pairs, [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])], strict=False)
    ):
        try:
            ia, ib = ANCESTRY_ORDER.index(pair[0]), ANCESTRY_ORDER.index(pair[1])
        except ValueError:
            continue
        allb = np.vstack([pooled[r] for r in rgs])
        if allb.shape[1] <= max(ia, ib):
            fs.no_data(ax, "ancestry column absent")
            continue
        pts = allb[:, [ia, ib]]
        show = pts if len(pts) <= 6000 else pts[rng.choice(len(pts), 6000, replace=False)]
        ax.scatter(show[:, 0], show[:, 1], s=9, alpha=0.25, lw=0, color=fs.BLUE_RAMP[2])
        lim = float(np.abs(show).max()) * 1.05
        ax.plot([-lim, lim], [-lim, lim], ls="--", lw=1.3, color=fs.INK2, zorder=5)
        r = float(np.corrcoef(pts[:, 0], pts[:, 1])[0, 1])
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"$\\beta$   {pair[0]}")
        ax.set_ylabel(f"$\\beta$   {pair[1]}")
        ax.text(
            0.035,
            0.965,
            f"r = {r:.2f}\nn = {len(pts):,}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color=fs.INK,
            bbox=dict(boxstyle="round,pad=0.3", facecolor=fs.SURFACE, edgecolor=fs.GRID),
        )
        fs.panel_letter(ax, "bc"[k], f"{pair[0]} vs {pair[1]}")

    fs.footnote(
        fig,
        f"Sampled {len(sub):,} of {len(man):,} instances. Panels b-c pool "
        "every target $r_g$, so their scatter is the design's spread across "
        "$r_g$ and not an estimation error.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# eda8 / eda9 -- the design itself
# --------------------------------------------------------------------------- #
def eda8_design_grid(pkg: dict, out: Path) -> Path | None:
    """The star design as a map. Worth one figure because 'extended, 15 points' is not
    a description anyone can hold in their head, and the two OFAT arms -- which every
    results figure slices along -- are only obvious when they are drawn."""
    man = pkg.get("manifest")
    if man is None:
        return None
    cells = man.groupby(["ncsl", "h2_target", "rg"]).size().rename("n").reset_index()
    h2s = sorted(cells["h2_target"].unique())
    rgs = sorted(cells["rg"].unique())
    ncsls = sorted(cells["ncsl"].unique())
    colors = dict(zip(ncsls, fs.ordinal_colors(len(ncsls)), strict=False))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.3))
    # Log scale FIRST: set_xscale rebuilds the locator and discards ticks set before it.
    ax1.set_xscale("log")
    ax1.minorticks_off()
    centre_rg, centre_h2 = max(rgs), (0.001 if 0.001 in h2s else h2s[len(h2s) // 2])

    for i, ncsl in enumerate(ncsls):
        arm = cells[(cells["ncsl"] == ncsl) & (cells["rg"] == centre_rg)]
        ax1.scatter(
            arm["h2_target"],
            np.full(len(arm), i),
            s=460,
            color=colors[ncsl],
            lw=1.0,
            edgecolor=fs.SURFACE,
            zorder=4,
        )
        for _, r in arm.iterrows():
            ax1.text(
                r["h2_target"],
                i,
                f"{int(r['n'])}",
                ha="center",
                va="center",
                fontsize=7.5,
                color=fs.SURFACE,
                fontweight="bold",
                zorder=6,
            )
        arm2 = cells[(cells["ncsl"] == ncsl) & (cells["h2_target"] == centre_h2)]
        ax2.scatter(
            arm2["rg"],
            np.full(len(arm2), i),
            s=460,
            color=colors[ncsl],
            lw=1.0,
            edgecolor=fs.SURFACE,
            zorder=4,
        )
        for _, r in arm2.iterrows():
            ax2.text(
                r["rg"],
                i,
                f"{int(r['n'])}",
                ha="center",
                va="center",
                fontsize=7.5,
                color=fs.SURFACE,
                fontweight="bold",
                zorder=6,
            )

    for ax, xs, xlab, title, letter in (
        (ax1, h2s, f"target per-locus $h^2$    ($r_g$ = {centre_rg:g})", "Heritability arm", "a"),
        (
            ax2,
            rgs,
            f"cross-ancestry $r_g$    ($h^2$ = {centre_h2:g})",
            "Genetic-correlation arm",
            "b",
        ),
    ):
        ax.set_yticks(range(len(ncsls)))
        ax.set_yticklabels([f"$n_{{csl}}$ = {int(n)}" for n in ncsls])
        ax.set_xticks(xs)
        ax.set_xticklabels([f"{v:g}" for v in xs])
        ax.set_xlabel(xlab)
        ax.set_ylim(-0.65, len(ncsls) - 0.35)
        ax.grid(axis="y", visible=False)
        fs.panel_letter(ax, letter, title)

    per_stratum = ""
    loci = pkg.get("loci")
    if loci is not None and "stratum" in loci.columns:
        counts = loci["stratum"].value_counts()
        per_stratum = (
            " The design is balanced within a cell and across strata up to "
            "the locus counts themselves — "
            + ", ".join(f"{s} {int(counts[s])} loci" for s in STRATUM_ORDER if s in counts)
            + f" x {int(man['replicate'].nunique())} replicates — so every "
            "design cell holds the same instances in the same proportion, "
            "and no comparison below is weighted by an accident of "
            "sampling."
        )
    fs.footnote(
        fig,
        "Numbers in each point are the instances simulated in that design "
        "cell (loci x replicates). The two arms share their centre cell, "
        "which is why the 15-point grid is not 3 x 3 x 3." + per_stratum,
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
FIGURES = [
    ("eda1_site_composition.png", eda1_site_composition),
    ("eda2_locus_selection.png", eda2_locus_selection),
    ("eda3_window_size.png", eda3_window_size),
    ("eda4_ld_divergence_maps.png", eda4_ld_divergence_maps),
    ("eda5_causal_maf.png", eda5_causal_maf),
    ("eda6_h2_calibration.png", eda6_h2_calibration),
    ("eda7_rg_recovery.png", eda7_rg_recovery),
    ("eda8_design_grid.png", eda8_design_grid),
]


def write_figures(data_root, out_dir, logger=None) -> list[Path]:
    """Draw every EDA figure the package supports. A figure that cannot be drawn is
    skipped with a note rather than failing the batch."""
    active = logger or log
    pkg = load_package(Path(data_root))
    fs.apply_style()
    out_dir = Path(out_dir)
    written = []
    for name, draw in FIGURES:
        try:
            p = draw(pkg, out_dir / name)
        except Exception as exc:  # noqa: BLE001 - one bad panel must not lose the batch
            active.warning("could not draw %s: %s", name, exc)
            continue
        if p is None:
            active.info("skipped %s (inputs absent)", name)
        else:
            written.append(p)
    active.info("wrote %d EDA figure(s) -> %s", len(written), out_dir)
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Render the fine-mapping EDA figures.")
    ap.add_argument(
        "--data-root",
        default="local/data/fine-mapping",
        help="the simulated package (holds loci/, ground_truth/, processed/)",
    )
    ap.add_argument(
        "--out-dir", default=None, help="default: <data-root>/../../output/fine-mapping/figures/eda"
    )
    args = ap.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )
    root = Path(args.data_root)
    if not root.exists():
        print(f"ERROR: {root} not found", file=sys.stderr)
        return 2
    default_out = root.parent.parent / "output/fine-mapping/figures/eda"
    out = Path(args.out_dir) if args.out_dir else default_out
    write_figures(root, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
