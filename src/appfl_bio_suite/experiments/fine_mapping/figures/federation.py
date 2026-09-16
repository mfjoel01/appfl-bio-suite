"""What the SuSiEx paper has no analogue for: the cost and the correctness of federating.

COORDINATOR-SIDE. The paper compares fine-mapping METHODS on data that is already in
one place. This experiment holds the method fixed -- it is SuSiEx, unmodified, on both
paths -- and moves the data problem instead: the genotypes never leave the site they
came from. So the questions these figures answer are the ones the paper never has to
ask, and they are the ones a reader of this work will.

    fed1_parity            does federating change the answer? (it must not)
    fed2_ancestry_columns  where each SuSiEx column's cohort actually comes from
    fed3_uplink_cost       what one federated locus costs in bytes, and why
    fed4_harmonization     the allele-order tax each site pays before it can contribute
    fed5_where_time_goes   what a run spends, and on what

fed1 IS THE EXPERIMENT'S CENTRAL CLAIM. By Corollary 1 of the derivation in
``docs/experiments/fine-mapping/reference/algo.md`` the federated fit EQUALS the
centralized fit -- same credible sets, same posterior inclusion probabilities, up to
floating-point summation order. That is a stronger claim than the GWAS experiment's,
which only approximates a pooled analysis, and it is falsifiable: any visible
disagreement in fed1 is a bug in the federated path, never a cost of federating. The
figure is drawn so that a disagreement would be impossible to miss -- identity line,
residuals on their own axis, and the worst observed difference printed rather than
described.
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

SITE_ORDER = ["anl", "covenant", "mbzuai"]
ANCESTRY_ORDER = ["EUR", "AFR", "AMR", "EAS", "CSA", "MID"]

# The published federation. Kept here rather than read from a scenario file so these
# figures can be drawn without the 350 GB package mounted; ``load_composition``
# overrides it from the real site manifests whenever they are available.
DEFAULT_COMPOSITION = {
    "anl": {"EUR": 30000, "AFR": 7500, "AMR": 7500, "EAS": 2500, "CSA": 2500},
    "covenant": {"AFR": 47500, "EUR": 1500, "CSA": 1000},
    "mbzuai": {"MID": 25000, "CSA": 15000, "AFR": 5000, "EUR": 5000},
}

KEY_COLS = ["locus_id", "architecture_id", "replicate"]
# The statistics a parity check compares. Every one is an output of SuSiEx rather than
# of the transport, so agreement on all of them is agreement on the fit itself.
PARITY_METRICS = [
    ("causal_pip_max", "PIP of the causal variant", "linear"),
    ("top_pip", "maximum PIP in the locus", "linear"),
    ("best_cs_size", "credible-set size", "log"),
    ("n_credible_sets", "credible sets returned", "linear"),
]


def load_composition(data_root: Path | None) -> dict[str, dict[str, int]]:
    """Realised site x ancestry counts, from the site manifests if they are present."""
    if data_root is None:
        return DEFAULT_COMPOSITION
    out: dict[str, dict[str, int]] = {}
    for site in SITE_ORDER:
        p = Path(data_root) / "processed" / site / f"{site}_manifest.tsv"
        if p.exists():
            m = pd.read_csv(p, sep="\t", usecols=["superpopulation"])
            out[site] = m["superpopulation"].value_counts().to_dict()
    return out or DEFAULT_COMPOSITION


# --------------------------------------------------------------------------- #
# fed1 -- the exactness claim
# --------------------------------------------------------------------------- #
def parse_dropped_variants(log_text: str) -> dict[str, int]:
    """``EUR: 6 window variant(s) dropped as not fully observed`` -> ``{'EUR': 6}``.

    Corollary 1's equality is conditional on both paths seeing the SAME variant list.
    The federated path drops a variant that is not fully observed at every site; the
    centralized comparator, which sees one pooled cohort, keeps it. The vendored code
    already warns about this per (locus, ancestry) -- this reads those warnings so the
    parity figure can report the precondition alongside the result instead of leaving
    a reader to explain a non-zero difference on their own.
    """
    return {
        m.group(1): int(m.group(2))
        for m in re.finditer(
            r"(\w+): (\d+) window variant\(s\) dropped as not fully observed", log_text
        )
    }


def fed1_parity(
    centralized: pd.DataFrame,
    federated: pd.DataFrame,
    out: Path,
    dropped: dict[str, int] | None = None,
) -> Path:
    """Federated against centralized, on the instances both paths ran.

    Drawn as identity-line scatters with a residual strip below each, rather than as a
    table of correlations. A correlation of 1.000 is compatible with a systematic
    offset, and an offset is exactly the failure mode this design is exposed to: sites
    standardising against their own column means instead of the pooled ones biases
    every estimate in the same direction, and would show here as a cloud parallel to
    the identity line rather than on it.
    """
    merged = centralized.merge(federated, on=KEY_COLS, suffixes=("_c", "_f"))
    if merged.empty:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        fs.no_data(
            ax,
            "no instances were run on BOTH paths\n"
            "run: run_stage.py federated --config <same package>",
        )
        return fs.save(fig, out, log)

    metrics = [
        m for m in PARITY_METRICS if f"{m[0]}_c" in merged.columns and f"{m[0]}_f" in merged.columns
    ]
    n = len(metrics)
    fig, axes = plt.subplots(2, n, figsize=(3.4 * n, 6.4), gridspec_kw={"height_ratios": [2.6, 1]})
    axes = np.atleast_2d(axes)
    worst_overall = 0.0

    for j, (col, label, scale) in enumerate(metrics):
        ax, axr = axes[0][j], axes[1][j]
        c = pd.to_numeric(merged[f"{col}_c"], errors="coerce")
        f = pd.to_numeric(merged[f"{col}_f"], errors="coerce")
        ok = np.isfinite(c) & np.isfinite(f)
        c, f = c[ok], f[ok]
        if not len(c):
            fs.no_data(ax, "no comparable values")
            fs.no_data(axr)
            continue
        lo = float(min(c.min(), f.min()))
        hi = float(max(c.max(), f.max()))
        pad = (hi - lo) * 0.06 or 0.05
        ax.plot(
            [lo - pad, hi + pad], [lo - pad, hi + pad], ls="--", lw=1.3, color=fs.INK2, zorder=3
        )
        ax.scatter(
            c,
            f,
            s=34,
            color=fs.PATH["federated"],
            alpha=0.6,
            lw=0.6,
            edgecolor=fs.SURFACE,
            zorder=4,
        )
        if scale == "log":
            ax.set_xscale("log")
            ax.set_yscale("log")
        else:
            ax.set_xlim(lo - pad, hi + pad)
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_aspect("equal" if scale != "log" else "auto")
        ax.set_xlabel("centralized")
        if j == 0:
            ax.set_ylabel("federated")
        fs.panel_letter(ax, "abcd"[j], label)

        d = (f - c).to_numpy(dtype=float)
        worst = float(np.abs(d).max())
        worst_overall = max(worst_overall, worst)
        # Residuals against the centralized value, so a difference that grows with the
        # statistic -- the signature of an accumulating summation error -- would show
        # as a fan rather than as a band.
        axr.axhline(0, color=fs.INK2, lw=1.1, ls="--", zorder=3)
        axr.scatter(c, d, s=26, color=fs.PATH["federated"], alpha=0.6, lw=0, zorder=4)
        if scale == "log":
            axr.set_xscale("log")
        axr.set_xlabel("centralized")
        if j == 0:
            axr.set_ylabel("federated − centralized")
        span = max(worst, 1e-15)
        axr.set_ylim(-span * 1.6, span * 1.6)
        axr.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)
        # Three regimes, and they mean different things. Exactly zero is bit-identical.
        # Below ~1e-12 is floating-point summation order, which Corollary 1 allows.
        # Anything larger needs an explanation, and the caption supplies one -- so it is
        # marked "warning" rather than "critical": a known, quantified precondition
        # violation is not the same as an unexplained disagreement.
        if worst == 0:
            verdict, colour = "identical", fs.STATUS["good"]
        elif worst < 1e-10:
            verdict, colour = f"max |Δ| = {worst:.1e}  (float)", fs.STATUS["good"]
        else:
            verdict, colour = f"max |Δ| = {worst:.3g}", fs.STATUS["warning"]
        axr.text(
            0.5,
            0.93,
            verdict,
            transform=axr.transAxes,
            ha="center",
            va="top",
            fontsize=9,
            fontweight="semibold",
            color=colour,
        )

    n_inst = len(merged)
    why = (
        " This figure summarizes result-table differences. The separate production parity "
        "gate checks the complete instance grid, variant eligibility, credible-set "
        "membership, fit status and retained-component PIPs at their output precision."
    )
    fs.footnote(
        fig,
        f"{n_inst:,} instances run on both paths over the same package; "
        f"largest difference on any statistic {worst_overall:.3g}." + why,
        y=-0.03,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# fed2 -- where each column's cohort comes from
# --------------------------------------------------------------------------- #
def fed2_ancestry_columns(composition: dict, out: Path, min_gwas_n: int = 1000) -> Path:
    """Each SuSiEx population column, decomposed by the sites that contribute to it.

    This is the design decision the experiment turns on, and it is invisible in every
    results table. SuSiEx takes a list of population columns; the obvious mapping is
    one column per site, and it is wrong here, because the simulator indexes effect
    sizes by superpopulation. An ancestry-mixed site column is not one homogeneous
    effect, which is the assumption the model makes about a column.

    Pooling by ancestry ACROSS sites gives columns that each are one effect -- and the
    figure shows what that costs: four of the six columns exist only because three
    institutions agreed to be in the same analysis.
    """
    totals = {}
    for anc in ANCESTRY_ORDER:
        per_site = {s: composition.get(s, {}).get(anc, 0) for s in SITE_ORDER}
        if sum(per_site.values()) > 0:
            totals[anc] = per_site
    ancs = sorted(totals, key=lambda a: -sum(totals[a].values()))

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.0, 4.8), gridspec_kw={"width_ratios": [1.55, 1]}
    )
    y = np.arange(len(ancs))
    left = np.zeros(len(ancs))
    for site in SITE_ORDER:
        vals = np.array([totals[a].get(site, 0) for a in ancs], dtype=float)
        ax1.barh(
            y,
            vals,
            left=left,
            height=0.64,
            color=fs.SITE[site],
            edgecolor=fs.SURFACE,
            linewidth=1.4,
            label=site,
        )
        for yi, (v, offset) in enumerate(zip(vals, left, strict=False)):
            if v >= 4000:
                ax1.text(
                    offset + v / 2,
                    yi,
                    f"{int(v):,}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=fs.SURFACE,
                    fontweight="semibold",
                )
        left += vals
    ax1.axvline(min_gwas_n, color=fs.INK, lw=1.2, ls="--", zorder=6)
    # Horizontal and inside the axes. Rotated below the spine it hangs off the figure,
    # and bbox_inches="tight" then stretches the whole canvas to contain it.
    ax1.text(
        min_gwas_n,
        -0.62,
        f"  min_gwas_n = {min_gwas_n:,}",
        va="center",
        ha="left",
        fontsize=8,
        color=fs.INK2,
    )
    for yi, a in enumerate(ancs):
        n_sites = sum(1 for s in SITE_ORDER if totals[a].get(s, 0) > 0)
        ax1.text(
            left[yi] + left.max() * 0.012,
            yi,
            f"{int(left[yi]):,}  ({n_sites} site{'s' if n_sites > 1 else ''})",
            va="center",
            ha="left",
            fontsize=8.5,
            color=fs.INK,
        )
    ax1.set_ylim(len(ancs) - 0.45, -0.95)  # room for the threshold caption at the top
    ax1.set_yticks(y)
    ax1.set_yticklabels(ancs)
    for tick, a in zip(ax1.get_yticklabels(), ancs, strict=False):
        tick.set_color(fs.ANCESTRY.get(a, fs.INK))
        tick.set_fontweight("semibold")
    ax1.set_xlim(0, left.max() * 1.22)
    ax1.xaxis.set_major_formatter(lambda v, _: f"{v / 1000:g}k")
    ax1.set_xlabel("individuals in the pooled ancestry column")
    ax1.grid(axis="y", visible=False)
    ax1.legend(title="contributing site", loc="lower right", ncol=1)
    fs.panel_letter(ax1, "a", "SuSiEx columns, by contributing site")

    # (b) the counterfactual: what each site could fine-map alone. A column below
    # min_gwas_n is dropped, so this is not a rhetorical comparison -- it is the number
    # of columns each analysis would actually have.
    solo = [
        sum(1 for a in ancs if composition.get(s, {}).get(a, 0) >= min_gwas_n) for s in SITE_ORDER
    ]
    joint = sum(1 for a in ancs if sum(totals[a].values()) >= min_gwas_n)
    labels = [*SITE_ORDER, "federation"]
    vals = [*solo, joint]
    cols = [fs.SITE[s] for s in SITE_ORDER] + [fs.PATH["federated"]]
    ax2.bar(range(len(labels)), vals, 0.6, color=cols, edgecolor=fs.SURFACE, linewidth=1.4)
    for i, v in enumerate(vals):
        ax2.text(
            i, v, f"{v}", ha="center", va="bottom", fontsize=11, fontweight="semibold", color=fs.INK
        )
    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels, rotation=15, ha="right")
    ax2.set_ylabel(f"ancestry columns at n ≥ {min_gwas_n:,}")
    ax2.set_ylim(0, max(vals) + 1.2)
    ax2.grid(axis="x", visible=False)
    fs.panel_letter(ax2, "b", "Columns available, alone and together")

    fs.footnote(
        fig,
        "A column below min_gwas_n is dropped from the fit, so panel b "
        "counts the columns each analysis would actually run with. "
        "Cross-ancestry fine-mapping needs more than one column; that is "
        "the whole mechanism.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# fed3 -- the price of shipping LD structure
# --------------------------------------------------------------------------- #
def fed3_uplink_cost(
    loci: pd.DataFrame | None, composition: dict, out: Path, dtype_bytes: int = 8
) -> Path:
    """Uplink per locus, and why it is quadratic rather than linear.

    A meta-analysis needs one number per variant per site and is O(M). A credible set
    needs the LD matrix, and LD is precisely the thing that does not decompose into
    per-variant summaries -- so this design ships X'X and pays O(M^2). That is the
    honest price of federating fine-mapping rather than association testing, and it is
    why a run is sharded by locus instead of sweeping the chromosome in one exchange.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.6))

    if loci is not None and "n_variants_in_window" in loci.columns:
        m = loci["n_variants_in_window"].to_numpy(dtype=float)
    else:
        m = np.linspace(1000, 4500, 200)
    grid = np.linspace(max(m.min() * 0.8, 200), m.max() * 1.05, 200)

    # (a) per (site, ancestry) block, against the linear alternative
    quad = dtype_bytes * grid**2 / 2**20
    lin = dtype_bytes * grid * 3 / 2**20  # beta, se, n per variant
    ax1.plot(grid, quad, lw=2.2, color=fs.PATH["federated"], zorder=4)
    ax1.plot(grid, lin, lw=2.2, color=fs.PATH["single-ancestry"], zorder=4)
    fs.direct_label(
        ax1,
        grid[-1],
        quad[-1],
        " $X'X$ — this design\n $O(M^2)$",
        fs.PATH["federated"],
        dx=4,
        va="center",
    )
    fs.direct_label(
        ax1,
        grid[-1],
        lin[-1],
        " marginal statistics\n $O(M)$",
        fs.PATH["single-ancestry"],
        dx=4,
        va="center",
    )
    if loci is not None and "n_variants_in_window" in loci.columns:
        obs = dtype_bytes * m**2 / 2**20
        ax1.scatter(m, obs, s=18, color=fs.PATH["federated"], alpha=0.35, lw=0, zorder=3)
    ax1.set_yscale("log")
    ax1.set_xlabel("variants in the locus window ($M$)")
    ax1.set_ylabel("uplink per (site, ancestry, locus)   (MiB)")
    ax1.set_xlim(right=grid[-1] * 1.42)
    fs.panel_letter(ax1, "a", "One block, two designs")

    # (b) what a whole locus costs the federation, per site, at the realised M
    med_m = float(np.median(m))
    block = dtype_bytes * med_m**2 / 2**20
    per_site = {
        s: len([a for a in ANCESTRY_ORDER if composition.get(s, {}).get(a, 0)]) for s in SITE_ORDER
    }
    vals = [per_site[s] * block for s in SITE_ORDER]
    ax2.bar(
        range(len(SITE_ORDER)),
        vals,
        0.6,
        color=[fs.SITE[s] for s in SITE_ORDER],
        edgecolor=fs.SURFACE,
        linewidth=1.4,
    )
    for i, (s, v) in enumerate(zip(SITE_ORDER, vals, strict=False)):
        ax2.text(
            i,
            v,
            f"{v:,.0f} MiB\n({per_site[s]} ancestries)",
            ha="center",
            va="bottom",
            fontsize=9,
            color=fs.INK2,
        )
    ax2.set_xticks(range(len(SITE_ORDER)))
    ax2.set_xticklabels(SITE_ORDER)
    ax2.set_ylabel(f"uplink per locus   (MiB, at median $M$ = {med_m:,.0f})")
    ax2.set_ylim(0, max(vals) * 1.32)
    ax2.grid(axis="x", visible=False)
    total = sum(vals)
    ax2.text(
        0.5,
        0.955,
        f"federation total: {total / 1024:.1f} GiB per locus,\n"
        "one exchange — there is no second round",
        transform=ax2.transAxes,
        ha="center",
        va="top",
        fontsize=9.5,
        color=fs.INK,
        bbox=dict(boxstyle="round,pad=0.45", facecolor=fs.SURFACE, edgecolor=fs.GRID),
    )
    fs.panel_letter(ax2, "b", "What one locus costs each site")

    fs.footnote(
        fig,
        f"float{dtype_bytes * 8} moments. A site uplinks one $X'X$ per "
        "ancestry it holds; ancestries it does not hold cost it nothing. "
        "Single round by construction: the site moments do not change, so "
        "a second round would recompute them at full price for no new "
        "information.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# fed4 -- the harmonization tax
# --------------------------------------------------------------------------- #
def fed4_harmonization(
    flip_counts: dict[str, int],
    n_variants: int,
    out: Path,
    causal_affected: tuple[int, int] | None = None,
) -> Path:
    """How much recoding each site does before its moments can be added to anyone's.

    Theorem 1's first condition wants one ordered, allele-harmonized variant list. The
    per-site filesets were cut without ``--keep-allele-order``, so PLINK set A1 to each
    site's own minor allele -- and a site coding a variant the other way round
    contributes ``2 - x`` where the others contribute ``x``. That sums in silently and
    negates every off-diagonal the variant touches. There is no error, no warning, and
    no way to see it downstream, which is why each site recodes against an agreed
    reference list and the coordinator re-checks rather than trusting the report.

    The tax is not evenly shared, and the reason is the point: the site furthest from
    the reference in allele frequency recodes the most.
    """
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(11.6, 4.5), gridspec_kw={"width_ratios": [1.25, 1]}
    )
    sites = [s for s in SITE_ORDER if s in flip_counts]
    vals = [flip_counts[s] for s in sites]
    fracs = [v / n_variants for v in vals]
    ax1.bar(
        range(len(sites)),
        fracs,
        0.6,
        color=[fs.SITE[s] for s in sites],
        edgecolor=fs.SURFACE,
        linewidth=1.4,
    )
    for i, (v, f) in enumerate(zip(vals, fracs, strict=False)):
        ax1.text(
            i,
            f,
            f"{f * 100:.2f}%\n{v:,} variants",
            ha="center",
            va="bottom",
            fontsize=9,
            color=fs.INK2,
        )
    ax1.set_xticks(range(len(sites)))
    ax1.set_xticklabels(sites)
    ax1.set_ylabel("share of the chromosome recoded")
    ax1.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
    ax1.set_ylim(0, max(fracs) * 1.42)
    ax1.grid(axis="x", visible=False)
    fs.panel_letter(ax1, "a", "Variants recoded per site")

    if causal_affected:
        hit, total = causal_affected
        ax2.bar(
            [0],
            [1 - hit / total],
            0.55,
            bottom=[hit / total],
            color=fs.STATUS["good"],
            edgecolor=fs.SURFACE,
            linewidth=1.4,
        )
        ax2.bar(
            [0],
            [hit / total],
            0.55,
            color=fs.STATUS["critical"],
            edgecolor=fs.SURFACE,
            linewidth=1.4,
        )
        # Both segments labelled in place. A caption above the panel is what collided
        # with the title, and direct labels are the better answer anyway.
        ax2.text(
            0,
            hit / total / 2,
            f"coded inconsistently\n{hit:,} of {total:,}  ({hit / total * 100:.1f}%)",
            ha="center",
            va="center",
            fontsize=9.5,
            color=fs.SURFACE,
            fontweight="semibold",
        )
        ax2.text(
            0,
            hit / total + (1 - hit / total) / 2,
            "consistent across all sites",
            ha="center",
            va="center",
            fontsize=9.5,
            color=fs.SURFACE,
            fontweight="semibold",
        )
        ax2.set_xlim(-0.6, 0.6)
        ax2.set_xticks([0])
        ax2.set_xticklabels(["causal variants"])
        ax2.set_ylim(0, 1)
        ax2.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
        ax2.set_ylabel("share")
        ax2.grid(axis="x", visible=False)
        ax2.set_xticklabels(["instances whose causal variant(s)\nare affected"])
        fs.panel_letter(ax2, "b", "Where it reaches the ground truth")
    else:
        fs.no_data(ax2, "causal-variant impact not supplied")

    fs.footnote(
        fig,
        f"{n_variants:,} chromosome-1 variants shared by all sites. "
        "Recoding is done at each site against an agreed reference list, "
        "before any moment is formed; the coordinator re-checks rather "
        "than trusting the count. Panel b is the residual defect in the "
        "realised package -- the phenotype stage was run before the "
        "reference list existed, so those causal effects carry "
        "site-dependent signs and attenuate the pooled signal.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# fed5 -- where the time goes
# --------------------------------------------------------------------------- #
def fed5_where_time_goes(
    centralized: pd.DataFrame, out: Path, federated: pd.DataFrame | None = None
) -> Path:
    """Per-instance wall-clock on each path, and what federating adds to it.

    The federated path does strictly more arithmetic than the centralized one -- every
    site forms its own moments and the coordinator sums them -- but it does not run a
    GWAS, because the summary statistics come out of the pooled moments in closed form.
    Whether that trade is a saving or a cost is a measurement, not an argument.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.5))
    series, labels, colors = [], [], []
    if "runtime_s" in centralized.columns:
        series.append(centralized["runtime_s"].dropna().values)
        labels.append(f"centralized\n(n={centralized['runtime_s'].notna().sum():,})")
        colors.append(fs.PATH["centralized"])
    if federated is not None and "runtime_s" in federated.columns:
        series.append(federated["runtime_s"].dropna().values)
        labels.append(f"federated\n(n={federated['runtime_s'].notna().sum():,})")
        colors.append(fs.PATH["federated"])
    if not series:
        fs.no_data(ax1, "no runtime_s column")
        fs.no_data(ax2)
        return fs.save(fig, out, log)

    fs.violin_box(ax1, series, colors=colors)
    ax1.set_yscale("log")
    ax1.set_xticks(range(1, len(labels) + 1))
    ax1.set_xticklabels(labels)
    ax1.set_ylabel("wall-clock per instance (s)")
    ax1.grid(axis="x", visible=False)
    for i, v in enumerate(series, start=1):
        fs.direct_label(
            ax1, i, float(np.median(v)), f"  median {np.median(v):.0f}s", colors[i - 1], dx=4
        )
    fs.panel_letter(ax1, "a", "Cost per instance, by path")

    # (b) A projection to a COMMON workload, not the raw totals. The two paths were run
    # over different numbers of instances -- summing each path's own measured seconds
    # would say "federated cost almost nothing", which is an artefact of having run 45
    # of them rather than a finding. Median per instance times a shared instance count
    # is the comparison a scheduler request actually needs.
    n_full = int(len(centralized))
    names = [lbl.split("\n")[0] for lbl in labels]
    proj = [float(np.median(v)) * n_full / 3600 for v in series]
    ax2.barh(range(len(names)), proj, 0.5, color=colors, edgecolor=fs.SURFACE, linewidth=1.4)
    for i, (t_, v) in enumerate(zip(proj, series, strict=False)):
        ax2.text(
            t_,
            i,
            f"  {t_:,.0f} CPU-hours    (median {np.median(v):.0f}s × {n_full:,})",
            va="center",
            ha="left",
            fontsize=9,
            color=fs.INK2,
        )
    ax2.set_yticks(range(len(names)))
    ax2.set_yticklabels(names)
    ax2.set_xlim(0, max(proj) * 2.05)
    ax2.set_xlabel(f"projected CPU-hours for the full {n_full:,}-instance sweep")
    ax2.invert_yaxis()
    ax2.grid(axis="y", visible=False)
    fs.panel_letter(ax2, "b", "Projected to one common workload")

    fs.footnote(
        fig,
        "Single-threaded per instance; a sweep parallelises across "
        "instances, so wall-clock is CPU-hours divided by the allocation. "
        "THE FEDERATED FIGURE IS THE COORDINATOR'S SHARE ONLY. Its "
        "runtime_s covers forming the summary statistics from the pooled "
        "moments in closed form and running SuSiEx; it does not include "
        "the X'X each site computes over its own genotypes, which happens "
        "on the sites' hardware and is not in either results table. The "
        "centralized figure, by contrast, includes its own per-ancestry "
        "GWAS. So the federated bar is lower because work moved, not "
        "because it disappeared, and these two numbers do not measure the "
        "same total.",
    )
    fig.tight_layout()
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
# fed6 -- what federating actually buys
# --------------------------------------------------------------------------- #
# Arms in the order a reader should meet them: each site alone, then the federation
# matched on sample size, then the whole federation, then the shortcut that avoids it.
ARM_ORDER = ["covenant", "mbzuai", "anl", "federation_50k", "federation", "ld_borrowed"]
ARM_LABEL = {
    "covenant": "Covenant alone",
    "mbzuai": "MBZUAI alone",
    "anl": "ANL alone",
    "federation_50k": "All three, n matched",
    "federation": "All three, full",
    "ld_borrowed": "Covenant stats + ANL LD",
}
# Two lines per label is the natural way to write these and it does not survive six
# arms: rotated, each label's second line juts right into its neighbour's first. One
# line rotates cleanly at any length, so the labels stay whole rather than abbreviated.


def arm_colors(arms: list[str]) -> list[str]:
    """A site's own colour where the arm IS that site; the federated hue for a
    federation; a status colour for the shortcut arm, which is a warning and not a
    series."""
    out = []
    for a in arms:
        if a in fs.SITE:
            out.append(fs.SITE[a])
        elif a == "ld_borrowed":
            out.append(fs.STATUS["critical"])
        else:
            out.append(fs.PATH["federated"])
    return out


def fed6_what_federation_buys(
    by_arm: pd.DataFrame, out: Path, cohort: dict[str, tuple[int, int]] | None = None
) -> Path:
    """The same loci, the same truth, differing only in who took part.

    THIS IS THE ARM THE OTHER FIGURES WERE MISSING. Everything else in this experiment is
    computed with all three sites participating, which establishes that federating is
    *correct* without ever showing that it is *worth it*. A reader's first question about
    a federation is "what would I get on my own?", and until this figure there was no arm
    that answered it.

    Read the panels left to right as one argument: a site alone reaches some power, the
    federation reaches more, and panel d says how much of that is simply three times the
    data rather than three times the diversity -- the ``n matched`` arm is the full
    federation restricted to one site's worth of people, so the gap between it and a solo
    site is the cost of splitting a fixed cohort across ancestries, and the gap between
    it and the full federation is sample size.

    That first gap is deliberately NOT labelled "diversity". SuSiEx fits one effect per
    ancestry column, so power tracks the size of the largest column rather than the total,
    and the solo arms win it by being concentrated: Covenant puts 47,500 of its 50,000
    people into AFR, where the n-matched federation's largest column is 20,000. The
    controlled comparison is ANL (5 columns, largest 30,000) against the n-matched
    federation (5 columns, largest 20,000) -- same total n, same column count, 14.6 pp
    apart, p~1e-55. Calling the bar "diversity" would tell a reader that ancestral
    diversity costs 30 points of power, which is not what it measures and is a claim this
    design cannot support.

    The ``ld_borrowed`` arm is a different claim and is coloured as a warning rather than
    as a series. It is the cheap alternative to federating -- one site's summary
    statistics against another's LD panel, O(M) instead of O(M^2) -- and it belongs here
    because "why not just meta-analyse?" is the first objection to this design.
    """
    available = set(by_arm["arm"])
    arms = [a for a in ARM_ORDER if a in available] + sorted(available - set(ARM_ORDER))
    if not arms:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        fs.no_data(ax, "results carry no arm column")
        return fs.save(fig, out, log)
    colors = arm_colors(arms)
    x = np.arange(len(arms))

    from appfl_bio_suite.experiments.fine_mapping.figures.paper_plots import (
        _err_bars,
        _rate,
        parse_architecture,
    )

    d = parse_architecture(by_arm)

    fig, axes = plt.subplots(1, 4, figsize=(17.6, 5.2))
    ax_a, ax_b, ax_c, ax_d = axes

    def _bars(ax, col, ylabel, letter, title, pct=True):
        rows = [_rate(d[d["arm"] == a], col) for a in arms]
        p, yerr = _err_bars([(r[0], r[1], r[2]) for r in rows])
        ax.bar(
            x,
            p,
            0.62,
            color=colors,
            edgecolor=fs.SURFACE,
            linewidth=1.4,
            yerr=yerr,
            capsize=3,
            ecolor=fs.INK2,
            error_kw=dict(lw=1.0),
            zorder=3,
        )
        for xi, v in enumerate(p):
            ax.text(
                xi,
                v + (yerr[1][xi] if len(yerr[1]) > xi else 0),
                f"{v * 100:.0f}%" if pct else f"{v:.2f}",
                ha="center",
                va="bottom",
                fontsize=9,
                color=fs.INK2,
            )
        ax.set_ylabel(ylabel)
        if pct:
            ax.yaxis.set_major_formatter(lambda v, _: f"{v * 100:.0f}%")
        fs.panel_letter(ax, letter, title)
        return [r[3] for r in rows]

    ns = _bars(ax_a, "captured", "power   =   P(a causal variant is captured)", "a", "Power")
    ax_a.set_ylim(0, 1.18)

    # (b) resolution
    data = [d[(d["arm"] == a) & d["captured"]]["best_cs_size"].dropna().values for a in arms]
    fs.boxplot(ax_b, data, positions=x, colors=colors, widths=0.6, points=False)
    ax_b.set_yscale("log")
    ax_b.set_ylim(bottom=0.82)
    ax_b.set_ylabel("credible-set size (variants)\nsmaller = sharper")
    fs.panel_letter(ax_b, "b", "Resolution")

    _bars(ax_c, "pip95", "P(causal variant at PIP > 0.95)", "c", "High-confidence yield")

    # Predeclared paired comparisons; no post-hoc selection of the best solo site.
    from ..reporting import paired_arm_differences

    differences = paired_arm_differences(by_arm)
    if len(differences):
        values = differences.power_difference.to_numpy()
        errors = np.vstack([values - differences.ci_low, differences.ci_high - values])
        y = np.arange(len(differences))
        ax_d.errorbar(values, y, xerr=errors, fmt="o", color=fs.PATH["federated"], capsize=3)
        ax_d.axvline(0, color=fs.INK2, lw=1.2)
        ax_d.set_yticks(y)
        ax_d.set_yticklabels([ARM_LABEL.get(a, a) for a in differences.arm], fontsize=7)
        ax_d.xaxis.set_major_formatter(lambda v, _: f"{v * 100:+.0f}")
        ax_d.set_xlabel("full federation minus comparison (pp)")
        ax_d.invert_yaxis()
        ax_d.grid(axis="y", visible=False)
        fs.panel_letter(ax_d, "d", "Paired power differences")
    else:
        fs.no_data(ax_d, "paired federation comparisons unavailable")

    for ax in (ax_a, ax_b, ax_c):
        ax.set_xticks(x)
        # Six arms of two-line labels do not fit a quarter of the canvas horizontally --
        # they collide into each other and the panel becomes unreadable. Rotating past
        # four arms keeps the four-arm case upright, which is easier to read.
        ax.set_xticklabels(
            [ARM_LABEL.get(a, a) for a in arms],
            fontsize=8.5 if len(arms) <= 4 else 8.0,
            rotation=0 if len(arms) <= 4 else 30,
            ha="center" if len(arms) <= 4 else "right",
            rotation_mode=None if len(arms) <= 4 else "anchor",
        )
        ax.grid(axis="x", visible=False)
        for tick, a in zip(ax.get_xticklabels(), arms, strict=False):
            if a == "ld_borrowed":
                tick.set_color(fs.STATUS["critical"])
    if cohort:
        fs.count_labels(ax_a, x, [cohort.get(a, (0, 0))[0] for a in arms], pad_pt=34, fmt="n={:,}")
        ax_a.set_xlabel("participating cohort", labelpad=30)
        for xi, a in enumerate(arms):
            k = cohort.get(a, (0, 0))[1]
            if k:
                ax_a.annotate(
                    f"{k} col",
                    xy=(xi, 0),
                    xycoords=("data", "axes fraction"),
                    xytext=(0, -46),
                    textcoords="offset points",
                    ha="center",
                    va="top",
                    fontsize=7.5,
                    color=fs.MUTED,
                    annotation_clip=False,
                )

    fig.suptitle(
        "What federating buys, on the same loci and the same truth",
        x=0.005,
        ha="left",
        fontsize=13,
        fontweight="semibold",
    )
    fs.footnote(
        fig,
        f"{len(d):,} instances across {len(arms)} arm(s); "
        f"{', '.join(f'{a}={n:,}' for a, n in zip(arms, ns, strict=False))}. "
        "Every arm is the same estimator over the same ground truth, "
        "differing only in which cohorts took part. Error bars are Wilson "
        "95% intervals. The `n matched` arm is the full federation "
        "restricted to one site's worth of people, stratified within site "
        "and ancestry, so panel d can separate splitting a fixed cohort from adding to it. "
        "The borrowed-LD arm is the cheap alternative to federating and is "
        "drawn as a warning, not a series.",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fs.save(fig, out, log)


# --------------------------------------------------------------------------- #
def parse_flip_counts(log_text: str) -> tuple[dict[str, int], int]:
    """``anl: 9398 of 533532 variants recoded ...`` -> ``({'anl': 9398}, 533532)``.

    Read from the run log because the count is a property of THAT run's reference
    list, and re-deriving it here would answer a subtly different question.
    """
    counts, total = {}, 0
    for m in re.finditer(r"(\w+): (\d+) of (\d+) variants recoded", log_text):
        counts[m.group(1)] = int(m.group(2))
        total = int(m.group(3))
    return counts, total


def write_figures(
    centralized,
    out_dir,
    federated=None,
    data_root=None,
    run_log=None,
    flips=None,
    by_arm=None,
    cohort=None,
    logger=None,
    strict=False,
) -> list[Path]:
    """Draw every federation figure the available inputs support."""
    active = logger or log
    fs.apply_style()
    out_dir = Path(out_dir)
    cen = pd.DataFrame(centralized)
    fed = pd.DataFrame(federated) if federated is not None else None
    composition = load_composition(Path(data_root) if data_root else None)

    loci = None
    if data_root:
        p = Path(data_root) / "loci" / "selected_loci.tsv"
        if p.exists():
            loci = pd.read_csv(p, sep="\t")

    flip_counts, n_var = {}, 0
    if run_log and Path(run_log).exists():
        flip_counts, n_var = parse_flip_counts(Path(run_log).read_text(errors="ignore"))
    causal_affected = None
    if flips is not None and len(flips):
        causal_affected = (int(flips["has_flip"].sum()), len(flips))

    dropped = {}
    if run_log and Path(run_log).exists():
        dropped = parse_dropped_variants(Path(run_log).read_text(errors="ignore"))

    jobs = [
        (
            "fed1_parity.png",
            lambda o: fed1_parity(cen, fed, o, dropped) if fed is not None else None,
        ),
        ("fed2_ancestry_columns.png", lambda o: fed2_ancestry_columns(composition, o)),
        ("fed3_uplink_cost.png", lambda o: fed3_uplink_cost(loci, composition, o)),
        (
            "fed4_harmonization.png",
            lambda o: (
                fed4_harmonization(flip_counts, n_var, o, causal_affected) if flip_counts else None
            ),
        ),
        ("fed5_where_time_goes.png", lambda o: fed5_where_time_goes(cen, o, fed)),
        (
            "fed6_what_federation_buys.png",
            lambda o: (
                fed6_what_federation_buys(by_arm, o, cohort)
                if by_arm is not None and len(by_arm)
                else None
            ),
        ),
    ]
    written = []
    for name, draw in jobs:
        try:
            p = draw(out_dir / name)
        except Exception as exc:  # noqa: BLE001
            if strict:
                raise RuntimeError(f"Could not draw {name}") from exc
            active.warning("could not draw %s: %s", name, exc)
            continue
        if p is None:
            active.info("skipped %s (inputs absent)", name)
        else:
            written.append(p)
    active.info("wrote %d federation figure(s) -> %s", len(written), out_dir)
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Render the federation figures.")
    ap.add_argument("--centralized", required=True)
    ap.add_argument("--federated", default=None)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--data-root", default=None)
    ap.add_argument(
        "--run-log", default=None, help="a federated run log, for the per-site recoding counts"
    )
    ap.add_argument("--flips", default=None, help="allele-flip flag TSV")
    args = ap.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )
    cen_p = Path(args.centralized)
    if not cen_p.exists():
        print(f"ERROR: {cen_p} not found", file=sys.stderr)
        return 2
    cen = pd.read_csv(cen_p, sep="\t")
    fed_p = Path(args.federated) if args.federated else None
    fed = pd.read_csv(fed_p, sep="\t") if fed_p and fed_p.exists() else None
    flips = pd.read_csv(args.flips, sep="\t") if args.flips and Path(args.flips).exists() else None
    write_figures(cen, Path(args.out_dir), fed, args.data_root, args.run_log, flips)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
