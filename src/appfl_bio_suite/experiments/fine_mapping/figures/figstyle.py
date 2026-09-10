"""The one visual system every fine-mapping figure is drawn in. COORDINATOR-SIDE.

WHY A SHARED MODULE RATHER THAN RCPARAMS PER FILE
--------------------------------------------------
These figures are read side by side -- a paper-analogue panel next to a federation
cost curve next to an EDA histogram -- so a colour that means ``AFR`` in one must not
mean ``covenant`` in the next. Every palette below is assigned to *one* entity and
imported rather than retyped, which is the only mechanism that actually holds across
four modules and thirty figures.

THE PALETTES ARE VALIDATED, NOT CHOSEN BY EYE
----------------------------------------------
Each set below was run through the data-viz palette validator (lightness band, chroma
floor, colour-vision-deficiency separation, normal-vision separation, surface contrast).
The recorded results, on the light surface these figures render on:

    ANCESTRY (6)   adjacent-pair gates PASS  (worst CVD dE 9.1, normal-vision dE 19.6)
                   ALL-PAIR gates FAIL       (green<->orange dE 3.2 protan)
    SITE (3)       all-pair gates PASS       (worst CVD dE 15.3, normal-vision dE 20.8)
    PATH (2)       all-pair gates PASS       (worst CVD dE 24.7, normal-vision dE 33.6)
    ordinal ramps  monotone lightness, adjacent dL >= 0.06, light end clears 2:1

The ancestry all-pair failure is load-bearing and is why :func:`ancestry_colors` is
documented for bars and stacks only. Six ancestries in one scatter cannot be told apart
by colour; the SuSiEx paper hits the same wall and solves it the same way, by faceting
one ancestry pair per panel (its Figure 4f/4g). Do that, do not add a seventh hue.

Three slots -- aqua, yellow, magenta -- sit below 3:1 against the surface. The relief
rule applies: anything drawn in them carries a visible direct label or appears in a
table, never colour alone. :func:`direct_label` is the tool for that.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # before pyplot: a headless node has no display to reach

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

__all__ = [
    "ANCESTRY",
    "SITE",
    "PATH",
    "STATUS",
    "INK",
    "INK2",
    "MUTED",
    "GRID",
    "AXIS",
    "SURFACE",
    "BLUE_RAMP",
    "ORANGE_RAMP",
    "DIVERGING",
    "apply_style",
    "ancestry_colors",
    "ordinal_colors",
    "panel_letter",
    "direct_label",
    "wilson",
    "binom_summary",
    "save",
    "legend_swatches",
    "boxplot",
    "violin_box",
    "footnote",
    "no_data",
    "count_labels",
    "title_with_letter",
    "VIOLET_RAMP",
    "spread_labels",
]

# --------------------------------------------------------------------------- #
# Chrome
# --------------------------------------------------------------------------- #
SURFACE = "#fcfcfb"
INK = "#0b0b0b"  # primary ink: titles, medians, values
INK2 = "#52514e"  # secondary ink: annotations
MUTED = "#898781"  # axis labels, tick marks
GRID = "#e1e0d9"  # hairline gridlines
AXIS = "#c3c2b7"  # baselines and spines

# --------------------------------------------------------------------------- #
# Categorical: ancestry. FIXED ORDER, never cycled, never re-sorted by value.
#
# The order is the one the config declares and the results table's min_p_* columns
# follow, so a legend here and a column there name the same thing in the same place.
# BARS AND STACKS ONLY -- see the module docstring on the all-pairs failure.
# --------------------------------------------------------------------------- #
ANCESTRY = {
    "EUR": "#2a78d6",  # blue
    "AFR": "#eb6834",  # orange
    "AMR": "#1baf7a",  # aqua      (contrast 2.74 -- direct-label it)
    "EAS": "#eda100",  # yellow    (contrast 2.11 -- direct-label it)
    "CSA": "#e87ba4",  # magenta   (contrast 2.62 -- direct-label it)
    "MID": "#008300",  # green
}

# Categorical: site. Deliberately a different family from ancestry, because a stacked
# bar of sites within an ancestry column puts both dimensions on one chart.
SITE = {
    "anl": "#4a3aa7",  # violet
    "covenant": "#e34948",  # red
    "mbzuai": "#eda100",  # yellow (direct-label it)
}

# Categorical: which analysis produced the number. Two slots, maximum contrast.
PATH = {
    "centralized": "#2a78d6",
    "federated": "#eb6834",
    "single-ancestry": "#898781",  # the weaker comparator reads as recessive on purpose
}

# Reserved. Never reused as "series 4"; always shipped with a label, never colour alone.
STATUS = {
    "good": "#0ca30c",
    "warning": "#fab219",
    "serious": "#ec835a",
    "critical": "#d03b3b",
}

# --------------------------------------------------------------------------- #
# Sequential. One hue, light -> dark, for any ORDERED factor: h2, rg, ncsl, stratum.
#
# Blue is the default and carries whichever ordered factor is the colour in a given
# chart. Orange is the second, used only when two ordered factors are coloured at once.
# Both are ordinal-validated: the light end clears 2:1 against the surface, so a step
# never dissolves into the background the way a true sequential ramp's 100 step may.
# --------------------------------------------------------------------------- #
BLUE_RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
ORANGE_RAMP = ["#f4a07a", "#eb6834", "#a8410f"]
VIOLET_RAMP = ["#a79fe8", "#6a5bc9", "#413094"]

# Diverging, for signed quantities only (a PIP difference, a beta difference). Two
# hues with a NEUTRAL grey midpoint -- never a hue at zero, or zero reads as a value.
DIVERGING = ("#1c5cab", "#f0efec", "#c0392b")


def apply_style() -> None:
    """Install the shared rcParams. Idempotent; call once per process."""
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.facecolor": SURFACE,
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans"],
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "semibold",
            "axes.labelsize": 10,
            "axes.edgecolor": AXIS,
            "axes.linewidth": 0.8,
            "axes.labelcolor": INK2,
            "text.color": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.axisbelow": True,  # grid behind the data, always
            "grid.color": GRID,
            "grid.linewidth": 0.7,
            "legend.frameon": False,
            "legend.fontsize": 9,
            "lines.linewidth": 2.0,  # 2px lines
            "lines.markersize": 6,  # >= 8px on the page at 300dpi
            "patch.linewidth": 0,
        }
    )


def ancestry_colors(names) -> list[str]:
    """Colours for a list of ancestry names, in the caller's order.

    Unknown names fall back to muted grey rather than a generated hue: a new ancestry
    is a decision about the palette, not something a plotting call invents.
    """
    return [ANCESTRY.get(str(n), MUTED) for n in names]


def ordinal_colors(n: int, ramp: list[str] | None = None) -> list[str]:
    """``n`` steps of an ordinal ramp, light -> dark, evenly spread over the ramp.

    For ``n`` at or below the ramp length this returns the ramp's own validated steps
    rather than interpolating, so the published contrast numbers still hold.
    """
    ramp = ramp or BLUE_RAMP
    if n <= 0:
        return []
    if n == 1:
        return [ramp[len(ramp) // 2]]
    if n <= len(ramp):
        idx = np.linspace(0, len(ramp) - 1, n).round().astype(int)
        return [ramp[i] for i in idx]
    # More steps than the ramp has: interpolate in RGB between its ends. Rare, and the
    # ordinal contrast guarantee degrades to the ramp's endpoints, which still hold.
    import matplotlib.colors as mcolors

    cmap = mcolors.LinearSegmentedColormap.from_list("ramp", ramp)
    return [mcolors.to_hex(cmap(i / (n - 1))) for i in range(n)]


def panel_letter(ax, letter: str, title: str = "", dx: float = -0.085, dy: float = 1.06) -> None:
    """Bold panel letter above-left of the axes, in the Nature/SuSiEx convention.

    Placed in axes coordinates rather than figure coordinates so it tracks the panel
    through ``tight_layout``, which is what stops letters drifting when a figure is
    resized.
    """
    ax.text(
        dx,
        dy,
        letter,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
        color=INK,
    )
    if title:
        ax.set_title(title, loc="left", pad=8)


def direct_label(
    ax,
    x,
    y,
    text: str,
    color: str,
    dx: float = 6,
    dy: float = 0,
    ha: str = "left",
    va: str = "center",
    fontsize: float = 8.5,
    weight: str = "semibold",
) -> None:
    """Label a series at its own end, in ITS OWN colour, offset in points.

    This is the relief mechanism for the three palette slots below 3:1 against the
    surface, and it is also what lets a reader identify a line without crossing to a
    legend box. Text elsewhere wears ink tokens; a series' own end label is the one
    place colour and text coincide, because the label IS the series' identity.
    """
    ax.annotate(
        text,
        xy=(x, y),
        xytext=(dx, dy),
        textcoords="offset points",
        color=color,
        fontsize=fontsize,
        fontweight=weight,
        ha=ha,
        va=va,
        zorder=6,
    )


def wilson(k: int, n: int, z: float = 1.959963985) -> tuple[float, float, float]:
    """Wilson score interval for a binomial proportion -> ``(point, lo, hi)``.

    Wilson rather than the normal approximation because coverage and power both live
    near 1.0 here, where the normal interval runs past 100% and stops being readable.
    """
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    lo = max(0.0, centre - half)
    hi = min(1.0, centre + half)
    # The Wilson interval is centred on a SHRUNK estimate, so at p = 1 with small n the
    # clamped upper bound lands a hair below p itself (p=1, n=10 gives hi=0.99977).
    # That is correct as an interval and wrong as an error bar -- matplotlib rejects a
    # negative yerr, which is how it surfaces. An interval must contain its own point
    # estimate; clamp rather than letting each caller discover this separately.
    return p, min(lo, p), max(hi, p)


def binom_summary(flags) -> tuple[float, float, float, int]:
    """``(point, lo, hi, n)`` for a boolean-ish series. NaN counts as False."""
    import pandas as pd

    s = pd.Series(flags)
    n = int(s.notna().sum() + s.isna().sum())  # every row is an opportunity
    k = int(s.fillna(False).astype(bool).sum())
    p, lo, hi = wilson(k, n)
    return p, lo, hi, n


def legend_swatches(
    ax, mapping: dict[str, str], title: str | None = None, marker: str = "patch", **kw
):
    """A legend built from an explicit ``{label: colour}`` map, in insertion order.

    Built by hand rather than from artists because several figures draw a series in
    more than one artist (a box plus its jittered points), and matplotlib's automatic
    legend then either duplicates the entry or picks the wrong one.
    """
    if marker == "patch":
        handles = [Patch(facecolor=c, edgecolor="none", label=k) for k, c in mapping.items()]
    else:
        handles = [
            Line2D([], [], color=c, marker=marker, ls="none", label=k, markersize=7)
            for k, c in mapping.items()
        ]
    return ax.legend(handles=handles, title=title, **kw)


def boxplot(
    ax,
    data,
    positions=None,
    colors=None,
    widths=0.62,
    points: bool = True,
    point_alpha: float = 0.16,
    seed: int = 0,
    showmeans: bool = False,
):
    """House boxplot: filled boxes, ink median, muted whiskers, jittered raw points.

    Outliers are suppressed and the raw points drawn instead. A box plot's flier dots
    and a strip plot's points look identical but mean different things, and showing
    both is the ambiguity; showing the strip alone is honest about the distribution's
    real shape, which matters here because credible-set sizes are heavily skewed.
    """
    positions = list(range(1, len(data) + 1)) if positions is None else list(positions)
    colors = colors or ordinal_colors(len(data))
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=widths,
        patch_artist=True,
        showfliers=False,
        showmeans=showmeans,
        medianprops=dict(color=INK, lw=1.6),
        whiskerprops=dict(color=MUTED, lw=1.0),
        capprops=dict(color=MUTED, lw=1.0),
        boxprops=dict(edgecolor=SURFACE, lw=1.4),  # 2px surface gap between fills
        meanprops=dict(marker="D", markerfacecolor=INK, markeredgecolor="none", ms=4),
    )
    for patch, c in zip(bp["boxes"], colors, strict=False):
        patch.set_facecolor(c)
        patch.set_alpha(0.9)
    if points:
        rng = np.random.default_rng(seed)
        for pos, v in zip(positions, data, strict=False):
            v = np.asarray(v, dtype=float)
            v = v[np.isfinite(v)]
            if not len(v):
                continue
            # Cap the strip so a 4,000-point cell does not become a solid block that
            # hides the box it is meant to qualify.
            if len(v) > 800:
                v = rng.choice(v, 800, replace=False)
            ax.scatter(
                pos + rng.uniform(-0.15, 0.15, size=len(v)),
                v,
                s=6,
                color=INK,
                alpha=point_alpha,
                lw=0,
                zorder=3,
            )
    return bp


def violin_box(ax, data, positions=None, colors=None, widths=0.8, showmean=True):
    """Violin with the box drawn inside it -- the SuSiEx paper's Figure 4a/4b form.

    The violin carries the shape, the box the quartiles, a red dot the mean. Used where
    the distribution's *shape* is the finding (bimodality in a PIP distribution is the
    whole point) rather than only its centre and spread.
    """
    positions = list(range(1, len(data) + 1)) if positions is None else list(positions)
    colors = colors or ordinal_colors(len(data))
    clean = [np.asarray(d, float)[np.isfinite(np.asarray(d, float))] for d in data]
    keep = [(p, d, c) for p, d, c in zip(positions, clean, colors, strict=False) if len(d) > 1]
    if keep:
        vp = ax.violinplot(
            [d for _, d, _ in keep],
            positions=[p for p, _, _ in keep],
            widths=widths,
            showextrema=False,
            showmedians=False,
        )
        for body, (_, _, c) in zip(vp["bodies"], keep, strict=False):
            body.set_facecolor(c)
            body.set_alpha(0.35)
            body.set_edgecolor(c)
            body.set_linewidth(1.0)
    ax.boxplot(
        clean,
        positions=positions,
        widths=widths * 0.22,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color=SURFACE, lw=1.4),
        whiskerprops=dict(color=INK2, lw=1.0),
        capprops=dict(color="none"),
        boxprops=dict(facecolor=INK2, edgecolor="none"),
    )
    if showmean:
        for p, d in zip(positions, clean, strict=False):
            if len(d):
                ax.scatter(
                    [p],
                    [d.mean()],
                    s=22,
                    color=STATUS["critical"],
                    zorder=6,
                    lw=0.8,
                    edgecolor=SURFACE,
                )


def footnote(fig, text: str, y: float = -0.02) -> None:
    """Provenance under a figure: what n, what filter, what caveat.

    WRAPPED TO THE FIGURE'S OWN WIDTH, which is not fussiness. ``savefig`` runs with
    ``bbox_inches="tight"``, so the saved canvas grows to contain every artist -- and a
    single unwrapped sentence of provenance is wider than any of these figures, which
    silently stretches the output to two or three times its intended aspect and leaves
    the plot marooned in the corner. Wrapping first is what keeps the panel geometry
    the one the figsize asked for.
    """
    import textwrap

    # ~= characters that fit across the figure at 8pt in this face. Empirical, and only
    # needs to be close: being slightly narrow costs a line, being wide costs the aspect.
    width = max(60, int(fig.get_figwidth() * 15.5))
    wrapped = "\n".join(textwrap.wrap(" ".join(text.split()), width=width))
    fig.text(0.0, y, wrapped, ha="left", va="top", fontsize=8, color=MUTED, linespacing=1.45)


def spread_labels(
    ax, items, x, min_gap_frac: float = 0.055, dx: float = 8, fontsize: float = 8.5
) -> None:
    """Direct-label several series at a shared x, nudged apart so none is unreadable.

    Series that converge -- which is what these do, since power saturates -- put their
    end labels within a few pixels of each other and the text overprints. Labels are
    laid out in data order and pushed up only as far as the minimum gap requires, so a
    label still sits nearest its own line and the reading stays honest.

    ``items`` is ``[(y, text, colour), ...]``; ``min_gap_frac`` is a fraction of the
    current y-range.
    """
    lo, hi = ax.get_ylim()
    gap = (hi - lo) * min_gap_frac
    ordered = sorted(items, key=lambda it: it[0])
    placed: list[float] = []
    for y, _, _ in ordered:
        y = max(y, placed[-1] + gap) if placed else y
        placed.append(y)
    # Re-centre the run so the nudge is shared rather than pushing the top label off
    # the axes, which is what a purely upward pass does when every series converges.
    overshoot = placed[-1] - hi + gap * 0.5
    if overshoot > 0:
        placed = [p - overshoot for p in placed]
    for (y, text, color), yy in zip(ordered, placed, strict=False):
        # Anchor the text at the NUDGED position, and draw a leader back to the series'
        # true end whenever the nudge is big enough to be noticed. Without the leader a
        # spread label silently misreports where its line finished.
        if abs(yy - y) > gap * 0.25:
            ax.annotate(
                text,
                xy=(x, y),
                xytext=(x, yy),
                xycoords="data",
                textcoords="data",
                color=color,
                fontsize=fontsize,
                fontweight="semibold",
                ha="left",
                va="center",
                zorder=6,
                annotation_clip=False,
                arrowprops=dict(
                    arrowstyle="-", color=color, lw=0.8, shrinkA=1, shrinkB=3, alpha=0.7
                ),
            )
            # annotate() places the text AT xytext; nudge it clear of the leader.
            ax.texts[-1].set_position((x, yy))
            ax.texts[-1].set_ha("left")
        else:
            ax.annotate(
                text,
                xy=(x, y),
                xytext=(dx, 0),
                textcoords="offset points",
                color=color,
                fontsize=fontsize,
                fontweight="semibold",
                ha="left",
                va="center",
                zorder=6,
                annotation_clip=False,
            )


def count_labels(ax, positions, counts, pad_pt: float = 24, fmt: str = "n={:,}") -> None:
    """Group sizes below the tick labels, offset in POINTS rather than axes fraction.

    An axes-fraction offset has to guess how tall the tick labels are, and guesses
    wrong whenever they wrap, rotate, or carry a superscript -- which is how these
    collide. Anchoring to the axes floor in data-x and stepping down a fixed number of
    points puts them under the labels whatever the labels turn out to be. The caller
    still has to make room: pass ``labelpad`` to ``set_xlabel``.
    """
    for pos, n in zip(positions, counts, strict=False):
        ax.annotate(
            fmt.format(int(n)),
            xy=(pos, 0),
            xycoords=("data", "axes fraction"),
            xytext=(0, -pad_pt),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=7.5,
            color=MUTED,
            annotation_clip=False,
        )


def title_with_letter(ax, letter: str, title: str, dx: float = -0.10, dy: float = 1.03) -> None:
    """Panel letter and title on ONE line above the axes, letter first.

    :func:`panel_letter` puts the title through ``set_title``, which sits below
    anything already anchored above the axes -- a legend placed at the top, typically --
    and then the letter and the title stack in the wrong order. This draws both as a
    single text anchored above-left, so a top legend cannot get between them.
    """
    ax.text(
        dx,
        dy,
        letter,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
        color=INK,
    )
    ax.text(
        dx + 0.045,
        dy,
        title,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="semibold",
        va="bottom",
        ha="left",
        color=INK,
    )


def no_data(ax, msg: str = "no data for this panel") -> None:
    """Say so in the panel rather than drawing empty axes that read as a real zero."""
    ax.text(
        0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", color=MUTED, fontsize=9.5
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)


def save(fig, path: Path, log=None) -> Path:
    """Write a figure and close it. Closing matters: thirty figures in one process
    otherwise accumulate and matplotlib starts warning about open figures."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    if log:
        log.info("wrote %s", path)
    return path
