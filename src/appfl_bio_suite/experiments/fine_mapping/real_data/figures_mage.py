"""Publication-quality figures for real-data MAGE + 1000G FedFM validation."""

from __future__ import annotations

import json
import logging
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from appfl_bio_suite.experiments.fine_mapping.figures import figstyle as fs
from appfl_bio_suite.experiments.fine_mapping.figures.figstyle import plt

from .benchmark import extract_locus_regional_profiles, normalize_variant_id
from .federation_builder import SITE_ORDER
from .mage import POPULATION_ORDER
from .panel import STRATUM_ORDER

log = logging.getLogger(__name__)

# Map 1000G/MAGE labels onto the shared visual system ("SAS" shares South Asian "CSA" hue).
MAGE_ANCESTRY_COLORS: dict[str, str] = {
    "EUR": fs.ANCESTRY["EUR"],
    "AFR": fs.ANCESTRY["AFR"],
    "AMR": fs.ANCESTRY["AMR"],
    "EAS": fs.ANCESTRY["EAS"],
    "SAS": fs.ANCESTRY["CSA"],
}

SITE_LABELS: dict[str, str] = {
    "site_atlantic": "Atlantic Biobank",
    "site_africa_eurasia": "Africa-Eurasia Biobank",
    "site_asiapacific": "Asia-Pacific Biobank",
}

STRATUM_LABELS: dict[str, str] = {
    "page_coloc": "PAGE GWAS coloc.",
    "high_fst": "High Fst",
    "multi_cs": "Multi-CS (>=4)",
    "moderate_cs": "Moderate single-CS",
}

STRATUM_COLORS: dict[str, str] = {
    "page_coloc": fs.BLUE_RAMP[2],
    "high_fst": fs.RED_RAMP[1],
    "multi_cs": fs.GREEN_RAMP[2],
    "moderate_cs": fs.MUTED,
}

ARM_LABELS: dict[str, str] = {
    "federated_susiex": "Federated SuSiEx (5-anc)",
    "centralized_susiex": "Pooled SuSiEx (5-anc)",
    "federated_multisusie": "Federated MultiSuSiE (5-anc)",
    "centralized_multisusie": "Pooled MultiSuSiE (5-anc)",
    "site_atlantic": "Site: Atlantic only",
    "site_africa_eurasia": "Site: Africa-Eurasia only",
    "site_asiapacific": "Site: Asia-Pacific only",
    "anc_EUR": "Single-anc: EUR",
    "anc_AFR": "Single-anc: AFR",
    "anc_AMR": "Single-anc: AMR",
    "anc_EAS": "Single-anc: EAS",
    "anc_SAS": "Single-anc: SAS",
}


def _save_pair(fig, out_stem: Path) -> list[Path]:
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = out_stem.with_suffix(".png")
    pdf_path = out_stem.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return [png_path, pdf_path]


def fig1_cohort_and_design(federation_root: Path, reports_dir: Path, out_stem: Path) -> list[Path]:
    """Fig 1: 3-site MAGE cohort structure, swap audit, and 16-locus panel overview."""
    fs.apply_style()
    meta_dir = Path(federation_root) / "meta"
    cohort = pd.read_csv(meta_dir / "cohort_metadata.tsv", sep="\t")
    panel = pd.read_csv(meta_dir / "panel_manifest.tsv", sep="\t")
    audit = json.loads((meta_dir / "sample_swap_audit.json").read_text(encoding="utf-8"))
    parity_path = Path(reports_dir) / "parity_summary.tsv"
    parity_df = pd.read_csv(parity_path, sep="\t") if parity_path.is_file() else pd.DataFrame()

    fig = plt.figure(figsize=(15.4, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.05, 1.45], wspace=0.34)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    # Panel a: stacked site x continental ancestry bars
    sites_present = [s for s in SITE_ORDER if s in set(cohort["site_id"])]
    counts = (
        cohort.groupby(["site_id", "continentalGroup"])
        .size()
        .unstack(fill_value=0)
        .reindex(sites_present, fill_value=0)
    )
    y_pos = np.arange(len(sites_present))
    left = np.zeros(len(sites_present), dtype=float)
    for pop in POPULATION_ORDER:
        vals = (
            counts[pop].to_numpy(dtype=float)
            if pop in counts.columns
            else np.zeros(len(sites_present))
        )
        ax_a.barh(
            y_pos,
            vals,
            left=left,
            height=0.56,
            color=MAGE_ANCESTRY_COLORS[pop],
            edgecolor=fs.SURFACE,
            linewidth=1.0,
            label=pop,
        )
        for i, v in enumerate(vals):
            if v >= 22:
                ax_a.text(
                    left[i] + v / 2.0,
                    y_pos[i],
                    f"{int(v)}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    fontweight="semibold",
                    color=fs.SURFACE,
                )
        left += vals

    n_pops_by_site = cohort.groupby("site_id")["population"].nunique()
    for i, site in enumerate(sites_present):
        ax_a.text(
            left[i] + 5,
            y_pos[i],
            f"n={int(left[i])} ({int(n_pops_by_site.get(site, 0))} pops)",
            ha="left",
            va="center",
            fontsize=8.5,
            fontweight="semibold",
            color=fs.INK,
        )
    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels([SITE_LABELS.get(s, s) for s in sites_present])
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, max(left.max() * 1.38, 100))
    ax_a.set_xlabel(f"Enrolled individuals (clean cohort, N = {len(cohort):,})")
    ax_a.grid(axis="y", visible=False)
    ax_a.legend(
        loc="upper center",
        bbox_to_anchor=(0.45, -0.18),
        ncol=5,
        frameon=False,
        fontsize=8.5,
        handlelength=1.1,
        columnspacing=0.9,
    )
    fs.title_with_letter(ax_a, "a", "Three-site 1000G / MAGE partition", dx=-0.22, dy=1.05)

    # Panel b: Sample-swap audit card
    ax_b.axis("off")
    swap_rec = audit.get("swapped_record", {})
    true_rec = audit.get("true_donor_record", {})
    swap_row = (
        parity_df[parity_df["comparison"] == "sample_swap_with_vs_clean"].iloc[0]
        if not parity_df.empty and "sample_swap_with_vs_clean" in set(parity_df["comparison"])
        else None
    )
    lines = [
        (
            "Raw Zenodo MAGE v1.0",
            f"N = {audit.get('n_raw_samples', 731)} expression columns",
            fs.INK,
        ),
        (
            "Detected duplicate donor:",
            f"{swap_rec.get('sample_kgpID', 'NA11919')} "
            f"({swap_rec.get('sra_accession', 'SRR19762653')})\n"
            f"batch {swap_rec.get('batch', 3)}, {swap_rec.get('annotated_population', 'CEU')}, "
            f"{swap_rec.get('annotated_sex', 'XY')}\n"
            f"-> true donor: {true_rec.get('sample_kgpID', 'HG00237')}\n"
            f"   ({true_rec.get('retained_sra_accession', 'SRR19762247')}, "
            f"batch {true_rec.get('retained_batch', 7)})",
            fs.STATUS["critical"],
        ),
        (
            "Remediation applied:",
            f"Exclude {','.join(audit.get('excluded_samples', [])) or 'NA11919'}\n"
            f"Clean cohort: N = {audit.get('n_clean_samples', len(cohort))} unique\n"
            f"donors across {audit.get('n_populations', 26)} populations",
            fs.STATUS["good"],
        ),
    ]
    if swap_row is not None:
        lines.append(
            (
                "Swap sensitivity on 16 loci:",
                f"max |ΔPIP| = {swap_row['worst_max_pip_abs_diff']:.4g}\n"
                f"CS Jaccard = {swap_row['mean_cs_union_jaccard']:.4f} | "
                f"Top SNP = {swap_row['frac_top_snp_identical'] * 100:.0f}%",
                fs.INK2,
            )
        )
    ax_b.set_xlim(0.0, 1.0)
    ax_b.set_ylim(0.0, 1.0)
    y_cursor = 0.96
    for hdr, body, col in lines:
        ax_b.text(
            0.0,
            y_cursor,
            hdr,
            fontsize=8.6,
            fontweight="bold",
            color=col,
            va="top",
            transform=ax_b.transAxes,
        )
        ax_b.text(
            0.03,
            y_cursor - 0.058,
            body,
            fontsize=7.6,
            color=fs.INK2,
            va="top",
            linespacing=1.24,
            transform=ax_b.transAxes,
        )
        n_body_lines = body.count("\n") + 1
        y_cursor -= 0.115 + 0.048 * n_body_lines
    ax_b.text(
        0.0,
        1.05,
        "b",
        fontsize=12,
        fontweight="bold",
        color=fs.INK,
        transform=ax_b.transAxes,
        va="bottom",
    )
    ax_b.text(
        0.08,
        1.05,
        "MAGE sample-swap audit",
        fontsize=11,
        fontweight="semibold",
        color=fs.INK,
        transform=ax_b.transAxes,
        va="bottom",
    )

    # Panel c: Selected 16-locus validation panel across 4 strata
    p_ord = panel.copy()
    p_ord["stratum_rank"] = p_ord["stratum"].map({s: i for i, s in enumerate(STRATUM_ORDER)})
    p_ord = p_ord.sort_values(
        ["stratum_rank", "chrom", "tss"], ascending=[True, True, True]
    ).reset_index(drop=True)
    yy = np.arange(len(p_ord))
    cs_sizes = p_ord["min_cs_size"].to_numpy(dtype=float)
    colors_c = [STRATUM_COLORS.get(s, fs.MUTED) for s in p_ord["stratum"]]
    ax_c.hlines(yy, 1.0, cs_sizes, color=fs.GRID, lw=1.2, zorder=2)
    ax_c.scatter(
        cs_sizes,
        yy,
        s=32 + 16 * np.clip(p_ord["n_published_cs"].to_numpy(dtype=float), 1, 6),
        c=colors_c,
        edgecolor=fs.SURFACE,
        linewidth=0.8,
        zorder=4,
    )
    for idx, row in p_ord.iterrows():
        tag = f"CS={int(row['n_published_cs'])}"
        if row["page_max_pph4"] >= 0.85:
            tag += f", H4={row['page_max_pph4']:.2f}"
        elif row["max_wc_fst"] >= 0.45:
            tag += f", Fst={row['max_wc_fst']:.2f}"
        ax_c.text(
            max(float(row["min_cs_size"]) * 1.16, 1.45),
            idx,
            tag,
            va="center",
            ha="left",
            fontsize=7.3,
            color=fs.INK2,
        )
    ax_c.set_xscale("log")
    ax_c.set_xlim(0.85, max(cs_sizes.max() * 4.5, 40))
    ax_c.set_yticks(yy)
    ax_c.set_yticklabels(
        [f"{r.gene_symbol} (chr{r.chrom})" for r in p_ord.itertuples()], fontsize=8
    )
    ax_c.invert_yaxis()
    ax_c.set_xlabel("Published MAGE minimum CS size (variants, log scale)")
    ax_c.grid(axis="y", visible=False)
    fs.legend_swatches(
        ax_c,
        {STRATUM_LABELS[s]: STRATUM_COLORS[s] for s in STRATUM_ORDER if s in set(p_ord["stratum"])},
        marker="o",
        loc="upper center",
        bbox_to_anchor=(0.45, -0.16),
        ncol=2,
        fontsize=8,
    )
    fs.title_with_letter(ax_c, "c", "16-locus validation panel across 4 strata", dx=-0.26, dy=1.05)

    fs.footnote(
        fig,
        "Clean MAGE v1.0 cohort (N=730 distinct 1000G NYGC 30x donors across 26 populations) "
        "partitioned into 3 geographically structured biobanks so EUR/AFR/AMR/SAS span sites "
        "while EAS concentrates at Asia-Pacific. Marker size in c scales with published MAGE "
        "credible-set count (+/-250 kb TSS cis-windows, excluding extended MHC).",
        y=-0.16,
    )
    return _save_pair(fig, out_stem)


def _collect_variant_pips_pair(
    ref_dir: Path, tst_dir: Path
) -> tuple[np.ndarray, np.ndarray, list[float]]:
    ref_insts = {
        p.parent.name: p.parent for p in (ref_dir / "data" / "instances").glob("*/variants.tsv")
    }
    tst_insts = {
        p.parent.name: p.parent for p in (tst_dir / "data" / "instances").glob("*/variants.tsv")
    }
    xs, ys, max_diffs = [], [], []
    for inst in sorted(set(ref_insts) & set(tst_insts)):
        va = pd.read_csv(
            ref_insts[inst] / "variants.tsv", sep="\t", usecols=["snp_id", "pip"]
        ).dropna()
        vb = pd.read_csv(
            tst_insts[inst] / "variants.tsv", sep="\t", usecols=["snp_id", "pip"]
        ).dropna()
        m = va.merge(vb, on="snp_id", suffixes=("_r", "_t"))
        if not m.empty:
            xr = m["pip_r"].to_numpy(dtype=float)
            yt = m["pip_t"].to_numpy(dtype=float)
            xs.append(xr)
            ys.append(yt)
            max_diffs.append(float(np.max(np.abs(xr - yt))))
    if not xs:
        return np.array([]), np.array([]), []
    return np.concatenate(xs), np.concatenate(ys), max_diffs


def _plot_pip_scatter(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    color: str,
    xlabel: str,
    ylabel: str,
    letter: str,
    title: str,
) -> None:
    if len(x) == 0:
        fs.no_data(ax, "no comparable variants")
        fs.title_with_letter(ax, letter, title)
        return
    ax.plot([0, 1], [0, 1], ls="--", lw=1.2, color=fs.INK2, zorder=2)
    # Downsample near-zero background points while plotting every informative PIP > 0.01
    hi = (x > 0.01) | (y > 0.01)
    rng = np.random.default_rng(7)
    lo_idx = np.flatnonzero(~hi)
    if len(lo_idx) > 2500:
        lo_idx = rng.choice(lo_idx, 2500, replace=False)
    keep_idx = np.concatenate([np.flatnonzero(hi), lo_idx])
    ax.scatter(
        x[keep_idx],
        y[keep_idx],
        s=18,
        color=color,
        alpha=0.58,
        lw=0.4,
        edgecolor=fs.SURFACE,
        zorder=4,
    )
    r = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else 1.0
    mdiff = float(np.max(np.abs(x - y)))
    ax.text(
        0.04,
        0.95,
        f"n = {len(x):,} variants\nr = {r:.6f}\nmax |ΔPIP| = {mdiff:.2e}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.2,
        color=fs.INK,
        bbox=dict(boxstyle="round,pad=0.3", facecolor=fs.SURFACE, edgecolor=fs.GRID, lw=0.8),
    )
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fs.title_with_letter(ax, letter, title, dx=-0.16, dy=1.05)


def fig2_federated_parity_and_controls(
    arms_dir: Path,
    reports_dir: Path,
    out_stem: Path,
) -> list[Path]:
    """Fig 2: Exact Federated vs Pooled parity contrasted against controls."""
    fs.apply_style()
    arms_dir = Path(arms_dir)
    reports_dir = Path(reports_dir)
    paired_df = pd.read_csv(reports_dir / "paired_arm_comparisons.tsv", sep="\t")
    parity_df = pd.read_csv(reports_dir / "parity_summary.tsv", sep="\t")
    ld_dist_path = reports_dir / "naive_ld_distortion.tsv"
    ld_dist = pd.read_csv(ld_dist_path, sep="\t") if ld_dist_path.is_file() else pd.DataFrame()

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.6))
    plt.subplots_adjust(wspace=0.32, hspace=0.36)

    # a: SuSiEx Federated vs Centralized
    x_sx, y_sx, _ = _collect_variant_pips_pair(
        arms_dir / "centralized_susiex", arms_dir / "federated_susiex"
    )
    _plot_pip_scatter(
        axes[0, 0],
        x_sx,
        y_sx,
        fs.PATH["federated"],
        "Centralized Pooled SuSiEx PIP",
        "Federated SuSiEx PIP",
        "a",
        "SuSiEx exact parity",
    )

    # b: MultiSuSiE Federated vs Centralized
    x_ms, y_ms, _ = _collect_variant_pips_pair(
        arms_dir / "centralized_multisusie", arms_dir / "federated_multisusie"
    )
    _plot_pip_scatter(
        axes[0, 1],
        x_ms,
        y_ms,
        fs.BLUE_RAMP[2],
        "Centralized Pooled MultiSuSiE PIP",
        "Federated MultiSuSiE PIP",
        "b",
        "MultiSuSiE exact parity",
    )

    # c: Naive Site-Local Standardization vs Centralized
    x_nv, y_nv, _ = _collect_variant_pips_pair(
        arms_dir / "centralized_susiex", arms_dir / "naive_site_standardized"
    )
    _plot_pip_scatter(
        axes[0, 2],
        x_nv,
        y_nv,
        fs.STATUS["critical"],
        "Centralized Pooled SuSiEx PIP",
        "Naive site-standardized PIP",
        "c",
        "Failure mode: site-local standardization",
    )
    if not ld_dist.empty:
        max_ld = ld_dist["max_ld_abs_diff"].max()
        max_z = ld_dist["max_z_abs_diff"].max()
        axes[0, 2].text(
            0.96,
            0.05,
            f"max |ΔLD| = {max_ld:.3f}\nmax |ΔZ| = {max_z:.2f}",
            transform=axes[0, 2].transAxes,
            va="bottom",
            ha="right",
            fontsize=8.0,
            color=fs.STATUS["critical"],
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor=fs.SURFACE,
                edgecolor=fs.STATUS["critical"],
                lw=0.8,
            ),
        )

    # d: Log-scale per-locus max |ΔPIP| across all methodological comparisons
    comp_order = [
        ("susiex_fed_vs_central", "SuSiEx\nFed vs\nPool", fs.PATH["federated"]),
        ("multisusie_fed_vs_central", "Multi-\nSuSiE\nFed/Pool", fs.BLUE_RAMP[2]),
        ("sample_swap_with_vs_clean", "Swap\nN=731\nvs 730", fs.GREEN_RAMP[1]),
        ("susiex_fwl_vs_yresid", "Covar.\nFWL vs\nY_resid", fs.INK2),
        ("susiex_naive_std_vs_central", "Naive\nstd. vs\nPool", fs.STATUS["critical"]),
    ]
    present_comps = [c for c in comp_order if c[0] in set(paired_df["comparison"])]
    data_d = [
        np.clip(
            paired_df.loc[paired_df["comparison"] == key, "max_pip_abs_diff"].to_numpy(dtype=float),
            1e-9,
            1.0,
        )
        for key, _, _ in present_comps
    ]
    fs.boxplot(
        axes[1, 0],
        data_d,
        colors=[c[2] for c in present_comps],
        points=True,
        point_alpha=0.55,
    )
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_ylim(3e-10, 2.5)
    axes[1, 0].set_xticks(range(1, len(present_comps) + 1))
    axes[1, 0].set_xticklabels([c[1] for c in present_comps], fontsize=7.2, linespacing=1.12)
    axes[1, 0].set_ylabel("Per-locus max |ΔPIP| (log scale)")
    fs.title_with_letter(axes[1, 0], "d", "Max PIP divergence per locus", dx=-0.16, dy=1.05)

    # e: FWL (X_perp, Y_perp) vs Y_resid + Raw X scatter
    x_fwl, y_fwl, _ = _collect_variant_pips_pair(
        arms_dir / "centralized_susiex_fwl", arms_dir / "federated_susiex"
    )
    _plot_pip_scatter(
        axes[1, 1],
        x_fwl,
        y_fwl,
        fs.INK2,
        "Exact FWL (X_perp, Y_perp) PIP",
        "Standard FedFM (Raw X, Y_resid) PIP",
        "e",
        "Covariate orthogonalization check",
    )

    # f: CS concordance across exactness & sensitivity checks
    if present_comps:
        pmap = parity_df.set_index("comparison")
        xpos = np.arange(len(present_comps))
        jacc = [float(pmap.loc[k, "mean_cs_union_jaccard"]) for k, _, _ in present_comps]
        exact = [float(pmap.loc[k, "frac_exact_cs_identical"]) for k, _, _ in present_comps]
        w = 0.40
        b1 = axes[1, 2].bar(
            xpos - w / 2, jacc, width=w, color=fs.BLUE_RAMP[2], label="Mean CS Jaccard"
        )
        b2 = axes[1, 2].bar(
            xpos + w / 2, exact, width=w, color=fs.GREEN_RAMP[2], label="Exact CS match rate"
        )
        for bar in [*b1, *b2]:
            h = bar.get_height()
            axes[1, 2].text(
                bar.get_x() + bar.get_width() / 2,
                h + 0.02,
                f"{h:.2f}",
                ha="center",
                va="bottom",
                fontsize=6.5,
                color=fs.INK,
            )
        axes[1, 2].set_xticks(xpos)
        axes[1, 2].set_xticklabels([c[1] for c in present_comps], fontsize=7.2, linespacing=1.12)
        axes[1, 2].set_ylim(0, 1.15)
        axes[1, 2].set_ylabel("Credible-set concordance")
        axes[1, 2].legend(
            loc="lower left",
            fontsize=8,
            frameon=True,
            facecolor=fs.SURFACE,
            edgecolor=fs.GRID,
        )
    fs.title_with_letter(axes[1, 2], "f", "CS agreement across controls", dx=-0.16, dy=1.05)

    fs.footnote(
        fig,
        "Comparison across all 16 autosomal MAGE loci (MAF >= 0.005). Exact federated pooling "
        "(a, b) matches pooled genotypes to machine/file precision whereas site-local "
        "standardization (c) distorts within-ancestry LD and marginal Z-scores before inference.",
        y=-0.04,
    )
    return _save_pair(fig, out_stem)


def _annotate_scatter_loci(
    ax,
    points: list[tuple[float, float, str]],
    placements: dict[str, tuple[float, float, str, str]],
    cluster_tol: float = 0.025,
) -> None:
    """Annotate scatter points cleanly, grouping nearly coincident loci with leader lines."""
    clusters: list[list[tuple[float, float, str]]] = []
    for pt in points:
        x, y, _ = pt
        placed = False
        for cl in clusters:
            cx = float(np.mean([p[0] for p in cl]))
            cy = float(np.mean([p[1] for p in cl]))
            if float(np.hypot(x - cx, y - cy)) <= cluster_tol:
                cl.append(pt)
                placed = True
                break
        if not placed:
            clusters.append([pt])

    for cl in clusters:
        cx = float(np.mean([p[0] for p in cl]))
        cy = float(np.mean([p[1] for p in cl]))
        genes = sorted(p[2] for p in cl)
        pos = next((placements[g] for g in genes if g in placements), None)
        if len(genes) >= 3:
            lines = [", ".join(genes[i : i + 2]) for i in range(0, len(genes), 2)]
            label = ",\n".join(lines)
        else:
            label = ", ".join(genes)
        if pos is None:
            ax.text(cx + 0.015, cy - 0.015, label, fontsize=7.3, color=fs.INK2)
            continue
        tx, ty, ha, va = pos
        need_arrow = float(np.hypot(tx - cx, ty - cy)) > 0.045
        ax.annotate(
            label,
            xy=(cx, cy),
            xytext=(tx, ty),
            ha=ha,
            va=va,
            fontsize=7.3,
            color=fs.INK2,
            linespacing=1.15,
            arrowprops=(
                {"arrowstyle": "-", "color": fs.MUTED, "lw": 0.7, "shrinkA": 2, "shrinkB": 3}
                if need_arrow
                else None
            ),
            zorder=5,
        )


def fig3_resolution_and_arms(reports_dir: Path, out_stem: Path) -> list[Path]:
    """Fig 3: Multi-ancestry federation vs solo biobanks, single ancestries, and MultiSuSiE."""
    fs.apply_style()
    reports_dir = Path(reports_dir)
    arm_loci = pd.read_csv(reports_dir / "arm_locus_metrics.tsv", sep="\t")
    arm_summ = pd.read_csv(reports_dir / "arm_summary.tsv", sep="\t")

    desired_arms = [
        "federated_susiex",
        "federated_multisusie",
        "site_atlantic",
        "site_africa_eurasia",
        "site_asiapacific",
        "anc_AFR",
        "anc_EUR",
        "anc_EAS",
        "anc_SAS",
        "anc_AMR",
    ]
    arms_present = [a for a in desired_arms if a in set(arm_summ["arm_id"])]
    summ = arm_summ.set_index("arm_id").loc[arms_present].reset_index()

    fig, axes = plt.subplots(2, 2, figsize=(14.6, 9.8))
    plt.subplots_adjust(wspace=0.36, hspace=0.48)

    def _arm_color(a: str) -> str:
        if a == "federated_susiex":
            return fs.PATH["federated"]
        if a == "federated_multisusie":
            return fs.BLUE_RAMP[2]
        if a.startswith("anc_"):
            return MAGE_ANCESTRY_COLORS.get(a.removeprefix("anc_"), fs.MUTED)
        return fs.INK2

    # a: Published MAGE CS capture rate across arms
    ax_a = axes[0, 0]
    yy = np.arange(len(summ))
    cap_rate = summ["total_pub_cs_captured"] / summ["total_pub_cs"].replace(0, np.nan)
    cols_a = [_arm_color(a) for a in summ["arm_id"]]
    ax_a.barh(yy, cap_rate, height=0.62, color=cols_a, edgecolor=fs.SURFACE, lw=1.0)
    for i, row in summ.iterrows():
        rate = row["total_pub_cs_captured"] / max(row["total_pub_cs"], 1)
        ax_a.text(
            rate + 0.015,
            i,
            f"{int(row['total_pub_cs_captured'])}/{int(row['total_pub_cs'])} ({rate * 100:.0f}%)",
            va="center",
            ha="left",
            fontsize=8,
            color=fs.INK,
        )
    ax_a.set_yticks(yy)
    ax_a.set_yticklabels([ARM_LABELS.get(a, a) for a in summ["arm_id"]], fontsize=8.5)
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, 1.18)
    ax_a.set_xlabel("Published MAGE credible sets captured")
    ax_a.grid(axis="y", visible=False)
    fs.title_with_letter(ax_a, "a", "Published CS capture across cohorts", dx=-0.34, dy=1.05)

    # b: Minimum CS size distribution per locus (smaller = finer resolution)
    ax_b = axes[0, 1]
    cs_data = []
    cs_counts = []
    for a in arms_present:
        vals = arm_loci.loc[arm_loci["arm_id"] == a, "min_cs_size"].dropna().to_numpy(dtype=float)
        cs_data.append(vals if len(vals) else np.array([np.nan]))
        cs_counts.append(len(vals))
    fs.boxplot(
        ax_b,
        cs_data,
        colors=cols_a,
        points=True,
        point_alpha=0.55,
    )
    ax_b.set_yscale("log")
    finite_vals = [vals[np.isfinite(vals)] for vals in cs_data if np.any(np.isfinite(vals))]
    max_cs = max((float(np.max(v)) for v in finite_vals if len(v) > 0), default=100.0)
    ax_b.set_ylim(0.75, max(max_cs * 3.2, 20.0))
    ax_b.set_xticks(range(1, len(arms_present) + 1))
    arm_short_labels = [
        ARM_LABELS.get(a, a)
        .replace("Federated ", "Fed ")
        .replace("Single-anc: ", "")
        .replace("Site: ", "")
        for a in arms_present
    ]
    ax_b.set_xticklabels(arm_short_labels, rotation=35, ha="right", fontsize=8)
    ax_b.set_ylabel("Minimum CS size per locus (log)")
    for pos, n in zip(range(1, len(arms_present) + 1), cs_counts, strict=True):
        ax_b.annotate(
            f"n={int(n)}",
            xy=(pos, 0.97),
            xycoords=("data", "axes fraction"),
            ha="center",
            va="top",
            fontsize=7.2,
            color=fs.MUTED,
        )
    fs.title_with_letter(ax_b, "b", "Credible-set resolution by cohort", dx=-0.16, dy=1.05)

    # c: Tight CSs (<=5 variants) and singleton CSs + PIP > 0.95 variants
    ax_c = axes[1, 0]
    xpos = np.arange(len(arms_present))
    w = 0.27
    b_sing = ax_c.bar(
        xpos - w, summ["total_singleton_cs"], width=w, color=fs.GREEN_RAMP[2], label="Singleton CS"
    )
    b_le5 = ax_c.bar(
        xpos, summ["total_cs_le5"], width=w, color=fs.PATH["federated"], label="Tight CS (<=5)"
    )
    ax_c.bar(
        xpos + w,
        summ["total_pip95_variants"],
        width=w,
        color=fs.BLUE_RAMP[1],
        label="Variants PIP > 0.95",
    )
    for b in [*b_sing, *b_le5]:
        h = b.get_height()
        if h > 0:
            ax_c.text(
                b.get_x() + b.get_width() / 2,
                h + 0.3,
                f"{int(h)}",
                ha="center",
                va="bottom",
                fontsize=7.2,
                color=fs.INK2,
            )
    ax_c.set_xticks(xpos)
    ax_c.set_xticklabels(arm_short_labels, rotation=35, ha="right", fontsize=8)
    ax_c.set_ylabel("Count across 16 MAGE loci")
    ax_c.legend(loc="upper right", fontsize=8)
    fs.title_with_letter(ax_c, "c", "High-confidence & singleton yield", dx=-0.16, dy=1.05)

    # d: Federated SuSiEx vs Federated MultiSuSiE per-locus top PIP comparison
    ax_d = axes[1, 1]
    sx = arm_loci[arm_loci["arm_id"] == "federated_susiex"].set_index("locus_id")
    ms = arm_loci[arm_loci["arm_id"] == "federated_multisusie"].set_index("locus_id")
    common_lids = sorted(set(sx.index) & set(ms.index))
    if common_lids:
        xv = sx.loc[common_lids, "top_pip"].to_numpy(dtype=float)
        yv = ms.loc[common_lids, "top_pip"].to_numpy(dtype=float)
        st = sx.loc[common_lids, "stratum"].tolist()
        ax_d.plot([0, 1], [0, 1], ls="--", lw=1.2, color=fs.INK2, zorder=2)
        pts_d: list[tuple[float, float, str]] = []
        for xi, yi, si, lid in zip(xv, yv, st, common_lids, strict=True):
            if np.isfinite(xi) and np.isfinite(yi):
                ax_d.scatter(
                    [xi],
                    [yi],
                    s=58,
                    color=STRATUM_COLORS.get(si, fs.MUTED),
                    edgecolor=fs.SURFACE,
                    lw=0.9,
                    zorder=4,
                )
                gene = str(sx.loc[lid, "gene_symbol"])
                if abs(xi - yi) > 0.12 or (xi > 0.90 and yi > 0.90):
                    pts_d.append((float(xi), float(yi), gene))
        _annotate_scatter_loci(
            ax_d,
            pts_d,
            placements={
                "NFATC1": (0.05, 1.01, "left", "center"),
                "GS1-124K5.4": (0.12, 0.90, "left", "center"),
                "ERN1": (0.06, 0.79, "left", "center"),
                "ILRUN": (0.06, 0.49, "left", "center"),
                "PCMTD1": (-0.02, 0.23, "left", "top"),
                "POT1": (0.15, 0.08, "left", "center"),
                "NPNT": (0.26, 0.42, "left", "bottom"),
                "ACOT1": (0.79, 0.97, "right", "top"),
            },
        )
        fs.legend_swatches(
            ax_d,
            {STRATUM_LABELS[s]: STRATUM_COLORS[s] for s in STRATUM_ORDER if s in set(st)},
            marker="o",
            loc="lower right",
            fontsize=8,
        )
        ax_d.set_xlim(-0.04, 1.06)
        ax_d.set_ylim(-0.04, 1.06)
        ax_d.set_xlabel("Federated SuSiEx locus top PIP")
        ax_d.set_ylabel("Federated MultiSuSiE locus top PIP")
    else:
        fs.no_data(ax_d, "MultiSuSiE results unavailable")
    fs.title_with_letter(ax_d, "d", "SuSiEx vs MultiSuSiE per-locus top PIP", dx=-0.16, dy=1.05)

    fs.footnote(
        fig,
        "Comparison across 16 MAGE autosomal eQTL loci: full 3-site 5-ancestry federation (N=730) "
        "vs solo biobanks and single-ancestry analyses across pooled sites.",
        y=-0.03,
    )
    return _save_pair(fig, out_stem)


def fig4_mage_agreement_and_ancestry(reports_dir: Path, out_stem: Path) -> list[Path]:
    """Fig 4: Agreement with published MAGE findings and ancestry-specific POST-HOC_PROB."""
    fs.apply_style()
    reports_dir = Path(reports_dir)
    arm_loci = pd.read_csv(reports_dir / "arm_locus_metrics.tsv", sep="\t")
    posthoc = pd.read_csv(reports_dir / "ancestry_posthoc_summary.tsv", sep="\t")

    fed_sx = arm_loci[arm_loci["arm_id"] == "federated_susiex"].copy()
    fed_ms = arm_loci[arm_loci["arm_id"] == "federated_multisusie"].copy()

    fig, axes = plt.subplots(2, 2, figsize=(15.0, 9.6))
    plt.subplots_adjust(wspace=0.35, hspace=0.38)

    # a: Published MAGE lead PIP vs Federated SuSiEx PIP on same variant
    ax_a = axes[0, 0]
    ok_a = fed_sx.dropna(subset=["published_lead_pip", "pub_lead_arm_pip"])
    if not ok_a.empty:
        ax_a.plot([0, 1], [0, 1], ls="--", lw=1.2, color=fs.INK2, zorder=2)
        pts_a: list[tuple[float, float, str]] = []
        for _, r in ok_a.iterrows():
            col = STRATUM_COLORS.get(str(r["stratum"]), fs.MUTED)
            ax_a.scatter(
                [r["published_lead_pip"]],
                [r["pub_lead_arm_pip"]],
                s=62,
                color=col,
                edgecolor=fs.SURFACE,
                lw=0.9,
                zorder=4,
            )
            gap = abs(r["published_lead_pip"] - r["pub_lead_arm_pip"])
            if gap > 0.18 or r["published_lead_pip"] > 0.90:
                pts_a.append(
                    (
                        float(r["published_lead_pip"]),
                        float(r["pub_lead_arm_pip"]),
                        str(r["gene_symbol"]),
                    )
                )
        _annotate_scatter_loci(
            ax_a,
            pts_a,
            placements={
                "ARIH1": (0.36, 0.05, "left", "center"),
                "ILRUN": (0.58, 0.08, "center", "bottom"),
                "ERN1": (0.68, 0.15, "center", "bottom"),
                "PCMTD1": (0.76, 0.24, "center", "bottom"),
                "NPNT": (0.83, 0.35, "center", "bottom"),
                "NFATC1": (0.945, 0.24, "center", "bottom"),
                "ACOT1": (1.06, 0.36, "right", "bottom"),
                "CETP": (0.93, 0.95, "right", "top"),
            },
        )
        r_lead = (
            float(np.corrcoef(ok_a["published_lead_pip"], ok_a["pub_lead_arm_pip"])[0, 1])
            if len(ok_a) >= 2
            else np.nan
        )
        in_cs_pct = ok_a["pub_lead_in_cs"].mean() * 100
        ax_a.text(
            0.04,
            0.94,
            f"n = {len(ok_a)} published leads\nr = {r_lead:.3f}\nin-CS rate = {in_cs_pct:.1f}%",
            transform=ax_a.transAxes,
            va="top",
            ha="left",
            fontsize=8.2,
            bbox=dict(boxstyle="round,pad=0.3", facecolor=fs.SURFACE, edgecolor=fs.GRID, lw=0.8),
        )
        ax_a.set_xlim(-0.04, 1.08)
        ax_a.set_ylim(-0.06, 1.08)
        ax_a.set_xlabel("Published MAGE SuSiE PIP (lead variant)")
        ax_a.set_ylabel("Federated 5-ancestry SuSiEx PIP")
        fs.legend_swatches(
            ax_a,
            {STRATUM_LABELS[s]: STRATUM_COLORS[s] for s in STRATUM_ORDER},
            marker="o",
            loc="center left",
            bbox_to_anchor=(0.02, 0.56),
            fontsize=8,
        )
    else:
        fs.no_data(ax_a)
    fs.title_with_letter(ax_a, "a", "Published MAGE lead PIP concordance", dx=-0.16, dy=1.05)

    # b: Minimum CS size by stratum: Published MAGE vs Federated SuSiEx vs Federated MultiSuSiE
    ax_b = axes[0, 1]
    strata = [s for s in STRATUM_ORDER if s in set(fed_sx["stratum"])]
    xpos = np.arange(len(strata))
    w = 0.25
    pub_med = [
        float(fed_sx.loc[fed_sx["stratum"] == s, "published_min_cs_size"].median()) for s in strata
    ]
    sx_med = [float(fed_sx.loc[fed_sx["stratum"] == s, "min_cs_size"].median()) for s in strata]
    ms_med = [
        float(fed_ms.loc[fed_ms["stratum"] == s, "min_cs_size"].median())
        if not fed_ms.empty
        else np.nan
        for s in strata
    ]
    ax_b.bar(xpos - w, pub_med, width=w, color=fs.MUTED, label="Published MAGE (SuSiE)")
    ax_b.bar(xpos, sx_med, width=w, color=fs.PATH["federated"], label="Federated SuSiEx")
    ax_b.bar(xpos + w, ms_med, width=w, color=fs.BLUE_RAMP[2], label="Federated MultiSuSiE")
    ax_b.set_xticks(xpos)
    ax_b.set_xticklabels([STRATUM_LABELS[s] for s in strata], rotation=18, ha="right", fontsize=8.5)
    ax_b.set_ylabel("Median minimum CS size per stratum")
    ax_b.legend(loc="upper left", fontsize=8)
    fs.title_with_letter(ax_b, "b", "CS resolution vs published MAGE by stratum", dx=-0.16, dy=1.05)

    # c: Heatmap of POST-HOC_PROB across EUR, AFR, AMR, EAS, SAS
    ax_c = axes[1, 0]
    pop_cols = [f"prob_{p}" for p in POPULATION_ORDER if f"prob_{p}" in posthoc.columns]
    if not posthoc.empty and pop_cols:
        # Sort by Fst descending and active pops ascending to showcase shared vs differentiated
        ph_show = (
            posthoc.sort_values(
                ["max_wc_fst", "n_active_pops_prob80", "gene_symbol", "cs_id"],
                ascending=[False, True, True, True],
            )
            .head(18)
            .copy()
        )
        mat = ph_show[pop_cols].to_numpy(dtype=float)
        cmap = plt.get_cmap("argonne_sequential")
        im = ax_c.imshow(np.nan_to_num(mat, nan=0.0), aspect="auto", vmin=0.0, vmax=1.0, cmap=cmap)
        for r_i in range(mat.shape[0]):
            for c_i in range(mat.shape[1]):
                v = mat[r_i, c_i]
                txt = "NA" if not np.isfinite(v) else f"{v:.2f}"
                col_txt = fs.SURFACE if (np.isfinite(v) and v >= 0.55) else fs.INK
                weight = "bold" if (np.isfinite(v) and v >= 0.80) else "normal"
                ax_c.text(
                    c_i,
                    r_i,
                    txt,
                    ha="center",
                    va="center",
                    fontsize=7.0,
                    color=col_txt,
                    fontweight=weight,
                )
        ax_c.set_xticks(range(len(pop_cols)))
        ax_c.set_xticklabels([c.removeprefix("prob_") for c in pop_cols], fontsize=8.5)
        ax_c.set_yticks(range(len(ph_show)))
        ax_c.set_yticklabels(
            [
                f"{r.gene_symbol} CS{int(r.cs_id)} (Fst={r.max_wc_fst:.2f})"
                for r in ph_show.itertuples()
            ],
            fontsize=7.3,
        )
        ax_c.grid(False)
        cbar = fig.colorbar(im, ax=ax_c, fraction=0.038, pad=0.03)
        cbar.set_label("SuSiEx POST-HOC_PROB", fontsize=8.5)
    else:
        fs.no_data(ax_c, "No post-hoc probability table")
    fs.title_with_letter(
        ax_c, "c", "Ancestry-specific CS probabilities (bold >= 0.80)", dx=-0.34, dy=1.05
    )

    # d: Cross-ancestry sharing (number of active ancestries at prob >= 0.80) vs Fst & MAGE overlap
    ax_d = axes[1, 1]
    if not posthoc.empty:
        for overlap, grp in posthoc.groupby("overlaps_published_mage_cs"):
            lbl = "Overlaps published MAGE CS" if overlap else "Additional federated CS"
            col = fs.PATH["federated"] if overlap else fs.RED_RAMP[1]
            rng = np.random.default_rng(11 if overlap else 19)
            xj = grp["n_active_pops_prob80"].to_numpy(dtype=float) + rng.uniform(
                -0.15, 0.15, size=len(grp)
            )
            ax_d.scatter(
                xj,
                grp["max_wc_fst"].to_numpy(dtype=float),
                s=45,
                color=col,
                alpha=0.80,
                edgecolor=fs.SURFACE,
                lw=0.8,
                label=f"{lbl} (n={len(grp)})",
            )
        ax_d.set_xticks(range(0, len(POPULATION_ORDER) + 1))
        ax_d.set_xlabel("Active ancestries per CS (POST-HOC_PROB >= 0.80)")
        ax_d.set_ylabel("Locus max cross-ancestry Fst")
        ax_d.legend(loc="upper right", fontsize=8)
    else:
        fs.no_data(ax_d)
    fs.title_with_letter(ax_d, "d", "Ancestry sharing vs Fst differentiation", dx=-0.16, dy=1.05)

    ms_ok = fed_ms.dropna(subset=["published_lead_pip", "pub_lead_arm_pip"])
    ms_note = ""
    if len(ms_ok) >= 2:
        ms_r = float(np.corrcoef(ms_ok["published_lead_pip"], ms_ok["pub_lead_arm_pip"])[0, 1])
        ms_in_cs = ms_ok["pub_lead_in_cs"].mean() * 100
        ms_jacc = float(fed_ms["best_cs_jaccard"].mean())
        ms_note = (
            f" Fed MultiSuSiE ({len(ms_ok)}/{len(fed_ms)} loci): "
            f"in-CS = {ms_in_cs:.0f}%, mean CS Jaccard = {ms_jacc:.3f}, r = {ms_r:.3f}."
        )
    fs.footnote(
        fig,
        "Comparison against MAGE v1.0 published single-cohort SuSiE credible sets "
        "(eQTL_finemapping.*) and Finemapped_Fst.csv. POST-HOC_PROB >= 0.80 follows the "
        f"ancestry-active threshold; 'NA' denotes a monomorphic/MAF-filtered variant.{ms_note}",
        y=-0.03,
    )
    return _save_pair(fig, out_stem)


def fig5_regional_exemplars(
    federation_root: Path,
    arms_dir: Path,
    reports_dir: Path,
    out_stem: Path,
    exemplar_loci: Sequence[str] | None = None,
) -> list[Path]:
    """Fig 5: Multi-track regional LocusZoom-style profiles for representative MAGE exemplars."""
    fs.apply_style()
    federation_root = Path(federation_root)
    arms_dir = Path(arms_dir)
    reports_dir = Path(reports_dir)
    arm_loci = pd.read_csv(reports_dir / "arm_locus_metrics.tsv", sep="\t")
    panel = pd.read_csv(federation_root / "meta" / "panel_manifest.tsv", sep="\t").set_index(
        "locus_id"
    )
    fed_sx = arm_loci[arm_loci["arm_id"] == "federated_susiex"].set_index("locus_id")

    if exemplar_loci is None:
        chosen: list[str] = []
        for s in ("page_coloc", "high_fst", "multi_cs"):
            sub = fed_sx[(fed_sx["stratum"] == s) & (fed_sx["pub_cs_captured"] > 0)].sort_values(
                ["top_pip", "n_cs", "pub_cs_captured"], ascending=[False, False, False]
            )
            if sub.empty:
                sub = fed_sx[fed_sx["stratum"] == s].sort_values("top_pip", ascending=False)
            if not sub.empty:
                chosen.append(str(sub.index[0]))
        if not chosen:
            chosen = fed_sx.index.tolist()[:3]
    else:
        chosen = list(exemplar_loci)[:3]

    track_arms = {
        "fed_sx": arms_dir / "federated_susiex",
        "fed_ms": arms_dir / "federated_multisusie",
        "eur": arms_dir / "anc_EUR_susiex",
        "afr": arms_dir / "anc_AFR_susiex",
    }
    track_arms = {k: v for k, v in track_arms.items() if v.exists()}

    n_cols = len(chosen)
    fig, axes = plt.subplots(4, n_cols, figsize=(5.2 * n_cols, 10.6), sharex="col")
    if n_cols == 1:
        axes = np.expand_dims(axes, axis=1)
    plt.subplots_adjust(wspace=0.22, hspace=0.20)

    for col_i, lid in enumerate(chosen):
        locus = panel.loc[lid]
        prof = extract_locus_regional_profiles(federation_root, track_arms, lid)
        mb = prof["bp"].to_numpy(dtype=float) / 1e6
        tss_mb = float(locus["tss"]) / 1e6
        pub_lead = normalize_variant_id(str(locus["published_lead_kgp"]))
        lead_rows = prof[prof["snp_id"] == pub_lead]
        lead_mb = float(lead_rows["bp"].iloc[0]) / 1e6 if not lead_rows.empty else None

        # Row 0: Marginal -log10(P) across 5 ancestries
        ax0 = axes[0, col_i]
        p_col = "max_neglog10p_fed_sx" if "max_neglog10p_fed_sx" in prof else None
        if p_col:
            yv = prof[p_col].to_numpy(dtype=float)
            ok = np.isfinite(yv)
            ax0.scatter(mb[ok], yv[ok], s=10, color=fs.MUTED, alpha=0.6, lw=0)
            in_cs = (
                prof["cs_fed_sx"].notna().to_numpy()
                if "cs_fed_sx" in prof
                else np.zeros(len(prof), dtype=bool)
            )
            if in_cs.any():
                ax0.scatter(
                    mb[ok & in_cs],
                    yv[ok & in_cs],
                    s=30,
                    color=fs.PATH["federated"],
                    edgecolor=fs.SURFACE,
                    lw=0.6,
                    zorder=4,
                )
        ax0.axvline(tss_mb, ls=":", lw=1.0, color=fs.INK2, alpha=0.7)
        if lead_mb is not None:
            ax0.axvline(lead_mb, ls="--", lw=0.9, color=fs.STATUS["good"], alpha=0.7)
        stratum_lbl = STRATUM_LABELS.get(locus["stratum"], locus["stratum"])
        ax0.set_title(
            f"{locus['gene_symbol']} (chr{int(locus['chrom'])}, {stratum_lbl})",
            fontsize=10.5,
            fontweight="semibold",
            loc="left",
        )
        if col_i == 0:
            ax0.set_ylabel("Max ancestry\n$-\\log_{10}(P)$", fontsize=9)

        # Row 1: Published MAGE PIP
        ax1 = axes[1, col_i]
        mp = prof["mage_pip"].to_numpy(dtype=float)
        ok_m = np.isfinite(mp)
        ax1.scatter(mb[~ok_m], np.zeros((~ok_m).sum()), s=6, color=fs.GRID, alpha=0.4, lw=0)
        if ok_m.any():
            ax1.scatter(
                mb[ok_m], mp[ok_m], s=28, color=fs.INK2, edgecolor=fs.SURFACE, lw=0.6, zorder=4
            )
        ax1.set_ylim(-0.05, 1.08)
        if col_i == 0:
            ax1.set_ylabel("Published MAGE\nSuSiE PIP", fontsize=9)

        # Row 2: Federated SuSiEx PIP + CS highlights
        ax2 = axes[2, col_i]
        if "pip_fed_sx" in prof:
            yv = prof["pip_fed_sx"].to_numpy(dtype=float)
            in_cs = prof["cs_fed_sx"].notna().to_numpy()
            ax2.scatter(
                mb[~in_cs],
                np.nan_to_num(yv[~in_cs], nan=0.0),
                s=8,
                color=fs.MUTED,
                alpha=0.45,
                lw=0,
            )
            if in_cs.any():
                ax2.scatter(
                    mb[in_cs],
                    yv[in_cs],
                    s=34,
                    color=fs.PATH["federated"],
                    edgecolor=fs.SURFACE,
                    lw=0.7,
                    zorder=4,
                )
            sx_row = fed_sx.loc[lid] if lid in fed_sx.index else None
            if sx_row is not None:
                min_sz = int(sx_row["min_cs_size"]) if pd.notna(sx_row["min_cs_size"]) else "NA"
                ax2.text(
                    0.97,
                    0.90,
                    f"n_CS={int(sx_row['n_cs'])}, min={min_sz}",
                    transform=ax2.transAxes,
                    ha="right",
                    va="top",
                    fontsize=7.6,
                    color=fs.PATH["federated"],
                    bbox=dict(
                        boxstyle="round,pad=0.2",
                        facecolor=fs.SURFACE,
                        edgecolor=fs.GRID,
                        lw=0.7,
                    ),
                )
        ax2.set_ylim(-0.05, 1.08)
        if col_i == 0:
            ax2.set_ylabel("Federated SuSiEx\nPIP (5-ancestry)", fontsize=9)

        # Row 3: Federated MultiSuSiE vs Single-Ancestry (EUR, AFR)
        ax3 = axes[3, col_i]
        for arm_k, lbl, col, mk, sz in [
            ("pip_fed_ms", "Fed MultiSuSiE", fs.BLUE_RAMP[2], "o", 26),
            ("pip_afr", "AFR-only SuSiEx", MAGE_ANCESTRY_COLORS["AFR"], "^", 24),
            ("pip_eur", "EUR-only SuSiEx", MAGE_ANCESTRY_COLORS["EUR"], "s", 22),
        ]:
            if arm_k in prof:
                yv = prof[arm_k].to_numpy(dtype=float)
                hi = np.isfinite(yv) & (yv > 0.02)
                if hi.any():
                    ax3.scatter(
                        mb[hi],
                        yv[hi],
                        s=sz,
                        marker=mk,
                        color=col,
                        alpha=0.85,
                        edgecolor=fs.SURFACE,
                        lw=0.5,
                        label=lbl,
                        zorder=4,
                    )
        ax3.set_ylim(-0.05, 1.08)
        ax3.set_xlabel(f"Chromosome {int(locus['chrom'])} position (Mb, GRCh38)")
        if col_i == 0:
            ax3.set_ylabel("MultiSuSiE &\nsingle-anc PIP", fontsize=9)
            ax3.legend(
                loc="upper left",
                fontsize=7.5,
                frameon=True,
                facecolor=fs.SURFACE,
                edgecolor=fs.GRID,
            )

    fs.footnote(
        fig,
        "Regional profiles across representative MAGE loci (+/-250 kb around dotted TSS; "
        "dashed green marks published MAGE lead). Row 1 shows peak marginal -log10(P); "
        "row 2 shows published MAGE CS PIPs; rows 3-4 compare Federated SuSiEx and MultiSuSiE.",
        y=-0.02,
    )
    return _save_pair(fig, out_stem)


def render_all_mage_figures(
    federation_root: Path,
    out_root: Path,
) -> list[Path]:
    """Render all 5 real-data MAGE validation figures into ``out_root / 'figures'``."""
    federation_root = Path(federation_root)
    out_root = Path(out_root)
    arms_dir = out_root / "arms"
    reports_dir = out_root / "reports"
    fig_dir = out_root / "figures"
    written: list[Path] = []
    written.extend(
        fig1_cohort_and_design(federation_root, reports_dir, fig_dir / "fig1_cohort_and_design")
    )
    written.extend(
        fig2_federated_parity_and_controls(
            arms_dir, reports_dir, fig_dir / "fig2_federated_parity_and_controls"
        )
    )
    written.extend(fig3_resolution_and_arms(reports_dir, fig_dir / "fig3_resolution_and_arms"))
    written.extend(
        fig4_mage_agreement_and_ancestry(reports_dir, fig_dir / "fig4_mage_agreement_and_ancestry")
    )
    written.extend(
        fig5_regional_exemplars(
            federation_root, arms_dir, reports_dir, fig_dir / "fig5_regional_locus_exemplars"
        )
    )
    return written
