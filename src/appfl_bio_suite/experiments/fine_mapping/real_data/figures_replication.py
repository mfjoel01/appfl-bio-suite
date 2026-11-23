"""Figures for the MAGE replication design (``benchmark.run_replication_suite``).

1. ``figR1_replication`` -- the pooled analysis against the published credible sets.
2. ``figR2_sites_and_federation`` -- each site alone against pooled; federated against pooled.
3. ``figR3_regional`` -- PIP tracks at two high-Fst loci, where covariate handling matters.

Colours come from the shared, validated figure system (``figures/figstyle.py``): the two
``PATH`` hues for pooled and federated, the three ``SITE`` hues for the three sites, and the
muted ink for the published reference. Each panel is single-series or directly labelled.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import (
    read_bed_variants,
    read_bim,
    read_fam,
)
from appfl_bio_suite.experiments.fine_mapping.figures import figstyle as fs
from appfl_bio_suite.experiments.fine_mapping.figures.figstyle import plt

from .benchmark import _instance_fit, normalize_variant_id
from .federation_builder import SITE_ORDER
from .mage import residualize_design

POOLED = fs.PATH["centralized"]
FEDERATED = fs.PATH["federated"]
REFERENCE = fs.MUTED
SITE_COLOR = dict(zip(SITE_ORDER, fs.brand.CATEGORICAL[:3], strict=True))
SITE_NAME = {
    "site_atlantic": "Atlantic",
    "site_africa_eurasia": "Africa-Eurasia",
    "site_asiapacific": "Asia-Pacific",
}

# (arm, label) in display order; arms missing from a run are skipped.
REPLICATION_BARS = (
    ("cen_susie_fwl_731", "SuSiE, published setup (N = 731)"),
    ("cen_susie_fwl_731_1mb", "  same, ±1 Mb window"),
    ("cen_susie_sitecov_731", "  same, site-derived covariates"),
    ("cen_susie_fwl_730", "SuSiE, swapped library removed (730)"),
    ("cen_susie_yresid_730", "SuSiE, expression-only adjustment (730)"),
    ("cen_multisusie_fwl_730", "MultiSuSiE, 5 ancestries (730)"),
    ("cen_susiex1pop_fwl_731", "SuSiEx, one population (731)"),
    ("cen_susiex_fwl_730", "SuSiEx, 5 ancestries (730)"),
)


def _save(fig, stem: Path) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = [stem.with_suffix(".png"), stem.with_suffix(".pdf")]
    fig.savefig(paths[0], dpi=300, bbox_inches="tight", facecolor=fs.SURFACE)
    fig.savefig(paths[1], bbox_inches="tight", facecolor=fs.SURFACE)
    plt.close(fig)
    return paths


def _published(federation_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = Path(federation_root) / "meta"
    panel = pd.read_csv(meta / "panel_manifest.tsv", sep="\t")
    pub = pd.read_csv(meta / "published_cs_variants.tsv", sep="\t")
    pub["snp"] = pub["variant_kgpID"].map(normalize_variant_id)
    return panel, pub


def _published_pip_pairs(arm_dir: Path, panel: pd.DataFrame, pub: pd.DataFrame):
    ours, theirs = [], []
    for loc in panel.itertuples(index=False):
        pips, _ = _instance_fit(arm_dir, str(loc.locus_id))
        pv = pub[pub["ensemblID"] == loc.ensembl_id]
        ours += [pips.get(s, 0.0) for s in pv["snp"]]
        theirs += pv["variantPIP"].astype(float).tolist()
    return np.asarray(theirs), np.asarray(ours)


def lead_pvalues(federation_root: Path) -> pd.DataFrame:
    """Pooled p-value of each published lead, adjusting expression only vs both (FWL)."""
    meta = Path(federation_root) / "meta"
    pre = meta / "pooled_plink" / "pooled_genotypes"
    bim, fam = read_bim(f"{pre}.bim"), read_fam(f"{pre}.fam")
    panel, _ = _published(federation_root)
    npz = np.load(meta / "expression_panel.npz", allow_pickle=True)
    samples = [str(s) for s in npz["samples"]]
    perm = [samples.index(i) for i in fam["IID"]]
    cov = pd.read_csv(meta / "covariates.tsv", sep="\t", dtype={"sample_id": str})
    Z = cov.set_index("sample_id").loc[fam["IID"].tolist()].to_numpy(dtype=np.float64)
    locus_ids = [str(x) for x in npz["locus_ids"]]

    def pval(x, y, extra_df):
        x, y = x - x.mean(), y - y.mean()
        b = x @ y / (x @ x)
        df = len(y) - 2 - extra_df
        se = np.sqrt(((y - b * x) @ (y - b * x)) / df / (x @ x))
        return 2 * stats.t.sf(abs(b / se), df)

    rows = []
    for loc in panel.itertuples(index=False):
        lead = normalize_variant_id(loc.published_lead_kgp)
        hit = np.flatnonzero(bim["snp_id"].astype(str).to_numpy() == lead)
        if hit.size != 1:
            continue
        x = read_bed_variants(pre, hit, len(fam))[:, 0].astype(np.float64)
        k = locus_ids.index(str(loc.locus_id))
        y_resid = npz["expr_resid"][k, perm].astype(np.float64)
        x_perp, y_perp = residualize_design(x[:, None], npz["expr_raw"][k, perm], Z)
        rows.append(
            {
                "locus_id": loc.locus_id,
                "gene": str(loc.gene_symbol),
                "stratum": str(loc.stratum),
                "p_expression_only": pval(x, y_resid, 0),
                "p_fwl": pval(x_perp[:, 0], y_perp, Z.shape[1]),
            }
        )
    return pd.DataFrame(rows)


def figR1_replication(federation_root: Path, out_root: Path, stem: Path) -> list[Path]:
    """Pooled SuSiE reproduces the published sets; adjusting genotypes too is what does it."""
    fs.apply_style()
    arms = Path(out_root) / "arms"
    reports = Path(out_root) / "reports"
    panel, pub = _published(federation_root)
    vs_mage = pd.read_csv(reports / "replication_vs_mage.tsv", sep="\t").set_index("arm_id")

    fig = plt.figure(figsize=(12.0, 9.6))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[1.0, 1.25], height_ratios=[1.0, 0.82], wspace=0.42, hspace=0.42
    )

    # a. published vs our PIPs over published credible-set variants
    ax = fig.add_subplot(gs[0, 0])
    main = "cen_susie_fwl_731" if "cen_susie_fwl_731" in vs_mage.index else "cen_susie_fwl_730"
    x, y = _published_pip_pairs(arms / main, panel, pub)
    ax.plot([0, 1], [0, 1], color=fs.GRID, lw=1, zorder=1)
    ax.scatter(x, y, s=16, color=POOLED, edgecolor=fs.SURFACE, linewidth=0.6, zorder=3)
    r = float(np.corrcoef(x, y)[0, 1])
    row = vs_mage.loc[main]
    ax.text(
        0.97,
        0.44,
        f"r = {r:.2f}\n{int(row['matched'])}/{int(row['published_cs'])} published sets matched"
        f"\n{len(x)} credible-set variants",
        transform=ax.transAxes,
        va="top",
        ha="right",
        fontsize=8.5,
        color=fs.INK,
    )
    ax.set(
        xlim=(-0.03, 1.03),
        ylim=(-0.03, 1.03),
        xlabel="Published MAGE PIP",
        ylabel="Pooled SuSiE PIP (N = 731)",
    )
    fs.panel_letter(ax, "a", "Replication of published PIPs", dx=-0.16)

    # b. lead-variant evidence, expression-only vs genotype+expression adjustment
    ax = fig.add_subplot(gs[0, 1])
    lp = lead_pvalues(federation_root).sort_values("p_fwl", ascending=False)
    ypos = np.arange(len(lp))
    a = -np.log10(lp["p_expression_only"].to_numpy())
    b = -np.log10(lp["p_fwl"].to_numpy())
    ax.hlines(ypos, np.minimum(a, b), np.maximum(a, b), color=fs.GRID, lw=1.5, zorder=1)
    ax.scatter(a, ypos, s=30, color=REFERENCE, zorder=3, label="Expression only")
    ax.scatter(b, ypos, s=30, color=POOLED, zorder=3, label="Genotypes and expression (FWL)")
    labels = [
        f"{g}{'  (high Fst)' if s == 'high_fst' else ''}"
        for g, s in zip(lp["gene"], lp["stratum"], strict=True)
    ]
    ax.set_yticks(ypos, labels, fontsize=8)
    ax.set_xscale("log")
    ticks = [0.1, 1, 3, 10, 30, 100]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel("−log10 p of the published lead (N = 730, log scale)")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    fs.panel_letter(ax, "b", "Covariate adjustment and the lead variant", dx=-0.30)

    # c. published sets matched, per pooled analysis (value in the label: static figure)
    ax = fig.add_subplot(gs[1, :])
    bars = [(a, lab) for a, lab in REPLICATION_BARS if a in vs_mage.index]
    vals = [int(vs_mage.loc[a, "matched"]) for a, _ in bars]
    total = int(vs_mage.iloc[0]["published_cs"])
    ypos = np.arange(len(bars))[::-1]
    colors = [POOLED if a == main else REFERENCE for a, _ in bars]
    ax.barh(ypos, vals, height=0.62, color=colors, edgecolor=fs.SURFACE, linewidth=1.5)
    ax.set_yticks(
        ypos, [f"{lab}   {v}/{total}" for (_, lab), v in zip(bars, vals, strict=True)], fontsize=8.5
    )
    ax.set_xlim(0, total)
    ax.set_xlabel(f"Published credible sets matched (Jaccard ≥ 0.5), of {total}")
    fs.panel_letter(ax, "c", "Which pooled analysis reproduces MAGE", dx=-0.28)
    fs.footnote(
        fig,
        "Published MAGE: susieR on 731 samples, covariates (sex, 5 genotype PCs, 60 PEER) "
        "regressed from genotypes and expression, L = 10. Here: ±250 kb windows unless noted; "
        "SuSiE run as MultiSuSiE with one population. 46 published sets across 16 eGenes.",
        y=0.02,
    )
    return _save(fig, stem)


def figR2_sites_and_federation(federation_root: Path, out_root: Path, stem: Path) -> list[Path]:
    """Each site alone recovers a fraction of the pooled result; the federation all of it."""
    fs.apply_style()
    arms = Path(out_root) / "arms"
    sites = pd.read_csv(Path(out_root) / "reports" / "sites_vs_pooled.tsv", sep="\t")
    fig, grid = plt.subplots(
        2,
        2,
        figsize=(12.0, 9.4),
        gridspec_kw={"height_ratios": [0.8, 1.0], "wspace": 0.55, "hspace": 0.45},
    )
    axes = grid.ravel()

    for ax, method, title, letter in (
        (axes[0], "susie", "SuSiE", "a"),
        (axes[1], "multisusie", "MultiSuSiE", "b"),
    ):
        rows = sites[sites["method"] == method].set_index("arm_id")
        order = [f"{s}_{method}_fwl" for s in SITE_ORDER] + [f"fed_{method}_fwl_730"]
        order = [o for o in order if o in rows.index]
        ref = int(rows.iloc[0]["reference_cs"])
        vals = [int(rows.loc[o, "recovered"]) for o in order]
        names = [SITE_NAME.get(o.rsplit(f"_{method}", 1)[0], "Federated, 3 sites") for o in order]
        colors = [SITE_COLOR.get(o.rsplit(f"_{method}", 1)[0], FEDERATED) for o in order]
        ypos = np.arange(len(order))[::-1]
        ax.barh(ypos, vals, height=0.62, color=colors, edgecolor=fs.SURFACE, linewidth=1.5)
        ax.set_yticks(
            ypos, [f"{n}   {v}/{ref}" for n, v in zip(names, vals, strict=True)], fontsize=8.5
        )
        ax.set_xlim(0, ref * 1.02)
        ax.set_xlabel(f"Pooled credible sets recovered (of {ref})")
        fs.panel_letter(ax, letter, f"{title}: each site alone vs federated", dx=-0.55)

    for ax, method, title, letter in (
        (axes[2], "susie", "SuSiE", "c"),
        (axes[3], "multisusie", "MultiSuSiE", "d"),
    ):
        cen, fed = arms / f"cen_{method}_fwl_730", arms / f"fed_{method}_fwl_730"
        xs, ys = [], []
        for inst in sorted((cen / "data" / "instances").glob("*")):
            lid = inst.name.split("_")[0]
            pc, _ = _instance_fit(cen, lid)
            pf, _ = _instance_fit(fed, lid)
            common = [k for k in pc if k in pf]
            xs += [pc[k] for k in common]
            ys += [pf[k] for k in common]
        xs, ys = np.asarray(xs), np.asarray(ys)
        ax.plot([0, 1], [0, 1], color=fs.GRID, lw=1, zorder=1)
        big = xs > 0.01
        ax.scatter(
            xs[big], ys[big], s=14, color=FEDERATED, edgecolor=fs.SURFACE, linewidth=0.5, zorder=3
        )
        ax.text(
            0.03,
            0.97,
            f"max |ΔPIP| = {np.max(np.abs(xs - ys)):.1e}\n{len(xs):,} variants, 16 loci",
            transform=ax.transAxes,
            va="top",
            fontsize=8.5,
            color=fs.INK,
        )
        ax.set(xlim=(-0.03, 1.03), ylim=(-0.03, 1.03), xlabel="Pooled PIP", ylabel="Federated PIP")
        fs.panel_letter(ax, letter, f"{title}: federated = pooled", dx=-0.16)
    fs.footnote(
        fig,
        "Clean cohort (N = 730), covariates regressed from genotypes and expression (FWL); "
        "in the federation the adjustment is computed from the sites' covariate "
        "cross-moments. A pooled set is recovered at Jaccard ≥ 0.5. Panels c-d plot variants "
        "with pooled PIP > 0.01; the maximum is over all variants. MultiSuSiE fits five "
        "ancestry columns.",
        y=0.02,
    )
    return _save(fig, stem)


def figR3_regional(
    federation_root: Path, out_root: Path, stem: Path, loci: tuple[str, ...] = ("ERN1", "NPNT")
) -> list[Path]:
    """PIP tracks at high-Fst loci: published, pooled FWL, expression-only, one site."""
    fs.apply_style()
    arms = Path(out_root) / "arms"
    panel, pub = _published(federation_root)
    tracks = (
        ("Published MAGE", None, REFERENCE),
        ("Pooled = federated, FWL", "cen_susie_fwl_730", POOLED),
        ("Pooled, expression-only adjustment", "cen_susie_yresid_730", REFERENCE),
        ("Asia-Pacific site alone", "site_asiapacific_susie_fwl", SITE_COLOR["site_asiapacific"]),
    )
    chosen = [panel[panel["gene_symbol"] == g].iloc[0] for g in loci if g in set(panel.gene_symbol)]
    if not chosen:  # another panel: fall back to its high-Fst loci, then to any
        pool = panel[panel["stratum"] == "high_fst"]
        chosen = [r for _, r in (pool if len(pool) else panel).head(2).iterrows()]
    fig, axes = plt.subplots(
        len(tracks), len(chosen), figsize=(6.2 * len(chosen), 7.2), sharex="col", squeeze=False
    )
    for j, loc in enumerate(chosen):
        pv = pub[pub["ensemblID"] == loc.ensembl_id]
        pub_pip = dict(zip(pv["snp"], pv["variantPIP"].astype(float), strict=True))
        pub_cs = set(pv["snp"])
        lead = normalize_variant_id(loc.published_lead_kgp)
        pos_of = {}
        base, _ = _instance_fit(arms / "cen_susie_fwl_730", str(loc.locus_id))
        for snp in base:
            pos_of[snp] = int(snp.split(":")[1])
        for i, (label, arm, color) in enumerate(tracks):
            ax = axes[i, j]
            if arm is None:
                pips, sets = pub_pip, []
            else:
                pips, sets = _instance_fit(arms / arm, str(loc.locus_id))
            in_cs = set().union(*sets) if sets else set()
            snps = [s for s in pips if s in pos_of]
            x = np.array([pos_of[s] for s in snps]) / 1e6
            y = np.array([pips[s] for s in snps])
            cs_mask = np.array([s in (pub_cs if arm is None else in_cs) for s in snps])
            ax.scatter(x[~cs_mask], y[~cs_mask], s=6, color=fs.GRID, zorder=2)
            ax.scatter(x[cs_mask], y[cs_mask], s=18, color=color, zorder=3)
            if lead in pos_of:
                ax.axvline(pos_of[lead] / 1e6, color=fs.INK, lw=0.6, alpha=0.35, zorder=1)
            ax.set_ylim(-0.04, 1.08)
            ax.set_yticks([0, 0.5, 1])
            ax.text(0.01, 0.96, label, transform=ax.transAxes, va="top", fontsize=8.5, color=fs.INK)
            if j == 0:
                ax.set_ylabel("PIP")
            if i == 0:
                fs.panel_letter(
                    ax, "ab"[j], f"{loc.gene_symbol} (chr{loc.chrom}, max Fst {loc.max_wc_fst:.2f})"
                )
        axes[-1, j].set_xlabel(f"chr{loc.chrom} position (Mb)")
    fs.footnote(
        fig,
        "Darker points are credible-set members, grey points other variants; MAGE releases "
        "PIPs only for its credible-set variants, so the top row has no background. The "
        "vertical line marks the published lead. Pooled and federated PIPs agree to 1e-13, "
        "so one track shows both. Clean cohort (N = 730) except the published row (731).",
        y=-0.03,
    )
    return _save(fig, stem)


def render_replication_figures(federation_root: Path, out_root: Path) -> list[Path]:
    """Render all replication figures into ``<out_root>/figures``."""
    figs = Path(out_root) / "figures"
    paths = []
    paths += figR1_replication(federation_root, out_root, figs / "figR1_replication")
    paths += figR2_sites_and_federation(
        federation_root, out_root, figs / "figR2_sites_and_federation"
    )
    paths += figR3_regional(federation_root, out_root, figs / "figR3_regional")
    return paths
