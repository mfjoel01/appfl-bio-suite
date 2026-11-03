"""Render accepted matched-estimator results without refitting or pooling cohorts.

Requires the coordinator plotting dependencies. Input paths in historical acceptance
records may predate relocation; resolve final tables by basename in --root and verify
against their original checksums before plotting.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

METRICS = {
    "power": ("Any-causal capture", "%", 100),
    "causal_recall": ("Causal-variant recall", "%", 100),
    "cs_coverage": ("Credible-set coverage", "%", 100),
    "mean_cs_size": ("Mean credible-set size", "variants", 1),
    "fdr_pip95": ("FDR at PIP > 0.95", "%", 100),
    "high_confidence_yield": ("Any true causal with PIP > 0.95", "%", 100),
}
POLICIES = {"common_purity_0.5": "Common purity ≥ 0.5", "unfiltered": "Unfiltered sets"}
METHODS = {"multisusie": ("MultiSuSiE", "#75468c"), "susiex": ("SuSiEx", "#087f8c")}
PHASES = {"benchmark": "Benchmark · 79 loci", "evaluation": "Held-out evaluation · 54 loci"}
ARMS = {
    "federation": "Full federation (150k)",
    "anl": "ANL (50k)",
    "covenant": "Covenant (50k)",
    "mbzuai": "MBZUAI (50k)",
    "federation_50k": "Matched 50k · draw 1",
    "federation_50k_seed2": "Matched 50k · draw 2",
    "federation_50k_seed3": "Matched 50k · draw 3",
    **{f"csa_to_anl_afr_1000_seed202609{d}": f"Selected 50k · draw {d - 27}" for d in (28, 29, 30)},
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.2)
    ax.set_axisbelow(True)


def interval(ax, row, y, scale, color, label=None):
    # Percentile intervals need not straddle the point estimate: draw endpoints directly.
    x, low, high = (float(row[k]) * scale for k in ("estimate", "ci_low", "ci_high"))
    if np.isfinite([low, high]).all():
        ax.plot([low, high], [y, y], color=color, linewidth=1.8)
    if np.isfinite(x):
        ax.plot(x, y, "o", color=color, markersize=5, label=label)


def save(fig, out, name, title, caption, figures):
    fig.suptitle(
        title + "\nPoints and 95% locus-bootstrap intervals · 2,000 replicates", fontsize=13
    )
    for extension in ("png", "pdf", "svg"):
        fig.savefig(out / "graphs" / f"{name}.{extension}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    figures.append(dict(name=name, title=title, caption=caption))


def render(root: Path, out: Path):
    accepted = json.loads((root / "accepted.json").read_text())
    if not accepted["pass_"]:
        raise ValueError("Comparison has not passed acceptance")
    for name, digest in accepted["sha256"].items():
        path = root / Path(name).name
        if sha(path) != digest:
            raise ValueError(f"Accepted numerical output changed: {path}")
    out.mkdir(parents=True, exist_ok=False)
    (out / "graphs").mkdir()
    (out / "tables").mkdir()
    tables = {}
    for name in (
        "summary",
        "paired_differences",
        "summary_stratified",
        "paired_differences_stratified",
    ):
        tables[name] = pd.read_csv(root / f"{name}.tsv", sep="\t")
        shutil.copy2(root / f"{name}.tsv", out / "tables" / f"{name}.tsv")
    shutil.copy2(root / "RESULTS.md", out / "RESULTS.md")
    for name in ("accepted.json", "design.json"):
        shutil.copy2(root / name, out / name)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    figures = []
    for phase, phase_label in PHASES.items():
        frame = tables["summary"]
        frame = frame[(frame.phase == phase) & (frame.arm == "federation")]
        fig, axes = plt.subplots(3, 2, figsize=(12, 9), layout="constrained")
        for ax, (metric, (label, unit, scale)) in zip(axes.flat, METRICS.items(), strict=True):
            for index, policy in enumerate(POLICIES):
                for method, (method_label, color) in METHODS.items():
                    row = frame[
                        (frame.metric == metric)
                        & (frame.policy == policy)
                        & (frame.method == method)
                    ].squeeze()
                    offset = -0.13 if method == "multisusie" else 0.13
                    interval(
                        ax, row, index + offset, scale, color, method_label if index == 0 else None
                    )
            ax.set(
                yticks=[0, 1],
                yticklabels=list(POLICIES.values()),
                ylim=(1.5, -0.5),
                title=label,
                xlabel=unit,
            )
            if metric == "cs_coverage":
                ax.axvline(95, color="#777777", linestyle="--", label="Nominal 95%")
            if metric == "power":
                ax.legend(loc="best", fontsize=8)
            style(ax)
        title = f"MultiSuSiE vs SuSiEx | {phase_label} | Full 150k cohort"
        fig.suptitle(title, fontsize=15)
        save(
            fig,
            out,
            f"{phase}_full_cohort",
            title,
            "Points and 95% locus-bootstrap intervals. Each policy uses the same fitted PIPs. "
            "Coverage targets 95%; mean size includes all retained sets. "
            "Axis scales differ by metric; these are estimator levels, not paired differences.",
            figures,
        )
        for policy, policy_label in POLICIES.items():
            frame = tables["paired_differences"]
            frame = frame[(frame.phase == phase) & (frame.policy == policy)]
            arms = sorted(frame.arm.unique(), key=lambda a: (a != "federation", a))
            fig, axes = plt.subplots(2, 3, figsize=(16, 8), layout="constrained")
            for ax, (metric, (label, unit, scale)) in zip(axes.flat, METRICS.items(), strict=True):
                for y, arm in enumerate(arms):
                    row = frame[(frame.arm == arm) & (frame.metric == metric)].squeeze()
                    interval(ax, row.rename({"difference": "estimate"}), y, scale, "#75468c")
                ax.axvline(0, color="#666666", linestyle="--", linewidth=1)
                ax.set(
                    yticks=range(len(arms)),
                    yticklabels=[ARMS[a] for a in arms],
                    ylim=(len(arms) - 0.5, -0.5),
                    title=label,
                    xlabel="Difference (variants)"
                    if unit == "variants"
                    else "Difference (percentage points)",
                )
                style(ax)
            title = f"{phase_label} | {policy_label} | MultiSuSiE − SuSiEx"
            fig.suptitle(title, fontsize=15)
            save(
                fig,
                out,
                f"{phase}_{policy}_cohorts",
                title,
                "Paired differences with marginal 95% locus-bootstrap intervals. Positive means "
                "a larger metric, not necessarily better. Cohort draws overlap and are not "
                "independent replications. Native PIP metrics do not change with the set filter.",
                figures,
            )
            frame = tables["paired_differences_stratified"]
            frame = frame[
                (frame.phase == phase) & (frame.arm == "federation") & (frame.policy == policy)
            ]
            keys = sorted(set(zip(frame.architecture_id, frame.causal_mode, strict=True)))
            fig, axes = plt.subplots(
                1, 3, figsize=(16, max(8, len(keys) * 0.38)), layout="constrained"
            )
            for ax, metric in zip(axes, ("causal_recall", "cs_coverage", "fdr_pip95"), strict=True):
                for y, (arch, mode) in enumerate(keys):
                    row = frame[
                        (frame.architecture_id == arch)
                        & (frame.causal_mode == mode)
                        & (frame.metric == metric)
                    ].squeeze()
                    interval(ax, row.rename({"difference": "estimate"}), y, 100, "#75468c")
                ax.axvline(0, color="#666666", linestyle="--", linewidth=1)
                ax.set(
                    yticks=range(len(keys)),
                    yticklabels=[f"{a} · {m}" for a, m in keys],
                    ylim=(len(keys) - 0.5, -0.5),
                    title=METRICS[metric][0],
                    xlabel="Difference (percentage points)",
                )
                if ax is not axes[0]:
                    ax.set_yticklabels([])
                style(ax)
            title = (
                f"{phase_label} | Full 150k | {policy_label} | MultiSuSiE − SuSiEx by architecture"
            )
            fig.suptitle(title, fontsize=15)
            save(
                fig,
                out,
                f"{phase}_{policy}_architectures",
                title,
                "MultiSuSiE − SuSiEx, stratified by architecture and causal mode. "
                "Intervals are marginal, with no multiplicity adjustment. Missing estimates or "
                "intervals are omitted, not plotted as zero. ncsl: causal count; h2: heritability; "
                "rg: simulated effect correlation (not MultiSuSiE's fixed rho=0.75 prior).",
                figures,
            )
    level = tables["summary"]
    level = level[
        (level.phase == "evaluation")
        & (level.arm == "federation")
        & (level.policy == "common_purity_0.5")
    ]
    display = level.pivot(index="metric", columns="method", values="estimate").reindex(METRICS)
    for metric, (_label, _unit, scale) in METRICS.items():
        display.loc[metric] *= scale
    display.index = [f"{METRICS[m][0]} ({METRICS[m][1]})" for m in display.index]
    display.columns = [METHODS[m][0] for m in display.columns]
    introduction = (
        "All 139,650 paired instances across 931 cohort/locus tasks passed acceptance. "
        "With common purity filtering in the full-cohort held-out evaluation, MultiSuSiE has "
        "slightly higher causal recall; SuSiEx has coverage closer to 95%, lower FDR and smaller "
        "credible sets. No universal winner emerges across policies and cohorts."
    )
    page = [
        "<!doctype html><html lang='en'><meta charset='utf-8'>",
        "<title>SuSiEx vs MultiSuSiE — completed matched comparison</title>",
        "<style>body{font:17px/1.5 system-ui;margin:40px auto;max-width:1200px;padding:0 24px;"
        "color:#243342}h1,h2{line-height:1.2}img{width:100%;height:auto}table{border-collapse:collapse}"
        "td,th{padding:8px 14px;border-bottom:1px solid #ddd;text-align:right}a{color:#087f8c}"
        "section{margin:48px 0}p{max-width:1000px}</style>",
        "<h1>SuSiEx vs MultiSuSiE</h1><p>Completed matched-filter comparison</p>",
        f"<p>{introduction}</p><h2>Full-cohort held-out results · common purity ≥ 0.5</h2>",
        display.to_html(float_format=lambda x: f"{x:.2f}"),
        "<p>95% intervals use 2,000 paired locus-bootstrap replicates. Native priors, PIP "
        "definitions and precision remain method-specific. The selected cohort was fixed "
        "using MultiSuSiE development results. Archived runtimes cannot support speed claims.</p>",
        "<p><a href='RESULTS.md'>Detailed report and intervals</a> · "
        "<a href='tables/summary.tsv'>All estimates</a> · "
        "<a href='tables/paired_differences.tsv'>All paired differences</a> · "
        "<a href='tables/summary_stratified.tsv'>Stratified estimates</a> · "
        "<a href='tables/paired_differences_stratified.tsv'>Stratified differences</a></p>",
    ]
    for figure in figures:
        name, title, caption = (figure[k] for k in ("name", "title", "caption"))
        page += [
            f"<section><h2>{html.escape(title)}</h2><p>{html.escape(caption)}</p>",
            f"<a href='graphs/{name}.pdf'>PDF</a> · <a href='graphs/{name}.svg'>SVG</a> · "
            f"<a href='graphs/{name}.png'>PNG</a>",
            f"<img src='graphs/{name}.png' alt='{html.escape(title, quote=True)}'></section>",
        ]
    (out / "index.html").write_text("\n".join(page) + "\n</html>\n")
    (out / "README.md").write_text(
        "# SuSiEx versus MultiSuSiE\n\nOpen index.html for the complete comparison gallery.\n\n"
        + introduction
        + "\n\nTen figures are supplied as PNG, PDF and SVG. Four tables contain all "
        "cohorts, policies and strata. RESULTS.md retains the detailed numerical "
        "interpretation.\n\n"
        "All numerical inputs were verified against the original acceptance checksums before "
        "rendering. report-manifest.json records input/output hashes and figure captions.\n"
    )
    manifest = dict(
        source=str(root.resolve()),
        n_paired_instances=accepted["n_paired_instances"],
        n_rows=accepted["n_rows"],
        generator_sha256=sha(__file__),
        source_sha256=accepted["sha256"],
        figures=figures,
        output_sha256={
            str(p.relative_to(out)): sha(p) for p in sorted(out.rglob("*")) if p.is_file()
        },
    )
    (out / "report-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {len(figures)} figures (PNG/PDF/SVG), tables and report to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    render(args.root, args.out)
