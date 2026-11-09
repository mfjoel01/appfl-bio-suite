#!/usr/bin/env python3
"""Build the fine-mapping results shown in the project site's Results tab.

    python scripts/fine-mapping/site_results.py \\
        --parity-run local/output/fine-mapping-full-rerun-20260914 \\
        --evaluation-run local/output/fine-mapping-final-20260928

The published results come from two runs:

* the PARITY run fine-mapped 11,850 simulated instances on both the federated and the
  centralized path over one package. It supplies the exploratory figures, the parity
  check, accuracy across the simulation design, and the cost of federating.
* the EVALUATION run searched 20 candidate 50,000-person cohorts on development draws,
  froze the winner, then re-ran every arm on fresh held-out draws. It supplies the
  comparison of the full federation against each site alone.

Figures are redrawn from the runs by the current plotting code, into a temporary
directory; the run directories are only read. The viewer embeds every result in the
map's metadata, so each figure's bytes are part of every load of the map. They are
scaled to the width the viewer displays and stored as lossless WebP, about half the
size of the equivalent PNG. Tables are reduced to readable columns.

Everything is written to website/results/fine-mapping/, which this script owns and
replaces. website/partners.json selects and captions the files, and
scripts/build_site.sh embeds them through `watch export`.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

import pandas as pd  # noqa: E402

from appfl_bio_suite.experiments.fine_mapping import composition  # noqa: E402
from appfl_bio_suite.experiments.fine_mapping.figures import (  # noqa: E402
    eda,
    federation,
    paper_plots,
)
from appfl_bio_suite.experiments.fine_mapping.reporting import clustered_ratio  # noqa: E402

DEFAULT_OUT = REPO / "website/results/fine-mapping"

# The viewer shows a figure at most 1,200 CSS pixels wide (its expanded dialog), so
# 1,800 covers a 1.5x display. The height cap keeps the tall locus zoom in proportion.
MAX_WIDTH, MAX_HEIGHT = 1800, 2000

# (run, path under the rendered figures, published name). The locus zoom is the
# parity example: one locus at the weakest heritability, resolved to its causal variant.
FIGURES = [
    ("parity", "eda/eda1_site_composition.png", "eda1_site_composition"),
    ("parity", "eda/eda4_ld_divergence_maps.png", "eda4_ld_divergence_maps"),
    ("parity", "eda/eda5_causal_maf.png", "eda5_causal_maf"),
    ("parity", "eda/eda6_h2_calibration.png", "eda6_h2_calibration"),
    ("parity", "eda/eda7_rg_recovery.png", "eda7_rg_recovery"),
    ("parity", "eda/eda8_design_grid.png", "eda8_design_grid"),
    ("parity", "federation/shared/fed1_parity.png", "fed1_parity"),
    ("parity", "federation/shared/fed2_ancestry_columns.png", "fed2_ancestry_columns"),
    ("parity", "federation/shared/fed3_uplink_cost.png", "fed3_uplink_cost"),
    ("parity", "federation/shared/fed5_where_time_goes.png", "fed5_where_time_goes"),
    ("parity", "federated/shared/pap1_recall_grid.png", "pap1_recall_grid"),
    ("parity", "federated/shared/pap2_diversity_panels.png", "pap2_diversity_panels"),
    ("parity", "federated/shared/pap3_convergence.png", "pap3_convergence"),
    ("parity", "federated/shared/pap4_power_coverage.png", "pap4_power_coverage"),
    ("parity", "federated/shared/pap5_resolution.png", "pap5_resolution"),
    ("parity", "federated/shared/pap6_high_confidence.png", "pap6_high_confidence"),
    ("parity", "federated/shared/pap7_ld_divergence.png", "pap7_ld_divergence"),
    ("parity", "federated/shared/pap8_scalability.png", "pap8_scalability"),
    ("parity", "federated/shared/pap10_population_probability.png", "pap10_population_probability"),
    ("parity", "federated/shared/pap12_effect_concordance.png", "pap12_effect_concordance"),
    (
        "parity",
        "federated/shared/pap11_locuszoom_L0000_ncsl1_h2-0.0005_rg1_rep0.png",
        "pap11_locuszoom_L0000",
    ),
    ("evaluation", "federation/shared/fed6_what_federation_buys.png", "fed6_what_federation_buys"),
    ("evaluation", "federation/harmonization.png", "fed4_harmonization"),
]

SITE_LABELS = {"anl": "ANL", "covenant": "Covenant", "mbzuai": "MBZUAI"}
# The five headline arms, in display order.
ARM_LABELS = {
    "federation": "Full federation",
    composition.SMART_ARM: "Smart composition",
    "covenant": "Covenant alone",
    "mbzuai": "MBZUAI alone",
    "anl": "ANL alone",
}

log = logging.getLogger("site_results")


# --------------------------------------------------------------------------- figures


def render_parity(root: Path, out: Path) -> None:
    """The parity run's figures, over its shared-causal instances.

    The paper figures are drawn from the federated results; parity makes them identical
    to the centralized ones. fed6 is not drawn here: this run predates the composition
    search, and the evaluation run's fed6 supersedes it.
    """
    eda.write_figures(root, out / "eda", logger=log, strict=True)
    frames = {}
    for name, sub, filename in [
        ("centralized", "fine_mapping", "fm_results.tsv"),
        ("federated", "fed_fine_mapping", "fed_fm_results.tsv"),
    ]:
        frame = pd.read_csv(root / "reports" / sub / filename, sep="\t")
        frames[name] = frame.loc[frame.causal_mode.eq("shared")].copy()
    paper_plots.write_figures(
        frames["federated"],
        out / "federated/shared",
        data_root=root,
        detail_dir=root / "reports/fed_fine_mapping/detail",
        logger=log,
        strict=True,
    )
    federation.write_figures(
        frames["centralized"],
        out / "federation/shared",
        federated=frames["federated"],
        data_root=root,
        logger=log,
        strict=True,
    )


def render_evaluation(root: Path, out: Path) -> None:
    headline = pd.read_csv(root / "fm_results_by_arm.tsv", sep="\t", low_memory=False)
    composition.evaluation_figures(root, headline, out)


def web_image(source: Path, destination: Path) -> None:
    from PIL import Image

    with Image.open(source) as image:
        image = image.convert("RGB")
        scale = min(1.0, MAX_WIDTH / image.width, MAX_HEIGHT / image.height)
        if scale < 1.0:
            size = (round(image.width * scale), round(image.height * scale))
            image = image.resize(size, Image.Resampling.LANCZOS)
        image.save(destination, "WEBP", lossless=True, quality=100, method=6)


# --------------------------------------------------------------------------- tables


def pct(value: float, digits: int = 1) -> str:
    return f"{value * 100:.{digits}f}%"


def interval(estimate: float, low: float, high: float, signed: bool = False) -> str:
    """``94.7% (94.3–95.1)``, or ``+16.6 (15.6–17.5)`` in percentage points."""
    if signed:
        return f"{estimate * 100:+.1f} ({low * 100:.1f}–{high * 100:.1f})"
    return f"{pct(estimate)} ({low * 100:.1f}–{high * 100:.1f})"


def describe(site_counts: dict[str, dict[str, int]]) -> str:
    """``Covenant 46,500 AFR · 1,500 EUR · 1,000 CSA; ANL 1,000 AFR``."""
    return "; ".join(
        f"{SITE_LABELS.get(site, site)} "
        + " · ".join(f"{n:,} {pop}" for pop, n in counts.items() if n)
        for site, counts in site_counts.items()
    )


def heldout_evaluation(root: Path) -> pd.DataFrame:
    summary = pd.read_csv(root / "evaluation_summary.tsv", sep="\t")
    raw = pd.read_csv(root / "fm_results_by_arm.tsv", sep="\t", low_memory=False)
    cohorts = composition.evaluation_cohorts(root)
    rows = []
    for arm, label in ARM_LABELS.items():
        metric = summary.loc[summary.arm.eq(arm)].set_index("metric")

        def cell(name, metric=metric):
            r = metric.loc[name]
            return interval(r.estimate, r.ci_low, r.ci_high)

        # Smallest credible set holding a causal variant, on the fits that have one.
        cs_size = raw.loc[raw.arm.eq(arm)].best_cs_size.mean()
        rows.append(
            {
                "Cohort": label,
                "Individuals": f"{cohorts[arm][0]:,}",
                "Power": cell("power"),
                "Yield, PIP > 0.95": cell("high_confidence_yield"),
                "Coverage": cell("cs_coverage"),
                "False discoveries": cell("fdr_pip95"),
                "Causal set size": f"{cs_size:.1f}",
            }
        )
    return pd.DataFrame(rows)


def paired_power_gain(root: Path) -> pd.DataFrame:
    paired = pd.read_csv(root / "paired_vs_full.tsv", sep="\t").set_index("arm")
    rows = []
    for arm, label in ARM_LABELS.items():
        if arm not in paired.index:
            continue
        r = paired.loc[arm]
        rows.append(
            {
                "Full federation compared with": label,
                "Power gain (pp)": interval(r.power_difference, r.ci_low, r.ci_high, signed=True),
                "Loci": int(r.n_loci),
                "Paired instances": f"{int(r.n_instances):,}",
            }
        )
    return pd.DataFrame(rows)


def composition_search(root: Path) -> pd.DataFrame:
    ranking = pd.read_csv(root / "development_ranking.tsv", sep="\t")
    design = json.loads((root / "design.json").read_text())
    selected = json.loads((root / "selection.json").read_text())["name"]
    counts = {c["name"]: c["site_counts"] for c in design["candidates"]}
    rows = []
    for rank, r in enumerate(ranking.itertuples(index=False), start=1):
        rows.append(
            {
                "Rank": rank,
                "Candidate": r.candidate,
                "Cohort of 50,000": describe(counts[r.candidate]),
                "Power": pct(r.power),
                "Yield, PIP > 0.95": pct(r.high_confidence_yield),
                "Coverage": pct(r.cs_coverage),
                "False discoveries": pct(r.fdr_pip95),
                "Outcome": "selected" if r.candidate == selected else "",
            }
        )
    return pd.DataFrame(rows)


def parity_benchmark(root: Path) -> pd.DataFrame:
    parity = json.loads((root / "parity.json").read_text())
    members = root / "reports/fed_fine_mapping/detail/fm_cs_members.tsv"
    with members.open() as handle:
        n_members = sum(1 for _ in handle) - 1
    rows = [
        ("Instances fine-mapped on both paths", f"{parity['n_instances']:,}"),
        ("Loci", f"{parity['n_loci']:,}"),
        ("Converged on both paths", f"{parity['n_converged_pairs']:,}"),
        ("Converged on one path only", f"{parity['n_nonconverged_pairs']:,}"),
        ("Instances that differ", f"{parity['n_mismatches']:,}"),
        ("Largest difference in any PIP", f"{parity['max_absolute_pip_difference']:g}"),
        ("Credible-set members compared", f"{n_members:,}"),
        ("Tolerance", f"rtol {parity['rtol']:g}, atol {parity['atol']:g}"),
        ("Verdict", "identical" if parity["pass"] else "differs"),
    ]
    return pd.DataFrame(rows, columns=["Check", "Federated vs pooled"])


def parity_arms(root: Path) -> pd.DataFrame:
    """The parity run's participation arms, with the design read from architecture_id."""
    arms = pd.read_csv(root / "arms/fm_results_by_arm.tsv", sep="\t", low_memory=False)
    arms = arms.loc[arms.causal_mode.eq("shared")].copy()
    design = arms.architecture_id.str.extract(r"ncsl(\d)_h2-([\d.]+)_rg([\d.]+)")
    arms["ncsl"] = design[0].astype(int)
    arms["h2"] = design[1].astype(float)
    arms["rg"] = design[2].astype(float)
    arms["opportunities"] = 1
    return arms


def power_by_architecture(root: Path) -> pd.DataFrame:
    arms = parity_arms(root)
    full = arms.loc[arms.arm.eq("federation")]
    rows = []
    for (h2, rg), cell in full.groupby(["h2", "rg"]):
        row = {"Per-locus h²": f"{h2:g}", "Cross-ancestry r_g": f"{rg:g}"}
        for ncsl, group in cell.groupby("ncsl"):
            row[f"{ncsl} causal"] = pct(group.any_causal_captured.mean())
        # Fewer than the design's 790 where withdrawn divergent instances were dropped.
        sizes = cell.groupby("ncsl").size()
        low, high = int(sizes.min()), int(sizes.max())
        row["Instances per cell"] = f"{low:,}" if low == high else f"{low:,}–{high:,}"
        rows.append(row)
    return pd.DataFrame(rows)


def weakest_signal_coverage(root: Path) -> pd.DataFrame:
    """Credible-set coverage at the faintest heritability, per participation arm."""
    arms = parity_arms(root)
    weak = arms.loc[arms.h2.eq(arms.h2.min()) & arms.rg.eq(1.0)]
    rows = []
    for arm in ("federation", "covenant", "mbzuai", "anl"):
        group = weak.loc[weak.arm.eq(arm)]
        rows.append(
            {
                "Cohort": ARM_LABELS[arm],
                f"Coverage at h² = {weak.h2.min():g}": interval(
                    *clustered_ratio(group, "n_cs_containing_causal", "n_credible_sets")
                ),
                "Credible sets": f"{int(group.n_credible_sets.sum()):,}",
                "Loci": group.locus_id.nunique(),
            }
        )
    return pd.DataFrame(rows)


TABLES = {
    "heldout_evaluation.tsv": ("evaluation", heldout_evaluation),
    "paired_power_gain.tsv": ("evaluation", paired_power_gain),
    "composition_search.tsv": ("evaluation", composition_search),
    "parity_benchmark.tsv": ("parity", parity_benchmark),
    "power_by_architecture.tsv": ("parity", power_by_architecture),
    "weakest_signal_coverage.tsv": ("parity", weakest_signal_coverage),
}


# --------------------------------------------------------------------------- main


def build(parity: Path, evaluation: Path, out: Path) -> list[Path]:
    runs = {"parity": parity, "evaluation": evaluation}
    for name, root in runs.items():
        if not root.is_dir():
            raise SystemExit(f"--{name}-run {root} is not a directory")
    with tempfile.TemporaryDirectory(prefix="site-results-") as scratch:
        rendered = {name: Path(scratch) / name for name in runs}
        render_parity(parity, rendered["parity"])
        render_evaluation(evaluation, rendered["evaluation"])
        missing = [p for run, p, _ in FIGURES if not (rendered[run] / p).is_file()]
        if missing:
            raise SystemExit(f"figures not drawn: {', '.join(missing)}")
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        written = []
        for run, path, name in FIGURES:
            destination = out / f"{name}.webp"
            web_image(rendered[run] / path, destination)
            written.append(destination)
    for filename, (run, make) in TABLES.items():
        destination = out / filename
        make(runs[run]).to_csv(destination, sep="\t", index=False)
        written.append(destination)
    return written


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--parity-run", type=Path, required=True)
    parser.add_argument("--evaluation-run", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    written = build(args.parity_run.resolve(), args.evaluation_run.resolve(), args.out)
    total = sum(p.stat().st_size for p in written)
    print(f"{len(written)} files, {total / 2**20:.1f} MiB -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
