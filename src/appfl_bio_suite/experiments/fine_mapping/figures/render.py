"""Render every fine-mapping figure from whatever inputs are on disk. COORDINATOR-SIDE.

One command rather than four, because the four modules want overlapping inputs and
working out which ones a given run can support is a job for code, not for the person
running it. Anything unavailable is reported and skipped; nothing here fails a batch
because one input is missing.

    python -m appfl_bio_suite.experiments.fine_mapping.figures.render \\
        --results local/papers/fine-mapping/results/fine_mapping/fm_results.tsv \\
        --data-root local/data/fine-mapping \\
        --out-dir local/output/fine-mapping/figures

Optional, and each unlocks specific figures:

    --detail-dir     output of ``figures.detail`` -- true per-credible-set coverage
                     (pap2c, pap4), the population-specific causal probabilities
                     (pap10), the estimated effect-size panels (pap12), and the
                     LocusZoom-style exemplars (pap11)
    --federated      a federated results table over the SAME package -- the parity
                     figure (fed1), which is this experiment's central claim
    --run-log        a federated run log -- per-site recoding counts (fed4) and the
                     dropped-variant precondition fed1 reports alongside its result
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from appfl_bio_suite.experiments.fine_mapping.figures import (
    eda as eda_mod,
)
from appfl_bio_suite.experiments.fine_mapping.figures import (
    federation as fed_mod,
)
from appfl_bio_suite.experiments.fine_mapping.figures import (
    paper_plots as paper_mod,
)

log = logging.getLogger(__name__)


def arm_cohort_sizes(data_root, by_arm) -> dict[str, tuple[int, int]]:
    """``{arm: (individuals, ancestry columns)}``, for annotating the arm axis.

    Read back from each arm's own config rather than assumed, because the down-sampled
    arm's cohort is a draw and its column count can differ from the full federation's --
    at one third of the package EAS falls below ``min_gwas_n`` and the arm has five
    columns, not six. Annotating it as six would misdescribe the run.
    """
    from appfl_bio_suite.experiments.fine_mapping.arms import ARMS_BY_NAME

    out: dict[str, tuple[int, int]] = {}
    for arm_name in by_arm["arm"].unique():
        arm = ARMS_BY_NAME.get(arm_name)
        if not data_root:
            continue
        try:
            import yaml

            cfg_p = Path(data_root)
            # arms write their configs beside the merged table; fall back to counting
            # nothing rather than guessing.
            for cand in (
                cfg_p.parent / "arms" / f"{arm_name}.yaml",
                cfg_p / "arms" / f"{arm_name}.yaml",
            ):
                if cand.exists():
                    raw = yaml.safe_load(cand.read_text())
                    sites = arm.sites if arm and arm.ld_from else tuple(raw["sites"])
                    pooled: dict[str, int] = {}
                    for s in sites:
                        for pop, n in raw["sites"][s]["composition"].items():
                            pooled[pop] = pooled.get(pop, 0) + n
                    floor = raw["fine_mapping"]["min_gwas_n"]
                    out[arm_name] = (
                        sum(v for v in pooled.values() if v >= floor),
                        sum(1 for v in pooled.values() if v >= floor),
                    )
                    break
        except Exception:  # noqa: BLE001 - an annotation must not lose the figure
            continue
    return out


def render_all(
    results: Path,
    out_dir: Path,
    data_root: Path | None = None,
    detail_dir: Path | None = None,
    federated: Path | None = None,
    run_log: Path | None = None,
    by_arm_results: Path | None = None,
    skip: tuple[str, ...] = (),
    logger=None,
) -> dict[str, list[Path]]:
    """Draw every group of figures the inputs support; return what was written."""
    active = logger or log
    out_dir = Path(out_dir)
    written: dict[str, list[Path]] = {}

    df = pd.read_csv(results, sep="\t")
    active.info(
        "results: %d instance(s), %d architecture(s), %d loci",
        len(df),
        df["architecture_id"].nunique(),
        df["locus_id"].nunique(),
    )

    if "eda" not in skip and data_root:
        active.info("--- exploratory figures ---")
        written["eda"] = eda_mod.write_figures(data_root, out_dir / "eda", active)
    elif "eda" not in skip:
        active.info("skipping the EDA group: --data-root not given")

    if "paper" not in skip:
        active.info("--- SuSiEx-paper analogues ---")
        written["paper"] = paper_mod.write_figures(
            df, out_dir / "paper", data_root, detail_dir, logger=active
        )

    if "federation" not in skip:
        active.info("--- federation figures ---")
        fed = pd.read_csv(federated, sep="\t") if federated and Path(federated).exists() else None
        by_arm = cohort = None
        if by_arm_results and Path(by_arm_results).exists():
            by_arm = pd.read_csv(by_arm_results, sep="\t")
            active.info(
                "arms: %s",
                ", ".join(f"{a} ({n:,})" for a, n in by_arm["arm"].value_counts().items()),
            )
            cohort = arm_cohort_sizes(data_root, by_arm)
        flips = None
        cache = out_dir / "cache" / "allele_flips.tsv"
        if cache.exists():
            flips = pd.read_csv(cache, sep="\t")
        written["federation"] = fed_mod.write_figures(
            df, out_dir / "federation", fed, data_root, run_log, flips, by_arm, cohort, active
        )

    total = sum(len(v) for v in written.values())
    active.info("=== %d figure(s) in total -> %s ===", total, out_dir)
    for group, paths in written.items():
        active.info("    %-11s %d", group, len(paths))
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Render every fine-mapping figure.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--results", required=True, help="per-instance results TSV")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--data-root", default=None, help="the simulated package")
    ap.add_argument("--detail-dir", default=None, help="output of figures.detail")
    ap.add_argument("--federated", default=None, help="federated results over the same package")
    ap.add_argument("--run-log", default=None, help="a federated run log")
    ap.add_argument(
        "--by-arm", default=None, help="fm_results_by_arm.tsv from the participation arms"
    )
    ap.add_argument("--skip", default="", help="comma-separated: eda,paper,federation")
    args = ap.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )
    p = Path(args.results)
    if not p.exists():
        print(
            f"ERROR: {p} not found -- no fine-mapping results to plot yet.\n"
            "Produce some with:\n"
            "    scripts/fine-mapping/run_stage.py centralized --config <pipeline_config.yaml>\n"
            "or, for a federated run:\n"
            "    appfl-bio-suite run fine-mapping --config loopback --driver serial",
            file=sys.stderr,
        )
        return 2
    render_all(
        p,
        Path(args.out_dir),
        Path(args.data_root) if args.data_root else None,
        Path(args.detail_dir) if args.detail_dir else None,
        Path(args.federated) if args.federated else None,
        Path(args.run_log) if args.run_log else None,
        Path(args.by_arm) if args.by_arm else None,
        tuple(s.strip() for s in args.skip.split(",") if s.strip()),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
