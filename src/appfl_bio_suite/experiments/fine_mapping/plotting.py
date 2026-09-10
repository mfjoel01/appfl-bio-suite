"""The aggregator's in-process figure entry point. COORDINATOR-SIDE.

WHY THIS IS NOW A THIN WRAPPER
-------------------------------
This module used to hold four figures of its own, ported from the standalone
repository's ``scripts/plot_fine_mapping.py``. They have moved into
``figures/paper_plots.py`` and grown into a full set that follows the SuSiEx paper's
own forms, alongside two more groups -- ``figures/eda.py`` for the simulated package
and ``figures/federation.py`` for the cost and correctness of federating.

The entry point stays here, and stays named ``write_figures``, because the aggregator
imports it by that name at the end of a run: a federated sweep that took hours should
produce its diagnostics in the same process, not leave someone to remember a second
command. Nothing about that call site changes.

The old figure names map onto the new ones as follows, for anyone holding a link to
an earlier run's output:

    fig1_power_grid.png     -> pap1_recall_grid.png
    fig2_cs_size.png        -> pap5_resolution.png       (panel b)
    fig3_causal_pip.png     -> pap5_resolution.png       (panel a)
    fig4_stratum_power.png  -> pap7_ld_divergence.png    (panel a)

For the full set -- including the EDA and federation groups, and the figures that need
harvested per-credible-set detail -- use the batch renderer instead::

    python -m appfl_bio_suite.experiments.fine_mapping.figures.render \\
        --results <fm_results.tsv> --data-root <package> --out-dir <figures/>
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from appfl_bio_suite.experiments.fine_mapping.figures.paper_plots import (
    parse_architecture,
)
from appfl_bio_suite.experiments.fine_mapping.figures.paper_plots import (
    write_figures as _write_paper_figures,
)

__all__ = ["write_figures", "parse_architecture", "main"]

log = logging.getLogger(__name__)


def write_figures(results, out_dir, logger=None, data_root=None, detail_dir=None):
    """Render the paper-analogue figures from a results DataFrame.

    Kept deliberately narrow -- this is the *aggregator's* call, made at the end of a
    federated run where the coordinator holds the results table and usually nothing
    else. The EDA group needs the simulated package and the federation group needs a
    second results table, neither of which the coordinator necessarily has, so both are
    left to the batch renderer.

    A figure that cannot be drawn is skipped with a note rather than failing the run:
    the numbers are the deliverable, and a plotting problem must never lose a
    fine-mapping sweep that took hours.
    """
    return _write_paper_figures(
        results, out_dir, data_root=data_root, detail_dir=detail_dir, logger=logger or log
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Render the fine-mapping figures from a results table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--results", default="local/output/fine-mapping/data/fed_fm_results.tsv")
    ap.add_argument("--out-dir", default=None, help="default: <results dir>/../graphs")
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--detail-dir", default=None)
    args = ap.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
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
        "loaded %d instance(s), %d architecture(s)", len(frame), frame["architecture_id"].nunique()
    )
    write_figures(frame, out_dir, data_root=args.data_root, detail_dir=args.detail_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
