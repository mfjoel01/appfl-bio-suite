#!/usr/bin/env python3
"""Build the GWAS and FLamby results shown in the project site's Results tab.

    python scripts/site_results.py gwas   --run local/runs/gwas-5site-fullscale
    python scripts/site_results.py flamby --run local/output/flamby-heart-disease-loopback \\
        --stamp 2026-09-14-17-56

Both runs are single-machine runs, and the captions in website/partners.json say so:

* GWAS: the five-site, full-variant-scale meta-analysis of 21 July 2026. Five simulated
  sites' summary statistics were combined by APPFL's MetaAnalysisAggregator on one node,
  not over Globus Compute. The aggregation is the federation's; the network is not.
* FLamby: a loopback run of Fed-Heart-Disease -- FedAvg over three of the dataset's four
  centres, five rounds, every site a process on one machine. No cross-institution FLamby
  run has completed yet; when one does, point this script at it.

Figures are scaled to the width the viewer displays and stored as lossless WebP, as
scripts/fine-mapping/site_results.py does: the viewer embeds every result in the map's
metadata, so each figure's bytes are part of every load of the map. Tables are reduced to
readable columns. No path, host name or account from the run reaches the output.

Writes website/results/gwas/ and website/results/flamby-heart-disease/, which this script
owns and replaces. website/partners.json selects and captions the files.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

MAX_WIDTH, MAX_HEIGHT = 1800, 2000
GENOME_WIDE = 5e-8
TOP_HITS = 15


def web_image(source: Path, destination: Path) -> None:
    from PIL import Image

    with Image.open(source) as image:
        image = image.convert("RGB")
        scale = min(1.0, MAX_WIDTH / image.width, MAX_HEIGHT / image.height)
        if scale < 1.0:
            size = (round(image.width * scale), round(image.height * scale))
            image = image.resize(size, Image.Resampling.LANCZOS)
        image.save(destination, "WEBP", lossless=True, quality=100, method=6)


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def replace_dir(target: Path, build) -> None:
    """Build into a scratch directory, then swap it in: a failed build changes nothing."""
    with tempfile.TemporaryDirectory(dir=target.parent) as scratch:
        staged = Path(scratch) / target.name
        staged.mkdir()
        build(staged)
        if target.exists():
            shutil.rmtree(target)
        shutil.move(str(staged), target)


# --------------------------------------------------------------------------- GWAS


def gwas(run: Path, out: Path) -> None:
    graphs, data = run / "server_output_5site" / "graphs", run / "server_output_5site" / "data"
    for trait in ("bmi", "t2d"):
        for kind in ("manhattan", "qq"):
            web_image(graphs / f"appfl_meta_gwas_{trait}_{kind}.png", out / f"{trait}_{kind}.webp")

    with (data / "appfl_site_pgs_metrics.csv").open(encoding="utf-8") as handle:
        sites = list(csv.DictReader(handle))
    with (data / "appfl_meta_summary.csv").open(encoding="utf-8") as handle:
        total = next(csv.DictReader(handle))
    rows = [
        [s["CLIENT_ID"], f"{int(s['GWAS_N']):,}", f"{float(s['LOCAL_BMI_R2']):.3f}",
         f"{float(s['LOCAL_T2D_AUROC']):.3f}"]
        for s in sites
    ]  # fmt: skip
    rows.append(
        ["All sites (N-weighted)", f"{int(total['TOTAL_GWAS_N']):,}",
         f"{float(total['WEIGHTED_LOCAL_BMI_R2']):.3f}",
         f"{float(total['WEIGHTED_LOCAL_T2D_AUROC']):.3f}"]
    )  # fmt: skip
    write_tsv(out / "site_cohorts.tsv", ["Site", "Individuals", "BMI R²", "T2D AUROC"], rows)

    with (data / "appfl_meta_gwas_hits.csv").open(encoding="utf-8") as handle:
        hits = [h for h in csv.DictReader(handle) if float(h["P"]) < GENOME_WIDE]
    rows = []
    for trait in ("BMI", "T2D"):
        mine = sorted((h for h in hits if h["TRAIT"] == trait), key=lambda h: float(h["P"]))
        for h in mine[:TOP_HITS]:
            effect = f"OR {float(h['OR']):.3f}" if h.get("OR") else f"β {float(h['BETA']):+.3f}"
            rows.append(
                [trait, h["CHR"], f"{int(h['BP']):,}", f"{h['EA']}/{h['NEA']}", effect,
                 f"{float(h['SE']):.4f}", f"{float(h['P']):.1e}", f"{float(h['MAF']):.3f}"]
            )  # fmt: skip
    write_tsv(
        out / "top_hits.tsv",
        ["Trait", "Chr", "Position (GRCh37)", "Effect/other allele", "Effect", "SE", "P", "MAF"],
        rows,
    )


# --------------------------------------------------------------------------- FLamby

_ANSI = re.compile(r"\x1b\[[0-9;]*m")
_ROW = re.compile(r"\]:\s+(\d+)\s+([YN])\s+([\d.\s]+)$")
_SIZE = re.compile(r"\[loopback\] (Site\d+): n = (\d+)")


def flamby_rounds(run: Path, stamp: str) -> tuple[dict[str, int], dict[str, dict]]:
    """Per site: training-set size, and per-round metrics from the site's own log.

    A ``Y`` row is the global model the site received, evaluated before local training;
    an ``N`` row is the same round after local training.
    """
    sizes: dict[str, int] = {}
    for server in run.parent.glob(f"loopback/*Server_{stamp}*.txt"):
        sizes.update({s: int(n) for s, n in _SIZE.findall(_ANSI.sub("", server.read_text()))})
    sites: dict[str, dict] = {}
    for log in sorted(run.glob(f"sites/*/appfl_*_{stamp}-*.txt")):
        site = log.parent.name
        received, trained = {}, {}
        for line in _ANSI.sub("", log.read_text()).splitlines():
            match = _ROW.search(line)
            if not match:
                continue
            round_no, pre, values = int(match[1]), match[2], [float(v) for v in match[3].split()]
            if pre == "Y":
                received[round_no] = {"loss": values[0], "acc": values[1]}
            else:
                trained[round_no] = {"train_loss": values[1], "train_acc": values[2],
                                     "loss": values[3], "acc": values[4]}  # fmt: skip
        sites[site] = {"received": received, "trained": trained}
    if not sites:
        raise SystemExit(f"no site logs for run {stamp} under {run}/sites")
    return sizes, sites


def flamby(run: Path, stamp: str, out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from appfl_bio_suite.core.plot_style import CATEGORICAL, apply_style

    apply_style()
    sizes, sites = flamby_rounds(run, stamp)
    names = sorted(sites)
    rounds = sorted(next(iter(sites.values()))["trained"])

    # Every site validates on the same pooled test set, so the global model each one
    # received scores the same everywhere; take it from the first site and say so.
    first = sites[names[0]]["received"]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), dpi=200)
    for axis, key, label in ((axes[0], "acc", "Accuracy"), (axes[1], "loss", "Loss")):
        received = sorted(first)
        axis.plot(
            [r + 1 for r in received],
            [first[r][key] for r in received],
            color="#222222",
            lw=2.2,
            marker="o",
            label="Global model (start of round)",
        )
        for colour, site in zip(CATEGORICAL, names, strict=False):
            trained = sites[site]["trained"]
            n = f", n = {sizes[site]}" if site in sizes else ""
            axis.plot(
                [r + 1 for r in rounds],
                [trained[r][key] for r in rounds],
                color=colour,
                lw=1.4,
                ls="--",
                marker=".",
                label=f"{site} after local training{n}",
            )
        axis.set_xlabel("Round")
        axis.set_ylabel(f"{label} on the pooled test set")
        axis.set_xticks([r + 1 for r in rounds])
    axes[0].legend(loc="lower right", fontsize=7, frameon=False)
    fig.suptitle(
        "Fed-Heart-Disease, FedAvg over three centres (loopback, one machine)", fontsize=10
    )
    fig.tight_layout()
    with tempfile.TemporaryDirectory() as scratch:
        png = Path(scratch) / "rounds.png"
        fig.savefig(png)
        plt.close(fig)
        web_image(png, out / "fedavg_rounds.webp")

    rows = []
    for r in rounds:
        row = [str(r + 1), f"{first[r]['acc']:.3f}" if r in first else ""]
        row += [f"{sites[s]['trained'][r]['acc']:.3f}" for s in names]
        rows.append(row)
    header = ["Round", "Global model accuracy (start of round)"]
    header += [
        f"{s} after local training" + (f" (n = {sizes[s]})" if s in sizes else "") for s in names
    ]
    write_tsv(out / "fedavg_rounds.tsv", header, rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    g = commands.add_parser("gwas")
    g.add_argument("--run", type=Path, required=True)
    f = commands.add_parser("flamby")
    f.add_argument("--run", type=Path, required=True)
    f.add_argument("--stamp", required=True, help="the run's timestamp, e.g. 2026-09-14-17-56")
    args = parser.parse_args(argv)

    results = REPO / "website" / "results"
    if args.command == "gwas":
        replace_dir(results / "gwas", lambda out: gwas(args.run, out))
        print(f"wrote {results / 'gwas'}")
    else:
        target = results / "flamby-heart-disease"
        replace_dir(target, lambda out: flamby(args.run, args.stamp, out))
        print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
