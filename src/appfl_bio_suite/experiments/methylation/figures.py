"""Analysis figures for a completed synthetic run. COORDINATOR-SIDE ONLY.

    python -m appfl_bio_suite.experiments.methylation.figures --run DIR [--curves DIR] --out DIR

``--run`` is a completed test run of a synthetic scenario. ``--curves`` is an optional
validation-only run of the same scenario with ``track_validation`` set, such as the
packaged ``uneven-curves`` config for ``uneven``. It must share the run's cohort, splits
and sites, and its training history must match the run's round for round, so that its
scores at the run's last round describe exactly the models the run tested.

Figures, each as PNG and SVG:

    01_who_gains         every site alone, against the federated and centralized models
    02_training_rounds   validation score after every round              (needs --curves)
    03_sparse_input      validation score against observed CpGs           (needs --curves)
    04_class_recall      which classes every model gets right
    05_site_composition  each site's share of every class
    A1_data_overview     training patients on two principal components, one class a panel
    A2_cpg_methylation   class-mean methylation at the most contrasting probes

Intervals are 95% bootstrap intervals computed from each confusion matrix, resampling
outcomes within each true class, as the stratified split was drawn. Rebuilding from the
same runs reproduces every file byte for byte, given the same fonts and package versions:
the PCA runs in float64, and SVG output carries no timestamp or random identifiers.
``figure_data.json`` records every plotted value and every file's hash.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

# Set before pyplot is imported: a headless node has no display.
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from appfl_bio_suite.core import plot_style as brand  # noqa: E402
from appfl_bio_suite.core.simulation import file_checksum  # noqa: E402

from .dataset import MethylationData  # noqa: E402

# Two identity hues, validated together against the white surface (lightness band,
# chroma floor, CVD and normal-vision separation all pass; cyan's 2.96:1 contrast is
# why every series also carries a legend entry or a direct label). Local sites are
# context and stay gray; the smallest site is emphasized by weight, not by a third hue.
CENTRALIZED = brand.COLORS["anl-blue-brand"]
FEDERATED = brand.COLORS["anl-cyan"]
FOCUS = brand.COLORS["anl-gray-700"]
CONTEXT = brand.COLORS["anl-gray-500"]
FAINT = brand.COLORS["anl-gray-300"]
INK = brand.COLORS["anl-ink"]
NOTE = brand.COLORS["anl-gray-500"]
CHECKED = ("metrics.json", "preprocessing.json", "cohort.npz", "config.json", "history.json")
# The only settings a curves run may change: it trains the same models for longer and
# scores them on validation patients at more coverage levels.
CURVES_ONLY = {"rounds", "coverage_counts", "track_validation", "validation_only"}
RC = {
    "font.size": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    # A fixed salt makes SVG element IDs deterministic instead of random.
    "svg.hashsalt": "appfl-bio-suite",
}
NUMBER_WORDS = dict(enumerate("zero one two three four five six seven eight nine ten".split()))


def interval(confusion, *, draws: int = 4000, seed: int = 0) -> tuple[float, float]:
    """95% interval for balanced accuracy, resampling each true class's outcomes."""
    matrix = np.asarray(confusion, dtype=np.int64)
    rng = np.random.default_rng(seed)
    recalls = [
        rng.binomial(row.sum(), row[c] / row.sum(), size=draws) / row.sum()
        for c, row in enumerate(matrix)
        if row.sum()
    ]
    low, high = np.percentile(np.mean(recalls, axis=0), [2.5, 97.5])
    return float(low), float(high)


def _verified(run: Path, names) -> dict:
    manifest = json.loads((run / "manifest.json").read_text())
    for name in names:
        if file_checksum(run / name) != manifest["outputs"][name]:
            raise ValueError(f"{run / name} does not match its run manifest")
    return manifest


def _check_curves(curves: Path, run: Path, manifest: dict, config: dict) -> dict:
    """Prove the curves run trained the tested models before using its scores."""
    # What kind of run it is comes first: a test run has no per-round scores to verify.
    kind = json.loads((curves / "manifest.json").read_text())
    if kind.get("evaluation_split") != "validation":
        raise ValueError(f"{curves} predicted test patients; a curves run must be validation-only")
    if "validation_history.json" not in kind["outputs"]:
        raise ValueError(f"{curves} has no per-round validation scores; set track_validation")
    curves_manifest = _verified(curves, (*CHECKED, "validation_history.json"))
    for name in ("cohort.npz", "preprocessing.json"):
        if curves_manifest["outputs"][name] != manifest["outputs"][name]:
            raise ValueError(f"{curves} has a different {name} from {run}")
    curves_config = json.loads((curves / "config.json").read_text())
    changed = {
        key
        for key in set(config) | set(curves_config)
        if config.get(key) != curves_config.get(key) and key not in CURVES_ONLY
    }
    if changed:
        raise ValueError(f"{curves} trains differently from {run}: {', '.join(sorted(changed))}")
    history = json.loads((run / "history.json").read_text())
    curves_history = json.loads((curves / "history.json").read_text())
    if curves_history[: len(history)] != history:
        raise ValueError(f"{curves} did not follow {run}'s training round for round")
    if max(config["coverage_counts"]) not in curves_config["coverage_counts"]:
        raise ValueError(f"{curves} does not score the run's highest coverage")
    return json.loads((curves / "validation_history.json").read_text())


def build(run: Path, out: Path, curves: Path | None = None) -> dict:
    """Write every figure for ``run`` into the new directory ``out``; return the record."""
    kind = json.loads((run / "manifest.json").read_text())
    if kind.get("input_kind") != "synthetic":
        raise ValueError("analysis figures are drawn only for synthetic runs")
    if kind.get("evaluation_split") != "test":
        raise ValueError("run evaluated validation patients only; there is no test result")
    manifest = _verified(run, CHECKED)
    config = json.loads((run / "config.json").read_text())
    tracked = _check_curves(curves, run, manifest, config) if curves else None
    out.mkdir(parents=True, exist_ok=False)

    data = MethylationData.load(run / "cohort.npz")
    prep = json.loads((run / "preprocessing.json").read_text())
    metrics = json.loads((run / "metrics.json").read_text())
    highest = max(row["requested_cpgs"] for row in metrics)
    scores = {row["model"]: row for row in metrics if row["requested_cpgs"] == highest}

    font = brand.apply_style()
    with plt.rc_context(RC):
        figure = _Figures(out, data, prep, config, scores, highest, tracked)
        record = figure.draw_all()
    record.update(
        {
            "run": str(run),
            "run_commit": manifest["suite_commit"],
            "curves_run": str(curves) if curves else None,
            "input_kind": manifest["input_kind"],
            "config": config,
            "intervals": "95% bootstrap, 4000 draws, outcomes resampled within each true class",
            "generator_source_sha256": manifest["source_sha256"]["dataset.py"],
            "figures_source_sha256": file_checksum(Path(__file__)),
            "font": font,
        }
    )
    record["files"] = {p.name: file_checksum(p) for p in sorted(out.iterdir()) if p.is_file()}
    (out / "figure_data.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


class _Figures:
    def __init__(self, out, data, prep, config, scores, highest, tracked):
        self.out, self.prep, self.config = out, prep, config
        self.scores, self.highest, self.tracked = scores, highest, tracked
        self.sites = list(prep["sites"])
        self.small = min(self.sites, key=lambda site: len(prep["sites"][site]))
        self.labels = prep["classes"]
        self.letters = [chr(ord("A") + i) for i in range(len(self.labels))]
        self.chance = 100 / len(self.labels)
        train = np.array(prep["splits"]["train"])
        self.panel = np.array(prep["panel_indices"])
        self.cpg_ids = data.cpg_ids
        self.X = data.X[np.ix_(train, self.panel)]
        self.y = data.y[train]
        self.names, self.record = [], {}

    # -- shared anatomy -----------------------------------------------------------------

    def save(self, fig, name, note):
        fig.text(0.06, 0.025, note, fontsize=9, color=NOTE)
        fig.savefig(self.out / f"{name}.png", dpi=190)
        svg = self.out / f"{name}.svg"
        fig.savefig(svg, metadata={"Date": None})
        svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
        self.names.append(name)
        plt.close(fig)

    @staticmethod
    def heading(fig, title, subtitle, y=0.955):
        fig.suptitle(title, x=0.06, ha="left", y=y, fontsize=18, weight="bold")
        below = y - 0.48 / fig.get_figheight()
        fig.text(0.06, below, subtitle, fontsize=11.5, color=brand.COLORS["anl-gray-700"])

    def chance_line(self, ax, orientation="h"):
        line = ax.axhline if orientation == "h" else ax.axvline
        line(self.chance, color=FAINT, linewidth=1, zorder=0)

    def chance_words(self):
        count = len(self.labels)
        return f"{NUMBER_WORDS.get(count, count)}-class chance, {self.chance:.1f}%".capitalize()

    def site_label(self, site):
        return f"{site} alone ({len(self.prep['sites'][site]):,} patients)"

    def tracked_scores(self, round_number, k):
        entry = next(e for e in self.tracked if e["round"] == round_number)
        return {s["model"]: s for s in entry["scores"] if s["requested_cpgs"] == k}

    def draw_all(self):
        self.who_gains()
        if self.tracked:
            self.training_rounds()
            self.sparse_input()
        self.class_recall()
        self.site_composition()
        self.data_overview()
        self.cpg_methylation()
        return {"figures": self.names, **self.record}

    # -- 01: who gains from joining -----------------------------------------------------

    def who_gains(self):
        values = {
            name: (
                100 * row["balanced_accuracy"],
                *(100 * v for v in interval(row["confusion_matrix"])),
            )
            for name, row in self.scores.items()
        }
        self.record["who_gains"] = {
            name: {"balanced_accuracy": v[0], "interval": [v[1], v[2]]}
            for name, v in values.items()
        }
        fig, ax = plt.subplots(figsize=(11, 6.4))
        fig.subplots_adjust(left=0.25, right=0.86, bottom=0.17, top=0.72)
        for name, color in (("centralized", CENTRALIZED), ("federated", FEDERATED)):
            score, low, high = values[name]
            ax.axvspan(low, high, color=color, alpha=0.12, linewidth=0, zorder=0)
            ax.axvline(score, color=color, linewidth=2, zorder=1)
        self.chance_line(ax, "v")
        federated = values["federated"][0]
        rows = list(reversed(self.sites))
        for y, site in enumerate(rows):
            score, low, high = values[f"local_{site}"]
            ax.plot([score, federated], [y, y], color=FAINT, linewidth=2, zorder=2)
            ax.plot([low, high], [y, y], color=CONTEXT, linewidth=2, zorder=3)
            ax.scatter(score, y, s=80, color=CONTEXT, edgecolors="white", linewidths=2, zorder=4)
            # Beside the whisker when a label above would sit on a reference line.
            crowded = min(abs(score - values[n][0]) for n in ("federated", "centralized")) < 12
            if crowded:
                ax.text(low - 1.2, y, f"{score:.1f}%", ha="right", va="center", fontsize=10.5)
            else:
                ax.text(score, y + 0.2, f"{score:.1f}%", ha="center", va="bottom", fontsize=10.5)
            gain = federated - score
            ax.text(
                1.03,
                y,
                f"{gain:+.1f}".replace("-", "\N{MINUS SIGN}"),
                transform=ax.get_yaxis_transform(),
                va="center",
                fontsize=13,
                weight="bold",
                color=INK,
            )
        ax.text(
            1.03,
            len(rows) - 0.35,
            "Federated\nminus alone",
            transform=ax.get_yaxis_transform(),
            va="bottom",
            fontsize=10,
            color=NOTE,
        )
        ax.set(
            xlim=(0, 100),
            ylim=(-0.6, len(rows) - 0.4),
            yticks=range(len(rows)),
            yticklabels=[self.site_label(site) for site in rows],
            xlabel="Balanced accuracy on the shared test set (%)",
        )
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        centralized = values["centralized"]
        handles = [
            Line2D(
                [],
                [],
                color=CENTRALIZED,
                linewidth=2,
                label=f"Centralized model {centralized[0]:.1f}%",
            ),
            Line2D([], [], color=FEDERATED, linewidth=2, label=f"Federated model {federated:.1f}%"),
            Line2D(
                [], [], color=CONTEXT, marker="o", linestyle="", markersize=8, label="Site alone"
            ),
            Line2D([], [], color=FAINT, linewidth=1, label=self.chance_words()),
        ]
        ax.legend(
            handles=handles,
            ncol=4,
            loc="lower left",
            bbox_to_anchor=(-0.33, 1.03),
            frameon=False,
            fontsize=10.5,
            handlelength=1.6,
            columnspacing=1.4,
        )
        self.heading(
            fig,
            "Which sites gain from joining the federation",
            f"Same {len(self.prep['splits']['test'])} held-out patients · "
            f"{self.highest:,} CpGs · {self.config['rounds']} rounds · "
            "shaded and whiskers: 95% intervals",
        )
        self.save(
            fig,
            "01_who_gains",
            "Synthetic data. One centralized and one federated model serve every site, "
            "so each is a line.",
        )

    # -- 02: score after every round ----------------------------------------------------

    def training_rounds(self):
        k = self.highest
        rounds = [entry["round"] for entry in self.tracked]
        series = {name: [] for name in self.scores}
        bands = {"centralized": [], "federated": []}
        for r in rounds:
            at = self.tracked_scores(r, k)
            for name in series:
                series[name].append(100 * at[name]["balanced_accuracy"])
            for name in bands:
                bands[name].append([100 * v for v in interval(at[name]["confusion_matrix"])])
        self.record["training_rounds"] = {
            "requested_cpgs": k,
            "rounds": rounds,
            "balanced_accuracy": series,
            "intervals": bands,
        }
        fig, ax = plt.subplots(figsize=(11, 6.4))
        fig.subplots_adjust(left=0.08, right=0.80, bottom=0.17, top=0.72)
        self.chance_line(ax)
        tested = self.config["rounds"]
        ax.axvline(tested, color=FAINT, linewidth=1, zorder=0)
        ax.text(
            tested + 0.3,
            3,
            f"{tested} rounds: the\ntested configuration",
            fontsize=10,
            color=NOTE,
            va="bottom",
        )
        for site in self.sites:
            focus = site == self.small
            ax.plot(
                rounds,
                series[f"local_{site}"],
                color=FOCUS if focus else CONTEXT,
                linewidth=2 if focus else 1,
                zorder=2,
            )
        for name, color in (("centralized", CENTRALIZED), ("federated", FEDERATED)):
            low, high = np.array(bands[name]).T
            ax.fill_between(rounds, low, high, color=color, alpha=0.12, linewidth=0, zorder=1)
            ax.plot(rounds, series[name], color=color, linewidth=2, zorder=3)
        ends = sorted(((series[f"local_{s}"][-1], s) for s in self.sites), reverse=True)
        self._end_labels(ax, rounds[-1], [(v, self.site_label(s)) for v, s in ends])
        ax.set(
            xlim=(1, rounds[-1]),
            ylim=(0, 102),
            xticks=[1, 5, 10, 15, 20, 25, 30][: len(rounds)],
            xlabel="Communication round",
            ylabel="Balanced accuracy (%)",
        )
        handles = [
            Line2D([], [], color=CENTRALIZED, linewidth=2, label="Centralized"),
            Line2D([], [], color=FEDERATED, linewidth=2, label="Federated"),
            Line2D([], [], color=FOCUS, linewidth=2, label=f"Smallest site alone ({self.small})"),
            Line2D([], [], color=CONTEXT, linewidth=1, label="Other sites alone"),
            Line2D([], [], color=FAINT, linewidth=1, label=self.chance_words()),
        ]
        ax.legend(
            handles=handles,
            ncol=5,
            loc="lower left",
            bbox_to_anchor=(0, 1.03),
            frameon=False,
            fontsize=10.5,
            handlelength=1.6,
            columnspacing=1.2,
        )
        final = {name: values[-1] for name, values in series.items()}
        swing = series[f"local_{self.small}"][1:]
        self.heading(
            fig,
            "With more rounds, federated catches up with centralized",
            f"{self.validation_size()} validation patients, scored after every round · "
            f"{k:,} CpGs · test patients never scored · shaded: 95% intervals",
        )
        self.save(
            fig,
            "02_training_rounds",
            f"Synthetic data. After round {rounds[-1]}: centralized {final['centralized']:.1f}%, "
            f"federated {final['federated']:.1f}%. From round 2 on, {self.small} alone ranges "
            f"{min(swing):.0f}–{max(swing):.0f}% from one round to the next.",
        )

    def validation_size(self):
        return len(self.prep["splits"]["validation"])

    @staticmethod
    def _end_labels(ax, x, items, gap=4.5):
        """Right-edge labels with leader lines, spaced so converging lines stay legible."""
        placed = []
        for value, text in items:  # descending by value
            y = value if not placed else min(value, placed[-1] - gap)
            placed.append(y)
            ax.annotate(
                text,
                xy=(x, value),
                xytext=(x + 0.8, y),
                textcoords="data",
                va="center",
                fontsize=10,
                color=INK,
                annotation_clip=False,
                arrowprops={
                    "arrowstyle": "-",
                    "color": FAINT,
                    "linewidth": 0.8,
                    "shrinkA": 0,
                    "shrinkB": 0,
                },
            )

    # -- 03: score against observed CpGs ------------------------------------------------

    def sparse_input(self):
        counts = sorted({s["requested_cpgs"] for s in self.tracked[0]["scores"]})
        tested, last = self.config["rounds"], self.tracked[-1]["round"]
        floor = self.config["min_coverage"] * len(self.panel)
        fig, axes = plt.subplots(1, 2, figsize=(12, 6.2), sharey=True)
        fig.subplots_adjust(left=0.08, right=0.97, bottom=0.17, top=0.70, wspace=0.08)
        record = {"requested_cpgs": counts, "training_minimum_cpgs": floor}
        for ax, round_number in zip(axes, (tested, last), strict=True):
            at = {k: self.tracked_scores(round_number, k) for k in counts}
            series = {
                name: [100 * at[k][name]["balanced_accuracy"] for k in counts]
                for name in self.scores
            }
            record[f"round_{round_number}"] = series
            if floor > counts[0]:
                ax.axvspan(
                    counts[0] / 1.3,
                    floor,
                    color=brand.COLORS["anl-gray-100"],
                    linewidth=0,
                    zorder=0,
                )
                ax.text(
                    np.sqrt(counts[0] * floor) / 1.15,
                    97,
                    "Below training\ncoverage",
                    ha="center",
                    va="top",
                    fontsize=9.5,
                    color=NOTE,
                )
            self.chance_line(ax)
            for site in self.sites:
                focus = site == self.small
                ax.plot(
                    counts,
                    series[f"local_{site}"],
                    color=FOCUS if focus else CONTEXT,
                    linewidth=2 if focus else 1,
                    marker="o" if focus else None,
                    markersize=4,
                    zorder=2,
                )
            for name, color in (("centralized", CENTRALIZED), ("federated", FEDERATED)):
                low, high = np.array(
                    [[100 * v for v in interval(at[k][name]["confusion_matrix"])] for k in counts]
                ).T
                ax.fill_between(counts, low, high, color=color, alpha=0.12, linewidth=0, zorder=1)
                ax.plot(
                    counts,
                    series[name],
                    color=color,
                    linewidth=2,
                    marker="o",
                    markersize=5,
                    markeredgecolor="white",
                    markeredgewidth=1.5,
                    zorder=3,
                )
            ax.set_xscale("log", base=2)
            ax.set(
                xlim=(counts[0] / 1.3, counts[-1] * 1.3),
                ylim=(0, 102),
                xticks=counts,
                xticklabels=[f"{k:,}" for k in counts],
                xlabel="Observed CpGs per patient (log scale)",
            )
            ax.minorticks_off()
            when = "the tested configuration" if round_number == tested else "trained longer"
            ax.set_title(f"After {round_number} rounds · {when}", loc="left", fontsize=12.5, pad=10)
        axes[0].set_ylabel("Balanced accuracy (%)")
        self.record["sparse_input"] = record
        handles = [
            Line2D([], [], color=CENTRALIZED, linewidth=2, marker="o", label="Centralized"),
            Line2D([], [], color=FEDERATED, linewidth=2, marker="o", label="Federated"),
            Line2D(
                [],
                [],
                color=FOCUS,
                linewidth=2,
                marker="o",
                markersize=4,
                label=f"Smallest site alone ({self.small})",
            ),
            Line2D([], [], color=CONTEXT, linewidth=1, label="Other sites alone"),
            Line2D([], [], color=FAINT, linewidth=1, label=self.chance_words()),
        ]
        fig.legend(
            handles=handles,
            ncol=5,
            loc="lower left",
            bbox_to_anchor=(0.06, 0.765),
            frameon=False,
            fontsize=10.5,
            handlelength=1.6,
            columnspacing=1.2,
        )
        self.heading(
            fig,
            "Accuracy falls steeply as fewer CpGs are observed",
            f"{self.validation_size()} validation patients · fewer observed CpGs mimic a "
            "shallower nanopore run · shaded: 95% intervals",
        )
        self.save(
            fig,
            "03_sparse_input",
            f"Synthetic data. Training drew between {floor:,.0f} and {len(self.panel):,} "
            "observed CpGs per patient; the gray region lies below that.",
        )

    # -- 04: what each model gets right -------------------------------------------------

    def class_recall(self):
        models = ["centralized", "federated"] + [f"local_{site}" for site in self.sites]
        rows = ["Centralized", "Federated"] + [self.site_label(site) for site in self.sites]
        recall = np.array(
            [[100 * (v or 0) for v in self.scores[m]["per_class_recall"]] for m in models]
        )
        test_counts = np.array(self.scores["centralized"]["confusion_matrix"]).sum(axis=1)
        balanced = [100 * self.scores[m]["balanced_accuracy"] for m in models]
        self.record["class_recall"] = {
            "models": models,
            "recall": recall.tolist(),
            "test_patients": test_counts.tolist(),
        }
        fig, ax = plt.subplots(figsize=(11, 6.6))
        fig.subplots_adjust(left=0.25, right=0.86, bottom=0.12, top=0.74)
        image = ax.imshow(recall, cmap="argonne_sequential", vmin=0, vmax=100, aspect="auto")
        colormap = plt.get_cmap("argonne_sequential")
        for (i, j), value in np.ndenumerate(recall):
            ax.text(
                j,
                i,
                f"{value:.0f}",
                ha="center",
                va="center",
                fontsize=11,
                color="white" if _dark(colormap(value / 100)) else INK,
            )
        for i, value in enumerate(balanced):
            ax.text(
                len(self.labels) - 0.35, i, f"{value:.1f}%", va="center", fontsize=11.5, color=INK
            )
        ax.text(
            len(self.labels) - 0.35,
            -0.75,
            "Balanced\naccuracy",
            va="bottom",
            fontsize=10,
            color=NOTE,
        )
        ax.set(
            xticks=range(len(self.labels)),
            yticks=range(len(rows)),
            yticklabels=rows,
            xticklabels=[
                f"Class {c}\nn = {n}" for c, n in zip(self.letters, test_counts, strict=True)
            ],
        )
        ax.tick_params(length=0)
        ax.xaxis.tick_top()
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.axhline(1.5, color="white", linewidth=4)
        self.heading(
            fig,
            "Which classes each model gets right",
            f"Recall (%) on the same {len(self.prep['splits']['test'])} held-out patients, "
            f"by true class · {self.highest:,} CpGs · {self.config['rounds']} rounds",
            y=0.965,
        )
        fig.colorbar(image, ax=ax, fraction=0.025, pad=0.13, label="Recall (%)")
        self.save(
            fig,
            "04_class_recall",
            "Synthetic data. Balanced accuracy is the mean of each row: "
            "every class counts equally.",
        )

    # -- 05: what each site holds -------------------------------------------------------

    def site_composition(self):
        counts = np.array(
            [
                [self.prep["site_class_counts"][site][label] for label in self.labels]
                for site in self.sites
            ]
        )
        share = 100 * counts / counts.sum(axis=1, keepdims=True)
        self.record["site_composition"] = {"sites": self.sites, "counts": counts.tolist()}
        fig, ax = plt.subplots(figsize=(11, 5.4))
        fig.subplots_adjust(left=0.2, right=0.9, bottom=0.12, top=0.68)
        image = ax.imshow(share, cmap="argonne_sequential", vmin=0, vmax=100, aspect="auto")
        colormap = plt.get_cmap("argonne_sequential")
        for (i, j), value in np.ndenumerate(counts):
            ax.text(
                j,
                i,
                f"{value}",
                ha="center",
                va="center",
                fontsize=12,
                color="white" if _dark(colormap(share[i, j] / 100)) else INK,
            )
        totals = counts.sum(axis=0)
        ax.set(
            xticks=range(len(self.labels)),
            yticks=range(len(self.sites)),
            xticklabels=[
                f"Class {c}\n{n} in total" for c, n in zip(self.letters, totals, strict=True)
            ],
            yticklabels=[
                f"{site} · {len(self.prep['sites'][site]):,} patients" for site in self.sites
            ],
        )
        ax.tick_params(length=0)
        ax.xaxis.tick_top()
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.colorbar(
            image, ax=ax, fraction=0.025, pad=0.03, label="Share of the site's patients (%)"
        )
        sparsest = int(counts.min())
        self.heading(
            fig,
            "Each site holds a skewed slice of the classes",
            f"Training patients per site and class · {counts.sum():,} in total · "
            f"the sparsest cell holds {sparsest} patient{'' if sparsest == 1 else 's'}",
            y=0.965,
        )
        self.save(
            fig,
            "05_site_composition",
            "Synthetic data. Dirichlet class skew and unequal capacity; every site is scored on "
            "the same test set.",
        )

    # -- appendix ----------------------------------------------------------------------

    def data_overview(self):
        """One class per panel, over every patient in gray: a single hue never has to
        separate six classes, which no palette can do under colour-vision deficiency."""
        from sklearn.decomposition import PCA

        filled = np.where(np.isnan(self.X), np.nanmean(self.X, axis=0), self.X).astype(np.float64)
        pca = PCA(n_components=2, svd_solver="randomized", random_state=0)
        embedding = pca.fit_transform(filled)
        explained = pca.explained_variance_ratio_
        self.record["pca_explained_variance"] = explained.tolist()
        columns = 3
        rows = -(-len(self.labels) // columns)
        fig, axes = plt.subplots(
            rows, columns, figsize=(12, 4 * rows + 1.6), sharex=True, sharey=True
        )
        fig.subplots_adjust(left=0.07, right=0.98, bottom=0.1, top=0.8, wspace=0.06, hspace=0.22)
        for index, ax in enumerate(np.ravel(axes)):
            if index >= len(self.labels):
                ax.set_visible(False)
                continue
            selected = self.y == self.labels[index]
            ax.scatter(
                *embedding[~selected].T, s=5, color=FAINT, edgecolors="none", rasterized=True
            )
            ax.scatter(
                *embedding[selected].T,
                s=9,
                color=CENTRALIZED,
                alpha=0.7,
                edgecolors="none",
                rasterized=True,
            )
            ax.set_title(
                f"Class {self.letters[index]} · {selected.sum()} patients",
                loc="left",
                fontsize=11.5,
            )
            ax.tick_params(labelsize=9.5)
        for ax in np.ravel(axes)[-columns:]:
            ax.set_xlabel(f"PC1 ({explained[0]:.1%} of variance)")
        for ax in np.ravel(axes)[::columns]:
            ax.set_ylabel(f"PC2 ({explained[1]:.1%})")
        self.heading(
            fig,
            "Appendix: training patients on the first two principal components",
            f"{len(self.y):,} training patients · beta values with training-mean imputation · "
            f"together the components explain {explained.sum():.1%} of variance",
            y=0.96,
        )
        self.save(
            fig,
            "A1_data_overview",
            "Synthetic data. Blue points away from their class's cluster are the generator's "
            "deliberate profile/label discordance.",
        )

    def cpg_methylation(self):
        class_means = np.array(
            [np.nanmean(self.X[self.y == label], axis=0) for label in self.labels]
        )
        overall = class_means.mean(axis=0)
        chosen = []
        for row in class_means:
            ranked = np.argsort(-(row - overall), kind="stable")
            chosen.extend([int(i) for i in ranked if i not in chosen][:10])
        controls = [int(i) for i in np.argsort(np.ptp(class_means, axis=0)) if i not in chosen][:10]
        chosen.extend(controls)
        self.record["display_cpg_ids"] = self.cpg_ids[self.panel[chosen]].tolist()
        groups = len(self.labels) + 1
        fig, ax = plt.subplots(figsize=(12, 6))
        fig.subplots_adjust(left=0.12, right=0.90, bottom=0.21, top=0.76)
        image = ax.imshow(
            class_means[:, chosen], aspect="auto", cmap="argonne_sequential", vmin=0, vmax=1
        )
        ax.set(
            yticks=range(len(self.labels)),
            yticklabels=[f"Class {letter}" for letter in self.letters],
            xticks=np.arange(groups) * 10 + 4.5,
            xticklabels=[f"{letter}-associated" for letter in self.letters] + ["Low contrast"],
        )
        ax.tick_params(axis="x", labelsize=10, length=0)
        ax.tick_params(axis="y", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for position in np.arange(1, groups) * 10 - 0.5:
            ax.axvline(position, color="white", linewidth=2)
        fig.colorbar(image, ax=ax, fraction=0.03, pad=0.025, label="Mean methylation beta")
        ax.set_xlabel("CpGs grouped by the class they contrast most", labelpad=13)
        self.heading(
            fig,
            "Appendix: class-mean methylation at the most contrasting CpGs",
            "Ten CpGs per class with the largest training-mean contrast, "
            "plus ten low-contrast controls",
            y=0.96,
        )
        self.save(
            fig,
            "A2_cpg_methylation",
            f"Synthetic data, training patients only. These {len(chosen)} CpGs were picked "
            "for contrast, so the blocks illustrate the generator; they are not a finding.",
        )


def _dark(rgba) -> bool:
    """Whether a fill is dark enough to need white text (relative luminance < 0.4)."""
    channels = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in rgba[:3]]
    return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2] < 0.4


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--run", type=Path, required=True, help="Completed test run of a synthetic scenario"
    )
    parser.add_argument("--curves", type=Path, help="Validation-only run with track_validation")
    parser.add_argument("--out", type=Path, required=True, help="New directory for figures")
    args = parser.parse_args(argv)
    try:
        record = build(args.run, args.out, args.curves)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"error: {exc}\n")
    print(f"Created {len(record['figures'])} figures in {args.out}")


if __name__ == "__main__":
    main()
