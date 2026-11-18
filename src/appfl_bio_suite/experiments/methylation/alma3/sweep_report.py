"""Presentation figures with paired uncertainty and an explicit best-local control."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def report_sweep(out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out)
    spec = json.loads((out / "specification.json").read_text())
    records = json.loads((out / "results.json").read_text())
    seeds, levels = spec["seeds"], spec["levels"]
    scores = {}
    detail = []
    for r in records:
        rows = [m for m in r["metrics"] if not m["unshifted"]]
        value = np.mean([m["balanced_accuracy"] for m in rows])
        scores[r["seed"], r["convergence"], r["arm"]] = value
        detail.append(
            dict(
                seed=r["seed"],
                convergence=r["convergence"],
                arm=r["arm"],
                balanced_accuracy=float(value),
            )
        )
        if r["arm"] == "zero-shot":
            scores[r["seed"], r["convergence"], "known-support"] = np.mean(
                [m["support_restricted_balanced_accuracy"] for m in rows]
            )
    methods = ["zero-shot", "known-support", "pooled", "fedavg"] + [f"local-S{i}" for i in range(5)]
    arrays = {
        arm: np.array([[scores[s, level, arm] for level in levels] for s in seeds])
        for arm in methods
    }
    locals_ = np.stack([arrays[f"local-S{i}"] for i in range(5)])
    arrays["mean-local"] = locals_.mean(0)
    # Best is chosen separately per seed: a deliberately strong oracle comparator.
    arrays["best-local"] = locals_.max(0)
    rng = np.random.default_rng(91041)
    indices = rng.integers(len(seeds), size=(4000, len(seeds)))

    def interval(values):
        return np.quantile(values[indices].mean(1), [0.025, 0.975], axis=0)

    summary, gaps = [], []
    for arm, values in arrays.items():
        lo, hi = interval(values)
        for j, level in enumerate(levels):
            summary.append(
                dict(
                    convergence=level,
                    method=arm,
                    mean=float(values[:, j].mean()),
                    ci_low=float(lo[j]) if len(seeds) > 1 else None,
                    ci_high=float(hi[j]) if len(seeds) > 1 else None,
                    n_seeds=len(seeds),
                )
            )
    for baseline in ["mean-local", "best-local", *[f"local-S{i}" for i in range(5)]]:
        values = arrays["fedavg"] - arrays[baseline]
        lo, hi = interval(values)
        for j, level in enumerate(levels):
            gaps.append(
                dict(
                    convergence=level,
                    comparator=baseline,
                    fedavg_minus_local=float(values[:, j].mean()),
                    ci_low=float(lo[j]) if len(seeds) > 1 else None,
                    ci_high=float(hi[j]) if len(seeds) > 1 else None,
                )
            )
    for name, rows in [("scores", detail), ("summary", summary), ("paired-gaps", gaps)]:
        (out / f"{name}.json").write_text(json.dumps(rows, indent=2) + "\n")
        with (out / f"{name}.csv").open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    colors = {
        "fedavg": "#007f86",
        "pooled": "#333333",
        "mean-local": "#d77915",
        "best-local": "#8053a6",
        "known-support": "#888888",
    }
    labels = {
        "fedavg": "Federated (FedAvg)",
        "pooled": "Pooled",
        "mean-local": "Mean local",
        "best-local": "Best local (per-seed oracle)",
        "known-support": "Zero-shot, known classes",
    }
    if np.array_equal(arrays["fedavg"], arrays["pooled"]):
        colors.pop("pooled")
        labels["fedavg"] = "FedAvg / pooled (same scores)"
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for arm in colors:
        values = arrays[arm]
        axes[0].plot(
            levels,
            100 * values.mean(0),
            marker="o",
            color=colors[arm],
            linestyle="--" if arm in {"pooled", "known-support"} else "-",
            label=labels[arm],
        )
        if arm in {"fedavg", "mean-local", "best-local"} and len(seeds) > 1:
            lo, hi = interval(values)
            axes[0].fill_between(levels, 100 * lo, 100 * hi, color=colors[arm], alpha=0.12)
    axes[0].set(ylabel="Balanced accuracy (%)", ylim=(0, 104), title="Do local models catch up?")
    axes[0].legend(fontsize=8, loc="lower right")
    for arm in ("mean-local", "best-local"):
        values = arrays["fedavg"] - arrays[arm]
        lo, hi = interval(values)
        axes[1].plot(levels, 100 * values.mean(0), marker="o", color=colors[arm], label=labels[arm])
        if len(seeds) > 1:
            axes[1].fill_between(levels, 100 * lo, 100 * hi, color=colors[arm], alpha=0.15)
    axes[1].axhline(0, color="#555555", linewidth=1)
    axes[1].set(ylabel="FedAvg − local (percentage points)", title="Paired federated advantage")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set(
            xlabel="Class-mix convergence: 0 = extreme skew; 1 = identical ratios",
            xticks=levels,
            xticklabels=[f"{v:g}" for v in levels],
            xlim=(min(levels) - 0.02, max(levels) + 0.02),
        )
        ax.grid(alpha=0.2)
    fig.suptitle(
        f"Synthetic five-site class-mix sweep · {len(seeds)} data seeds\n"
        "Fixed pooled class totals, training budget and balanced evaluation cohorts"
    )
    fig.text(
        0.01,
        0.008,
        (
            "Validation only; seed-bootstrap 95% bands. "
            if len(seeds) > 1
            else "One-seed development check; no uncertainty estimate. "
        )
        + "Synthetic mechanism study, not clinical validation. Analysis performed using ALMA3.",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"sweep.{ext}", dpi=220)
    if len(levels) > 2 and 0.25 in levels:
        for ax in axes:
            ax.set(
                xlim=(-0.005, 0.255),
                xticks=[v for v in levels if v <= 0.25],
                xlabel="Class-mix convergence (early region; identical endpoint = 1)",
            )
        fig.suptitle(
            f"Early convergence detail · {len(seeds)} data seeds\n"
            "Same prespecified sweep; 0–25% of movement toward shared class ratios"
        )
        for ext in ("png", "svg", "pdf"):
            fig.savefig(out / f"sweep-early.{ext}", dpi=220)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for site in range(5):
        arm = f"local-S{site}"
        ax.plot(
            levels,
            100 * arrays[arm].mean(0),
            marker=["o", "s", "^", "D", "v"][site],
            label=f"Local S{site}",
        )
    ax.plot(levels, 100 * arrays["fedavg"].mean(0), color="black", linewidth=3, label="FedAvg")
    ax.set(
        xlabel="Class-mix convergence (0: extreme; 1: identical)",
        ylabel="Balanced accuracy (%)",
        xticks=levels,
        xticklabels=[f"{v:g}" for v in levels],
        ylim=(0, 104),
        title="Every local model — no local site omitted",
    )
    ax.grid(alpha=0.2)
    ax.legend(ncol=3, fontsize=9)
    fig.text(
        0.01,
        0.01,
        "Mean across data seeds; synthetic validation only. Analysis performed using ALMA3.",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out / "all-locals.png", dpi=220)
    fig.savefig(out / "all-locals.svg")
    plt.close(fig)
    counts = {r["convergence"]: r["class_counts"] for r in records}
    columns = min(4, len(levels))
    row_count = (len(levels) + columns - 1) // columns
    fig, axes = plt.subplots(
        row_count, columns, figsize=(3.5 * columns, 3.4 * row_count), squeeze=False
    )
    for ax in axes.flat:
        ax.set_visible(False)
    for ax, level in zip(list(axes.flat)[: len(levels)], levels, strict=True):
        ax.set_visible(True)
        matrix = np.array(counts[level])
        bottom = np.zeros(5)
        for j, label in enumerate(spec["class_indices"]):
            ax.bar(range(5), matrix[:, j], bottom=bottom, label=f"Class {label}")
            bottom += matrix[:, j]
        ax.set(
            xticks=range(5),
            xticklabels=[f"S{i}" for i in range(5)],
            title=f"Convergence {level:g}",
            ylim=(0, spec["patients_per_site"]),
        )
    axes[0, 0].set_ylabel("Training patients per site")
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, names, loc="lower center", ncol=len(names))
    fig.suptitle("Class mix changes; site sizes and pooled totals stay fixed")
    fig.tight_layout(rect=(0, 0.1, 1, 0.9))
    fig.savefig(out / "site-mixes.png", dpi=220)
    plt.close(fig)
    lines = [
        "# Five-site synthetic class-mix sweep",
        "",
        "Analysis performed using ALMA3.",
        "",
        "Validation-only mechanism study. This is not clinical evidence or a DP release.",
        "",
        "| Convergence | FedAvg | Mean local | Best local | FedAvg minus best (pp) |",
        "|---|---:|---:|---:|---:|",
    ]
    for j, level in enumerate(levels):
        fed, mean, best = [
            100 * arrays[a][:, j].mean() for a in ("fedavg", "mean-local", "best-local")
        ]
        lines.append(f"| {level:g} | {fed:.2f}% | {mean:.2f}% | {best:.2f}% | {fed - best:.2f} |")
    lines += [
        "",
        "Best local is the maximum per seed over all five local models, evaluated "
        "on the same balanced cohorts at all sites. It is an oracle comparator, "
        "not a selected deployable model.",
        "",
        "Interpolation preserves all site counts and pooled class totals. Patient "
        "reservoirs and evaluation cohorts are paired across sweep points. Site "
        "measurement shifts remain fixed; only class ratios converge.",
        "",
        "Federated/pooled training uses 500 patients; each local model uses 100. "
        "Learning rate, clipping, number of rounds and local steps are fixed. No "
        "early stopping or per-point retuning. Zero-shot known-class support is "
        "an explicit no-training control.",
        "",
        "Intervals resample data seeds jointly across convergence points and "
        "compare models within seed. With few seeds these intervals are "
        "descriptive, not a clinical-population uncertainty estimate. No "
        "calibrated decision policy is fitted in this sweep.",
        "",
        "See sweep.png/svg/pdf, all-locals.png, site-mixes.png, summary.csv and paired-gaps.csv.",
    ]
    lines += [
        "",
        "Local training uses ordinary unweighted examples; class-balanced local "
        "resampling or reweighting is not evaluated. This does not establish superiority "
        "over all possible local training strategies.",
        "",
        "Own-site comparisons and four-class unshifted retention are supplementary "
        "views of the same recorded predictions, not additional training or test data.",
    ]
    (out / "REPORT.md").write_text("\n".join(lines) + "\n")
    from .sweep_diagnostics import diagnostics

    diagnostics(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    report_sweep(parser.parse_args().directory)
