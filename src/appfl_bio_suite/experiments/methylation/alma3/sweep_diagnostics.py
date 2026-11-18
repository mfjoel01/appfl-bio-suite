"""Supplementary own-site comparisons and retention, using existing sweep scores."""

import csv
import json
from pathlib import Path

import numpy as np


def diagnostics(out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out)
    spec = json.loads((out / "specification.json").read_text())
    records = json.loads((out / "results.json").read_text())
    seeds, levels = spec["seeds"], spec["levels"]
    lookup = {(r["seed"], r["convergence"], r["arm"]): r for r in records}
    own = np.empty((len(seeds), len(levels), 5))
    fed = np.empty_like(own)
    rows = []
    for i, seed in enumerate(seeds):
        for j, level in enumerate(levels):
            for site in range(5):
                local = next(
                    m
                    for m in lookup[seed, level, f"local-S{site}"]["metrics"]
                    if m["site"] == site and not m["unshifted"]
                )
                global_ = next(
                    m
                    for m in lookup[seed, level, "fedavg"]["metrics"]
                    if m["site"] == site and not m["unshifted"]
                )
                own[i, j, site] = local["balanced_accuracy"]
                fed[i, j, site] = global_["balanced_accuracy"]
                rows.append(
                    dict(
                        seed=seed,
                        convergence=level,
                        site=site,
                        local_balanced_accuracy=own[i, j, site],
                        fedavg_balanced_accuracy=fed[i, j, site],
                        paired_gap=fed[i, j, site] - own[i, j, site],
                    )
                )
    with (out / "own-site-scores.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    rng = np.random.default_rng(51071)
    indices = rng.integers(len(seeds), size=(4000, len(seeds)))
    gap = fed - own
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for site in range(5):
        axes[0].plot(
            levels, 100 * own[:, :, site].mean(0), marker="o", label=f"Local S{site} at S{site}"
        )
    axes[0].plot(
        levels, 100 * fed.mean(2).mean(0), color="black", linewidth=3, label="FedAvg (site mean)"
    )
    axes[0].set(
        ylabel="Balanced accuracy (%)", ylim=(0, 104), title="Each local model at its own site"
    )
    axes[0].legend(fontsize=8, ncol=2, loc="lower right")
    summaries = []
    for name, values, color in (
        ("Mean own-site gain", gap.mean(2), "#d77915"),
        ("Smallest own-site gain", gap.min(2), "#8053a6"),
    ):
        lo, hi = np.quantile(values[indices].mean(1), [0.025, 0.975], axis=0)
        axes[1].plot(levels, 100 * values.mean(0), marker="o", color=color, label=name)
        if len(seeds) > 1:
            axes[1].fill_between(levels, 100 * lo, 100 * hi, color=color, alpha=0.15)
        for j, level in enumerate(levels):
            summaries.append(
                dict(
                    convergence=level,
                    contrast=name,
                    mean=float(values[:, j].mean()),
                    ci_low=float(lo[j]) if len(seeds) > 1 else None,
                    ci_high=float(hi[j]) if len(seeds) > 1 else None,
                )
            )
    (out / "own-site-gaps.json").write_text(json.dumps(summaries, indent=2) + "\n")
    axes[1].axhline(0, color="#555555", linewidth=1)
    axes[1].set(
        ylabel="FedAvg − own-site local (percentage points)", title="Paired within each site"
    )
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set(
            xlabel="Class-mix convergence (0: extreme; 1: identical)",
            xticks=levels,
            xticklabels=[f"{v:g}" for v in levels],
        )
        ax.grid(alpha=0.2)
    fig.suptitle(f"Own-site evaluation control · {len(seeds)} seeds · fixed balanced validation")
    fig.text(
        0.01,
        0.01,
        "Synthetic mechanism study; no clinical claims. Analysis performed using ALMA3.",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"own-site.{ext}", dpi=220)
    plt.close(fig)
    retention = {}
    for arm in ["zero-shot", "pooled", "fedavg", *[f"local-S{i}" for i in range(5)]]:
        retention[arm] = np.array(
            [
                [
                    np.mean(
                        [
                            m["balanced_accuracy"]
                            for m in lookup[seed, level, arm]["metrics"]
                            if m["unshifted"]
                        ]
                    )
                    for level in levels
                ]
                for seed in seeds
            ]
        )
    local_retention = np.stack([retention[f"local-S{i}"] for i in range(5)])
    retention["mean-local"] = local_retention.mean(0)
    retention["worst-local"] = local_retention.min(0)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    retention_rows = []
    for arm in ["zero-shot", "pooled", "fedavg", "mean-local", "worst-local"]:
        values = retention[arm]
        ax.plot(levels, 100 * values.mean(0), marker="o", label=arm)
        for i, seed in enumerate(seeds):
            for j, level in enumerate(levels):
                retention_rows.append(
                    dict(
                        seed=seed,
                        convergence=level,
                        method=arm,
                        unshifted_balanced_accuracy=float(values[i, j]),
                    )
                )
    with (out / "retention.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(retention_rows[0]))
        writer.writeheader()
        writer.writerows(retention_rows)
    ax.set(
        xlabel="Class-mix convergence",
        ylabel="Unshifted balanced accuracy (%)",
        xticks=levels,
        xticklabels=[f"{v:g}" for v in levels],
        ylim=(0, 104),
        title="Retention on the four studied classes",
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=9, ncol=2)
    fig.text(
        0.01,
        0.01,
        "Other 105 classes not tested. Synthetic validation. Analysis performed using ALMA3.",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out / "retention.png", dpi=220)
    fig.savefig(out / "retention.svg")
    plt.close(fig)
