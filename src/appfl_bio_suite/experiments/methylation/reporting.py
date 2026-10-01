"""Export self-contained figures and an evidence-based run summary."""

from __future__ import annotations

import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from appfl_bio_suite.core import plot_style as brand  # noqa: E402

COLORS = {
    "local": brand.COLORS["anl-gray-500"],
    "federated": brand.COLORS["anl-cyan"],
    "centralized": brand.COLORS["anl-navy"],
}


def report(out, results, preprocessing):
    brand.apply_style()
    counts = sorted({row["requested_cpgs"] for row in results})
    full = {r["model"]: r for r in results if r["requested_cpgs"] == counts[-1]}
    fig, ax = plt.subplots(figsize=(7, 4))
    for i, name in enumerate(("local", "federated", "centralized")):
        scores = [
            full[f"local_S{s + 1}" if name == "local" else name]["balanced_accuracy"]
            for s in range(4)
        ]
        ax.bar(np.arange(4) + (i - 1) * 0.25, scores, width=0.25, label=name, color=COLORS[name])
    ax.set(
        xticks=range(4),
        xticklabels=[f"S{i + 1}" for i in range(4)],
        ylim=(0, 1),
        ylabel="Balanced accuracy on shared test set",
        title="Training-site comparison",
    )
    fig.suptitle(f"{preprocessing['input_kind']} · serial simulation")
    fig.legend(loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    fig.savefig(out / "site_comparison.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4))
    for name in ("centralized", "federated"):
        rows = sorted((r for r in results if r["model"] == name), key=lambda r: r["requested_cpgs"])
        ax.plot(
            [r["requested_cpgs"] for r in rows],
            [r["balanced_accuracy"] for r in rows],
            marker="o",
            label=name,
            color=COLORS[name],
        )
    ax.set(xlabel="Requested observed CpGs", ylabel="Balanced accuracy", ylim=(0, 1))
    fig.suptitle(f"{preprocessing['input_kind']} · serial simulation")
    fig.legend(loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    fig.savefig(out / "coverage_curve.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 5))
    matrix = np.array(full["federated"]["confusion_matrix"])
    ax.imshow(matrix, cmap="argonne_sequential")
    for (i, j), count in np.ndenumerate(matrix):
        # Light text on the dark half of the scale, dark text on the light half.
        color = "white" if count > matrix.max() / 2 else brand.COLORS["anl-ink"]
        ax.text(j, i, str(count), ha="center", va="center", color=color)
    labels = preprocessing["classes"]
    ax.set(
        xticks=range(len(labels)),
        yticks=range(len(labels)),
        xticklabels=labels,
        yticklabels=labels,
        xlabel="Predicted",
        ylabel="True",
        title="Federated model",
    )
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.suptitle(f"{preprocessing['input_kind']} · serial simulation")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out / "confusion_matrix.png", dpi=160)
    plt.close(fig)
    with (out / "site_class_counts.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["site", *labels])
        for name, values in preprocessing["site_class_counts"].items():
            writer.writerow([name, *[values[label] for label in labels]])
    gap = full["federated"]["balanced_accuracy"] - full["centralized"]["balanced_accuracy"]
    gains = {
        f"S{i + 1}": full["federated"]["balanced_accuracy"]
        - full[f"local_S{i + 1}"]["balanced_accuracy"]
        for i in range(4)
    }
    lines = [
        "# Methylation experiment results",
        "",
        f"Single-node research simulation. Input: {preprocessing['input_kind']}.",
        "",
        f"Highest evaluated coverage: {counts[-1]} CpGs.",
        f"Federated minus centralized balanced accuracy: {gap:+.4f}.",
        "",
        "| Training site | Federated minus local balanced accuracy |",
        "| --- | --- |",
    ]
    lines.extend(f"| {site} | {gain:+.4f} |" for site, gain in gains.items())
    lines.extend(
        [
            "",
            "All arms use the same test samples. The repeated federated and pooled",
            "bars are the same shared-test score, not site-specific test cohorts.",
            "",
            "See metrics.json for calibration, class recall, macro-F1 and confidence results.",
            "A null confident_accuracy means no predictions reached 0.90 confidence.",
            "",
            "## From PoC to a network",
            "",
            "Validate an authorized real cohort and the probe/label schema; benchmark real",
            "nanopore inputs; then add partner execution, signed containers and metadata.",
            "Privacy mechanisms and clinical validation remain separate work.",
        ]
    )
    (out / "REPORT.md").write_text("\n".join(lines) + "\n")
