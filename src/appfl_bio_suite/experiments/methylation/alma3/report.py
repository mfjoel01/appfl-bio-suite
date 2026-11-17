"""Meeting figures with paired seed uncertainty; never claims ALMA3 validation."""

import csv
import json
from collections import defaultdict

import numpy as np

from ..pipeline import write_json


def summarize(records):
    groups = defaultdict(list)
    for record in records:
        key = record["depth"], record["sigma"], record["arm"]
        rows = [m for m in record["metrics"] if not m["unshifted"]]
        max_coverage = max(m["coverage"] for m in rows)
        rows = [m for m in rows if m["coverage"] == max_coverage]
        # Unweighted across sites, deliberately not dominated by the largest hospital.
        value = float(np.mean([m["balanced_accuracy"] for m in rows]))
        groups[key].append((record["seed"], value))
        diagnostic = "support_restricted_balanced_accuracy"
        if record["arm"] == "zero-shot" and all(diagnostic in m for m in rows):
            control = record["depth"], record["sigma"], "zero-shot-known-support"
            groups[control].append((record["seed"], float(np.mean([m[diagnostic] for m in rows]))))
    output = []
    rng = np.random.default_rng(9041)
    for (depth, sigma, arm), pairs in sorted(groups.items()):
        values = np.array([p[1] for p in pairs])
        boot = rng.choice(values, (2000, len(values)), replace=True).mean(1)
        output.append(
            {
                "depth": depth,
                "sigma": sigma,
                "arm": arm,
                "n_seeds": len(values),
                "mean": float(values.mean()),
                "sd": float(values.std(ddof=1)) if len(values) > 1 else None,
                "ci_low": float(np.quantile(boot, 0.025)) if len(values) > 1 else None,
                "ci_high": float(np.quantile(boot, 0.975)) if len(values) > 1 else None,
            }
        )
    # Paired differences use the same data seed within each configuration.
    differences = []
    for (depth, sigma, arm), pairs in sorted(groups.items()):
        if arm == "zero-shot":
            continue
        for reference in ("zero-shot", "zero-shot-known-support", "pooled", "scratch-fedavg"):
            baseline = dict(groups.get((depth, sigma, reference), []))
            values = np.array([value - baseline[seed] for seed, value in pairs if seed in baseline])
            if len(values):
                boot = rng.choice(values, (2000, len(values)), replace=True).mean(1)
                differences.append(
                    {
                        "depth": depth,
                        "sigma": sigma,
                        "arm": arm,
                        "reference": reference,
                        "mean_difference": float(values.mean()),
                        "ci_low": float(np.quantile(boot, 0.025)) if len(values) > 1 else None,
                        "ci_high": float(np.quantile(boot, 0.975)) if len(values) > 1 else None,
                    }
                )
    return output, differences


def report(out, records, identity):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary, differences = summarize(records)
    write_json(out / "summary.json", summary)
    write_json(out / "paired-differences.json", differences)
    with (out / "summary.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    fixture = identity["kind"] == "tiny-stand-in"
    attribution = (
        "Tiny stand-in; no ALMA3 weights." if fixture else "Analysis performed using ALMA3."
    )
    (out / "REPORT.md").write_text(
        "# Synthetic methylation adaptation study\n\n" + attribution + "\n\n"
        "Mechanism illustration on synthetic data, not clinical validation. "
        "The manifest identifies validation-only versus the single frozen test evaluation.\n\n"
        "See summary.csv, paired-differences.json and the two PNG figures. "
        "Intervals bootstrap data seeds (only three by default), not hospitals or "
        "clinical populations. When present, zero-shot-known-support suppresses "
        "classes outside the configured synthetic support without updating weights. "
        "Nonprivate reference arms and synthetic evaluation "
        "mean this whole study is not a DP release.\n"
    )
    heading = (
        "Tiny stand-in — pipeline demonstration"
        if fixture
        else "ALMA3 member adaptation — synthetic data"
    )
    split = records[0]["metrics"][0]["split"]
    seed_count = len({r["seed"] for r in records})
    heading += f" ({split}; {seed_count} seed{'s' if seed_count != 1 else ''})"
    uncertainty_label = (
        "seed-bootstrap 95% intervals"
        if seed_count > 1
        else "single-seed pilot; no uncertainty estimate"
    )
    depths = list(dict.fromkeys(r["depth"] for r in records))
    fig, axes = plt.subplots(1, len(depths), figsize=(7 * len(depths), 5), squeeze=False)
    for ax, depth in zip(axes[0], depths, strict=True):
        subset = [
            r
            for r in summary
            if r["depth"] == depth and r["sigma"] == min(s["sigma"] for s in summary)
        ]
        x = np.arange(len(subset))
        ax.bar(x, [s["mean"] for s in subset], color="#247a86")
        for i, s in enumerate(subset):
            if s["ci_low"] is not None:
                ax.vlines(i, s["ci_low"], s["ci_high"], color="black", lw=2)
        ax.set_xticks(x, [s["arm"] for s in subset], rotation=40, ha="right")
        ax.set(title=f"Depth: {depth}", ylabel="Site-mean balanced accuracy", ylim=(0, 1))
    fig.suptitle(
        heading + f"\nNoise multiplier {min(r['sigma'] for r in records):g}; "
        "maximum requested coverage; " + uncertainty_label
    )
    fig.text(0.01, 0.005, attribution, fontsize=9)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(out / "comparison.png", dpi=180)
    fig.savefig(out / "comparison.svg")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for depth in depths:
        rows = [r for r in records if r["depth"] == depth and r["arm"] == "fedavg"]
        groups = defaultdict(list)
        calls = defaultdict(list)
        for r in rows:
            if r["privacy"]["epsilon"] is not None:
                score = np.mean(
                    [m["balanced_accuracy"] for m in r["metrics"] if not m["unshifted"]]
                )
                groups[r["privacy"]["epsilon"]].append(score)
                calls[r["privacy"]["epsilon"]].append(
                    np.mean(
                        [
                            m["hierarchy"][-1]["call_rate"]
                            for m in r["metrics"]
                            if not m["unshifted"]
                        ]
                    )
                )
        if groups:
            x = sorted(groups)
            axes[0].errorbar(
                x,
                [np.mean(groups[v]) for v in x],
                yerr=[np.std(groups[v]) for v in x],
                marker="o",
                label=f"{depth} accuracy",
            )
            axes[0].plot(
                x,
                [np.mean(calls[v]) for v in x],
                linestyle="--",
                marker="s",
                label=f"{depth} call rate",
            )
        for arm, marker in (
            ("zero-shot", "o"),
            ("fedavg", "s"),
            (f"local-S{identity.get('smallest_site', 3)}", "^"),
        ):
            values = [
                r for r in records if r["depth"] == depth and r["arm"] == arm and r["sigma"] == 0
            ]
            if not values:
                continue
            unshifted = [
                np.mean([m["balanced_accuracy"] for m in r["metrics"] if m["unshifted"]])
                for r in values
            ]
            # Smallest site's index is passed in identity by the driver.
            site = identity.get("smallest_site", 3)
            smallest = [
                np.mean(
                    [
                        m["balanced_accuracy"]
                        for m in r["metrics"]
                        if not m["unshifted"] and m["site"] == site
                    ]
                )
                for r in values
            ]
            label = f"{depth}/{arm}"
            axes[1].errorbar(label, np.mean(unshifted), yerr=np.std(unshifted), fmt=marker)
            axes[2].errorbar(label, np.mean(smallest), yerr=np.std(smallest), fmt=marker)
    axes[0].set(
        xlabel="Per-arm epsilon (training + calibration)",
        ylabel="Balanced accuracy / call rate",
        title="Privacy cost (accuracy mean ± seed SD)",
    )
    if axes[0].lines:
        axes[0].legend()
    for ax, title in zip(
        axes[1:], ("Unshifted holdout / forgetting", "Smallest-site holdout"), strict=True
    ):
        ax.set(title=title, ylim=(0, 1))
        ax.tick_params(axis="x", rotation=65)
    fig.suptitle(heading)
    fig.text(0.01, 0.005, attribution, fontsize=9)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(out / "privacy-forgetting-smallest.png", dpi=180)
    plt.close(fig)
    manifest = json.loads((out / "manifest.json").read_text())
    config = manifest["settings"]["config"]
    if config.get("site_class_proportions"):
        site_comparison(out, records, config, attribution)


def site_comparison(out, records, config, attribution):
    """Show declared synthetic class mixes and each model's cross-site performance."""
    import matplotlib.pyplot as plt

    proportions = np.asarray(config["site_class_proportions"])
    sites = list(range(len(proportions)))
    arms = ["zero-shot", "zero-shot-known-support"] + [
        a for a in dict.fromkeys(r["arm"] for r in records) if a != "zero-shot"
    ]
    # One panel per depth/noise pair so settings are never silently pooled.
    for depth, sigma in dict.fromkeys((r["depth"], r["sigma"]) for r in records):
        matrix = np.full((len(arms), len(sites)), np.nan)
        rows = []
        for i, arm in enumerate(arms):
            for site in sites:
                control = arm == "zero-shot-known-support"
                key = "support_restricted_balanced_accuracy" if control else "balanced_accuracy"
                values = []
                for r in records:
                    if (r["depth"], r["sigma"], r["arm"]) != (
                        depth,
                        sigma,
                        "zero-shot" if control else arm,
                    ):
                        continue
                    for m in r["metrics"]:
                        if (
                            m["site"] == site
                            and not m["unshifted"]
                            and m["coverage"] == max(config["coverage"])
                            and key in m
                        ):
                            values.append(m[key])
                            rows.append(
                                dict(
                                    arm=arm,
                                    site=site,
                                    seed=r["seed"],
                                    balanced_accuracy=m[key],
                                    accuracy=None if control else m["accuracy"],
                                )
                            )
                if values:
                    matrix[i, site] = np.mean(values)
        stem = f"sites-{depth}-noise-{sigma}"
        write_json(out / f"{stem}.json", rows)
        with (out / "REPORT.md").open("a") as stream:
            stream.write(
                f"\n## Local models evaluated at their own site ({depth}, noise {sigma:g})\n\n"
                "Mean across seeds on shifted cohorts at maximum coverage. "
                "Compare this with the cross-site balanced-accuracy heatmap.\n\n"
                "| Site | Accuracy | Balanced accuracy |\n|---|---:|---:|\n"
            )
            for site in sites:
                own = [r for r in rows if r["site"] == site and r["arm"] == f"local-S{site}"]
                if own:
                    stream.write(
                        f"| S{site} | {np.mean([r['accuracy'] for r in own]):.4f} | "
                        f"{np.mean([r['balanced_accuracy'] for r in own]):.4f} |\n"
                    )
        fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.5]})
        bottom = np.zeros(len(sites))
        for j, label in enumerate(config["class_indices"]):
            axes[0].bar(sites, proportions[:, j] * 100, bottom=bottom, label=f"Class {label}")
            bottom += proportions[:, j] * 100
        axes[0].set(
            xticks=sites,
            xticklabels=[f"S{i}" for i in sites],
            ylabel="Patients (%)",
            title="Configured class proportions",
            ylim=(0, 100),
        )
        axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2)
        plot = axes[1].imshow(matrix, vmin=0, vmax=1, cmap="YlGnBu", aspect="auto")
        for i in range(len(arms)):
            for j in sites:
                if np.isfinite(matrix[i, j]):
                    axes[1].text(
                        j,
                        i,
                        f"{matrix[i, j]:.2f}",
                        ha="center",
                        va="center",
                        color="white" if matrix[i, j] > 0.6 else "black",
                    )
        axes[1].set(
            xticks=sites,
            xticklabels=[f"S{i}" for i in sites],
            yticks=list(range(len(arms))),
            yticklabels=arms,
            xlabel="Evaluation site",
            title="Shifted balanced accuracy (mean across seeds)",
        )
        fig.colorbar(plot, ax=axes[1], fraction=0.04, pad=0.03)
        fig.suptitle(
            f"Synthetic site heterogeneity — {depth}; noise {sigma:g}; "
            f"{len(set(r['seed'] for r in records))} seeds; "
            f"{records[0]['metrics'][0]['split']}"
        )
        fig.text(0.01, 0.005, attribution, fontsize=9)
        fig.tight_layout(rect=(0, 0.04, 1, 0.96))
        fig.savefig(out / f"{stem}.png", dpi=180)
        fig.savefig(out / f"{stem}.svg")
        plt.close(fig)
