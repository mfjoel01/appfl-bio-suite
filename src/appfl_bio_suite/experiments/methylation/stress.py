"""Paired synthetic sensitivity runs with and without a minimum class count per site."""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from .pipeline import load_config, run_experiment, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 43, 44])
    args = parser.parse_args()
    if len(set(args.seeds)) != len(args.seeds):
        parser.error("seeds must be unique")
    base = load_config("stress")
    if base.batch_shift != 0:
        raise ValueError("paired class-support study requires batch_shift=0")
    for seed in args.seeds:
        replace(base, seed=seed).validate()
    args.out.mkdir(parents=True, exist_ok=False)
    write_json(args.out / "status.json", {"complete": False, "completed_runs": []})
    rows, completed = [], []
    for seed in args.seeds:
        for condition, minimum in (("unconstrained", 0), ("all_classes_present", 2)):
            directory = args.out / f"seed-{seed}-{condition}"
            config = replace(base, seed=seed, min_site_class_patients=minimum)
            run_experiment(directory, config)
            results = json.loads((directory / "metrics.json").read_text())
            rows.extend({"seed": seed, "condition": condition, **result} for result in results)
            completed.append(directory.name)
            write_json(args.out / "status.json", {"complete": False, "completed_runs": completed})
    for seed in args.seeds:
        pooled = [
            [
                r
                for r in rows
                if r["seed"] == seed and r["condition"] == condition and r["model"] == "centralized"
            ]
            for condition in ("unconstrained", "all_classes_present")
        ]
        comparable = [
            {
                r["requested_cpgs"]: {k: v for k, v in r.items() if k != "condition"}
                for r in condition_rows
            }
            for condition_rows in pooled
        ]
        if comparable[0] != comparable[1]:
            raise RuntimeError("paired conditions changed the centralized reference results")
    summary = []
    for condition in ("unconstrained", "all_classes_present"):
        for coverage in base.coverage_counts:
            group = [
                r for r in rows if r["condition"] == condition and r["requested_cpgs"] == coverage
            ]
            entry = {"condition": condition, "cpgs": coverage, "seeds": args.seeds}
            for model in ("centralized", "federated", "local_S4"):
                values = [r["balanced_accuracy"] for r in group if r["model"] == model]
                entry[model] = {"mean": float(np.mean(values)), "std": float(np.std(values))}
            summary.append(entry)
    write_json(args.out / "summary.json", summary)
    lines = [
        "# Synthetic stress comparison",
        "",
        "Only 2% of CpGs carry a class-specific signal. This is a synthetic sensitivity",
        "test, not a biological benchmark. Paired conditions keep data, holdouts and site",
        "capacities fixed; reserving class examples changes site label composition.",
        "Batch shifts are disabled. Pooled metrics must match exactly within each seed.",
        "",
        f"Seeds: {', '.join(map(str, args.seeds))}. Values are mean ± population SD",
        "across seeds, not confidence intervals.",
        "",
        "| Site constraint | CpGs | Pooled | Federated | Small-site local |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in summary:
        values = [
            f"{row[name]['mean']:.4f} ± {row[name]['std']:.4f}"
            for name in ("centralized", "federated", "local_S4")
        ]
        lines.append(f"| {row['condition']} | {row['cpgs']} | " + " | ".join(values) + " |")
    (args.out / "REPORT.md").write_text("\n".join(lines) + "\n")
    write_json(args.out / "status.json", {"complete": True, "completed_runs": completed})
    print(f"Stress comparison saved to {args.out}")


if __name__ == "__main__":
    main()
