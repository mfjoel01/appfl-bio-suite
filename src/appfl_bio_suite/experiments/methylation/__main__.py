"""Generate a development fixture: python -m ...methylation --out DIR."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .dataset import synthetic_data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=120)
    parser.add_argument("--cpgs", type=int, default=256)
    parser.add_argument("--classes", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--run", action="store_true", help="Run the complete E1–E4 comparison")
    parser.add_argument("--config", default="ci-tiny", help="Scenario name or YAML path for --run")
    parser.add_argument("--input", type=Path, help="Supplied NPZ cohort for --run")
    parser.add_argument(
        "--validation-only",
        action="store_true",
        help="Evaluate validation patients only; never predict test patients",
    )
    args = parser.parse_args()
    if args.run:
        from .pipeline import load_config, run_experiment

        print(
            run_experiment(
                args.out, load_config(args.config), args.input, validation_only=args.validation_only
            )
        )
        return
    if args.input is not None or args.config != "ci-tiny" or args.validation_only:
        parser.error("--input, --config and --validation-only require --run")
    parameters = {
        key: value
        for key, value in vars(args).items()
        if key in ("samples", "cpgs", "classes", "seed")
    }
    data = synthetic_data(**parameters)
    args.out.mkdir(parents=True, exist_ok=False)
    dataset_path = args.out / "cohort.npz"
    data.save(dataset_path)
    manifest = {
        "schema_version": 1,
        "experiment": "methylation",
        "stage": "synthetic-fixture-only",
        "parameters": parameters,
        "numpy_version": np.__version__,
        "shape": list(data.X.shape),
        "class_counts": dict(zip(*np.unique(data.y, return_counts=True), strict=True)),
        "dataset_sha256": hashlib.sha256(dataset_path.read_bytes()).hexdigest(),
        "limitations": "No patient split, site assignment, training or biological validation.",
    }
    manifest["class_counts"] = {k: int(v) for k, v in manifest["class_counts"].items()}
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Synthetic fixture saved to {args.out}")


if __name__ == "__main__":
    main()
