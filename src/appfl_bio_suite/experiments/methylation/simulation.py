"""Suite simulate command: create an unsplit synthetic cohort and provenance."""

import importlib.metadata
import json
import platform
import subprocess
import sys
from pathlib import Path

from appfl_bio_suite.core.simulation import file_checksum

from .dataset import synthetic_data
from .pipeline import load_config, write_json


def cli_entry(*, scenario, out_dir, verify_manifest, list_scenarios):
    """Back ``appfl-bio-suite simulate methylation``. Returns an exit code."""
    if list_scenarios:
        print("\n".join(p.stem for p in sorted((Path(__file__).parent / "configs").glob("*.yaml"))))
        return 0
    try:
        if verify_manifest:
            return _verify(Path(verify_manifest))
        if out_dir is None:
            raise ValueError("--out is required")
        return _simulate(scenario or "ci-tiny", out_dir)
    except (ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


def _verify(path):
    expected = json.loads(path.read_text()).get("dataset_sha256")
    if expected is None:
        raise ValueError(f"{path} is not a methylation cohort manifest")
    if expected != file_checksum(path.parent / "cohort.npz"):
        print(f"cohort checksum mismatch: {path.parent / 'cohort.npz'}", file=sys.stderr)
        return 1
    print("Cohort checksum verified")
    return 0


def _simulate(scenario, out_dir):
    config = load_config(scenario)
    data = synthetic_data(
        samples=config.samples,
        cpgs=config.cpgs,
        classes=config.classes,
        seed=config.seed,
        signal_fraction=config.signal_fraction,
        class_weights=config.class_weights,
        profile_noise=config.profile_noise,
        missing_fraction=config.missing_fraction,
    )
    out_dir.mkdir(parents=True, exist_ok=False)
    data.save(out_dir / "cohort.npz")
    write_json(
        out_dir / "manifest.json",
        {
            "experiment": "methylation",
            "stage": "synthetic-cohort",
            "scenario": scenario,
            "seed": config.seed,
            "samples": config.samples,
            "cpgs": config.cpgs,
            "classes": config.classes,
            "signal_fraction": config.signal_fraction,
            "class_weights": config.class_weights,
            "profile_noise": config.profile_noise,
            "missing_fraction": config.missing_fraction,
            "dataset_sha256": file_checksum(out_dir / "cohort.npz"),
            "python": platform.python_version(),
            "numpy": importlib.metadata.version("numpy"),
            "suite_commit": subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=Path(__file__).parent,
                capture_output=True,
                text=True,
                check=False,
            ).stdout.strip(),
            "generator_sha256": file_checksum(Path(__file__).parent / "dataset.py"),
        },
    )
    print(f"Cohort saved to {out_dir}; patient holdouts and sites are assigned by run")
    return 0
