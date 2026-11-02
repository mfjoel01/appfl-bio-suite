"""Adapter for the suite's serial experiment entry point."""

from dataclasses import asdict
from pathlib import Path

from .pipeline import load_config, run_experiment, write_json


def run(*, variant, driver, out_dir, dry_run, watch, data_root, federation_path):
    if driver != "serial":
        raise ValueError("the methylation experiment currently requires --driver serial")
    if watch or federation_path:
        raise ValueError("the methylation simulation does not accept --watch or --federation")
    config = load_config("ci-tiny" if variant in ("default", "loopback") else variant)
    out = out_dir or Path("local/output/methylation")
    input_path = None
    if data_root:
        input_path = data_root / "cohort.npz" if data_root.is_dir() else data_root
        if not input_path.is_file():
            raise ValueError(f"cohort not found: {input_path}")
    if dry_run:
        out.mkdir(parents=True, exist_ok=False)
        write_json(out / "config.json", asdict(config))
        return {"dry_run": True, "output": str(out)}
    return run_experiment(out, config, input_path)
