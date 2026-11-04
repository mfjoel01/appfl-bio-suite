"""The run the GWAS golden checksums describe, defined once.

The parity tests use it in two ways. In-process, it is the pipeline run that every
committed checksum is compared against. As a script, it is the same run in a fresh
interpreter started with CPU-feature overrides in its environment::

    python golden_pipeline.py WORK_DIR

That second use has to be a separate process because OpenBLAS and numpy read those
overrides once, when they load. The script prints one line of JSON: the checksum of every
output the golden file records, and what numpy and BLAS actually ran on, so the caller can
tell whether the override took effect.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

GOLDEN_PATH = Path(__file__).with_name("gwas_legacy_golden.json")


def load_golden() -> dict:
    return json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))


def run(root: Path, golden: dict) -> Path:
    """Build the fixture under ``root`` and run simulation steps 1 and 2 on it."""
    from make_fixture import build

    from appfl_bio_suite.experiments.gwas.simulation.bundler import bundle_sites
    from appfl_bio_suite.experiments.gwas.simulation.phenotypes import simulate_phenotypes

    fixture = golden["fixture"]
    info = build(
        root,
        n_samples=fixture["n_samples"],
        n_variants=fixture["n_variants"],
        seed=fixture["seed"],
    )
    simulate_phenotypes(
        plink_prefix=info["plink_prefix"],
        pgs_t2d_path=info["pgs_t2d"],
        pgs_bmi_path=info["pgs_bmi"],
        out_dir=root / "out",
    )
    bundle_sites(
        plink_prefix=info["plink_prefix"],
        pooled_dir=root / "out",
        sites_dir=root / "sites",
        site_sizes=golden["site_sizes"],
        split_seed=42,
    )
    return root


def resolve(root: Path, key: str) -> Path:
    """Where the output a golden key names lives under a run's ``root``."""
    if key.startswith("pooled/"):
        return root / "out" / key.split("/", 1)[1]
    site, name = key.split("/", 1)
    return root / "sites" / site / "data" / name


def checksums(root: Path, keys) -> dict[str, str | None]:
    """Checksum each named output. None for one the run did not produce."""
    from appfl_bio_suite.core.simulation import file_checksum

    out: dict[str, str | None] = {}
    for key in keys:
        path = resolve(root, key)
        out[key] = file_checksum(path) if path.is_file() else None
    return out


def cpu_state() -> dict:
    """The CPU features numpy dispatches to and the kernel BLAS selected, in this process.

    Both are chosen at load time from the CPU and from environment overrides. Recording
    them is what lets a caller confirm that an override changed anything, rather than
    passing a comparison of two identical configurations.
    """
    import numpy as np  # noqa: F401  (BLAS is loaded, and so visible, only after this)

    try:
        from numpy._core._multiarray_umath import __cpu_dispatch__, __cpu_features__
    except ImportError:  # numpy < 2
        from numpy.core._multiarray_umath import __cpu_dispatch__, __cpu_features__

    state: dict = {
        "numpy_dispatch": sorted(f for f in __cpu_dispatch__ if __cpu_features__.get(f)),
        "blas": None,
    }
    try:
        from threadpoolctl import threadpool_info
    except ImportError:
        return state
    # A set: numpy and scipy each bundle an OpenBLAS, and whether scipy's is loaded yet
    # says nothing about which kernel runs.
    state["blas"] = sorted(
        {
            f"{lib['internal_api']}:{lib.get('architecture')}"
            for lib in threadpool_info()
            if lib["user_api"] == "blas"
        }
    )
    return state


def main(argv: list[str]) -> int:
    root = Path(argv[1])
    golden = load_golden()
    run(root, golden)
    report = {"outputs": checksums(root, golden["outputs"]), "cpu": cpu_state()}
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
