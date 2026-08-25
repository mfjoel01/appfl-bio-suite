"""Run the PRE-MIGRATION simulation scripts, unmodified, to capture a golden reference.

WHY THIS EXISTS
---------------
The GWAS simulation is scientifically load-bearing: it reproduces the design of a
published paper, and its outputs are the data behind figures that have been reported. Any
refactor of it therefore has exactly one acceptance criterion -- **the output must be
byte-identical for a fixed seed** -- and that criterion is worthless unless it is checked
against the original code rather than against a reimplementation of what the original was
believed to do.

So this harness runs the legacy scripts *as they are*, with no edits, and checksums
everything they produce. The port is then held against that.

This must be built and green BEFORE any numerical code is touched. A characterization
test written after a refactor characterizes the refactor.

HOW THE LEGACY SCRIPTS ARE INVOKED
----------------------------------
They cannot be imported: they are top-level scripts that do their work at module scope,
and they locate their inputs by walking up from ``__file__``::

    BASE = Path(__file__).resolve().parent.parent        # local/
    IN_DIR = BASE / "data_sim" / "input"
    sys.path.insert(0, str(BASE.parent / "examples" / "resources" / "configs" / "gwas"))

So the harness builds a scratch tree with that exact shape, copies the scripts into it,
drops the fixture into the input directory, and runs them as subprocesses. Nothing in the
real source trees is read for writing or modified -- they are live, and one of them backs
a partner experiment currently in flight.

WHAT IT PINS DOWN
-----------------
Two determinism details that a well-meaning modernization would silently break, and which
the checksums exist to catch:

* ``data_bundler.py`` uses the LEGACY global RNG -- ``np.random.seed(42)`` followed by
  ``np.random.shuffle(indices)``. Switching to ``np.random.default_rng(42)`` is the
  obvious modernization and produces a completely different sample split, with no error
  and no obvious symptom.
* ``standardise()`` uses ``ndarray.std()``, which is the population standard deviation
  (ddof=0). ``pandas.Series.std()`` defaults to ddof=1. Swapping one for the other
  changes every polygenic score by a factor of sqrt(n/(n-1)).
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

# Where the pre-migration APPFL_GWAS checkout lives.
#
# Supplied by the environment rather than hardcoded: those trees are one person's working
# directories on one machine, and a path to them has no business in a public repository.
# Anyone else running these tests points the variable at their own copy -- and on a
# machine without one, the legacy-tree tests skip.
#
#     export APPFL_BIO_SUITE_LEGACY_GWAS=/path/to/APPFL_GWAS/APPFL
#
# The committed golden checksums in tests/fixtures/ are the durable record; this harness
# only exists to regenerate or re-confirm them while the original trees still exist.
_LEGACY_ENV = "APPFL_BIO_SUITE_LEGACY_GWAS"

LEGACY_GWAS_REPO = Path(os.environ.get(_LEGACY_ENV, "/nonexistent"))
LEGACY_SIM_DIR = LEGACY_GWAS_REPO / "local" / "data_sim"
LEGACY_CONFIG_DIR = LEGACY_GWAS_REPO / "examples" / "resources" / "configs" / "gwas"

# Outputs the legacy pipeline produces, in the order the pipeline writes them. The port
# must produce all of these, with identical content.
PHENOTYPE_OUTPUTS = (
    "covariates.csv",
    "phenotypes_gwas.csv",
    "phenotypes_pgs_eval.csv",
    "pgs_scores.csv",
    "score_t2d.txt",
    "score_bmi.txt",
    "simulation_summary.csv",
    "variant_counts.csv",
)

# Per-site files the bundler writes. The GWAS client loader requires the first six; the
# rest are simulation and central-baseline artifacts.
SITE_OUTPUTS = (
    "EUR.synthetic.100k.ld.maf.bed",
    "EUR.synthetic.100k.ld.maf.bim",
    "EUR.synthetic.100k.ld.maf.fam",
    "phenotypes_gwas.csv",
    "phenotypes_pgs_eval.csv",
    "covariates.csv",
    "pgs_scores.csv",
    "score_t2d.txt",
    "score_bmi.txt",
)


def legacy_tree_available() -> bool:
    """Whether the pre-migration source trees are reachable.

    They are not part of this repository, so tests that need them skip when
    ``APPFL_BIO_SUITE_LEGACY_GWAS`` is unset or points somewhere without them. Once the
    migration is complete and the trees are retired, the golden checksums remain as the
    record and these tests stop being runnable -- which is exactly why the checksums are
    committed rather than regenerated on demand.
    """
    return (LEGACY_SIM_DIR / "simulate_phenotypes.py").is_file()


def build_scratch_tree(root: Path, site_sizes: list[int] | None = None) -> Path:
    """Recreate the directory shape the legacy scripts expect, under ``root``.

    Returns the path standing in for the original ``local/`` directory.
    """
    local = root / "local"
    (local / "data_sim" / "input").mkdir(parents=True, exist_ok=True)
    (local / "data_sim" / "output" / "data").mkdir(parents=True, exist_ok=True)
    (local / "sites").mkdir(parents=True, exist_ok=True)

    # The scripts fall back to <local>/../examples/resources/configs/gwas for gwas_config.
    config_dir = root / "examples" / "resources" / "configs" / "gwas"
    config_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(LEGACY_CONFIG_DIR / "gwas_config.py", config_dir / "gwas_config.py")

    # Data_Sim_Scaling=1.0 so the harness exercises the full variant set. The live
    # environment ships 0.02, a smoke fraction, which would leave most of the pipeline
    # untested.
    (config_dir / "gwas_env.env").write_text(
        "Use_cuML=false\nVariant_Scaling=1.0\nHit_P_Threshold=5e-8\nData_Sim_Scaling=1.0\n",
        encoding="utf-8",
    )

    shutil.copy2(LEGACY_SIM_DIR / "simulate_phenotypes.py", local / "data_sim")

    bundler_source = _bundler_source()
    bundler_text = bundler_source.read_text(encoding="utf-8")
    if site_sizes is not None:
        # The ONLY modification the harness makes, and only because the fixture has 240
        # samples rather than 100,000. The legacy bundler asserts its five hardcoded
        # sizes sum to exactly 100,000 -- which is itself the defect the port fixes by
        # moving sizes into a scenario config.
        #
        # Substituting the literal keeps every other line, including the RNG calls and
        # the slicing, byte-identical to the original.
        bundler_text = bundler_text.replace(
            "SITE_SIZES = [18032, 9237, 25028, 40752, 6951]",
            f"SITE_SIZES = {site_sizes}",
        ).replace(
            'assert sum(SITE_SIZES) == 100000, "Site sizes must sum to 100,000"',
            f"assert sum(SITE_SIZES) == {sum(site_sizes)}",
        )
    (local / "data_sim" / "data_bundler.py").write_text(bundler_text, encoding="utf-8")

    return local


def _bundler_source() -> Path:
    """Locate data_bundler.py.

    It exists in two places in the pre-migration layout. They are byte-identical
    (verified), and collapsing them to one copy is part of the migration.
    """
    candidates = [
        LEGACY_SIM_DIR / "data_bundler.py",
        LEGACY_GWAS_REPO.parent / "TEMP_5site_fullscale" / "data_sim" / "data_bundler.py",
    ]
    if not LEGACY_GWAS_REPO.is_dir():
        raise FileNotFoundError(
            f"the pre-migration tree is not available. Set ${_LEGACY_ENV} to point at "
            "an APPFL_GWAS checkout, or rely on the committed golden checksums in "
            "tests/fixtures/gwas_legacy_golden.json instead."
        )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"no data_bundler.py in any of: {candidates}")


def run_legacy_stage(local: Path, script: str, timeout: int = 900) -> subprocess.CompletedProcess:
    """Run one legacy script in the scratch tree."""
    env = dict(os.environ)
    env["PYTHONNOUSERSITE"] = "1"
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env[var] = "1"
    # Deliberately NOT set: the scripts' own fallback path discovery is part of what is
    # being characterized.
    env.pop("GWAS_PROJECT_DIR", None)

    return subprocess.run(
        [sys.executable, str(local / "data_sim" / script)],
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
        cwd=str(local),
    )


def run_legacy_pipeline(root: Path, site_sizes: list[int]) -> Path:
    """Run both legacy stages end to end. Returns the scratch ``local/`` directory."""
    from make_fixture import build

    local = build_scratch_tree(root, site_sizes=site_sizes)
    build(local, n_samples=sum(site_sizes))

    for script in ("simulate_phenotypes.py", "data_bundler.py"):
        result = run_legacy_stage(local, script)
        if result.returncode != 0:
            raise RuntimeError(
                f"legacy {script} failed (exit {result.returncode})\n"
                f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
            )
    return local


def collect_outputs(local: Path, n_sites: int) -> dict[str, str]:
    """Checksum every pipeline output, keyed by a stable relative name."""
    from appfl_bio_suite.core.simulation import file_checksum

    out: dict[str, str] = {}
    data_dir = local / "data_sim" / "output" / "data"
    for name in PHENOTYPE_OUTPUTS:
        path = data_dir / name
        if path.is_file():
            out[f"pooled/{name}"] = file_checksum(path)

    for index in range(1, n_sites + 1):
        site_dir = local / "sites" / f"Site{index}" / "data"
        for name in SITE_OUTPUTS:
            path = site_dir / name
            if path.is_file():
                out[f"Site{index}/{name}"] = file_checksum(path)
    return out
