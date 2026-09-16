#!/usr/bin/env python3
"""Run one vendored fine-mapping stage. A thin entry point, and thin on purpose.

WHY THIS FILE EXISTS AT ALL, RATHER THAN ``python -m``
-------------------------------------------------------
Three of these stages parallelize with joblib's loky backend, which pickles the
per-instance function **by qualified module name** and re-imports it in each worker.
Launched as ``python -m appfl_bio_suite...fine_mapping``, the module's ``__name__``
becomes ``__main__``; the workers then try to import ``__main__``, get the launcher
instead of the stage, and fail in a way that names neither.

Dispatching through a separate file keeps every stage module a normal import. The
standalone repository learned this the same way and kept the same shape.

EVERYTHING HERE IS COORDINATOR-SIDE. NONE OF IT IS THE FEDERATION.
-------------------------------------------------------------------
The federation is ``appfl-bio-suite run fine-mapping``. Nothing dispatched from this file
contacts a partner or a Globus endpoint; these are the stages that build the data and the
stages that produce the comparators a federated result is judged against.

Two groups:

**Building the data package.** ``appfl-bio-suite simulate fine-mapping`` does all of this
in one process, and for anything smaller than the published scenario that is the command
to use. The stages are exposed separately because two of them shard over loci, and the
published scenario -- 100 loci x 15 architectures x 10 replicates over 150,000
individuals -- is not a single-process job. See ``simulation/stages.py``.

**Analysing it.** ``centralized`` is the SuSiEx baseline: one column per ancestry, pooling
that ancestry's individuals across every site, which is the fit you would get with all the
data in one place. ``federated`` is the same statistics the APPFL aggregator runs, driven
from one process instead of over the wire -- useful for reproducing a federated result
locally, and for sharding a run larger than one exchange can carry. ``validation`` and
``qc`` check the simulated package itself.

USAGE
-----
Analysis stages take ``--config``, the ``pipeline_config.yaml`` that ``simulate`` (or
``prepare``) writes into the run directory. Arguments are passed through to the stage
unchanged, so ``--n-shards``/``--shard-index``/``--merge``/``--limit`` work as upstream::

    appfl-bio-suite simulate fine-mapping --scenario ci-tiny --out /scratch/fm

    scripts/fine-mapping/run_stage.py centralized --config /scratch/fm/pipeline_config.yaml
    scripts/fine-mapping/run_stage.py federated   --config /scratch/fm/pipeline_config.yaml \\
        --n-shards 4 --shard-index 0
    scripts/fine-mapping/run_stage.py federated   --config /scratch/fm/pipeline_config.yaml \\
        --merge --n-shards 4
    scripts/fine-mapping/run_stage.py figures \\
        --results /scratch/fm/reports/fine_mapping/fm_results.tsv

The sharded build, for the published scenario, takes ``--scenario``/``--out`` at the ends
and ``--config`` in the middle. scripts/fine-mapping/submit_polaris_simulate.pbs runs it.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# BLAS reads these at import time, so they have to be set before numpy is reachable.
# Every stage fans out with joblib and pins one thread per worker; without these a
# 32-worker shard tries to start 32x32 threads and a busy node refuses.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

# Importable from a checkout without installing. Harmless when the package is installed.
_REPO_ROOT = Path(__file__).resolve().parents[2]
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

# stage -> (module providing main(argv), one-line description).
#
# Ordered as a run goes. The simulation stages are only needed for a sharded run; a
# single-process one is `appfl-bio-suite simulate fine-mapping`, which does all five.
_STAGES = {
    # -- building the data package, stage by stage (sharded runs only) --------
    "prepare": (
        "appfl_bio_suite.experiments.fine_mapping.simulation.stages",
        "materialize the genotype pool, write pipeline_config.yaml",
        "prepare_main",
    ),
    "sampling": (
        "appfl_bio_suite.experiments.fine_mapping.fedfm.sampling",
        "cut disjoint per-site cohorts out of the pool",
        "main",
    ),
    "locus-selection": (
        "appfl_bio_suite.experiments.fine_mapping.fedfm.locus_selection",
        "score LD divergence and pick stratified loci  [shardable]",
        "main",
    ),
    "phenotypes": (
        "appfl_bio_suite.experiments.fine_mapping.fedfm.phenotype_sim",
        "causal variants, effect sizes, phenotypes  [shardable]",
        "main",
    ),
    "bundle": (
        "appfl_bio_suite.experiments.fine_mapping.simulation.stages",
        "assemble per-site bundles, check disjointness, write the manifest",
        "bundle_main",
    ),
    # -- analysing it ---------------------------------------------------------
    "centralized": (
        "appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping",
        "the centralized SuSiEx baseline  [shardable]",
        "main",
    ),
    "federated": (
        "appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping",
        "the standalone federated path, no APPFL  [shardable]",
        "main",
    ),
    "validation": (
        "appfl_bio_suite.experiments.fine_mapping.fedfm.validation",
        "the five-check self-validation harness",
        "main",
    ),
    "qc": (
        "appfl_bio_suite.experiments.fine_mapping.fedfm.qc",
        "per-site QC and HTML reports",
        "main",
    ),
    "figures": (
        "appfl_bio_suite.experiments.fine_mapping.plotting",
        "the four result figures",
        "main",
    ),
}


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        width = max(len(name) for name in _STAGES)
        print(__doc__)
        print("Stages:")
        for name, (_module, description, _entry) in _STAGES.items():
            print(f"  {name:<{width}}  {description}")
        print("\n  [shardable] stages take --n-shards N --shard-index I, then")
        print("  --merge --n-shards N to reduce the parts.")
        return 0 if argv else 2

    stage, rest = argv[0], argv[1:]
    if stage not in _STAGES:
        print(
            f"unknown stage '{stage}'. Known: {', '.join(_STAGES)}.",
            file=sys.stderr,
        )
        return 2

    import importlib
    import inspect

    module_name, _description, entry = _STAGES[stage]
    module = importlib.import_module(module_name)
    run = getattr(module, entry)

    # The vendored stages are not uniform about this: most take `argv`, but
    # validation.main() reads sys.argv directly. Set sys.argv either way -- it is what a
    # no-argument main() will parse, and it makes argparse's usage line and error
    # messages name the stage instead of this launcher.
    sys.argv = [f"run_stage.py {stage}", *rest]
    if inspect.signature(run).parameters:
        return int(run(rest) or 0)
    return int(run() or 0)


if __name__ == "__main__":
    raise SystemExit(main())
