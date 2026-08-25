"""appfl-bio-suite: cross-silo federated learning experiments in computational biology.

A downstream extension package for APPFL. Custom trainers, aggregators, models, and
dataset loaders live here and are referenced by path from configuration; nothing needs to
be merged into APPFL for an experiment in this suite to run.

Layout:
    core/         experiment-agnostic infrastructure (endpoints, identity, preflight,
                  launch, partner bundles, simulation contract, config)
    experiments/  one package per experiment; the science lives here

Start at docs/coordinator/new-federation.md.
"""

from __future__ import annotations

__version__ = "0.1.0"

# Applied at import time, before anything can reach APPFL's Globus Compute communicator.
#
# appfl 1.10.0 imports a globus-compute-sdk module that 4.9.0 does not ship, and this
# suite pins both. A no-op on any version pair where that is not true. See core/compat.py
# for the full explanation and the deletion criteria.
from appfl_bio_suite.core.compat import ensure_appfl_globus_compute_importable

ensure_appfl_globus_compute_importable()

__all__ = ["__version__"]
