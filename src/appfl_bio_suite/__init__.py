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

import os
from pathlib import Path

from appfl_bio_suite.core.compat import ensure_appfl_globus_compute_importable

__version__ = "0.1.0"

# Vendored command-line tools, put on PATH at import time.
#
# scripts/fine-mapping/install_plink.sh and install_susiex.sh drop static binaries into
# vendor/bin/, and preflight resolves them as "on PATH *or* in vendor/bin". The code that
# actually shells out -- run_plink in the vendored fedfm tree -- only ever looks at PATH.
# Those two disagreed, so a green `preflight --check env` was followed by
# `PlinkError: PLINK binary not found on PATH: plink` on the very next command.
#
# Prepending once here is what makes them agree, and it is done at import rather than in
# the CLI so that pytest, the drivers, and programmatic callers get it too. In an
# installed-only environment there is no checkout and nothing vendored, so this is a
# no-op. Prepended, not appended: a vendored build is chosen deliberately and should win
# over whatever a cluster module happens to expose.
_vendor_bin = Path(__file__).resolve().parents[2] / "vendor" / "bin"
if _vendor_bin.is_dir():
    _entries = os.environ.get("PATH", "").split(os.pathsep)
    if str(_vendor_bin) not in _entries:
        os.environ["PATH"] = os.pathsep.join([str(_vendor_bin), *filter(None, _entries)])

# Applied at import time, before anything can reach APPFL's Globus Compute communicator.
#
# appfl 1.10.0 imports a globus-compute-sdk module that 4.9.0 does not ship, and this
# suite pins both. A no-op on any version pair where that is not true. See core/compat.py
# for the full explanation and the deletion criteria.
ensure_appfl_globus_compute_importable()

__all__ = ["__version__"]
