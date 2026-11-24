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

__version__ = "0.3.2"

# Where partners install from: PyPI, at exactly this version.
#
# An exact pin, never a range. A federation whose sites resolved `>=` on different days
# is running different code with no way to say which. PyPI never lets a published
# version's files be replaced, so every site installing `==X.Y.Z` runs the same wheel,
# for the same reason the Globus Compute stack is pinned in constraints.txt: skew here
# does not fail at install time, it fails as a deserialization error partway into
# someone else's run.
#
# The matching git tag (INSTALL_TAG) is still cut for every release. The publish
# workflow builds the PyPI distributions from it, and coordinators check it out for the
# scripts and example files that are not in the wheel. The release checklist in
# docs/coordinator/releasing.md keeps the version, the tag and the PyPI upload in step.
REPO_URL = "https://github.com/mfjoel01/appfl-bio-suite"
INSTALL_TAG = f"v{__version__}"


def install_spec(extras: str = "") -> str:
    """The pip requirement a partner installs: this exact release, from PyPI.

    ``extras`` is a comma-separated string ("gwas", "flamby"), matching the shape the
    experiment registry already stores. Empty means no extras -- ``pkg[]`` is not valid
    requirement syntax, so it must be omitted rather than rendered empty.
    """
    name = f"appfl-bio-suite[{extras}]" if extras else "appfl-bio-suite"
    return f"{name}=={__version__}"


def release_file_url(path: str) -> str:
    """A repository file as of this release's tag, for installs that have no checkout.

    pip accepts URLs for ``-c`` and ``-r``, so a PyPI install can apply this release's
    ``constraints.txt`` or ``requirements-multisusie.txt`` without cloning anything.
    """
    raw = REPO_URL.replace("https://github.com/", "https://raw.githubusercontent.com/")
    return f"{raw}/{INSTALL_TAG}/{path}"


def multisusie_install_hint() -> str:
    """How to install the optional MultiSuSiE backend, which PyPI metadata cannot carry."""
    return (
        "Install the coordinator extra and the pinned MultiSuSiE backend:\n"
        f"    pip install '{install_spec('multisusie')}' "
        f"-r {release_file_url('requirements-multisusie.txt')}"
    )


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

__all__ = [
    "INSTALL_TAG",
    "REPO_URL",
    "__version__",
    "install_spec",
    "multisusie_install_hint",
    "release_file_url",
]
