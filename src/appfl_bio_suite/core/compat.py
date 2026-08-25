"""Bridge one incompatibility between the two versions this suite pins.

THIS MODULE IS MEANT TO BE DELETED. `tests/test_upstream_shim.py` fails once it is no
longer needed, so it cannot quietly outlive its purpose.

The problem
-----------
The suite pins ``appfl==1.10.0`` and ``globus-compute-sdk==4.9.0``. Both pins are
correct in isolation and both are required:

* globus-compute-sdk/-endpoint must be exactly 4.9.0 across the whole federation. 4.12
  redesigned ``endpoint start`` to always launch the multi-user manager, which rejects
  the single-user ``engine``-in-config layout that every coordinator endpoint uses.
* appfl 1.10.0 is the current release and the one this suite is written against.

They do not work together. ``appfl-1.10.0-py3-none-any.whl`` contains, at module scope in
``appfl/comm/globus_compute/globus_compute_server_communicator.py`` (line 18)::

    from globus_compute_sdk.sdk.login_manager import AuthorizerLoginManager

``globus_compute_sdk.sdk.login_manager`` was removed in the 4.9 line. So::

    pip install "appfl==1.10.0" "globus-compute-sdk==4.9.0"
    python -c "from appfl.comm.globus_compute import GlobusComputeServerCommunicator"
    # ModuleNotFoundError: No module named 'globus_compute_sdk.sdk.login_manager'

which takes out the Globus Compute driver path -- the production path for every
experiment in this suite -- before any user code runs.

Why a shim rather than a fork or a vendored copy
------------------------------------------------
The import is dead weight on our path. ``AuthorizerLoginManager`` is used in exactly one
branch of APPFL's communicator, guarded by ``if "compute_token" in kwargs and
"openid_token" in kwargs`` -- the hosted APPFLx backend's delegated-token flow. A
coordinator driving their own endpoints never reaches it. Upstream's own fix for this is
to move the import inside that branch, which is a five-line change.

So the entire defect is a module-level import of something that only a path we do not use
would need. Registering a placeholder is proportionate to that. Forking APPFL to change
one import line, or vendoring a ~700-line communicator subclass to shadow it, would not
be.

The placeholder is not silently wrong: if anything ever does instantiate it -- meaning
someone genuinely took the hosted-token path -- it raises with an explanation rather than
failing obscurely later.

Removing this module
--------------------
When APPFL ships a release with the import moved into the branch that uses it:

1. Bump the ``appfl`` pin in ``constraints.txt`` and ``pyproject.toml``.
2. ``tests/test_upstream_shim.py`` starts failing. That is the signal, not a regression.
3. Delete this file, its call in ``appfl_bio_suite/__init__.py``, and that test.
"""

from __future__ import annotations

import importlib
import importlib.util
import sys
import types
from pathlib import Path

__all__ = [
    "shim_is_needed",
    "ensure_appfl_globus_compute_importable",
    "appfl_install_is_editable",
]

_SHIMMED_MODULE = "globus_compute_sdk.sdk.login_manager"

_COMMUNICATOR_RELPATH = Path("comm/globus_compute/globus_compute_server_communicator.py")

_HOSTED_TOKEN_PATH_MESSAGE = (
    "AuthorizerLoginManager was instantiated, but this environment has "
    "globus-compute-sdk 4.9.x, which does not ship it.\n"
    "\n"
    "This means something took APPFL's hosted-backend token-auth path -- the branch "
    "guarded by `compute_token`/`openid_token` in GlobusComputeServerCommunicator. "
    "appfl_bio_suite does not use that path: a coordinator authenticates as themselves "
    "and dispatches to endpoints directly.\n"
    "\n"
    "If you need the hosted APPFLx backend, you need a globus-compute-sdk that still "
    "provides this module -- which conflicts with the 4.9.0 pin the rest of the "
    "federation depends on. Resolve that deliberately; do not paper over it here.\n"
    "\n"
    "See appfl_bio_suite/core/compat.py."
)


def _real_module_exists() -> bool:
    """True when globus-compute-sdk genuinely provides the module."""
    if _SHIMMED_MODULE in sys.modules:
        return not getattr(sys.modules[_SHIMMED_MODULE], "__appfl_bio_suite_shim__", False)
    try:
        return importlib.util.find_spec(_SHIMMED_MODULE) is not None
    except (ImportError, AttributeError, ValueError):
        # A parent package that itself fails to import means the module is not available,
        # which is the answer the caller needs.
        return False


def shim_is_needed() -> bool:
    """True when the installed appfl/globus-compute-sdk pair needs the workaround.

    Both halves have to hold: APPFL must still import the module at module scope, and the
    SDK must not provide it. Used by ``tests/test_upstream_shim.py`` to fail loudly once
    this file has outlived its purpose.
    """
    if _real_module_exists():
        return False
    return _appfl_imports_login_manager_at_module_scope()


def _appfl_package_root() -> Path | None:
    """Locate the installed ``appfl`` package directory without importing it.

    ``find_spec`` is called on the top-level package only. Passing a dotted submodule
    path would import every parent package on the way -- including ``appfl.misc.utils``,
    which imports torch -- and importing is precisely what this module cannot rely on.
    """
    try:
        spec = importlib.util.find_spec("appfl")
    except (ImportError, AttributeError, ValueError):
        return None
    if spec is None or not spec.submodule_search_locations:
        return None
    return Path(next(iter(spec.submodule_search_locations)))


def appfl_install_is_editable() -> bool:
    """True when ``appfl`` resolves to a source checkout rather than a normal install.

    Matters because an editable checkout may carry local patches. In that case
    :func:`shim_is_needed` describes *that working tree*, not the pinned release the suite
    declares a dependency on -- so it is not a valid signal for whether the shim can be
    deleted. Tests use this to skip rather than report a misleading answer.
    """
    root = _appfl_package_root()
    if root is None:
        return False
    return "site-packages" not in root.parts and "dist-packages" not in root.parts


def _appfl_imports_login_manager_at_module_scope() -> bool:
    """Check whether the installed APPFL still has the top-level import.

    Read as source rather than imported, because importing is the thing that fails.
    """
    root = _appfl_package_root()
    if root is None:
        return False
    source = root / _COMMUNICATOR_RELPATH
    if not source.is_file():
        return False

    try:
        with open(source, encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                # A module-scope import has no leading whitespace. Once the import moves
                # inside the branch that uses it, it becomes indented and stops matching.
                if line[:1].isspace():
                    continue
                if stripped.startswith(f"from {_SHIMMED_MODULE} import"):
                    return True
    except OSError:
        return False
    return False


def ensure_appfl_globus_compute_importable() -> bool:
    """Make ``appfl.comm.globus_compute`` importable on globus-compute-sdk 4.9.x.

    Idempotent, and a no-op whenever the real module is available -- so it stays harmless
    on any future version pair where the problem is gone.

    Returns True if a placeholder was installed, False if nothing was needed.
    """
    if _real_module_exists():
        return False
    if sys.modules.get(_SHIMMED_MODULE) is not None:
        return True  # already shimmed by an earlier call

    class AuthorizerLoginManager:  # noqa: N801 - name is fixed by APPFL's import
        """Placeholder for a class globus-compute-sdk 4.9 does not ship.

        Exists to satisfy an unused module-level import. Raises if actually used.
        """

        def __init__(self, *args, **kwargs):
            raise RuntimeError(_HOSTED_TOKEN_PATH_MESSAGE)

    module = types.ModuleType(_SHIMMED_MODULE)
    module.AuthorizerLoginManager = AuthorizerLoginManager
    module.__appfl_bio_suite_shim__ = True
    module.__doc__ = (
        "Placeholder installed by appfl_bio_suite.core.compat. Not part of "
        "globus-compute-sdk. See that module for why it exists and when to delete it."
    )

    sys.modules[_SHIMMED_MODULE] = module

    # `from globus_compute_sdk.sdk.login_manager import X` resolves the parent package
    # first, so the attribute has to be reachable there too.
    try:
        parent = importlib.import_module("globus_compute_sdk.sdk")
        if not hasattr(parent, "login_manager"):
            parent.login_manager = module
    except ImportError:
        # globus-compute-sdk is not installed at all. Nothing more to do here; the caller
        # gets a clear ImportError naming the real missing package.
        pass

    return True
