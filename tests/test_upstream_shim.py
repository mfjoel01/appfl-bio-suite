"""Guard the compat shim so it gets deleted rather than rotting.

`core/compat.py` works around appfl 1.10.0 importing a globus-compute-sdk module that
4.9.0 removed. That workaround is correct today and should not outlive the defect.

The load-bearing test here is `test_shim_is_still_necessary`, which FAILS once upstream
fixes the import. A failure is the signal to delete compat.py, not a regression to
investigate.
"""

from __future__ import annotations

import importlib
import sys

import pytest

from appfl_bio_suite.core import compat


def test_shim_is_still_necessary():
    """Fails when the workaround can be removed. That is the point of this test.

    If this fails, upstream APPFL no longer imports
    globus_compute_sdk.sdk.login_manager at module scope, or the SDK ships it again.
    Either way:

      1. Delete src/appfl_bio_suite/core/compat.py
      2. Delete its call in src/appfl_bio_suite/__init__.py
      3. Delete this file
      4. Note the removal in docs/coordinator/releasing.md

    Skipped against an editable appfl checkout, which may carry local patches and so says
    nothing about the pinned release. CI installs from the pins, where this always runs.
    """
    if compat.appfl_install_is_editable():
        pytest.skip(
            "appfl is an editable checkout, which may already be patched locally. "
            "Whether the shim can be dropped is a question about the pinned release, "
            "so this check is only meaningful against a normal install."
        )
    assert compat.shim_is_needed(), (
        "The compat shim is no longer needed -- delete it.\n"
        "See this test's docstring for the four-step removal."
    )


def test_appfl_globus_compute_communicator_imports():
    """The whole point: APPFL's Globus Compute driver path must import.

    This is the production path for every experiment in the suite. Without the shim it
    raises ModuleNotFoundError before any user code runs.

    Skips when torch cannot load at all. That is an environment problem -- a GPU
    transport library missing on the current host, say -- and it is unrelated to the
    shim. Letting it fail here would report a broken workaround when the workaround is
    fine, which is worse than not running the check.
    """
    pytest.importorskip("appfl", reason="appfl not installed")
    try:
        importlib.import_module("torch")
    except Exception as exc:  # noqa: BLE001 - torch raises OSError, not ImportError
        pytest.skip(f"torch cannot load in this environment, so appfl cannot: {exc}")

    module = importlib.import_module("appfl.comm.globus_compute")
    assert hasattr(module, "GlobusComputeServerCommunicator")


def test_the_shimmed_module_itself_imports():
    """The narrow thing the shim is responsible for, checkable without torch.

    Separated from the test above so that a torch problem cannot mask a shim problem,
    and vice versa.
    """
    compat.ensure_appfl_globus_compute_importable()
    from globus_compute_sdk.sdk.login_manager import AuthorizerLoginManager

    assert AuthorizerLoginManager is not None


def test_shim_is_idempotent():
    """Calling it repeatedly must not stack placeholders or raise."""
    first = compat.ensure_appfl_globus_compute_importable()
    second = compat.ensure_appfl_globus_compute_importable()
    assert first == second


def test_placeholder_raises_a_useful_error_if_actually_used():
    """A silent no-op class would be worse than the original crash.

    If someone genuinely takes APPFL's hosted-token path, they must get a message that
    explains the version conflict, not an obscure failure three frames deeper.
    """
    compat.ensure_appfl_globus_compute_importable()
    module = sys.modules.get("globus_compute_sdk.sdk.login_manager")
    if module is None or not getattr(module, "__appfl_bio_suite_shim__", False):
        pytest.skip("real globus-compute-sdk module present; nothing was shimmed")

    with pytest.raises(RuntimeError, match="hosted-backend token-auth path"):
        module.AuthorizerLoginManager()


def test_shim_marks_itself_as_not_the_real_module():
    """The marker attribute is how shim_is_needed() avoids fooling itself."""
    compat.ensure_appfl_globus_compute_importable()
    module = sys.modules.get("globus_compute_sdk.sdk.login_manager")
    if module is None:
        pytest.skip("globus-compute-sdk not installed")
    if getattr(module, "__appfl_bio_suite_shim__", False):
        assert compat._real_module_exists() is False
