"""The vendored binaries must be findable the same way everywhere.

`preflight --check env` resolved PLINK and SuSiEx as "on PATH *or* in vendor/bin", while
the code that shells out to them only ever consulted PATH. A green preflight was
therefore followed by `PlinkError: PLINK binary not found on PATH: plink` on the next
command -- which is exactly how CI failed, after preflight had passed.

Importing the package now prepends vendor/bin, so there is one answer instead of two.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

import appfl_bio_suite  # noqa: F401  -- imported for its PATH side effect
from appfl_bio_suite.core.experiments import REGISTRY, repo_root

VENDOR_BIN = repo_root() / "vendor" / "bin"

needs_vendor = pytest.mark.skipif(
    not VENDOR_BIN.is_dir(),
    reason="nothing vendored here; run scripts/fine-mapping/install_{plink,susiex}.sh",
)


@needs_vendor
def test_importing_the_package_puts_vendor_bin_on_path():
    assert str(VENDOR_BIN) in os.environ.get("PATH", "").split(os.pathsep)


@needs_vendor
def test_vendored_binaries_are_reachable_by_plain_which():
    """`shutil.which` is what the vendored fedfm tree uses, so it is the one that counts."""
    wanted = {name for spec in REGISTRY.values() for name in spec.required_binaries}
    vendored = {p.name for p in VENDOR_BIN.iterdir() if p.is_file() and os.access(p, os.X_OK)}
    for name in sorted(wanted & vendored):
        assert shutil.which(name), (
            f"{name} is in {VENDOR_BIN} but shutil.which cannot find it; "
            "preflight would pass and the run would then fail"
        )


def test_the_path_entry_is_not_duplicated_on_reimport():
    """Import is not guaranteed to happen once per process; prepending twice is a smell."""
    import importlib

    importlib.reload(appfl_bio_suite)
    entries = os.environ.get("PATH", "").split(os.pathsep)
    assert entries.count(str(VENDOR_BIN)) <= 1


def test_an_installed_only_tree_is_a_no_op(monkeypatch, tmp_path: Path):
    """No checkout, nothing vendored: the import must not invent a PATH entry."""
    before = os.environ.get("PATH", "")
    monkeypatch.setattr(Path, "is_dir", lambda self: False)
    import importlib

    importlib.reload(appfl_bio_suite)
    assert os.environ.get("PATH", "") == before
