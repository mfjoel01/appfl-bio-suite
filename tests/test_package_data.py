"""Every file pyproject.toml declares as package data must actually be in the repo.

THE FAILURE MODE
----------------
Package data is the one kind of source file that a developer's own checkout will
happily supply from disk while git has never heard of it. Nothing imports it, so no
import fails; it is read at run time by path, and that path resolves locally. The
result passes every test on the machine that wrote it and fails on every machine that
clones it.

This is not hypothetical. ``.gitignore`` carried a bare ``data/`` rule meaning "no
experiment inputs or outputs in git". A bare directory pattern matches at any depth, so
it also matched ``src/appfl_bio_suite/core/ga4gh/data/`` -- the vendored DUO release.
The snapshot sat untracked for the life of the branch. Locally everything passed. In CI
the checkout had no ``duo.json``, and because consent evaluation is on the path of both
the fine-mapping simulation and config resolution, it surfaced as forty-odd unrelated
looking failures and a loopback that died in step 5 of 6.

The lesson that generalizes: a negation for such a rule must re-include the DIRECTORY.
git does not descend into an excluded directory, so ``!.../data/duo.json`` on its own is
silently a no-op -- the naive fix looks right and changes nothing.

WHY A TEST AND NOT A CONVENTION
-------------------------------
The author cannot see this bug by construction: their working tree is the one place the
file exists. It costs milliseconds to check, and the check is exact -- ask git, do not
infer from ``.gitignore``, because the rule interactions are precisely what went wrong.
"""

from __future__ import annotations

import subprocess
import tomllib

import pytest

from appfl_bio_suite.core.experiments import repo_root

ROOT = repo_root()
PACKAGE_ROOT = ROOT / "src" / "appfl_bio_suite"


def declared_globs() -> list[str]:
    """The package-data patterns, read from pyproject.toml rather than duplicated here."""
    with (ROOT / "pyproject.toml").open("rb") as fh:
        config = tomllib.load(fh)
    data = config["tool"]["setuptools"]["package-data"]
    return list(data["appfl_bio_suite"])


def tracked_files() -> set[str]:
    """Paths git knows about, relative to the repo root. Empty set if this is not a checkout."""
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return set()
    return {p for p in result.stdout.split("\0") if p}


requires_checkout = pytest.mark.skipif(
    not (ROOT / ".git").exists(),
    reason="not a git checkout -- nothing to compare against",
)


@pytest.mark.parametrize("pattern", declared_globs())
def test_every_declared_pattern_matches_something(pattern: str):
    """A glob that matches nothing ships nothing, and does it without complaining."""
    assert list(PACKAGE_ROOT.glob(pattern)), (
        f"package-data pattern {pattern!r} matches no file under {PACKAGE_ROOT}. "
        "Either the file moved and the pattern did not, or it was never added."
    )


@requires_checkout
@pytest.mark.parametrize("pattern", declared_globs())
def test_declared_package_data_is_tracked_by_git(pattern: str):
    """On disk is not enough. A file git does not carry does not survive a clone."""
    tracked = tracked_files()
    missing = sorted(
        str(path.relative_to(ROOT))
        for path in PACKAGE_ROOT.glob(pattern)
        if str(path.relative_to(ROOT)) not in tracked
    )
    assert not missing, (
        f"declared as package data but untracked: {missing}. These exist in this working "
        "tree and in no clone of it. Check .gitignore -- and note that re-including a "
        "file under an excluded directory requires negating the directory, not the file."
    )
