"""The three places that carry the version must agree.

Partners install from a git tag, and `INSTALL_TAG` -- which is what the generated setup
documents tell them to use -- is derived from `__version__`. If `pyproject.toml` says
0.2.0 and `__version__` still says 0.1.0, the wheel is one version and every bundle
points at another; if `CITATION.cff` lags, the citation records software nobody ran.

docs/coordinator/releasing.md lists keeping these in step as a checklist item. This is
that item, done mechanically, because a checklist is the thing you skip at the end of a
release.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import appfl_bio_suite

ROOT = Path(__file__).resolve().parent.parent


def test_pyproject_citation_and_dunder_version_agree():
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    packaged = pyproject["project"]["version"]

    citation = re.search(
        r"^version:\s*(\S+)\s*$", (ROOT / "CITATION.cff").read_text(encoding="utf-8"), re.M
    )
    assert citation, "CITATION.cff has no top-level `version:` line"

    assert packaged == appfl_bio_suite.__version__ == citation.group(1), (
        "the version is not the same in all three places:\n"
        f"  pyproject.toml            {packaged}\n"
        f"  appfl_bio_suite.__version__ {appfl_bio_suite.__version__}\n"
        f"  CITATION.cff              {citation.group(1)}\n"
        "\nSee the release checklist in docs/coordinator/releasing.md."
    )


def test_install_tag_is_the_taggable_form_of_the_version():
    """What a partner is told to clone must be what the tag will be called."""
    assert appfl_bio_suite.INSTALL_TAG == f"v{appfl_bio_suite.__version__}"
