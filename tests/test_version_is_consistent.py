"""The three places that carry the version must agree.

Partners install an exact version from PyPI, and `install_spec()` -- which is what the
generated setup documents tell them to use -- is derived from `__version__`. If
`pyproject.toml` says 0.2.0 and `__version__` still says 0.1.0, the wheel is one version
and every bundle points at another; if `CITATION.cff` lags, the citation records software
nobody ran.

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


def test_install_spec_pins_this_release_on_pypi():
    """A partner's step one: an exact PyPI pin, never a range, a branch or a URL."""
    version = appfl_bio_suite.__version__
    assert appfl_bio_suite.install_spec() == f"appfl-bio-suite=={version}"
    assert appfl_bio_suite.install_spec("gwas") == f"appfl-bio-suite[gwas]=={version}"


# Where a document names a concrete release: a PyPI pin, a wheel filename, a checkout of
# the tag, or a raw file fetched from the tag. Placeholders such as X.Y.Z do not match.
_RELEASE_REFERENCES = (
    re.compile(r"appfl-bio-suite(?:\[[^\]]*\])?==(\d+\.\d+\.\d+)"),
    re.compile(r"appfl_bio_suite-(\d+\.\d+\.\d+)-py3"),
    re.compile(r"--branch v(\d+\.\d+\.\d+)"),
    re.compile(r"appfl-bio-suite/v(\d+\.\d+\.\d+)/"),
)


def test_documented_install_commands_name_this_release():
    """Every copy-pasteable install command must install the version this tree is.

    Release notes describe their own release and are exempt; README and the guides are
    what someone follows today, so a stale `==0.2.1` there installs the wrong federation.
    """
    documents = [ROOT / "README.md", *sorted((ROOT / "docs").rglob("*.md"))]
    documents += sorted((ROOT / "examples").rglob("*.md"))
    documents += sorted((ROOT / "website").glob("*.html"))
    stale = [
        f"{path.relative_to(ROOT)}: {match.group(0)}"
        for path in documents
        for pattern in _RELEASE_REFERENCES
        for match in pattern.finditer(path.read_text(encoding="utf-8"))
        if match.group(1) != appfl_bio_suite.__version__
    ]
    assert not stale, (
        f"install commands that do not name {appfl_bio_suite.__version__}:\n  " + "\n  ".join(stale)
    )
