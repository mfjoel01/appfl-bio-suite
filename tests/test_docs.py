"""Documentation structure, drift, and the partner-bundle closed set.

Three things are checked, each guarding a failure that has actually happened on this
project:

1. **Every experiment has all four documents, non-empty.** Uniformity is the point: no
   experiment gets a special extra document and none is omitted, so a reader always knows
   where to look.

2. **Documented values match shipped configuration.** Documentation quotes version pins
   and expected counts. When those drift, the person who finds out is a partner whose
   cluster does not match the guide.

3. **Bundles read only from docs/partner/.** That directory is a closed set. Enforcing it
   in code means "did I accidentally send my internal notes" stops being a judgment call.
"""

from __future__ import annotations

import re

import pytest

from appfl_bio_suite.core.config import EXPERIMENT_NAMES, load_federation
from appfl_bio_suite.core.experiments import REGISTRY, repo_root

DOCS = repo_root() / "docs"
REQUIRED_EXPERIMENT_DOCS = ("ABOUT.md", "DATA.md", "RUNBOOK.md")

pytestmark = pytest.mark.skipif(
    not DOCS.is_dir(), reason="docs/ is not packaged; these tests need a source checkout"
)


# ---------------------------------------------------------------------------
# structure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("experiment", EXPERIMENT_NAMES)
@pytest.mark.parametrize("doc", REQUIRED_EXPERIMENT_DOCS)
def test_experiment_has_all_four_docs(experiment: str, doc: str):
    """Uniform across every experiment, implemented or not."""
    path = DOCS / "experiments" / experiment / doc
    assert path.is_file(), (
        f"{path} is missing.\n"
        "Every experiment has exactly these documents: ABOUT.md (what it is and the "
        "record of runs), DATA.md (how the data comes into existence), RUNBOOK.md (how "
        "to run it), plus docs/partner/experiments/<name>.md. None is optional -- an "
        "unimplemented experiment gets one that says so."
    )
    assert path.stat().st_size > 200, f"{path} exists but is essentially empty"


@pytest.mark.parametrize("experiment", EXPERIMENT_NAMES)
def test_experiment_has_a_partner_doc(experiment: str):
    """The fourth document. Part 2 of what a partner receives."""
    path = DOCS / "partner" / "experiments" / f"{experiment}.md"
    assert path.is_file(), f"{path} is missing"
    assert path.stat().st_size > 200, f"{path} exists but is essentially empty"


def test_shared_partner_doc_exists():
    """Part 1 is shared by every experiment, which is why it is not per-experiment."""
    path = DOCS / "partner" / "endpoint-setup.md"
    assert path.is_file()
    assert path.stat().st_size > 1000


@pytest.mark.parametrize(
    "name",
    ["new-federation.md", "architecture.md", "troubleshooting.md",
     "reference-deployment.md", "releasing.md"],
)
def test_coordinator_docs_exist(name: str):
    path = DOCS / "coordinator" / name
    assert path.is_file(), f"{path} is missing"
    assert path.stat().st_size > 500


def test_docs_index_exists():
    assert (DOCS / "README.md").is_file()


def test_unimplemented_experiments_say_so():
    """A planned experiment's docs must not read as though it works."""
    for name, spec in REGISTRY.items():
        if spec.implemented:
            continue
        about = (DOCS / "experiments" / name / "ABOUT.md").read_text(encoding="utf-8")
        assert re.search(r"not implemented", about, re.IGNORECASE), (
            f"{name} is registered as unimplemented, but its ABOUT.md does not say so. "
            "A reader must not have to infer it."
        )


# ---------------------------------------------------------------------------
# drift
# ---------------------------------------------------------------------------


def _constraints() -> dict[str, str]:
    pins = {}
    path = repo_root() / "constraints.txt"
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if "==" in line:
            name, _, version = line.partition("==")
            pins[name.strip().lower()] = version.strip()
    return pins


def test_partner_doc_quotes_the_real_globus_compute_pin():
    """The version a partner is told to verify must be the one we actually pin.

    Every site in a federation must match. A guide quoting a stale version means a
    partner installs the wrong one and the mismatch surfaces mid-run.
    """
    from appfl_bio_suite.core.partner import render_context

    federation = load_federation(repo_root() / "federation.yaml.example")
    context = render_context(federation, "gwas", "site-north")
    pinned = _constraints()["globus-compute-endpoint"]

    assert context["globus_compute_version"] == pinned, (
        f"partner documents tell sites to expect globus-compute-endpoint "
        f"{context['globus_compute_version']}, but constraints.txt pins {pinned}."
    )


def test_python_version_in_partner_context_matches_pyproject():
    from appfl_bio_suite.core.partner import render_context

    federation = load_federation(repo_root() / "federation.yaml.example")
    context = render_context(federation, "gwas", "site-north")
    pyproject = (repo_root() / "pyproject.toml").read_text(encoding="utf-8")

    assert f'>={context["python_version"]}' in pyproject, (
        f"partner documents say Python {context['python_version']}, which does not match "
        "requires-python in pyproject.toml"
    )


def test_troubleshooting_covers_the_expensive_failures():
    """These cost days each. If one is dropped, the knowledge is gone."""
    text = (DOCS / "coordinator" / "troubleshooting.md").read_text(encoding="utf-8").lower()
    for topic in [
        "anchor",            # ^ $ in identity mapping -- the most expensive one
        "af_unix",           # TMPDIR overflow
        "can't start new thread",  # BLAS thread fan-out
        "idle_heartbeats",   # block teardown between rounds
        "unknown opcode",    # Python minor mismatch
        "login_manager",     # the appfl/SDK conflict
        "403",               # unprivileged endpoint start
        "422",               # mapping did not resolve
    ]:
        assert topic in text, f"troubleshooting.md no longer covers {topic!r}"


# ---------------------------------------------------------------------------
# the closed set
# ---------------------------------------------------------------------------


def test_bundle_reads_only_partner_docs(tmp_path):
    """A bundle may contain documentation ONLY from docs/partner/.

    This is the structural guarantee that replaces remembering. Everything under
    docs/partner/ is safe to send; nothing outside it is.
    """
    from appfl_bio_suite.core.partner import generate_bundle, partner_docs_root

    federation = load_federation(repo_root() / "federation.yaml.example")
    destination = generate_bundle(federation, "gwas", "site-north", tmp_path)

    partner_root = partner_docs_root().resolve()
    sources = {p.resolve() for p in partner_root.rglob("*.md")}
    source_texts = {p.read_text(encoding="utf-8") for p in sources}

    for produced in destination.glob("*.md"):
        rendered = produced.read_text(encoding="utf-8")
        # Every rendered doc must derive from a docs/partner/ source. Compare on the
        # unrendered skeleton: template variables differ after rendering, so a structural
        # anchor from the source is the check that survives.
        assert any(
            _shares_structure(rendered, text) for text in source_texts
        ), (
            f"{produced.name} in the bundle does not correspond to any file under "
            f"{partner_root}. Bundle content must come from docs/partner/ and nowhere "
            "else."
        )


def _shares_structure(rendered: str, source: str) -> bool:
    """Do these share enough literal headings to be the same document?"""
    def headings(text: str) -> set[str]:
        return {
            line.strip()
            for line in text.splitlines()
            if line.startswith("#") and "{{" not in line
        }

    source_headings = headings(source)
    if not source_headings:
        return False
    return len(source_headings & headings(rendered)) >= max(1, len(source_headings) // 2)


def test_bundle_refuses_content_outside_partner_docs(tmp_path):
    """The guard must actually reject, not merely be documented."""
    from appfl_bio_suite.core.partner import BundleSecurityError, _assert_inside_partner_docs

    with pytest.raises(BundleSecurityError, match="may only come from"):
        _assert_inside_partner_docs(DOCS / "coordinator" / "troubleshooting.md")

    with pytest.raises(BundleSecurityError):
        _assert_inside_partner_docs(repo_root() / "local" / "notes.md")


def test_generated_bundle_has_no_placeholders(tmp_path):
    """Every placeholder left in a bundle is a question a partner has to ask."""
    from appfl_bio_suite.core.partner import generate_bundle

    federation = load_federation(repo_root() / "federation.yaml.example")
    for experiment, site in (("gwas", "site-north"), ("flamby-heart-disease", "site-east")):
        destination = generate_bundle(federation, experiment, site, tmp_path / experiment)
        for path in destination.rglob("*"):
            if path.suffix not in (".md", ".yaml", ".j2", ".json"):
                continue
            text = path.read_text(encoding="utf-8")
            if path.suffix == ".j2":
                continue  # a template by design; rendered when installed
            leftovers = re.findall(r"\{\{[^}]+\}\}", text)
            assert not leftovers, f"{path} still contains {leftovers}"
