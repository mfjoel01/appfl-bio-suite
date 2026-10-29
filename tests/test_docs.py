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
    [
        "new-federation.md",
        "architecture.md",
        "troubleshooting.md",
        "reference-deployment.md",
        "releasing.md",
        "network-map.md",
        "ga4gh.md",
    ],
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

    assert f">={context['python_version']}" in pyproject, (
        f"partner documents say Python {context['python_version']}, which does not match "
        "requires-python in pyproject.toml"
    )


def test_troubleshooting_covers_the_expensive_failures():
    """These cost days each. If one is dropped, the knowledge is gone."""
    text = (DOCS / "coordinator" / "troubleshooting.md").read_text(encoding="utf-8").lower()
    for topic in [
        "anchor",  # ^ $ in identity mapping -- the most expensive one
        "af_unix",  # TMPDIR overflow
        "can't start new thread",  # BLAS thread fan-out
        "idle_heartbeats",  # block teardown between rounds
        "unknown opcode",  # Python minor mismatch
        "login_manager",  # the appfl/SDK conflict
        "403",  # unprivileged endpoint start
        "422",  # mapping did not resolve
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
        assert any(_shares_structure(rendered, text) for text in source_texts), (
            f"{produced.name} in the bundle does not correspond to any file under "
            f"{partner_root}. Bundle content must come from docs/partner/ and nowhere "
            "else."
        )


def _shares_structure(rendered: str, source: str) -> bool:
    """Do these share enough literal headings to be the same document?"""

    def headings(text: str) -> set[str]:
        return {
            line.strip() for line in text.splitlines() if line.startswith("#") and "{{" not in line
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


def _partner_federation(tmp_path, **experiment_fields):
    """A one-partner FLamby federation that knows nothing about the partner's scheduler.

    The shipped example fills every field in, so tests about what a bundle does with a
    value the coordinator lacks build this instead. Skipping would make them agree with
    whatever the example happened to say, which is the one thing they must not do.
    """
    import yaml

    from appfl_bio_suite.core.config import load_federation as _load

    path = tmp_path / "federation.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "coordinator": {"identity": "a.researcher@example.edu"},
                # No `account`, no `partition`.
                "sites": [{"id": "partner", "name": "Partner University", "scheduler": "slurm"}],
                "experiments": {
                    "flamby-heart-disease": {
                        "enabled": True,
                        "service_account": "flamby_svc",
                        "endpoint_name": "flamby-ep",
                        "dataset": "HeartDisease",
                        "num_clients": 4,
                        "rounds": 2,
                        "sites": [
                            {
                                "site": "partner",
                                "client_id": "Site1",
                                "endpoint_uuid": "bbbbbbbb-0000-0000-0000-000000000001",
                                "center": 0,
                                "output_dir": "/home/flamby_svc/outputs",
                                "expected_train_samples": 199,
                            }
                        ],
                        **experiment_fields,
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    return _load(path)


def test_unknown_scheduler_fields_are_marked_not_left_blank(tmp_path):
    """`account:` with nothing after it parses as null and fails on the partner's cluster.

    federation.yaml legitimately may not carry a partner's scheduler account or queue --
    the coordinator often does not know them before onboarding. Rendering them empty made
    the endpoint template look complete and install cleanly, then fail at job submission
    with an error naming neither the field nor the file. An explicit marker is the honest
    form, and the coordinator is told which fields carry one.

    The setup guide is checked too. It used to show its own hand-kept copy of the
    template, which kept rendering the blanks after the file itself was fixed.
    """
    import yaml

    from appfl_bio_suite.core.partner import generate_bundle, unresolved_fields

    federation = _partner_federation(tmp_path)

    assert unresolved_fields(federation, "flamby-heart-disease", "partner") == [
        "account",
        "partition",
    ]

    destination = generate_bundle(federation, "flamby-heart-disease", "partner", tmp_path / "out")
    template = (destination / "user_config_template.yaml.j2").read_text(encoding="utf-8")
    provider = yaml.safe_load(template)["engine"]["provider"]

    for field in ("account", "partition"):
        assert provider[field], f"{field} rendered empty rather than as an explicit marker"
        assert str(provider[field]).startswith("REPLACE_WITH_"), provider[field]

    guide = (destination / "1-endpoint-setup.md").read_text(encoding="utf-8")
    assert template.strip() in guide, "the setup guide shows a template other than the file"
    assert not re.search(r"^\s*(account|partition):\s*$", guide, re.M), "guide has a blank"
    assert not re.search(r"^\s*module load\s*$", guide, re.M), "guide loads an empty module"


def test_setup_guide_creates_the_environment_the_workers_activate(tmp_path):
    """Step 1 must create, by name, the environment worker_init activates.

    They used to disagree: Step 1 led with a venv at ~/appfl-env while worker_init ran
    `conda activate appfl_env`, so a partner following the first instruction built an
    endpoint whose workers could not start. And the name was fixed, so a partner already
    hosting an earlier experiment's endpoint had to collide with that endpoint's
    environment or install over it.
    """
    from appfl_bio_suite.core.partner import generate_bundle

    federation = _partner_federation(tmp_path, conda_env="finemap_env")
    destination = generate_bundle(federation, "flamby-heart-disease", "partner", tmp_path / "out")
    guide = (destination / "1-endpoint-setup.md").read_text(encoding="utf-8")
    template = (destination / "user_config_template.yaml.j2").read_text(encoding="utf-8")

    assert "conda create -n finemap_env " in guide
    assert "conda activate finemap_env\n" in template
    assert "conda activate finemap_env &&" in guide, "the service-account check tests another env"
    for stale in ("venv", "appfl_env", "appfl-env"):
        assert stale not in guide, f"the guide still mentions {stale!r}"

    # Part 2 installs FLamby into that same environment; naming another one puts the
    # dataset where the workers cannot import it.
    part2 = (destination / "2-flamby-heart-disease.md").read_text(encoding="utf-8")
    assert "conda activate finemap_env " in part2
    assert "appfl_env" not in part2, "Part 2 activates a different environment than Part 1"


def test_setup_guide_mapping_is_the_installable_file(tmp_path):
    """The JSON Step 4b shows must parse, and must be exactly the file the bundle ships.

    It used to be assembled by pasting the regex between quotes, which rendered
    `"a\\.researcher@example\\.edu"` with single backslashes -- an invalid JSON escape. A
    partner who followed "replace the entire contents with the block below" installed a
    file the mapper could not read, while the bundle's own copy was fine.
    """
    import json
    import re

    from appfl_bio_suite.core.partner import generate_bundle

    destination = generate_bundle(
        _partner_federation(tmp_path), "flamby-heart-disease", "partner", tmp_path / "out"
    )
    guide = (destination / "1-endpoint-setup.md").read_text(encoding="utf-8")
    shipped = (destination / "example_identity_mapping_config.json").read_text(encoding="utf-8")

    blocks = re.findall(r"```json\n(.*?)```", guide, re.S)
    assert len(blocks) == 1, "expected exactly one JSON block in Part 1"
    assert json.loads(blocks[0]) == json.loads(shipped)
    assert blocks[0] == shipped

    assert _partner_federation(tmp_path).experiment("flamby-heart-disease").conda_env == (
        "appfl_env"
    ), "the default changed, which would re-point every existing partner's workers"
