"""Validate a partner's Globus Compute identity-mapping file, offline, in about a second.

WHY THIS EXISTS
---------------
A multi-user endpoint decides whether to accept your tasks by matching your Globus
identity against an expression in a JSON file, then running the task as a local account.
Get that expression wrong and every submission fails with::

    422 SEMANTICALLY_INVALID -- Request payload failed validation:
    Identity failed to map to a local user name. (LookupError)

The failure is entirely on the partner's side of the wire, so the natural debugging loop
is: you submit, it fails, you email them, they change something, they wait for you to
re-test. That loop has taken over a week on this project, and it produced one confident
wrong diagnosis that a partner acted on before anyone noticed.

This module closes the loop locally. The partner runs one command against their own file
and gets PASS or FAIL immediately, with the compiled pattern printed, without needing
anyone at the coordinating site to be awake.

THE BUG IT PRIMARILY CATCHES
----------------------------
``match`` is not a full regex. ``ExpressionIdentityMapping._compile_match`` escapes
``^ $ + { } [ ]`` into *literal* characters and then wraps the expression in its own
``^...$``. So a mapping that reads::

    { "source": "{username}", "match": "^you@example\\.org$", "output": "svc" }

compiles to ``^\\^you@example\\.org\\$$`` -- a pattern matching only the literal string
including the caret and dollar sign. It can never match anything, with either source
field. Removing the anchors does not loosen the match; the mapper anchors it for you.

Two secondary traps this also catches, both of which produce a file that looks correct:

* ``map_identities`` silently skips any identity record whose ``status`` is not ``used``
  or ``private`` (expression.py:307). A hand-built test record without a status matches
  nothing and looks exactly like a bad regex.
* Escapes other than ``\\. \\? \\* \\| \\( \\) \\\\ \\0-\\9`` are rejected outright, so a
  pattern built with Python's ``re.escape`` (which escapes ``-``) makes the whole mapping
  document invalid.

WHAT IT CANNOT SEE
------------------
A ``403 ENDPOINT_ACCESS_FORBIDDEN`` rather than a 422. That is a different failure: the
endpoint was started by an unprivileged user, so identity mapping was skipped entirely
and only the starting identity may submit. The endpoint starts successfully anyway, with
no error, which is what makes it confusing. No amount of editing the JSON fixes it.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from appfl_bio_suite.core.config import FederationError

__all__ = ["MappingRule", "MappingReport", "validate_mapping_file", "build_mapping_document"]

# The mapper skips any identity record whose status is not one of these, so a test record
# must carry one or it matches nothing for reasons that have nothing to do with the regex.
_ACTIVE_STATUS = "used"


@dataclass
class MappingRule:
    """One rule from a mapping document, with what it actually compiles to."""

    source: str
    match: str
    output: str
    compiled: str | None = None
    error: str | None = None

    @property
    def has_anchors(self) -> bool:
        return "^" in self.match or "$" in self.match


@dataclass
class MappingReport:
    """The result of checking one mapping file against one identity."""

    path: Path
    identity: str
    expected_output: str
    rules: list[MappingRule] = field(default_factory=list)
    mapped_to: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.errors and self.expected_output in self.mapped_to

    def diagnosis(self) -> str:
        """The most likely cause, in the order these actually occur in practice."""
        if self.passed:
            return ""
        if self.errors:
            return self.errors[0]
        anchored = [r for r in self.rules if r.has_anchors]
        if anchored:
            return (
                "The `match` field contains ^ or $ anchors. This is the cause in nearly "
                "every case: the mapper escapes them into literal characters and then "
                "anchors the pattern itself, so the rule can never match. Remove them -- "
                "this does not loosen the match."
            )
        if self.mapped_to:
            return (
                f"The identity maps to {self.mapped_to}, but tasks must run as "
                f"'{self.expected_output}'. Check the `output` field."
            )
        if not self.rules:
            return (
                "The file contains no mapping rules. Check that this is the file that "
                "`identity_mapping_config_path` in config.yaml actually points at."
            )
        return (
            f"No rule matched '{self.identity}'. Check for a typo, an unescaped dot "
            "(write `\\.`), or that you are editing the file "
            "`identity_mapping_config_path` names -- after switching to a privileged "
            "user, that path may still point into a personal home directory."
        )

    def render(self) -> str:
        lines = [f"Mapping file: {self.path}", f"Identity:     {self.identity}", ""]
        for rule in self.rules:
            lines.append(f"  source={rule.source!r}  match={rule.match!r}")
            if rule.compiled:
                lines.append(f"    compiles to: {rule.compiled}")
            if rule.error:
                lines.append(f"    !! {rule.error}")
            if rule.has_anchors:
                lines.append("    !! contains ^ or $ -- these become LITERAL characters.")
                lines.append("       Remove them; matching is already anchored for you.")
            lines.append(f"    output: {rule.output}")
            lines.append("")
        for err in self.errors:
            lines.append(f"  !! {err}")
        if self.errors:
            lines.append("")
        if self.passed:
            lines.append(f"PASS: {self.identity} maps to '{self.expected_output}'.")
            lines.append("")
            lines.append(
                "No restart is needed -- the endpoint polls this file and reloads within "
                "about five seconds. Tell the coordinator and they can re-test right away."
            )
        else:
            lines.append(f"FAIL: {self.identity} does not map to '{self.expected_output}'.")
            lines.append("")
            lines.append(self.diagnosis())
        return "\n".join(lines)


def build_mapping_document(identity_match: str, output_account: str) -> list[dict]:
    """The mapping document a partner should install, as a Python object.

    ``identity_match`` comes from ``Coordinator.identity_regex`` -- already correctly
    escaped and deliberately unanchored.
    """
    return [
        {
            "DATA_TYPE": "expression_identity_mapping#1.0.0",
            "mappings": [
                {
                    "source": "{username}",
                    "match": identity_match,
                    "output": output_account,
                }
            ],
        }
    ]


def validate_mapping_file(
    path: str | Path,
    identity: str,
    expected_output: str,
    identity_id: str | None = None,
) -> MappingReport:
    """Check that ``path`` maps ``identity`` to ``expected_output``.

    Uses the same ``globus-identity-mapping`` library the endpoint uses, so the answer is
    authoritative rather than a reimplementation of its rules.
    """
    path = Path(path)
    report = MappingReport(path=path, identity=identity, expected_output=expected_output)

    try:
        from globus_identity_mapping.loader import load_mappers
    except ImportError:
        report.errors.append(
            "globus-identity-mapping is not installed. It ships with "
            "globus-compute-endpoint, so run this in the same environment as the "
            "endpoint -- typically the conda env named in your worker_init."
        )
        return report

    if not path.is_file():
        report.errors.append(
            f"no such file: {path}\n"
            "     The authoritative path is whatever `identity_mapping_config_path` in "
            "the endpoint's config.yaml points at -- check there first."
        )
        return report

    try:
        document = json.loads(path.read_bytes())
    except json.JSONDecodeError as exc:
        report.errors.append(f"{path} is not valid JSON: {exc}")
        return report

    # The identity record as the Compute service presents it to the endpoint. `status`
    # matters: map_identities silently skips anything not "used" or "private", which
    # looks identical to a bad regex.
    record = [
        {
            "id": identity_id or "00000000-0000-0000-0000-000000000000",
            "sub": identity_id or "00000000-0000-0000-0000-000000000000",
            "username": identity,
            "status": _ACTIVE_STATUS,
        }
    ]

    entries = document if isinstance(document, list) else [document]
    for entry in entries:
        try:
            mappers = load_mappers([entry], None, None)
        except Exception as exc:
            report.errors.append(
                f"could not load mapping document: {type(exc).__name__}: {exc}"
            )
            continue

        for mapper in mappers:
            for mapping in getattr(mapper, "mappings", []):
                compiled = mapping.get("compiled_match")
                report.rules.append(
                    MappingRule(
                        source=mapping.get("source", "?"),
                        match=mapping.get("match", "?"),
                        output=mapping.get("output", "?"),
                        compiled=compiled.pattern if compiled is not None else None,
                    )
                )
            try:
                for result in mapper.map_identities(record):
                    for names in result.values():
                        report.mapped_to.extend(names)
            except Exception as exc:
                report.errors.append(f"mapping raised {type(exc).__name__}: {exc}")

    return report


def validate_for_federation(path: str | Path, federation, experiment: str) -> MappingReport:
    """Validate a mapping file against a federation's coordinator and service account.

    The identity and the expected local account both come from ``federation.yaml``, never
    from a literal -- a different coordinator running this gets the right answer for their
    own federation with no code change.
    """
    exp = federation.experiment(experiment)
    try:
        identity = federation.coordinator.identity
    except AttributeError as exc:  # pragma: no cover - defensive
        raise FederationError("federation config has no coordinator identity") from exc
    return validate_mapping_file(
        path,
        identity=identity,
        expected_output=exp.service_account,
        identity_id=federation.coordinator.identity_id,
    )
