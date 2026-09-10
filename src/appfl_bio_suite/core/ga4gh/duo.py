"""DUO -- the Data Use Ontology, as an access decision a machine can make.

A DUO consent code says what a dataset may be used for. This module turns one into a
decision: given a site's ``DATA_USE.json`` and the study a coordinator declared, is this
site permitted to compute for this run, and if not, which term refused?

THE ONTOLOGY IS VENDORED, NOT FETCHED
-------------------------------------
``data/duo.json`` is a snapshot of the DUO release named in its ``version_iri``,
extracted from the ontology's own OWL. Two reasons it is not resolved over the network:

* A partner's compute node has no outbound web access. The site-side check has to work
  with the standard library and a file on disk, or it does not run where it matters.
* A consent decision that changes because an ontology release moved under a running
  federation is a governance incident, not a feature. The snapshot's version IRI is
  recorded in every decision, so an upgrade is visible in the provenance.

One quirk of that release is preserved deliberately: the IRI for *population origins or
ancestry research prohibited* is ``.../DUO_00000044`` -- one zero too many -- while its
canonical ``id`` is ``DUO:0000044``. Ids are what this module keys on and IRIs are what
it reports, so both are carried and neither is "corrected".

THE THREE VOCABULARIES, WHICH ARE EASY TO CONFLATE
--------------------------------------------------
    permission  what the dataset allows          (subclasses of DUO:0000001)
    modifier    extra conditions on that         (subclasses of DUO:0000017)
    purpose     what a STUDY is trying to do     (DUO:0000031 .. DUO:0000040)

A profile carries a permission and modifiers. A request carries purposes and
attestations. Matching one against the other is the whole job, and the direction matters:
the study must fit inside the permission, not the other way round.

DENY, AND THE THIRD OUTCOME
---------------------------
A decision is ``permitted``, ``denied``, or ``undetermined``. The third is not a hedge --
it is the honest result for terms that carry free text a program cannot evaluate
(*research specific restrictions*, DUO:0000012, whose value might be "no use in studies
of X"). ``undetermined`` blocks dispatch exactly as ``denied`` does, and is cleared only
by a human recording that they read it, in the request's ``acknowledged`` list. Silently
treating an unevaluable restriction as satisfied is the failure mode this exists to
prevent.

THIS MODULE IS PURE
-------------------
:func:`evaluate` takes and returns plain dicts, with no pydantic and no suite imports.
That is what lets ``experiments/fine_mapping/dataset.py`` -- which is shipped to a
partner's worker and may import neither -- carry a faithful copy of it, held to this one
by ``tests/test_ga4gh_shipped_parity.py``.
"""

from __future__ import annotations

import json
from datetime import date
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

__all__ = [
    "DUO_SNAPSHOT_PATH",
    "PERMISSIONS",
    "MODIFIERS",
    "PURPOSES",
    "DataUseProfile",
    "DataUseRequest",
    "Decision",
    "Reason",
    "DuoError",
    "evaluate",
    "load_profile",
    "write_profile",
    "ontology",
    "term",
    "term_label",
    "describe_term",
]

DUO_SNAPSHOT_PATH = Path(__file__).resolve().parent / "data" / "duo.json"

# The five data use permissions of DUO:0000001. Exactly one appears on a profile: a
# dataset that "is GRU and also HMB" is a dataset whose consent code was never agreed.
PERMISSIONS = {
    "DUO:0000004": "no restriction",
    "DUO:0000042": "general research use",
    "DUO:0000006": "health or medical or biomedical research",
    "DUO:0000007": "disease specific research",
    "DUO:0000011": "population origins or ancestry research only",
}

# Research purposes a REQUEST declares. Subclasses of OBI:0000066 in DUO, not of
# DUO:0000001 -- which is the formal statement of the direction above.
PURPOSES = {
    "DUO:0000031": "method development",
    "DUO:0000032": "population research",
    "DUO:0000033": "ancestry research",
    "DUO:0000034": "age category research",
    "DUO:0000035": "gender category research",
    "DUO:0000036": "research control",
    "DUO:0000037": "biomedical research",
    "DUO:0000038": "genetic research",
    "DUO:0000039": "drug development research",
    "DUO:0000040": "disease category research",
}

# Purposes that are biomedical: DUO:0000037 and its three subclasses. HMB is satisfied by
# these and by nothing else.
BIOMEDICAL_PURPOSES = frozenset({"DUO:0000037", "DUO:0000038", "DUO:0000039", "DUO:0000040"})

# Purposes that are population/ancestry work. POA permits only these; DUO:0000044
# prohibits exactly these.
ANCESTRY_PURPOSES = frozenset({"DUO:0000032", "DUO:0000033"})

# Every modifier this matcher knows how to evaluate, and what satisfies it. The values
# are the request keys consulted; a modifier absent from here is reported as
# `undetermined` rather than quietly ignored, which is what keeps an unrecognized term
# from reading as an approval.
MODIFIERS = {
    "DUO:0000012": "research specific restrictions",
    "DUO:0000015": "no general methods research",
    "DUO:0000016": "genetic studies only",
    "DUO:0000018": "not for profit, non commercial use only",
    "DUO:0000019": "publication required",
    "DUO:0000020": "collaboration required",
    "DUO:0000021": "ethics approval required",
    "DUO:0000022": "geographical restriction",
    "DUO:0000024": "publication moratorium",
    "DUO:0000025": "time limit on use",
    "DUO:0000026": "user specific restriction",
    "DUO:0000027": "project specific restriction",
    "DUO:0000028": "institution specific restriction",
    "DUO:0000029": "return to database or resource",
    "DUO:0000043": "clinical care use",
    "DUO:0000044": "population origins or ancestry research prohibited",
    "DUO:0000045": "not for profit organisation use only",
    "DUO:0000046": "non-commercial use only",
}

PERMITTED = "permitted"
DENIED = "denied"
UNDETERMINED = "undetermined"


class DuoError(ValueError):
    """A data use profile or request is malformed."""


# ---------------------------------------------------------------------------
# the vendored ontology
# ---------------------------------------------------------------------------


@lru_cache(maxsize=1)
def ontology() -> dict[str, Any]:
    """The vendored DUO release, as a dict. Cached; the file never changes at run time."""
    try:
        return json.loads(DUO_SNAPSHOT_PATH.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:  # pragma: no cover - packaging failure
        raise DuoError(
            f"the DUO snapshot is missing from the install ({DUO_SNAPSHOT_PATH}). It is "
            "package data; reinstall the suite, or check "
            "[tool.setuptools.package-data] in pyproject.toml."
        ) from exc


def term(term_id: str) -> dict[str, Any]:
    """One DUO term, or raise naming what was expected."""
    terms = ontology()["terms"]
    try:
        return terms[term_id]
    except KeyError:
        raise DuoError(
            f"'{term_id}' is not a DUO term in {ontology()['version_iri']}. Ids look "
            "like 'DUO:0000042'; see `appfl-bio-suite ga4gh duo terms`."
        ) from None


def term_label(term_id: str) -> str:
    return str(term(term_id)["label"])


def describe_term(term_id: str) -> str:
    entry = term(term_id)
    return f"{entry['id']} {entry['label']} -- {entry['definition']}"


# ---------------------------------------------------------------------------
# the documents
# ---------------------------------------------------------------------------


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ModifierTerm(_Strict):
    """One DUO modifier on a dataset, with the value some modifiers carry.

    ``value`` is a list because every value-carrying modifier in DUO is naturally
    plural -- three approved institutions, two permitted countries -- and a schema that
    allows one and then needs two is a schema people work around with commas.
    """

    id: str
    value: list[str] = Field(default_factory=list)

    @field_validator("id")
    @classmethod
    def _known(cls, value: str) -> str:
        value = value.strip()
        term(value)  # raises with the right message if unknown
        if value in PERMISSIONS:
            raise ValueError(
                f"{value} ({term_label(value)}) is a data use PERMISSION, not a modifier. "
                "It belongs in `permission`, and only one permission may be declared."
            )
        return value

    @property
    def label(self) -> str:
        return term_label(self.id)


class DataUseProfile(_Strict):
    """What one dataset permits: the consent code, machine-readable.

    This is what a site ships inside its own bundle as ``DATA_USE.json`` and what its
    worker enforces. The dataset it describes is named by ``drs_uri`` where one exists,
    which is the join between this standard and DRS: a consent code that does not say
    which bytes it governs governs nothing in particular.
    """

    # Free-form identity of the dataset these terms cover.
    dataset_id: str
    # DRS URI of the object these terms govern. Optional only because a profile can be
    # written before the bundle is cut; preflight warns when it stays empty.
    drs_uri: str | None = None
    site: str | None = None
    description: str = ""
    contact: str | None = None
    # Who asserts these terms. A consent code with no steward is a claim with no author.
    steward: str | None = None
    issued: date | None = None

    permission: str
    # The DUO release these term ids were written against, recorded by write_profile.
    # Carried rather than ignored because a decision is only reproducible if the
    # vocabulary it was made in is known -- and because a profile written under a much
    # older release is a thing a steward should be told about rather than a thing that
    # silently still parses.
    duo_version: str | None = None
    # Disease terms for DUO:0000007. Ontology ids (MONDO/DOID/...) or labels; matched
    # case-insensitively against the request's `disease`.
    permission_value: list[str] = Field(default_factory=list)
    modifiers: list[ModifierTerm] = Field(default_factory=list)

    @field_validator("permission")
    @classmethod
    def _is_a_permission(cls, value: str) -> str:
        value = value.strip()
        if value not in PERMISSIONS:
            known = "\n  ".join(f"{k}  {v}" for k, v in PERMISSIONS.items())
            raise ValueError(
                f"'{value}' is not a DUO data use permission. Exactly one of:\n  {known}"
            )
        return value

    @model_validator(mode="after")
    def _disease_specific_names_a_disease(self) -> DataUseProfile:
        if self.permission == "DUO:0000007" and not self.permission_value:
            raise ValueError(
                "DUO:0000007 (disease specific research) requires `permission_value` -- "
                "the disease(s) the data may be used to study. Without it the term "
                "permits nothing in particular and no request can be matched against it."
            )
        return self

    @model_validator(mode="after")
    def _no_duplicate_modifiers(self) -> DataUseProfile:
        seen = [m.id for m in self.modifiers]
        dupes = sorted({m for m in seen if seen.count(m) > 1})
        if dupes:
            raise ValueError(
                f"duplicate modifier(s) {dupes}. Merge their values into one entry: two "
                "entries for the same term with different values is ambiguous about "
                "whether it means 'and' or 'or'."
            )
        return self

    def modifier(self, term_id: str) -> ModifierTerm | None:
        for entry in self.modifiers:
            if entry.id == term_id:
                return entry
        return None

    def to_dict(self) -> dict[str, Any]:
        return json.loads(self.model_dump_json(exclude_none=True))

    def render(self) -> str:
        lines = [
            f"dataset  {self.dataset_id}" + (f"  ({self.site})" if self.site else ""),
            f"terms    {self.permission}  {term_label(self.permission)}",
        ]
        if self.permission_value:
            lines.append(f"         for: {', '.join(self.permission_value)}")
        for entry in self.modifiers:
            suffix = f": {', '.join(entry.value)}" if entry.value else ""
            lines.append(f"         + {entry.id}  {entry.label}{suffix}")
        if self.drs_uri:
            lines.append(f"governs  {self.drs_uri}")
        if self.steward:
            lines.append(f"steward  {self.steward}")
        return "\n".join(lines)


class DataUseRequest(_Strict):
    """What this study intends to do, and what its requester attests to.

    Written once by the coordinator in ``federation.yaml`` and evaluated against every
    site's profile. It is deliberately not per-site: a study has one purpose, and a
    coordinator who can vary the declared purpose per site to get past a refusal has a
    gate that gates nothing.
    """

    # Who is asking. `requester` is the coordinator's Globus identity -- the same string
    # partners already authorize in their identity mapping, which is what makes
    # DUO:0000026 (user specific restriction) checkable rather than aspirational.
    #
    # Defaulted rather than required so that `federation.yaml` need not repeat the
    # identity it already declares; `Federation.data_use_request()` fills it in. A
    # request that reaches a site with this still empty fails DUO:0000026 and is
    # reported as such, which is the right outcome -- an anonymous requester is not an
    # approved user.
    requester: str = ""
    project: str | None = None
    institution: str | None = None
    # ISO 3166 country/region code or free text, matched case-insensitively against
    # DUO:0000022 values.
    geography: str | None = None

    purposes: list[str] = Field(default_factory=list)
    # For matching DUO:0000007 (disease specific research).
    disease: str | None = None

    # Attestations. Each corresponds to a modifier that a requester can satisfy by
    # agreeing to something rather than by being something.
    non_commercial: bool = False
    not_for_profit_organisation: bool = False
    ethics_approval: str | None = None
    collaboration_agreed: bool = False
    publication_agreed: bool = False
    moratorium_accepted: bool = False
    return_to_resource_agreed: bool = False

    # Terms a human read and accepted, by id. The ONLY way to clear `undetermined`.
    acknowledged: list[str] = Field(default_factory=list)

    @field_validator("purposes")
    @classmethod
    def _known_purposes(cls, value: list[str]) -> list[str]:
        for entry in value:
            if entry not in PURPOSES:
                known = "\n  ".join(f"{k}  {v}" for k, v in PURPOSES.items())
                raise ValueError(
                    f"'{entry}' is not a DUO research purpose. A request declares what "
                    f"the STUDY does, from:\n  {known}\n"
                    "(Permissions like DUO:0000042 describe a dataset, not a study.)"
                )
        return value

    @model_validator(mode="after")
    def _purposes_present(self) -> DataUseRequest:
        if not self.purposes:
            raise ValueError(
                "a data use request must declare at least one research purpose. A "
                "request with no purpose cannot be matched against any permission, and "
                "defaulting it to 'general' would make the whole check ceremonial."
            )
        return self

    def to_dict(self) -> dict[str, Any]:
        return json.loads(self.model_dump_json(exclude_none=True))

    def render(self) -> str:
        lines = [f"requester {self.requester}"]
        if self.project:
            lines.append(f"project   {self.project}")
        if self.institution:
            lines.append(f"institute {self.institution}")
        lines.append("purposes  " + ", ".join(f"{p} {PURPOSES[p]}" for p in self.purposes))
        if self.disease:
            lines.append(f"disease   {self.disease}")
        attested = [
            name
            for name, value in (
                ("non-commercial", self.non_commercial),
                ("not-for-profit organisation", self.not_for_profit_organisation),
                ("collaboration", self.collaboration_agreed),
                ("publication", self.publication_agreed),
                ("moratorium", self.moratorium_accepted),
                ("return to resource", self.return_to_resource_agreed),
            )
            if value
        ]
        if self.ethics_approval:
            attested.append(f"ethics approval {self.ethics_approval}")
        if attested:
            lines.append("attests   " + ", ".join(attested))
        if self.acknowledged:
            lines.append("acknowledged " + ", ".join(self.acknowledged))
        return "\n".join(lines)


class Reason(_Strict):
    """One term's contribution to a decision."""

    term_id: str
    term_label: str
    outcome: Literal["permitted", "denied", "undetermined"]
    detail: str = ""

    def render(self) -> str:
        mark = {PERMITTED: "ok", DENIED: "DENY", UNDETERMINED: "????"}[self.outcome]
        return f"  [{mark:>4}] {self.term_id} {self.term_label}: {self.detail}"


class Decision(_Strict):
    """The answer, with every term that produced it."""

    outcome: Literal["permitted", "denied", "undetermined"]
    dataset_id: str
    requester: str
    reasons: list[Reason] = Field(default_factory=list)
    evaluated_at: str = ""
    duo_version: str = ""

    @property
    def permitted(self) -> bool:
        return self.outcome == PERMITTED

    @property
    def blocking(self) -> list[Reason]:
        return [r for r in self.reasons if r.outcome != PERMITTED]

    def render(self) -> str:
        head = f"{self.outcome.upper()}: {self.dataset_id} for {self.requester}"
        return "\n".join([head, *(r.render() for r in self.reasons)])

    def to_dict(self) -> dict[str, Any]:
        return json.loads(self.model_dump_json())


# ---------------------------------------------------------------------------
# the matcher
# ---------------------------------------------------------------------------


def evaluate(
    profile: dict[str, Any],
    request: dict[str, Any],
    as_of: str | None = None,
) -> dict[str, Any]:
    """Match one profile against one request. Plain dicts in, plain dict out.

    Kept free of pydantic and of every suite import so that the shipped site loader can
    carry a faithful copy (see this module's docstring). :func:`decide` is the typed
    wrapper most callers want.

    ``as_of`` is an ISO date; it defaults to today and exists so that DUO:0000025 (time
    limit on use) is testable and so a decision recorded in provenance can be re-derived.
    """
    today = as_of or date.today().isoformat()
    reasons: list[dict[str, Any]] = []

    permission = str(profile.get("permission", "")).strip()
    purposes = [str(p) for p in request.get("purposes", [])]
    acknowledged = {str(a) for a in request.get("acknowledged", [])}

    def add(term_id: str, outcome: str, detail: str) -> None:
        labels = {**PERMISSIONS, **MODIFIERS, **PURPOSES}
        reasons.append(
            {
                "term_id": term_id,
                "term_label": labels.get(term_id, ""),
                "outcome": outcome,
                "detail": detail,
            }
        )

    # -- the permission ---------------------------------------------------
    if permission not in PERMISSIONS:
        add(permission or "(none)", DENIED, "not a DUO data use permission")
    elif permission == "DUO:0000004":
        add(permission, PERMITTED, "no restriction on use")
    elif permission == "DUO:0000042":
        add(permission, PERMITTED, "general research use permits any research purpose")
    elif permission == "DUO:0000006":
        outside = [p for p in purposes if p not in BIOMEDICAL_PURPOSES]
        if outside:
            add(
                permission,
                DENIED,
                "permits health/medical/biomedical research only, but the study declares "
                + ", ".join(f"{p} ({PURPOSES.get(p, '?')})" for p in outside),
            )
        else:
            add(permission, PERMITTED, "every declared purpose is biomedical research")
    elif permission == "DUO:0000007":
        allowed = [str(v).strip().lower() for v in profile.get("permission_value", [])]
        disease = str(request.get("disease") or "").strip().lower()
        if not disease:
            add(
                permission,
                DENIED,
                "disease specific research: the study declares no disease. Permitted "
                "for: " + (", ".join(allowed) or "(unspecified)"),
            )
        elif disease not in allowed:
            add(
                permission,
                DENIED,
                f"disease specific research for {', '.join(allowed) or '(unspecified)'}; "
                f"the study declares '{disease}'",
            )
        else:
            add(permission, PERMITTED, f"disease specific research, and the study is on {disease}")
    elif permission == "DUO:0000011":
        outside = [p for p in purposes if p not in ANCESTRY_PURPOSES]
        if outside:
            add(
                permission,
                DENIED,
                "permits population origins or ancestry research only, but the study "
                "declares " + ", ".join(f"{p} ({PURPOSES.get(p, '?')})" for p in outside),
            )
        else:
            add(permission, PERMITTED, "the study is population/ancestry research")

    # -- the modifiers ----------------------------------------------------
    for entry in profile.get("modifiers", []):
        term_id = str(entry.get("id", "")).strip()
        values = [str(v) for v in entry.get("value", [])]
        lowered = [v.strip().lower() for v in values]
        outcome, detail = _evaluate_modifier(term_id, values, lowered, request, purposes, today)
        if outcome == UNDETERMINED and term_id in acknowledged:
            outcome = PERMITTED
            detail = f"acknowledged by the requester: {detail}"
        add(term_id, outcome, detail)

    outcomes = {r["outcome"] for r in reasons}
    if DENIED in outcomes:
        overall = DENIED
    elif UNDETERMINED in outcomes:
        overall = UNDETERMINED
    else:
        overall = PERMITTED

    return {
        "outcome": overall,
        "dataset_id": str(profile.get("dataset_id", "")),
        "requester": str(request.get("requester", "")),
        "reasons": reasons,
        "evaluated_at": today,
        "duo_version": str(profile.get("duo_version", "")) or _snapshot_version(),
    }


def _snapshot_version() -> str:
    """The vendored release IRI, or "" where the snapshot is unreadable.

    Never raises: a decision that is otherwise complete should not fail because the
    version string could not be read, and the empty value is visible in the provenance.
    """
    try:
        return str(ontology().get("version_iri", ""))
    except Exception:  # noqa: BLE001 - provenance detail, never a reason to refuse
        return ""


def _evaluate_modifier(
    term_id: str,
    values: list[str],
    lowered: list[str],
    request: dict[str, Any],
    purposes: list[str],
    today: str,
) -> tuple[str, str]:
    """One modifier against the request. Returns ``(outcome, detail)``.

    Split out so both this module and its shipped copy read as one table of rules rather
    than a long branch inside a loop.
    """
    ask = request.get

    if term_id == "DUO:0000015":  # no general methods research
        if "DUO:0000031" in purposes:
            return DENIED, "the study declares method development (DUO:0000031)"
        return PERMITTED, "the study declares no methods development"

    if term_id == "DUO:0000016":  # genetic studies only
        outside = [p for p in purposes if p != "DUO:0000038"]
        if outside:
            return DENIED, "genetic studies only; the study also declares " + ", ".join(outside)
        return PERMITTED, "the study is genetic research"

    if term_id == "DUO:0000044":  # population origins or ancestry research prohibited
        offending = [p for p in purposes if p in ANCESTRY_PURPOSES]
        if offending:
            return DENIED, "population/ancestry research is prohibited: " + ", ".join(offending)
        return PERMITTED, "the study declares no population or ancestry research"

    if term_id in ("DUO:0000018", "DUO:0000046"):  # non-commercial (and not-for-profit)
        if not ask("non_commercial", False):
            return DENIED, "the requester does not attest to non-commercial use"
        if term_id == "DUO:0000018" and not ask("not_for_profit_organisation", False):
            return DENIED, "the requester does not attest to being a not-for-profit organisation"
        return PERMITTED, "attested"

    if term_id == "DUO:0000045":  # not for profit organisation use only
        if not ask("not_for_profit_organisation", False):
            return DENIED, "the requester does not attest to being a not-for-profit organisation"
        return PERMITTED, "attested"

    if term_id == "DUO:0000019":  # publication required
        if not ask("publication_agreed", False):
            return DENIED, "the requester has not agreed to publish results"
        return PERMITTED, "attested"

    if term_id == "DUO:0000020":  # collaboration required
        if not ask("collaboration_agreed", False):
            return DENIED, "the requester has not agreed to collaborate" + (
                f" with {', '.join(values)}" if values else ""
            )
        return PERMITTED, "attested" + (f" ({', '.join(values)})" if values else "")

    if term_id == "DUO:0000021":  # ethics approval required
        approval = str(ask("ethics_approval") or "").strip()
        if not approval:
            return DENIED, "no ethics/IRB approval reference was supplied"
        return PERMITTED, f"approval {approval}"

    if term_id == "DUO:0000022":  # geographical restriction
        where = str(ask("geography") or "").strip().lower()
        if not values:
            return UNDETERMINED, "restricted geographically, but no region was specified"
        if not where:
            return DENIED, "restricted to " + ", ".join(values) + "; the request names no region"
        if where not in lowered:
            return DENIED, f"restricted to {', '.join(values)}; the request names '{where}'"
        return PERMITTED, f"within {where}"

    if term_id == "DUO:0000024":  # publication moratorium
        if not ask("moratorium_accepted", False):
            return DENIED, "the requester has not accepted the publication moratorium" + (
                f" (until {values[0]})" if values else ""
            )
        return PERMITTED, "accepted" + (f" until {values[0]}" if values else "")

    if term_id == "DUO:0000025":  # time limit on use
        if not values:
            return UNDETERMINED, "a time limit applies, but no date was specified"
        limit = values[0].strip()
        if today > limit:
            return DENIED, f"use expired on {limit} (evaluating as of {today})"
        return PERMITTED, f"within the limit, which ends {limit}"

    if term_id == "DUO:0000026":  # user specific restriction
        who = str(ask("requester") or "").strip().lower()
        if not values:
            return UNDETERMINED, "restricted to approved users, but none were listed"
        if who not in lowered:
            return DENIED, f"'{who}' is not among the approved users"
        return PERMITTED, f"'{who}' is an approved user"

    if term_id == "DUO:0000027":  # project specific restriction
        project = str(ask("project") or "").strip().lower()
        if not values:
            return UNDETERMINED, "restricted to approved projects, but none were listed"
        if project not in lowered:
            return DENIED, (
                f"restricted to {', '.join(values)}; the request names "
                f"'{project or '(no project)'}'"
            )
        return PERMITTED, f"project '{project}' is approved"

    if term_id == "DUO:0000028":  # institution specific restriction
        institution = str(ask("institution") or "").strip().lower()
        if not values:
            return UNDETERMINED, "restricted to approved institutions, but none were listed"
        if institution not in lowered:
            return DENIED, (
                f"restricted to {', '.join(values)}; the request names "
                f"'{institution or '(no institution)'}'"
            )
        return PERMITTED, f"institution '{institution}' is approved"

    if term_id == "DUO:0000029":  # return to database or resource
        if not ask("return_to_resource_agreed", False):
            return DENIED, "the requester has not agreed to return derived data to the resource"
        return PERMITTED, "attested"

    if term_id == "DUO:0000043":  # clinical care use
        # A permissive modifier: it widens what the data may be used for rather than
        # narrowing it, so a research study is unaffected by its presence.
        return PERMITTED, "clinical care use is additionally permitted; no effect on research use"

    if term_id == "DUO:0000012":  # research specific restrictions
        return UNDETERMINED, (
            "free-text research restriction"
            + (f": {'; '.join(values)}" if values else "")
            + ". A person must read this and record the decision in `acknowledged`."
        )

    return UNDETERMINED, (
        "this matcher has no rule for that term, so it cannot be treated as satisfied"
    )


def decide(profile: DataUseProfile, request: DataUseRequest, as_of: str | None = None) -> Decision:
    """Typed wrapper over :func:`evaluate`."""
    return Decision.model_validate(evaluate(profile.to_dict(), request.to_dict(), as_of))


# ---------------------------------------------------------------------------
# profiles on disk
# ---------------------------------------------------------------------------

PROFILE_FILENAME = "DATA_USE.json"


def load_profile(path: str | Path) -> DataUseProfile:
    """Read a ``DATA_USE.json``, with the message naming the file when it is wrong."""
    path = Path(path)
    if path.is_dir():
        path = path / PROFILE_FILENAME
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise DuoError(
            f"no data use profile at {path}. Every dataset in a federation declares its "
            f"terms in {PROFILE_FILENAME}; the simulation writes one into each bundle, "
            "and a real site writes its own."
        ) from exc
    except json.JSONDecodeError as exc:
        raise DuoError(f"{path} is not valid JSON: {exc}") from exc
    try:
        return DataUseProfile.model_validate(raw)
    except Exception as exc:
        raise DuoError(f"{path} is not a valid data use profile:\n\n{exc}") from exc


def write_profile(profile: DataUseProfile, path: str | Path) -> Path:
    """Write a profile as ``DATA_USE.json``, creating parents."""
    path = Path(path)
    if path.is_dir():
        path = path / PROFILE_FILENAME
    path.parent.mkdir(parents=True, exist_ok=True)
    if not profile.duo_version:
        profile = profile.model_copy(update={"duo_version": _snapshot_version()})
    path.write_text(json.dumps(profile.to_dict(), indent=2) + "\n", encoding="utf-8")
    return path
