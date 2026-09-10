"""DUO: the matcher, the vendored ontology, and the decisions that gate a run.

The cases here are the ones a federation actually hits. Each is a shape a real consent
code takes, not a synthetic permutation: general research use, biomedical-only, disease
specific, ancestry-only, and the modifiers a requester satisfies by attesting to
something versus the ones only a person can clear.
"""

from __future__ import annotations

import json

import pytest

from appfl_bio_suite.core.ga4gh.duo import (
    DENIED,
    MODIFIERS,
    PERMISSIONS,
    PERMITTED,
    PURPOSES,
    UNDETERMINED,
    DataUseProfile,
    DataUseRequest,
    DuoError,
    decide,
    evaluate,
    load_profile,
    ontology,
    write_profile,
)

# ---------------------------------------------------------------------------
# the vendored ontology
# ---------------------------------------------------------------------------


def test_snapshot_is_a_real_duo_release():
    """The snapshot must name the release it came from, or a decision is unattributable."""
    snapshot = ontology()
    assert snapshot["ontology"] == "DUO"
    assert snapshot["version_iri"].startswith("http://purl.obolibrary.org/obo/duo/releases/")
    assert len(snapshot["terms"]) >= 35


@pytest.mark.parametrize("term_id", sorted({*PERMISSIONS, *MODIFIERS, *PURPOSES}))
def test_every_term_this_matcher_uses_exists_in_the_release(term_id: str):
    """A rule keyed on an id the ontology does not define is a rule that never fires."""
    terms = ontology()["terms"]
    assert term_id in terms, f"{term_id} is not in the vendored DUO release"
    assert not terms[term_id]["deprecated"], f"{term_id} is deprecated upstream"


@pytest.mark.parametrize("term_id,label", sorted(PERMISSIONS.items()))
def test_permission_labels_match_the_ontology(term_id: str, label: str):
    assert ontology()["terms"][term_id]["label"] == label


@pytest.mark.parametrize("term_id,label", sorted(MODIFIERS.items()))
def test_modifier_labels_match_the_ontology(term_id: str, label: str):
    assert ontology()["terms"][term_id]["label"] == label


def test_permissions_and_modifiers_are_disjoint():
    """Confusing the two is the mistake this vocabulary invites; keep it impossible."""
    assert not set(PERMISSIONS) & set(MODIFIERS)
    assert not set(PERMISSIONS) & set(PURPOSES)


def test_the_ancestry_id_keeps_the_releases_own_iri_quirk():
    """DUO:0000044's IRI carries an extra zero upstream. Preserved, not 'corrected'."""
    entry = ontology()["terms"]["DUO:0000044"]
    assert entry["iri"].endswith("DUO_00000044")
    assert entry["id"] == "DUO:0000044"


# ---------------------------------------------------------------------------
# the matcher
# ---------------------------------------------------------------------------


def profile(permission="DUO:0000042", modifiers=(), values=()):
    return {
        "dataset_id": "test/site",
        "permission": permission,
        "permission_value": list(values),
        "modifiers": [dict(m) for m in modifiers],
    }


def request(**overrides):
    base = {"requester": "me@example.org", "purposes": ["DUO:0000038"]}
    base.update(overrides)
    return base


def test_general_research_use_permits_any_purpose():
    for purpose in PURPOSES:
        decision = evaluate(profile(), request(purposes=[purpose]))
        assert decision["outcome"] == PERMITTED, purpose


def test_no_restriction_permits_everything():
    decision = evaluate(profile("DUO:0000004"), request(purposes=list(PURPOSES)))
    assert decision["outcome"] == PERMITTED


def test_biomedical_only_refuses_a_population_origins_study():
    """The case the shipped scenario is built around: HMB versus ancestry research."""
    decision = evaluate(profile("DUO:0000006"), request(purposes=["DUO:0000032"]))
    assert decision["outcome"] == DENIED
    assert "biomedical" in decision["reasons"][0]["detail"]


def test_biomedical_only_permits_genetic_research():
    assert evaluate(profile("DUO:0000006"), request())["outcome"] == PERMITTED


def test_a_mixed_study_fails_the_narrower_permission():
    """Every declared purpose must fit. One that does not is a refusal, not a warning."""
    decision = evaluate(profile("DUO:0000006"), request(purposes=["DUO:0000038", "DUO:0000032"]))
    assert decision["outcome"] == DENIED


def test_disease_specific_needs_a_matching_disease():
    terms = profile("DUO:0000007", values=["type 2 diabetes"])
    assert evaluate(terms, request(disease="Type 2 Diabetes"))["outcome"] == PERMITTED
    assert evaluate(terms, request(disease="asthma"))["outcome"] == DENIED
    assert evaluate(terms, request())["outcome"] == DENIED


def test_population_origins_only_permits_only_that():
    terms = profile("DUO:0000011")
    assert evaluate(terms, request(purposes=["DUO:0000033"]))["outcome"] == PERMITTED
    assert evaluate(terms, request(purposes=["DUO:0000038"]))["outcome"] == DENIED


def test_ancestry_prohibited_is_the_mirror_image():
    terms = profile(modifiers=[{"id": "DUO:0000044"}])
    assert evaluate(terms, request(purposes=["DUO:0000032"]))["outcome"] == DENIED
    assert evaluate(terms, request(purposes=["DUO:0000038"]))["outcome"] == PERMITTED


def test_no_methods_research():
    terms = profile(modifiers=[{"id": "DUO:0000015"}])
    assert evaluate(terms, request(purposes=["DUO:0000031"]))["outcome"] == DENIED
    assert evaluate(terms, request())["outcome"] == PERMITTED


@pytest.mark.parametrize(
    "term_id,attestation",
    [
        ("DUO:0000019", "publication_agreed"),
        ("DUO:0000020", "collaboration_agreed"),
        ("DUO:0000024", "moratorium_accepted"),
        ("DUO:0000029", "return_to_resource_agreed"),
        ("DUO:0000046", "non_commercial"),
        ("DUO:0000045", "not_for_profit_organisation"),
    ],
)
def test_attestation_modifiers_deny_until_attested(term_id: str, attestation: str):
    terms = profile(modifiers=[{"id": term_id}])
    assert evaluate(terms, request())["outcome"] == DENIED
    assert evaluate(terms, request(**{attestation: True}))["outcome"] == PERMITTED


def test_not_for_profit_non_commercial_needs_both():
    """DUO:0000018 is the conjunction, and half of it is not enough."""
    terms = profile(modifiers=[{"id": "DUO:0000018"}])
    assert evaluate(terms, request(non_commercial=True))["outcome"] == DENIED
    assert (
        evaluate(terms, request(non_commercial=True, not_for_profit_organisation=True))["outcome"]
        == PERMITTED
    )


def test_ethics_approval_requires_a_reference_not_a_boolean():
    terms = profile(modifiers=[{"id": "DUO:0000021"}])
    assert evaluate(terms, request())["outcome"] == DENIED
    assert evaluate(terms, request(ethics_approval="IRB-2026-114"))["outcome"] == PERMITTED


def test_user_specific_restriction_matches_the_requesting_identity():
    """The identity partners already authorize is the one checked here."""
    terms = profile(modifiers=[{"id": "DUO:0000026", "value": ["me@example.org"]}])
    assert evaluate(terms, request())["outcome"] == PERMITTED
    assert evaluate(terms, request(requester="someone@else.org"))["outcome"] == DENIED


def test_an_empty_requester_is_not_an_approved_user():
    terms = profile(modifiers=[{"id": "DUO:0000026", "value": ["me@example.org"]}])
    assert evaluate(terms, request(requester=""))["outcome"] == DENIED


@pytest.mark.parametrize(
    "term_id,field,ok,bad",
    [
        ("DUO:0000027", "project", "fedfm", "other"),
        ("DUO:0000028", "institution", "example university", "elsewhere"),
        ("DUO:0000022", "geography", "us", "cn"),
    ],
)
def test_value_matching_is_case_insensitive(term_id, field, ok, bad):
    terms = profile(modifiers=[{"id": term_id, "value": [ok.upper()]}])
    assert evaluate(terms, request(**{field: ok}))["outcome"] == PERMITTED
    assert evaluate(terms, request(**{field: bad}))["outcome"] == DENIED


def test_a_value_modifier_with_no_values_is_undetermined_not_permitted():
    """The failure mode this whole third outcome exists to prevent."""
    terms = profile(modifiers=[{"id": "DUO:0000028"}])
    assert evaluate(terms, request(institution="anywhere"))["outcome"] == UNDETERMINED


def test_time_limit_expires():
    terms = profile(modifiers=[{"id": "DUO:0000025", "value": ["2026-12-31"]}])
    assert evaluate(terms, request(), as_of="2026-06-01")["outcome"] == PERMITTED
    assert evaluate(terms, request(), as_of="2027-01-01")["outcome"] == DENIED


def test_free_text_restrictions_are_undetermined_until_acknowledged():
    terms = profile(modifiers=[{"id": "DUO:0000012", "value": ["no use in studies of X"]}])
    assert evaluate(terms, request())["outcome"] == UNDETERMINED
    cleared = evaluate(terms, request(acknowledged=["DUO:0000012"]))
    assert cleared["outcome"] == PERMITTED
    assert "acknowledged" in cleared["reasons"][-1]["detail"]


def test_an_unknown_modifier_is_undetermined_rather_than_ignored():
    """Silently passing a term nobody implemented is how a gate stops gating."""
    terms = profile(modifiers=[{"id": "DUO:0000043"}, {"id": "DUO:0000001"}])
    assert evaluate(terms, request())["outcome"] == UNDETERMINED


def test_clinical_care_use_does_not_restrict_research():
    terms = profile(modifiers=[{"id": "DUO:0000043"}])
    assert evaluate(terms, request())["outcome"] == PERMITTED


def test_acknowledging_cannot_clear_a_denial():
    """`acknowledged` is for terms a program cannot evaluate, not for overriding one."""
    terms = profile(modifiers=[{"id": "DUO:0000021"}])
    assert evaluate(terms, request(acknowledged=["DUO:0000021"]))["outcome"] == DENIED


def test_a_decision_records_the_ontology_release_it_was_made_under():
    decision = evaluate(profile(), request())
    assert "duo/releases/" in decision["duo_version"]


# ---------------------------------------------------------------------------
# the documents
# ---------------------------------------------------------------------------


def test_a_permission_in_the_modifier_slot_is_refused_by_name():
    with pytest.raises(ValueError, match="PERMISSION, not a modifier"):
        DataUseProfile(dataset_id="x", permission="DUO:0000042", modifiers=[{"id": "DUO:0000006"}])


def test_disease_specific_without_a_disease_is_refused():
    with pytest.raises(ValueError, match="requires `permission_value`"):
        DataUseProfile(dataset_id="x", permission="DUO:0000007")


def test_duplicate_modifiers_are_refused():
    with pytest.raises(ValueError, match="duplicate modifier"):
        DataUseProfile(
            dataset_id="x",
            permission="DUO:0000042",
            modifiers=[{"id": "DUO:0000019"}, {"id": "DUO:0000019"}],
        )


def test_a_request_must_declare_a_purpose():
    with pytest.raises(ValueError, match="at least one research purpose"):
        DataUseRequest(requester="me@example.org")


def test_a_permission_is_not_a_valid_purpose():
    """The direction of the vocabulary, enforced where someone would get it wrong."""
    with pytest.raises(ValueError, match="not a DUO research purpose"):
        DataUseRequest(requester="me@example.org", purposes=["DUO:0000042"])


def test_profile_round_trips_through_disk(tmp_path):
    original = DataUseProfile(
        dataset_id="fine-mapping/anl",
        permission="DUO:0000006",
        modifiers=[{"id": "DUO:0000018"}, {"id": "DUO:0000026", "value": ["me@example.org"]}],
    )
    path = write_profile(original, tmp_path)
    reloaded = load_profile(path)
    assert reloaded.permission == original.permission
    assert [m.id for m in reloaded.modifiers] == [m.id for m in original.modifiers]
    # The writer stamps the release; a profile that does not know its vocabulary cannot
    # be re-evaluated later with confidence.
    assert "duo/releases/" in (reloaded.duo_version or "")


def test_load_profile_names_the_file_when_it_is_malformed(tmp_path):
    path = tmp_path / "DATA_USE.json"
    path.write_text(json.dumps({"dataset_id": "x", "permission": "not-a-term"}), encoding="utf-8")
    with pytest.raises(DuoError) as exc:
        load_profile(path)
    assert str(path) in str(exc.value)


def test_decide_wraps_evaluate_identically():
    terms = DataUseProfile(dataset_id="x", permission="DUO:0000006")
    ask = DataUseRequest(requester="me@example.org", purposes=["DUO:0000038"])
    assert decide(terms, ask).outcome == evaluate(terms.to_dict(), ask.to_dict())["outcome"]
