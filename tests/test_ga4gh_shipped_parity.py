"""The DUO matcher exists twice. This is what keeps the two copies the same.

``core/ga4gh/duo.py`` is the coordinator's; ``experiments/fine_mapping/dataset.py``
carries an inlined copy, because its source is shipped to a worker that has neither this
suite nor any way to import it. The arrangement is the same one ``trainer.py`` and
``fedfm/fed_fine_mapping.py`` are in for the site aggregates, and it is acceptable for
the same reason: the duplication is *tested*, not trusted.

WHY THIS MATTERS MORE THAN THE USUAL DUPLICATION
------------------------------------------------
A drift between two copies of a statistic shows up as a number that disagrees. A drift
between two copies of a consent matcher shows up as a coordinator being told a site will
accept a study that the site then refuses -- or, far worse, the reverse: a preflight that
says 'denied' while the shipped copy says 'permitted' and computes anyway.

So the case table below is the full cross product of the shapes that decide anything, and
both implementations must return the same outcome AND the same per-term reasons.
"""

from __future__ import annotations

import itertools

import pytest

from appfl_bio_suite.core.ga4gh import duo as core_duo
from appfl_bio_suite.experiments.fine_mapping import dataset as shipped

PROFILES = [
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": []},
    {"dataset_id": "d", "permission": "DUO:0000004", "modifiers": []},
    {"dataset_id": "d", "permission": "DUO:0000006", "modifiers": []},
    {
        "dataset_id": "d",
        "permission": "DUO:0000007",
        "permission_value": ["type 2 diabetes"],
        "modifiers": [],
    },
    {"dataset_id": "d", "permission": "DUO:0000011", "modifiers": []},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000015"}]},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000016"}]},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000018"}]},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000019"}]},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000020"}]},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000021"}]},
    {
        "dataset_id": "d",
        "permission": "DUO:0000042",
        "modifiers": [{"id": "DUO:0000022", "value": ["US", "CA"]}],
    },
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000022"}]},
    {
        "dataset_id": "d",
        "permission": "DUO:0000042",
        "modifiers": [{"id": "DUO:0000024", "value": ["2027-01-01"]}],
    },
    {
        "dataset_id": "d",
        "permission": "DUO:0000042",
        "modifiers": [{"id": "DUO:0000025", "value": ["2026-06-30"]}],
    },
    {
        "dataset_id": "d",
        "permission": "DUO:0000042",
        "modifiers": [{"id": "DUO:0000026", "value": ["me@example.org"]}],
    },
    {
        "dataset_id": "d",
        "permission": "DUO:0000042",
        "modifiers": [{"id": "DUO:0000027", "value": ["fedfm"]}],
    },
    {
        "dataset_id": "d",
        "permission": "DUO:0000042",
        "modifiers": [{"id": "DUO:0000028", "value": ["Example University"]}],
    },
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000029"}]},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000043"}]},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000044"}]},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000045"}]},
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000046"}]},
    {
        "dataset_id": "d",
        "permission": "DUO:0000042",
        "modifiers": [{"id": "DUO:0000012", "value": ["ask first"]}],
    },
    # Several modifiers at once, which is what a real consent code looks like.
    {
        "dataset_id": "d",
        "permission": "DUO:0000006",
        "modifiers": [
            {"id": "DUO:0000018"},
            {"id": "DUO:0000021"},
            {"id": "DUO:0000026", "value": ["me@example.org"]},
        ],
    },
    # A term neither implementation has a rule for.
    {"dataset_id": "d", "permission": "DUO:0000042", "modifiers": [{"id": "DUO:0000001"}]},
    # A malformed permission.
    {"dataset_id": "d", "permission": "not-a-term", "modifiers": []},
]

REQUESTS = [
    {"requester": "me@example.org", "purposes": ["DUO:0000038"]},
    {"requester": "SOMEONE@else.org", "purposes": ["DUO:0000032", "DUO:0000033"]},
    {"requester": "me@example.org", "purposes": ["DUO:0000031"]},
    {
        "requester": "me@example.org",
        "purposes": ["DUO:0000038"],
        "disease": "Type 2 Diabetes",
        "project": "FedFM",
        "institution": "example university",
        "geography": "us",
        "non_commercial": True,
        "not_for_profit_organisation": True,
        "ethics_approval": "IRB-1",
        "collaboration_agreed": True,
        "publication_agreed": True,
        "moratorium_accepted": True,
        "return_to_resource_agreed": True,
    },
    {"requester": "", "purposes": ["DUO:0000038"], "acknowledged": ["DUO:0000012"]},
]

# Fixed, so a case that depends on the date (DUO:0000025) is decided identically on both
# sides and does not change meaning next year.
AS_OF = "2026-06-01"

CASES = list(itertools.product(range(len(PROFILES)), range(len(REQUESTS))))


@pytest.mark.parametrize(("profile_index", "request_index"), CASES)
def test_the_two_matchers_agree(profile_index: int, request_index: int):
    profile = PROFILES[profile_index]
    request = REQUESTS[request_index]

    mine = core_duo.evaluate(profile, request, as_of=AS_OF)
    theirs = shipped.evaluate_data_use(profile, request, as_of=AS_OF)

    assert mine["outcome"] == theirs["outcome"], (
        f"outcome differs for profile {profile_index} / request {request_index}:\n"
        f"  core/ga4gh/duo.py: {mine['outcome']}\n"
        f"  dataset.py:        {theirs['outcome']}\n"
        "The shipped copy of the matcher has drifted from the coordinator's. Both "
        "decide whether a site's data may be used; a difference means the preflight "
        "and the site can disagree about that."
    )
    assert [(r["term_id"], r["outcome"]) for r in mine["reasons"]] == [
        (r["term_id"], r["outcome"]) for r in theirs["reasons"]
    ]
    assert [r["detail"] for r in mine["reasons"]] == [r["detail"] for r in theirs["reasons"]]


def test_the_two_label_tables_agree_where_they_overlap():
    """The shipped copy carries its own labels; a stale one misreports a refusal."""
    core_labels = {**core_duo.PERMISSIONS, **core_duo.MODIFIERS, **core_duo.PURPOSES}
    for term_id, label in core_labels.items():
        assert shipped.DUO_LABELS.get(term_id) == label, term_id


def test_the_shipped_copy_knows_the_same_permission_set():
    assert shipped.PERMISSION_IDS == frozenset(core_duo.PERMISSIONS)


def test_the_shipped_copy_uses_the_same_purpose_groupings():
    assert shipped.BIOMEDICAL_PURPOSES == core_duo.BIOMEDICAL_PURPOSES
    assert shipped.ANCESTRY_PURPOSES == core_duo.ANCESTRY_PURPOSES


def test_the_shipped_outcome_strings_are_the_same_three():
    assert (shipped.PERMITTED, shipped.DENIED, shipped.UNDETERMINED) == (
        core_duo.PERMITTED,
        core_duo.DENIED,
        core_duo.UNDETERMINED,
    )
