"""The smoke probe must not require a partner to have installed this package.

An experiment's dataset and model are shipped to a worker as source text precisely so a
partner need not install `appfl_bio_suite` -- see experiments/*/dataset.py. A readiness
check that is stricter than the run it clears is a check that fails on correctly
configured endpoints, which is the one thing it must never do.
"""

from __future__ import annotations

import base64

import pytest

from appfl_bio_suite.core.endpoint import explain_submit_failure, smoke_test, where_am_i


def _decoded(payload: str) -> bytes:
    """Everything in a serialized payload, base64 chunks decoded where they decode.

    The wire format is an strategy-id header followed by base64 lines, and a combined
    payload carries several. Searching the decoded bytes rather than matching the format
    keeps this test about what reaches the worker, not about how the SDK frames it.
    """
    out = [payload.encode()]
    for chunk in payload.replace(":", "\n").split("\n"):
        try:
            out.append(base64.b64decode(chunk, validate=True))
        except Exception:  # noqa: BLE001 - a header line is not meant to decode
            continue
    return b"".join(out)


def test_default_strategy_ships_only_a_module_reference():
    """The bug, pinned: 85 bytes naming a module the partner does not have.

    Verified against a live partner endpoint -- the default strategy came back
    `ModuleNotFoundError: No module named 'appfl_bio_suite'`, raised in the unpickler
    before the probe body, whose imports are all local precisely to avoid this, was ever
    reached.
    """
    from globus_compute_sdk.serialize import DillCode

    payload = _decoded(DillCode().serialize(where_am_i))
    assert b"appfl_bio_suite" in payload
    assert b"def where_am_i" not in payload


def test_probe_is_shipped_as_source():
    """The fix: the function's own text travels with it, so no import is needed."""
    from globus_compute_sdk.serialize import CombinedCode

    assert b"def where_am_i" in _decoded(CombinedCode().serialize(where_am_i))


def test_smoke_test_installs_the_source_shipping_serializer(monkeypatch):
    """Guard the wiring, not just the strategy: the Executor must actually be told."""
    import globus_compute_sdk
    from globus_compute_sdk.serialize import CombinedCode

    executors = []

    class _FakeFuture:
        def result(self, timeout=None):
            return {
                "node": "example-node",
                "user": "example-user",
                "cwd": "/example",
                "python": "3.12",
            }

    class _FakeExecutor:
        def __init__(self, endpoint_id=None):
            self.endpoint_id = endpoint_id
            self.serializer = None
            executors.append(self)

        def submit(self, fn):
            self.submitted = fn
            return _FakeFuture()

        def shutdown(self, **kwargs):
            pass

    monkeypatch.setattr(globus_compute_sdk, "Executor", _FakeExecutor)

    result = smoke_test("bbbbbbbb-0000-0000-0000-000000000001", label="Site1")

    assert result.ok
    assert result.payload["user"] == "example-user"

    executor = executors[0]
    assert executor.submitted is where_am_i
    assert isinstance(executor.serializer.code_serializer, CombinedCode)


@pytest.mark.parametrize(
    ("message", "diagnosis"),
    [
        ("422 SEMANTICALLY_INVALID -- Identity failed to map to a local user name", "422 --"),
        ("('POST', 'https://compute/v3/x', None, 403, 'ENDPOINT_ACCESS_FORBIDDEN')", "403 --"),
        ("the service returned 404: no such endpoint", "404 --"),
    ],
)
def test_a_status_code_selects_its_diagnosis(message: str, diagnosis: str):
    assert explain_submit_failure(RuntimeError(message)).startswith(diagnosis)


@pytest.mark.parametrize("digits", ["403", "404", "422"])
def test_digits_inside_an_endpoint_uuid_are_not_a_status_code(digits: str):
    """These messages quote the endpoint UUID, and a UUID is hex.

    About one in fifty contains one of the three codes, and that used to pick the
    diagnosis: a scheduler timeout at such a site was explained as an identity-mapping
    failure, sending the partner to edit a file that was correct.
    """
    uuid = f"0a1b2c3d-{digits}0-4e5f-8a9b-0c1d2e3f{digits}a"
    timed_out = TimeoutError(f"no result from endpoint {uuid} within 300 seconds")
    assert explain_submit_failure(timed_out).startswith("Timed out waiting")
    unknown = RuntimeError(f"endpoint {uuid} rejected the task")
    assert explain_submit_failure(unknown).startswith("Unrecognized failure")


def test_the_status_attribute_is_used_when_the_message_does_not_carry_it():
    class ApiError(Exception):
        http_status = 403

    assert explain_submit_failure(ApiError("forbidden")).startswith("403 --")
