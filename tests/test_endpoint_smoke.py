"""The smoke probe must not require a partner to have installed this package.

An experiment's dataset and model are shipped to a worker as source text precisely so a
partner need not install `appfl_bio_suite` -- see experiments/*/dataset.py. A readiness
check that is stricter than the run it clears is a check that fails on correctly
configured endpoints, which is the one thing it must never do.
"""

from __future__ import annotations

import base64

from appfl_bio_suite.core.endpoint import smoke_test, where_am_i


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
            return {"node": "test-worker", "user": "flamby_svc", "cwd": "/", "python": "3.12.13"}

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

    result = smoke_test("00000000-0000-0000-0000-000000000001", label="Test site")

    assert result.ok
    assert result.payload["user"] == "flamby_svc"

    executor = executors[0]
    assert executor.submitted is where_am_i
    assert isinstance(executor.serializer.code_serializer, CombinedCode)
