"""The result-spool handle is written by a site. The coordinator must not trust its shape.

The round trip itself is covered where the trainer is (test_fine_mapping_shipped_parity);
this is the other half: what the coordinator does with a handle that is not what the
trainer would have written.
"""

from __future__ import annotations

import hashlib
import json

import pytest
import torch

from appfl_bio_suite.core.result_spool import SPOOL_KEY, fetch_spooled, spool_info


def handle(**overrides):
    info = {
        "format": "torch-save-v1",
        "path": "/work/spool/site-abc.pt",
        "size": 10,
        "sha256": "0" * 64,
        "client_id": "site",
        **overrides,
    }
    blob = json.dumps(info).encode()
    return {SPOOL_KEY: torch.frombuffer(bytearray(blob), dtype=torch.uint8)}


def test_a_well_formed_handle_is_accepted():
    assert spool_info(handle())["client_id"] == "site"
    assert spool_info({"manifest": torch.zeros(1)}) is None


@pytest.mark.parametrize(
    "overrides",
    [
        {"sha256": "/../../" + "a" * 57},  # 64 characters, and a path
        {"sha256": "Z" * 64},
        {"sha256": 0},
        {"size": True},
        {"size": 0},
        {"path": "relative/spool/x.pt"},
        {"path": ["/work/spool/x.pt"]},
        {"client_id": ["site"]},
        {"format": "pickle"},
    ],
)
def test_a_malformed_handle_is_refused(overrides):
    with pytest.raises(ValueError, match="malformed result-spool handle"):
        spool_info(handle(**overrides))


def test_a_handle_that_is_not_json_is_refused_as_a_handle():
    blob = torch.frombuffer(bytearray(b"not json"), dtype=torch.uint8)
    with pytest.raises(ValueError, match="malformed result-spool handle"):
        spool_info({SPOOL_KEY: blob})


class _Done:
    def __init__(self, value):
        self._value = value

    def result(self, timeout=None):
        return self._value


def test_the_fetched_copy_is_written_inside_the_spool_directory(tmp_path):
    """A site-supplied client id must not choose where the coordinator writes.

    The fetched file is created, filled with bytes the site sends, and then deleted. With
    the id used as written, "../../<name>" did all three to a path outside the directory.
    """
    victim = tmp_path / "victim-0000000000000000.pt"
    victim.write_bytes(b"the coordinator's own file")

    payload = tmp_path / "payload.pt"
    torch.save({"x": torch.arange(4)}, payload)
    data = payload.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    # Aim the traversal at `victim` exactly: same directory, same sha-derived suffix.
    hostile = victim.with_name(f"victim-{digest[:16]}.pt")
    hostile.write_bytes(b"the coordinator's own file")

    seen = []

    def submit(function, path, offset, length):
        seen.append((offset, length))
        return _Done(data[offset : offset + length])

    info = spool_info(handle(client_id="../victim", size=len(data), sha256=digest))
    spool_dir = tmp_path / "spool"
    fetched = fetch_spooled(info, submit, spool_dir, log=lambda message: None)

    assert torch.equal(fetched["x"], torch.arange(4))
    assert hostile.read_bytes() == b"the coordinator's own file"
    assert list(spool_dir.iterdir()) == []  # its own copy is removed once loaded
    assert seen == [(0, len(data))]
