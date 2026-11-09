"""Masked wire sums preserve fine-mapping moments and reject partial/mixed rounds."""

import copy
import json

import numpy as np
import pytest
import torch

from appfl_bio_suite.experiments.fine_mapping.aggregator import FineMappingAggregator
from appfl_bio_suite.experiments.fine_mapping.secure_aggregation import pool_secure_payloads
from appfl_bio_suite.experiments.fine_mapping.trainer import (
    KEY_MANIFEST,
    SECURE_SCALE,
    mask_moments,
    mask_site_payload,
)

CONFIG = {"run_id": "test-run-1", "participants": ["a", "b", "c"]}
KEYS = {
    "a": {"b": "ab" * 32, "c": "ac" * 32},
    "b": {"a": "ab" * 32, "c": "bc" * 32},
    "c": {"a": "ac" * 32, "b": "bc" * 32},
}


def pack(manifest):
    return torch.frombuffer(bytearray(json.dumps(manifest).encode()), dtype=torch.uint8)


def inputs(cid, index):
    rng = np.random.default_rng(index)
    x = rng.integers(0, 3, (8 + index, 3)).astype(float)
    y = rng.normal(size=len(x))
    manifest = {
        "client_id": cid,
        "loci": [{"locus_id": "L0", "chrom": "1", "start_bp": 1, "end_bp": 9}],
        "geno_blocks": [
            {
                "locus_id": "L0",
                "pop": "EUR",
                "n": len(x),
                "n_incomplete_variants": 0,
                "variants": {
                    "snp_id": ["r1", "r2", "r3"],
                    "bp": [1, 2, 3],
                    "chrom": ["1"] * 3,
                    "a1": ["A"] * 3,
                    "a2": ["G"] * 3,
                },
            }
        ],
        "pheno_blocks": [{"instance": "L0_arch_rep0", "locus_id": "L0", "pop": "EUR", "n": len(x)}],
    }
    payload = {
        "G::L0::EUR": torch.from_numpy(x.T @ x),
        "u::L0::EUR": torch.from_numpy(x.sum(axis=0)),
        "c::L0_arch_rep0::EUR": torch.from_numpy(x.T @ y),
        "s::L0_arch_rep0::EUR": torch.tensor([y.sum(), y @ y, len(x)], dtype=torch.float64),
    }
    payload[KEY_MANIFEST] = pack(manifest)
    return payload, manifest


@pytest.fixture
def submissions():
    plain, masked = {}, {}
    for index, cid in enumerate(CONFIG["participants"]):
        payload, manifest = inputs(cid, index)
        plain[cid] = payload
        masked[cid] = mask_moments(payload, manifest, cid, CONFIG, KEYS[cid])
    return plain, masked


def test_round_trip_and_decode(submissions):
    plain, masked = submissions
    pooled, manifests = pool_secure_payloads(masked, CONFIG)
    assert set(manifests) == set(CONFIG["participants"])
    for name in pooled:
        if name == KEY_MANIFEST:
            continue
        expected = sum(p[name].numpy() for p in plain.values())
        np.testing.assert_allclose(pooled[name], expected, atol=3 / (2 * SECURE_SCALE), rtol=0)
        assert masked["a"][name].dtype == torch.int64
        assert not torch.equal(masked["a"][name], plain["a"][name])
    aggregator = object.__new__(FineMappingAggregator)
    aggregator.aggregator_configs = {"secure_aggregation": CONFIG}
    geno, pheno, _, originals = aggregator._decode(list(masked), masked)
    assert len(geno["L0"]) == 1
    assert geno["L0"][0].n == 27
    assert pheno["L0_arch_rep0"]["EUR"][0].n == 27
    assert originals == manifests
    np.testing.assert_array_equal(geno["L0"][0].G, sum(p["G::L0::EUR"] for p in plain.values()))


@pytest.mark.parametrize("failure", ["dropout", "plaintext", "run", "layout", "dtype", "key"])
def test_fail_closed(submissions, failure):
    plain, masked = submissions
    if failure == "dropout":
        del masked["c"]
    elif failure == "plaintext":
        masked["a"] = plain["a"]
    elif failure == "dtype":
        masked["a"]["G::L0::EUR"] = masked["a"]["G::L0::EUR"].double()
    else:
        payload, manifest = inputs("a", 0)
        config = copy.deepcopy(CONFIG)
        keys = dict(KEYS["a"])
        if failure == "run":
            config["run_id"] = "another-run"
        elif failure == "layout":
            manifest["geno_blocks"][0]["variants"]["a1"][0] = "T"
        else:
            keys["b"] = "ff" * 32
        masked["a"] = mask_moments(payload, manifest, "a", config, keys)
    with pytest.raises(ValueError, match="secure aggregation"):
        pool_secure_payloads(masked, CONFIG)


def test_disabled_mode_rejects_masks(submissions):
    plain, masked = submissions
    assert pool_secure_payloads(plain, None)[0] is None
    with pytest.raises(ValueError, match="disabled"):
        pool_secure_payloads(masked, None)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), 2**60, -(2**60)])
def test_invalid_moments(value):
    payload, manifest = inputs("a", 0)
    payload["G::L0::EUR"][0, 0] = value
    with pytest.raises(ValueError, match="range"):
        mask_moments(payload, manifest, "a", CONFIG, KEYS["a"])


def test_site_authorization_and_reuse_guard(tmp_path, monkeypatch):
    path = tmp_path / "keys.json"
    path.write_text(json.dumps({**CONFIG, "client_id": "a", "keys": KEYS["a"]}))
    monkeypatch.setenv("APPFL_FINEMAP_SECAGG_KEYS", str(path))
    payload, _ = inputs("a", 0)
    with pytest.raises(ValueError, match="authorization"):
        mask_site_payload(payload, "a", {**CONFIG, "run_id": "unapproved"})
    masked = mask_site_payload(payload, "a", CONFIG)
    assert masked["G::L0::EUR"].dtype == torch.int64
    with pytest.raises(FileExistsError):
        mask_site_payload(payload, "a", CONFIG)


def test_run_id_changes_masks():
    payload, manifest = inputs("a", 0)
    one = mask_moments(payload, copy.deepcopy(manifest), "a", CONFIG, KEYS["a"])
    two = mask_moments(payload, manifest, "a", {**CONFIG, "run_id": "new"}, KEYS["a"])
    assert not torch.equal(one["G::L0::EUR"], two["G::L0::EUR"])


@pytest.mark.parametrize(
    "update",
    [
        {"participants": ["a"]},
        {"participants": ["a", "a"]},
        {"run_id": ""},
        {"keys": {"a": "must-never-ship"}},
        {"scale": 1},
        {"protocol": "unknown"},
    ],
)
def test_invalid_public_configuration(update):
    from appfl_bio_suite.experiments.fine_mapping.trainer import secure_config

    with pytest.raises(ValueError):
        secure_config({**CONFIG, **update})


def test_native_appfl_cannot_preserve_moment_precision():
    """Compatibility probe: a legitimate large count is rounded before masking."""
    from appfl.privacy.secure_aggregator import SecureAggregator

    from appfl_bio_suite.experiments.fine_mapping.trainer import secure_config

    original = torch.tensor([2**40 + 1], dtype=torch.float64)
    masked, _ = SecureAggregator("a", ["a", "b"], b"probe-key", torch.device("cpu")).mask_update(
        {"moment": original}, round_id=0
    )
    assert masked.dtype == torch.float32
    assert masked.double().item() == 2**40  # unit-scale mask disappears at this magnitude
    assert masked.double().item() != original.item()
    with pytest.raises(ValueError, match="float32"):
        secure_config({**CONFIG, "protocol": "appfl"})
