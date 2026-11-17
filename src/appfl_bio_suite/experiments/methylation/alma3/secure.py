"""Adapt the suite's fixed-roster fine-mapping protocol to arbitrary named tensors.

This serial simulation provisions secrets locally. The coordinator sum helper
receives masked payloads only. No dropout recovery or host isolation is claimed.
"""

import hashlib
import hmac
import json
import secrets
import uuid

import numpy as np
import torch

KEY_MANIFEST = "manifest"
SCALE = 2**20


def provision(roster, context):
    config = {"run_id": f"alma3:{context}:{uuid.uuid4()}", "participants": list(roster)}
    keys = {cid: {} for cid in roster}
    for i, cid in enumerate(roster):
        for peer in roster[i + 1 :]:
            key = secrets.token_hex(32)
            keys[cid][peer] = key
            keys[peer][cid] = key
    return config, keys


def mask(values, cid, config, keys):
    # Wire-compatible specialization of fine_mapping.trainer.mask_moments.
    # Kept dependency-light: the upstream runtime must not import APPFL.
    roster = sorted(config["participants"])
    if len(roster) < 2 or len(set(roster)) != len(roster) or cid not in roster:
        raise ValueError("invalid secure roster")
    if set(keys) != set(roster) - {cid}:
        raise ValueError("missing peer secrets")
    public = {
        "protocol": "pairwise-u64-v1",
        "run_id": config["run_id"],
        "participants": roster,
        "scale": SCALE,
    }
    if not public["run_id"]:
        raise ValueError("empty secure context")
    blocks = {
        "geno_blocks": [],
        "pheno_blocks": [],
        "loci": [],
        "tensors": {k: list(v.shape) for k, v in sorted(values.items())},
    }
    layout = hashlib.sha256(json.dumps(blocks, sort_keys=True).encode()).hexdigest()
    context = json.dumps([public, layout], sort_keys=True).encode()
    derived, confirmations = {}, {}
    for peer, value in keys.items():
        key = bytes.fromhex(value)
        if len(key) != 32:
            raise ValueError("pairwise keys must have 32 bytes")
        pair = json.dumps(sorted([cid, peer])).encode()
        derived[peer] = hmac.digest(key, context + pair, "sha256")
        confirmations[peer] = hmac.digest(derived[peer], b"confirmation", "sha256").hex()
    result = {}
    for name, tensor in values.items():
        scaled = tensor.detach().cpu().numpy().astype(np.float64)
        scaled *= SCALE
        np.rint(scaled, out=scaled)
        if not np.isfinite(scaled).all() or (np.abs(scaled) >= 2**62 // len(roster)).any():
            raise ValueError("secure tensor outside fixed-point range")
        encoded = np.asarray(scaled.astype(np.int64)).view(np.uint64)
        del scaled
        for peer, key in derived.items():
            stream = hmac.digest(key, name.encode(), "sha256")
            pad = np.frombuffer(
                hashlib.shake_256(stream).digest(encoded.size * 8), dtype="<u8"
            ).reshape(encoded.shape)
            if cid < peer:
                np.add(encoded, pad, out=encoded)
            else:
                np.subtract(encoded, pad, out=encoded)
            del pad
        result[name] = torch.from_numpy(encoded.view(np.int64))
    manifest = {
        "client_id": cid,
        "geno_blocks": [],
        "pheno_blocks": [],
        "loci": [],
        "secure_aggregation": {**public, "layout": layout, "confirmations": confirmations},
    }
    result[KEY_MANIFEST] = torch.frombuffer(
        bytearray(json.dumps(manifest).encode()), dtype=torch.uint8
    )
    return result


def pool(payloads, config):
    from appfl_bio_suite.experiments.fine_mapping.secure_aggregation import pool_secure_payloads

    values, _ = pool_secure_payloads(payloads, config)
    return {k: v for k, v in values.items() if k != KEY_MANIFEST}


def wire_bytes(values):
    return sum(v.numel() * v.element_size() for v in values.values())


def sum_values(values, context, secure=True):
    """Convenience for in-process tests; production RPC masks on the worker."""
    if not secure:
        return {k: sum(v[k].double() for v in values.values()) for k in next(iter(values.values()))}
    config, keys = provision(list(values), context)
    return pool({cid: mask(v, cid, config, keys[cid]) for cid, v in values.items()}, config)


def tensor_state(model):
    return {k: p.detach().cpu().clone() for k, p in model.named_parameters() if p.requires_grad}


def load_state(model, state):
    parameters = dict(model.named_parameters())
    expected = {k for k, p in parameters.items() if p.requires_grad}
    if set(state) != expected:
        raise ValueError("checkpoint does not match selected trainable parameters")
    with torch.no_grad():
        for k, value in state.items():
            if value.shape != parameters[k].shape or not torch.isfinite(value).all():
                raise ValueError(f"invalid parameter {k}")
            parameters[k].copy_(value)
