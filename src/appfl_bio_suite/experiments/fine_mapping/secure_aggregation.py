"""Coordinator-side fixed-roster masked sum; never receives pairwise secrets."""

import copy
import json

import numpy as np
import torch

from appfl_bio_suite.experiments.fine_mapping.trainer import (
    KEY_MANIFEST,
    SECURE_SCALE,
    secure_config,
    secure_layout,
)


def pool_secure_payloads(local_models, config):
    """Validate all submissions before summing; return pooled payload and manifests."""
    for cid, payload in local_models.items():
        if KEY_MANIFEST not in payload:
            raise ValueError(f"site {cid} returned no '{KEY_MANIFEST}' block")
    manifests = {
        cid: json.loads(payload[KEY_MANIFEST].detach().cpu().numpy().tobytes())
        for cid, payload in local_models.items()
    }
    if config is None:
        if any("secure_aggregation" in m for m in manifests.values()):
            raise ValueError("masked payload received with secure aggregation disabled")
        return None, manifests
    public = secure_config(config)
    roster = public["participants"]
    if set(local_models) != set(roster):
        raise ValueError("secure aggregation requires the complete configured roster")
    first = local_models[roster[0]]
    layout = secure_layout(manifests[roster[0]], first)
    for cid in roster:
        payload, manifest = local_models[cid], manifests[cid]
        reported = manifest.get("secure_aggregation", {})
        if manifest.get("client_id") != cid:
            raise ValueError("secure aggregation client identity mismatch")
        if any(reported.get(k) != v for k, v in public.items()):
            raise ValueError("secure aggregation configuration mismatch or plaintext payload")
        if reported.get("layout") != layout or secure_layout(manifest, payload) != layout:
            raise ValueError("secure aggregation requires identical ordered block/variant layouts")
        confirmations = reported.get("confirmations", {})
        if set(confirmations) != set(roster) - {cid}:
            raise ValueError("secure aggregation missing pairwise confirmations")
        for peer, confirmation in confirmations.items():
            other = manifests[peer].get("secure_aggregation", {}).get("confirmations", {})
            if not isinstance(confirmation, str) or other.get(cid) != confirmation:
                raise ValueError("secure aggregation pairwise key/context mismatch")
        for name, tensor in payload.items():
            if name != KEY_MANIFEST and tensor.dtype != torch.int64:
                raise ValueError("secure aggregation requires int64 wire tensors")
    pooled = {}
    for name in first:
        if name == KEY_MANIFEST:
            continue
        total = np.zeros(tuple(first[name].shape), dtype=np.uint64)
        for cid in roster:
            total += local_models[cid][name].detach().cpu().numpy().view(np.uint64)
        pooled[name] = torch.from_numpy(total.view(np.int64).astype(np.float64) / SECURE_SCALE)
    manifest = copy.deepcopy(manifests[roster[0]])
    del manifest["secure_aggregation"]
    manifest["client_id"] = "secure-pooled"
    for kind in ("geno_blocks", "pheno_blocks"):
        for i, block in enumerate(manifest[kind]):
            for field in ("n", "n_incomplete_variants"):
                if field in block:
                    block[field] = sum(m[kind][i][field] for m in manifests.values())
    pooled[KEY_MANIFEST] = torch.frombuffer(
        bytearray(json.dumps(manifest).encode()), dtype=torch.uint8
    )
    return pooled, manifests
