"""Optional ALMA3 adapter and deterministic tiny contract stand-in."""

import importlib.metadata
import json
from pathlib import Path

import numpy as np
import torch
from torch import nn

from . import REVISION, VERSION


class Tiny(nn.Module):
    """Analytic prototype classifier. Never reported as ALMA3 pretrained weights."""

    def __init__(self):
        super().__init__()
        self.direct_cpg = nn.Linear(24, 12)
        self.trunk = nn.Linear(12, 12)
        self.heads = nn.Linear(12, 3)
        with torch.no_grad():
            self.direct_cpg.weight.zero_()
            self.direct_cpg.bias.zero_()
            self.direct_cpg.weight[:, :12].copy_(torch.eye(12))
            self.trunk.weight.copy_(torch.eye(12))
            self.trunk.bias.zero_()
            self.heads.weight.zero_()
            self.heads.bias.zero_()
            for j in range(3):
                self.heads.weight[j, j * 4 : (j + 1) * 4] = 3

    def forward(self, beta, observed, uncertainty, chrom, positions):
        x = torch.where(observed, 2 * beta - 1, 0)
        return self.heads(self.trunk(self.direct_cpg(x))).log_softmax(-1)


class Member(nn.Module):
    """One released neural member, including its released probability transforms.

    This is a single-member adaptation, explicitly distinguished from the frozen
    ten-member ensemble. No upstream files or inference wrappers are patched.
    """

    def __init__(self, root, metadata, key, device):
        super().__init__()
        from alma3.dx import load_dx
        from alma3.final_ensemble import DiagnosisResidual

        rows = metadata["ensemble"]["members"]
        self.row = next((r for r in rows if r["key"] == key), None)
        if self.row is None or self.row["kind"] == "sparse_templates":
            raise ValueError("select a released neural member key")
        self.taxonomy = metadata["taxonomy"]
        self.neural = load_dx(root / self.row["export"], device=device)
        if "residual" in self.row:
            state = torch.load(root / self.row["residual"], weights_only=True, map_location="cpu")[
                "state"
            ]
            self.residual = DiagnosisResidual(state).to(device)
        if "partition" in self.row:
            self.register_buffer(
                "partition", torch.from_numpy(np.load(root / self.row["partition"]))
            )

    def raw(self, inputs):
        logits, embedding = self.neural(*inputs, _return_embedding=True)
        return self.probabilities(logits, embedding)

    def probabilities(self, logits, embedding):
        from alma3.final_ensemble import mix_stages
        from alma3.terminal import terminal_log_probabilities, terminal_paths

        logits = {k: v.float() for k, v in logits.items()}
        native = terminal_log_probabilities(
            logits, self.taxonomy, output_mode=self.neural.config.output_mode
        )
        if hasattr(self, "residual"):
            return self.residual(embedding.double(), native.double())
        if "temperatures" in self.row:
            calibrated = terminal_log_probabilities(
                {k: v / self.row["temperatures"][j] for j, (k, v) in enumerate(logits.items())},
                self.taxonomy,
                output_mode="conditional_heads_v1",
            )
            if "native_roles" in self.row:
                return mix_stages(
                    native,
                    calibrated,
                    self.row["native_roles"],
                    torch.tensor(terminal_paths(self.taxonomy), device=native.device),
                )
            return calibrated
        return native

    def forward(self, *inputs):
        native = self.raw(inputs)
        if not hasattr(self, "partition"):
            return native
        beta, observed, uncertainty, chrom, positions = inputs
        masks = [observed & (self.partition != j)[None] for j in range(2)]
        eligible = torch.stack([m.sum(1) >= 5000 for m in masks], 1).all(1)
        views = [native] + [
            self.raw(
                (beta, torch.where(eligible[:, None], m, observed), uncertainty, chrom, positions)
            )
            for m in masks
        ]
        average = torch.logsumexp(torch.stack(views), 0) - np.log(3)
        return torch.where(eligible[:, None], average, native)


class CachedHeads(nn.Module):
    """Worker-local frozen embeddings; identical trainable parameters and transforms."""

    def __init__(self, member):
        super().__init__()
        self.member = member

    def forward(self, embedding):
        logits = self.member.neural.logits_from_embedding(embedding)
        return self.member.probabilities(logits, embedding)


def select_depth(model, depth):
    if depth not in {"heads", "trunk", "full"}:
        raise ValueError("unknown fine-tuning depth")
    for name, parameter in model.named_parameters():
        parts = name.split(".")
        parameter.requires_grad_(
            depth == "full"
            or "heads" in parts
            or "residual" in parts
            or (depth == "trunk" and "trunk" in parts)
        )
    if not any(p.requires_grad for p in model.parameters()):
        raise ValueError("depth selected no trainable parameters")


def reset_neural(model, seed):
    # Reset every learnable parameter, including custom attention/token parameters.
    # Only taxonomy and genomic-reference buffers survive the random-start arm.
    torch.manual_seed(seed)
    with torch.no_grad():
        for module in model.modules():
            if hasattr(module, "reset_parameters"):
                module.reset_parameters()
            elif isinstance(module, nn.MultiheadAttention):
                module._reset_parameters()
            else:
                for parameter in module.parameters(recurse=False):
                    if parameter.ndim >= 2:
                        nn.init.xavier_uniform_(parameter)
                    else:
                        parameter.zero_()
        for name, p in model.named_parameters():
            if name.endswith("log_scale"):
                p.fill_(np.log(20.0))
        if hasattr(model, "residual"):
            model.residual.mean.zero_()
            model.residual.scale.fill_(1)


def source(root=None, member=None, device="cpu", fixture=False):
    if fixture:
        profiles = np.full((3, 24), 0.15, dtype=np.float32)
        for j in range(3):
            profiles[j, j * 4 : (j + 1) * 4] = 0.85
        return (
            Tiny().to(device),
            profiles,
            torch.zeros(24, dtype=torch.long),
            torch.linspace(0, 1, 24),
            np.array([[0], [1], [2]]),
            {"kind": "tiny-stand-in", "pretrained": False, "version": "fixture-v1"},
        )
    if importlib.metadata.version("alma3") != VERSION:
        raise ValueError(f"ALMA3 runtime must be exactly {VERSION}")
    if torch.__version__.split("+")[0] != "2.14.1":
        raise ValueError("ALMA3 must run in its isolated Torch 2.14.1 runtime")
    from alma3.release import validate_release
    from alma3.terminal import terminal_paths

    root = Path(root).resolve()
    metadata = validate_release(root, verify_hashes=True)
    row = next(r for r in metadata["ensemble"]["members"] if r["kind"] == "sparse_templates")
    profiles = torch.load(root / row["templates"], weights_only=True, map_location="cpu")["means"]
    if profiles.shape != (109, 387849) or not torch.isfinite(profiles).all():
        raise ValueError("invalid released class profiles")
    return (
        Member(root, metadata, member, device).to(device),
        profiles.numpy(),
        metadata["cpg"].chr_id,
        metadata["cpg"].pos,
        np.asarray(terminal_paths(metadata["taxonomy"])),
        {
            "kind": "alma3-neural-member",
            "pretrained": True,
            "version": VERSION,
            "source_revision": REVISION,
            "member": member,
            "release_sha256": metadata["manifest_sha256"],
            "release_metadata": metadata["release"],
            "precision": "fp32-research",
            "license": (root / "LICENSE").read_text(),
            "member_config": json.loads(
                (
                    root
                    / next(
                        r["export"] for r in metadata["ensemble"]["members"] if r["key"] == member
                    )
                    / "config.json"
                ).read_text()
            ),
        },
    )
