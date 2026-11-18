"""Isolated ALMA3 runtime RPC worker (numpy/torch/sklearn; no APPFL dependency)."""

import argparse
import copy
import json
import sys
import time
import traceback
from dataclasses import replace
from pathlib import Path

import numpy as np
import torch

from ..evaluation import metrics
from .calibration import accepted, histogram, path_scores
from .cohort import cohort
from .config import Config
from .model import CachedHeads, Member, reset_neural, select_depth, source
from .privacy import train
from .secure import load_state, mask, tensor_state


def save_tensors(path, tensors):
    np.savez(path, **{k: v.detach().cpu().numpy() for k, v in tensors.items()})


def read_tensors(path):
    with np.load(path, allow_pickle=False) as archive:
        return {k: torch.from_numpy(archive[k].copy()) for k in archive.files}


class Worker:
    def __init__(self, settings):
        self.settings = settings
        self.config = Config(**settings["config"])
        self.config.validate()
        torch.set_num_threads(self.config.threads)
        self.device = settings.get("device", "cpu")
        self.model, self.profiles, self.chrom, self.positions, self.paths, self.identity = source(
            root=settings.get("release"),
            member=settings.get("member"),
            device=self.device,
            fixture=settings.get("fixture", False),
        )
        self.base = copy.deepcopy(self.model.state_dict())
        self.active = None
        self.cache = {}
        self.baselines = {}
        self.input_cache = {}
        self.input_cache_bytes = 0
        self.feature_cache = {}

    def prepare(self, seed, depth, scratch):
        key = seed, depth, scratch
        if self.active != key:
            self.model.load_state_dict(self.base)
            if scratch:
                reset_neural(self.model, seed)
            select_depth(self.model, depth)
            self.active = key
            self.model.zero_grad(set_to_none=True)
        self.model.eval()

    def data(self, seed, site, split, coverage, unshifted=False):
        key = seed, site, split, coverage, unshifted
        if key in self.input_cache:
            return self.input_cache[key][0]
        inputs, y = cohort(
            self.profiles,
            self.chrom,
            self.positions,
            self.config,
            seed,
            site,
            split,
            coverage,
            unshifted=unshifted,
        )
        value = tuple(x.to(self.device) for x in inputs), y.to(self.device)
        size = sum(v.numel() * v.element_size() for v in (*value[0], value[1]))
        # Bound even expanded reference tensors conservatively. No data is written
        # to coordinator storage. Large cohorts are regenerated, not retained.
        budget = 256 * 1024**2
        if size <= budget:
            while self.input_cache and self.input_cache_bytes + size > budget:
                _, removed = self.input_cache.pop(next(iter(self.input_cache)))
                self.input_cache_bytes -= removed
            self.input_cache[key] = value, size
            self.input_cache_bytes += size
        return value

    def can_cache_features(self, request):
        return (
            isinstance(self.model, Member)
            and request["depth"] == "heads"
            and not request.get("scratch", False)
            and not hasattr(self.model, "partition")
            and not self.model.neural.config.direct_logit_correction
        )

    def training_data(self, request):
        if "class_counts" in request:
            return self.sweep_training_data(request)
        cached = self.can_cache_features(request)
        all_data = []
        for site in request["sites"]:
            key = request["seed"], site, self.config.coverage[-1]
            if cached and key in self.feature_cache:
                value = self.feature_cache[key]
            else:
                inputs, labels = self.data(*key[:2], "train", key[2])
                if cached:
                    with torch.no_grad():
                        embedding = torch.cat(
                            [
                                self.model.neural.embed(*(v[i : i + 8] for v in inputs))
                                for i in range(0, len(labels), 8)
                            ]
                        )
                    value = (embedding,), labels
                    # At most 64 site cohorts; embeddings stay entirely in the worker.
                    if len(self.feature_cache) >= 64:
                        self.feature_cache.pop(next(iter(self.feature_cache)))
                    self.feature_cache[key] = value
                else:
                    value = inputs, labels
            all_data.append(value)
        inputs = tuple(torch.cat([d[0][j] for d in all_data]) for j in range(len(all_data[0][0])))
        labels = torch.cat([d[1] for d in all_data])
        return CachedHeads(self.model) if cached else self.model, inputs, labels

    def sweep_training_data(self, request):
        """Paired class-wise patient prefixes; only embeddings stay in worker memory."""
        if not self.can_cache_features(request) or request.get("sigma", 0) != 0:
            raise ValueError("sweep reservoirs require nonprivate pretrained heads")
        support = self.config.class_indices
        counts = np.asarray(request["class_counts"])
        if (
            counts.shape != (len(self.config.site_sizes), len(support))
            or counts.dtype.kind not in "iu"
            or (counts < 1).any()
            or not np.array_equal(counts.sum(1), self.config.site_sizes)
        ):
            raise ValueError("invalid sweep class counts")
        if any(n % len(support) for n in self.config.site_sizes):
            raise ValueError("reservoir cohort sizes must divide evenly by class count")
        pool_config = replace(
            self.config,
            site_class_proportions=tuple(
                (1 / len(support),) * len(support) for _ in self.config.site_sizes
            ),
        )
        all_features, all_labels = [], []
        for site in request["sites"]:
            key = "sweep-reservoir", request["seed"], site
            if key not in self.feature_cache:
                features, labels = [], []
                for replicate in range(1, len(support) + 1):
                    inputs, y = cohort(
                        self.profiles,
                        self.chrom,
                        self.positions,
                        pool_config,
                        request["seed"],
                        site,
                        "train",
                        self.config.coverage[-1],
                        replicate=replicate,
                    )
                    inputs = tuple(v.to(self.device) for v in inputs)
                    with torch.no_grad():
                        features.append(
                            torch.cat(
                                [
                                    self.model.neural.embed(*(v[i : i + 8] for v in inputs))
                                    for i in range(0, len(y), 8)
                                ]
                            )
                        )
                    labels.append(y.to(self.device))
                if len(self.feature_cache) >= 64:
                    self.feature_cache.pop(next(iter(self.feature_cache)))
                self.feature_cache[key] = torch.cat(features), torch.cat(labels)
            features, y = self.feature_cache[key]
            indices = torch.cat(
                [
                    torch.where(y == label)[0][: int(counts[site, j])]
                    for j, label in enumerate(support)
                ]
            )
            if len(indices) != self.config.site_sizes[site]:
                raise ValueError("reservoir cannot supply requested class counts")
            all_features.append(features[indices])
            all_labels.append(y[indices])
        return CachedHeads(self.model), (torch.cat(all_features),), torch.cat(all_labels)

    def evaluation_data(self, request):
        key = (
            request["seed"],
            request["site"],
            request["split"],
            request["coverage"],
            request.get("unshifted", False),
        )
        embedded = self.can_cache_features(request) and not request.get("baseline")
        if embedded and key in self.feature_cache:
            inputs, labels = self.feature_cache[key]
        else:
            inputs, labels = self.data(*key)
            if embedded:
                with torch.no_grad():
                    inputs = (
                        torch.cat(
                            [
                                self.model.neural.embed(*(v[i : i + 8] for v in inputs))
                                for i in range(0, len(labels), 8)
                            ]
                        ),
                    )
                if len(self.feature_cache) >= 64:
                    self.feature_cache.pop(next(iter(self.feature_cache)))
                self.feature_cache[key] = inputs, labels
        return inputs, labels, embedded

    def predict(self, inputs, baseline=None, *, embedded=False):
        if baseline and self.settings.get("fixture"):
            raise ValueError("release baselines require real ALMA3 weights")
        if baseline and baseline not in self.baselines:
            from alma3.final_ensemble import FinalEnsemble, SparseTemplates
            from alma3.release import validate_release

            root = Path(self.settings["release"])
            metadata = validate_release(root)
            if baseline == "template":
                row = next(
                    r for r in metadata["ensemble"]["members"] if r["kind"] == "sparse_templates"
                )
                state = torch.load(root / row["templates"], weights_only=True, map_location="cpu")
                self.baselines[baseline] = SparseTemplates(
                    state["means"], state["group_index"], row["sites"], row["temperature"]
                ).to(self.device)
            elif baseline == "ensemble":
                self.baselines[baseline] = FinalEnsemble(root, device=self.device)
            else:
                raise ValueError("invalid release baseline")
        with torch.no_grad():
            # Bounded forward memory. Deliberately research FP32, not release certification.
            outputs = []
            for i in range(0, len(inputs[0]), 8):
                batch = tuple(v[i : i + 8] for v in inputs)
                if baseline == "template":
                    result = self.baselines[baseline](*batch[:2])
                elif baseline == "ensemble":
                    result = self.baselines[baseline].predict(*batch)
                elif embedded:
                    result = CachedHeads(self.model)(*batch)
                else:
                    result = self.model(*batch)
                if not torch.isfinite(result).all() or not torch.allclose(
                    result.exp().sum(1),
                    torch.ones(len(result), device=result.device, dtype=result.dtype),
                    atol=3e-5,
                    rtol=0,
                ):
                    raise ValueError("nonfinite or unnormalized terminal probabilities")
                outputs.append(result.cpu())
            return torch.cat(outputs)

    def handle(self, request):
        started = time.monotonic()
        op = request["op"]
        if op == "identity":
            return {
                **self.identity,
                "levels": self.paths.shape[1],
                "classes": len(self.paths),
                "cpgs": self.profiles.shape[1],
                "torch": torch.__version__,
            }
        self.prepare(request["seed"], request["depth"], request.get("scratch", False))
        if request.get("state"):
            load_state(self.model, read_tensors(request["state"]))
        if op == "initial":
            # Reinitialize even when the preceding request used this same mode.
            self.active = None
            self.prepare(request["seed"], request["depth"], request.get("scratch", False))
            result = tensor_state(self.model)
        elif op == "train":
            training_model, inputs, labels = self.training_data(request)
            train(
                training_model,
                inputs,
                labels,
                steps=self.config.local_steps,
                lr=self.config.learning_rate,
                clip=self.config.clip_norm,
                sigma=request["sigma"],
                mu=request.get("mu", 0.0),
                microbatch_size=request.get("gradient_batch_size", 1),
            )
            result = tensor_state(self.model)
            result = {k: v * request.get("weight", 1.0) for k, v in result.items()}
        elif op in {"calibrate", "evaluate"}:
            split = request["split"]
            if op == "calibrate" and split != "calibration":
                raise ValueError("calibration must use the separate calibration fold")
            if op == "evaluate" and split not in {"validation", "test"}:
                raise ValueError("evaluation split must be validation or test")
            # The worker also enforces test authorization, independent of the driver.
            if split == "test" and not self.settings.get("allow_test", False):
                raise ValueError("test evaluation requires a frozen study")
            cache_key = (
                request.get("baseline"),
                request.get("state"),
                request["seed"],
                request["depth"],
                request.get("scratch", False),
                request["site"],
                split,
                request["coverage"],
                request.get("unshifted", False),
            )
            if cache_key not in self.cache:
                inputs, y, embedded = self.evaluation_data(request)
                logp = self.predict(inputs, request.get("baseline"), embedded=embedded)
                taxonomy = getattr(self.model, "taxonomy", None)
                chosen, scores = path_scores(logp, self.paths, taxonomy)
                # This cache holds predictions; inputs use a separate bounded cache.
                # Keep sequential hierarchy queries from repeating expensive forwards.
                if len(self.cache) >= 128:
                    self.cache.pop(next(iter(self.cache)))
                self.cache[cache_key] = (logp, y.cpu().numpy(), chosen, scores)
            logp, y, chosen, scores = self.cache[cache_key]
            if op == "calibrate":
                values = histogram(
                    chosen,
                    scores,
                    self.paths[y],
                    level=request["level"],
                    prior=request["prior"],
                    bins=self.config.bins,
                    sigma=request["sigma"],
                )
                result = {"histogram": torch.from_numpy(values)}
            else:
                report = metrics(logp.numpy(), y)
                if self.config.class_indices:
                    # Diagnostic only: known synthetic support, not a tuned model.
                    # Reveals gains achievable by suppressing unsupported classes.
                    support = np.asarray(self.config.class_indices)
                    restricted = support[logp[:, support].argmax(1).numpy()]
                    report["support_restricted_balanced_accuracy"] = float(
                        np.mean(
                            [(restricted[y == label] == label).mean() for label in np.unique(y)]
                        )
                    )
                decisions = accepted(chosen, scores, request["thresholds"])
                report["hierarchy"] = [
                    {
                        "call_rate": float(decisions[:, j].mean()),
                        "accuracy": float(
                            (chosen[decisions[:, j], j] == self.paths[y][decisions[:, j], j]).mean()
                        )
                        if decisions[:, j].any()
                        else None,
                    }
                    for j in range(self.paths.shape[1])
                ]
                report["seconds"] = time.monotonic() - started
                return report
        else:
            raise ValueError(f"unknown worker operation {op}")
        if request.get("secure"):
            result = mask(result, request["cid"], request["secure"], request["keys"])
        save_tensors(request["output"], result)
        return {
            "seconds": time.monotonic() - started,
            "tensor_bytes": sum(v.numel() * v.element_size() for v in result.values()),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", required=True, type=Path)
    args = parser.parse_args()
    worker = Worker(json.loads(args.settings.read_text()))
    for line in sys.stdin:
        try:
            request = json.loads(line)
            result = worker.handle(request)
            print(json.dumps({"ok": True, "result": result}, allow_nan=False), flush=True)
        except Exception as error:
            traceback.print_exc(file=sys.stderr)
            print(json.dumps({"ok": False, "error": str(error)}), flush=True)


if __name__ == "__main__":
    main()
