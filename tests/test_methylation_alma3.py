"""ALMA3 research protocol tests: no upstream package, container or weights."""

import json
from dataclasses import replace

import numpy as np
import pytest
import torch

from appfl_bio_suite.experiments.fine_mapping.trainer import mask_moments
from appfl_bio_suite.experiments.methylation.alma3.calibration import (
    accepted,
    histogram,
    select_threshold,
)
from appfl_bio_suite.experiments.methylation.alma3.cohort import cohort
from appfl_bio_suite.experiments.methylation.alma3.config import Config
from appfl_bio_suite.experiments.methylation.alma3.model import Tiny, select_depth, source
from appfl_bio_suite.experiments.methylation.alma3.privacy import accounting, epsilon, train
from appfl_bio_suite.experiments.methylation.alma3.secure import mask, pool, provision, tensor_state
from appfl_bio_suite.experiments.methylation.alma3.study import (
    aggregate,
    freeze,
    guard_output,
    run,
)
from appfl_bio_suite.experiments.methylation.alma3.study import (
    test_frozen as evaluate_frozen,
)
from appfl_bio_suite.experiments.methylation.alma3.worker import Worker


@pytest.fixture(autouse=True)
def threads():
    previous = torch.get_num_threads()
    torch.set_num_threads(1)
    yield
    torch.set_num_threads(previous)


def tiny_config():
    return Config(
        seeds=(4,),
        depths=("heads",),
        noise_multipliers=(0.0,),
        site_sizes=(5, 3),
        validation_size=4,
        test_size=5,
        rounds=1,
        local_steps=1,
        coverage=(0.75,),
        bins=5,
        min_calibration=2,
    )


def test_secure_protocol_parity_and_weighted_appfl():
    protocol, keys = provision(["a", "b"], "test")
    values = {
        "weights": torch.tensor([1.1, -2.5], dtype=torch.float64),
        "log_scale": torch.tensor(1.2, dtype=torch.float64),
    }
    manifest = {"client_id": "a", "geno_blocks": [], "pheno_blocks": [], "loci": []}
    old = mask_moments(values, manifest, "a", protocol, keys["a"])
    new = mask(values, "a", protocol, keys["a"])
    for name in old:
        torch.testing.assert_close(old[name], new[name])
    a = {"w": torch.tensor([1.0]), "scale": torch.tensor(2.0)}
    b = {"w": torch.tensor([3.0]), "scale": torch.tensor(0.5)}
    plain = aggregate({"a": a, "b": b}, {"a": 1, "b": 3}, secure=False, context=None)
    masked = {
        "a": mask({k: v * 0.25 for k, v in a.items()}, "a", protocol, keys["a"]),
        "b": mask({k: v * 0.75 for k, v in b.items()}, "b", protocol, keys["b"]),
    }
    secured = aggregate(masked, {"a": 1, "b": 3}, secure=True, context=protocol)
    torch.testing.assert_close(plain["w"], torch.tensor([2.5]))
    torch.testing.assert_close(secured, plain)
    torch.testing.assert_close(secured["scale"], torch.tensor(0.875))
    double = aggregate(
        masked,
        {"a": 1, "b": 3},
        secure=True,
        context=protocol,
        reference={k: v.double() for k, v in a.items()},
    )
    assert all(v.dtype == torch.float64 for v in double.values())
    with pytest.raises(ValueError, match="complete"):
        pool({"a": masked["a"]}, protocol)
    with pytest.raises(ValueError, match="range"):
        mask({"w": torch.tensor([float("inf")])}, "a", protocol, keys["a"])
    changed = {**protocol, "run_id": "different"}
    with pytest.raises(ValueError, match="configuration"):
        pool(masked, changed)


def test_gradients_depth_proximal_and_privacy():
    config = tiny_config()
    model, profiles, chrom, pos, _, _ = source(fixture=True)
    select_depth(model, "heads")
    initial = {k: v.clone() for k, v in model.state_dict().items()}
    inputs, y = cohort(profiles, chrom, pos, config, 4, 0, "train", 0.75)
    train(model, inputs, y, steps=2, lr=0.2, clip=1, sigma=0)
    assert not torch.equal(initial["heads.weight"], model.heads.weight)
    torch.testing.assert_close(initial["direct_cpg.weight"], model.direct_cpg.weight)
    other = Tiny()
    other.load_state_dict(initial)
    select_depth(other, "heads")
    train(other, inputs, y, steps=2, lr=0.2, clip=1, sigma=0, mu=2)
    assert not torch.equal(other.heads.weight, model.heads.weight)
    a = accounting(2, 8, 5, 1e-5)
    assert a["rho"] == 1
    assert a["epsilon"] == pytest.approx(epsilon(1, 1e-5))
    assert accounting(8, 8, 5, 1e-5)["epsilon"] < a["epsilon"]
    assert accounting(0, 8, 5, 1e-5)["epsilon"] is None


def test_calibration_aggregates_gates_and_abstains():
    chosen = np.array([[1, 2], [1, 2], [0, -1]])
    scores = np.array([[0.99, 0.95], [0.2, 0.95], [0.9, 0]])
    truth = np.array([[1, 2], [1, 1], [0, -1]])
    counts = histogram(chosen, scores, truth, level=1, prior=[0.8], bins=10, sigma=0)
    assert counts.sum() == 1
    assert counts[1, 9] == 1
    assert not accepted(chosen, scores, [0.8, 0.9])[1, 1]
    assert select_threshold(np.zeros((2, 10)), floor=0.9, minimum=20) == 1.01
    counts = np.zeros((2, 10))
    counts[1, 9] = 100
    assert select_threshold(counts, floor=0.9, minimum=20) <= 0.9
    assert select_threshold(counts, floor=0.9, minimum=20, noise_std=100) == 1.01


def test_validation_stream_independent_of_test_and_test_gate(tmp_path):
    config = tiny_config()
    _, profiles, chrom, pos, _, _ = source(fixture=True)
    inputs, labels = cohort(profiles, chrom, pos, config, 4, 0, "validation", 0.75)
    other, other_labels = cohort(
        profiles, chrom, pos, replace(config, test_size=999), 4, 0, "validation", 0.75
    )
    for a, b in zip(inputs, other, strict=True):
        torch.testing.assert_close(a, b)
    torch.testing.assert_close(labels, other_labels)
    worker = Worker({"config": config.json(), "fixture": True})
    with pytest.raises(ValueError, match="frozen"):
        worker.handle({"op": "evaluate", "seed": 4, "depth": "heads", "split": "test"})
    with pytest.raises(ValueError, match="immutable"):
        guard_output(tmp_path / "release" / "output", tmp_path / "release")


def test_end_to_end_freeze_and_one_test(tmp_path, monkeypatch):
    from appfl_bio_suite.experiments.methylation.alma3 import worker as worker_module

    original = worker_module.cohort
    splits = []

    def checked(*args, **kwargs):
        splits.append(args[6])
        return original(*args, **kwargs)

    monkeypatch.setattr(worker_module, "cohort", checked)
    config = tiny_config()
    out = tmp_path / "validation"
    run(out, config, fixture=True)
    assert "test" not in splits
    records = json.loads((out / "results.json").read_text())
    assert {r["arm"] for r in records} == {
        "zero-shot",
        "pooled",
        "fedavg",
        "fedprox",
        "scratch-fedavg",
        "local-S0",
        "local-S1",
    }
    assert next(r for r in records if r["arm"] == "fedavg")["training_wire_bytes"] > 0
    locked = freeze(out)
    evaluate_frozen(locked, tmp_path / "test")
    assert "test" in splits
    with pytest.raises(FileExistsError):
        evaluate_frozen(locked, tmp_path / "test-again")
    assert (out / "comparison.png").is_file()
    # Tampering is rejected before any test access, even when the test marker exists.
    (out / records[0]["calibration"]).write_text("{}")
    with pytest.raises(ValueError, match="changed"):
        evaluate_frozen(locked, tmp_path / "tampered")


def test_private_noise_not_replayed_by_public_seed():
    a, b = Tiny(), Tiny()
    select_depth(a, "heads")
    select_depth(b, "heads")
    _, profiles, chrom, pos, _, _ = source(fixture=True)
    inputs, labels = cohort(profiles, chrom, pos, tiny_config(), 4, 0, "train", 0.75)
    torch.manual_seed(4)
    train(a, inputs, labels, steps=1, lr=0.1, clip=1, sigma=2)
    torch.manual_seed(4)
    train(b, inputs, labels, steps=1, lr=0.1, clip=1, sigma=2)
    assert not torch.equal(tensor_state(a)["heads.weight"], tensor_state(b)["heads.weight"])


def test_calibration_reuses_site_predictions_without_reusing_noise(tmp_path, monkeypatch):
    worker = Worker({"config": tiny_config().json(), "fixture": True})
    calls = []
    original = worker.predict

    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(worker, "predict", counted)
    for index, site in enumerate((0, 1, 0)):
        worker.handle(
            {
                "op": "calibrate",
                "seed": 4,
                "depth": "heads",
                "split": "calibration",
                "site": site,
                "coverage": 0.75,
                "level": 0,
                "prior": [],
                "sigma": 2,
                "output": str(tmp_path / f"{index}.npz"),
            }
        )
    assert len(calls) == 2
    with np.load(tmp_path / "0.npz") as first, np.load(tmp_path / "2.npz") as second:
        assert not np.array_equal(first["histogram"], second["histogram"])


def test_invalid_prediction_cannot_be_reported_as_accuracy(monkeypatch):
    worker = Worker({"config": tiny_config().json(), "fixture": True})
    inputs, _ = worker.data(4, 0, "validation", 0.75)
    monkeypatch.setattr(
        worker.model, "forward", lambda *batch: torch.full((len(batch[0]), 3), float("nan"))
    )
    with pytest.raises(ValueError, match="terminal probabilities"):
        worker.predict(inputs)


def test_shared_class_support_and_selected_arms(tmp_path):
    config = replace(
        tiny_config(),
        class_indices=(0, 2),
        minimum_per_class=1,
        arms=("zero-shot", "fedavg"),
        site_class_proportions=((0.8, 0.2), (0.2, 0.8)),
    )
    config.validate()
    _, profiles, chrom, pos, _, _ = source(fixture=True)
    for site in range(len(config.site_sizes)):
        for split in ("train", "calibration", "validation", "test"):
            _, labels = cohort(profiles, chrom, pos, config, 4, site, split, 0.75)
            assert set(labels.tolist()) == {0, 2}
    with pytest.raises(ValueError, match="accommodate"):
        replace(config, minimum_per_class=2).validate()
    with pytest.raises(ValueError, match="taxonomy"):
        cohort(profiles, chrom, pos, replace(config, class_indices=(0, 3)), 4, 0, "train", 0.75)
    run(tmp_path / "selected", config, fixture=True)
    records = json.loads((tmp_path / "selected" / "results.json").read_text())
    assert [r["arm"] for r in records] == ["zero-shot", "fedavg"]
    for record in records:
        for row in record["metrics"]:
            assert row["support_restricted_balanced_accuracy"] >= row["balanced_accuracy"] - 1e-12
    assert not list((tmp_path / "selected").glob("initial-*-True.npz"))
    summary = json.loads((tmp_path / "selected" / "summary.json").read_text())
    assert "zero-shot-known-support" in {row["arm"] for row in summary}
    assert (tmp_path / "selected" / "sites-heads-noise-0.0.png").is_file()


def test_cached_heads_preserve_gradients_and_probability_transform():
    from appfl_bio_suite.experiments.methylation.alma3.model import CachedHeads

    class Neural(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = torch.nn.Linear(5, 8)
            self.heads = torch.nn.Linear(8, 3)

        def logits_from_embedding(self, embedding):
            return self.heads(embedding)

    class MemberLike(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.neural = Neural()

        def probabilities(self, logits, embedding):
            return (logits / 1.7).log_softmax(-1)

        def forward(self, inputs):
            embedding = self.neural.encoder(inputs)
            return self.probabilities(self.neural.logits_from_embedding(embedding), embedding)

    model = MemberLike()
    select_depth(model, "heads")
    inputs, labels = torch.randn(4, 5), torch.tensor([0, 1, 1, 2])
    initial = tensor_state(model)
    features = model.neural.encoder(inputs).detach()
    cached = CachedHeads(model)
    torch.testing.assert_close(model(inputs), cached(features))
    train(model, (inputs,), labels, steps=3, lr=0.3, clip=1, sigma=0, mu=0.1)
    expected = tensor_state(model)
    from appfl_bio_suite.experiments.methylation.alma3.secure import load_state

    load_state(model, initial)
    train(cached, (features,), labels, steps=3, lr=0.3, clip=1, sigma=0, mu=0.1)
    for name, value in tensor_state(model).items():
        torch.testing.assert_close(value, expected[name])


def test_signature_rotation_preserves_information_and_is_label_blind():
    from appfl_bio_suite.experiments.methylation.alma3.cohort import rotate_signature

    _, profiles, chrom, pos, _, _ = source(fixture=True)
    support = np.arange(3)
    logits = np.log(profiles / (1 - profiles))
    rotated = rotate_signature(logits, profiles, support, 75)
    recovered = rotate_signature(rotated, profiles, support, -75)
    np.testing.assert_allclose(recovered, logits, atol=3e-6)
    before = np.linalg.norm(logits[:, None] - logits[None, :], axis=2)
    after = np.linalg.norm(rotated[:, None] - rotated[None, :], axis=2)
    np.testing.assert_allclose(after, before, atol=3e-6)
    # A patient's transform cannot depend on which other patients share its batch.
    np.testing.assert_allclose(
        rotate_signature(logits[:1], profiles, support, 75), rotated[:1], atol=3e-6
    )
    config = replace(
        tiny_config(), class_indices=(0, 1, 2), minimum_per_class=1, signature_rotation_degrees=75
    )
    config.validate()
    original = replace(config, signature_rotation_degrees=0)
    shifted, labels = cohort(profiles, chrom, pos, config, 4, 0, "validation", 0.75)
    unshifted, unshifted_labels = cohort(
        profiles, chrom, pos, config, 4, 0, "validation", 0.75, unshifted=True
    )
    baseline, baseline_labels = cohort(
        profiles, chrom, pos, original, 4, 0, "validation", 0.75, unshifted=True
    )
    assert torch.equal(labels, unshifted_labels)
    assert torch.equal(unshifted_labels, baseline_labels)
    for a, b in zip(unshifted, baseline, strict=True):
        assert torch.equal(a, b)
    assert not torch.equal(shifted[0], baseline[0])
    with pytest.raises(ValueError, match="three explicit classes"):
        replace(config, class_indices=(), minimum_per_class=0).validate()


def test_five_sites_exact_skew_preserves_all_classes():
    from appfl_bio_suite.experiments.methylation.alma3.cohort import class_counts

    ratios = (
        (0.90, 0.06, 0.03, 0.01),
        (0.01, 0.90, 0.06, 0.03),
        (0.03, 0.01, 0.90, 0.06),
        (0.06, 0.03, 0.01, 0.90),
        (0.47, 0.01, 0.01, 0.51),
    )
    config = replace(
        tiny_config(),
        site_sizes=(100,) * 5,
        validation_size=100,
        test_size=100,
        class_indices=(0, 1, 2, 3),
        minimum_per_class=1,
        site_class_proportions=ratios,
    )
    config.validate()
    profiles = np.full((4, 24), 0.5, dtype=np.float32)
    chrom, pos = torch.zeros(24, dtype=torch.long), torch.zeros(24)
    for seed in (42, 43):
        for site, proportions in enumerate(ratios):
            for split in ("train", "calibration", "validation", "test"):
                _, labels = cohort(profiles, chrom, pos, config, seed, site, split, 0.75)
                np.testing.assert_array_equal(
                    np.bincount(labels, minlength=4),
                    np.rint(np.array(proportions) * 100).astype(int),
                )
    counts = class_counts(7, ratios[0], 1)
    assert counts.sum() == 7 and counts.min() >= 1
    for change in (
        {"site_class_proportions": ratios[:2]},
        {"site_class_proportions": ((1.0, 0.0, 0.0, 0.0),) * 5},
        {"site_class_proportions": ((0.5, 0.5, 0.5, 0.5),) * 5},
        {"minimum_per_class": 0},
    ):
        with pytest.raises(ValueError):
            replace(config, **change).validate()


def test_evaluation_feature_cache_recomputes_mutable_heads(monkeypatch):
    from types import SimpleNamespace

    from appfl_bio_suite.experiments.methylation.alma3 import worker as worker_module

    class Neural(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = torch.nn.Linear(24, 8)
            self.heads = torch.nn.Linear(8, 3)
            self.config = SimpleNamespace(direct_logit_correction=None)

        def embed(self, beta, *unused):
            return self.encoder(beta)

        def logits_from_embedding(self, embedding):
            return self.heads(embedding)

    class MemberLike(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.neural = Neural()

        def probabilities(self, logits, embedding):
            return logits.log_softmax(-1)

        def forward(self, *inputs):
            embedding = self.neural.embed(*inputs)
            return self.probabilities(self.neural.logits_from_embedding(embedding), embedding)

    worker = Worker({"config": tiny_config().json(), "fixture": True})
    worker.model = MemberLike()
    select_depth(worker.model, "heads")
    monkeypatch.setattr(worker_module, "Member", MemberLike)
    request = dict(seed=4, site=0, split="validation", coverage=0.75, depth="heads")
    raw, labels = worker.data(4, 0, "validation", 0.75)
    features, cached_labels, embedded = worker.evaluation_data(request)
    assert embedded and torch.equal(labels, cached_labels)
    before = worker.predict(raw)
    torch.testing.assert_close(before, worker.predict(features, embedded=True))
    with torch.no_grad():
        worker.model.neural.heads.bias[0] += 1
    reused, _, _ = worker.evaluation_data(request)
    assert reused[0].data_ptr() == features[0].data_ptr()
    after = worker.predict(reused, embedded=True)
    torch.testing.assert_close(worker.predict(raw), after)
    assert not torch.equal(before, after)
    assert not worker.can_cache_features({**request, "scratch": True})
    assert not worker.can_cache_features({**request, "depth": "trunk"})
