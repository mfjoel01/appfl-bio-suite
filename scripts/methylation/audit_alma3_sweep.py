#!/usr/bin/env python3
"""Audit completed sweep records, balanced evaluation and saved model hashes."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path


def audit(out):
    out = Path(out).resolve()
    spec = json.loads((out / "specification.json").read_text())
    records = json.loads((out / "results.json").read_text())
    assert (out / "COMPLETE.json").exists(), "sweep is incomplete"
    sites, classes = len(spec["extreme_counts"]), len(spec["class_indices"])
    arms = ["zero-shot", "pooled", "fedavg"] + [f"local-S{i}" for i in range(sites)]
    expected = set(itertools.product(spec["seeds"], spec["levels"], arms))
    actual = {(r["seed"], r["convergence"], r["arm"]) for r in records}
    assert actual == expected and len(records) == len(expected), "missing or duplicate arms"
    zero, state_hashes = {}, {}
    evaluated = 0
    for r in records:
        counts, t = r["class_counts"], r["convergence"]
        assert all(sum(row) == spec["patients_per_site"] and min(row) >= 1 for row in counts)
        assert [sum(row[j] for row in counts) for j in range(classes)] == [
            sites * v for v in spec["target_counts"]
        ]
        for i, row in enumerate(counts):
            for j, n in enumerate(row):
                desired = (1 - t) * spec["extreme_counts"][i][j] + t * spec["target_counts"][j]
                assert abs(n - desired) <= 1 + 1e-8
        if t == 0:
            assert counts == spec["extreme_counts"]
        if t == 1:
            assert counts == [spec["target_counts"]] * sites
        assert len(r["metrics"]) == sites * 2
        assert {(m["site"], m["unshifted"]) for m in r["metrics"]} == set(
            itertools.product(range(sites), [False, True])
        )
        for m in r["metrics"]:
            assert m["class_counts"] == [spec["validation_per_class"]] * classes
            assert m["n_test"] == classes * spec["validation_per_class"]
            assert m["split"] == "validation"
            evaluated += 1
        if r["arm"] == "zero-shot":
            if r["seed"] in zero:
                assert zero[r["seed"]] == r["metrics"], "zero-shot changed along sweep"
            zero[r["seed"]] = r["metrics"]
        path = (out / r["state"]).resolve()
        assert path.is_relative_to(out), "checkpoint outside output"
        if str(path) not in state_hashes:
            state_hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        assert state_hashes[str(path)] == r["state_sha256"], "checkpoint changed"
    assert not list(out.rglob("TEST_STARTED"))
    return dict(
        model_records=len(records),
        evaluation_records=evaluated,
        unique_checkpoints_verified=len(state_hashes),
        site_counts_preserved=True,
        pooled_counts_preserved=True,
        evaluation_counts=[spec["validation_per_class"]] * classes,
        zero_shot_identical_across_levels=True,
        test_access=False,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    result = audit(args.directory)
    (args.directory / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
