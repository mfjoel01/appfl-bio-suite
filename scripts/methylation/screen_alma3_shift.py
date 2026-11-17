#!/usr/bin/env python3
"""Screen synthetic shift severity in the isolated ALMA3 worker environment.

Uses only the first configured seed and validation split. This is generator
selection, not an independent evaluation. Run with NUMPY_MADVISE_HUGEPAGE=0.
"""

import argparse
import json
from dataclasses import replace
from pathlib import Path

from appfl_bio_suite.experiments.methylation.alma3.worker import Worker


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--release", required=True, type=Path)
    parser.add_argument("--member", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--parameter", choices=["batch_shift", "signature_rotation_degrees"], default="batch_shift"
    )
    parser.add_argument("--shifts", nargs="+", type=float, default=[0, 2, 4, 8, 12])
    args = parser.parse_args()
    out = args.out.resolve()
    if out.is_relative_to(args.release.resolve()):
        raise ValueError("screen outputs must be outside the immutable release")
    config = json.loads(args.config.read_text())
    worker = Worker(
        dict(
            config=config,
            fixture=False,
            release=str(args.release.resolve()),
            member=args.member,
            device="cpu",
            allow_test=False,
        )
    )
    for shift in args.shifts:
        replace(worker.config, **{args.parameter: shift}).validate()
    out.mkdir(parents=True, exist_ok=False)
    identity = worker.handle(dict(op="identity"))
    (out / "LICENSE.ALMA3").write_text(identity.pop("license"))
    (out / "ATTRIBUTION.txt").write_text(
        "Analysis performed using ALMA3. Based on ALMA3 Version 3.1.2.\n"
    )
    (out / "identity.json").write_text(json.dumps(identity, indent=2))
    (out / "config.json").write_text(json.dumps(config, indent=2))
    rows = []
    for shift in args.shifts:
        worker.config = replace(worker.config, **{args.parameter: shift})
        worker.cache.clear()
        worker.feature_cache.clear()
        worker.input_cache.clear()
        worker.input_cache_bytes = 0
        for site in range(len(worker.config.site_sizes)):
            for unshifted in [True, False] if shift == args.shifts[0] else [False]:
                result = worker.handle(
                    dict(
                        op="evaluate",
                        seed=worker.config.seeds[0],
                        depth="heads",
                        site=site,
                        split="validation",
                        coverage=max(worker.config.coverage),
                        unshifted=unshifted,
                        thresholds=[1.01] * identity["levels"],
                    )
                )
                rows.append(
                    dict(
                        parameter=args.parameter,
                        shift=shift,
                        site=site,
                        unshifted=unshifted,
                        metrics=result,
                    )
                )
                print(
                    shift,
                    site,
                    unshifted,
                    result["balanced_accuracy"],
                    result.get("support_restricted_balanced_accuracy"),
                    flush=True,
                )
        (out / "screen.json").write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
