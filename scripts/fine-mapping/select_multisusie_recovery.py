"""Preserve the frozen selection rule, including the outcome of no eligible candidate."""

import argparse
from pathlib import Path

from appfl_bio_suite.experiments.fine_mapping import multisusie_pipeline as pipeline


def select(root):
    pipeline.verify(root)
    try:
        pipeline.select(root)
        return
    except ValueError as exc:
        if str(exc) != "No candidate passes the prespecified calibration screen":
            raise
    # Complete, converged development results were collected and ranked by the
    # original code. A failed screen is a scientific outcome, not an exception
    # permitting us to pick the highest-power ineligible candidate.
    if (root / "selection.json").exists():
        raise RuntimeError("Cannot replace an existing selection")
    design = pipeline.read(root / "design.json")
    selection = dict(
        name=None,
        solo=None,
        site_counts={},
        ancestry_counts={},
        eligible=False,
        reason="No candidate passed the prespecified calibration screen",
        evaluation_used_for_selection=False,
        development_results_sha256=pipeline.sha(root / "development_results.tsv"),
    )
    pipeline.write(root / "selection.json", selection)
    source = Path(design["source_composition"])
    arms = []
    for name in ["federation", "anl", "covenant", "mbzuai"]:
        cfg = pipeline.freeze_config(
            root, "evaluation", name, source / "evaluation/arms" / f"{name}.yaml"
        )
        arms.append(
            dict(
                name=name,
                config=str(cfg),
                candidate=name,
                sampling_seed=design["sampling_seeds"][0],
                central_check=False,
            )
        )
    pipeline.make_tasks(root, "evaluation", arms)
    files = [
        root / "selection.json",
        root / "evaluation_tasks.json",
        *list((root / "configs/evaluation").glob("*.yaml")),
    ]
    pipeline.write(root / "selection.lock.json", {str(p): pipeline.sha(p) for p in files})


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    args = p.parse_args()
    select(args.root.resolve())
