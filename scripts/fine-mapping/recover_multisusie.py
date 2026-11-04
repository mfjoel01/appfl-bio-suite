"""Recovery overlay: preserve frozen source and successful fits; retry capped fits.

The only estimator change is an iteration ceiling of 10,000 after failure at 1,000.
Both use the same convergence tolerance, starting values and model. Numerical PIP
parity remains mandatory. Exact CS membership is a separately reported endpoint;
strict_pass is false for CS differences even when numerical PIP parity passes.
"""

import argparse
import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from appfl_bio_suite.experiments.fine_mapping import multisusie as backend
from appfl_bio_suite.experiments.fine_mapping import multisusie_pipeline as pipeline
from appfl_bio_suite.experiments.fine_mapping.multisusie_compare import compare_fits

BASE_FIT = backend.fit_multisusie
BASE_RUN = backend.run_multisusie


def run_with_trace(geno, pheno, destination, options):
    captured = []

    def capture(inputs, opts):
        fit = BASE_FIT(inputs, opts)
        captured.append(fit)
        return fit

    backend.fit_multisusie = capture
    try:
        record = BASE_RUN(geno, pheno, destination, options)
    finally:
        backend.fit_multisusie = BASE_FIT
    if captured:
        fit = captured[0]
        np.savez_compressed(destination / "iteration_trace.npz", elbo=fit.elbo)
    return record


def recover_fit(geno, pheno, destination, options):
    record = run_with_trace(geno, pheno, destination, options)
    if record["converged"]:
        return record
    archive = destination.with_name(destination.name + f".iteration-cap-{time.time_ns()}")
    destination.rename(archive)
    record = run_with_trace(geno, pheno, destination, replace(options, max_iter=10000))
    record["recovery"] = dict(
        previous_attempt=str(archive), change="max_iter 1000 to 10000; same model and tolerance"
    )
    pipeline.write(destination / "fit.json", record)
    return record


def recovery_parity(federated, central, atol=1e-7, rtol=1e-5):
    records = [pipeline.read(p / "fit.json") for p in [federated, central]]
    opts = [{k: v for k, v in r["options"].items() if k != "max_iter"} for r in records]
    frames = [pd.read_csv(p / "variants.tsv", sep="\t") for p in [federated, central]]
    a, b = frames
    cols = [c for c in a if c != "pip"]
    aligned = set(a.columns) == set(b.columns) and a[cols].equals(b[cols])

    def sets(path):
        f = pd.read_csv(path / "credible_sets.tsv", sep="\t")
        return sorted(tuple(sorted(g.snp_id)) for _, g in f.groupby("component"))

    budgets = [r["options"]["max_iter"] for r in records]
    passed = (
        opts[0] == opts[1]
        and set(budgets) <= {1000, 10000}
        and aligned
        and all(r["converged"] for r in records)
        and records[0]["populations"] == records[1]["populations"]
        and records[0]["sample_sizes"] == records[1]["sample_sizes"]
        and np.allclose(a.pip, b.pip, atol=atol, rtol=rtol)
    )
    exact_cs = sets(federated) == sets(central)
    return {
        "pass": bool(passed),
        "acceptance_criterion": "converged numerical PIP parity; CS membership reported separately",
        "strict_pass": bool(passed and exact_cs),
        "credible_sets_exact": exact_cs,
        "credible_sets_federated": sets(federated) if not exact_cs else None,
        "credible_sets_centralized": sets(central) if not exact_cs else None,
        "variants_aligned": aligned,
        "max_absolute_pip_difference": float(np.max(np.abs(a.pip - b.pip))) if aligned else None,
        "atol": atol,
        "rtol": rtol,
        "iteration_budgets": budgets,
        "comparison": "identical model/tolerance; both converged; allowed iteration ceilings 1000/10000",
    }


def work(root, task):
    pipeline.run_multisusie = recover_fit
    pipeline.compare_fits = recovery_parity
    return pipeline.work_item(root, task)


def diagnostic(root, recovery, cases=None):
    cases = (
        cases
        if cases is not None
        else [
            ("benchmark", "federation_50k_seed3", "L0071", "L0071_ncsl2_h2-0.005_rg1_rep1"),
            (
                "development",
                "afr_to_anl_eur_2500_seed20260930",
                "L0014",
                "L0014_ncsl1_h2-0.0005_rg1_rep4",
            ),
        ]
    )
    for phase, name, lid, inst in cases:
        cfg, loci, jobs = pipeline.jobs_for(root / f"configs/{phase}/{name}.yaml")
        group = jobs[jobs.locus_id == lid]
        instances = [
            f"{r.locus_id}_{r.architecture_id}_rep{r.replicate}" for r in group.itertuples()
        ]
        locus = loci.set_index("locus_id", drop=False).loc[lid]
        columns = pipeline.plan_columns(
            {s: dict(c.composition) for s, c in cfg.sites.items()},
            cfg.superpopulations,
            cfg.fine_mapping.min_gwas_n,
        )
        pops = sorted(c.name for c in columns)
        ref = pipeline.reference_variants(cfg)
        geno, pheno = [], []
        for site in cfg.sites:
            g, p = pipeline.site_aggregates(
                cfg, site, locus, pops, instances, index=pipeline.SiteIndex.load(cfg, site, ref)
            )
            geno.extend(g)
            pheno.extend(p)
        pheno = [p for p in pheno if p.instance == inst]
        original = pipeline.read(root / phase / name / lid / inst / "multisusie/fit.json")
        inputs = backend.prepare_inputs(geno, pheno, pipeline.MultiSuSiEOptions())
        fingerprint = hashlib.sha256(inputs.variants.to_csv(index=False).encode())
        for pop, n, x, c, w in zip(
            inputs.populations, inputs.sample_sizes, inputs.xtx, inputs.xty, inputs.yty, strict=True
        ):
            fingerprint.update(json.dumps([pop, n, w]).encode())
            fingerprint.update(x.tobytes())
            fingerprint.update(c.tobytes())
        if fingerprint.hexdigest() != original["moments_sha256"]:
            raise ValueError("Diagnostic staging does not reproduce original moments")
        print(f"Verified original batched moments: {inst}", flush=True)
        target = recovery / "diagnostics" / phase / inst
        if (target / "fit.json").exists():
            record = pipeline.read(target / "fit.json")
        else:
            record = run_with_trace(geno, pheno, target, pipeline.MultiSuSiEOptions(max_iter=10000))
        if original["moments_sha256"] != record["moments_sha256"]:
            raise ValueError("Diagnostic changed input moments")
        print(json.dumps(record), flush=True)
        if not record["converged"]:
            raise RuntimeError(f"Extended fit did not converge: {inst}")
    pipeline.write(recovery / "diagnostics.accepted.json", {"pass": True, "n_cases": len(cases)})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=["diagnostic", "run"])
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--recovery", type=Path, required=True)
    p.add_argument("--phase", default="benchmark")
    p.add_argument("--index", type=int, default=0)
    p.add_argument("--workers", type=int, default=2)
    args = p.parse_args()
    pipeline.verify(args.root)
    if args.stage == "diagnostic":
        diagnostic(args.root, args.recovery)
        return
    if args.phase == "evaluation":
        for name, digest in pipeline.read(args.root / "selection.lock.json").items():
            if pipeline.sha(name) != digest:
                raise ValueError(f"Frozen selection changed: {name}")
        tasks = pipeline.read(args.root / "evaluation_tasks.json")[args.index :: 64]
    else:
        tasks = pipeline.read(args.recovery / f"{args.phase}.batches.json")[args.index]
    failures = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        fs = {pool.submit(work, args.root, t): t for t in tasks}
        for f in as_completed(fs):
            try:
                print("Completed", f.result(), flush=True)
            except Exception as exc:
                failures.append({"task": fs[f], "error": repr(exc)})
                print("FAILED", fs[f], repr(exc), flush=True)
    pipeline.write(args.recovery / f"{args.phase}.failures.{args.index}.json", failures)
    if failures:
        raise RuntimeError(f"{len(failures)} tasks failed")


if __name__ == "__main__":
    main()
