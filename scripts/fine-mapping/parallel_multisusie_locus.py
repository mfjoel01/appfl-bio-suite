"""Prefit independent instances in parallel, then use the original locus validator.

All phenotypes are staged together exactly as in the frozen pipeline. Forked workers
share read-only staged moments. Each worker owns a distinct fit directory. Use only
when no other job is writing this locus. Scientific settings are unchanged.
"""

import argparse
import multiprocessing
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import recover_multisusie as recovery

p = recovery.pipeline
_GENO = None
_PHENO = None
_OUT = None


def prefit(instance):
    dest = _OUT / instance / "multisusie"
    if dest.exists():
        ok = (dest / "fit.json").exists() and p.read(dest / "fit.json").get("converged")
        if ok and all(
            (dest / f).exists() for f in ["variants.tsv", "credible_sets.tsv", "posterior.npz"]
        ):
            return instance
        dest.rename(dest.with_name(dest.name + f".interrupted-{time.time_ns()}"))
    record = recovery.recover_fit(
        _GENO, [v for v in _PHENO if v.instance == instance], dest, p.MultiSuSiEOptions()
    )
    if not record["converged"]:
        raise RuntimeError(f"Nonconverged fit: {dest}")
    return instance


def complete_locus(root, task, workers):
    global _GENO, _PHENO, _OUT
    if task["central_check"]:
        raise ValueError("This runner is for evaluation without centralized checks")
    cfg, loci, jobs = p.jobs_for(task["config"])
    jobs = jobs[jobs.locus_id == task["locus_id"]]
    instances = [f"{r.locus_id}_{r.architecture_id}_rep{r.replicate}" for r in jobs.itertuples()]
    _OUT = root / task["phase"] / task["name"] / task["locus_id"]
    provenance = {
        str(f): p.sha(f)
        for site in cfg.sites
        for inst in instances
        for f in [cfg.resolved_path("ground_truth_dir") / "phenotypes" / site / f"{inst}.pheno"]
    }
    saved = _OUT / "phenotypes.sha256.json"
    if saved.exists() and p.read(saved) != provenance:
        raise ValueError("Phenotypes changed during restart")
    locus = loci.set_index("locus_id", drop=False).loc[task["locus_id"]]
    columns = p.plan_columns(
        {s: dict(c.composition) for s, c in cfg.sites.items()},
        cfg.superpopulations,
        cfg.fine_mapping.min_gwas_n,
    )
    pops = sorted(c.name for c in columns)
    ref = p.reference_variants(cfg)
    _GENO, _PHENO = [], []
    for site in cfg.sites:
        g, h = p.site_aggregates(
            cfg, site, locus, pops, instances, index=p.SiteIndex.load(cfg, site, ref)
        )
        _GENO.extend(g)
        _PHENO.extend(h)
    if {pop: sum(g.n for g in _GENO if g.pop == pop) for pop in pops} != {
        c.name: c.n for c in columns
    }:
        raise ValueError("Observed sample sizes differ from design")
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=multiprocessing.get_context("fork")
    ) as pool:
        futures = [pool.submit(prefit, inst) for inst in instances]
        for future in as_completed(futures):
            print("Fit checkpoint", future.result(), flush=True)
    # Original validation, scoring, phenotype hashes and acceptance checks remain authoritative.
    _GENO, _PHENO = None, None
    result = recovery.work(root, task)
    print("Accepted", result, flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    p.verify(args.root)
    for name, digest in p.read(args.root / "selection.lock.json").items():
        if p.sha(name) != digest:
            raise ValueError(f"Frozen selection changed: {name}")
    tasks = [
        t
        for t in p.read(args.root / "evaluation_tasks.json")
        if t["name"] == "federation" and t["locus_id"] == "L0063"
    ]
    if len(tasks) != 1:
        raise ValueError("Expected exactly one final locus")
    complete_locus(args.root, tasks[0], args.workers)
