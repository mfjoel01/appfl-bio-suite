"""Run a bounded MultiSuSiE pilot on an existing synthetic pipeline configuration.

This standalone runner invokes the existing site stage locally, as the SuSiEx
standalone driver does. It does not dispatch remote APPFL tasks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from .fedfm.fed_fine_mapping import SiteIndex, reference_variants, site_aggregates
from .fedfm.fine_mapping import plan_columns
from .fedfm.utils import load_config
from .multisusie import MultiSuSiEOptions, run_multisusie


def run_pilot(config: Path, out: Path, options: MultiSuSiEOptions, limit=1, locus_id=None):
    if limit < 1:
        raise ValueError("limit must be positive")
    cfg = load_config(config)
    loci = pd.read_csv(cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t")
    manifest = pd.read_csv(cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv", sep="\t")
    manifest = manifest[manifest.replicate < cfg.architecture.replicates]
    if locus_id is not None:
        loci = loci[loci.locus_id == locus_id]
    jobs = (
        manifest[manifest.locus_id.isin(loci.locus_id)]
        .sort_values(["locus_id", "architecture_id", "replicate"])
        .head(limit)
    )
    if jobs.empty or jobs.duplicated(["locus_id", "architecture_id", "replicate"]).any():
        raise ValueError("No matching instances, or duplicate manifest keys")
    out.mkdir(parents=True, exist_ok=False)
    (out / "run.json").write_text(
        json.dumps(
            {
                "config": str(config.resolve()),
                "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
                "options": asdict(options),
                "n_instances": len(jobs),
                "execution": "local site stages; coordinator inference; synthetic data only",
            },
            indent=2,
        )
    )
    planned = plan_columns(
        {s: dict(sc.composition) for s, sc in cfg.sites.items()},
        cfg.superpopulations,
        cfg.fine_mapping.min_gwas_n,
    )
    pops = [c.name for c in planned]
    expected_n = {c.name: c.n for c in planned}
    reference = reference_variants(cfg)
    indexes = {site: SiteIndex.load(cfg, site, reference) for site in cfg.sites}
    results = []
    for lid, group in jobs.groupby("locus_id", sort=True):
        locus = loci.set_index("locus_id", drop=False).loc[lid]
        instances = [f"{lid}_{r.architecture_id}_rep{r.replicate}" for r in group.itertuples()]
        geno, pheno = [], []
        for site in cfg.sites:
            g, p = site_aggregates(cfg, site, locus, pops, instances, index=indexes[site])
            geno.extend(g)
            pheno.extend(p)
        actual_n = {pop: sum(g.n for g in geno if g.pop == pop) for pop in pops}
        if actual_n != expected_n:
            raise ValueError("Observed ancestry sample sizes differ from planned cohorts")
        for row, inst in zip(group.itertuples(), instances, strict=True):
            fit = run_multisusie(
                geno, [p for p in pheno if p.instance == inst], out / inst, options
            )
            # Ground truth is used only after fitting, never to set rho or select variants.
            truth = set(str(row.causal_snp_ids).split(","))
            variants = pd.read_csv(out / inst / "variants.tsv", sep="\t")
            cs = pd.read_csv(out / inst / "credible_sets.tsv", sep="\t")
            metrics = {
                "instance": inst,
                "locus_id": lid,
                "architecture_id": row.architecture_id,
                "replicate": row.replicate,
                "causal_mode": getattr(row, "causal_mode", "shared"),
                "fit_status": fit["fit_status"],
                "converged": fit["converged"],
                "n_causal": len(truth),
                "n_eligible_causal": len(truth & set(variants.snp_id)),
                "n_credible_sets": fit["n_credible_sets"],
                "runtime_s": fit["runtime_s"],
            }
            if fit["converged"]:
                captured = truth & set(cs.snp_id)
                metrics.update(
                    n_causal_captured=len(captured),
                    any_causal_captured=bool(captured),
                    n_cs_containing_causal=sum(
                        bool(truth & set(g.snp_id)) for _, g in cs.groupby("component")
                    ),
                )
                for threshold in (0.5, 0.9, 0.95):
                    selected = set(variants.loc[variants.pip > threshold, "snp_id"])
                    metrics[f"n_true_pip{int(100 * threshold)}"] = len(selected & truth)
                    metrics[f"n_false_pip{int(100 * threshold)}"] = len(selected - truth)
            results.append(metrics)
            pd.DataFrame(results).to_csv(out / "results.tsv", sep="\t", index=False)
    if not all(r["converged"] for r in results):
        raise RuntimeError("Pilot contains nonconverged fits; inspect fit.json and results.tsv")
    return pd.DataFrame(results)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=1)
    parser.add_argument("--locus")
    parser.add_argument("--rho", type=float, default=0.75)
    parser.add_argument("--maf", type=float, default=0.005)
    parser.add_argument("--min-mac", type=float, default=0)
    parser.add_argument("--n-signals", type=int, default=10)
    parser.add_argument("--max-iter", type=int, default=1000)
    parser.add_argument("--tol", type=float, default=1e-6)
    parser.add_argument("--coverage", type=float, default=0.95)
    parser.add_argument("--min-abs-corr", type=float, default=0)
    args = parser.parse_args(argv)
    options = MultiSuSiEOptions(
        **{key: getattr(args, key) for key in MultiSuSiEOptions.__dataclass_fields__}
    )
    run_pilot(args.config, args.out, options, args.limit, args.locus)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
