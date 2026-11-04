"""Paired three-method pilot with an independent centralized genotype comparator."""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from .fedfm.fed_fine_mapping import (
    GenoAggregate,
    PhenoAggregate,
    SiteIndex,
    build_columns,
    pool_pheno,
    pooled_sumstats,
    reference_variants,
    site_aggregates,
    write_sumstats,
)
from .fedfm.fine_mapping import plan_columns, run_susiex
from .fedfm.utils import load_config, read_bed_variants, read_bim, read_fam, read_site_manifest
from .inference import parse_outputs, pip_by_snp, read_table
from .multisusie import MultiSuSiEOptions, run_multisusie


def centralized_blocks(cfg, locus, pops, instances):
    """Read the source cohort directly and form moments over pooled individuals.

    This deliberately does not add site moments. Missing variants are excluded over
    the entire enrolled ancestry, just as the federated complete-variant intersection.
    """
    prefix = cfg.resolved_path("hapnest_dir") / f"chr{cfg.chromosome}"
    bim, fam = read_bim(prefix.with_suffix(".bim")), read_fam(prefix.with_suffix(".fam"))
    variants = bim.loc[
        (bim.chrom.astype(str) == str(locus.chrom)) & bim.bp.between(locus.start_bp, locus.end_bp)
    ]
    x = read_bed_variants(prefix, variants.index.to_numpy(), len(fam))
    variants = variants[["chrom", "snp_id", "bp", "a1", "a2"]].reset_index(drop=True)
    enrolled = pd.concat(
        [read_site_manifest(cfg.site_dir(site) / f"{site}_manifest.tsv") for site in cfg.sites],
        ignore_index=True,
    )
    keys = ["FID", "IID"]
    if enrolled.duplicated(keys).any() or fam.duplicated(keys).any():
        raise ValueError("Duplicated pooled individual identifiers")
    source_index = pd.MultiIndex.from_frame(fam[keys])
    rows = source_index.get_indexer(pd.MultiIndex.from_frame(enrolled[keys]))
    if (rows < 0).any():
        raise ValueError("An enrolled individual is absent from the source cohort")
    ys = []
    for inst in instances:
        frames = [
            pd.read_csv(
                cfg.resolved_path("ground_truth_dir") / "phenotypes" / site / f"{inst}.pheno",
                sep=r"\s+",
                dtype={"FID": str, "IID": str},
            )
            for site in cfg.sites
        ]
        frame = pd.concat(frames, ignore_index=True)
        if frame.duplicated(keys).any():
            raise ValueError("Duplicated pooled phenotype identifiers")
        # Matched-N arms may reuse phenotype files for the full original site.
        y = frame.set_index(keys).reindex(pd.MultiIndex.from_frame(enrolled[keys])).y.to_numpy()
        if not np.isfinite(y).all():
            raise ValueError("Missing or nonfinite pooled phenotypes")
        ys.append(y)
    geno, pheno = [], []
    for pop in pops:
        members = np.flatnonzero(enrolled.superpopulation.to_numpy() == pop)
        xp = x[rows[members]]
        complete = np.isfinite(xp).all(axis=0)
        xp = xp[:, complete].astype(np.float64)
        vp = variants.loc[complete].reset_index(drop=True)
        n = len(xp)
        geno.append(GenoAggregate("central", pop, str(locus.locus_id), vp, xp.T @ xp, xp.sum(0), n))
        for inst, y in zip(instances, ys, strict=True):
            yp = y[members]
            pheno.append(
                PhenoAggregate(
                    "central",
                    pop,
                    inst,
                    vp.snp_id.to_numpy(),
                    xp.T @ yp,
                    float(yp.sum()),
                    float(yp @ yp),
                    n,
                )
            )
    return geno, pheno


def compare_fits(federated: Path, central: Path, atol=1e-7, rtol=1e-5):
    """Fail on missing variants, eligibility, nonconvergence, PIPs or set membership."""
    a, b = [pd.read_csv(p / "variants.tsv", sep="\t") for p in (federated, central)]
    records = [json.loads((p / "fit.json").read_text()) for p in (federated, central)]
    columns = [c for c in a if c != "pip"]
    aligned = set(a.columns) == set(b.columns) and a[columns].equals(b[columns])
    delta = float(np.max(np.abs(a.pip - b.pip))) if aligned and len(a) else None

    def sets(path):
        cs = pd.read_csv(path / "credible_sets.tsv", sep="\t")
        return sorted(tuple(sorted(g.snp_id)) for _, g in cs.groupby("component"))

    passed = bool(
        aligned
        and all(r["converged"] for r in records)
        and records[0]["populations"] == records[1]["populations"]
        and records[0]["sample_sizes"] == records[1]["sample_sizes"]
        and records[0]["options"] == records[1]["options"]
        and np.allclose(a.pip, b.pip, atol=atol, rtol=rtol)
        and sets(federated) == sets(central)
    )
    return {
        "pass": passed,
        "variants_aligned": aligned,
        "max_absolute_pip_difference": delta,
        "atol": atol,
        "rtol": rtol,
    }


def score(variants, sets, truth, converged):
    """Common metric definitions; failed fits never receive discovery scores."""
    truth = set(truth)
    result = {"n_causal": len(truth), "n_eligible_causal": len(truth & set(variants))}
    if not converged:
        return result
    captured = truth & set().union(*sets)
    result.update(
        n_credible_sets=len(sets),
        n_causal_captured=len(captured),
        any_causal_captured=bool(captured),
        n_cs_containing_causal=sum(bool(s & truth) for s in sets),
        total_cs_snps=sum(map(len, sets)),
    )
    for threshold in (0.5, 0.9, 0.95):
        selected = {v for v, pip in variants.items() if pip > threshold}
        result[f"n_true_pip{int(100 * threshold)}"] = len(selected & truth)
        result[f"n_false_pip{int(100 * threshold)}"] = len(selected - truth)
    return result


def run_comparison(config, panel, out, susiex, options, central_check=True):
    """Run one arm on explicit panel keys. The caller supplies a fresh output path."""
    if options.min_mac != 0:
        raise ValueError("Paired pilot currently requires min_mac=0 to match SuSiEx eligibility")
    cfg = load_config(config)
    keys = ["locus_id", "architecture_id", "replicate"]
    jobs = pd.read_csv(panel, sep="\t")[keys]
    if jobs.empty or jobs.duplicated(keys).any():
        raise ValueError("Panel must contain unique, nonempty instance keys")
    manifest = pd.read_csv(cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv", sep="\t")
    jobs = jobs.merge(manifest, on=keys, how="left", validate="one_to_one", indicator=True)
    if not jobs._merge.eq("both").all() or (jobs.replicate >= cfg.architecture.replicates).any():
        raise ValueError("Panel includes instances unavailable in this arm")
    loci = pd.read_csv(cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t")
    if not set(jobs.locus_id) <= set(loci.locus_id) or loci.locus_id.duplicated().any():
        raise ValueError("Panel has missing or ambiguous loci")
    out.mkdir(parents=True, exist_ok=False)
    (out / "panel.tsv").write_bytes(Path(panel).read_bytes())
    planned = plan_columns(
        {s: dict(c.composition) for s, c in cfg.sites.items()},
        cfg.superpopulations,
        cfg.fine_mapping.min_gwas_n,
    )
    pops = sorted(c.name for c in planned)
    expected_n = {c.name: c.n for c in planned}
    ref = reference_variants(cfg)
    indexes = {s: SiteIndex.load(cfg, s, ref) for s in cfg.sites}
    results, parities = [], []
    for lid, group in jobs.groupby("locus_id", sort=True):
        locus = loci.set_index("locus_id", drop=False).loc[lid]
        instances = [f"{lid}_{r.architecture_id}_rep{r.replicate}" for r in group.itertuples()]
        geno, pheno = [], []
        start = time.perf_counter()
        for site in cfg.sites:
            g, p = site_aggregates(cfg, site, locus, pops, instances, index=indexes[site])
            geno.extend(g)
            pheno.extend(p)
        if {p: sum(g.n for g in geno if g.pop == p) for p in pops} != expected_n:
            raise ValueError("Planned and observed sample sizes differ")
        aggregation_s = time.perf_counter() - start
        columns = build_columns(geno, pops, out / f"{lid}_ld", options.maf)
        if central_check:
            cg, cp = centralized_blocks(cfg, locus, pops, instances)
        for row, inst in zip(group.itertuples(), instances, strict=True):
            ps = [p for p in pheno if p.instance == inst]
            target = out / inst
            target.mkdir()
            fit = run_multisusie(geno, ps, target / "multisusie", options)
            if central_check:
                run_multisusie(
                    cg, [p for p in cp if p.instance == inst], target / "centralized", options
                )
                parity = compare_fits(target / "multisusie", target / "centralized")
                parities.append({"instance": inst, **parity})
                (out / "parity.json").write_text(json.dumps(parities, indent=2))
                if not parity["pass"]:
                    raise RuntimeError(f"Centralized-federated parity failed: {inst}")
            truth = str(row.causal_snp_ids).split(",")
            base = {
                "instance": inst,
                "locus_id": lid,
                "architecture_id": row.architecture_id,
                "replicate": row.replicate,
                "stratum": locus.stratum,
                "causal_mode": getattr(row, "causal_mode", "shared"),
                "analyzed_n": sum(expected_n.values()),
                "site_stage_s": aggregation_s,
            }
            variants = pd.read_csv(target / "multisusie/variants.tsv", sep="\t")
            cs = pd.read_csv(target / "multisusie/credible_sets.tsv", sep="\t")
            results.append(
                {
                    **base,
                    "method": "multisusie",
                    "fit_status": fit["fit_status"],
                    "runtime_s": fit["runtime_s"],
                    "pip_definition": fit["pip_definition"],
                    **score(
                        dict(zip(variants.snp_id, variants.pip, strict=True)),
                        [set(g.snp_id) for _, g in cs.groupby("component")],
                        truth,
                        fit["converged"],
                    ),
                }
            )
            ss, ld, ns = [], [], []
            for col in columns:
                m = pool_pheno([p for p in ps if p.pop == col.pop], col.variants.snp_id.to_numpy())
                frame = pooled_sumstats(col.variants, col.g_diag, col.u, m.c, m.q, m.w, col.n)
                ss.append(write_sumstats(frame, target / f"{col.pop}.sumstats"))
                ld.append(col.ld_prefix)
                ns.append(col.n)
            for method, purity, pvalue in [("susiex", 0.5, 1e-5), ("susiex_unfiltered", 0, 1)]:
                directory = target / method
                start = time.perf_counter()
                run_susiex(
                    ss,
                    ld,
                    ns,
                    int(locus.chrom),
                    int(locus.start_bp),
                    int(locus.end_bp),
                    directory,
                    "cs",
                    str(susiex),
                    "plink",
                    options.coverage,
                    pvalue,
                    options={
                        "min_purity": purity,
                        "n_signals": options.n_signals,
                        "max_iter": options.max_iter,
                        "tol": options.tol,
                    },
                )
                elapsed = time.perf_counter() - start
                parsed = parse_outputs(directory, "cs", truth)
                snps = read_table(directory / "cs.snp")
                found = pip_by_snp(directory / "cs.snp")
                if set(snps.SNP.astype(str)) != set(variants.snp_id):
                    raise RuntimeError(f"{method}: variant eligibility differs from MultiSuSiE")
                cs = read_table(directory / "cs.cs")
                sets = (
                    [set(g.SNP.astype(str)) for _, g in cs.groupby("CS_ID")]
                    if "CS_ID" in cs
                    else []
                )
                results.append(
                    {
                        **base,
                        "method": method,
                        "fit_status": parsed["fit_status"],
                        "runtime_s": elapsed,
                        "pip_definition": parsed["pip_definition"],
                        **score(found, sets, truth, parsed["converged"]),
                    }
                )
            pd.DataFrame(results).to_csv(out / "results.tsv", sep="\t", index=False)
            print(f"Completed {inst}: MultiSuSiE and two SuSiEx fits", flush=True)
    frame = pd.DataFrame(results)
    if not frame.fit_status.str.startswith("converged_").all():
        raise RuntimeError("Comparison contains failed or nonconverged fits")
    (out / "accepted.json").write_text(
        json.dumps(
            {
                "pass": True,
                "n_instances": len(jobs),
                "n_fits": len(frame),
                "centralized_parity_checked": central_check,
                "options": asdict(options),
                "comparison": "native MultiSuSiE moments; SuSiEx six-digit sumstats and float32 LD",
            },
            indent=2,
        )
    )
    return frame


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--susiex", type=Path, required=True)
    parser.add_argument("--rho", type=float, default=0.75)
    parser.add_argument("--skip-centralized", action="store_true")
    args = parser.parse_args(argv)
    run_comparison(
        args.config,
        args.panel,
        args.out,
        args.susiex,
        MultiSuSiEOptions(rho=args.rho),
        not args.skip_centralized,
    )


if __name__ == "__main__":
    main()
