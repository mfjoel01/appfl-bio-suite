#!/usr/bin/env python3
"""End-to-end real-data validation pipeline for FedFM on MAGE + 1000 Genomes NYGC 30x GRCh38."""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--data-base",
        type=Path,
        default=Path("local/data/mage"),
        help="Base directory for MAGE staged inputs, cache, and site bundles.",
    )
    p.add_argument(
        "--out-root",
        type=Path,
        default=Path("local/output/mage-validation"),
        help="Directory for validation arms, summary tables, and publication figures.",
    )
    p.add_argument("--n-per-stratum", type=int, default=4)
    p.add_argument("--window-half-bp", type=int, default=250_000)
    p.add_argument("--min-gap-bp", type=int, default=1_000_000)
    p.add_argument("--maf-build", type=float, default=0.01)
    p.add_argument("--maf-analysis", type=float, default=0.005)
    p.add_argument("--n-signals", type=int, default=10)
    p.add_argument("--max-iter", type=int, default=1000)
    p.add_argument("--tol", type=float, default=1e-6)
    p.add_argument("--multisusie-max-iter", type=int, default=1000)
    p.add_argument("--multisusie-tol", type=float, default=1e-5)
    p.add_argument("--multisusie-min-abs-corr", type=float, default=0.5)
    p.add_argument("--max-workers", type=int, default=4)
    p.add_argument("--locus-limit", type=int, default=None)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--no-swap-control", action="store_true")
    p.add_argument("--stage", action="store_true", help="Stage MAGE raw files and build cache.")
    p.add_argument(
        "--build",
        action="store_true",
        help="Build clean (and swap) 3-site federations.",
    )
    p.add_argument("--loopback", action="store_true", help="Run APPFL loopback driver check.")
    p.add_argument("--run", action="store_true", help="Run multi-arm validation benchmark.")
    p.add_argument("--plot", action="store_true", help="Render PNG + PDF validation figures.")
    p.add_argument(
        "--replication",
        action="store_true",
        help=(
            "Build raw-expression federations and run the replication design: pooled SuSiE "
            "against the published MAGE credible sets, single sites against pooled, and "
            "federated against pooled, all with federated FWL covariate adjustment."
        ),
    )
    p.add_argument(
        "--replication-out",
        type=Path,
        default=Path("local/output/mage-replication"),
        help="Output directory for the replication design.",
    )
    p.add_argument("--no-susiex", action="store_true", help="Skip SuSiEx replication arms.")
    p.add_argument(
        "--no-wide-window",
        action="store_true",
        help="Skip the sensitivity run at the published ±1 Mb window.",
    )
    p.add_argument(
        "--no-site-covariates",
        action="store_true",
        help="Skip the sensitivity run with covariates each site derives from its own data.",
    )
    p.add_argument("--all", action="store_true", help="Run all stages end to end.")
    return p


def _write_loopback_parity(loop_results: Path, fed_results: Path, out_path: Path) -> dict | None:
    """Compare APPFL serial loopback output against the federated SuSiEx arm on overlapping loci."""
    if not (loop_results.is_file() and fed_results.is_file()):
        return None
    import pandas as pd

    lp = pd.read_csv(loop_results, sep="\t").set_index("locus_id")
    fd = pd.read_csv(fed_results, sep="\t").set_index("locus_id")
    common = sorted(set(lp.index) & set(fd.index))
    if not common:
        return None

    def cs_memberships(results, lid):
        paths = sorted((results.parent / "instances").glob(f"{lid}_*/credible_sets.tsv"))
        if not paths:
            return None
        sets = []
        for path in paths:
            frame = pd.read_csv(path, sep="\t")
            sets.extend(
                tuple(sorted(g["snp_id"].astype(str))) for _, g in frame.groupby("component")
            )
        return sorted(sets)

    records = []
    max_top_pip_diff = 0.0
    for lid in common:
        row_l = lp.loc[lid]
        row_f = fd.loc[lid]
        pip_l = float(row_l["top_pip"]) if pd.notna(row_l["top_pip"]) else float("nan")
        pip_f = float(row_f["top_pip"]) if pd.notna(row_f["top_pip"]) else float("nan")
        diff = abs(pip_l - pip_f)
        if pd.notna(diff):
            max_top_pip_diff = max(max_top_pip_diff, diff)
        cs_l = cs_memberships(loop_results, lid)
        cs_f = cs_memberships(fed_results, lid)
        records.append(
            {
                "locus_id": lid,
                "loopback_cs_sizes": str(row_l["cs_sizes_json"]),
                "federated_cs_sizes": str(row_f["cs_sizes_json"]),
                "loopback_top_pip": pip_l,
                "federated_top_pip": pip_f,
                "top_pip_abs_diff": diff,
                "exact_cs_match": cs_l is not None and cs_f is not None and cs_l == cs_f,
            }
        )
    summary = {
        "n_common_loci": len(common),
        "all_exact_cs_match": all(r["exact_cs_match"] for r in records),
        "max_top_pip_abs_diff": max_top_pip_diff,
        "loci": records,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    args = _build_parser().parse_args(argv)
    if not any(
        [args.stage, args.build, args.loopback, args.run, args.plot, args.replication, args.all]
    ):
        args.all = True

    from appfl_bio_suite.experiments.fine_mapping.real_data.benchmark import run_validation_suite
    from appfl_bio_suite.experiments.fine_mapping.real_data.federation_builder import (
        build_mage_federation,
    )
    from appfl_bio_suite.experiments.fine_mapping.real_data.figures_mage import (
        render_all_mage_figures,
    )
    from appfl_bio_suite.experiments.fine_mapping.real_data.mage import (
        prepare_mage_cache,
        stage_mage_data,
    )

    data_base = args.data_base.resolve()
    raw_dir = data_base / "raw"
    cache_dir = data_base / "cache"
    clean_root = data_base / "federation_clean"
    swap_root = data_base / "federation_with_swap"
    out_root = args.out_root.resolve()

    if args.stage or args.all:
        stage_mage_data(raw_dir)
        prepare_mage_cache(raw_dir, cache_dir, exclude_swapped=True, overwrite=args.overwrite)
        if not args.no_swap_control:
            prepare_mage_cache(raw_dir, cache_dir, exclude_swapped=False, overwrite=args.overwrite)

    if args.build or args.all:
        clean_manifest = build_mage_federation(
            out_root=clean_root,
            raw_dir=raw_dir,
            cache_dir=cache_dir,
            n_per_stratum=args.n_per_stratum,
            window_half_bp=args.window_half_bp,
            min_gap_bp=args.min_gap_bp,
            maf=args.maf_build,
            exclude_swapped=True,
            overwrite=args.overwrite,
        )
        print("Clean federation built:", json.dumps(clean_manifest["sites"]))
        if not args.no_swap_control:
            swap_manifest = build_mage_federation(
                out_root=swap_root,
                raw_dir=raw_dir,
                cache_dir=cache_dir,
                n_per_stratum=args.n_per_stratum,
                window_half_bp=args.window_half_bp,
                min_gap_bp=args.min_gap_bp,
                maf=args.maf_build,
                exclude_swapped=False,
                overwrite=args.overwrite,
            )
            print("Swap-control federation built:", json.dumps(swap_manifest["sites"]))

    if args.loopback or args.all:
        from appfl_bio_suite.core.launch import launch
        from appfl_bio_suite.core.loopback import loopback_federation

        loop_out = out_root / "appfl_loopback"
        loop_out.mkdir(parents=True, exist_ok=True)
        fed = loopback_federation(
            experiment="fine-mapping",
            data_root=clean_root,
            out_root=loop_out / "sites",
        )
        exp = fed.experiment("fine-mapping")
        object.__setattr__(exp, "locus_limit", 1 if args.locus_limit is None else args.locus_limit)
        object.__setattr__(exp, "fine_mapping_output_dir", str((loop_out / "results").resolve()))
        object.__setattr__(exp, "maf", args.maf_analysis)
        object.__setattr__(exp, "n_signals", args.n_signals)
        object.__setattr__(exp, "max_iter", args.max_iter)
        object.__setattr__(exp, "tol", args.tol)
        rc = launch(
            fed,
            "fine-mapping",
            driver="serial",
            out_dir=loop_out / "configs",
        )
        if rc != 0:
            raise RuntimeError(f"APPFL loopback serial check failed with exit code {rc}")
        print("APPFL loopback serial check completed ->", loop_out / "results")

    if args.run or args.all:
        report = run_validation_suite(
            clean_root=clean_root,
            out_root=out_root,
            swap_root=None if args.no_swap_control else swap_root,
            maf=args.maf_analysis,
            n_signals=args.n_signals,
            max_iter=args.max_iter,
            tol=args.tol,
            multisusie_max_iter=args.multisusie_max_iter,
            multisusie_tol=args.multisusie_tol,
            multisusie_min_abs_corr=args.multisusie_min_abs_corr,
            max_workers=args.max_workers,
            locus_limit=args.locus_limit,
            skip_existing=not args.overwrite,
        )
        print(
            "Validation report complete:",
            json.dumps(report.get("parity_summary", []), indent=2),
        )

    loop_parity = _write_loopback_parity(
        out_root / "appfl_loopback" / "results" / "data" / "fed_fm_results.tsv",
        out_root / "arms" / "federated_susiex" / "data" / "fed_fm_results.tsv",
        out_root / "reports" / "appfl_loopback_parity.json",
    )
    if loop_parity is not None:
        print("APPFL loopback parity:", json.dumps(loop_parity))

    if args.plot or args.all:
        figs = render_all_mage_figures(clean_root, out_root)
        print(f"Rendered {len(figs)} figure files under {out_root / 'figures'}")

    if args.replication or args.all:
        from appfl_bio_suite.experiments.fine_mapping.real_data.benchmark import (
            run_replication_suite,
        )

        fwl_roots = {True: data_base / "federation_clean_fwl"}
        if not args.no_swap_control:
            fwl_roots[False] = data_base / "federation_with_swap_fwl"
        common = dict(
            raw_dir=raw_dir,
            cache_dir=cache_dir,
            n_per_stratum=args.n_per_stratum,
            min_gap_bp=args.min_gap_bp,
            maf=args.maf_build,
            overwrite=args.overwrite,
            phenotype="raw",
            selection_window_half_bp=args.window_half_bp,
        )
        tags = {True: "730", False: "731"}
        wide_roots, site_local_roots = {}, {}
        for exclude_swapped, root in fwl_roots.items():
            build_mage_federation(
                out_root=root,
                window_half_bp=args.window_half_bp,
                exclude_swapped=exclude_swapped,
                **common,
            )
            tag = tags[exclude_swapped]
            if not args.no_wide_window:
                wide_roots[tag] = root.with_name(f"{root.name}_1mb")
                build_mage_federation(
                    out_root=wide_roots[tag],
                    window_half_bp=1_000_000,
                    exclude_swapped=exclude_swapped,
                    **common,
                )
            if not args.no_site_covariates:
                site_local_roots[tag] = root.with_name(f"{root.name}_sitecov")
                build_mage_federation(
                    out_root=site_local_roots[tag],
                    window_half_bp=args.window_half_bp,
                    exclude_swapped=exclude_swapped,
                    covariate_source="site_local",
                    **common,
                )
        report = run_replication_suite(
            clean_root=fwl_roots[True],
            swap_root=fwl_roots.get(False),
            out_root=args.replication_out.resolve(),
            maf=args.maf_analysis,
            n_signals=args.n_signals,
            max_workers=args.max_workers,
            run_susiex=not args.no_susiex,
            locus_limit=args.locus_limit,
            wide_roots=wide_roots or None,
            site_local_roots=site_local_roots or None,
        )
        print("Replication vs MAGE:", json.dumps(report["vs_mage"], indent=2, default=str))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
