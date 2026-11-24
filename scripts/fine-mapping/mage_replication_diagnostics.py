#!/usr/bin/env python3
"""Diagnostics behind MAGE_VALIDATION.md that sit outside the replication suite.

Run from the repository root with ``PYTHONPATH=src`` after
``run_mage_validation.py --replication`` has built the federations. Each command writes an
arm directory under ``<out>/diagnostics`` that ``benchmark.score_against_published`` scores.

``covariates VARIANT``
    Pooled SuSiE (N = 731, FWL) with one preprocessing step done per site instead of on the
    pooled cohort, to locate what site-derived covariates cost:
    ``siteINT_publishedcov``, ``pooledINT_site15``, ``siteINT_site30``, ``siteINT_site45``,
    ``pooledINT_pooledPC60``.
``susiex VARIANT``
    SuSiEx with one population (where it should reduce to SuSiE) or five ancestries:
    ``1pop``, ``1pop_nopfilter``, ``1pop_adjusted``, ``1pop_adjusted_nopfilter``,
    ``5anc_adjusted``, ``5anc_adjusted_nopfilter``. ``adjusted`` passes ``susiex_finite_sample_adjust=True``.
``fixed-prior``
    SuSiE with prior estimation off and SuSiEx's prior (the largest squared standardized
    marginal effect, ``data.cpp:482``) instead.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

COVARIATE_VARIANTS = {
    # name: (expression normalization, covariates)
    "siteINT_publishedcov": ("site", "published"),
    "pooledINT_site15": ("pooled", 15),
    "siteINT_site30": ("site", 30),
    "siteINT_site45": ("site", 45),
    "pooledINT_pooledPC60": ("pooled", "pooled60"),
}
SUSIEX_VARIANTS = {
    # name: (single population, pval_thresh, finite-sample adjustment, N = 731 root)
    "1pop": (True, 1e-5, False, True),
    "1pop_nopfilter": (True, 1.0, False, True),
    "1pop_adjusted": (True, 1e-5, True, True),
    "1pop_adjusted_nopfilter": (True, 1.0, True, True),
    "5anc_adjusted": (False, 1e-5, True, False),
    "5anc_adjusted_nopfilter": (False, 1.0, True, False),
}


def _covariate_root(name: str, data_base: Path, out: Path) -> Path:
    import numpy as np
    import pandas as pd

    from appfl_bio_suite.experiments.fine_mapping.real_data.federation_builder import SITE_ORDER
    from appfl_bio_suite.experiments.fine_mapping.real_data.mage import (
        expression_pcs,
        site_local_expression_and_covariates,
    )

    norm, covs = COVARIATE_VARIANTS[name]
    src = data_base / "federation_with_swap_fwl"
    root = out / "roots" / name
    if root.exists():
        shutil.rmtree(root)
    (root / "meta").mkdir(parents=True)
    for f in ("panel_manifest.tsv", "cohort_metadata.tsv"):
        shutil.copy2(src / "meta" / f, root / "meta" / f)
    (root / "meta" / "pooled_plink").symlink_to(src / "meta" / "pooled_plink")
    shutil.copy2(src / "mage_manifest.json", root / "mage_manifest.json")
    npz = dict(np.load(src / "meta" / "expression_panel.npz", allow_pickle=True))
    cohort = pd.read_csv(src / "meta" / "cohort_metadata.tsv", sep="\t", dtype={"IID": str})
    site_expr, site_cov = site_local_expression_and_covariates(
        data_base / "cache",
        npz["ensembl_ids"].tolist(),
        cohort.set_index("IID")["site_id"],
        SITE_ORDER,
        exclude_swapped=False,
        n_factors=covs if isinstance(covs, int) else 15,
    )
    if norm == "site":
        npz["expr_raw"] = site_expr
    published = pd.read_csv(src / "meta" / "covariates.tsv", sep="\t", dtype={"sample_id": str})
    if covs == "published":
        table = published
    elif covs == "pooled60":
        cache = data_base / "cache" / "raw_with_swap"
        order = pd.read_csv(cache / "covariates.tsv", sep="\t", dtype={"sample_id": str})
        if order["sample_id"].tolist() != published["sample_id"].tolist():
            raise RuntimeError("expression cache and covariate table disagree on sample order")
        scores = expression_pcs(np.load(cache / "expression_autosomal.npy").astype(float).T, 60)
        table = published[["sample_id", "sex", "PC1", "PC2", "PC3", "PC4", "PC5"]].copy()
        for j in range(60):
            table[f"ePC{j + 1}"] = scores[:, j]
    else:
        table = site_cov
    table.to_csv(root / "meta" / "covariates.tsv", sep="\t", index=False)
    np.savez(root / "meta" / "expression_panel.npz", **npz)
    return root


def _fixed_prior(data_base: Path, out: Path) -> None:
    import numpy as np
    import pandas as pd
    from MultiSuSiE.susiepy_ss import susie_multi_ss

    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import (
        read_bed_variants,
        read_bim,
        read_fam,
    )
    from appfl_bio_suite.experiments.fine_mapping.multisusie import (
        MomentInputs,
        credible_set_purity,
    )
    from appfl_bio_suite.experiments.fine_mapping.real_data.mage import residualize_design

    meta = data_base / "federation_with_swap_fwl" / "meta"
    pre = meta / "pooled_plink" / "pooled_genotypes"
    bim, fam = read_bim(f"{pre}.bim"), read_fam(f"{pre}.fam")
    panel = pd.read_csv(meta / "panel_manifest.tsv", sep="\t")
    npz = np.load(meta / "expression_panel.npz", allow_pickle=True)
    samples = [str(s) for s in npz["samples"]]
    perm = [samples.index(i) for i in fam["IID"]]
    cov = pd.read_csv(meta / "covariates.tsv", sep="\t", dtype={"sample_id": str})
    Z = cov.set_index("sample_id").loc[fam["IID"].tolist()].to_numpy(float)
    arm = out / "susie_fixedprior_731"
    for loc in panel.itertuples():
        mask = (bim["chrom"].astype(str) == str(loc.chrom)) & bim["bp"].between(
            loc.start_bp, loc.end_bp
        )
        v = bim.loc[mask, ["chrom", "snp_id", "bp", "a1", "a2"]].reset_index(drop=True)
        X = read_bed_variants(pre, bim.loc[mask].index.to_numpy(), len(fam)).astype(float)
        mac = np.minimum(X.sum(0), 2 * len(fam) - X.sum(0))
        keep = mac >= 2 * len(fam) * 0.005 - 1e-10
        k = list(npz["locus_ids"]).index(loc.locus_id)
        Xf, yf = residualize_design(X[:, keep], npz["expr_raw"][k, perm], Z)
        v = v.loc[keep].reset_index(drop=True)
        ok = Xf.var(0) > 1e-12
        Xf, v = Xf[:, ok], v.loc[ok].reset_index(drop=True)
        n = len(yf)
        Xs = (Xf - Xf.mean(0)) / Xf.std(0)
        ys = (yf - yf.mean()) / yf.std()
        xtx, xty = Xs.T @ Xs, Xs.T @ ys
        fit = susie_multi_ss(
            XTX_list=[xtx.copy()],
            XTY_list=[xty.copy()],
            YTY_list=[float(ys @ ys)],
            population_sizes=[n],
            rho=np.ones((1, 1)),
            L=10,
            scaled_prior_variance=float(np.max((xty / n) ** 2)),
            standardize=False,
            estimate_residual_variance=True,
            estimate_prior_variance=False,
            pop_spec_effect_priors=True,
            max_iter=1000,
            tol=1e-5,
            coverage=0.95,
            min_abs_corr=0.5,
            float_type=np.float64,
            low_memory_mode=False,
            R_list=[xtx / n],
            variant_ids=None,
        )
        inputs = MomentInputs(
            v, ["ALL"], [n], [xtx], [xty], [float(ys @ ys)], [np.ones(len(v), bool)]
        )
        dest = arm / "data" / "instances" / f"{loc.locus_id}_mage_eqtl_rep0"
        dest.mkdir(parents=True, exist_ok=True)
        v.assign(pip=fit.pip).to_csv(dest / "variants.tsv", sep="\t", index=False)
        rows = [
            {"component": c + 1, "snp_id": v.snp_id[i]}
            for c, keep_cs in enumerate(fit.sets[3])
            if keep_cs and credible_set_purity(inputs, fit.sets[0][c]) >= 0.5
            for i in fit.sets[0][c]
        ]
        pd.DataFrame(rows, columns=["component", "snp_id"]).to_csv(
            dest / "credible_sets.tsv", sep="\t", index=False
        )
    (arm / "data" / "fed_fm_results.tsv").write_text("diagnostic\n")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("command", choices=["covariates", "susiex", "fixed-prior"])
    p.add_argument("variant", nargs="?")
    p.add_argument("--data-base", type=Path, default=Path("local/data/mage"))
    p.add_argument("--out", type=Path, default=Path("local/output/mage-replication/diagnostics"))
    p.add_argument("--n-workers", type=int, default=6)
    args = p.parse_args()
    from appfl_bio_suite.experiments.fine_mapping.real_data.benchmark import (
        run_centralized_arm,
    )

    data_base, out = args.data_base.resolve(), args.out.resolve()
    if args.command == "fixed-prior":
        _fixed_prior(data_base, out)
        return 0
    if args.command == "covariates":
        root = _covariate_root(args.variant, data_base, out)
        rec = run_centralized_arm(
            root,
            out / args.variant,
            args.variant,
            backend="multisusie",
            fwl=True,
            single_population=True,
            maf=0.005,
            n_signals=10,
            max_iter=1000,
            tol=1e-5,
            multisusie_min_abs_corr=0.5,
            n_workers=args.n_workers,
        )
    else:
        single, pth, adjust, swap = SUSIEX_VARIANTS[args.variant]
        root = data_base / ("federation_with_swap_fwl" if swap else "federation_clean_fwl")
        name = f"susiex_{args.variant}"
        rec = run_centralized_arm(
            root,
            out / name,
            name,
            backend="susiex",
            fwl=True,
            single_population=single,
            pval_thresh=pth,
            maf=0.005,
            n_signals=10,
            max_iter=1000,
            tol=1e-6,
            n_workers=args.n_workers,
            susiex_finite_sample_adjust=adjust,
        )
    print(json.dumps(rec, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
