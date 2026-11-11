"""Real-data FedFM validation suite on MAGE v1.0 + 1000 Genomes NYGC GRCh38."""

from __future__ import annotations

import concurrent.futures
import gzip
import json
import logging
import multiprocessing
import shutil
import time
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.ga4gh.drs import DrsRegistry
from appfl_bio_suite.core.ga4gh.resolve import trim_drs_object
from appfl_bio_suite.experiments.fine_mapping.aggregator import FineMappingAggregator
from appfl_bio_suite.experiments.fine_mapping.backend_outputs import export_susiex
from appfl_bio_suite.experiments.fine_mapping.dataset import get_dataset
from appfl_bio_suite.experiments.fine_mapping.fedfm.fed_fine_mapping import (
    write_ld_panel,
    write_sumstats,
)
from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import (
    _resolve_binary,
    run_susiex,
)
from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import (
    read_bed_variants,
    read_bim,
    read_fam,
)
from appfl_bio_suite.experiments.fine_mapping.inference import (
    parse_outputs,
)
from appfl_bio_suite.experiments.fine_mapping.multisusie import MultiSuSiEOptions
from appfl_bio_suite.experiments.fine_mapping.trainer import (
    SiteFineMappingTrainer,
    harmonize_bim,
)
from appfl_bio_suite.experiments.fine_mapping.trainer import (
    read_bim as trainer_read_bim,
)
from appfl_bio_suite.experiments.fine_mapping.trainer import (
    read_fam as trainer_read_fam,
)

from .mage import POPULATION_ORDER, residualize_design

log = logging.getLogger(__name__)

# DUO research request matching MAGE general research use (DUO:0000042).
MAGE_DATA_USE_REQUEST: dict = {
    "requester": "mage-validation@localhost",
    "purposes": ["DUO:0000038"],
    "non_commercial": True,
    "not_for_profit_organisation": True,
    "publication_agreed": True,
}


def normalize_variant_id(variant_id: str) -> str:
    """Normalize ``chr16:56973441:C:T`` and ``16:56973441:C:T`` to bare ``chrom:pos:ref:alt``."""
    return str(variant_id).removeprefix("chr")


def _active_site_ids(federation_root: Path, sites: Sequence[str] | None = None) -> list[str]:
    cfg_sites = sorted(
        {
            p.parent.parent.name
            for p in federation_root.glob("*/data/selected_loci.tsv")
            if p.parent.parent.name != "meta"
        }
    )
    if sites is None:
        return cfg_sites
    missing = [s for s in sites if s not in set(cfg_sites)]
    if missing:
        raise FileNotFoundError(f"Requested site(s) missing in {federation_root}: {missing}")
    return list(sites)


def _site_pops_for(
    federation_root: Path,
    sites: Sequence[str],
    pops: Sequence[str] | None,
    population_order: Sequence[str] = POPULATION_ORDER,
) -> list[str]:
    present: set[str] = set()
    for site in sites:
        manifest = pd.read_csv(
            federation_root / site / "data" / "site_manifest.tsv",
            sep="\t",
            usecols=["superpopulation"],
        )
        present.update(manifest["superpopulation"].astype(str))
    if pops is not None:
        present &= set(pops)
    ordered = [p for p in population_order if p in present]
    ordered.extend(sorted(present - set(ordered)))
    return ordered


def _arm_complete(out_dir: Path, expected_loci: Sequence[str]) -> bool:
    res_path = out_dir / "data" / "fed_fm_results.tsv"
    if not res_path.is_file():
        return False
    df = pd.read_csv(res_path, sep="\t")
    if len(df) < len(expected_loci):
        return False
    for lid in expected_loci:
        if not list((out_dir / "data" / "instances").glob(f"{lid}_*/variants.tsv")):
            return False
    return True


def _resolve_arm_pops_from_variants(
    variants: pd.DataFrame,
    population_order: Sequence[str] = POPULATION_ORDER,
) -> list[str]:
    present = [c.removeprefix("eligible_") for c in variants.columns if c.startswith("eligible_")]
    ordered = [p for p in population_order if p in present]
    ordered.extend(sorted(set(present) - set(ordered)))
    return ordered


def run_federated_arm(
    federation_root: Path,
    out_dir: Path,
    arm_id: str,
    backend: str = "susiex",
    sites: Sequence[str] | None = None,
    pops: Sequence[str] | None = None,
    population_order: Sequence[str] = POPULATION_ORDER,
    maf: float = 0.005,
    pval_thresh: float = 1e-5,
    level: float = 0.95,
    n_signals: int = 10,
    max_iter: int = 1000,
    tol: float = 1e-6,
    multisusie_rho: float = 0.75,
    multisusie_min_abs_corr: float = 0.5,
    locus_limit: int | None = None,
    n_workers: int = 1,
    skip_existing: bool = True,
    overwrite: bool = False,
) -> dict:
    """Run shipped ``get_dataset -> SiteFineMappingTrainer -> FineMappingAggregator``."""
    federation_root = Path(federation_root).resolve()
    out_dir = Path(out_dir).resolve()
    panel = pd.read_csv(federation_root / "meta" / "panel_manifest.tsv", sep="\t")
    all_loci = panel["locus_id"].tolist()
    expected_loci = all_loci[:locus_limit] if locus_limit else all_loci

    if skip_existing and not overwrite and _arm_complete(out_dir, expected_loci):
        res_df = pd.read_csv(out_dir / "data" / "fed_fm_results.tsv", sep="\t")
        return {
            "arm_id": arm_id,
            "mode": "federated",
            "backend": backend,
            "out_dir": str(out_dir),
            "n_loci": len(res_df),
            "n_converged": int(res_df["converged"].fillna(False).sum()),
            "resumed": True,
        }

    if out_dir.exists() and overwrite:
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    active_sites = _active_site_ids(federation_root, sites=sites)
    requested_pops = list(pops) if pops is not None else list(population_order)
    present_pops = _site_pops_for(federation_root, active_sites, requested_pops, population_order)
    if not present_pops:
        raise ValueError(f"{arm_id}: no populations present across sites={active_sites}")

    registry_path = federation_root / "drs_registry.json"
    registry = DrsRegistry.load(registry_path) if registry_path.is_file() else None
    bundles = registry.site_bundles() if registry is not None else {}

    t0 = time.perf_counter()
    local_models = {}
    ga4gh_objects: dict[str, str] = {}
    for site in active_sites:
        bundle_obj = bundles.get(site)
        drs_payload = trim_drs_object(bundle_obj, registry) if bundle_obj is not None else None
        if bundle_obj is not None:
            ga4gh_objects[site] = bundle_obj.self_uri
        dataset, _ = get_dataset(
            data_dir=federation_root / site / "data",
            site_id=site,
            data_use_request=MAGE_DATA_USE_REQUEST,
            enforce_data_use=True,
            drs_object=drs_payload,
            verify_bundles="metadata" if drs_payload is not None else "off",
        )
        trainer = SiteFineMappingTrainer(
            train_dataset=dataset,
            train_configs={
                "trainer_output_dirname": str(out_dir / "sites" / site),
                "pops": present_pops,
                "locus_limit": locus_limit,
                "uplink_gram_dtype": "float64",
            },
            logger=log,
            client_id=site,
        )
        trainer.train()
        local_models[site] = trainer.get_parameters()

    root = repo_root()
    susiex_bin = _resolve_binary("SuSiEx", root)
    plink_bin = _resolve_binary("plink", root)
    ga4gh_cfg = {}
    if registry is not None:
        ga4gh_cfg = {
            "data_use_request": MAGE_DATA_USE_REQUEST,
            "drs": {"hostname": registry.hostname, "objects": ga4gh_objects},
        }

    aggregator = FineMappingAggregator(
        aggregator_configs={
            "backend": backend,
            "population_order": list(population_order),
            "level": float(level),
            "pval_thresh": float(pval_thresh),
            "maf": float(maf),
            "n_signals": int(n_signals),
            "max_iter": int(max_iter),
            "tol": float(tol),
            "multisusie_rho": float(multisusie_rho),
            "multisusie_min_abs_corr": float(multisusie_min_abs_corr),
            "n_workers": max(1, int(n_workers)),
            "keep_work": False,
            "make_figures": False,
            "output_dir": str(out_dir),
            "vendor_root": str(root),
            "susiex_binary": susiex_bin,
            "plink_binary": plink_bin,
            "causal_manifest": None,
            "ga4gh": ga4gh_cfg,
        },
        logger=log,
    )
    nonconverged_error = ""
    try:
        aggregator.aggregate(local_models)
    except RuntimeError as exc:
        if "did not converge" not in str(exc):
            raise
        nonconverged_error = str(exc)

    res_df = pd.read_csv(out_dir / "data" / "fed_fm_results.tsv", sep="\t")
    return {
        "arm_id": arm_id,
        "mode": "federated",
        "backend": backend,
        "sites": active_sites,
        "pops": present_pops,
        "out_dir": str(out_dir),
        "n_loci": len(res_df),
        "n_converged": int(res_df["converged"].fillna(False).sum()),
        "runtime_s": round(time.perf_counter() - t0, 3),
        "warning": nonconverged_error,
        "resumed": False,
    }


def _ols_ld_and_sumstats(
    Xp: np.ndarray,
    yp: np.ndarray,
    variants: pd.DataFrame,
) -> tuple[np.ndarray, pd.DataFrame, np.ndarray, float]:
    """Independent centralized covariance LD + single-variant OLS on centered (Xp, yp)."""
    n = Xp.shape[0]
    if n < 3:
        raise ValueError(f"Need n >= 3 for centralized OLS, got {n}")
    x_mean = Xp.mean(axis=0)
    Xc = Xp - x_mean
    sxx = np.einsum("ij,ij->j", Xc, Xc)
    yc = yp - float(yp.mean())
    syy = float(yc @ yc)
    sxy = Xc.T @ yc

    cov = (Xc.T @ Xc) / float(n)
    sd = np.sqrt(np.clip(np.diag(cov).copy(), 0.0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        R = cov / np.outer(sd, sd)
    R[~np.isfinite(R)] = 0.0
    np.fill_diagonal(R, 1.0)

    ok = sxx > 0.0
    beta = np.divide(sxy, sxx, out=np.zeros_like(sxy), where=ok)
    sse = np.maximum(syy - beta * sxy, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(sse / float(n - 2) / sxx)
        stat = beta / se
    pval = 2.0 * stats.t.sf(np.abs(stat), n - 2)

    ss = pd.DataFrame(
        {
            "chr": variants["chrom"].astype(str).to_numpy(),
            "snp": variants["snp_id"].to_numpy(),
            "bp": variants["bp"].to_numpy(),
            "A1": variants["a1"].to_numpy(),
            "A2": variants["a2"].to_numpy(),
            "beta": beta,
            "se": se,
            "stat": stat,
            "p": pval,
        }
    ).loc[ok & np.isfinite(se) & (se > 0.0)]
    return R, ss.reset_index(drop=True), x_mean / 2.0, float(np.var(yp, ddof=1))


def _fit_multisusie_central_direct(
    Xp_by_pop: Mapping[str, np.ndarray],
    yp_by_pop: Mapping[str, np.ndarray],
    variants_by_pop: Mapping[str, pd.DataFrame],
    eligible_masks: Mapping[str, np.ndarray],
    pops: Sequence[str],
    destination: Path,
    options: MultiSuSiEOptions,
) -> dict:
    """Fit MultiSuSiE directly from centered pooled genotypes/phenotypes (no moment-pooling)."""
    from MultiSuSiE.susiepy_ss import susie_multi_ss

    destination.mkdir(parents=True, exist_ok=True)
    metadata = ["chrom", "snp_id", "bp", "a1", "a2"]
    frames = [
        variants_by_pop[pop]
        .loc[eligible_masks[pop], metadata]
        .assign(chrom=lambda d: d["chrom"].astype(str))
        for pop in pops
    ]
    union = (
        pd.concat(frames, ignore_index=True)
        .drop_duplicates()
        .sort_values(["chrom", "bp", "snp_id"])
        .reset_index(drop=True)
    )
    index = pd.Index(union["snp_id"])

    xtx_list, xty_list, yty_list, n_list, elig_list = [], [], [], [], []
    # Match multisusie.py population sort order internally so array ordering is canonical.
    ms_pops = sorted(pops)
    for pop in ms_pops:
        mask = eligible_masks[pop]
        Xp = Xp_by_pop[pop][:, mask]
        yp = yp_by_pop[pop]
        Xc = Xp - Xp.mean(axis=0)
        yc = yp - float(yp.mean())
        gram_c = Xc.T @ Xc
        cross_c = Xc.T @ yc
        frame_ids = variants_by_pop[pop].loc[mask, "snp_id"]
        idx = index.get_indexer(frame_ids)
        mat = np.zeros((len(index), len(index)), dtype=np.float64)
        vec = np.zeros(len(index), dtype=np.float64)
        emask = np.zeros(len(index), dtype=bool)
        mat[np.ix_(idx, idx)], vec[idx], emask[idx] = gram_c, cross_c, True
        xtx_list.append(mat)
        xty_list.append(vec)
        yty_list.append(float(yc @ yc))
        n_list.append(int(Xp.shape[0]))
        elig_list.append(emask)

    rho = np.full((len(ms_pops), len(ms_pops)), options.rho, dtype=np.float64)
    np.fill_diagonal(rho, 1.0)
    t0 = time.perf_counter()
    fit = susie_multi_ss(
        XTX_list=xtx_list,
        XTY_list=xty_list,
        YTY_list=yty_list,
        population_sizes=n_list,
        rho=rho,
        L=options.n_signals,
        standardize=False,
        estimate_residual_variance=True,
        estimate_prior_variance=True,
        estimate_prior_method="early_EM",
        pop_spec_effect_priors=True,
        max_iter=options.max_iter,
        tol=options.tol,
        coverage=options.coverage,
        min_abs_corr=options.min_abs_corr,
        float_type=np.float64,
        low_memory_mode=True,
        variant_ids=None,
    )
    runtime_s = time.perf_counter() - t0
    variants = union.copy()
    variants["pip"] = fit.pip
    for pop, emask in zip(ms_pops, elig_list, strict=True):
        variants[f"eligible_{pop}"] = emask
    variants.to_csv(destination / "variants.tsv", sep="\t", index=False)

    members = []
    for comp, keep in enumerate(fit.sets[3]):
        if keep:
            for i in fit.sets[0][comp]:
                members.append(
                    {
                        "component": comp + 1,
                        "snp_id": variants.iloc[i]["snp_id"],
                        "component_pip": float(fit.alpha[comp, i]),
                        "pip": float(fit.pip[i]),
                        "coverage": float(fit.sets[2][comp]),
                    }
                )
    pd.DataFrame(
        members,
        columns=["component", "snp_id", "component_pip", "pip", "coverage"],
    ).to_csv(destination / "credible_sets.tsv", sep="\t", index=False)
    fit_status = (
        ("converged_cs" if members else "converged_no_cs") if fit.converged else "nonconverged"
    )
    rec = {
        "method": "MultiSuSiE",
        "options": asdict(options),
        "populations": ms_pops,
        "sample_sizes": n_list,
        "converged": bool(fit.converged),
        "niter": int(fit.niter),
        "n_credible_sets": int(np.sum(fit.sets[3])),
        "fit_status": fit_status,
        "runtime_s": runtime_s,
    }
    (destination / "fit.json").write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    sizes = [len(fit.sets[0][c]) for c, k in enumerate(fit.sets[3]) if k]
    return {
        "converged": bool(fit.converged),
        "fit_status": fit_status,
        "n_credible_sets": len(sizes),
        "cs_sizes_json": json.dumps(sorted(sizes)),
        "total_cs_snps": int(sum(sizes)),
        "top_pip": float(np.max(fit.pip)) if len(fit.pip) else np.nan,
        "runtime_s": round(runtime_s, 3),
        "error": "",
    }


def _run_centralized_locus(task: dict) -> dict:
    """Execute independent centralized inference for one MAGE locus."""
    lid = task["lid"]
    arch_id = task["arch_id"]
    rep_id = task["rep_id"]
    inst = f"{lid}_{arch_id}_rep{rep_id}"
    pops = task["pops"]
    fwl = task["fwl"]
    backend = task["backend"]
    maf = task["maf"]
    v_win = task["v_win"]

    X_all = read_bed_variants(task["pooled_prefix"], task["v_indices"], task["n_samples"]).astype(
        np.float64
    )
    if fwl:
        X_use, y_use = residualize_design(
            X_all, task["expr_raw"], task["Z_cov"], include_intercept=True
        )
    else:
        X_use = X_all
        y_use = task["expr_resid"]

    Xp_by_pop: dict[str, np.ndarray] = {}
    yp_by_pop: dict[str, np.ndarray] = {}
    vars_by_pop: dict[str, pd.DataFrame] = {}
    elig_by_pop: dict[str, np.ndarray] = {}
    for pop in pops:
        p_rows = task["pop_rows"][pop]
        Xp = X_use[p_rows]
        yp = y_use[p_rows]
        # Use raw allele sum for MAF filter even when FWL residuals are centered
        u_raw = X_all[p_rows].sum(axis=0)
        mac = np.minimum(u_raw, 2.0 * len(p_rows) - u_raw)
        keep = (mac >= 2.0 * len(p_rows) * maf - 1e-10) & (np.var(Xp, axis=0) > 1e-12)
        Xp_by_pop[pop] = Xp
        yp_by_pop[pop] = yp
        vars_by_pop[pop] = v_win
        elig_by_pop[pop] = keep

    inst_dest = task["data_dir"] / "instances" / inst
    row = {
        "locus_id": lid,
        "architecture_id": arch_id,
        "replicate": rep_id,
        "stratum": task["stratum"],
        "backend": backend,
    }
    if backend == "multisusie":
        ms_res = _fit_multisusie_central_direct(
            Xp_by_pop, yp_by_pop, vars_by_pop, elig_by_pop, pops, inst_dest, task["ms_opts"]
        )
        row.update(ms_res)
        return row

    priv = task["work_root"] / inst
    ld_dir = task["work_root"] / f"{lid}_ld"
    ld_dir.mkdir(parents=True, exist_ok=True)
    priv.mkdir(parents=True, exist_ok=True)
    sst_files, ld_prefixes, n_list, col_like = [], [], [], []
    for pop in pops:
        keep = elig_by_pop[pop]
        v_ok = vars_by_pop[pop].loc[keep].reset_index(drop=True)
        Xp_ok = Xp_by_pop[pop][:, keep]
        yp = yp_by_pop[pop]
        R, ss, _freq, _ = _ols_ld_and_sumstats(Xp_ok, yp, v_ok)
        row[f"min_p_{pop}"] = float(ss["p"].min()) if len(ss) else np.nan
        sst_files.append(write_sumstats(ss, priv / f"{pop}.sumstats"))
        # Preserve exact raw allele frequency in .frq even in FWL mode
        raw_frq = X_all[task["pop_rows"][pop]][:, keep].mean(0) / 2.0
        ld_prefixes.append(write_ld_panel(ld_dir / f"{pop}_ld", v_ok, R, raw_frq, len(yp)))
        n_list.append(len(yp))
        col_like.append(type("Col", (), {"pop": pop, "variants": v_ok})())

    t_sx = time.perf_counter()
    run_susiex(
        sst_files,
        ld_prefixes,
        n_list,
        task["chrom"],
        task["start_bp"],
        task["end_bp"],
        priv,
        "cs",
        task["susiex_bin"],
        task["plink_bin"],
        task["level"],
        task["pval_thresh"],
        threads=1,
        options={
            "n_signals": task["n_signals"],
            "max_iter": task["max_iter"],
            "tol": task["tol"],
        },
    )
    metrics = parse_outputs(priv, "cs", [])
    metrics["runtime_s"] = round(time.perf_counter() - t_sx, 3)
    row.update(metrics)
    row["error"] = ""
    # Archive cs.* so POST-HOC_PROB harvesting works symmetrically
    art_dir = task["art_root"] / inst
    art_dir.mkdir(parents=True, exist_ok=True)
    for src in [*priv.glob("cs.*"), *priv.glob("*.sumstats")]:
        if src.is_file():
            with (
                src.open("rb") as fh_in,
                gzip.open(art_dir / f"{src.name}.gz", "wb") as fh_out,
            ):
                shutil.copyfileobj(fh_in, fh_out)
    export_susiex(
        art_dir,
        inst_dest,
        col_like,
        row,
        {
            "n_signals": task["n_signals"],
            "max_iter": task["max_iter"],
            "tol": task["tol"],
            "maf": maf,
            "fwl": fwl,
        },
    )
    shutil.rmtree(priv, ignore_errors=True)
    shutil.rmtree(ld_dir, ignore_errors=True)
    return row


def run_centralized_arm(
    federation_root: Path,
    out_dir: Path,
    arm_id: str,
    backend: str = "susiex",
    fwl: bool = False,
    population_order: Sequence[str] = POPULATION_ORDER,
    maf: float = 0.005,
    pval_thresh: float = 1e-5,
    level: float = 0.95,
    n_signals: int = 10,
    max_iter: int = 1000,
    tol: float = 1e-6,
    multisusie_rho: float = 0.75,
    multisusie_min_abs_corr: float = 0.5,
    locus_limit: int | None = None,
    n_workers: int = 1,
    skip_existing: bool = True,
    overwrite: bool = False,
) -> dict:
    """Independent centralized pooled comparator reading directly from ``meta/pooled_plink``."""
    federation_root = Path(federation_root).resolve()
    out_dir = Path(out_dir).resolve()
    meta_dir = federation_root / "meta"
    panel = pd.read_csv(meta_dir / "panel_manifest.tsv", sep="\t")
    if locus_limit is not None:
        panel = panel.iloc[: int(locus_limit)].copy()
    expected_loci = panel["locus_id"].tolist()

    if skip_existing and not overwrite and _arm_complete(out_dir, expected_loci):
        res_df = pd.read_csv(out_dir / "data" / "fed_fm_results.tsv", sep="\t")
        return {
            "arm_id": arm_id,
            "mode": "centralized_fwl" if fwl else "centralized",
            "backend": backend,
            "out_dir": str(out_dir),
            "n_loci": len(res_df),
            "n_converged": int(res_df["converged"].fillna(False).sum()),
            "resumed": True,
        }

    if out_dir.exists() and overwrite:
        shutil.rmtree(out_dir)
    data_dir = out_dir / "data"
    art_root = out_dir / "artifacts"
    work_root = out_dir / "work"
    for d in (data_dir / "instances", art_root, work_root):
        d.mkdir(parents=True, exist_ok=True)

    pooled_prefix = meta_dir / "pooled_plink" / "pooled_genotypes"
    bim = read_bim(f"{pooled_prefix}.bim")
    fam = read_fam(f"{pooled_prefix}.fam")
    cohort = pd.read_csv(meta_dir / "cohort_metadata.tsv", sep="\t", dtype={"IID": str})
    cohort_idx = cohort.set_index("IID").loc[fam["IID"].tolist()]
    pops = [p for p in population_order if p in set(cohort_idx["continentalGroup"])]
    pop_rows = {
        pop: np.flatnonzero(cohort_idx["continentalGroup"].to_numpy() == pop) for pop in pops
    }

    cov_df = pd.read_csv(meta_dir / "covariates.tsv", sep="\t", dtype={"sample_id": str})
    Z_cov = cov_df.set_index("sample_id").loc[fam["IID"].tolist()].to_numpy(dtype=np.float64)

    npz = np.load(meta_dir / "expression_panel.npz", allow_pickle=True)
    expr_samples = [str(s) for s in npz["samples"].tolist()]
    col_perm = [expr_samples.index(iid) for iid in fam["IID"].tolist()]
    expr_by_locus_raw = {
        str(lid): npz["expr_raw"][idx, col_perm].astype(np.float64)
        for idx, lid in enumerate(npz["locus_ids"].tolist())
    }
    expr_by_locus_resid = {
        str(lid): npz["expr_resid"][idx, col_perm].astype(np.float64)
        for idx, lid in enumerate(npz["locus_ids"].tolist())
    }

    manifest = json.loads((federation_root / "mage_manifest.json").read_text(encoding="utf-8"))
    arch_id = str(manifest["architecture_id"])
    rep_id = int(manifest["replicate_id"])
    root = repo_root()
    susiex_bin = _resolve_binary("SuSiEx", root)
    plink_bin = _resolve_binary("plink", root)

    ms_opts = MultiSuSiEOptions(
        rho=float(multisusie_rho),
        n_signals=int(n_signals),
        max_iter=int(max_iter),
        tol=float(tol),
        coverage=float(level),
        min_abs_corr=float(multisusie_min_abs_corr),
        maf=float(maf),
        min_mac=0.0,
    )

    t0 = time.perf_counter()
    tasks: list[dict] = []
    for locus in panel.itertuples(index=False):
        lid = str(locus.locus_id)
        mask = (bim["chrom"].astype(str) == str(locus.chrom)) & bim["bp"].between(
            int(locus.start_bp), int(locus.end_bp)
        )
        v_win = bim.loc[mask, ["chrom", "snp_id", "bp", "a1", "a2"]].reset_index(drop=True)
        tasks.append(
            {
                "lid": lid,
                "chrom": int(locus.chrom),
                "start_bp": int(locus.start_bp),
                "end_bp": int(locus.end_bp),
                "stratum": str(locus.stratum),
                "arch_id": arch_id,
                "rep_id": rep_id,
                "pooled_prefix": pooled_prefix,
                "v_indices": bim.loc[mask].index.to_numpy(),
                "n_samples": len(fam),
                "v_win": v_win,
                "pop_rows": pop_rows,
                "pops": pops,
                "Z_cov": Z_cov,
                "expr_raw": expr_by_locus_raw[lid],
                "expr_resid": expr_by_locus_resid[lid],
                "backend": backend,
                "fwl": fwl,
                "maf": maf,
                "pval_thresh": pval_thresh,
                "level": level,
                "n_signals": n_signals,
                "max_iter": max_iter,
                "tol": tol,
                "ms_opts": ms_opts,
                "data_dir": data_dir,
                "art_root": art_root,
                "work_root": work_root,
                "susiex_bin": susiex_bin,
                "plink_bin": plink_bin,
            }
        )

    if int(n_workers) > 1 and len(tasks) > 1:
        n_proc = min(int(n_workers), len(tasks))
        ctx = multiprocessing.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_proc, mp_context=ctx) as pool:
            rows = list(pool.map(_run_centralized_locus, tasks))
    else:
        rows = [_run_centralized_locus(task) for task in tasks]

    res_df = pd.DataFrame(rows)
    res_df.to_csv(data_dir / "fed_fm_results.tsv", sep="\t", index=False)
    return {
        "arm_id": arm_id,
        "mode": "centralized_fwl" if fwl else "centralized",
        "backend": backend,
        "out_dir": str(out_dir),
        "n_loci": len(res_df),
        "n_converged": int(res_df["converged"].fillna(False).sum()),
        "runtime_s": round(time.perf_counter() - t0, 3),
        "resumed": False,
    }


def run_naive_site_standardized_arm(
    federation_root: Path,
    out_dir: Path,
    arm_id: str = "naive_site_standardized",
    population_order: Sequence[str] = POPULATION_ORDER,
    maf: float = 0.005,
    pval_thresh: float = 1e-5,
    level: float = 0.95,
    n_signals: int = 10,
    max_iter: int = 1000,
    tol: float = 1e-6,
    locus_limit: int | None = None,
    skip_existing: bool = True,
    overwrite: bool = False,
) -> dict:
    """Deliberate anti-pattern: standardize dosages locally at each site before combining."""
    federation_root = Path(federation_root).resolve()
    out_dir = Path(out_dir).resolve()
    panel = pd.read_csv(federation_root / "meta" / "panel_manifest.tsv", sep="\t")
    if locus_limit is not None:
        panel = panel.iloc[: int(locus_limit)].copy()
    expected_loci = panel["locus_id"].tolist()

    if skip_existing and not overwrite and _arm_complete(out_dir, expected_loci):
        res_df = pd.read_csv(out_dir / "data" / "fed_fm_results.tsv", sep="\t")
        return {
            "arm_id": arm_id,
            "mode": "naive_site_standardized",
            "backend": "susiex",
            "out_dir": str(out_dir),
            "n_loci": len(res_df),
            "n_converged": int(res_df["converged"].fillna(False).sum()),
            "resumed": True,
        }

    if out_dir.exists() and overwrite:
        shutil.rmtree(out_dir)
    data_dir = out_dir / "data"
    art_root = out_dir / "artifacts"
    work_root = out_dir / "work"
    for d in (data_dir / "instances", art_root, work_root):
        d.mkdir(parents=True, exist_ok=True)

    active_sites = _active_site_ids(federation_root)
    site_state: dict[str, dict] = {}
    for site in active_sites:
        sdata = federation_root / site / "data"
        prefix = sdata / "site_genotypes"
        bim = trainer_read_bim(f"{prefix}.bim")
        fam = trainer_read_fam(f"{prefix}.fam")
        ref = pd.read_csv(
            sdata / "reference_variants.tsv",
            sep="\t",
            dtype={"snp_id": str, "chrom": str, "bp": np.int64, "a1": str, "a2": str},
        ).set_index("snp_id")
        bim_h, flip_all = harmonize_bim(bim, ref, site)
        s_man = pd.read_csv(sdata / "site_manifest.tsv", sep="\t", dtype=str)
        site_state[site] = {
            "prefix": prefix,
            "bim": bim_h,
            "fam": fam,
            "flip_all": flip_all,
            "manifest": s_man,
            "pheno_dir": sdata / "phenotypes",
        }

    manifest_meta = json.loads((federation_root / "mage_manifest.json").read_text(encoding="utf-8"))
    arch_id = str(manifest_meta["architecture_id"])
    rep_id = int(manifest_meta["replicate_id"])
    root = repo_root()
    susiex_bin = _resolve_binary("SuSiEx", root)
    plink_bin = _resolve_binary("plink", root)

    # Also read pooled cohort for exact R_pooled vs R_naive distortion diagnostic
    pooled_prefix = federation_root / "meta" / "pooled_plink" / "pooled_genotypes"
    pooled_bim = read_bim(f"{pooled_prefix}.bim")
    pooled_fam = read_fam(f"{pooled_prefix}.fam")
    cohort = (
        pd.read_csv(federation_root / "meta" / "cohort_metadata.tsv", sep="\t", dtype={"IID": str})
        .set_index("IID")
        .loc[pooled_fam["IID"].tolist()]
    )

    t0 = time.perf_counter()
    rows: list[dict] = []
    for locus in panel.itertuples(index=False):
        lid = str(locus.locus_id)
        inst = f"{lid}_{arch_id}_rep{rep_id}"
        priv = work_root / inst
        ld_dir = work_root / f"{lid}_ld"
        priv.mkdir(parents=True, exist_ok=True)
        ld_dir.mkdir(parents=True, exist_ok=True)

        pops = _site_pops_for(federation_root, active_sites, None, population_order)
        # Read pooled truth window to compute pooled MAF filter + LD distortion diagnostic
        p_mask = (pooled_bim["chrom"].astype(str) == str(locus.chrom)) & pooled_bim["bp"].between(
            int(locus.start_bp), int(locus.end_bp)
        )
        p_vars = pooled_bim.loc[p_mask, ["chrom", "snp_id", "bp", "a1", "a2"]].reset_index(
            drop=True
        )
        X_pool_win = read_bed_variants(
            pooled_prefix, pooled_bim.loc[p_mask].index.to_numpy(), len(pooled_fam)
        ).astype(np.float64)

        sst_files, ld_prefixes, n_list, col_like = [], [], [], []
        ld_max_diffs, ld_mean_diffs, z_max_diffs = [], [], []
        row = {
            "locus_id": lid,
            "architecture_id": arch_id,
            "replicate": rep_id,
            "stratum": str(locus.stratum),
            "backend": "susiex",
        }

        for pop in pops:
            # Determine pooled MAF mask identically so variants align
            p_rows = np.flatnonzero(cohort["continentalGroup"].to_numpy() == pop)
            Xp_pool = X_pool_win[p_rows]
            n_tot = len(Xp_pool)
            af_pool = Xp_pool.mean(axis=0) / 2.0
            keep = np.minimum(af_pool, 1.0 - af_pool) >= maf - 1e-10
            v_ok = p_vars.loc[keep].reset_index(drop=True)
            Xc_p = Xp_pool[:, keep] - Xp_pool[:, keep].mean(axis=0)
            sd_p = np.sqrt(np.clip((Xc_p**2).mean(axis=0), 0.0, None))
            with np.errstate(divide="ignore", invalid="ignore"):
                R_exact = (Xc_p.T @ Xc_p / n_tot) / np.outer(sd_p, sd_p)
            R_exact[~np.isfinite(R_exact)] = 0.0
            np.fill_diagonal(R_exact, 1.0)

            R_naive = np.zeros((len(v_ok), len(v_ok)), dtype=np.float64)
            z_naive = np.zeros(len(v_ok), dtype=np.float64)
            w_sum = 0.0
            for site in active_sites:
                st = site_state[site]
                members = st["manifest"].loc[st["manifest"]["superpopulation"] == pop, "IID"]
                if members.empty:
                    continue
                s_bim = st["bim"]
                s_mask = (s_bim["chrom"].astype(str) == str(locus.chrom)) & s_bim["bp"].between(
                    int(locus.start_bp), int(locus.end_bp)
                )
                s_sub = s_bim.loc[s_mask]
                X_s = read_bed_variants(
                    st["prefix"], s_sub.index.to_numpy(), len(st["fam"])
                ).astype(np.float64)
                flip_s = st["flip_all"][s_sub.index.to_numpy()]
                if flip_s.any():
                    X_s[:, flip_s] = 2.0 - X_s[:, flip_s]
                row_of = pd.Series(np.arange(len(st["fam"])), index=st["fam"]["IID"].values)
                s_rows = np.sort(row_of.reindex(members.values).to_numpy(dtype=np.int64))
                # Map site variants in exact pooled v_ok order
                col_idx = pd.Index(s_sub["snp_id"]).get_indexer(v_ok["snp_id"])
                X_sp = X_s[s_rows][:, col_idx]
                ph = pd.read_csv(st["pheno_dir"] / f"{inst}.pheno", sep="\t", dtype={"IID": str})
                iids_s = st["fam"].iloc[s_rows]["IID"].values
                y_sp = ph.set_index("IID")["y"].reindex(iids_s).to_numpy(dtype=np.float64)

                ns = float(len(y_sp))
                Xc_s = X_sp - X_sp.mean(axis=0)
                yc_s = y_sp - y_sp.mean()
                sd_s = np.sqrt(np.clip((Xc_s**2).mean(axis=0), 0.0, None))
                with np.errstate(divide="ignore", invalid="ignore"):
                    R_s = (Xc_s.T @ Xc_s / ns) / np.outer(sd_s, sd_s)
                    sxx_s = np.einsum("ij,ij->j", Xc_s, Xc_s)
                    syy_s = float(yc_s @ yc_s)
                    sxy_s = Xc_s.T @ yc_s
                    b_s = np.where(sxx_s > 0, sxy_s / np.where(sxx_s > 0, sxx_s, 1.0), 0.0)
                    sse_s = np.maximum(syy_s - b_s * sxy_s, 0.0)
                    se_s = np.sqrt(sse_s / max(ns - 2.0, 1.0) / np.where(sxx_s > 0, sxx_s, 1.0))
                    z_s = np.where((sxx_s > 0) & (se_s > 0), b_s / se_s, 0.0)
                R_s[~np.isfinite(R_s)] = 0.0
                np.fill_diagonal(R_s, 1.0)
                R_naive += ns * R_s
                z_naive += np.sqrt(ns) * z_s
                w_sum += ns

            R_naive /= max(w_sum, 1.0)
            np.fill_diagonal(R_naive, 1.0)
            z_comb = z_naive / np.sqrt(max(w_sum, 1.0))
            # Convert combined Stouffer z-score back to pseudo-beta/se on unit variance scale
            b_pseudo = z_comb / np.sqrt(float(n_tot))
            se_pseudo = np.full(len(v_ok), 1.0 / np.sqrt(float(n_tot)), dtype=np.float64)
            p_naive = 2.0 * stats.norm.sf(np.abs(z_comb))
            ss = pd.DataFrame(
                {
                    "chr": v_ok["chrom"].astype(str).to_numpy(),
                    "snp": v_ok["snp_id"].to_numpy(),
                    "bp": v_ok["bp"].to_numpy(),
                    "A1": v_ok["a1"].to_numpy(),
                    "A2": v_ok["a2"].to_numpy(),
                    "beta": b_pseudo,
                    "se": se_pseudo,
                    "stat": z_comb,
                    "p": np.clip(p_naive, 1e-300, 1.0),
                }
            )
            ld_diff = np.abs(R_naive - R_exact)
            ld_max_diffs.append(float(np.max(ld_diff)))
            ld_mean_diffs.append(float(np.mean(ld_diff)))

            # Exact pooled Z-scores for distortion tracking across the same individuals
            ph_map = {}
            for s in active_sites:
                df_ph = pd.read_csv(
                    site_state[s]["pheno_dir"] / f"{inst}.pheno",
                    sep="\t",
                    dtype={"IID": str},
                )
                ph_map.update(dict(zip(df_ph["IID"], df_ph["y"], strict=True)))
            ye = np.array([ph_map[iid] for iid in cohort.index[p_rows]], dtype=np.float64)
            _, ss_exact, _, _ = _ols_ld_and_sumstats(Xp_pool[:, keep], ye, v_ok)
            merged_z = pd.merge(
                ss[["snp", "stat"]],
                ss_exact[["snp", "stat"]],
                on="snp",
                suffixes=("_naive", "_exact"),
            )
            z_max_diffs.append(
                float(np.max(np.abs(merged_z["stat_naive"] - merged_z["stat_exact"])))
            )

            row[f"min_p_{pop}"] = float(ss["p"].min()) if len(ss) else np.nan
            sst_files.append(write_sumstats(ss, priv / f"{pop}.sumstats"))
            ld_prefixes.append(
                write_ld_panel(ld_dir / f"{pop}_ld", v_ok, R_naive, af_pool[keep], n_tot)
            )
            n_list.append(n_tot)
            col_like.append(type("Col", (), {"pop": pop, "variants": v_ok})())

        t_sx = time.perf_counter()
        run_susiex(
            sst_files,
            ld_prefixes,
            n_list,
            int(locus.chrom),
            int(locus.start_bp),
            int(locus.end_bp),
            priv,
            "cs",
            susiex_bin,
            plink_bin,
            level,
            pval_thresh,
            threads=1,
            options={"n_signals": n_signals, "max_iter": max_iter, "tol": tol},
        )
        metrics = parse_outputs(priv, "cs", [])
        metrics["runtime_s"] = round(time.perf_counter() - t_sx, 3)
        row.update(metrics)
        row["max_ld_abs_diff"] = float(max(ld_max_diffs)) if ld_max_diffs else np.nan
        row["mean_ld_abs_diff"] = float(np.mean(ld_mean_diffs)) if ld_mean_diffs else np.nan
        row["max_z_abs_diff"] = float(max(z_max_diffs)) if z_max_diffs else np.nan
        row["error"] = ""
        art_dir = art_root / inst
        art_dir.mkdir(parents=True, exist_ok=True)
        for src in [*priv.glob("cs.*"), *priv.glob("*.sumstats")]:
            if src.is_file():
                with src.open("rb") as fh_in, gzip.open(art_dir / f"{src.name}.gz", "wb") as fh_out:
                    shutil.copyfileobj(fh_in, fh_out)
        export_susiex(
            art_dir,
            data_dir / "instances" / inst,
            col_like,
            row,
            {"n_signals": n_signals, "max_iter": max_iter, "tol": tol, "maf": maf, "naive": True},
        )
        shutil.rmtree(priv, ignore_errors=True)
        shutil.rmtree(ld_dir, ignore_errors=True)
        rows.append(row)

    res_df = pd.DataFrame(rows)
    res_df.to_csv(data_dir / "fed_fm_results.tsv", sep="\t", index=False)
    return {
        "arm_id": arm_id,
        "mode": "naive_site_standardized",
        "backend": "susiex",
        "out_dir": str(out_dir),
        "n_loci": len(res_df),
        "n_converged": int(res_df["converged"].fillna(False).sum()),
        "runtime_s": round(time.perf_counter() - t0, 3),
        "resumed": False,
    }


def _load_instance_cs_sets(inst_dir: Path) -> list[set[str]]:
    cs_path = inst_dir / "credible_sets.tsv"
    if not cs_path.is_file():
        return []
    df = pd.read_csv(cs_path, sep="\t")
    if df.empty or "component" not in df.columns:
        return []
    return [set(g["snp_id"].astype(str)) for _, g in df.groupby("component")]


def compare_arm_pairs(
    ref_dir: Path,
    test_dir: Path,
    comparison_name: str,
) -> tuple[pd.DataFrame, dict]:
    """Compare per-locus PIPs and credible sets between ``ref_dir`` and ``test_dir``."""
    ref_insts = {
        p.parent.name: p.parent
        for p in (Path(ref_dir) / "data" / "instances").glob("*/variants.tsv")
    }
    tst_insts = {
        p.parent.name: p.parent
        for p in (Path(test_dir) / "data" / "instances").glob("*/variants.tsv")
    }
    common = sorted(set(ref_insts) & set(tst_insts))
    rows: list[dict] = []
    for inst in common:
        va = pd.read_csv(ref_insts[inst] / "variants.tsv", sep="\t")
        vb = pd.read_csv(tst_insts[inst] / "variants.tsv", sep="\t")
        merged = pd.merge(
            va[["snp_id", "pip"]].dropna(subset=["pip"]),
            vb[["snp_id", "pip"]].dropna(subset=["pip"]),
            on="snp_id",
            suffixes=("_ref", "_test"),
        )
        if merged.empty:
            continue
        diffs = np.abs(merged["pip_ref"].to_numpy() - merged["pip_test"].to_numpy())
        max_diff = float(np.max(diffs))
        mean_diff = float(np.mean(diffs))
        top_ref = str(merged.loc[merged["pip_ref"].idxmax(), "snp_id"])
        top_tst = str(merged.loc[merged["pip_test"].idxmax(), "snp_id"])
        pearson = (
            float(np.corrcoef(merged["pip_ref"], merged["pip_test"])[0, 1])
            if merged["pip_ref"].std() > 0 and merged["pip_test"].std() > 0
            else 1.0
        )

        cs_a = _load_instance_cs_sets(ref_insts[inst])
        cs_b = _load_instance_cs_sets(tst_insts[inst])
        exact_cs = sorted(tuple(sorted(s)) for s in cs_a) == sorted(tuple(sorted(s)) for s in cs_b)
        union_a = set().union(*cs_a) if cs_a else set()
        union_b = set().union(*cs_b) if cs_b else set()
        jaccard = len(union_a & union_b) / len(union_a | union_b) if (union_a or union_b) else 1.0
        rows.append(
            {
                "comparison": comparison_name,
                "instance": inst,
                "locus_id": inst.split("_")[0],
                "n_variants_common": len(merged),
                "n_variants_ref": int(va["pip"].notna().sum()),
                "n_variants_test": int(vb["pip"].notna().sum()),
                "max_pip_abs_diff": max_diff,
                "mean_pip_abs_diff": mean_diff,
                "pip_pearson_r": pearson,
                "top_snp_ref": top_ref,
                "top_snp_test": top_tst,
                "top_snp_identical": top_ref == top_tst,
                "n_cs_ref": len(cs_a),
                "n_cs_test": len(cs_b),
                "exact_cs_identical": bool(exact_cs),
                "cs_union_jaccard": float(jaccard),
                "max_pip_ref": float(merged["pip_ref"].max()),
                "max_pip_test": float(merged["pip_test"].max()),
            }
        )
    df = pd.DataFrame(rows)
    summary = {
        "comparison": comparison_name,
        "n_loci": len(df),
        "worst_max_pip_abs_diff": float(df["max_pip_abs_diff"].max()) if len(df) else np.nan,
        "median_max_pip_abs_diff": float(df["max_pip_abs_diff"].median()) if len(df) else np.nan,
        "mean_pip_pearson_r": float(df["pip_pearson_r"].mean()) if len(df) else np.nan,
        "frac_exact_cs_identical": float(df["exact_cs_identical"].mean()) if len(df) else np.nan,
        "frac_top_snp_identical": float(df["top_snp_identical"].mean()) if len(df) else np.nan,
        "mean_cs_union_jaccard": float(df["cs_union_jaccard"].mean()) if len(df) else np.nan,
    }
    return df, summary


def evaluate_arm_against_mage(
    arm_dir: Path,
    federation_root: Path,
    arm_id: str,
    backend: str = "susiex",
) -> pd.DataFrame:
    """Evaluate one arm's resolution and concordance with published MAGE fine-mapping findings."""
    arm_dir = Path(arm_dir)
    meta_dir = Path(federation_root) / "meta"
    panel = pd.read_csv(meta_dir / "panel_manifest.tsv", sep="\t")
    pub_cs = pd.read_csv(meta_dir / "published_cs_summary.tsv", sep="\t")
    pub_vars = pd.read_csv(meta_dir / "published_cs_variants.tsv", sep="\t")
    pub_vars["norm_id"] = pub_vars["variant_kgpID"].map(normalize_variant_id)

    manifest_meta = json.loads(
        (Path(federation_root) / "mage_manifest.json").read_text(encoding="utf-8")
    )
    arch_id = str(manifest_meta["architecture_id"])
    rep_id = int(manifest_meta["replicate_id"])
    res_df = pd.read_csv(arm_dir / "data" / "fed_fm_results.tsv", sep="\t").set_index("locus_id")

    rows: list[dict] = []
    for locus in panel.itertuples(index=False):
        lid = str(locus.locus_id)
        eid = str(locus.ensembl_id)
        inst = f"{lid}_{arch_id}_rep{rep_id}"
        inst_dir = arm_dir / "data" / "instances" / inst
        va_path = inst_dir / "variants.tsv"
        va = pd.read_csv(va_path, sep="\t") if va_path.is_file() else pd.DataFrame()
        sets = _load_instance_cs_sets(inst_dir)
        sizes = sorted(len(s) for s in sets)

        pips = (
            dict(
                zip(
                    va["snp_id"].astype(str),
                    pd.to_numeric(va["pip"], errors="coerce"),
                    strict=True,
                )
            )
            if not va.empty
            else {}
        )
        valid_pips = {k: v for k, v in pips.items() if np.isfinite(v)}
        top_snp = max(valid_pips, key=valid_pips.get) if valid_pips else ""
        top_pip = float(valid_pips[top_snp]) if top_snp else np.nan

        pv = pub_vars[pub_vars["ensemblID"] == eid].copy()
        pc = pub_cs[pub_cs["ensemblID"] == eid].copy()
        pub_sets = (
            [set(g["norm_id"]) for _, g in pv.groupby("variantCredibleSet")] if not pv.empty else []
        )
        pub_union = set().union(*pub_sets) if pub_sets else set()
        pub_lead_norm = normalize_variant_id(locus.published_lead_kgp)

        arm_union = set().union(*sets) if sets else set()
        pub_sets_captured = sum(bool(ps & arm_union) for ps in pub_sets)
        arm_sets_matching_pub = sum(bool(s & pub_union) for s in sets)
        best_cs_matching = min((len(s) for s in sets if s & pub_union), default=np.nan)
        pub_lead_in_cs = pub_lead_norm in arm_union
        pub_lead_pip_arm = float(valid_pips.get(pub_lead_norm, np.nan)) if valid_pips else np.nan

        best_jaccard = 0.0
        for s in sets:
            for ps in pub_sets:
                if s | ps:
                    best_jaccard = max(best_jaccard, len(s & ps) / len(s | ps))

        # PIP correlation across published CS variants tested in arm
        sub_corr = pd.DataFrame(
            {
                "norm_id": pv["norm_id"],
                "mage_pip": pd.to_numeric(pv["variantPIP"], errors="coerce"),
            }
        )
        sub_corr["arm_pip"] = sub_corr["norm_id"].map(valid_pips)
        sub_corr = sub_corr.dropna(subset=["mage_pip", "arm_pip"])
        if len(sub_corr) >= 3 and sub_corr["mage_pip"].std() > 0 and sub_corr["arm_pip"].std() > 0:
            r_pip = float(np.corrcoef(sub_corr["mage_pip"], sub_corr["arm_pip"])[0, 1])
            rho_pip = float(stats.spearmanr(sub_corr["mage_pip"], sub_corr["arm_pip"]).statistic)
        else:
            r_pip, rho_pip = np.nan, np.nan

        res_row = res_df.loc[lid] if lid in res_df.index else pd.Series(dtype=object)
        rows.append(
            {
                "arm_id": arm_id,
                "backend": backend,
                "locus_id": lid,
                "ensembl_id": eid,
                "gene_symbol": str(locus.gene_symbol),
                "chrom": int(locus.chrom),
                "tss": int(locus.tss),
                "stratum": str(locus.stratum),
                "max_wc_fst": float(locus.max_wc_fst),
                "max_fst_pop": str(locus.max_fst_pop),
                "page_max_pph4": float(locus.page_max_pph4),
                "page_trait": str(locus.page_trait) if pd.notna(locus.page_trait) else "",
                "published_n_cs": int(len(pc)),
                "published_min_cs_size": int(locus.min_cs_size),
                "published_median_cs_size": float(locus.median_cs_size),
                "published_lead_snp": pub_lead_norm,
                "published_lead_pip": float(locus.published_lead_pip),
                "converged": bool(res_row.get("converged", False)),
                "fit_status": str(res_row.get("fit_status", "")),
                "n_cs": len(sets),
                "min_cs_size": int(sizes[0]) if sizes else np.nan,
                "median_cs_size": float(np.median(sizes)) if sizes else np.nan,
                "best_matching_cs_size": best_cs_matching,
                "total_cs_snps": int(sum(sizes)),
                "n_singleton_cs": int(sum(1 for sz in sizes if sz == 1)),
                "n_cs_le5": int(sum(1 for sz in sizes if sz <= 5)),
                "top_snp": top_snp,
                "top_pip": top_pip,
                "n_pip50": int(sum(v > 0.5 for v in valid_pips.values())),
                "n_pip95": int(sum(v > 0.95 for v in valid_pips.values())),
                "pub_cs_captured": int(pub_sets_captured),
                "any_pub_cs_captured": bool(pub_sets_captured > 0),
                "frac_pub_cs_captured": (
                    float(pub_sets_captured / len(pub_sets)) if pub_sets else np.nan
                ),
                "arm_cs_matching_pub": int(arm_sets_matching_pub),
                "pub_lead_in_cs": bool(pub_lead_in_cs),
                "pub_lead_arm_pip": pub_lead_pip_arm,
                "best_cs_jaccard": float(best_jaccard),
                "mage_cs_pip_pearson": r_pip,
                "mage_cs_pip_spearman": rho_pip,
                "cs_size_delta_vs_mage_min": (
                    float(sizes[0] - int(locus.min_cs_size)) if sizes else np.nan
                ),
                "runtime_s": float(res_row.get("runtime_s", np.nan)),
            }
        )
    return pd.DataFrame(rows)


def summarize_arms(arm_loci_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-locus MAGE validation metrics by arm."""
    rows: list[dict] = []
    for arm_id, grp in arm_loci_df.groupby("arm_id", sort=False):
        with_cs = grp[grp["n_cs"] > 0]
        rows.append(
            {
                "arm_id": arm_id,
                "backend": grp["backend"].iloc[0],
                "n_loci": len(grp),
                "n_converged": int(grp["converged"].sum()),
                "loci_with_cs": int(len(with_cs)),
                "total_cs": int(grp["n_cs"].sum()),
                "total_singleton_cs": int(grp["n_singleton_cs"].sum()),
                "total_cs_le5": int(grp["n_cs_le5"].sum()),
                "median_min_cs_size": (
                    float(with_cs["min_cs_size"].median()) if len(with_cs) else np.nan
                ),
                "mean_min_cs_size": (
                    float(with_cs["min_cs_size"].mean()) if len(with_cs) else np.nan
                ),
                "frac_loci_pub_cs_captured": float(grp["any_pub_cs_captured"].mean()),
                "total_pub_cs_captured": int(grp["pub_cs_captured"].sum()),
                "total_pub_cs": int(grp["published_n_cs"].sum()),
                "frac_pub_lead_in_cs": float(grp["pub_lead_in_cs"].mean()),
                "mean_pub_lead_arm_pip": float(grp["pub_lead_arm_pip"].mean()),
                "median_top_pip": float(grp["top_pip"].median()),
                "mean_best_cs_jaccard": float(grp["best_cs_jaccard"].mean()),
                "total_pip95_variants": int(grp["n_pip95"].sum()),
            }
        )
    return pd.DataFrame(rows)


def harvest_posthoc_probs(
    arm_dir: Path,
    federation_root: Path,
    population_order: Sequence[str] = POPULATION_ORDER,
) -> pd.DataFrame:
    """Harvest SuSiEx ``POST-HOC_PROB_POP*`` and CS overlap with published MAGE annotations."""
    arm_dir = Path(arm_dir)
    meta_dir = Path(federation_root) / "meta"
    panel = pd.read_csv(meta_dir / "panel_manifest.tsv", sep="\t").set_index("locus_id")
    pub_vars = pd.read_csv(meta_dir / "published_cs_variants.tsv", sep="\t")
    pub_vars["norm_id"] = pub_vars["variant_kgpID"].map(normalize_variant_id)

    rows: list[dict] = []
    for summ_gz in sorted((arm_dir / "artifacts").glob("*/cs.summary.gz")):
        inst = summ_gz.parent.name
        lid = inst.split("_")[0]
        if lid not in panel.index:
            continue
        locus = panel.loc[lid]
        va = pd.read_csv(arm_dir / "data" / "instances" / inst / "variants.tsv", sep="\t")
        pops = _resolve_arm_pops_from_variants(va, population_order)
        cs_df = pd.read_csv(arm_dir / "data" / "instances" / inst / "credible_sets.tsv", sep="\t")
        summ = pd.read_csv(summ_gz, sep="\t", comment="#")
        if summ.empty or "CS_ID" not in summ.columns:
            continue
        pv = set(pub_vars.loc[pub_vars["ensemblID"] == locus["ensembl_id"], "norm_id"])
        for _, r in summ.iterrows():
            cid = int(r["CS_ID"])
            members = set(cs_df.loc[cs_df["component"] == cid, "snp_id"].astype(str))
            rec = {
                "locus_id": lid,
                "gene_symbol": str(locus["gene_symbol"]),
                "stratum": str(locus["stratum"]),
                "max_wc_fst": float(locus["max_wc_fst"]),
                "max_fst_pop": str(locus["max_fst_pop"]),
                "cs_id": cid,
                "cs_length": int(r["CS_LENGTH"]),
                "cs_purity": float(r["CS_PURITY"]),
                "max_pip": float(r["MAX_PIP"]),
                "max_pip_snp": str(r["MAX_PIP_SNP"]),
                "overlaps_published_mage_cs": bool(members & pv),
                "n_active_pops_prob80": 0,
            }
            n_80 = 0
            for idx, pop in enumerate(pops, start=1):
                col = f"POST-HOC_PROB_POP{idx}"
                val = pd.to_numeric(r.get(col, np.nan), errors="coerce")
                rec[f"prob_{pop}"] = float(val) if pd.notna(val) else np.nan
                if pd.notna(val) and float(val) >= 0.80:
                    n_80 += 1
            rec["n_active_pops_prob80"] = n_80
            rows.append(rec)
    return pd.DataFrame(rows)


def extract_locus_regional_profiles(
    federation_root: Path,
    arm_dirs: Mapping[str, Path],
    locus_id: str,
) -> pd.DataFrame:
    """Build unified per-variant table (position, multi-arm PIPs, min -log10P, MAGE PIP)."""
    meta_dir = Path(federation_root) / "meta"
    panel = pd.read_csv(meta_dir / "panel_manifest.tsv", sep="\t").set_index("locus_id")
    locus = panel.loc[locus_id]
    manifest_meta = json.loads(
        (Path(federation_root) / "mage_manifest.json").read_text(encoding="utf-8")
    )
    inst = f"{locus_id}_{manifest_meta['architecture_id']}_rep{manifest_meta['replicate_id']}"

    pub_vars = pd.read_csv(meta_dir / "published_cs_variants.tsv", sep="\t")
    pv = pub_vars[pub_vars["ensemblID"] == locus["ensembl_id"]].copy()
    pv["snp_id"] = pv["variant_kgpID"].map(normalize_variant_id)

    base: pd.DataFrame | None = None
    for arm_name, adir in arm_dirs.items():
        v_path = Path(adir) / "data" / "instances" / inst / "variants.tsv"
        cs_path = Path(adir) / "data" / "instances" / inst / "credible_sets.tsv"
        if not v_path.is_file():
            continue
        vdf = pd.read_csv(v_path, sep="\t")
        cdf = (
            pd.read_csv(cs_path, sep="\t")
            if cs_path.is_file()
            else pd.DataFrame(columns=["snp_id", "component"])
        )
        cs_map = (
            dict(zip(cdf["snp_id"].astype(str), cdf["component"], strict=True))
            if not cdf.empty
            else {}
        )
        sub = vdf[["chrom", "snp_id", "bp", "a1", "a2"]].copy()
        sub[f"pip_{arm_name}"] = pd.to_numeric(vdf["pip"], errors="coerce")
        sub[f"cs_{arm_name}"] = sub["snp_id"].map(cs_map)

        # Attach min(-log10 P) across ancestries if SuSiEx sumstats exist
        art_dir = Path(adir) / "artifacts" / inst
        minp_series: pd.Series | None = None
        for gz in art_dir.glob("*.sumstats.gz"):
            sdf = pd.read_csv(gz, sep="\t", usecols=["snp", "p"])
            vals_lp = -np.log10(np.clip(sdf["p"].to_numpy(dtype=float), 1e-300, 1.0))
            lp = pd.Series(vals_lp, index=sdf["snp"].astype(str))
            minp_series = (
                lp
                if minp_series is None
                else np.maximum(minp_series.reindex(lp.index, fill_value=0.0), lp)
            )
        if minp_series is not None:
            sub[f"max_neglog10p_{arm_name}"] = sub["snp_id"].map(minp_series).astype(float)

        if base is None:
            base = sub
        else:
            merge_cols = [c for c in sub.columns if c not in ("chrom", "bp", "a1", "a2")]
            base = base.merge(sub[merge_cols], on="snp_id", how="outer")

    if base is None:
        raise FileNotFoundError(f"No variant outputs found for {locus_id}")
    mage_pip = dict(
        zip(pv["snp_id"], pd.to_numeric(pv["variantPIP"], errors="coerce"), strict=True)
    )
    mage_cs = dict(zip(pv["snp_id"], pv["variantCredibleSet"].astype(str), strict=True))
    base["mage_pip"] = base["snp_id"].map(mage_pip)
    base["mage_cs"] = base["snp_id"].map(mage_cs)
    return base.sort_values("bp").reset_index(drop=True)


def _execute_arm_task(spec: dict) -> dict:
    """Dispatch a single validation arm specification inside a worker process."""
    kind = spec["kind"]
    kwargs = spec["kwargs"]
    if kind == "federated":
        return run_federated_arm(**kwargs)
    if kind == "centralized":
        return run_centralized_arm(**kwargs)
    if kind == "naive":
        return run_naive_site_standardized_arm(**kwargs)
    raise ValueError(f"Unknown arm kind: {kind}")


def run_validation_suite(
    clean_root: Path,
    out_root: Path,
    swap_root: Path | None = None,
    maf: float = 0.005,
    n_signals: int = 10,
    max_iter: int = 1000,
    tol: float = 1e-6,
    multisusie_max_iter: int = 1000,
    multisusie_tol: float = 1e-5,
    multisusie_min_abs_corr: float = 0.5,
    max_workers: int = 4,
    run_single_ancestries: bool = True,
    run_solo_sites: bool = True,
    locus_limit: int | None = None,
    skip_existing: bool = True,
) -> dict:
    """Execute the full MAGE real-data FedFM validation suite and write summary tables."""
    clean_root = Path(clean_root).resolve()
    out_root = Path(out_root).resolve()
    arms_dir = out_root / "arms"
    reports_dir = out_root / "reports"
    arms_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Order longest arms first so worker pools schedule MultiSuSiE without tail stragglers.
    arm_specs: list[dict] = [
        {
            "kind": "federated",
            "kwargs": {
                "federation_root": clean_root,
                "out_dir": arms_dir / "federated_multisusie",
                "arm_id": "federated_multisusie",
                "backend": "multisusie",
                "maf": maf,
                "n_signals": n_signals,
                "max_iter": multisusie_max_iter,
                "tol": multisusie_tol,
                "multisusie_min_abs_corr": multisusie_min_abs_corr,
                "locus_limit": locus_limit,
                "n_workers": max(1, int(max_workers)),
                "skip_existing": skip_existing,
            },
        },
        {
            "kind": "centralized",
            "kwargs": {
                "federation_root": clean_root,
                "out_dir": arms_dir / "centralized_multisusie",
                "arm_id": "centralized_multisusie",
                "backend": "multisusie",
                "fwl": False,
                "maf": maf,
                "n_signals": n_signals,
                "max_iter": multisusie_max_iter,
                "tol": multisusie_tol,
                "multisusie_min_abs_corr": multisusie_min_abs_corr,
                "locus_limit": locus_limit,
                "n_workers": max(1, int(max_workers)),
                "skip_existing": skip_existing,
            },
        },
        {
            "kind": "federated",
            "kwargs": {
                "federation_root": clean_root,
                "out_dir": arms_dir / "federated_susiex",
                "arm_id": "federated_susiex",
                "backend": "susiex",
                "maf": maf,
                "n_signals": n_signals,
                "max_iter": max_iter,
                "tol": tol,
                "locus_limit": locus_limit,
                "skip_existing": skip_existing,
            },
        },
        {
            "kind": "centralized",
            "kwargs": {
                "federation_root": clean_root,
                "out_dir": arms_dir / "centralized_susiex",
                "arm_id": "centralized_susiex",
                "backend": "susiex",
                "fwl": False,
                "maf": maf,
                "n_signals": n_signals,
                "max_iter": max_iter,
                "tol": tol,
                "locus_limit": locus_limit,
                "skip_existing": skip_existing,
            },
        },
        {
            "kind": "centralized",
            "kwargs": {
                "federation_root": clean_root,
                "out_dir": arms_dir / "centralized_susiex_fwl",
                "arm_id": "centralized_susiex_fwl",
                "backend": "susiex",
                "fwl": True,
                "maf": maf,
                "n_signals": n_signals,
                "max_iter": max_iter,
                "tol": tol,
                "locus_limit": locus_limit,
                "skip_existing": skip_existing,
            },
        },
        {
            "kind": "naive",
            "kwargs": {
                "federation_root": clean_root,
                "out_dir": arms_dir / "naive_site_standardized",
                "arm_id": "naive_site_standardized",
                "maf": maf,
                "n_signals": n_signals,
                "max_iter": max_iter,
                "tol": tol,
                "locus_limit": locus_limit,
                "skip_existing": skip_existing,
            },
        },
    ]

    if swap_root is not None:
        arm_specs.append(
            {
                "kind": "federated",
                "kwargs": {
                    "federation_root": Path(swap_root).resolve(),
                    "out_dir": arms_dir / "with_swap_federated_susiex",
                    "arm_id": "with_swap_federated_susiex",
                    "backend": "susiex",
                    "maf": maf,
                    "n_signals": n_signals,
                    "max_iter": max_iter,
                    "tol": tol,
                    "locus_limit": locus_limit,
                    "skip_existing": skip_existing,
                },
            }
        )

    if run_solo_sites:
        for site in _active_site_ids(clean_root):
            arm_specs.append(
                {
                    "kind": "federated",
                    "kwargs": {
                        "federation_root": clean_root,
                        "out_dir": arms_dir / f"{site}_susiex",
                        "arm_id": site,
                        "backend": "susiex",
                        "sites": [site],
                        "maf": maf,
                        "n_signals": n_signals,
                        "max_iter": max_iter,
                        "tol": tol,
                        "locus_limit": locus_limit,
                        "skip_existing": skip_existing,
                    },
                }
            )

    if run_single_ancestries:
        for pop in POPULATION_ORDER:
            arm_specs.append(
                {
                    "kind": "federated",
                    "kwargs": {
                        "federation_root": clean_root,
                        "out_dir": arms_dir / f"anc_{pop}_susiex",
                        "arm_id": f"anc_{pop}",
                        "backend": "susiex",
                        "pops": [pop],
                        "maf": maf,
                        "n_signals": n_signals,
                        "max_iter": max_iter,
                        "tol": tol,
                        "locus_limit": locus_limit,
                        "skip_existing": skip_existing,
                    },
                }
            )

    # Run fast single-process arms in parallel, then fan MultiSuSiE loci across workers cleanly.
    results_by_idx: dict[int, dict] = {}
    ms_indexed = [
        (i, s) for i, s in enumerate(arm_specs) if s["kwargs"].get("backend") == "multisusie"
    ]
    other_indexed = [
        (i, s) for i, s in enumerate(arm_specs) if s["kwargs"].get("backend") != "multisusie"
    ]
    if int(max_workers) > 1 and len(other_indexed) > 1:
        n_proc = min(int(max_workers), len(other_indexed))
        ctx = multiprocessing.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_proc, mp_context=ctx) as pool:
            other_res = list(pool.map(_execute_arm_task, [s for _, s in other_indexed]))
        for (i, _), res in zip(other_indexed, other_res, strict=True):
            results_by_idx[i] = res
    else:
        for i, spec in other_indexed:
            results_by_idx[i] = _execute_arm_task(spec)
    for i, spec in ms_indexed:
        results_by_idx[i] = _execute_arm_task(spec)
    run_records = [results_by_idx[i] for i in range(len(arm_specs))]

    # Pairwise exactness & control comparisons
    pair_specs = [
        (
            "susiex_fed_vs_central",
            arms_dir / "centralized_susiex",
            arms_dir / "federated_susiex",
        ),
        (
            "multisusie_fed_vs_central",
            arms_dir / "centralized_multisusie",
            arms_dir / "federated_multisusie",
        ),
        (
            "susiex_fwl_vs_yresid",
            arms_dir / "centralized_susiex_fwl",
            arms_dir / "federated_susiex",
        ),
        (
            "susiex_naive_std_vs_central",
            arms_dir / "centralized_susiex",
            arms_dir / "naive_site_standardized",
        ),
        (
            "fed_multisusie_vs_fed_susiex",
            arms_dir / "federated_susiex",
            arms_dir / "federated_multisusie",
        ),
    ]
    if swap_root is not None:
        pair_specs.append(
            (
                "sample_swap_with_vs_clean",
                arms_dir / "federated_susiex",
                arms_dir / "with_swap_federated_susiex",
            )
        )

    pair_frames, pair_summaries = [], []
    for name, ref_p, tst_p in pair_specs:
        if ref_p.exists() and tst_p.exists():
            df_p, summ_p = compare_arm_pairs(ref_p, tst_p, name)
            pair_frames.append(df_p)
            pair_summaries.append(summ_p)

    paired_loci_df = pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()
    parity_summary_df = pd.DataFrame(pair_summaries)
    paired_loci_df.to_csv(reports_dir / "paired_arm_comparisons.tsv", sep="\t", index=False)
    parity_summary_df.to_csv(reports_dir / "parity_summary.tsv", sep="\t", index=False)

    # MAGE evaluation per arm
    eval_frames = []
    for rec in run_records:
        eval_frames.append(
            evaluate_arm_against_mage(
                Path(rec["out_dir"]), clean_root, arm_id=rec["arm_id"], backend=rec["backend"]
            )
        )
    arm_loci_df = pd.concat(eval_frames, ignore_index=True)
    arm_summary_df = summarize_arms(arm_loci_df)
    arm_loci_df.to_csv(reports_dir / "arm_locus_metrics.tsv", sep="\t", index=False)
    arm_summary_df.to_csv(reports_dir / "arm_summary.tsv", sep="\t", index=False)

    # Ancestry posthoc probabilities for primary federated SuSiEx
    posthoc_df = harvest_posthoc_probs(arms_dir / "federated_susiex", clean_root)
    posthoc_df.to_csv(reports_dir / "ancestry_posthoc_summary.tsv", sep="\t", index=False)

    # Naive distortion metrics per locus
    naive_res_path = arms_dir / "naive_site_standardized" / "data" / "fed_fm_results.tsv"
    if naive_res_path.is_file():
        ndf = pd.read_csv(naive_res_path, sep="\t")
        want_cols = [
            "locus_id",
            "stratum",
            "max_ld_abs_diff",
            "mean_ld_abs_diff",
            "max_z_abs_diff",
        ]
        cols = [c for c in want_cols if c in ndf]
        if cols:
            ndf[cols].to_csv(reports_dir / "naive_ld_distortion.tsv", sep="\t", index=False)

    report = {
        "clean_root": str(clean_root),
        "out_root": str(out_root),
        "runs": run_records,
        "parity_summary": pair_summaries,
        "arm_summary": arm_summary_df.to_dict(orient="records"),
    }
    (reports_dir / "validation_report.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8"
    )
    return report
