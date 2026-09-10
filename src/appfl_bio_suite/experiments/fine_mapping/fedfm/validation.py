"""Exhaustive validation harness for the realised FedFM data package.

Runs every validation described in `paper/design.md` §10 against the on-disk
outputs and writes machine-readable results (CSV + a human-readable SUMMARY.txt)
under `reports/validation/`. Designed to run unattended on a Polaris compute
node (see `scripts/submit_polaris_validation.pbs`) so the full ~35k-instance
sweep can be done with no time/compute bound; the chat then just reads the
emitted files.

Checks implemented (all run by default):

  1. Bit-exact phenotype re-derivation (design §10.5 / handoff §7.1) — NEW.
     Independently recompute g_i = Σ_k g_{i,k}·β_{k,pop(i)} from raw .bed +
     β.npy + per-individual ancestry; confirm Var(g)/Var(y_stored) reproduces
     the manifest's recorded empirical_h2_<site>. Also runs a deliberately
     WRONG site-level (dominant-pop) β model to show the per-ancestry lookup
     is actually exercised. Folds in the structural per-instance facts
     (phenotype row count, NaN check) since it already reads every .pheno.

  2. Heritability calibration (design §10.1) — empirical_h2 vs h2_target.

  3. Cross-ancestry effect-size correlation (design §10.2) — per-rg 6×6 sample
     correlation of the per-variant β rows; mean off-diagonal should ≈ rg.

  4. Structural completeness (design §10.3) — file inventories, ncsl ==
     |causal_snp_ids|, causal bp within locus window, phenotype rows == site n,
     no NaN.

  5. Tag-SNP r²/LD matrix sanity (design §10.4) — diagonal ≈ 1, |R| ≤ 1, no NaN
     across every cached .npz.

Run:
    python -m src.validation --config config/simulation_config.yaml \
        --n-workers 32 --out reports/validation
    # quick spot-check on a sample:
    python -m src.validation --sample 200
"""
from __future__ import annotations

import argparse
import functools
import glob
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .phenotype_sim import build_architecture_grid
from .utils import (
    SimulationConfig,
    load_config,
    read_bed_variants,
    read_bim,
    read_fam,
    read_site_manifest,
    setup_logging,
)

# ---------------------------------------------------------------------------
# 1. Bit-exact phenotype re-derivation
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=2)
def _site_io(config_path: str) -> tuple[dict, tuple[str, ...], str]:
    """Per-worker cache of the per-site I/O needed for re-derivation.

    Returns (io_by_site, site_names, ground_truth_dir). lru_cache means each
    loky worker process loads this exactly once and reuses it across rows.
    Only the snp_id→variant-index map is retained from each (large) .bim.
    """
    cfg = load_config(config_path)
    pop_to_idx = {p: i for i, p in enumerate(cfg.superpopulations)}
    io: dict[str, dict] = {}
    for site, site_cfg in cfg.sites.items():
        site_dir = cfg.site_dir(site)
        prefix = site_dir / f"{site}_chr{cfg.chromosome}"
        bim = read_bim(str(prefix) + ".bim")
        fam = read_fam(str(prefix) + ".fam")
        manifest = read_site_manifest(site_dir / f"{site}_manifest.tsv")
        # Align manifest to .fam (= .bed row) order by (FID, IID).
        manifest = (
            manifest.set_index(["FID", "IID"])
            .loc[list(zip(fam["FID"], fam["IID"]))]
            .reset_index()
        )
        pop_index = np.array(
            [pop_to_idx[p] for p in manifest["superpopulation"]], dtype=np.int64
        )
        io[site] = {
            "prefix": str(prefix),
            "fam_iid": fam["IID"].to_numpy(),
            "n": len(fam),
            "pop_index": pop_index,
            "snpid_to_vidx": pd.Series(bim.index.to_numpy(), index=bim["snp_id"].to_numpy()),
            "dom_idx": pop_to_idx[site_cfg.dominant],
        }
    return io, tuple(cfg.sites.keys()), str(cfg.resolved_path("ground_truth_dir"))


def _genetic_value(dosage: np.ndarray, beta: np.ndarray, pop_index: np.ndarray) -> np.ndarray:
    """Independent reimplementation of the per-individual β lookup (NOT importing
    the production helper). Mean-imputes missing dosages column-wise to match
    production so a present-day NaN cannot masquerade as a mismatch.
    """
    dosage = dosage.astype(np.float64, copy=True)
    if np.isnan(dosage).any():
        col_means = np.nanmean(dosage, axis=0)
        r, c = np.where(np.isnan(dosage))
        dosage[r, c] = col_means[c]
    beta_per_ind = beta.T[pop_index]  # (n_samples, ncsl)
    return np.einsum("ik,ik->i", dosage, beta_per_ind)


def _rederive_row(config_path: str, row: dict) -> list[dict]:
    io, sites, gt = _site_io(config_path)
    gt = Path(gt)
    base = f"{row['locus_id']}_{row['architecture_id']}_rep{int(row['replicate'])}"
    beta = np.load(gt / "effect_sizes" / f"{base}.npy").astype(np.float64)
    causal_ids = str(row["causal_snp_ids"]).split(",")

    recs: list[dict] = []
    for site in sites:
        s = io[site]
        common = {
            "locus_id": row["locus_id"],
            "architecture_id": row["architecture_id"],
            "replicate": int(row["replicate"]),
            "site": site,
            "ncsl": int(row["ncsl"]),
            "rg": float(row["rg"]),
            "h2_target": float(row["h2_target"]),
            "manifest_h2": float(row[f"empirical_h2_{site}"]),
        }
        pheno = gt / "phenotypes" / site / f"{base}.pheno"
        if not pheno.exists():
            recs.append({**common, "status": "missing_pheno"})
            continue
        try:
            vidx = s["snpid_to_vidx"].loc[causal_ids].to_numpy(dtype=np.int64)
        except KeyError:
            recs.append({**common, "status": "missing_snp"})
            continue

        dosage = read_bed_variants(s["prefix"], vidx, s["n"])
        ydf = pd.read_csv(pheno, sep="\t", dtype={"FID": str, "IID": str})
        n_rows = len(ydf)
        y_has_nan = bool(np.isnan(ydf["y"].to_numpy()).any())
        y = ydf.set_index("IID")["y"].loc[s["fam_iid"]].to_numpy(dtype=np.float64)
        y_var = float(np.var(y))

        # --- model under test: per-individual ancestry lookup ---
        g = _genetic_value(dosage, beta, s["pop_index"])
        g_var = float(np.var(g))
        rederived_h2 = g_var / y_var if y_var > 0 else 0.0
        mh2 = common["manifest_h2"]
        h2_relerr = abs(rederived_h2 - mh2) / mh2 if mh2 > 0 else abs(rederived_h2 - mh2)

        # --- contrast: WRONG single-site (dominant-pop) β model ---
        wrong_pop = np.full_like(s["pop_index"], s["dom_idx"])
        g_wrong = _genetic_value(dosage, beta, wrong_pop)
        wrong_h2 = float(np.var(g_wrong)) / y_var if y_var > 0 else 0.0
        wrong_relerr = abs(wrong_h2 - mh2) / mh2 if mh2 > 0 else abs(wrong_h2 - mh2)

        recs.append({
            **common,
            "status": "ok",
            "rederived_h2": rederived_h2,
            "h2_relerr": h2_relerr,
            "wrong_model_h2": wrong_h2,
            "wrong_relerr": wrong_relerr,
            "g_var": g_var,
            "y_var": y_var,
            "n_pheno_rows": n_rows,
            "y_has_nan": y_has_nan,
        })
    return recs


def run_rederivation(cfg_path: str, manifest: pd.DataFrame, n_workers: int, logger) -> pd.DataFrame:
    logger.info("rederivation: %d rows x sites, n_workers=%d", len(manifest), n_workers)
    rows = manifest.to_dict("records")
    results = Parallel(n_jobs=n_workers, backend="loky", verbose=5)(
        delayed(_rederive_row)(cfg_path, r) for r in rows
    )
    flat = [rec for sub in results for rec in sub]
    return pd.DataFrame(flat)


def summarise_rederivation(df: pd.DataFrame, rtol: float) -> dict:
    ok = df[df["status"] == "ok"]
    bad_status = df[df["status"] != "ok"]
    n = len(ok)
    n_pass = int((ok["h2_relerr"] <= rtol).sum())
    # Discrimination is INFORMATIONAL, not a pass gate: at rg=1 the per-ancestry
    # and site-level models coincide by construction, so the wrong model is not
    # expected to be rejected there. We just report how often it IS materially off
    # (the rg<1 ancestry-mixed instances), confirming the lookup is exercised.
    off = ok["wrong_relerr"] > 10 * rtol
    materially_off = int(off.sum())
    differ_by_rg = ok.assign(_off=off).groupby("rg")["_off"].mean().round(3).to_dict()
    return {
        "instances_ok": n,
        "instances_bad_status": int(len(bad_status)),
        "h2_match": n_pass,
        "h2_match_frac": n_pass / n if n else 0.0,
        "worst_h2_relerr": float(ok["h2_relerr"].max()) if n else None,
        "median_h2_relerr": float(ok["h2_relerr"].median()) if n else None,
        "wrong_model_materially_off": materially_off,
        "wrong_model_off_frac_by_rg": {str(k): v for k, v in differ_by_rg.items()},
        "rtol": rtol,
        "pass": bool(n > 0 and n_pass == n and len(bad_status) == 0),
    }


# ---------------------------------------------------------------------------
# 2. Heritability calibration (design §10.1)
# ---------------------------------------------------------------------------


def run_h2_calibration(cfg: SimulationConfig, manifest: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sites = list(cfg.sites.keys())
    recs = []
    for site in sites:
        col = f"empirical_h2_{site}"
        emp = manifest[col].to_numpy(dtype=float)
        tgt = manifest["h2_target"].to_numpy(dtype=float)
        relerr = np.abs(emp / tgt - 1.0)
        for i, (lid, aid, rep) in enumerate(zip(manifest["locus_id"], manifest["architecture_id"], manifest["replicate"])):
            recs.append({
                "locus_id": lid, "architecture_id": aid, "replicate": int(rep),
                "site": site, "h2_target": float(tgt[i]),
                "empirical_h2": float(emp[i]), "target_relerr": float(relerr[i]),
            })
    per_inst = pd.DataFrame(recs)
    tol = cfg.phenotype.h2_tolerance_relative
    summary = (
        per_inst.groupby("site")["target_relerr"]
        .agg(median="median", p95=lambda x: np.percentile(x, 95), max="max")
        .reset_index()
    )
    summary["tolerance"] = tol
    summary["pass"] = summary["max"] <= tol
    return per_inst, summary


# ---------------------------------------------------------------------------
# 3. Cross-ancestry effect-size correlation (design §10.2)
# ---------------------------------------------------------------------------


def run_effect_size_correlation(cfg: SimulationConfig, manifest: pd.DataFrame, gt_dir: Path) -> tuple[pd.DataFrame, dict]:
    pops = cfg.superpopulations
    recs = []
    matrices_out: dict[str, list[list[float]]] = {}
    for rg, grp in manifest.groupby("rg"):
        rows_beta = []
        for _, r in grp.iterrows():
            base = f"{r['locus_id']}_{r['architecture_id']}_rep{int(r['replicate'])}"
            beta = np.load(gt_dir / "effect_sizes" / f"{base}.npy")  # (ncsl, n_pops)
            rows_beta.append(np.asarray(beta, dtype=np.float64))
        stacked = np.vstack(rows_beta)  # (sum_ncsl, n_pops)
        # Ancestry-divergent instances carry private-variant β rows that are
        # structurally zero outside one pop; those would bias the cross-pop
        # correlation. Keep only fully-nonzero rows (shared variants and all
        # shared-mode rows), which are the ones that actually carry an rg tie.
        nz = np.all(stacked != 0.0, axis=1)
        stacked = stacked[nz]
        # Need ≥2 rows to estimate a correlation. Heavy divergent-mode use can
        # shrink the fully-nonzero pool below that; record NaN (inconclusive)
        # rather than a spurious value, and let the summary drop NaN groups.
        if stacked.shape[0] < 2:
            recs.append({
                "rg": float(rg), "n_beta_rows": int(stacked.shape[0]),
                "mean_offdiag_corr": float("nan"), "min_offdiag_corr": float("nan"),
                "max_offdiag_corr": float("nan"), "abs_err_vs_rg": float("nan"),
            })
            continue
        corr = np.corrcoef(stacked, rowvar=False)  # (n_pops, n_pops)
        offdiag = corr[~np.eye(len(pops), dtype=bool)]
        recs.append({
            "rg": float(rg),
            "n_beta_rows": int(stacked.shape[0]),
            "mean_offdiag_corr": float(np.mean(offdiag)),
            "min_offdiag_corr": float(np.min(offdiag)),
            "max_offdiag_corr": float(np.max(offdiag)),
            "abs_err_vs_rg": float(abs(np.mean(offdiag) - rg)),
        })
        matrices_out[f"rg={rg}"] = corr.tolist()
    return pd.DataFrame(recs).sort_values("rg").reset_index(drop=True), {"pops": pops, "matrices": matrices_out}


# ---------------------------------------------------------------------------
# 4. Structural completeness (design §10.3)
# ---------------------------------------------------------------------------


def run_structural(cfg: SimulationConfig, manifest: pd.DataFrame, gt_dir: Path,
                   loci: pd.DataFrame, rederiv: pd.DataFrame) -> pd.DataFrame:
    sites = list(cfg.sites.keys())
    n_rows = len(manifest)
    invariants = []

    def add(name, observed, expected, ok):
        invariants.append({"invariant": name, "observed": observed,
                           "expected": expected, "pass": bool(ok)})

    # Expected manifest size = n_loci × n_architectures × n_replicates, derived
    # from the realised config (NOT hardcoded), so the check survives a change in
    # the selected-loci count (e.g. 78 → 79 after a locus-selection re-run).
    n_arch = len(build_architecture_grid(
        ncsl=cfg.architecture.ncsl,
        h2=cfg.architecture.h2,
        rg=cfg.architecture.rg,
        mode=cfg.architecture.factorial_mode,
    ))
    n_rep = cfg.architecture.replicates
    expected_rows = len(loci) * n_arch * n_rep

    add("selected_loci", len(loci), "<=100 after dedup", len(loci) <= 100)
    n_eff = len(list((gt_dir / "effect_sizes").glob("*.npy")))
    add("effect_size_npy_files", n_eff, n_rows, n_eff == n_rows)
    add("manifest_rows", n_rows, f"{len(loci)}*{n_arch}*{n_rep}={expected_rows}",
        n_rows == expected_rows)

    for site in sites:
        n_ph = len(list((gt_dir / "phenotypes" / site).glob("*.pheno")))
        add(f"pheno_files_{site}", n_ph, n_rows, n_ph == n_rows)

    # |causal_snp_ids| matches ncsl (shared mode) or the union size ncsl +
    # (n_pops-1)*n_private (ancestry-divergent mode).
    n_pops = len(cfg.superpopulations)
    n_priv = cfg.architecture.ancestry_divergent_causal.n_private_per_pop
    if "causal_mode" in manifest.columns:
        is_div = (manifest["causal_mode"] == "divergent").to_numpy()
    else:
        is_div = np.zeros(n_rows, dtype=bool)
    csl_len = manifest["causal_snp_ids"].astype(str).str.split(",").apply(len)
    expected_len = manifest["ncsl"].to_numpy() + is_div * (n_pops - 1) * n_priv
    n_match = int((csl_len.to_numpy() == expected_len).sum())
    add("ncsl_matches_|causal_snp_ids|", f"{n_match}/{n_rows}", n_rows, n_match == n_rows)

    # Ancestry-divergent sparsity: each divergent instance's β must have shared
    # rows nonzero in all pops and private rows nonzero in exactly one, and the
    # per-pop causal sets must genuinely differ. (No-op when the mode is off.)
    if is_div.any():
        div_rows = manifest.loc[is_div]
        n_div = len(div_rows)
        n_sparse_ok = 0
        n_sets_differ = 0
        for _, r in div_rows.iterrows():
            base = f"{r['locus_id']}_{r['architecture_id']}_rep{int(r['replicate'])}"
            beta = np.load(gt_dir / "effect_sizes" / f"{base}.npy").astype(np.float64)
            nz_per_row = (beta != 0.0).sum(axis=1)  # pops each variant is causal in
            # exactly n_shared rows are causal in all pops; the rest in exactly one
            if np.all((nz_per_row == n_pops) | (nz_per_row == 1)) \
                    and (nz_per_row == n_pops).sum() == int(r["ncsl"]) - n_priv:
                n_sparse_ok += 1
            try:
                by_pop = json.loads(r["causal_snp_ids_by_pop_json"])
                sets = {frozenset(v) for v in by_pop.values()}
                if len(sets) > 1:
                    n_sets_differ += 1
            except (json.JSONDecodeError, TypeError):
                pass
        add("ancestry_divergent_sparsity", f"{n_sparse_ok}/{n_div}", n_div,
            n_sparse_ok == n_div)
        add("ancestry_divergent_sets_differ", f"{n_sets_differ}/{n_div}", n_div,
            n_sets_differ == n_div)

    # causal bp within locus window
    win = loci.set_index("locus_id")[["start_bp", "end_bp"]]
    m = manifest.copy()
    m["start_bp"] = m["locus_id"].map(win["start_bp"])
    m["end_bp"] = m["locus_id"].map(win["end_bp"])
    def _bp_in(row):
        bps = [int(x) for x in str(row["causal_bp"]).split(",")]
        return all(row["start_bp"] <= b <= row["end_bp"] for b in bps)
    n_in = int(m.apply(_bp_in, axis=1).sum())
    add("causal_bp_within_window", f"{n_in}/{n_rows}", n_rows, n_in == n_rows)

    # From the re-derivation pass (covers every produced .pheno):
    ok = rederiv[rederiv["status"] == "ok"]
    if len(ok):
        for site in sites:
            sub = ok[ok["site"] == site]
            site_n = cfg.sites[site].n
            n_rowsok = int((sub["n_pheno_rows"] == site_n).sum())
            add(f"pheno_rows_{site}==n", f"{n_rowsok}/{len(sub)}", site_n,
                n_rowsok == len(sub))
        n_nan = int(ok["y_has_nan"].sum())
        add("phenotype_y_no_nan", f"{n_nan} NaN files", 0, n_nan == 0)
    return pd.DataFrame(invariants)


# ---------------------------------------------------------------------------
# 5. LD / tag-SNP r² matrix sanity (design §10.4)
# ---------------------------------------------------------------------------


def _ld_one(path: str) -> dict:
    z = np.load(path, allow_pickle=True)
    r = np.asarray(z["r2"], dtype=np.float64)
    n = r.shape[0]
    diag = np.diag(r)
    off = r[~np.eye(n, dtype=bool)]
    return {
        "file": Path(path).name,
        "n_tag": n,
        "median_diag_dev": float(np.median(np.abs(diag - 1.0))),
        "max_diag_dev": float(np.max(np.abs(diag - 1.0))),
        "max_abs_value": float(np.max(np.abs(r))),
        "frac_abs_gt_1": float(np.mean(np.abs(off) > 1.0 + 1e-6)),
        "any_nan": bool(np.isnan(r).any()),
    }


def run_ld_sanity(cfg: SimulationConfig, n_workers: int, logger, ld_tol: float = 1e-2) -> tuple[pd.DataFrame, dict]:
    ld_dir = cfg.resolved_path("loci_dir") / "ld_matrices"
    files = sorted(glob.glob(str(ld_dir / "*.npz")))
    logger.info("ld sanity: %d .npz files", len(files))
    recs = Parallel(n_jobs=n_workers, backend="loky", verbose=2)(
        delayed(_ld_one)(f) for f in files
    )
    df = pd.DataFrame(recs)
    # Values may exceed 1 by a small floor because the r matrix is computed on
    # mean-imputed dosages (~0.4% missing in HAPNEST chr1), so the imputed columns
    # are not exactly unit-variance. We tolerate |R| <= 1 + ld_tol, commensurate
    # with the diagonal deviation, and report the realised overshoot honestly.
    max_abs = float(df["max_abs_value"].max()) if len(df) else 0.0
    summary = {
        "n_files": len(df),
        "median_diag_dev": float(df["median_diag_dev"].median()) if len(df) else None,
        "max_diag_dev": float(df["max_diag_dev"].max()) if len(df) else None,
        "max_abs_value": max_abs,
        "abs_tolerance": 1.0 + ld_tol,
        "files_with_nan": int(df["any_nan"].sum()) if len(df) else 0,
        "files_exceeding_tol": int((df["max_abs_value"] > 1.0 + ld_tol).sum()) if len(df) else 0,
        "files_with_any_overshoot": int((df["max_abs_value"] > 1.0 + 1e-6).sum()) if len(df) else 0,
        "pass": bool(len(df) > 0 and df["any_nan"].sum() == 0 and max_abs <= 1.0 + ld_tol),
    }
    return df, summary


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config/simulation_config.yaml")
    ap.add_argument("--out", default=None,
                    help="where to write the report; defaults to the config's "
                         "reports_dir/validation, so the output lands beside the run "
                         "it validates rather than in the current directory")
    ap.add_argument("--n-workers", type=int, default=8)
    ap.add_argument("--rtol", type=float, default=1e-5,
                    help="relative tolerance on re-derived vs manifest empirical h2")
    ap.add_argument("--ld-tol", type=float, default=1e-2,
                    help="tolerance on |R| overshoot above 1 from missing-data floor")
    ap.add_argument("--sample", type=int, default=0,
                    help="if >0, run re-derivation on a random sample of N manifest rows")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--skip", default="", help="comma list of checks to skip: "
                    "rederivation,h2,corr,structural,ld")
    ap.add_argument("--shard-index", type=int, default=0,
                    help="0-based instance shard for multi-node re-derivation")
    ap.add_argument("--n-shards", type=int, default=1,
                    help="total shards; >1 runs the re-derivation MAP step only")
    ap.add_argument("--merge", action="store_true",
                    help="REDUCE: merge re-derivation shard parts, then run the "
                         "remaining (cheap) checks and write SUMMARY")
    args = ap.parse_args()

    logger = setup_logging(name="validation")
    t0 = time.time()
    cfg = load_config(args.config)
    cfg_path = str(Path(args.config).resolve())
    gt_dir = cfg.resolved_path("ground_truth_dir")
    # Relative to the RUN, not to the caller. This used to default to the literal
    # "reports/validation", which is resolved against the working directory -- so the
    # Polaris job, which cds to the suite checkout before launching stages, wrote every
    # site's validation report into the source tree while its own completion banner
    # pointed at ${FM_RUN_DIR}/reports/validation/SUMMARY.txt, where nothing had been
    # written. An explicit --out still wins.
    out = Path(args.out) if args.out else cfg.resolved_path("reports_dir") / "validation"
    out.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(gt_dir / "causal_manifest.tsv", sep="\t")
    loci = pd.read_csv(cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t")
    skip = {s.strip() for s in args.skip.split(",") if s.strip()}

    # --- multi-node re-derivation MAP step: score this shard, write a part, exit.
    if args.n_shards > 1 and not args.merge:
        if not (0 <= args.shard_index < args.n_shards):
            ap.error("need 0 <= --shard-index < --n-shards")
        shard = manifest.iloc[args.shard_index::args.n_shards].reset_index(drop=True)
        logger.info("re-derivation shard %d/%d: %d instances",
                    args.shard_index, args.n_shards, len(shard))
        part = run_rederivation(cfg_path, shard, args.n_workers, logger)
        tag = f"part{args.shard_index:03d}of{args.n_shards:03d}"
        part.to_csv(out / f"rederivation_per_instance.{tag}.csv", index=False)
        logger.info("wrote shard %s", tag)
        return 0

    # --- REDUCE step: stitch the shard parts back into the full per-instance CSV.
    if args.merge:
        parts = []
        for r in range(args.n_shards):
            tag = f"part{r:03d}of{args.n_shards:03d}"
            p = out / f"rederivation_per_instance.{tag}.csv"
            if not p.exists():
                raise FileNotFoundError(f"Missing re-derivation shard part: {p}")
            parts.append(pd.read_csv(p))
        rederiv = pd.concat(parts, ignore_index=True)
        rederiv.to_csv(out / "rederivation_per_instance.csv", index=False)
        for r in range(args.n_shards):
            (out / f"rederivation_per_instance.part{r:03d}of{args.n_shards:03d}.csv").unlink()
        logger.info("merged %d re-derivation shards -> %d instances",
                    args.n_shards, len(rederiv))
        skip = skip | {"rederivation"}  # already done by the shards; reuse merged df
        return _finish_validation(cfg, cfg_path, gt_dir, out, manifest, loci, skip,
                                  args, logger, t0, rederiv=rederiv)

    rederiv_manifest = manifest
    if args.sample > 0:
        rederiv_manifest = manifest.sample(n=min(args.sample, len(manifest)),
                                           random_state=args.seed).reset_index(drop=True)
    return _finish_validation(cfg, cfg_path, gt_dir, out, manifest, loci, skip,
                              args, logger, t0, rederiv_manifest=rederiv_manifest)


def _finish_validation(cfg, cfg_path, gt_dir, out, manifest, loci, skip, args, logger,
                       t0, rederiv=None, rederiv_manifest=None) -> int:
    """Run checks 2–5, the (optional) re-derivation, and emit the summary.

    ``rederiv`` may be supplied pre-computed (the multi-node REDUCE path passes
    the merged shard results, with "rederivation" in ``skip``); otherwise it is
    computed here for the single-node path."""
    if rederiv is None:
        rederiv = pd.DataFrame()
    if rederiv_manifest is None:
        rederiv_manifest = manifest

    summary: dict = {"config": cfg_path, "n_manifest_rows": len(manifest)}
    if not rederiv.empty:
        summary["rederivation"] = summarise_rederivation(rederiv, args.rtol)

    # --- 1. re-derivation (also feeds structural) ---
    if "rederivation" not in skip:
        logger.info("=== [1/5] phenotype re-derivation ===")
        rederiv = run_rederivation(cfg_path, rederiv_manifest, args.n_workers, logger)
        rederiv.to_csv(out / "rederivation_per_instance.csv", index=False)
        summary["rederivation"] = summarise_rederivation(rederiv, args.rtol)
        logger.info("rederivation summary: %s", summary["rederivation"])

    # --- 2. h2 calibration ---
    if "h2" not in skip:
        logger.info("=== [2/5] heritability calibration ===")
        h2_inst, h2_sum = run_h2_calibration(cfg, manifest)
        h2_inst.to_csv(out / "h2_calibration_per_instance.csv", index=False)
        h2_sum.to_csv(out / "h2_calibration_summary.csv", index=False)
        summary["h2_calibration"] = {
            "per_site": h2_sum.to_dict("records"),
            "pass": bool(h2_sum["pass"].all()),
        }

    # --- 3. effect-size correlation ---
    if "corr" not in skip:
        logger.info("=== [3/5] cross-ancestry effect-size correlation ===")
        corr_df, corr_mats = run_effect_size_correlation(cfg, manifest, gt_dir)
        corr_df.to_csv(out / "effect_size_correlation.csv", index=False)
        (out / "effect_size_correlation_matrices.json").write_text(json.dumps(corr_mats, indent=2))
        summary["effect_size_correlation"] = {
            "per_rg": corr_df.to_dict("records"),
            # Drop inconclusive (NaN) groups — a group with <2 fully-nonzero β
            # rows cannot be scored and must not fail the check.
            "pass": bool((corr_df["abs_err_vs_rg"].dropna() < 0.1).all()),
        }

    # --- 4. structural completeness ---
    if "structural" not in skip and not rederiv.empty:
        logger.info("=== [4/5] structural completeness ===")
        struct = run_structural(cfg, manifest, gt_dir, loci, rederiv)
        struct.to_csv(out / "structural_completeness.csv", index=False)
        summary["structural"] = {
            "invariants": struct.to_dict("records"),
            "pass": bool(struct["pass"].all()),
        }

    # --- 5. LD matrix sanity ---
    if "ld" not in skip:
        logger.info("=== [5/5] LD / tag-SNP r2 matrix sanity ===")
        ld_df, ld_sum = run_ld_sanity(cfg, args.n_workers, logger, args.ld_tol)
        ld_df.to_csv(out / "ld_matrix_per_file.csv", index=False)
        summary["ld_sanity"] = ld_sum

    summary["wall_seconds"] = round(time.time() - t0, 1)
    summary["sampled"] = args.sample if args.sample > 0 else "all"

    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    _write_text_summary(out / "SUMMARY.txt", summary)
    logger.info("validation complete in %.1fs -> %s", summary["wall_seconds"], out)
    print((out / "SUMMARY.txt").read_text())

    checks = [v.get("pass") for k, v in summary.items() if isinstance(v, dict) and "pass" in v]
    return 0 if all(c for c in checks) else 1


def _write_text_summary(path: Path, s: dict) -> None:
    L = []
    L.append("=" * 72)
    L.append("FedFM EXHAUSTIVE VALIDATION SUMMARY")
    L.append("=" * 72)
    L.append(f"config            : {s.get('config')}")
    L.append(f"manifest rows     : {s.get('n_manifest_rows')}")
    L.append(f"re-derivation scope: {s.get('sampled')}")
    L.append(f"wall seconds      : {s.get('wall_seconds')}")
    L.append("")

    def status(d):
        return "PASS" if d.get("pass") else "FAIL"

    if "rederivation" in s:
        r = s["rederivation"]
        L.append(f"[1] PHENOTYPE RE-DERIVATION ............... {status(r)}")
        L.append(f"      instances checked       : {r['instances_ok']}")
        L.append(f"      h2 match within rtol={r['rtol']:g}: {r['h2_match']}/{r['instances_ok']}")
        L.append(f"      worst h2 rel error      : {r['worst_h2_relerr']:.3e}")
        L.append(f"      median h2 rel error     : {r['median_h2_relerr']:.3e}")
        L.append(f"      wrong site-level model materially off (info): "
                 f"{r['wrong_model_materially_off']}/{r['instances_ok']}")
        L.append(f"        fraction off by rg    : {r['wrong_model_off_frac_by_rg']}")
        L.append(f"        (rg=1.0 coincides by construction -> 0 expected)")
        if r["instances_bad_status"]:
            L.append(f"      !! instances with bad status : {r['instances_bad_status']}")
        L.append("")

    if "h2_calibration" in s:
        h = s["h2_calibration"]
        L.append(f"[2] HERITABILITY CALIBRATION ............. {status(h)}")
        for rec in h["per_site"]:
            L.append(f"      {rec['site']:<9} median={rec['median']:.4f} "
                     f"p95={rec['p95']:.4f} max={rec['max']:.4f} (tol={rec['tolerance']})")
        L.append("")

    if "effect_size_correlation" in s:
        c = s["effect_size_correlation"]
        L.append(f"[3] CROSS-ANCESTRY EFFECT-SIZE CORR ...... {status(c)}")
        for rec in c["per_rg"]:
            L.append(f"      rg={rec['rg']:<4} mean_offdiag={rec['mean_offdiag_corr']:.4f} "
                     f"range=[{rec['min_offdiag_corr']:.3f},{rec['max_offdiag_corr']:.3f}] "
                     f"|err|={rec['abs_err_vs_rg']:.4f}")
        L.append("")

    if "structural" in s:
        st = s["structural"]
        L.append(f"[4] STRUCTURAL COMPLETENESS .............. {status(st)}")
        for inv in st["invariants"]:
            mark = "ok " if inv["pass"] else "XX "
            L.append(f"      {mark}{inv['invariant']:<32} obs={inv['observed']} exp={inv['expected']}")
        L.append("")

    if "ld_sanity" in s:
        ld = s["ld_sanity"]
        L.append(f"[5] LD / r2 MATRIX SANITY ................ {status(ld)}")
        L.append(f"      files                   : {ld['n_files']}")
        L.append(f"      median diag deviation   : {ld['median_diag_dev']:.3e}")
        L.append(f"      max diag deviation      : {ld['max_diag_dev']:.3e}")
        L.append(f"      max |value|             : {ld['max_abs_value']:.4f} "
                 f"(tol {ld['abs_tolerance']:.3f})")
        L.append(f"      files w/ NaN            : {ld['files_with_nan']}")
        L.append(f"      files |.|>1 (floor)     : {ld['files_with_any_overshoot']} "
                 f"(exceeding tol: {ld['files_exceeding_tol']})")
        L.append("")

    checks = [(k, v.get("pass")) for k, v in s.items() if isinstance(v, dict) and "pass" in v]
    overall = all(p for _, p in checks)
    L.append("=" * 72)
    L.append(f"OVERALL: {'PASS' if overall else 'FAIL'}  "
             f"({sum(p for _, p in checks)}/{len(checks)} checks passed)")
    L.append("=" * 72)
    path.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
