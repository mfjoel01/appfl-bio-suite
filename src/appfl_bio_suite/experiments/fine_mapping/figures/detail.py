"""Per-credible-set and per-variant records, harvested from SuSiEx's own output.

COORDINATOR-SIDE. WHY THIS EXISTS AS A SEPARATE MODULE
-------------------------------------------------------
``fedfm/fine_mapping.py`` reduces each SuSiEx run to one row of summary metrics and
then deletes the working directory. That is the right shape for a 11,850-instance
sweep -- the alternative is a table nobody can open -- but it discards three things
several figures need, all of which SuSiEx already computed:

  ``cs.summary``  one row per credible set: length, purity, max PIP, and
                  POST-HOC_PROB_POP*, the population-specific causal probability.
                  That last one is the quantity the SuSiEx paper thresholds at 0.8 to
                  say "this signal is causal in this population", and it is the only
                  per-population statement the method makes. It is not recoverable
                  from anything in the results table.
  ``cs.cs``       one row per credible-set MEMBER, carrying per-population BETA, SE
                  and -log10 p as comma-joined fields -- the estimated effect sizes
                  the paper's Figure 4f/4g compares across ancestries.
  ``cs.snp``      PIP for every variant in the window, plus per-population log Bayes
                  factors. The PIP track of a LocusZoom-style panel.

The obvious fix -- widen ``parse_susiex`` -- is not available. ``fedfm/`` is vendored
byte-for-byte from the upstream repository and ``tests/test_fine_mapping_configs.py``
asserts that with ``cmp``; editing it to add a column would break the property that
makes the migration checkable. So this module reads the same files from outside
instead, and a sampled re-run with ``--keep-work`` supplies them.

That re-run is a sample rather than the full sweep on purpose. The per-CS
distributions and the population-specific probabilities are population quantities:
a few hundred instances estimate them to a precision far finer than the figures can
show, and the full sweep is ~56 CPU-hours for records that would not move a line.
The headline metrics in the results table stay the full 11,850.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

__all__ = ["harvest_instance", "harvest_tree", "sample_instances", "run_sample"]


# SuSiEx writes population columns positionally, as POP1..POPn in the order the
# --sst_file list was given. That order is the config's superpopulation list filtered
# by min_gwas_n, which is what plan_columns returns -- so the mapping back to names
# has to be passed in, never guessed. Getting it wrong silently mislabels every
# population-specific probability, which is the kind of error that survives review.
def _pop_names(columns) -> list[str]:
    return [getattr(c, "name", str(c)) for c in columns]


def _split_per_pop(value, n_pop: int) -> list[float]:
    """``"0.058,0.060,0.039"`` -> ``[0.058, 0.060, 0.039]``, NA-tolerant."""
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return [np.nan] * n_pop
    parts = str(value).split(",")
    out = []
    for p in parts:
        try:
            out.append(float(p))
        except ValueError:  # SuSiEx writes a literal "NA" for a missing column
            out.append(np.nan)
    while len(out) < n_pop:
        out.append(np.nan)
    return out[:n_pop]


def harvest_instance(
    work_dir: Path, pops: list[str], truth_snps: list[str], out_name: str = "cs"
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """``(credible_sets, members, variants)`` for one instance's SuSiEx output.

    Every frame carries the instance key columns so the three can be concatenated
    across instances and joined back to the results table without a second lookup.
    Missing or empty files give empty frames, not exceptions: a locus that converged
    to no credible set is a legitimate outcome and shows up as an absent ``.summary``.
    """
    work_dir = Path(work_dir)
    inst = work_dir.name
    truth = set(truth_snps)
    n_pop = len(pops)

    def _read(suffix):
        p = work_dir / f"{out_name}{suffix}"
        if not p.exists() or p.stat().st_size == 0:
            return None
        try:
            df = pd.read_csv(p, sep="\t", comment="#")
        except (pd.errors.EmptyDataError, OSError, ValueError):
            return None
        # SuSiEx writes a bare "NULL" line when nothing converged.
        return None if df.empty or "NULL" in df.columns else df

    # ---- credible sets (one row each) ------------------------------------- #
    summ = _read(".summary")
    members = _read(".cs")
    cs_rows: list[dict] = []
    if summ is not None and "CS_ID" in summ.columns:
        # Which credible sets contain a true causal variant -- the numerator of
        # coverage, and the whole reason this harvest exists. It has to come from the
        # MEMBER table: .summary names only each set's top variant, and the causal
        # variant is frequently in the set without being its maximum.
        contains = {}
        if members is not None and {"CS_ID", "SNP"}.issubset(members.columns):
            for cid, grp in members.groupby("CS_ID"):
                hit = set(grp["SNP"].astype(str)) & truth
                contains[cid] = (len(hit), sorted(hit))
        for row in summ.itertuples(index=False):
            cid = row.CS_ID
            n_hit, hits = contains.get(cid, (0, []))
            rec = {
                "instance": inst,
                "cs_id": int(cid),
                "cs_length": int(getattr(row, "CS_LENGTH", np.nan)),
                "cs_purity": float(getattr(row, "CS_PURITY", np.nan)),
                "max_pip": float(getattr(row, "MAX_PIP", np.nan)),
                "max_pip_snp": str(getattr(row, "MAX_PIP_SNP", "")),
                "bp": getattr(row, "BP", np.nan),
                "n_causal_in_cs": n_hit,
                "contains_causal": bool(n_hit > 0),
                "causal_snps_in_cs": ",".join(hits),
                # Is the set's own top variant the causal one? A set can contain the
                # truth and still point at a tag, and those are different successes.
                "top_is_causal": str(getattr(row, "MAX_PIP_SNP", "")) in truth,
            }
            for i, pop in enumerate(pops, start=1):
                col = f"POST-HOC_PROB_POP{i}"
                v = getattr(row, col.replace("-", "_"), None)
                if v is None and col in summ.columns:
                    v = summ.loc[summ["CS_ID"] == cid, col].iloc[0]
                rec[f"prob_causal_{pop}"] = float(v) if v is not None else np.nan
            for field, key in (
                ("BETA", "beta"),
                ("SE", "se"),
                ("_LOG10P", "log10p"),
                ("REF_FRQ", "frq"),
            ):
                raw = getattr(row, field, None)
                if raw is None:
                    src = {"_LOG10P": "-LOG10P"}.get(field, field)
                    raw = (
                        summ.loc[summ["CS_ID"] == cid, src].iloc[0] if src in summ.columns else None
                    )
                for pop, val in zip(pops, _split_per_pop(raw, n_pop), strict=False):
                    rec[f"{key}_{pop}"] = val
            cs_rows.append(rec)
    cs_df = pd.DataFrame(cs_rows)

    # ---- credible-set members --------------------------------------------- #
    mem_rows: list[dict] = []
    if members is not None and {"CS_ID", "SNP"}.issubset(members.columns):
        for row in members.itertuples(index=False):
            snp = str(row.SNP)
            rec = {
                "instance": inst,
                "cs_id": int(row.CS_ID),
                "snp": snp,
                "bp": getattr(row, "BP", np.nan),
                "cs_pip": float(getattr(row, "CS_PIP", np.nan)),
                "ovrl_pip": float(getattr(row, "OVRL_PIP", np.nan)),
                "is_causal": snp in truth,
            }
            for field, key in (
                ("BETA", "beta"),
                ("SE", "se"),
                ("_LOG10P", "log10p"),
                ("REF_FRQ", "frq"),
            ):
                raw = getattr(row, field, None)
                for pop, val in zip(pops, _split_per_pop(raw, n_pop), strict=False):
                    rec[f"{key}_{pop}"] = val
            mem_rows.append(rec)
    mem_df = pd.DataFrame(mem_rows)

    # ---- every variant's PIP ---------------------------------------------- #
    snp = _read(".snp")
    var_df = pd.DataFrame()
    if snp is not None and "SNP" in snp.columns:
        pip_cols = [c for c in snp.columns if c.startswith("PIP(")]
        if pip_cols:
            var_df = pd.DataFrame(
                {
                    "instance": inst,
                    "snp": snp["SNP"].astype(str),
                    "bp": snp["BP"] if "BP" in snp.columns else np.nan,
                    # Max across credible sets, matching how the results table's top_pip is
                    # defined, so a value here and a value there mean the same thing.
                    "pip": snp[pip_cols].apply(pd.to_numeric, errors="coerce").max(axis=1),
                }
            )
            var_df["is_causal"] = var_df["snp"].isin(truth)
    return cs_df, mem_df, var_df


def harvest_tree(
    work_root: Path, truth_map: dict, pops: list[str], out_dir: Path, logger=None
) -> dict[str, Path]:
    """Harvest every instance directory under ``work_root`` into three TSVs.

    Per-variant records are the bulky one -- ~2,000 rows per instance -- so they are
    written only for instances that produced a credible set. A window with no signal
    contributes 2,000 rows of near-zero PIP and nothing a figure would draw.
    """
    active = logger or log
    work_root, out_dir = Path(work_root), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cs_all, mem_all, var_all = [], [], []
    n_seen = 0
    for d in sorted(work_root.iterdir()):
        if not d.is_dir() or d.name.endswith("_refs") or d.name == "keeps":
            continue
        truth = truth_map.get(d.name)
        if truth is None:
            active.debug("no truth for %s; skipping", d.name)
            continue
        n_seen += 1
        cs, mem, var = harvest_instance(d, pops, truth)
        if len(cs):
            cs_all.append(cs)
        if len(mem):
            mem_all.append(mem)
        if len(var) and len(cs):
            var_all.append(var)
    written = {}
    for name, frames in (
        ("fm_credible_sets.tsv", cs_all),
        ("fm_cs_members.tsv", mem_all),
        ("fm_variant_pips.tsv", var_all),
    ):
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        p = out_dir / name
        df.to_csv(p, sep="\t", index=False)
        written[name] = p
        active.info("wrote %s (%d rows)", p, len(df))
    active.info("harvested %d instance director(ies)", n_seen)
    return written


def sample_instances(
    loci: pd.DataFrame,
    manifest: pd.DataFrame,
    loci_per_stratum: int = 3,
    reps: int = 3,
    seed: int = 20260601,
) -> pd.DataFrame:
    """A stratified sample of instances: every architecture, at loci spanning strata.

    Stratified rather than truncated because ``--limit`` takes the first N instances in
    iteration order, which is all one locus -- and one locus is one LD structure, which
    is the variable these figures are about.
    """
    rng = np.random.default_rng(seed)
    picked = []
    for _stratum, grp in loci.groupby("stratum"):
        take = grp.sample(min(loci_per_stratum, len(grp)), random_state=int(rng.integers(1e6)))
        picked.append(take)
    chosen = pd.concat(picked, ignore_index=True)
    want = manifest[
        manifest["locus_id"].isin(set(chosen["locus_id"])) & (manifest["replicate"] < reps)
    ].copy()
    return want.merge(
        chosen[
            ["locus_id", "stratum", "window_id", "chrom", "start_bp", "end_bp", "divergence_score"]
        ],
        on="locus_id",
        how="left",
    )


def run_sample(
    config_path: Path,
    out_dir: Path,
    loci_per_stratum: int = 3,
    reps: int = 3,
    n_workers: int = 8,
    level: float = 0.95,
    pval_thresh: float = 1e-5,
    keep_raw: bool = False,
    logger=None,
) -> dict[str, Path]:
    """Re-run a stratified sample with the working directories kept, then harvest.

    Every statistic still comes from the vendored code path -- this only chooses which
    instances to run and reads the files afterwards, so a number here and the same
    number in the full sweep are produced by identical code.
    """
    active = logger or log
    from joblib import Parallel, delayed

    from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import (
        _resolve_binary,
        extract_locus_window,
        finemap_instance,
        materialize_column_keeps,
        plan_columns,
        pooled_ids_path,
        precompute_ld,
    )
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import ensure_dir

    try:
        from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import load_config
    except ImportError:  # upstream names this differently across versions
        from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import (
            SimulationConfig as _SC,
        )

        load_config = _SC.load  # type: ignore[assignment]

    cfg = load_config(Path(config_path))
    repo = cfg.repo_root
    susiex = _resolve_binary("SuSiEx", repo)
    plink = _resolve_binary(cfg.tools.plink, repo)
    plink2 = _resolve_binary(cfg.tools.plink2, repo)

    out_dir = ensure_dir(Path(out_dir))
    work_root = ensure_dir(out_dir / "work")
    columns = plan_columns(
        {s: dict(sc.composition) for s, sc in cfg.sites.items()},
        cfg.superpopulations,
        cfg.fine_mapping.min_gwas_n,
    )
    pops = _pop_names(columns)
    keep_dir = ensure_dir(work_root / "keeps")
    materialize_column_keeps(columns, cfg, keep_dir)
    enrolled = pooled_ids_path(cfg, keep_dir)
    active.info("ancestry columns: %s", ", ".join(f"{c.name}(n={c.n})" for c in columns))

    loci = pd.read_csv(cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t")
    manifest = pd.read_csv(cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv", sep="\t")
    want = sample_instances(loci, manifest, loci_per_stratum, reps)
    truth_map = {
        f"{r.locus_id}_{r.architecture_id}_rep{int(r.replicate)}": str(r.causal_snp_ids).split(",")
        for r in want.itertuples(index=False)
    }
    active.info(
        "sampled %d instances over %d loci (%d arch x %d rep)",
        len(want),
        want["locus_id"].nunique(),
        want["architecture_id"].nunique(),
        reps,
    )

    rows = []
    for locus_id, grp in want.groupby("locus_id", sort=True):
        locus = loci.set_index("locus_id").loc[locus_id]
        locus = pd.Series({**locus.to_dict(), "locus_id": locus_id})
        ref_dir = ensure_dir(work_root / f"{locus_id}_refs")
        window = extract_locus_window(cfg, locus, enrolled, ref_dir, plink2)
        ld = precompute_ld(locus, columns, window, ref_dir, plink, 0.005)
        jobs = [
            (r.architecture_id, int(r.replicate), str(r.causal_snp_ids).split(","))
            for r in grp.itertuples(index=False)
        ]
        rows += Parallel(n_jobs=n_workers, backend="loky")(
            delayed(finemap_instance)(
                cfg,
                locus,
                aid,
                rep,
                truth,
                columns,
                window,
                ld,
                work_root,
                susiex,
                plink,
                plink2,
                level,
                pval_thresh,
                True,
            )
            for aid, rep, truth in jobs
        )
        active.info("  %s done (%d instances, cumulative %d)", locus_id, len(jobs), len(rows))

    res = pd.DataFrame(rows)
    res_path = out_dir / "fm_results_sample.tsv"
    res.to_csv(res_path, sep="\t", index=False)
    active.info("wrote %s (%d rows)", res_path, len(res))

    written = harvest_tree(work_root, truth_map, pops, out_dir, active)
    written["fm_results_sample.tsv"] = res_path
    (out_dir / "populations.txt").write_text("\n".join(pops) + "\n")

    if not keep_raw:
        # The harvest is complete and the working tree is ~10 GB of PLINK
        # intermediates. Keep only the ancestry sumstats for the loci an exemplar
        # LocusZoom panel is drawn from -- those cannot be regenerated cheaply.
        keep_inst = sorted(truth_map)[:6]
        for d in work_root.iterdir():
            if d.is_dir() and d.name not in keep_inst and not d.name.endswith("_refs"):
                shutil.rmtree(d, ignore_errors=True)
            elif d.is_dir() and d.name.endswith("_refs"):
                shutil.rmtree(d, ignore_errors=True)
        active.info("pruned working tree, kept %d exemplar instance(s)", len(keep_inst))
    return written


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Re-run a stratified instance sample and harvest SuSiEx detail."
    )
    ap.add_argument("--config", required=True, help="pipeline_config.yaml")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--loci-per-stratum", type=int, default=3)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--n-workers", type=int, default=8)
    ap.add_argument(
        "--keep-raw", action="store_true", help="do not prune the working tree after harvesting"
    )
    ap.add_argument(
        "--harvest-only",
        metavar="WORK_ROOT",
        help="skip the re-run; harvest an existing --keep-work tree",
    )
    args = ap.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )
    if args.harvest_only:
        print(
            "--harvest-only needs the same config to recover truth and populations", file=sys.stderr
        )
        return 2
    run_sample(
        Path(args.config),
        Path(args.out_dir),
        args.loci_per_stratum,
        args.reps,
        args.n_workers,
        keep_raw=args.keep_raw,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
