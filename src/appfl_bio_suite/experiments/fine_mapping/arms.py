"""What federating buys: fine-map the same loci with each site alone, then together.

COORDINATOR-SIDE. Never shipped to a worker.

THE QUESTION THIS ANSWERS, AND WHY IT NEEDED A NEW MODULE
----------------------------------------------------------
Every result in this experiment so far was produced with all three sites participating.
That measures whether federating is *correct* -- and ``fed1_parity`` shows it is, exactly
-- but it says nothing about whether federating is *worth it*, because there is no arm in
which a site works alone. A reader's first question about a federation is "what would I
get on my own?", and nothing in the results table could answer it.

An arm is one fine-mapping run over the same 79 loci, the same architectures, the same
ground truth, differing only in **who participated**:

    anl        50,000 people, 5 ancestry columns, EUR-dominant
    covenant   50,000 people, 3 ancestry columns, 95% AFR
    mbzuai     50,000 people, 4 ancestry columns, MID/CSA
    federation 150,000 people, 6 ancestry columns

That is a federation axis, not an ancestry axis. Ancestry breadth comes along for the
ride -- see ``downsampled``.

THE THREE SITE-ALONE ARMS NEED NO CODE, ONLY A CONFIG
------------------------------------------------------
``plan_columns``, ``materialize_column_keeps`` and ``pooled_phenotype`` are all scoped by
``cfg.sites``, so a config listing one site produces a genuine site-alone fit: that site's
ancestry columns, that site's individuals, that site's LD. Nothing is stubbed and no
vendored code is touched, which is what makes the arms comparable -- every one of them is
the same estimator over the same truth.

``downsampled`` IS THE CONFOUND CONTROL
----------------------------------------
Site-alone and federation differ in two things at once: 50,000 people vs 150,000, and 3-5
ancestry columns vs 6. For the headline claim that is fine and arguably the point -- a
site cannot obtain either on its own. But it leaves "is this just three times the sample
size?" unanswered, so ``downsampled`` runs the full federation restricted to 50,000 people
total, stratified to preserve each site's ancestry mix. Matched on n, any remaining gap is
diversity.

``materialize_column_keeps`` cross-checks the manifest against the config's declared
composition and raises if they disagree, so a down-sampled arm cannot be faked by editing
numbers -- it needs real subsetted manifests. :func:`build_downsampled_site_dir` writes
them into a shadow ``processed/`` tree whose PLINK filesets are symlinks, so the arm costs
a few megabytes rather than a copy of a 287 GB directory.

``ld_borrowed`` IS THE ARGUMENT AGAINST THE SHORTCUT
------------------------------------------------------
The cheap alternative to federating is "just send me your summary statistics and I will
use my own LD panel" -- it is O(M) instead of O(M^2) and needs no site-side compute. This
arm does exactly that: GWAS from one site's cohort, LD from another's. It is the
fine-mapping analogue of applying a model fitted on one cohort to another, and it is the
failure mode the SuSiEx paper studies as LD mismatch. If it produces well-calibrated
credible sets, the O(M^2) uplink is not justified; if it does not, that is the reason the
uplink exists.
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

log = logging.getLogger(__name__)

SITE_ORDER = ("anl", "covenant", "mbzuai")

__all__ = ["Arm", "ARMS", "build_configs", "run_arm", "merge_arms", "main"]


@dataclass(frozen=True)
class Arm:
    """One participation scenario. ``sites`` is who takes part; ``fraction`` subsets them."""

    name: str
    sites: tuple[str, ...]
    fraction: float = 1.0
    ld_from: str | None = None  # borrow LD from this site instead of using own
    label: str = ""

    @property
    def is_solo(self) -> bool:
        return len(self.sites) == 1 and self.fraction == 1.0 and self.ld_from is None


# The published set. `federation` is the arm the existing full-scale run already covers,
# and is listed so the comparison figure can name it even when it is not re-run.
ARMS: tuple[Arm, ...] = (
    Arm("ld_borrowed", ("covenant",), ld_from="anl", label="Covenant statistics, ANL's LD"),
    Arm("federation_50k", SITE_ORDER, fraction=1 / 3, label="All three, n matched to one site"),
    Arm("covenant", ("covenant",), label="Covenant alone"),
    Arm("anl", ("anl",), label="ANL alone"),
    Arm("mbzuai", ("mbzuai",), label="MBZUAI alone"),
    Arm("federation", SITE_ORDER, label="All three"),
)
ARMS_BY_NAME = {a.name: a for a in ARMS}


# --------------------------------------------------------------------------- #
# config construction
# --------------------------------------------------------------------------- #
def _load(config_path: Path) -> dict:
    return yaml.safe_load(Path(config_path).read_text())


def build_downsampled_site_dir(
    base_cfg: dict, out_root: Path, fraction: float, seed: int = 20260601
) -> dict[str, dict[str, int]]:
    """Shadow ``processed/`` tree: symlinked filesets, subsetted manifests.

    Subsets **within ancestry within site**, so a down-sampled federation keeps every
    site's ancestry mix and only shrinks. Sampling by individual across the pooled cohort
    instead would let one site dominate by luck and confound the very thing the arm
    exists to control.

    The PLINK filesets are symlinked rather than copied: ``--keep`` selects individuals at
    read time, so a subset needs no new genotype file, and copying would cost 287 GB.
    """
    out_root = Path(out_root)
    processed = Path(base_cfg["paths"]["processed_dir"])
    rng = np.random.default_rng(seed)
    compositions: dict[str, dict[str, int]] = {}

    for site in base_cfg["sites"]:
        src, dst = processed / site, out_root / site
        dst.mkdir(parents=True, exist_ok=True)
        for f in src.iterdir():
            if f.name == f"{site}_manifest.tsv":
                continue
            link = dst / f.name
            if not link.exists():
                link.symlink_to(f.resolve())

        man = pd.read_csv(src / f"{site}_manifest.tsv", sep="\t")
        keep = []
        for _pop, grp in man.groupby("superpopulation"):
            n = int(round(len(grp) * fraction))
            if n < 1:
                continue
            keep.append(grp.iloc[rng.permutation(len(grp))[:n]])
        sub = pd.concat(keep, ignore_index=True) if keep else man.iloc[:0]
        sub.to_csv(dst / f"{site}_manifest.tsv", sep="\t", index=False)
        compositions[site] = sub["superpopulation"].value_counts().to_dict()
        log.info("  %s: %d of %d individuals kept", site, len(sub), len(man))
    return compositions


def build_configs(
    base_config: Path, out_dir: Path, arms: tuple[Arm, ...] = ARMS, replicates: int | None = None
) -> dict[str, Path]:
    """Write one pipeline config per arm. Returns ``{arm name: config path}``.

    Everything except ``sites:`` and the output paths is copied from the base config
    verbatim, so an arm cannot drift from the published design by accident -- the diff
    between any two arm configs is exactly the participation being tested.
    """
    base = _load(base_config)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    for arm in arms:
        cfg = _load(base_config)
        # Cutting REPLICATES rather than loci is the honest way to make an arm cheaper.
        # The driver's inner loop is `range(cfg.architecture.replicates)`, so lowering it
        # keeps every locus and every architecture -- and therefore every LD-divergence
        # stratum and every design cell -- and only thins the repeats within a cell. The
        # obvious alternative, `--limit`, truncates locus by locus and would silently
        # drop whole strata off the end of the list.
        if replicates is not None:
            cfg["architecture"]["replicates"] = int(replicates)
        missing = [s for s in arm.sites if s not in base["sites"]]
        if missing:
            log.warning("arm %s wants unknown site(s) %s; skipping", arm.name, missing)
            continue

        if arm.fraction != 1.0:
            shadow = out_dir / f"{arm.name}_processed"
            log.info("arm %s: building shadow cohort at %.0f%%", arm.name, arm.fraction * 100)
            comps = build_downsampled_site_dir(base, shadow, arm.fraction)
            cfg["paths"]["processed_dir"] = str(shadow)
            cfg["sites"] = {
                s: {**base["sites"][s], "n": sum(comps[s].values()), "composition": comps[s]}
                for s in arm.sites
                if s in comps
            }
        elif arm.ld_from:
            # Both cohorts must stay described: the borrower supplies the summary
            # statistics and the lender supplies the LD, so scoping the config down to
            # one of them here would throw away the other's composition. The split is
            # done at run time instead, and recorded so the run cannot guess wrong.
            cfg["sites"] = {
                s: base["sites"][s] for s in (*arm.sites, arm.ld_from) if s in base["sites"]
            }
            cfg["arm"] = {"name": arm.name, "borrower": arm.sites[0], "lender": arm.ld_from}
        else:
            cfg["sites"] = {s: base["sites"][s] for s in arm.sites}

        run_dir = out_dir / arm.name
        cfg["paths"]["reports_dir"] = str(run_dir / "reports")
        cfg["paths"]["logs_dir"] = str(run_dir / "logs")
        p = out_dir / f"{arm.name}.yaml"
        p.write_text(yaml.safe_dump(cfg, sort_keys=False))
        written[arm.name] = p

        active = arm.sites if arm.ld_from else tuple(cfg["sites"])
        comp = {s: cfg["sites"][s]["composition"] for s in active}
        pooled: dict[str, int] = {}
        for c in comp.values():
            for pop, n in c.items():
                pooled[pop] = pooled.get(pop, 0) + n
        cols = {k: v for k, v in pooled.items() if v >= base["fine_mapping"]["min_gwas_n"]}
        log.info(
            "arm %-15s n=%7d  %d column(s): %s",
            arm.name,
            sum(pooled.values()),
            len(cols),
            ", ".join(f"{k}={v:,}" for k, v in sorted(cols.items(), key=lambda kv: -kv[1])),
        )
    return written


def _scoped_config(
    config_path: Path, sites: list[str], scratch: Path, reports_dir: Path | None = None
):
    """Load ``config_path`` with ``sites:`` narrowed to ``sites``.

    Written to disk rather than mutated in memory because the config object validates on
    construction -- ``materialize_column_keeps`` cross-checks each column's size against
    the declared composition, and that check is the reason a scoped config is trustworthy.
    """
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import ensure_dir, load_config

    raw = _load(config_path)
    missing = [s for s in sites if s not in raw["sites"]]
    if missing:
        raise ValueError(f"config {config_path} has no site(s) {missing}")
    raw["sites"] = {s: raw["sites"][s] for s in sites}
    scratch = ensure_dir(scratch)
    if reports_dir is not None:
        raw["paths"]["reports_dir"] = str(reports_dir)
    else:
        raw["paths"]["reports_dir"] = str(scratch / "reports")
    raw["paths"]["logs_dir"] = str(scratch / "logs")
    p = scratch / f"{'_'.join(sites)}.yaml"
    p.write_text(yaml.safe_dump(raw, sort_keys=False))
    return load_config(p)


# --------------------------------------------------------------------------- #
# running
# --------------------------------------------------------------------------- #
def run_arm(
    arm_name: str,
    config_path: Path,
    n_workers: int = 32,
    shard_index: int = 0,
    n_shards: int = 1,
    limit: int | None = None,
    level: float = 0.95,
    pval_thresh: float = 1e-5,
) -> Path:
    """Run one arm through the vendored driver, unchanged.

    Delegating rather than reimplementing is the whole point: an arm's numbers and the
    published run's numbers come from the same function, so a difference between arms is
    participation and can be nothing else.
    """
    from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import (
        run_fine_mapping,
    )
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import load_config

    arm = ARMS_BY_NAME.get(arm_name)
    if arm is not None and arm.ld_from:
        return run_ld_borrowed_arm(arm, config_path, n_workers, limit)

    cfg = load_config(Path(config_path))
    run_fine_mapping(
        cfg,
        shard_index=shard_index,
        n_shards=n_shards,
        n_workers=n_workers,
        level=level,
        pval_thresh=pval_thresh,
        keep_work=False,
        limit=limit,
    )
    return Path(cfg.resolved_path("reports_dir")) / "fine_mapping" / "fm_results.tsv"


def run_ld_borrowed_arm(
    arm: Arm,
    config_path: Path,
    n_workers: int = 32,
    limit: int | None = None,
    level: float = 0.95,
    pval_thresh: float = 1e-5,
) -> Path:
    """GWAS from one site's cohort, LD panel from another's.

    Written out here rather than delegated because it is the one arm the vendored driver
    cannot express: ``finemap_instance`` pairs each column's summary statistics with that
    same column's LD, which is exactly the pairing under test. Everything inside the loop
    is still the vendored function -- ``run_gwas``, ``precompute_ld``, ``run_susiex`` and
    ``parse_susiex`` are called unmodified; only which LD file reaches which column
    changes.

    Restricted to ancestries both sites hold, because borrowing a panel for an ancestry
    the lender does not have is not a shortcut anyone would take -- it is just a missing
    column.
    """
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
    from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import (
        ensure_dir,
        load_config,
    )

    cfg = load_config(Path(config_path))
    lender, borrower = arm.ld_from, arm.sites[0]
    repo = cfg.repo_root
    susiex = _resolve_binary("SuSiEx", repo)
    plink = _resolve_binary(cfg.tools.plink, repo)
    plink2 = _resolve_binary(cfg.tools.plink2, repo)

    fm_dir = ensure_dir(Path(cfg.resolved_path("reports_dir")) / "fine_mapping")
    work_root = ensure_dir(fm_dir / "work")

    # The borrower's own columns -- these produce the summary statistics. Scoped, because
    # the arm config deliberately carries the lender's cohort too.
    cfg = _scoped_config(
        config_path,
        [borrower],
        work_root / "borrower_cfg",
        reports_dir=Path(cfg.resolved_path("reports_dir")),
    )
    own = plan_columns(
        {s: dict(sc.composition) for s, sc in cfg.sites.items()},
        cfg.superpopulations,
        cfg.fine_mapping.min_gwas_n,
    )
    materialize_column_keeps(own, cfg, ensure_dir(work_root / "keeps_own"))

    # The lender's columns. Built by scoping the SAME arm config down to the lender,
    # which is why build_configs keeps both cohorts in it -- deriving the lender from the
    # borrower-scoped config is impossible, and guessing it from the base config would
    # let the two drift apart.
    lend_cfg = _scoped_config(config_path, [lender], work_root / "lender_cfg")
    lend_cols = plan_columns(
        {s: dict(sc.composition) for s, sc in lend_cfg.sites.items()},
        lend_cfg.superpopulations,
        lend_cfg.fine_mapping.min_gwas_n,
    )
    materialize_column_keeps(lend_cols, lend_cfg, ensure_dir(work_root / "keeps_lender"))

    shared = [c for c in own if c.pop in {x.pop for x in lend_cols}]
    if not shared:
        raise ValueError(f"{borrower} and {lender} share no ancestry column")
    lend_by_pop = {c.pop: c for c in lend_cols}
    log.info(
        "LD-borrow arm: %s statistics + %s LD, on %d shared column(s): %s",
        borrower,
        lender,
        len(shared),
        ", ".join(c.pop for c in shared),
    )

    enrolled = pooled_ids_path(cfg, work_root / "keeps_own")
    lend_enrolled = pooled_ids_path(lend_cfg, work_root / "keeps_lender")

    loci = pd.read_csv(cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t")
    manifest = pd.read_csv(cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv", sep="\t")
    truth = {
        (r.locus_id, r.architecture_id, int(r.replicate)): str(r.causal_snp_ids).split(",")
        for r in manifest.itertuples(index=False)
    }

    rows: list[dict] = []
    for _, locus in loci.iterrows():
        ref_dir = ensure_dir(work_root / f"{locus['locus_id']}_refs")
        # Two windows: the borrower's individuals for the GWAS, the lender's for the LD.
        win_own = extract_locus_window(cfg, locus, enrolled, ref_dir, plink2)
        lend_dir = ensure_dir(ref_dir / "lender")
        win_lend = extract_locus_window(lend_cfg, locus, lend_enrolled, lend_dir, plink2)
        ld = precompute_ld(
            locus, [lend_by_pop[c.pop] for c in shared], win_lend, lend_dir, plink, 0.005
        )

        # Replicate count comes from the CONFIG, not from the manifest. The manifest
        # always holds all 10; the config is what `--reps` lowers, and the vendored
        # driver's inner loop reads it. Taking it from the manifest here would silently
        # give this arm a different instance set from every other arm, which is exactly
        # the comparison the figure makes.
        jobs = [
            (a, r, truth[(locus["locus_id"], a, r)])
            for a in manifest["architecture_id"].unique()
            for r in range(int(cfg.architecture.replicates))
            if (locus["locus_id"], a, r) in truth
        ]
        if limit is not None:
            jobs = jobs[: max(0, limit - len(rows))]
        if jobs:
            rows += Parallel(n_jobs=n_workers, backend="loky")(
                delayed(finemap_instance)(
                    cfg,
                    locus,
                    aid,
                    rep,
                    tr,
                    shared,
                    win_own,
                    ld,
                    work_root,
                    susiex,
                    plink,
                    plink2,
                    level,
                    pval_thresh,
                    False,
                )
                for aid, rep, tr in jobs
            )
        shutil.rmtree(ref_dir, ignore_errors=True)
        log.info("  %s done (cumulative %d)", locus["locus_id"], len(rows))
        if limit is not None and len(rows) >= limit:
            break

    out = fm_dir / "fm_results.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    log.info("wrote %d row(s) -> %s", len(rows), out)
    return out


# --------------------------------------------------------------------------- #
KEY = ["locus_id", "architecture_id", "replicate"]


def merge_arms(arm_results: dict[str, Path], out_path: Path, pair: bool = True) -> pd.DataFrame:
    """One table, one ``arm`` column, restricted to the instances **every** arm ran.

    THE PAIRING IS NOT OPTIONAL BOOKKEEPING. The arms are compared to each other, and the
    published all-sites run carries ten replicates per architecture while a cheaper arm
    may carry three. Concatenating them as-is would compare each arm on a different
    instance set -- and since instances differ in difficulty, a difference between arms
    could then be a difference in which loci and replicates each happened to cover rather
    than in who participated. Intersecting on (locus, architecture, replicate) makes every
    arm answer the same questions, which is what turns the comparison into a paired one.

    Pass ``pair=False`` only to inspect what an arm produced in isolation.
    """
    frames = []
    for name, p in arm_results.items():
        p = Path(p)
        if not p.exists():
            log.warning("arm %s has no results at %s; skipping", name, p)
            continue
        df = pd.read_csv(p, sep="\t")
        if not set(KEY).issubset(df.columns):
            log.warning("arm %s lacks the key columns; skipping", name)
            continue
        df.insert(0, "arm", name)
        frames.append(df)
    if not frames:
        raise FileNotFoundError("no arm produced results")

    if pair and len(frames) > 1:
        common = None
        for df in frames:
            keys = set(map(tuple, df[KEY].itertuples(index=False, name=None)))
            common = keys if common is None else (common & keys)
        if not common:
            raise ValueError(
                "the arms share no (locus, architecture, replicate) -- they cannot be compared"
            )
        idx = pd.MultiIndex.from_tuples(sorted(common), names=KEY)
        trimmed = []
        for df in frames:
            before = len(df)
            # `arm` is an ordinary column and survives the index round-trip, so it must
            # not be re-inserted afterwards.
            kept = df.set_index(KEY).loc[idx].reset_index()
            trimmed.append(kept)
            log.info(
                "  %-15s %6d of %6d instance(s) kept by pairing",
                df["arm"].iloc[0],
                len(kept),
                before,
            )
        frames = trimmed
    else:
        for df in frames:
            log.info("  %-15s %6d instance(s)", df["arm"].iloc[0], len(df))

    merged = pd.concat(frames, ignore_index=True)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)
    log.info(
        "wrote %d row(s) across %d arm(s) -> %s", len(merged), merged["arm"].nunique(), out_path
    )
    return merged


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Build, run and merge the participation arms.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "--base-config", required=True, help="the published pipeline_config.yaml (all three sites)"
    )
    ap.add_argument("--out-dir", required=True)
    ap.add_argument(
        "--arms",
        default=",".join(a.name for a in ARMS if a.name != "federation"),
        help="comma-separated arm names (default: everything but "
        "`federation`, which the published run already covers)",
    )
    ap.add_argument("--build-only", action="store_true")
    ap.add_argument("--merge-only", action="store_true")
    ap.add_argument(
        "--no-pair",
        action="store_true",
        help="do not restrict arms to their common instances (diagnostic "
        "only -- it makes the arms non-comparable)",
    )
    ap.add_argument("--n-workers", type=int, default=32)
    ap.add_argument(
        "--limit",
        type=int,
        default=None,
        help="truncate each arm to N instances -- SMOKE TESTS ONLY; it "
        "truncates locus by locus and biases the stratum mix",
    )
    ap.add_argument(
        "--reps",
        type=int,
        default=None,
        help="replicates per architecture (default: the config's). Cuts "
        "cost while keeping every locus and architecture.",
    )
    ap.add_argument(
        "--federation-results",
        default=None,
        help="the published full-scale results TSV, folded in as the `federation` arm",
    )
    args = ap.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S"
    )

    out_dir = Path(args.out_dir)
    wanted = [
        ARMS_BY_NAME[n]
        for n in (s.strip() for s in args.arms.split(","))
        if n and n in ARMS_BY_NAME
    ]
    configs = build_configs(Path(args.base_config), out_dir, tuple(wanted), args.reps)
    if args.build_only:
        return 0

    results: dict[str, Path] = {}
    if not args.merge_only:
        for arm in wanted:
            if arm.name not in configs:
                continue
            log.info("=== arm %s: %s ===", arm.name, arm.label)
            try:
                results[arm.name] = run_arm(
                    arm.name, configs[arm.name], n_workers=args.n_workers, limit=args.limit
                )
            except Exception as exc:  # noqa: BLE001 - one arm must not lose the others
                log.error("arm %s failed: %s", arm.name, exc)
    else:
        for arm in wanted:
            p = out_dir / arm.name / "reports" / "fine_mapping" / "fm_results.tsv"
            if p.exists():
                results[arm.name] = p

    if args.federation_results and Path(args.federation_results).exists():
        results["federation"] = Path(args.federation_results)
    if results:
        merge_arms(results, out_dir / "fm_results_by_arm.tsv", pair=not args.no_pair)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
