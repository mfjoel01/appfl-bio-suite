"""Downstream SuSiEx cross-ancestry fine-mapping — the centralized baseline.

SuSiEx takes a *list of population columns*, one ``(sumstats, LD, n)`` triple
per column, and jointly fine-maps them in a single run. Here **one column = one
superpopulation, pooled across every site**: the AFR column is all 60,000 AFR
individuals (7,500 ANL + 47,500 Covenant + 5,000 MBZUAI) treated as a single
cohort. Six columns cover all 150,000 individuals exactly once.

This is the centralized reference point — it is precisely the fit you would get
if all three sites' data lived in one place, since pooling by ancestry with no
regard to site boundaries is what "centralized" means here. It is also the
estimand the simulator is built against: effect sizes are indexed by
superpopulation, not by site (``design.md`` §5.2), so an ancestry column is
exactly one homogeneous effect — which is the population unit SuSiEx's model
assumes. Maximal per-ancestry sample size and in-sample LD computed on the full
pooled cohort make this the upper bound against which federated variants (which
must reconstruct these columns from site-resident aggregates) are compared.

For each (locus, architecture, replicate) instance this stage:

  1. Extracts the locus window ONCE per locus from the HAPNEST source fileset,
     restricted to the 150,000 individuals enrolled at the three sites, then
     derives each ancestry's LD reference panel from it via ``plink --keep`` and
     computes the SuSiEx LD matrix. Both depend only on (ancestry, locus), so
     they are built once per locus and reused across every architecture x
     replicate instance there.
  2. Concatenates the three sites' phenotype files for the instance into one
     pooled phenotype, then runs a per-ancestry GWAS (``plink2 --glm``,
     ``--keep``-restricted to that ancestry) on the window and reformats the
     result into SuSiEx summary-statistics format.
  3. Runs SuSiEx across all ancestry columns' sumstats + precomputed LD to
     produce credible sets (``.cs`` / ``.snp`` / ``.summary``).
  4. Evaluates the credible sets against the ground-truth causal variants in
     ``causal_manifest.tsv`` and emits one standard-fine-mapping-metrics row.

Note that the per-site phenotype noise calibration (``design.md`` §5.3) sets
sigma^2_eps per *site*, so a pooled ancestry column concatenates phenotypes
whose noise was calibrated against different site cohorts. The marginal effect
estimates stay consistent; only the homoscedasticity assumption is mildly
violated, which costs a little efficiency and nothing in correctness.

SuSiEx is a C++ CLI (installed by ``scripts/install_susiex.sh`` into
``vendor/bin/SuSiEx``); we shell out to it exactly as we do for PLINK. It
consumes GWAS summary statistics + LD reference panels, NOT phenotypes — hence
the GWAS step in front.

Sharding mirrors ``phenotype_sim``: ``--n-shards N`` partitions loci by
``index % N`` so each node runs a disjoint subset; ``--merge`` reduces the
per-shard result parts into ``reports/fine_mapping/fm_results.tsv``.
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .phenotype_sim import build_architecture_grid
from .utils import (
    SimulationConfig,
    ensure_dir,
    get_logger,
    load_config,
    read_site_manifest,
    run_plink,
    setup_logging,
    write_ids_file,
)

# SuSiEx sumstats column order we always emit (1-indexed positions passed to
# SuSiEx below): chr=1 snp=2 bp=3 A1=4 A2=5 beta=6 se=7 stat=8 p=9.
_SUMSTATS_COLS = ["chr", "snp", "bp", "A1", "A2", "beta", "se", "stat", "p"]

# SuSiEx probes exactly these three files at the ``--ld_file`` prefix; if all
# are present it sets precmp=true, skips LD computation, and no longer needs
# ``--ref_file`` (main.cpp:204-222).
_LD_SUFFIXES = (".ld.bin", "_frq.frq", "_ref.bim")

# Metric columns every result row carries, failed instances included, so the
# result table has one schema regardless of how many instances errored.
# ``parse_susiex`` returns exactly these (asserted in tests).
_METRIC_COLS = (
    "n_causal", "n_credible_sets", "cs_sizes_json", "total_cs_snps",
    "n_causal_captured", "any_causal_captured", "best_cs_size",
    "causal_pip_max", "causal_pip_mean", "top_pip",
    "mean_cs_purity", "min_cs_purity", "converged",
)


# ---------------------------------------------------------------------------
# Tool resolution
# ---------------------------------------------------------------------------

def _resolve_binary(name: str, repo_root: Path) -> str:
    """Absolute path to a tool: PATH first, then repo ``vendor/bin``."""
    found = shutil.which(name)
    if found:
        return found
    vendored = repo_root / "vendor" / "bin" / name
    if vendored.exists():
        return str(vendored)
    raise FileNotFoundError(
        f"{name} not found on PATH or in vendor/bin. "
        f"For SuSiEx run scripts/install_susiex.sh; for plink run scripts/install_plink.sh."
    )


# ---------------------------------------------------------------------------
# SuSiEx population columns
# ---------------------------------------------------------------------------

@dataclass
class FMColumn:
    """One SuSiEx population column: a ``(sumstats, LD, n)`` triple to build.

    ``pop`` is the superpopulation this column covers, pooled across every site
    that holds it; ``n`` is that pooled cohort size; ``sites`` records which
    sites contributed (diagnostics only — the column is one undivided cohort).
    ``keep_path`` is the ``plink --keep`` ID file selecting the pooled
    individuals out of the locus window.
    """

    pop: str
    n: int
    sites: tuple[str, ...] = ()
    keep_path: Path | None = None

    @property
    def name(self) -> str:
        return self.pop


def plan_columns(
    site_compositions: dict[str, dict[str, int]],
    superpopulations: list[str], min_n: int,
) -> list[FMColumn]:
    """Plan one pooled column per superpopulation, in ``superpopulations`` order.

    Each ancestry's individuals are pooled across every site holding them, so a
    column's ``n`` is the sum over sites (AFR: 7500 + 47500 + 5000 = 60000).
    Ancestries whose pooled cohort is below ``min_n`` are dropped. Pure (no I/O):
    ``keep_path`` is filled later by ``materialize_column_keeps``.
    """
    cols: list[FMColumn] = []
    for pop in superpopulations:
        contributors = tuple(
            site for site, comp in site_compositions.items() if comp.get(pop)
        )
        n = sum(comp.get(pop, 0) for comp in site_compositions.values())
        if n >= min_n and contributors:
            cols.append(FMColumn(pop, int(n), contributors))
    if not cols:
        raise ValueError(
            f"No fine-mapping columns after planning (min_gwas_n={min_n}); "
            f"loosen min_gwas_n"
        )
    return cols


def materialize_column_keeps(
    columns: list[FMColumn], cfg: SimulationConfig, keep_dir: Path
) -> None:
    """Write each ancestry column's pooled ``--keep`` ID file; set ``keep_path``.

    Reads every site manifest once and concatenates each ancestry's individuals
    across sites into a single ID file — this is where the pooling happens.
    """
    ensure_dir(keep_dir)
    manifests = pd.concat(
        [read_site_manifest(cfg.site_dir(s) / f"{s}_manifest.tsv") for s in cfg.sites],
        ignore_index=True,
    )
    for col in columns:
        sub = manifests[manifests["superpopulation"] == col.pop]
        if len(sub) != col.n:
            raise RuntimeError(
                f"Pooled {col.pop} cohort is {len(sub)} individuals across the "
                f"site manifests but config composition says {col.n}"
            )
        path = keep_dir / f"{col.name}.keep"
        write_ids_file(sub, path)
        col.keep_path = path


def pooled_ids_path(cfg: SimulationConfig, keep_dir: Path) -> Path:
    """Write (once) the union of all enrolled individuals across the three sites.

    The locus window is cut from the HAPNEST source fileset, which holds all
    1,008,000 synthetic individuals; this restricts it to the 150,000 actually
    enrolled at a site. Sites are disjoint ``--keep`` subsets of HAPNEST sharing
    its variant list, so cutting the window centrally is equivalent to cutting it
    per site and merging, without the allele-order hazards of ``--bmerge``.
    """
    ensure_dir(keep_dir)
    path = keep_dir / "enrolled.keep"
    if not path.exists():
        man = pd.concat(
            [read_site_manifest(cfg.site_dir(s) / f"{s}_manifest.tsv") for s in cfg.sites],
            ignore_index=True,
        )
        write_ids_file(man, path)
    return path


# ---------------------------------------------------------------------------
# Per-locus reference panels + per-column GWAS
# ---------------------------------------------------------------------------

def extract_locus_window(
    cfg: SimulationConfig, locus: pd.Series, enrolled_keep: Path,
    out_dir: Path, plink2: str,
) -> Path:
    """Cut the locus [start,end] window as ONE pooled panel over all sites.

    Taken from the HAPNEST source fileset restricted to the enrolled 150,000
    (``enrolled_keep``) rather than merging three per-site windows: the sites are
    disjoint ``--keep`` subsets of HAPNEST sharing its variant list and order, so
    the result is identical while avoiding ``--bmerge``.

    Returns the bfile prefix. Idempotent: skips the cut if the .bed exists.
    """
    chrom, start, end = int(locus["chrom"]), int(locus["start_bp"]), int(locus["end_bp"])
    ensure_dir(out_dir)
    prefix = out_dir / "pooled_ref"
    if not prefix.with_suffix(".bed").exists():
        bf = cfg.resolved_path("hapnest_dir") / f"chr{cfg.chromosome}"
        run_plink(
            ["--bfile", str(bf), "--keep", str(enrolled_keep), "--chr", str(chrom),
             "--from-bp", str(start), "--to-bp", str(end),
             "--make-bed", "--out", str(prefix)],
            binary=plink2,
        )
    return prefix


def _ld_ready(prefix: Path) -> bool:
    """True if SuSiEx would find a usable precomputed LD matrix at ``prefix``."""
    return all(Path(f"{prefix}{suffix}").exists() for suffix in _LD_SUFFIXES)


def precompute_ld(
    locus: pd.Series, columns: list[FMColumn], window: Path,
    out_dir: Path, plink: str, maf: float,
) -> dict[str, Path]:
    """Compute each ancestry's SuSiEx LD matrix ONCE per (ancestry, locus).

    LD is a function of the window genotypes alone — it does not depend on the
    phenotype — so it is invariant across the ``n_arch x n_rep`` instances at a
    locus. Computing it here instead of letting each instance's SuSiEx call
    rebuild it collapses that many redundant PLINK passes into one.

    Each panel is the pooled window restricted to that ancestry's individuals
    across all sites via ``--keep``, so both the MAF filter and the correlations
    are computed *within* the ancestry over its full pooled cohort — the same LD
    reference a centralized analysis holding all the data would build.
    The PLINK invocations otherwise mirror the ones SuSiEx would issue itself
    (``data.cpp:65-122``), so the matrices match the un-cached path.

    Returns {ancestry: ld_prefix}. Idempotent: skips a column whose LD exists.
    """
    chrom, start, end = int(locus["chrom"]), int(locus["start_bp"]), int(locus["end_bp"])
    ensure_dir(out_dir)
    ld: dict[str, Path] = {}
    for col in columns:
        ref = window
        prefix = out_dir / f"{col.name}_ld"
        ld[col.name] = prefix
        if _ld_ready(prefix):
            continue
        # the variant list SuSiEx would extract: in-window variants of the panel
        bim = pd.read_csv(ref.with_suffix(".bim"), sep=r"\s+", header=None,
                          usecols=[0, 1, 3], names=["chrom", "snp", "bp"])
        keep = bim.loc[(bim["chrom"].astype(str) == str(chrom))
                       & bim["bp"].between(start, end), "snp"]
        snp_list = Path(f"{prefix}.snp")
        keep.to_csv(snp_list, index=False, header=False)

        extract = ["--bfile", str(ref), "--keep-allele-order", "--chr", str(chrom),
                   "--extract", str(snp_list), "--maf", str(maf),
                   "--make-bed", "--out", f"{prefix}_ref"]
        if col.keep_path is not None:
            extract += ["--keep", str(col.keep_path)]
        run_plink(extract, binary=plink)
        run_plink(["--bfile", f"{prefix}_ref", "--keep-allele-order",
                   "--r", "square", "bin4", "--out", str(prefix)], binary=plink)
        run_plink(["--bfile", f"{prefix}_ref", "--keep-allele-order",
                   "--freq", "--out", f"{prefix}_frq"], binary=plink)
        junk = [snp_list, Path(f"{prefix}_ref.bed"), Path(f"{prefix}_ref.fam")]
        junk += out_dir.glob(f"{prefix.name}*.log")
        junk += out_dir.glob(f"{prefix.name}*.nosex")
        for path in junk:
            path.unlink(missing_ok=True)
        if not _ld_ready(prefix):
            raise RuntimeError(
                f"LD precompute produced no usable matrix for column {col.name} "
                f"at {prefix} (expected {', '.join(_LD_SUFFIXES)})"
            )
    return ld


def run_gwas(
    ref_prefix: Path, pheno_path: Path, out_prefix: Path, plink2: str,
    keep: Path | None = None,
) -> Path:
    """Per-ancestry quantitative GWAS on the window; write a SuSiEx sumstats file.

    ``keep`` (an FID/IID file) restricts the GWAS to one ancestry's pooled
    individuals across all sites; ``None`` uses the whole window. The phenotype
    file is the pooled one covering all enrolled individuals, so plink2 aligns to
    the kept intersection.

    Returns the sumstats path. Reformats plink2 ``.glm.linear`` (ADD rows) into
    ``chr snp bp A1 A2 beta se stat p`` with A2 = the non-effect allele.
    """
    glm = ["--bfile", str(ref_prefix), "--pheno", str(pheno_path),
           "--pheno-name", "y", "--glm", "allow-no-covars",
           "--out", str(out_prefix)]
    if keep is not None:
        glm += ["--keep", str(keep)]
    run_plink(glm, binary=plink2)
    globbed = list(out_prefix.parent.glob(out_prefix.name + "*.glm.linear"))
    if not globbed:
        raise FileNotFoundError(f"plink2 --glm produced no .glm.linear for {out_prefix}")
    g = pd.read_csv(globbed[0], sep="\t")
    g = g[g["TEST"] == "ADD"].copy()
    a1, ref, alt = g["A1"], g["REF"], g["ALT"]
    g["A2"] = np.where(a1 == alt, ref, alt)
    out = g[["#CHROM", "ID", "POS", "A1", "A2", "BETA", "SE", "T_STAT", "P"]].copy()
    out.columns = _SUMSTATS_COLS
    out = out.dropna(subset=["beta", "se", "p"])
    ss_path = out_prefix.with_suffix(".sumstats")
    out.to_csv(ss_path, sep="\t", index=False)
    return ss_path


def pooled_phenotype(cfg: SimulationConfig, inst: str, out_path: Path) -> Path:
    """Concatenate the three sites' phenotype files for one instance.

    Phenotypes are simulated and stored per site, but a pooled ancestry column
    draws its individuals from every site, so the GWAS needs one phenotype table
    spanning all of them. Individuals are disjoint across sites, so this is a
    plain row concatenation; ``--keep`` then selects the ancestry.
    """
    pheno_root = cfg.resolved_path("ground_truth_dir") / "phenotypes"
    frames = [
        pd.read_csv(pheno_root / site / f"{inst}.pheno", sep="\t")
        for site in cfg.sites
    ]
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(out_path, sep="\t", index=False)
    return out_path


# ---------------------------------------------------------------------------
# SuSiEx invocation + output parsing
# ---------------------------------------------------------------------------

def run_susiex(
    sumstats: list[Path], ld_prefixes: list[Path],
    n_gwas: list[int], chrom: int, start: int, end: int,
    out_dir: Path, out_name: str, susiex: str, plink: str,
    level: float, pval_thresh: float, threads: int = 1,
) -> None:
    """Shell out to SuSiEx for one instance across all sites.

    No ``--ref_file``: ``precompute_ld`` has already written the LD matrices, so
    SuSiEx takes its precmp path (main.cpp:204-222) and reads them directly,
    never touching a genotype panel. That is also the call shape a genuinely
    federated deployment would use — sites ship LD, the coordinator has none.
    ``--plink`` still has to be passed because the argument validator demands it
    unconditionally (main.cpp:274), but on this path it is never invoked.
    """
    def _csv(xs) -> str:
        return ",".join(str(x) for x in xs)

    ncol = len(sumstats)
    run_plink(
        ["--sst_file", _csv(sumstats),
         "--n_gwas", _csv(n_gwas),
         "--ld_file", _csv(ld_prefixes),
         "--out_dir", str(out_dir), "--out_name", out_name,
         "--chr", str(chrom), "--bp", f"{start},{end}",
         "--chr_col", _csv(["1"] * ncol), "--snp_col", _csv(["2"] * ncol),
         "--bp_col", _csv(["3"] * ncol), "--a1_col", _csv(["4"] * ncol),
         "--a2_col", _csv(["5"] * ncol), "--eff_col", _csv(["6"] * ncol),
         "--se_col", _csv(["7"] * ncol), "--pval_col", _csv(["9"] * ncol),
         "--plink", plink, "--level", str(level),
         "--pval_thresh", str(pval_thresh), "--threads", str(threads)],
        binary=susiex,
    )


def _pip_by_snp(snp_file: Path) -> dict[str, float]:
    """Map SNP id -> max PIP across all PIP(CSk) columns, from the ``.snp`` file."""
    if not snp_file.exists():
        return {}
    try:
        df = pd.read_csv(snp_file, sep="\t")
    except (pd.errors.EmptyDataError, OSError):
        return {}
    pip_cols = [c for c in df.columns if c.startswith("PIP(")]
    if "SNP" not in df.columns or not pip_cols:
        return {}
    pip = df[pip_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)
    return dict(zip(df["SNP"].astype(str), pip.astype(float)))


def parse_susiex(out_dir: Path, out_name: str, truth_snps: list[str]) -> dict:
    """Standard fine-mapping metrics from SuSiEx ``.cs`` / ``.summary`` / ``.snp``."""
    cs_file = out_dir / f"{out_name}.cs"
    sum_file = out_dir / f"{out_name}.summary"
    snp_file = out_dir / f"{out_name}.snp"
    truth = set(truth_snps)

    cs = pd.DataFrame()
    if cs_file.exists():
        try:
            cs = pd.read_csv(cs_file, sep="\t", comment="#")
        except (pd.errors.EmptyDataError, OSError):
            cs = pd.DataFrame()
    # SuSiEx can write a headerless "no credible set" note; guard on the columns.
    has_cs = not cs.empty and {"CS_ID", "SNP"}.issubset(cs.columns)

    cs_sizes: list[int] = []
    total_cs_snps = 0
    n_causal_captured = 0
    best_cs_size = None
    if has_cs:
        for _cid, grp in cs.groupby("CS_ID"):
            members = set(grp["SNP"].astype(str))
            cs_sizes.append(len(members))
            total_cs_snps += len(members)
            if members & truth:
                if best_cs_size is None or len(members) < best_cs_size:
                    best_cs_size = len(members)
        captured = set(cs["SNP"].astype(str)) & truth
        n_causal_captured = len(captured)

    # per-CS purity/length from the .summary table
    mean_purity = min_purity = np.nan
    if sum_file.exists():
        try:
            s = pd.read_csv(sum_file, sep="\t", comment="#")
            if "CS_PURITY" in s.columns and len(s):
                pur = pd.to_numeric(s["CS_PURITY"], errors="coerce")
                mean_purity, min_purity = float(pur.mean()), float(pur.min())
        except (pd.errors.EmptyDataError, OSError):
            pass

    pip_map = _pip_by_snp(snp_file)
    causal_pips = [pip_map[s] for s in truth_snps if s in pip_map]
    top_pip = max(pip_map.values()) if pip_map else np.nan

    return {
        "n_causal": len(truth_snps),
        "n_credible_sets": len(cs_sizes),
        "cs_sizes_json": json.dumps(sorted(cs_sizes)),
        "total_cs_snps": total_cs_snps,
        "n_causal_captured": n_causal_captured,
        "any_causal_captured": bool(n_causal_captured > 0),
        "best_cs_size": best_cs_size if best_cs_size is not None else np.nan,
        "causal_pip_max": max(causal_pips) if causal_pips else np.nan,
        "causal_pip_mean": float(np.mean(causal_pips)) if causal_pips else np.nan,
        "top_pip": top_pip,
        "mean_cs_purity": mean_purity,
        "min_cs_purity": min_purity,
        "converged": bool(has_cs or bool(pip_map)),
    }


# ---------------------------------------------------------------------------
# One instance
# ---------------------------------------------------------------------------

def finemap_instance(
    cfg: SimulationConfig, locus: pd.Series, arch_id: str, rep: int,
    truth_snps: list[str], columns: list[FMColumn], window: Path,
    ld: dict[str, Path], work_root: Path, susiex: str, plink: str, plink2: str,
    level: float, pval_thresh: float, keep_work: bool,
) -> dict:
    """Run GWAS -> SuSiEx -> eval for a single instance; return a result row.

    The three sites' phenotypes are pooled once, then one GWAS per ancestry
    (``--keep``-restricted to that ancestry's pooled individuals), then a single
    SuSiEx run across all ancestry columns' sumstats + precomputed LD.
    """
    t0 = time.time()
    chrom, start, end = int(locus["chrom"]), int(locus["start_bp"]), int(locus["end_bp"])
    inst = f"{locus['locus_id']}_{arch_id}_rep{rep}"
    priv = ensure_dir(work_root / inst)
    row = {
        "locus_id": locus["locus_id"], "architecture_id": arch_id, "replicate": rep,
        "stratum": locus.get("stratum", ""),
    }
    min_p: dict[str, float] = {}
    try:
        pheno = pooled_phenotype(cfg, inst, priv / "pooled.pheno")
        sumstats, ld_list, n_list = [], [], []
        for col in columns:
            ss = run_gwas(window, pheno, priv / col.name, plink2,
                          keep=col.keep_path)
            ssdf = pd.read_csv(ss, sep="\t", usecols=["snp", "p"])
            min_p[col.name] = float(ssdf["p"].min()) if len(ssdf) else np.nan
            sumstats.append(ss)
            ld_list.append(ld[col.name])  # shared read-only, built once per locus
            n_list.append(col.n)
        run_susiex(sumstats, ld_list, n_list, chrom, start, end,
                   priv, "cs", susiex, plink, level, pval_thresh)
        metrics = parse_susiex(priv, "cs", truth_snps)
        row.update(metrics)
        row["error"] = ""
    except Exception as exc:  # keep the shard alive; record the failure
        get_logger().warning("fine-map FAILED %s: %s", inst, exc)
        # Fill the whole metric schema, not just the two flags: a shard in which
        # every instance fails must still emit the columns _write_rollup reads.
        row.update(dict.fromkeys(_METRIC_COLS, np.nan))
        row.update({"any_causal_captured": False, "converged": False,
                    "error": str(exc)[:200]})
    finally:
        for col in columns:
            row[f"min_p_{col.name}"] = min_p.get(col.name, np.nan)
        row["runtime_s"] = round(time.time() - t0, 2)
        if not keep_work:
            shutil.rmtree(priv, ignore_errors=True)
    return row


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def _fm_dir(cfg: SimulationConfig) -> Path:
    return cfg.resolved_path("reports_dir") / "fine_mapping"


def run_fine_mapping(
    cfg: SimulationConfig, shard_index: int = 0, n_shards: int = 1,
    n_workers: int = 8, level: float = 0.95, pval_thresh: float = 1e-5,
    keep_work: bool = False, limit: int | None = None, maf: float = 0.005,
) -> None:
    logger = get_logger()
    repo = cfg.repo_root
    susiex = _resolve_binary("SuSiEx", repo)
    plink = _resolve_binary(cfg.tools.plink, repo)
    plink2 = _resolve_binary(cfg.tools.plink2, repo)

    fm_dir = ensure_dir(_fm_dir(cfg))
    work_root = ensure_dir(fm_dir / "work")

    # SuSiEx columns are locus-independent: plan them once and write each
    # ancestry's pooled --keep ID file up front, then reuse across all loci.
    keep_dir = work_root / "keeps"
    columns = plan_columns(
        {s: dict(sc.composition) for s, sc in cfg.sites.items()},
        cfg.superpopulations, cfg.fine_mapping.min_gwas_n,
    )
    materialize_column_keeps(columns, cfg, keep_dir)
    enrolled_keep = pooled_ids_path(cfg, keep_dir)
    logger.info("SuSiEx ancestry columns (pooled across sites, min_gwas_n=%d): %s",
                cfg.fine_mapping.min_gwas_n,
                ", ".join(f"{c.name}(n={c.n} from {'+'.join(c.sites)})"
                          for c in columns))

    loci = pd.read_csv(
        cfg.resolved_path("loci_dir") / "selected_loci.tsv", sep="\t"
    ).reset_index(drop=True)
    manifest = pd.read_csv(cfg.resolved_path("ground_truth_dir") / "causal_manifest.tsv", sep="\t")
    truth_map = {
        (r.locus_id, r.architecture_id, int(r.replicate)): str(r.causal_snp_ids).split(",")
        for r in manifest.itertuples(index=False)
    }
    architectures = build_architecture_grid(
        ncsl=cfg.architecture.ncsl, h2=cfg.architecture.h2,
        rg=cfg.architecture.rg, mode=cfg.architecture.factorial_mode,
    )
    n_rep = cfg.architecture.replicates

    if n_shards > 1:
        loci = loci.loc[(loci.index % n_shards) == shard_index].copy()
    logger.info("Fine-mapping %d loci × %d arch × %d rep (shard %d/%d, %d workers)",
                len(loci), len(architectures), n_rep, shard_index, n_shards, n_workers)

    all_rows: list[dict] = []
    for _, locus in loci.iterrows():
        # window + LD built ONCE per locus, serially, then shared read-only
        # across instances: both are phenotype-independent, and building them
        # here also keeps the workers from racing to write the same LD files.
        locus_ref_dir = ensure_dir(work_root / f"{locus['locus_id']}_refs")
        window = extract_locus_window(cfg, locus, enrolled_keep, locus_ref_dir, plink2)
        ld = precompute_ld(locus, columns, window, locus_ref_dir, plink, maf)
        jobs = []
        for arch in architectures:
            for rep in range(n_rep):
                truth = truth_map.get((locus["locus_id"], arch.architecture_id, rep))
                if truth is None:
                    continue
                jobs.append((arch.architecture_id, rep, truth))
                if limit is not None and len(all_rows) + len(jobs) >= limit:
                    break
            if limit is not None and len(all_rows) + len(jobs) >= limit:
                break
        rows = Parallel(n_jobs=n_workers, backend="loky")(
            delayed(finemap_instance)(
                cfg, locus, aid, rep, truth, columns, window, ld, work_root,
                susiex, plink, plink2, level, pval_thresh, keep_work)
            for (aid, rep, truth) in jobs
        )
        all_rows.extend(rows)
        if not keep_work:
            shutil.rmtree(locus_ref_dir, ignore_errors=True)
        logger.info("  locus %s done (%d instances, cumulative %d)",
                    locus["locus_id"], len(jobs), len(all_rows))
        if limit is not None and len(all_rows) >= limit:
            break

    df = pd.DataFrame(all_rows)
    if n_shards > 1:
        tag = f"part{shard_index:03d}of{n_shards:03d}"
        df.to_csv(fm_dir / f"fm_results.{tag}.tsv", sep="\t", index=False)
        logger.info("Shard %d/%d wrote %d result rows (%s)", shard_index, n_shards, len(df), tag)
    else:
        df.to_csv(fm_dir / "fm_results.tsv", sep="\t", index=False)
        logger.info("Wrote %d result rows -> %s", len(df), fm_dir / "fm_results.tsv")
        _write_rollup(cfg, df)


def merge_fm_results(cfg: SimulationConfig, n_shards: int) -> None:
    """REDUCE: concat per-shard result parts, sort, write fm_results.tsv + rollup."""
    logger = get_logger()
    fm_dir = _fm_dir(cfg)
    parts = []
    for r in range(n_shards):
        p = fm_dir / f"fm_results.part{r:03d}of{n_shards:03d}.tsv"
        if not p.exists():
            raise FileNotFoundError(f"Missing fine-mapping result part: {p}")
        parts.append(pd.read_csv(p, sep="\t"))
    df = pd.concat(parts, ignore_index=True).sort_values(
        ["locus_id", "architecture_id", "replicate"]).reset_index(drop=True)
    df.to_csv(fm_dir / "fm_results.tsv", sep="\t", index=False)
    for r in range(n_shards):
        (fm_dir / f"fm_results.part{r:03d}of{n_shards:03d}.tsv").unlink()
    logger.info("Merged %d shards -> %d rows", n_shards, len(df))
    _write_rollup(cfg, df)


def _write_rollup(cfg: SimulationConfig, df: pd.DataFrame) -> None:
    """Aggregate power/credible-set metrics by architecture and by stratum."""
    if df.empty:
        return
    fm_dir = _fm_dir(cfg)
    arch = build_architecture_grid(
        ncsl=cfg.architecture.ncsl, h2=cfg.architecture.h2,
        rg=cfg.architecture.rg, mode=cfg.architecture.factorial_mode)
    meta = pd.DataFrame([{"architecture_id": a.architecture_id, "ncsl": a.ncsl,
                          "h2_target": a.h2, "rg": a.rg} for a in arch])
    d = df.merge(meta, on="architecture_id", how="left")

    def _agg(g: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "n_instances": len(g),
            "power_any_causal": g["any_causal_captured"].mean(),
            "mean_n_cs": g["n_credible_sets"].mean(),
            "median_best_cs_size": g["best_cs_size"].median(),
            "mean_causal_pip": g["causal_pip_max"].mean(),
        })

    for by, name in [(["ncsl", "h2_target", "rg"], "by_architecture"),
                     (["stratum", "rg"], "by_stratum_rg")]:
        cols = [c for c in by if c in d.columns]
        if cols:
            d.groupby(cols).apply(_agg, include_groups=False).reset_index().to_csv(
                fm_dir / f"fm_rollup_{name}.tsv", sep="\t", index=False)
    get_logger().info("Wrote rollups -> %s/fm_rollup_*.tsv", fm_dir)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="FedFM SuSiEx fine-mapping")
    parser.add_argument("--config", default="config/simulation_config.yaml")
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--n-shards", type=int, default=1)
    parser.add_argument("--merge", action="store_true",
                        help="REDUCE: merge per-shard result parts, then exit")
    parser.add_argument("--n-workers", type=int, default=8)
    parser.add_argument("--level", type=float, default=0.95,
                        help="credible-set coverage level (SuSiEx --level)")
    parser.add_argument("--pval-thresh", type=float, default=1e-5,
                        help="SuSiEx marginal p-value filter for credible sets")
    parser.add_argument("--maf", type=float, default=0.005,
                        help="MAF filter applied to the LD panel (SuSiEx default)")
    parser.add_argument("--keep-work", action="store_true",
                        help="keep per-instance GWAS/LD scratch (debug)")
    parser.add_argument("--limit", type=int, default=None,
                        help="stop after N instances (smoke test)")
    args = parser.parse_args(list(argv) if argv is not None else None)
    cfg = load_config(args.config)
    setup_logging(log_dir=cfg.resolved_path("logs_dir"))

    if args.merge:
        merge_fm_results(cfg, args.n_shards)
    else:
        if args.n_shards > 1 and not (0 <= args.shard_index < args.n_shards):
            parser.error("need 0 <= --shard-index < --n-shards")
        run_fine_mapping(cfg, args.shard_index, args.n_shards, args.n_workers,
                         args.level, args.pval_thresh, args.keep_work, args.limit,
                         args.maf)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
