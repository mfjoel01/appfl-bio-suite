"""Programmatic locus panel selection over real MAGE autosomal eGenes."""

from __future__ import annotations

import collections
import csv
import gzip
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

STRATUM_ORDER: tuple[str, ...] = ("page_coloc", "high_fst", "multi_cs", "moderate_cs")

# Extended MHC / HLA region on GRCh38 chromosome 6 (excluded by default to avoid
# hyper-polymorphic multi-megabase LD complexity during cross-site validation).
HLA_CHROM = 6
HLA_START_BP = 28_510_120
HLA_END_BP = 33_480_577


def _is_hla(chrom: int, tss: int) -> bool:
    return chrom == HLA_CHROM and HLA_START_BP <= tss <= HLA_END_BP


def build_mage_gene_catalog(
    raw_dir: Path,
    window_half_bp: int = 250_000,
) -> pd.DataFrame:
    """Assemble autosomal finemapped eGene metadata with Fst, CS, and PAGE annotations."""
    raw_dir = Path(raw_dir)
    genes: dict[str, dict] = {}
    with gzip.open(raw_dir / "inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz", "rt") as fh:
        fh.readline()
        for line in fh:
            parts = line.split("\t", 5)
            c = parts[0].removeprefix("chr")
            if c.isdigit() and 1 <= int(c) <= 22:
                genes[parts[3]] = {
                    "chrom": int(c),
                    "tss": int(parts[1]),
                    "gene_id_base": parts[3].split(".")[0],
                }

    summary: dict[str, dict] = {}
    with gzip.open(raw_dir / "eQTL_summary.MAGE.v1.0.txt.gz", "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            summary[row["ensemblID"]] = row

    cs_sum: dict[str, list[dict]] = collections.defaultdict(list)
    with gzip.open(raw_dir / "eQTL_finemapping.credibleSet_summary.MAGE.v1.0.txt.gz", "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cs_sum[row["ensemblID"]].append(row)

    max_dist: dict[str, int] = collections.defaultdict(int)
    top_pip: dict[str, tuple[float, str, str, bool]] = {}
    any_indel: dict[str, bool] = collections.defaultdict(bool)
    sig_path = raw_dir / "eQTL_finemapping.significantAssociations.MAGE.v1.0.txt.gz"
    with gzip.open(sig_path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            eid = row["ensemblID"]
            if eid not in genes:
                continue
            dist = abs(int(row["variantPosition"]) - genes[eid]["tss"])
            if dist > max_dist[eid]:
                max_dist[eid] = dist
            is_indel = len(row["variantRef"]) > 1 or len(row["variantAlt"]) > 1
            if is_indel:
                any_indel[eid] = True
            pip = float(row["variantPIP"])
            if eid not in top_pip or pip > top_pip[eid][0]:
                top_pip[eid] = (pip, row["variant_kgpID"], row["variant_rsID"], is_indel)

    fst_by_gene: dict[str, dict[str, float]] = collections.defaultdict(dict)
    with open(raw_dir / "Finemapped_Fst.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            fst_by_gene[row["geneID"]][row["continentalGroup"]] = float(row["mean_wcFst"])

    page_by_gene: dict[str, list[tuple[float, str, str]]] = collections.defaultdict(list)
    page_path = raw_dir / "PAGE.eQTL_colocalization.moderate_results.MAGE.v1.0.txt.gz"
    with gzip.open(page_path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            page_by_gene[row["eQTL_gene_ensemblID"]].append(
                (float(row["PP.H4"]), row["traitDescription"], row["eQTL_geneCS"])
            )

    rows: list[dict] = []
    for eid, g in genes.items():
        if eid not in cs_sum:
            continue
        s_row = summary.get(eid, {})
        sizes = sorted(int(c["numVariants"]) for c in cs_sum[eid])
        fsts = fst_by_gene.get(g["gene_id_base"], {})
        max_fst_pop = max(fsts, key=fsts.get) if fsts else ""
        max_fst = float(fsts[max_fst_pop]) if max_fst_pop else 0.0
        p_hits = sorted(page_by_gene.get(eid, []), reverse=True)
        lead = top_pip.get(eid, (0.0, "", "", False))
        rows.append(
            {
                "ensembl_id": eid,
                "gene_id_base": g["gene_id_base"],
                "gene_symbol": str(s_row.get("geneSymbol", eid)),
                "chrom": int(g["chrom"]),
                "tss": int(g["tss"]),
                "start_bp": max(1, int(g["tss"]) - int(window_half_bp)),
                "end_bp": int(g["tss"]) + int(window_half_bp),
                "variants_tested": int(s_row.get("variantsTested", 0)),
                "significant_assocs": int(s_row.get("significantAssocs", 0)),
                "n_published_cs": len(sizes),
                "min_cs_size": sizes[0],
                "median_cs_size": sizes[len(sizes) // 2],
                "max_cs_size": sizes[-1],
                "max_cs_dist_from_tss": int(max_dist.get(eid, 0)),
                "all_cs_within_window": int(max_dist.get(eid, 0)) <= int(window_half_bp),
                "max_wc_fst": max_fst,
                "max_fst_pop": max_fst_pop,
                "page_max_pph4": float(p_hits[0][0]) if p_hits else 0.0,
                "page_trait": str(p_hits[0][1]) if p_hits else "",
                "page_eqtl_cs": str(p_hits[0][2]) if p_hits else "",
                "published_lead_pip": float(lead[0]),
                "published_lead_kgp": str(lead[1]),
                "published_lead_rsid": str(lead[2]),
                "published_lead_is_indel": bool(lead[3]),
                "has_indel_cs_variant": bool(any_indel.get(eid, False)),
                "is_hla": _is_hla(int(g["chrom"]), int(g["tss"])),
            }
        )
    return pd.DataFrame(rows).sort_values(["chrom", "tss", "ensembl_id"]).reset_index(drop=True)


def _overlaps_existing(row: pd.Series, chosen: Sequence[dict], min_gap_bp: int) -> bool:
    c, t = int(row["chrom"]), int(row["tss"])
    return any(item["chrom"] == c and abs(item["tss"] - t) < min_gap_bp for item in chosen)


def _pick_top(
    pool: pd.DataFrame,
    chosen: list[dict],
    used_ids: set[str],
    stratum: str,
    quota: int,
    min_gap_bp: int,
) -> None:
    picked = 0
    for _, row in pool.iterrows():
        if picked >= quota:
            break
        eid = str(row["ensembl_id"])
        if eid in used_ids or _overlaps_existing(row, chosen, min_gap_bp):
            continue
        item = row.to_dict()
        item["stratum"] = stratum
        chosen.append(item)
        used_ids.add(eid)
        picked += 1


def select_mage_panel(
    catalog: pd.DataFrame,
    n_per_stratum: int = 4,
    min_gap_bp: int = 1_000_000,
    exclude_hla: bool = True,
) -> pd.DataFrame:
    """Select pairwise non-overlapping eGenes across four complementary validation strata."""
    eligible = catalog[catalog["all_cs_within_window"]].copy()
    if exclude_hla and "is_hla" in eligible:
        eligible = eligible[~eligible["is_hla"]].copy()

    chosen: list[dict] = []
    used: set[str] = set()

    # 1. PAGE GWAS-colocalized eGenes ranked by PP.H4 then significance strength
    page_pool = eligible[eligible["page_max_pph4"] >= 0.85].sort_values(
        ["page_max_pph4", "significant_assocs", "ensembl_id"],
        ascending=[False, False, True],
    )
    _pick_top(page_pool, chosen, used, "page_coloc", n_per_stratum, min_gap_bp)

    # 2. High cross-ancestry allele-frequency differentiation (Fst >= 0.50)
    fst_pool = eligible[eligible["max_wc_fst"] >= 0.50].sort_values(
        ["max_wc_fst", "significant_assocs", "ensembl_id"],
        ascending=[False, False, True],
    )
    _pick_top(fst_pool, chosen, used, "high_fst", n_per_stratum, min_gap_bp)

    # 3. High allelic heterogeneity (multiple independent eQTL credible sets)
    multi_pool = eligible[eligible["n_published_cs"] >= 4].sort_values(
        ["n_published_cs", "significant_assocs", "ensembl_id"],
        ascending=[False, False, True],
    )
    _pick_top(multi_pool, chosen, used, "multi_cs", n_per_stratum, min_gap_bp)

    # 4. Moderate single-CS eGenes with wider single-cohort uncertainty (size 8..35)
    mod_pool = eligible[
        (eligible["n_published_cs"] == 1)
        & (eligible["max_wc_fst"] < 0.25)
        & (eligible["max_cs_size"].between(8, 35))
    ].sort_values(
        ["significant_assocs", "published_lead_pip", "ensembl_id"],
        ascending=[False, False, True],
    )
    _pick_top(mod_pool, chosen, used, "moderate_cs", n_per_stratum, min_gap_bp)

    if len(chosen) < len(STRATUM_ORDER) * n_per_stratum:
        # Backfill any shortfall on small synthetic fixtures deterministically
        remaining = eligible.sort_values(
            ["significant_assocs", "ensembl_id"], ascending=[False, True]
        )
        for stratum in STRATUM_ORDER:
            have = sum(1 for c in chosen if c["stratum"] == stratum)
            if have < n_per_stratum:
                _pick_top(
                    remaining,
                    chosen,
                    used,
                    stratum,
                    n_per_stratum - have,
                    min_gap_bp,
                )

    panel = pd.DataFrame(chosen).sort_values(["chrom", "tss", "ensembl_id"]).reset_index(drop=True)
    panel["locus_id"] = [
        f"M{idx:02d}-{str(sym).replace('_', '-')}" for idx, sym in enumerate(panel["gene_symbol"])
    ]
    ordered_cols = [
        "locus_id",
        "stratum",
        "ensembl_id",
        "gene_symbol",
        "chrom",
        "tss",
        "start_bp",
        "end_bp",
        "n_published_cs",
        "min_cs_size",
        "median_cs_size",
        "max_cs_size",
        "max_wc_fst",
        "max_fst_pop",
        "page_max_pph4",
        "page_trait",
        "published_lead_kgp",
        "published_lead_rsid",
        "published_lead_pip",
        "significant_assocs",
    ]
    extra = [c for c in panel.columns if c not in ordered_cols]
    return panel[ordered_cols + extra]


def load_published_mage_cs(
    raw_dir: Path,
    ensembl_ids: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load published MAGE CS summaries and CS variant PIPs for selected ``ensembl_ids``."""
    raw_dir = Path(raw_dir)
    wanted = set(ensembl_ids)
    cs_sum = pd.read_csv(
        raw_dir / "eQTL_finemapping.credibleSet_summary.MAGE.v1.0.txt.gz", sep="\t"
    )
    cs_vars = pd.read_csv(
        raw_dir / "eQTL_finemapping.significantAssociations.MAGE.v1.0.txt.gz", sep="\t"
    )
    return (
        cs_sum[cs_sum["ensemblID"].isin(wanted)].copy().reset_index(drop=True),
        cs_vars[cs_vars["ensemblID"].isin(wanted)].copy().reset_index(drop=True),
    )
