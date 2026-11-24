"""Build realistic multi-site FedFM bundles from MAGE expression + 1000G GRCh38 genotypes."""

from __future__ import annotations

import json
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from appfl_bio_suite.core.experiments import repo_root
from appfl_bio_suite.core.ga4gh.drs import build_registry
from appfl_bio_suite.core.ga4gh.duo import DataUseProfile, write_profile

from ..fedfm.fine_mapping import _resolve_binary
from ..fedfm.utils import read_bed_variants, read_fam
from .mage import (
    load_cached_gene_expression,
    prepare_mage_cache,
    residualize_covariates,
    site_local_expression_and_covariates,
    stage_mage_data,
)
from .panel import build_mage_gene_catalog, load_published_mage_cs, select_mage_panel
from .vcf_slice import extract_locus_plink

SITE_ORDER: tuple[str, ...] = ("site_atlantic", "site_africa_eurasia", "site_asiapacific")

# Disjoint assignment of all 26 MAGE/1000G populations across three geographically
# structured clinical biobanks so EUR/AFR/AMR/SAS span multiple sites while EAS
# concentrates at site_asiapacific.
SITE_POPULATIONS: dict[str, tuple[str, ...]] = {
    "site_atlantic": ("CEU", "GBR", "IBS", "ASW", "ACB", "GWD", "PUR", "CLM", "GIH"),
    "site_africa_eurasia": ("TSI", "FIN", "YRI", "LWK", "MSL", "ESN", "PJL", "BEB", "MXL"),
    "site_asiapacific": ("CHB", "JPT", "CHS", "CDX", "KHV", "STU", "ITU", "PEL"),
}

ARCHITECTURE_ID = "mage_eqtl"
REPLICATE_ID = 0


def assign_mage_sites(
    metadata: pd.DataFrame,
    site_populations: Mapping[str, Sequence[str]] = SITE_POPULATIONS,
) -> pd.DataFrame:
    """Assign each MAGE donor to an institutional site via its 1000G population label."""
    pop_to_site: dict[str, str] = {}
    for site, pops in site_populations.items():
        for pop in pops:
            if pop in pop_to_site:
                raise ValueError(f"Population {pop!r} assigned to multiple sites")
            pop_to_site[pop] = site

    observed_pops = set(metadata["population"].astype(str))
    missing_pops = observed_pops - set(pop_to_site)
    if missing_pops:
        raise ValueError(f"Unassigned MAGE populations: {sorted(missing_pops)}")

    assigned = metadata.copy()
    assigned["site_id"] = assigned["population"].map(pop_to_site)
    return assigned


def _merge_locus_plink_sets(
    prefixes: Sequence[Path],
    merged_prefix: Path,
    plink_bin: str,
) -> Path:
    if not prefixes:
        raise ValueError("Cannot merge an empty list of locus PLINK filesets")
    merged_prefix.parent.mkdir(parents=True, exist_ok=True)
    if len(prefixes) == 1:
        for ext in (".bed", ".bim", ".fam"):
            shutil.copy2(f"{prefixes[0]}{ext}", f"{merged_prefix}{ext}")
        return merged_prefix

    merge_list = merged_prefix.parent / f"{merged_prefix.name}_merge_list.txt"
    merge_list.write_text(
        "\n".join(f"{p}.bed {p}.bim {p}.fam" for p in prefixes[1:]) + "\n",
        encoding="utf-8",
    )
    cmd = [
        plink_bin,
        "--bfile",
        str(prefixes[0]),
        "--merge-list",
        str(merge_list),
        "--make-bed",
        "--threads",
        "1",
        "--out",
        str(merged_prefix),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        raise RuntimeError(f"plink --merge-list failed:\n{res.stderr or res.stdout}")
    return merged_prefix


def _subset_site_plink(
    pooled_prefix: Path,
    keep_ids: pd.DataFrame,
    out_prefix: Path,
    plink2_bin: str,
) -> Path:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    keep_file = out_prefix.parent / f".{out_prefix.name}.keep"
    keep_ids[["FID", "IID"]].to_csv(keep_file, sep="\t", header=False, index=False)
    cmd = [
        plink2_bin,
        "--bfile",
        str(pooled_prefix),
        "--keep",
        str(keep_file),
        "--make-bed",
        "--threads",
        "1",
        "--out",
        str(out_prefix),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    keep_file.unlink(missing_ok=True)
    if res.returncode != 0:
        raise RuntimeError(f"plink2 --keep failed for {out_prefix}:\n{res.stderr or res.stdout}")
    return out_prefix


def build_mage_federation(
    out_root: Path,
    raw_dir: Path | None = None,
    cache_dir: Path | None = None,
    n_per_stratum: int = 4,
    window_half_bp: int = 250_000,
    min_gap_bp: int = 1_000_000,
    maf: float = 0.01,
    snps_only: bool = False,
    covariate_mode: str = "all",
    exclude_swapped: bool = True,
    site_populations: Mapping[str, Sequence[str]] = SITE_POPULATIONS,
    drs_hostname: str = "mage-local.appfl-bio-suite.org",
    overwrite: bool = False,
    phenotype: str = "residualized",
    selection_window_half_bp: int | None = None,
    covariate_source: str = "published",
    n_site_factors: int = 15,
    panel_gene_ids: Sequence[str] | None = None,
) -> dict:
    """Materialize a complete multi-site FedFM real-data validation bundle from MAGE + 1000G.

    ``phenotype="residualized"`` ships expression with the covariates already regressed out
    by the coordinator (the original design). ``phenotype="raw"`` ships the released
    inverse-normal TMM expression untouched, so no phenotype preprocessing happens centrally
    and covariate adjustment is left to the federated FWL step (``covariates.py``). Either
    way each site's bundle carries its own rows of the covariate table as ``covariates.tsv``.

    ``selection_window_half_bp`` (default: ``window_half_bp``) is the window the panel is
    chosen with; a wider ``window_half_bp`` then widens the analysis around the same loci.
    ``covariate_source="site_local"`` replaces the published pooled covariates with ones
    each site derives from its own data (``mage.site_local_expression_and_covariates``).
    """
    if phenotype not in ("residualized", "raw"):
        raise ValueError(f"phenotype must be 'residualized' or 'raw', got {phenotype!r}")
    if covariate_source not in ("published", "site_local"):
        raise ValueError("covariate_source must be 'published' or 'site_local'")
    root = repo_root()
    mage_base = root / "local" / "data" / "mage"
    raw_dir = Path(raw_dir) if raw_dir is not None else mage_base / "raw"
    cache_dir = Path(cache_dir) if cache_dir is not None else mage_base / "cache"
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    plink_bin = _resolve_binary("plink", root)
    plink2_bin = _resolve_binary("plink2", root)

    stage_mage_data(raw_dir)
    cache_paths = prepare_mage_cache(
        raw_dir, cache_dir, exclude_swapped=exclude_swapped, overwrite=False
    )
    audit = json.loads(cache_paths["audit"].read_text(encoding="utf-8"))
    samples_df = pd.read_csv(cache_paths["samples"], sep="\t", dtype={"sample_kgpID": str})
    samples_df = assign_mage_sites(samples_df, site_populations=site_populations)

    selection_half = int(selection_window_half_bp or window_half_bp)
    catalog = build_mage_gene_catalog(raw_dir, window_half_bp=selection_half)
    if panel_gene_ids is None:
        panel = select_mage_panel(
            catalog,
            n_per_stratum=n_per_stratum,
            min_gap_bp=min_gap_bp,
            exclude_hla=True,
        )
    else:
        if not len(panel_gene_ids) or len(set(panel_gene_ids)) != len(panel_gene_ids):
            raise ValueError("Expected a nonempty panel of distinct gene IDs")
        panel = catalog.set_index("ensembl_id").loc[list(panel_gene_ids)].reset_index()
        panel["stratum"] = "prespecified"
        panel["locus_id"] = "H-" + panel["ensembl_id"]
    panel["start_bp"] = (panel["tss"] - int(window_half_bp)).clip(lower=1)
    panel["end_bp"] = panel["tss"] + int(window_half_bp)
    ensembl_ids = panel["ensembl_id"].tolist()
    meta_sub, expr_raw, expr_samples, Z_all = load_cached_gene_expression(
        cache_dir, ensembl_ids, exclude_swapped=exclude_swapped
    )
    _ = meta_sub

    # Slice covariate matrix according to selected adjustment mode
    cov_table = pd.read_csv(cache_paths["cov"], sep="\t")
    if covariate_mode == "all":
        cov_cols = [c for c in cov_table.columns if c != "sample_id"]
    elif covariate_mode == "pcs_sex":
        cov_cols = [c for c in cov_table.columns if c == "sex" or c.startswith("PC")]
    elif covariate_mode == "peers_sex":
        cov_cols = [c for c in cov_table.columns if c == "sex" or c.startswith("PEER")]
    else:
        raise ValueError(f"Unsupported covariate_mode {covariate_mode!r}")
    if covariate_source == "site_local":
        site_of = samples_df.set_index("sample_kgpID")["site_id"]
        expr_raw, cov_table = site_local_expression_and_covariates(
            cache_dir,
            ensembl_ids,
            site_of,
            [s for s in site_populations if s in set(site_of)],
            exclude_swapped=exclude_swapped,
            n_factors=n_site_factors,
        )
        cov_cols = [c for c in cov_table.columns if c != "sample_id"]
    Z_use = cov_table[cov_cols].to_numpy(dtype=np.float64)
    _ = Z_all

    # Align residuals row-wise across panel genes (n_samples x n_loci -> transpose back)
    expr_resid = residualize_covariates(expr_raw.T, Z_use, include_intercept=True).T

    # Write shared keep file ("0" FID matches NYGC 1000G VCF -> PLINK2 default FID)
    meta_dir = out_root / "meta"
    meta_dir.mkdir(parents=True, exist_ok=True)
    keep_all = meta_dir / "clean_samples.keep"
    pd.DataFrame({"FID": "0", "IID": expr_samples}).to_csv(
        keep_all, sep="\t", header=False, index=False
    )

    # Extract clean PLINK1 fileset per selected locus and merge
    loci_plink_dir = cache_dir / "plink" / ("clean" if exclude_swapped else "with_swap")
    locus_prefixes: list[Path] = []
    for row in panel.itertuples():
        lp = loci_plink_dir / f"{row.locus_id}_{window_half_bp}"
        extract_locus_plink(
            chrom=int(row.chrom),
            start_bp=int(row.start_bp),
            end_bp=int(row.end_bp),
            out_prefix=lp,
            keep_samples_file=keep_all,
            cache_dir=cache_dir,
            plink2_bin=plink2_bin,
            maf=maf,
            snps_only=snps_only,
            overwrite=overwrite,
        )
        locus_prefixes.append(lp)

    pooled_dir = meta_dir / "pooled_plink"
    pooled_prefix = pooled_dir / "pooled_genotypes"
    _merge_locus_plink_sets(locus_prefixes, pooled_prefix, plink_bin=plink_bin)

    pooled_bim = pd.read_csv(
        f"{pooled_prefix}.bim",
        sep="\t",
        header=None,
        names=["chrom", "snp_id", "cm", "bp", "a1", "a2"],
        dtype={"chrom": str, "snp_id": str, "bp": np.int64, "a1": str, "a2": str},
    )
    pooled_fam = read_fam(f"{pooled_prefix}.fam")
    if set(pooled_fam["IID"]) != set(expr_samples):
        raise RuntimeError("Pooled PLINK .fam sample IDs do not match MAGE sample IDs")

    # Verify zero missing genotypes in pooled binary (FedFM _gram integrality invariant)
    n_check = min(len(pooled_bim), 50)
    if n_check > 0:
        G_chk = read_bed_variants(pooled_prefix, np.arange(n_check), len(pooled_fam))
        if np.isnan(G_chk).any():
            raise RuntimeError("Unexpected missing genotypes after --geno 0.0 extraction")

    ref_variants = pooled_bim[["snp_id", "chrom", "bp", "a1", "a2"]].copy()
    selected_loci = panel[["locus_id", "chrom", "start_bp", "end_bp", "stratum"]].copy()

    # Re-index samples_df to exact expression sample order
    samples_idx = samples_df.set_index("sample_kgpID").loc[expr_samples].reset_index()
    samples_idx["FID"] = "0"
    samples_idx["IID"] = samples_idx["sample_kgpID"]
    samples_idx["superpopulation"] = samples_idx["continentalGroup"]

    sample_to_expr_col = {sid: i for i, sid in enumerate(expr_samples)}

    site_names = [s for s in site_populations if s in set(samples_idx["site_id"])]
    site_counts: dict[str, dict[str, int]] = {}
    for site in site_names:
        site_data = out_root / site / "data"
        pheno_dir = site_data / "phenotypes"
        pheno_dir.mkdir(parents=True, exist_ok=True)

        sub = samples_idx[samples_idx["site_id"] == site].copy()
        site_prefix = site_data / "site_genotypes"
        _subset_site_plink(pooled_prefix, sub[["FID", "IID"]], site_prefix, plink2_bin=plink2_bin)
        site_fam = read_fam(f"{site_prefix}.fam")
        # Order site manifest & phenotypes to match site PLINK .fam order deterministically
        sub_ordered = sub.set_index("IID").loc[site_fam["IID"].tolist()].reset_index()
        sub_ordered[["FID", "IID", "superpopulation"]].to_csv(
            site_data / "site_manifest.tsv", sep="\t", index=False
        )
        ref_variants.to_csv(site_data / "reference_variants.tsv", sep="\t", index=False)
        selected_loci.to_csv(site_data / "selected_loci.tsv", sep="\t", index=False)

        col_indices = [sample_to_expr_col[sid] for sid in sub_ordered["IID"].tolist()]
        site_cov = (
            cov_table.astype({"sample_id": str})
            .set_index("sample_id")
            .loc[sub_ordered["IID"].tolist(), cov_cols]
            .reset_index(drop=True)
        )
        site_cov.insert(0, "IID", sub_ordered["IID"].tolist())
        site_cov.to_csv(site_data / "covariates.tsv", sep="\t", index=False)
        site_expr = expr_raw if phenotype == "raw" else expr_resid
        for g_idx, row in enumerate(panel.itertuples()):
            inst = f"{row.locus_id}_{ARCHITECTURE_ID}_rep{REPLICATE_ID}"
            pheno = pd.DataFrame(
                {
                    "FID": sub_ordered["FID"].tolist(),
                    "IID": sub_ordered["IID"].tolist(),
                    "y": site_expr[g_idx, col_indices],
                }
            )
            pheno.to_csv(pheno_dir / f"{inst}.pheno", sep="\t", index=False)

        site_counts[site] = {
            k: int(v) for k, v in sub_ordered["superpopulation"].value_counts().items()
        }

    # Build DRS registry over site bundles (DATA_USE.json is excluded from Merkle checksums)
    # and stamp each site's DATA_USE.json with DUO:0000042 + its bundle DRS URI.
    registry = build_registry(out_root, hostname=drs_hostname, experiment="fine-mapping")
    drs_path = out_root / "drs_registry.json"
    registry.save(drs_path)
    bundle_map = registry.site_bundles()
    for site in site_names:
        profile = DataUseProfile(
            dataset_id=f"mage-1000g-chr-panel-{site}",
            drs_uri=bundle_map[site].self_uri,
            site=site,
            permission="DUO:0000042",
            description=(
                "MAGE v1.0 LCL RNA-seq inverse-normal TMM expression paired with 1000 Genomes "
                "NYGC 30x GRCh38 genotypes (open-access general research use)."
            ),
        )
        write_profile(profile, out_root / site / "data" / "DATA_USE.json")

    # Save coordinator-side metadata and published MAGE benchmark truth references
    pub_cs_sum, pub_cs_vars = load_published_mage_cs(raw_dir, ensembl_ids)
    panel.to_csv(meta_dir / "panel_manifest.tsv", sep="\t", index=False)
    samples_idx[
        ["FID", "IID", "site_id", "continentalGroup", "population", "sex", "batch", "RIN"]
    ].to_csv(meta_dir / "cohort_metadata.tsv", sep="\t", index=False)
    cov_table.to_csv(meta_dir / "covariates.tsv", sep="\t", index=False)
    pub_cs_sum.to_csv(meta_dir / "published_cs_summary.tsv", sep="\t", index=False)
    pub_cs_vars.to_csv(meta_dir / "published_cs_variants.tsv", sep="\t", index=False)
    (meta_dir / "sample_swap_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    np.savez(
        meta_dir / "expression_panel.npz",
        ensembl_ids=np.asarray(ensembl_ids, dtype=object),
        locus_ids=np.asarray(panel["locus_id"].tolist(), dtype=object),
        samples=np.asarray(expr_samples, dtype=object),
        expr_raw=expr_raw,
        expr_resid=expr_resid,
        covariate_cols=np.asarray(cov_cols, dtype=object),
    )

    manifest = {
        "n_loci": int(len(panel)),
        "n_samples": int(len(expr_samples)),
        "n_variants_pooled": int(len(pooled_bim)),
        "exclude_swapped": bool(exclude_swapped),
        "covariate_mode": covariate_mode,
        "phenotype": phenotype,
        "covariate_source": covariate_source,
        "window_half_bp": int(window_half_bp),
        "selection_window_half_bp": selection_half,
        "maf": float(maf),
        "snps_only": bool(snps_only),
        "architecture_id": ARCHITECTURE_ID,
        "replicate_id": REPLICATE_ID,
        "sites": site_counts,
        "drs_registry": str(drs_path),
    }
    import yaml

    # Compatibility config consumed by `core.loopback._scenario_maf` so
    # `appfl-bio-suite run --experiment fine-mapping --loopback --data-root <out_root>`
    # picks up the exact MAF threshold and inference defaults automatically.
    pipeline_cfg = {
        "sites": {s: site_counts[s] for s in site_names},
        "fine_mapping": {
            "maf": float(maf),
            "min_gwas_n": 0,
            "keep_ambiguous": True,
            "n_signals": 10,
            "max_iter": 1000,
            "tol": 1e-6,
        },
        "mage": {
            "exclude_swapped": bool(exclude_swapped),
            "covariate_mode": covariate_mode,
            "phenotype": phenotype,
            "covariate_source": covariate_source,
            "window_half_bp": int(window_half_bp),
            "selection_window_half_bp": selection_half,
        },
    }
    (out_root / "pipeline_config.yaml").write_text(
        yaml.safe_dump(pipeline_cfg, sort_keys=False), encoding="utf-8"
    )
    (out_root / "mage_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
