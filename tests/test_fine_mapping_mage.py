"""Offline unit and integration tests for MAGE + 1000 Genomes real-data validation."""

from __future__ import annotations

import gzip
import json
import struct
import zlib
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

from appfl_bio_suite.experiments.fine_mapping.aggregator import FineMappingAggregator
from appfl_bio_suite.experiments.fine_mapping.dataset import get_dataset
from appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping import _resolve_binary
from appfl_bio_suite.experiments.fine_mapping.fedfm.utils import read_bim
from appfl_bio_suite.experiments.fine_mapping.figures.paper_plots import write_figures
from appfl_bio_suite.experiments.fine_mapping.real_data import (
    benchmark,
    federation_builder,
    figures_mage,
    mage,
    panel,
    vcf_slice,
)
from appfl_bio_suite.experiments.fine_mapping.simulation.cohort import write_plink1
from appfl_bio_suite.experiments.fine_mapping.trainer import SiteFineMappingTrainer

matplotlib.use("Agg", force=True)


def _bgzf_block(payload: bytes) -> bytes:
    """Wrap ``payload`` into one valid BGZF block with BSIZe subfield."""
    comp = zlib.compressobj(6, zlib.DEFLATED, -zlib.MAX_WBITS)
    cdata = comp.compress(payload) + comp.flush()
    crc = zlib.crc32(payload) & 0xFFFFFFFF
    isize = len(payload) & 0xFFFFFFFF
    extra = b"BC" + struct.pack("<HH", 2, 0)  # placeholder replaced once total length known
    bsize = 12 + len(extra) + len(cdata) + 8 - 1
    extra = b"BC" + struct.pack("<HH", 2, bsize)
    header = struct.pack("<4sIBB", b"\x1f\x8b\x08\x04", 0, 0, 255) + struct.pack("<H", len(extra))
    return header + extra + cdata + struct.pack("<II", crc, isize)


def _make_tbi_for_block(
    ref_name: str,
    bin_id: int,
    beg_voff: int,
    end_voff: int,
    linear_voff: int,
) -> bytes:
    names_raw = ref_name.encode("ascii") + b"\x00"
    header = struct.pack("<4s8i", b"TBI\x01", 1, 2, 1, 2, 2, ord("#"), 0, len(names_raw))
    bins_bytes = struct.pack("<i", 1) + struct.pack("<Ii2Q", bin_id, 1, beg_voff, end_voff)
    lin_bytes = struct.pack("<i", 1) + struct.pack("<Q", linear_voff)
    return header + names_raw + bins_bytes + lin_bytes


def test_vcf_slice_tbi_and_bgzf_slicing():
    """Tabix parsing, bin math, virtual offset merging, and BGZF extraction stay exact offline."""
    assert vcf_slice.normalize_chrom("chr4") == "4"
    assert vcf_slice.normalize_chrom("04") == "4"
    assert vcf_slice.normalize_chrom(22) == "22"
    with pytest.raises(ValueError, match="autosomal"):
        vcf_slice.normalize_chrom("X")

    bins = vcf_slice.reg2bins(1000, 2000)
    assert bins[0] == 0 and 4681 in bins
    with pytest.raises(ValueError, match="Invalid 0-based"):
        vcf_slice.reg2bins(2000, 2000)

    vcf_header = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    vcf_body = b"chr4\t1500\trs1\tA\tG\t.\tPASS\t.\n"
    b0 = _bgzf_block(vcf_header)
    b1 = _bgzf_block(vcf_body)
    eof = _bgzf_block(b"")
    stream = b0 + b1 + eof

    rec_voff_beg = len(b0) << 16
    rec_voff_end = (len(b0) + len(b1)) << 16
    tbi_bytes = _make_tbi_for_block("chr4", 4681, rec_voff_beg, rec_voff_end, rec_voff_beg)
    tbi = vcf_slice.parse_tbi(gzip.compress(tbi_bytes))
    assert tbi.resolve_ref_index(4) == 0
    with pytest.raises(KeyError):
        tbi.resolve_ref_index(19)

    ranges = vcf_slice.query_byte_ranges(tbi, "4", 1000, 2000, file_size=len(stream))
    assert len(ranges) == 1
    r0, r1 = ranges[0]
    assert r0 == len(b0)
    recovered = b"".join(
        payload for _pos, payload in vcf_slice.iter_bgzf_blocks(stream[r0 : r1 + 1])
    )
    assert recovered == vcf_body


def _write_synthetic_mage_raw(
    raw_dir: Path,
    *,
    include_swap: bool = True,
    n_per_pop: int = 2,
    n_variants: int = 24,
    rng_seed: int = 202611,
) -> dict:
    """Create a miniature 26-population MAGE v1.0 raw directory with planted eQTL signals."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(rng_seed)

    pop_to_group = {
        "CEU": "EUR",
        "GBR": "EUR",
        "IBS": "EUR",
        "TSI": "EUR",
        "FIN": "EUR",
        "ASW": "AFR",
        "ACB": "AFR",
        "GWD": "AFR",
        "YRI": "AFR",
        "LWK": "AFR",
        "MSL": "AFR",
        "ESN": "AFR",
        "PUR": "AMR",
        "CLM": "AMR",
        "MXL": "AMR",
        "PEL": "AMR",
        "CHB": "EAS",
        "JPT": "EAS",
        "CHS": "EAS",
        "CDX": "EAS",
        "KHV": "EAS",
        "GIH": "SAS",
        "PJL": "SAS",
        "BEB": "SAS",
        "STU": "SAS",
        "ITU": "SAS",
    }

    rows: list[dict] = []
    idx = 0
    for pop, grp in pop_to_group.items():
        for k in range(n_per_pop):
            if pop == "CEU" and k == 0 and include_swap:
                rows.append(
                    {
                        "SRA_accession": mage.SWAP_AFFECTED_SRA,
                        "internal_libraryID": "GM11919",
                        "sample_coriellID": "GM11919",
                        "sample_kgpID": mage.SWAP_AFFECTED_SAMPLE,
                        "batch": 3,
                        "continentalGroup": "EUR",
                        "population": "CEU",
                        "sex": "XY",
                        "numReads": 30_000_000,
                        "totalVolume_ul": 40,
                        "RNAQubitConc_ng-ul": 150.0,
                        "RNAQubitTotalAmount_ng": 6000.0,
                        "RIN": 9.5,
                    }
                )
            elif pop == "GBR" and k == 0 and include_swap:
                rows.append(
                    {
                        "SRA_accession": "SRR19762247",
                        "internal_libraryID": "HG00237-7-2",
                        "sample_coriellID": "HG00237",
                        "sample_kgpID": mage.SWAP_REPLACEMENT_SAMPLE,
                        "batch": 7,
                        "continentalGroup": "EUR",
                        "population": "GBR",
                        "sex": "XX",
                        "numReads": 31_000_000,
                        "totalVolume_ul": 42,
                        "RNAQubitConc_ng-ul": 160.0,
                        "RNAQubitTotalAmount_ng": 6200.0,
                        "RIN": 9.6,
                    }
                )
            else:
                sid = f"SAMP{idx:04d}"
                rows.append(
                    {
                        "SRA_accession": f"SRR900{idx:04d}",
                        "internal_libraryID": f"LIB{idx:04d}",
                        "sample_coriellID": sid,
                        "sample_kgpID": sid,
                        "batch": (idx % 5) + 1,
                        "continentalGroup": grp,
                        "population": pop,
                        "sex": "XX" if idx % 2 == 0 else "XY",
                        "numReads": 25_000_000,
                        "totalVolume_ul": 40,
                        "RNAQubitConc_ng-ul": 150.0,
                        "RNAQubitTotalAmount_ng": 6000.0,
                        "RIN": 9.3,
                    }
                )
                idx += 1

    sm = pd.DataFrame(rows)
    sm.to_csv(raw_dir / "sample.metadata.MAGE.v1.0.txt", sep="\t", index=False)
    seq_extra = sm.copy()
    if include_swap:
        dropped_row = seq_extra[seq_extra["sample_kgpID"] == mage.SWAP_REPLACEMENT_SAMPLE].iloc[0]
        dropped_dict = dropped_row.to_dict()
        dropped_dict["SRA_accession"] = mage.SWAP_DROPPED_SRA
        dropped_dict["batch"] = 3
        seq_extra = pd.concat([seq_extra, pd.DataFrame([dropped_dict])], ignore_index=True)
    seq_extra.to_csv(raw_dir / "sequencing_library.metadata.MAGE.v1.0.txt", sep="\t", index=False)

    sample_ids = sm["sample_kgpID"].tolist()
    n_samples = len(sample_ids)

    # Plant 5 autosomal eGenes (across all 4 panel strata + 1 HLA exclusion check)
    genes = [
        ("ENSG00000100001.1", "GENE_COLOC", 1, 10_000_000, 0.96, 0.12, "EUR", 1, 12, "page_coloc"),
        ("ENSG00000100002.1", "GENE_FST", 2, 20_000_000, 0.10, 0.72, "AFR", 1, 6, "high_fst"),
        ("ENSG00000100003.1", "GENE_MULTI", 3, 30_000_000, 0.10, 0.15, "EAS", 4, 4, "multi_cs"),
        ("ENSG00000100004.1", "GENE_MOD", 4, 40_000_000, 0.05, 0.12, "AMR", 1, 12, "moderate_cs"),
        ("ENSG00000100005.1", "GENE_HLA", 6, 30_000_000, 0.99, 0.85, "EUR", 2, 3, "hla_drop"),
    ]

    # Covariates: sex + PC1 + PEER1
    cov_rows = [
        ["sex", *sm["sex"].tolist()],
        ["PC1", *[f"{v:.6f}" for v in rng.normal(size=n_samples)]],
        ["PEER1", *[f"{v:.6f}" for v in rng.normal(size=n_samples)]],
    ]
    with gzip.open(raw_dir / "eQTL_covariates.tab.gz", "wt") as fh:
        fh.write("\t".join(["id", *sample_ids]) + "\n")
        for r in cov_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    # Genotypes & planted expression per non-HLA gene so phenotypes carry real genetic signal
    genos_by_gene: dict[str, tuple[np.ndarray, np.ndarray, list[tuple[str, str]]]] = {}
    Z_num = np.column_stack(
        [
            (sm["sex"].to_numpy() == "XY").astype(np.float64),
            np.asarray(cov_rows[1][1:], dtype=np.float64),
            np.asarray(cov_rows[2][1:], dtype=np.float64),
        ]
    )
    bed_lines = ["\t".join(["#Chr", "start", "end", "ID", *sample_ids]) + "\n"]
    eqtl_sum_rows = []
    cs_sum_rows = []
    sig_rows = []
    fst_rows = []
    page_rows = []

    for eid, sym, chrom, tss, pph4, fst, fst_pop, n_cs, cs_sz, _tag in genes:
        positions = tss - 60_000 + np.arange(n_variants) * 5_000
        alleles: list[tuple[str, str]] = [("A", "G")] * n_variants
        alleles[2] = ("A", "ATG")  # exercise biallelic indel preservation
        # Site-correlated MAF difference + strong shared signal at variant index 8
        G_chr = np.empty((n_samples, n_variants), dtype=np.int8)
        for j in range(n_variants):
            p_base = 0.25 if j != 8 else 0.35
            for r_i, grp in enumerate(sm["continentalGroup"]):
                freq = min(0.85, max(0.15, p_base + (0.20 if grp == "AFR" else -0.08)))
                G_chr[r_i, j] = rng.binomial(2, freq)
        genos_by_gene[eid] = (G_chr, positions, alleles)

        y = (
            0.85 * G_chr[:, 8]
            + 0.35 * Z_num[:, 0]
            + 0.25 * Z_num[:, 1]
            + rng.normal(0, 0.9, n_samples)
        )
        bed_lines.append(
            "\t".join([f"chr{chrom}", str(tss), str(tss + 1), eid, *[f"{v:.6f}" for v in y]]) + "\n"
        )
        eqtl_sum_rows.append(
            {
                "ensemblID": eid,
                "geneSymbol": sym,
                "variantsTested": n_variants,
                "significantAssocs": 100 - chrom,
                "is_eGene": True,
                "finemappedCredibleSets": n_cs,
            }
        )
        for c_idx in range(1, n_cs + 1):
            cs_sum_rows.append(
                {
                    "ensemblID": eid,
                    "geneSymbol": sym,
                    "credibleSet": f"L{c_idx}",
                    "numVariants": cs_sz,
                    "coverage": 0.96,
                    "minCorr": 0.85,
                    "meanCorr": 0.92,
                    "medianCorr": 0.91,
                }
            )
        lead_pos = int(positions[8])
        lead_kgp = f"{chrom}:{lead_pos}:{alleles[8][0]}:{alleles[8][1]}"
        sig_rows.append(
            {
                "variantChrom": f"chr{chrom}",
                "variantPosition": lead_pos,
                "variantRef": alleles[8][0],
                "variantAlt": alleles[8][1],
                "variant_kgpID": lead_kgp,
                "variant_rsID": f"rs{chrom}008",
                "ensemblID": eid,
                "geneSymbol": sym,
                "variantPIP": 0.91,
                "variantCredibleSet": "L1",
            }
        )
        fst_rows.append(
            {
                "geneID": eid.split(".")[0],
                "continentalGroup": fst_pop,
                "mean_wcFst": fst,
                "DEpadj": 0.01,
            }
        )
        if pph4 >= 0.8:
            page_rows.append(
                {
                    "GWAS_catalog_accessionID": "GCST000001",
                    "traitDescription": "Mock trait",
                    "GWAS_sentinelVariant": f"chr{chrom}:{lead_pos}:A:G",
                    "eQTL_gene_ensemblID": eid,
                    "eQTL_gene_geneSymbol": sym,
                    "eQTL_variantChrom": f"chr{chrom}",
                    "eQTL_variantPosition": lead_pos,
                    "eQTL_variantRef": "A",
                    "eQTL_variantAlt": "G",
                    "eQTL_variant_kgpID": lead_kgp,
                    "eQTL_variant_rsID": f"rs{chrom}008",
                    "eQTL_geneCS": "L1",
                    "GWAS_variantChrom": f"chr{chrom}",
                    "GWAS_variantPosition": lead_pos,
                    "GWAS_variantRef": "A",
                    "GWAS_variantAlt": "G",
                    "GWAS_variant_kgpID": lead_kgp,
                    "GWAS_variant_rsID": f"rs{chrom}008",
                    "GWAS_CS": "L1",
                    "PP.H0": 0.0,
                    "PP.H1": 0.0,
                    "PP.H2": 0.0,
                    "PP.H3": 1.0 - pph4,
                    "PP.H4": pph4,
                }
            )

    with gzip.open(raw_dir / "inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz", "wt") as fh:
        fh.writelines(bed_lines)
    pd.DataFrame(eqtl_sum_rows).to_csv(
        raw_dir / "eQTL_summary.MAGE.v1.0.txt.gz", sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(cs_sum_rows).to_csv(
        raw_dir / "eQTL_finemapping.credibleSet_summary.MAGE.v1.0.txt.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(sig_rows).to_csv(
        raw_dir / "eQTL_finemapping.significantAssociations.MAGE.v1.0.txt.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(fst_rows).to_csv(raw_dir / "Finemapped_Fst.csv", index=False)
    pd.DataFrame(page_rows).to_csv(
        raw_dir / "PAGE.eQTL_colocalization.moderate_results.MAGE.v1.0.txt.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    return {"samples": sample_ids, "genos_by_gene": genos_by_gene, "genes": genes}


def test_mage_sample_swap_and_fwl_projection(tmp_path):
    """Sample-swap audit excludes NA11919 -> N-1 donors, and FWL equals joint OLS adjustment."""
    raw = tmp_path / "raw"
    cache = tmp_path / "cache"
    _write_synthetic_mage_raw(raw, include_swap=True, n_per_pop=2)

    audit_clean = mage.audit_and_remediate_sample_swap(raw, exclude_swapped=True)
    audit_swap = mage.audit_and_remediate_sample_swap(raw, exclude_swapped=False)
    assert audit_clean["swap_detected_in_release"] is True
    assert audit_clean["excluded_samples"] == [mage.SWAP_AFFECTED_SAMPLE]
    assert audit_clean["n_clean_samples"] == audit_swap["n_clean_samples"] - 1
    assert audit_clean["n_populations"] == 26
    assert audit_clean["n_continental_groups"] == 5

    paths = mage.prepare_mage_cache(raw, cache, exclude_swapped=True)
    assert paths["audit"].is_file()
    gene_meta, expr, expr_samples, Z_cov = mage.load_cached_gene_expression(
        cache, [row["ensembl_id"] for row in gene_meta_list(raw)], exclude_swapped=True
    )
    assert mage.SWAP_AFFECTED_SAMPLE not in expr_samples
    assert expr.shape == (len(gene_meta), len(expr_samples))
    assert Z_cov.shape[0] == len(expr_samples)
    assert len(gene_meta) == 5

    # Exact Frisch-Waugh-Lovell algebraic identity: beta(y_perp ~ X_perp) == beta(y ~ 1 + X + Z).
    rng = np.random.default_rng(7)
    n, m, k = 120, 4, 5
    X = rng.normal(size=(n, m))
    Z = rng.normal(size=(n, k))
    y = 1.5 + 0.7 * X[:, 1] - 0.9 * Z[:, 0] + 0.4 * Z[:, 3] + rng.normal(0, 0.25, n)

    X_perp, y_perp = mage.residualize_design(X, y, Z, include_intercept=True)
    assert np.allclose(X_perp.T @ Z, 0.0, atol=1e-11)
    assert np.allclose(y_perp.T @ Z, 0.0, atol=1e-11)
    for j in range(m):
        beta_fwl = float((X_perp[:, j] @ y_perp) / (X_perp[:, j] @ X_perp[:, j]))
        design = np.column_stack([np.ones(n), X[:, j], Z])
        beta_full = float(np.linalg.lstsq(design, y, rcond=None)[0][1])
        assert abs(beta_fwl - beta_full) < 1e-11


def gene_meta_list(raw_dir: Path) -> list[dict]:
    df = pd.read_csv(raw_dir / "eQTL_summary.MAGE.v1.0.txt.gz", sep="\t")
    return [{"ensembl_id": eid} for eid in df["ensemblID"].tolist()]


def test_panel_selection_non_overlapping_and_strata(tmp_path):
    """HLA loci are excluded and all strata remain non-overlapping by min_gap_bp."""
    raw = tmp_path / "raw"
    _write_synthetic_mage_raw(raw, include_swap=True, n_per_pop=1)
    cat = panel.build_mage_gene_catalog(raw, window_half_bp=250_000)
    assert bool(cat.loc[cat["gene_symbol"] == "GENE_HLA", "is_hla"].iloc[0]) is True

    selected = panel.select_mage_panel(cat, n_per_stratum=1, min_gap_bp=1_000_000, exclude_hla=True)
    assert "GENE_HLA" not in selected["gene_symbol"].tolist()
    assert set(selected["stratum"]) == set(panel.STRATUM_ORDER)
    assert selected["all_cs_within_window"].all()

    cs_sum, cs_vars = panel.load_published_mage_cs(raw, selected["ensembl_id"].tolist())
    assert len(cs_sum) > 0 and len(cs_vars) == 4


def _mock_extract_from_synthetic(genos_by_gene: dict, sample_ids: list[str], genes: list[tuple]):
    """Build deterministic local extract_locus_plink stand-in preserving GRCh38 IDs and indels."""
    by_chr_pos = {(c, tss): (eid, *genos_by_gene[eid]) for eid, _s, c, tss, *_rest in genes}

    def _fake_extract(
        chrom,
        start_bp,
        end_bp,
        out_prefix,
        keep_samples_file,
        cache_dir,
        plink2_bin="plink2",
        maf=0.01,
        snps_only=True,
        rm_dup_mode="exclude-mismatch",
        overwrite=False,
        **_kwargs,
    ):
        out_prefix = Path(out_prefix)
        c = int(vcf_slice.normalize_chrom(chrom))
        match = next(
            (v for (mc, mtss), v in by_chr_pos.items() if mc == c and start_bp <= mtss <= end_bp),
            None,
        )
        assert match is not None
        _eid, G_chr, positions, alleles = match
        keep_df = pd.read_csv(
            keep_samples_file, sep="\t", header=None, names=["FID", "IID"], dtype=str
        )
        keep_idx = [sample_ids.index(iid) for iid in keep_df["IID"].tolist()]
        G_sub = G_chr[keep_idx]
        kept_iids = [sample_ids[i] for i in keep_idx]

        write_plink1(out_prefix, G_sub, c, positions, kept_iids, alleles=alleles)
        with Path(f"{out_prefix}.fam").open("w", encoding="utf-8") as fh:
            for iid in kept_iids:
                fh.write(f"0\t{iid}\t0\t0\t0\t-9\n")
        with Path(f"{out_prefix}.bim").open("w", encoding="utf-8") as fh:
            for pos, (a1, a2) in zip(positions, alleles, strict=True):
                fh.write(f"{c}\t{c}:{int(pos)}:{a1}:{a2}\t0\t{int(pos)}\t{a1}\t{a2}\n")
        return out_prefix

    return _fake_extract


def test_mage_federation_and_shipped_pipeline_round_trip(tmp_path, monkeypatch):
    """End-to-end bundle creation, APPFL site->aggregator round-trip, and benchmark parity."""
    pytest.importorskip("MultiSuSiE")
    try:
        plink_bin = _resolve_binary("plink", Path(__file__).resolve().parents[1])
        plink2_bin = _resolve_binary("plink2", Path(__file__).resolve().parents[1])
    except FileNotFoundError:
        pytest.skip("PLINK binaries unavailable")
    assert plink_bin and plink2_bin

    raw = tmp_path / "raw"
    cache = tmp_path / "cache"
    fed_root = tmp_path / "federation"
    meta_synth = _write_synthetic_mage_raw(raw, include_swap=True, n_per_pop=10, rng_seed=99)

    monkeypatch.setattr(
        "urllib.request.urlopen",
        lambda *a, **kw: pytest.fail("Unexpected network access in offline test"),
    )
    monkeypatch.setattr(federation_builder, "stage_mage_data", lambda *a, **kw: {})
    monkeypatch.setattr(
        federation_builder,
        "extract_locus_plink",
        _mock_extract_from_synthetic(
            meta_synth["genos_by_gene"], meta_synth["samples"], meta_synth["genes"]
        ),
    )
    manifest = federation_builder.build_mage_federation(
        out_root=fed_root,
        raw_dir=raw,
        cache_dir=cache,
        n_per_stratum=1,
        window_half_bp=250_000,
        min_gap_bp=1_000_000,
        maf=0.01,
        snps_only=False,
        exclude_swapped=True,
    )
    assert set(manifest["sites"]) == set(federation_builder.SITE_ORDER)
    assert manifest["n_loci"] == 4
    assert (fed_root / "drs_registry.json").is_file()
    assert (fed_root / "pipeline_config.yaml").is_file()

    # Strict single-arg APPFL logger to guard paper_plots formatting regression too
    class StrictLogger:
        def info(self, msg):
            assert isinstance(msg, str)

        def warning(self, msg):
            assert isinstance(msg, str)

        def debug(self, msg):
            assert isinstance(msg, str)

        def error(self, msg):
            assert isinstance(msg, str)

    strict_log = StrictLogger()

    # Shipped contracts & DUO/BIM indel verification on every site
    payloads = {}
    for site in federation_builder.SITE_ORDER:
        sdata = fed_root / site / "data"
        assert (sdata / "site_genotypes.bed").is_file()
        assert (sdata / "site_manifest.tsv").is_file()
        duo = json.loads((sdata / "DATA_USE.json").read_text(encoding="utf-8"))
        assert duo["permission"] == "DUO:0000042"
        bim = read_bim(f"{sdata / 'site_genotypes.bim'}")
        assert any(a == "ATG" for a in bim["a2"].astype(str)), "Expected preserved indel allele"

        dataset, _ = get_dataset(
            data_dir=sdata,
            site_id=site,
            data_use_request=benchmark.MAGE_DATA_USE_REQUEST,
            enforce_data_use=True,
        )
        trainer = SiteFineMappingTrainer(
            train_dataset=dataset,
            train_configs={
                "trainer_output_dirname": str(tmp_path / "trainer_out" / site),
                "pops": list(benchmark.POPULATION_ORDER),
                "locus_limit": 2,
                "uplink_gram_dtype": "float64",
            },
            logger=strict_log,
            client_id=site,
        )
        trainer.train()
        payloads[site] = trainer.get_parameters()

    ms_out = tmp_path / "agg_ms"
    agg = FineMappingAggregator(
        aggregator_configs={
            "backend": "multisusie",
            "population_order": list(benchmark.POPULATION_ORDER),
            "level": 0.95,
            "maf": 0.005,
            "n_signals": 3,
            "max_iter": 200,
            "tol": 1e-3,
            "multisusie_rho": 0.75,
            "multisusie_min_abs_corr": 0.5,
            "n_workers": 2,
            "make_figures": False,
            "output_dir": str(ms_out),
            "causal_manifest": None,
        },
        logger=strict_log,
    )
    state = agg.aggregate(payloads)
    assert state["fine_mapping_ready"].item() == 1
    res_ms = pd.read_csv(ms_out / "data" / "fed_fm_results.tsv", sep="\t")
    assert len(res_ms) == 2 and res_ms["converged"].all()
    assert res_ms["any_causal_captured"].isna().all()
    assert (
        pd.read_csv(ms_out / "data" / "fed_fm_site_summary.csv")["N_VARIANTS_FLIPPED"] == 0
    ).all()

    # Check paper_plots.write_figures with StrictLogger when truth columns are NaN
    written_pp = write_figures(
        res_ms,
        ms_out / "paper_figs",
        loci=pd.read_csv(fed_root / "site_atlantic" / "data" / "selected_loci.tsv", sep="\t"),
        logger=strict_log,
    )
    assert isinstance(written_pp, list)

    # Full SuSiEx + FWL + naive control + reporting + figures_mage suite when SuSiEx is available
    try:
        susiex_bin = _resolve_binary("SuSiEx", Path(__file__).resolve().parents[1])
    except FileNotFoundError:
        susiex_bin = ""
    if not susiex_bin:
        return

    val_out = tmp_path / "val_out"
    report = benchmark.run_validation_suite(
        clean_root=fed_root,
        out_root=val_out,
        swap_root=None,
        maf=0.005,
        n_signals=3,
        max_iter=200,
        tol=1e-4,
        multisusie_max_iter=200,
        multisusie_tol=1e-3,
        max_workers=1,
        run_single_ancestries=True,
        run_solo_sites=True,
        locus_limit=None,
    )
    parity = pd.read_csv(val_out / "reports" / "parity_summary.tsv", sep="\t").set_index(
        "comparison"
    )
    assert parity.loc["susiex_fed_vs_central", "worst_max_pip_abs_diff"] < 1e-8
    assert parity.loc["susiex_fed_vs_central", "frac_exact_cs_identical"] == 1.0
    assert parity.loc["multisusie_fed_vs_central", "worst_max_pip_abs_diff"] < 1e-10
    assert parity.loc["multisusie_fed_vs_central", "frac_exact_cs_identical"] == 1.0

    ndf = pd.read_csv(val_out / "reports" / "naive_ld_distortion.tsv", sep="\t")
    assert (ndf["max_ld_abs_diff"] > 0.0).all()
    assert (ndf["max_z_abs_diff"] > 0.0).all()
    assert len(report["arm_summary"]) >= 9

    figs = figures_mage.render_all_mage_figures(fed_root, val_out)
    assert len(figs) == 10
    for fig_path in figs:
        assert fig_path.is_file() and fig_path.stat().st_size > 1000
