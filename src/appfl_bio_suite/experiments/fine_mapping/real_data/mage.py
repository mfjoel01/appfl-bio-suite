"""MAGE v1.0 dataset staging, sample-swap remediation, caching, and covariate adjustment."""

from __future__ import annotations

import gzip
import hashlib
import json
import struct
import urllib.request
import zlib
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from .vcf_slice import http_content_length, http_get_bytes

MAGE_ZENODO_RECORD_ID = "10535719"
MAGE_ZENODO_API_URL = f"https://zenodo.org/api/records/{MAGE_ZENODO_RECORD_ID}"
MAGE_DATA_ZIP_URL = (
    f"https://zenodo.org/api/records/{MAGE_ZENODO_RECORD_ID}/files/MAGE.v1.0.data.zip/content"
)
MAGE_PUBLICATION_DOI = "10.1038/s41586-024-07708-2"

CONTINENTAL_GROUPS: tuple[str, ...] = ("AFR", "EUR", "EAS", "SAS", "AMR")
POPULATION_ORDER: tuple[str, ...] = ("EUR", "AFR", "AMR", "EAS", "SAS")

SWAP_AFFECTED_SAMPLE = "NA11919"
SWAP_AFFECTED_SRA = "SRR19762653"
SWAP_REPLACEMENT_SAMPLE = "HG00237"
SWAP_DROPPED_SRA = "SRR19762530"
SWAP_CLEAN_SRA = "SRR19762247"

# Verified file sizes and SHA-256 digests from Zenodo record 10535719 (MAGE v1.0).
MAGE_FILES: dict[str, tuple[int, str]] = {
    "sample.metadata.MAGE.v1.0.txt": (
        58645,
        "d3729189a1443e0b6e09caf6a13c3edaa4ff1a3bb552e147af4ffd385559c1e6",
    ),
    "sequencing_library.metadata.MAGE.v1.0.txt": (
        62631,
        "f03a87e5a42010e48e575d96208a486de0703c3f103af8e16c49725753c60fa4",
    ),
    "inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz": (
        49334484,
        "388e92200d1a7dd1bef24ef8b2e12b63f1558d7e72629f3f514def90df61209f",
    ),
    "eQTL_covariates.tab.gz": (
        184413,
        "c36e18f1068dc0497f64401737babdc3f170f93480de01aad8a1786a9f19748e",
    ),
    "eQTL_summary.MAGE.v1.0.txt.gz": (
        283864,
        "3597db30f22e262c755c7dfb687faa2def0e7045a0b0163db18d4c7b8099e110",
    ),
    "eQTL_finemapping.credibleSet_summary.MAGE.v1.0.txt.gz": (
        529257,
        "f474fc30588b87949a072ebc7df4c8cabfed4abcd8313e808ca0049c0a85732e",
    ),
    "eQTL_finemapping.significantAssociations.MAGE.v1.0.txt.gz": (
        4340911,
        "28673e882796907e12a96c9620c5fc663f81daee6cbeab88d66072b316e55fab",
    ),
    "eQTL_FastQTL_results.permutation_pass.MAGE.v1.0.txt.gz": (
        1725994,
        "80ee188ecc3a66a57ccc2a5eddacf3d657aaedc6e8eeb9e4c3ea8e69f26bc9f1",
    ),
    "PAGE.eQTL_colocalization.moderate_results.MAGE.v1.0.txt.gz": (
        5691,
        "82f9745778cf3be860e7463b3e82807447b326ce96659b3c57150285cb7740db",
    ),
    "PAGE.eQTL_colocalization.moderate_eQTL_CS.MAGE.v1.0.txt.gz": (
        41721,
        "8c38aeab954c4868ef60c2139588a4cddf7c485174bc6ddd00026ad0f3ce79d2",
    ),
    "Finemapped_Fst.csv": (
        2486104,
        "0c4916817acac86b43fd287e8c93c4aecf0389a3ec0a598c8ab494ff25dc3379",
    ),
}


def sha256_file(path: Path) -> str:
    """Compute hex SHA-256 checksum of ``path``."""
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def _parse_zip_central_directory(url: str, total_size: int) -> dict[str, tuple[int, int, int, int]]:
    """Locate member offsets in a remote ZIP archive via HTTP Range on the EOCD/CD."""
    tail_len = min(131_072, total_size)
    tail = http_get_bytes(url, total_size - tail_len, total_size - 1)
    eocd_sig = b"PK\x05\x06"
    idx = tail.rfind(eocd_sig)
    if idx < 0:
        raise RuntimeError(f"ZIP end-of-central-directory marker not found at {url}")
    eocd = tail[idx : idx + 22]
    _sig, _d1, _d2, _n1, n_entries, cd_size, cd_offset, _clen = struct.unpack("<IHHHHIIH", eocd)
    # Handle Zip64 EOCD locator if standard fields are sentinel 0xFFFFFFFF
    if cd_offset == 0xFFFFFFFF or cd_size == 0xFFFFFFFF:
        loc_idx = tail.rfind(b"PK\x06\x07", 0, idx)
        if loc_idx < 0:
            raise RuntimeError("Zip64 locator missing in remote ZIP")
        _, _d_eocd, zip64_off, _n_disks = struct.unpack("<IIQI", tail[loc_idx : loc_idx + 20])
        z64 = http_get_bytes(url, zip64_off, zip64_off + 55)
        cd_size, cd_offset = struct.unpack("<QQ", z64[40:56])

    cd_bytes = http_get_bytes(url, cd_offset, cd_offset + cd_size - 1)
    members: dict[str, tuple[int, int, int, int]] = {}
    pos = 0
    while pos + 46 <= len(cd_bytes):
        if cd_bytes[pos : pos + 4] != b"PK\x01\x02":
            break
        fields = struct.unpack("<IHHHHHHIIIHHHHHII", cd_bytes[pos : pos + 46])
        comp_method = fields[4]
        comp_size = fields[8]
        uncomp_size = fields[9]
        fname_len, extra_len, comm_len = fields[10], fields[11], fields[12]
        local_off = fields[16]
        fname = cd_bytes[pos + 46 : pos + 46 + fname_len].decode("utf-8")
        extra = cd_bytes[pos + 46 + fname_len : pos + 46 + fname_len + extra_len]
        # Resolve Zip64 extra field if needed
        if 0xFFFFFFFF in (comp_size, uncomp_size, local_off):
            ep = 0
            vals = []
            while ep + 4 <= len(extra):
                hid, hlen = struct.unpack("<HH", extra[ep : ep + 4])
                hdata = extra[ep + 4 : ep + 4 + hlen]
                if hid == 0x0001:
                    vals = list(struct.unpack(f"<{hlen // 8}Q", hdata))
                    vi = 0
                    if uncomp_size == 0xFFFFFFFF and vi < len(vals):
                        uncomp_size = vals[vi]
                        vi += 1
                    if comp_size == 0xFFFFFFFF and vi < len(vals):
                        comp_size = vals[vi]
                        vi += 1
                    if local_off == 0xFFFFFFFF and vi < len(vals):
                        local_off = vals[vi]
                    break
                ep += 4 + hlen
        members[fname] = (comp_method, local_off, comp_size, uncomp_size)
        pos += 46 + fname_len + extra_len + comm_len
    return members


def _fetch_zip_member(url: str, entry: tuple[int, int, int, int]) -> bytes:
    comp_method, local_off, comp_size, _uncomp_size = entry
    lhdr = http_get_bytes(url, local_off, local_off + 29)
    if lhdr[:4] != b"PK\x03\x04":
        raise RuntimeError(f"Bad local ZIP header at offset {local_off}")
    fn_len, ex_len = struct.unpack("<HH", lhdr[26:30])
    data_start = local_off + 30 + fn_len + ex_len
    raw = http_get_bytes(url, data_start, data_start + comp_size - 1)
    if comp_method == 0:
        return raw
    if comp_method == 8:
        return zlib.decompress(raw, -zlib.MAX_WBITS)
    raise RuntimeError(f"Unsupported ZIP compression method {comp_method}")


def stage_mage_data(
    raw_dir: Path,
    filenames: Iterable[str] | None = None,
    verify_checksums: bool = True,
) -> dict[str, Path]:
    """Ensure required MAGE v1.0 tables exist locally with verified SHA-256 digests."""
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    wanted = list(filenames) if filenames is not None else list(MAGE_FILES)

    missing: list[str] = []
    staged: dict[str, Path] = {}
    for name in wanted:
        target = raw_dir / name
        exp_size, exp_sha = MAGE_FILES[name]
        if target.is_file() and target.stat().st_size == exp_size:
            if not verify_checksums or sha256_file(target) == exp_sha:
                staged[name] = target
                continue
        missing.append(name)

    if not missing:
        return staged

    # Query remote ZIP central directory first, then fall back to direct record files
    zip_size = http_content_length(MAGE_DATA_ZIP_URL)
    zip_members = _parse_zip_central_directory(MAGE_DATA_ZIP_URL, zip_size)
    by_basename = {Path(k).name: k for k in zip_members}

    record_urls: dict[str, str] | None = None
    for name in missing:
        target = raw_dir / name
        exp_size, exp_sha = MAGE_FILES[name]
        if name in by_basename:
            payload = _fetch_zip_member(MAGE_DATA_ZIP_URL, zip_members[by_basename[name]])
        else:
            if record_urls is None:
                with urllib.request.urlopen(MAGE_ZENODO_API_URL, timeout=60) as resp:
                    meta = json.loads(resp.read().decode("utf-8"))
                record_urls = {f["key"]: f["links"]["self"] for f in meta.get("files", [])}
            if name not in record_urls:
                raise FileNotFoundError(
                    f"{name!r} not found in MAGE archive or Zenodo record {MAGE_ZENODO_RECORD_ID}"
                )
            payload = http_get_bytes(record_urls[name])

        if len(payload) != exp_size:
            raise ValueError(f"{name}: expected {exp_size} bytes, got {len(payload)}")
        got_sha = hashlib.sha256(payload).hexdigest()
        if verify_checksums and got_sha != exp_sha:
            raise ValueError(f"{name}: SHA-256 mismatch ({got_sha} != {exp_sha})")
        tmp = target.with_suffix(target.suffix + ".tmp")
        tmp.write_bytes(payload)
        tmp.replace(target)
        staged[name] = target

    return staged


def audit_and_remediate_sample_swap(
    raw_dir: Path,
    exclude_swapped: bool = True,
) -> dict:
    """Detect and optionally remediate the documented MAGE sample-swap (NA11919/HG00237).

    In MAGE v1.0, library ``GM11919`` (``SRR19762653``, annotated as CEU ``NA11919``,
    batch 3, XY) was derived in reality from ``HG00237``, while batch-3 library
    ``SRR19762530`` labeled ``HG00237`` was dropped and replaced by batch-7 library
    ``HG00237-7-2`` (``SRR19762247``, GBR, XX). Because Zenodo v1.0 retains the
    mislabeled ``NA11919`` column alongside ``HG00237``, keeping ``NA11919`` duplicates
    donor ``HG00237`` under two ancestries/sexes; default remediation removes
    ``NA11919`` leaving 730 genetically distinct individuals across 26 populations.
    """
    raw_dir = Path(raw_dir)
    sm = pd.read_csv(
        raw_dir / "sample.metadata.MAGE.v1.0.txt",
        sep="\t",
        dtype={"sample_kgpID": str},
    )
    seq_meta_path = raw_dir / "sequencing_library.metadata.MAGE.v1.0.txt"
    sl = (
        pd.read_csv(seq_meta_path, sep="\t", dtype={"sample_kgpID": str})
        if seq_meta_path.is_file()
        else pd.DataFrame()
    )

    with gzip.open(raw_dir / "inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz", "rt") as fh:
        bed_samples = fh.readline().rstrip("\n").split("\t")[4:]
    cov_header = (
        pd.read_csv(raw_dir / "eQTL_covariates.tab.gz", sep="\t", nrows=0).columns[1:].tolist()
    )
    if set(bed_samples) != set(cov_header) or set(bed_samples) != set(sm["sample_kgpID"]):
        raise ValueError("MAGE expression, covariate, and metadata sample sets do not align")

    na11919_rows = sm[sm["sample_kgpID"] == SWAP_AFFECTED_SAMPLE]
    hg00237_rows = sm[sm["sample_kgpID"] == SWAP_REPLACEMENT_SAMPLE]
    seq_dropped = (
        sl[sl["SRA_accession"] == SWAP_DROPPED_SRA]
        if not sl.empty and "SRA_accession" in sl
        else pd.DataFrame()
    )

    na11919_present = not na11919_rows.empty and (SWAP_AFFECTED_SAMPLE in bed_samples)
    excluded = [SWAP_AFFECTED_SAMPLE] if (exclude_swapped and na11919_present) else []
    clean_sm = sm[~sm["sample_kgpID"].isin(excluded)].copy()
    clean_samples = [sid for sid in bed_samples if sid not in set(excluded)]

    raw_pop_counts = sm["continentalGroup"].value_counts().to_dict()
    clean_pop_counts = clean_sm["continentalGroup"].value_counts().to_dict()
    clean_subpop_counts = clean_sm["population"].value_counts().to_dict()

    return {
        "exclude_swapped": bool(exclude_swapped),
        "n_raw_samples": int(len(sm)),
        "n_clean_samples": int(len(clean_sm)),
        "n_populations": int(clean_sm["population"].nunique()),
        "n_continental_groups": int(clean_sm["continentalGroup"].nunique()),
        "swap_detected_in_release": bool(na11919_present),
        "excluded_samples": excluded,
        "swapped_record": {
            "sample_kgpID": SWAP_AFFECTED_SAMPLE,
            "sra_accession": (
                str(na11919_rows.iloc[0]["SRA_accession"]) if not na11919_rows.empty else None
            ),
            "internal_libraryID": (
                str(na11919_rows.iloc[0]["internal_libraryID"]) if not na11919_rows.empty else None
            ),
            "batch": int(na11919_rows.iloc[0]["batch"]) if not na11919_rows.empty else None,
            "annotated_population": (
                str(na11919_rows.iloc[0]["population"]) if not na11919_rows.empty else None
            ),
            "annotated_sex": str(na11919_rows.iloc[0]["sex"]) if not na11919_rows.empty else None,
        },
        "true_donor_record": {
            "sample_kgpID": SWAP_REPLACEMENT_SAMPLE,
            "retained_sra_accession": (
                str(hg00237_rows.iloc[0]["SRA_accession"]) if not hg00237_rows.empty else None
            ),
            "retained_libraryID": (
                str(hg00237_rows.iloc[0]["internal_libraryID"]) if not hg00237_rows.empty else None
            ),
            "retained_batch": (
                int(hg00237_rows.iloc[0]["batch"]) if not hg00237_rows.empty else None
            ),
            "dropped_batch3_sra_present_in_eQTL_inputs": bool(
                not seq_dropped.empty and (SWAP_DROPPED_SRA in set(sm["SRA_accession"].astype(str)))
            ),
        },
        "raw_continental_counts": {k: int(v) for k, v in raw_pop_counts.items()},
        "clean_continental_counts": {k: int(v) for k, v in clean_pop_counts.items()},
        "clean_subpopulation_counts": {k: int(v) for k, v in clean_subpop_counts.items()},
        "clean_sample_ids": clean_samples,
    }


def load_clean_metadata(raw_dir: Path, exclude_swapped: bool = True) -> pd.DataFrame:
    """Return sample metadata table indexed by ``sample_kgpID`` after swap remediation."""
    audit = audit_and_remediate_sample_swap(raw_dir, exclude_swapped=exclude_swapped)
    sm = pd.read_csv(
        raw_dir / "sample.metadata.MAGE.v1.0.txt",
        sep="\t",
        dtype={"sample_kgpID": str},
    )
    keep = set(audit["clean_sample_ids"])
    sm = sm[sm["sample_kgpID"].isin(keep)].copy()
    return sm.set_index("sample_kgpID").loc[audit["clean_sample_ids"]]


def load_covariate_matrix(
    raw_dir: Path,
    sample_ids: Sequence[str],
    covariate_mode: str = "all",
) -> tuple[np.ndarray, list[str]]:
    """Load encoded MAGE covariates (``sex``, ``PC1..PC5``, ``PEER1..PEER60``) in sample order."""
    cov = pd.read_csv(Path(raw_dir) / "eQTL_covariates.tab.gz", sep="\t", index_col=0)
    cov_s = cov[list(sample_ids)].T.copy()
    if "sex" in cov_s.columns:
        cov_s["sex"] = (cov_s["sex"].astype(str) == "XY").astype(np.float64)
    if covariate_mode == "all":
        cols = list(cov_s.columns)
    elif covariate_mode == "pcs_sex":
        cols = [c for c in cov_s.columns if c == "sex" or c.startswith("PC")]
    elif covariate_mode == "peers_sex":
        cols = [c for c in cov_s.columns if c == "sex" or c.startswith("PEER")]
    else:
        raise ValueError(f"Unknown covariate_mode {covariate_mode!r}")
    return cov_s[cols].to_numpy(dtype=np.float64), cols


def residualize_covariates(
    values: np.ndarray,
    Z: np.ndarray,
    include_intercept: bool = True,
) -> np.ndarray:
    """Project out covariates ``Z`` via orthonormal QR decomposition."""
    arr = np.asarray(values, dtype=np.float64)
    Z_mat = np.asarray(Z, dtype=np.float64)
    if Z_mat.ndim != 2 or Z_mat.shape[0] != arr.shape[0]:
        raise ValueError(f"Incompatible shapes values={arr.shape}, Z={Z_mat.shape}")
    if include_intercept:
        Z_mat = np.column_stack([np.ones(arr.shape[0], dtype=np.float64), Z_mat])
    Q, _ = np.linalg.qr(Z_mat, mode="reduced")
    return arr - Q @ (Q.T @ arr)


def residualize_design(
    X: np.ndarray,
    y: np.ndarray,
    Z: np.ndarray,
    include_intercept: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Frisch-Waugh-Lovell orthogonalization of genotypes ``X`` and phenotype ``y``."""
    Z_mat = np.asarray(Z, dtype=np.float64)
    if include_intercept:
        Z_mat = np.column_stack([np.ones(Z_mat.shape[0], dtype=np.float64), Z_mat])
    Q, _ = np.linalg.qr(Z_mat, mode="reduced")
    X64 = np.asarray(X, dtype=np.float64)
    y64 = np.asarray(y, dtype=np.float64)
    return X64 - Q @ (Q.T @ X64), y64 - Q @ (Q.T @ y64)


def prepare_mage_cache(
    raw_dir: Path,
    cache_dir: Path,
    exclude_swapped: bool = True,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Stream the 49 MB MAGE BED once into fast autosomal arrays and annotation tables."""
    raw_dir = Path(raw_dir)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    tag = "clean" if exclude_swapped else "raw_with_swap"
    out_dir = cache_dir / tag
    out_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "audit": out_dir / "sample_swap_audit.json",
        "samples": out_dir / "samples.tsv",
        "genes": out_dir / "genes.tsv",
        "expr": out_dir / "expression_autosomal.npy",
        "cov": out_dir / "covariates.tsv",
    }
    if not overwrite and all(p.is_file() and p.stat().st_size > 0 for p in paths.values()):
        return paths

    audit = audit_and_remediate_sample_swap(raw_dir, exclude_swapped=exclude_swapped)
    paths["audit"].write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")

    samples: list[str] = audit["clean_sample_ids"]
    meta_df = load_clean_metadata(raw_dir, exclude_swapped=exclude_swapped).reset_index()
    meta_df.to_csv(paths["samples"], sep="\t", index=False)

    Z, cov_cols = load_covariate_matrix(raw_dir, samples, covariate_mode="all")
    cov_out = pd.DataFrame(Z, columns=cov_cols)
    cov_out.insert(0, "sample_id", samples)
    cov_out.to_csv(paths["cov"], sep="\t", index=False)

    # Load gene summary annotations for autosomal metadata enrichment
    eq_sum = pd.read_csv(raw_dir / "eQTL_summary.MAGE.v1.0.txt.gz", sep="\t").set_index("ensemblID")

    autosomal_chrs = {
        **{f"chr{i}": str(i) for i in range(1, 23)},
        **{str(i): str(i) for i in range(1, 23)},
    }
    gene_rows: list[dict] = []
    expr_rows: list[np.ndarray] = []

    with gzip.open(raw_dir / "inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col_idx = [header.index(sid) for sid in samples]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            chrom_norm = autosomal_chrs.get(parts[0])
            if chrom_norm is None:
                continue
            eid = parts[3]
            s_row = eq_sum.loc[eid] if eid in eq_sum.index else None
            gene_rows.append(
                {
                    "row_index": len(gene_rows),
                    "chrom": int(chrom_norm),
                    "tss": int(parts[1]),
                    "end": int(parts[2]),
                    "ensembl_id": eid,
                    "gene_id_base": eid.split(".")[0],
                    "gene_symbol": str(s_row["geneSymbol"]) if s_row is not None else eid,
                    "variants_tested": int(s_row["variantsTested"]) if s_row is not None else 0,
                    "significant_assocs": (
                        int(s_row["significantAssocs"]) if s_row is not None else 0
                    ),
                    "is_egene": bool(s_row["is_eGene"]) if s_row is not None else False,
                    "finemapped_cs": (
                        int(s_row["finemappedCredibleSets"])
                        if s_row is not None and pd.notna(s_row["finemappedCredibleSets"])
                        else 0
                    ),
                }
            )
            expr_rows.append(np.fromiter((float(parts[i]) for i in col_idx), dtype=np.float32))

    genes_df = pd.DataFrame(gene_rows)
    expr_mat = np.vstack(expr_rows).astype(np.float32)

    tmp_expr = paths["expr"].with_suffix(".npy.tmp")
    with tmp_expr.open("wb") as out_f:
        np.save(out_f, expr_mat)
    tmp_expr.replace(paths["expr"])
    genes_df.to_csv(paths["genes"], sep="\t", index=False)
    return paths


def load_cached_gene_expression(
    cache_dir: Path,
    ensembl_ids: Sequence[str],
    exclude_swapped: bool = True,
) -> tuple[pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    """Load cached autosomal gene rows, expression array, sample order, and covariates."""
    tag = "clean" if exclude_swapped else "raw_with_swap"
    root = Path(cache_dir) / tag
    genes = pd.read_csv(root / "genes.tsv", sep="\t")
    cov_df = pd.read_csv(root / "covariates.tsv", sep="\t")
    samples = cov_df["sample_id"].astype(str).tolist()
    Z = cov_df.drop(columns=["sample_id"]).to_numpy(dtype=np.float64)

    expr_mmap = np.load(root / "expression_autosomal.npy", mmap_mode="r")
    by_eid = genes.set_index("ensembl_id")
    wanted = genes[genes["ensembl_id"].isin(ensembl_ids)].copy()
    if len(wanted) != len(set(ensembl_ids)):
        missing = sorted(set(ensembl_ids) - set(wanted["ensembl_id"]))
        raise KeyError(f"Missing Ensembl IDs in MAGE cache: {missing}")
    ordered_idx = [int(by_eid.loc[eid, "row_index"]) for eid in ensembl_ids]
    expr_sub = np.asarray(expr_mmap[ordered_idx], dtype=np.float64)
    meta_sub = by_eid.loc[list(ensembl_ids)].reset_index()
    return meta_sub, expr_sub, samples, Z
