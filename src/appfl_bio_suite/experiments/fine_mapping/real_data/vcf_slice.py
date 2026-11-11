"""Pure-Python Tabix (.tbi) + HTTP Range BGZF slicer for NYGC 1000G GRCh38 VCFs.

Streams only the compressed genomic windows needed for selected cis-eQTL loci
over HTTPS without downloading multi-gigabyte whole-chromosome phased VCFs,
caches index/header/locus slices locally, and converts each locus to a clean
PLINK1 binary fileset (``{bed,bim,fam}``) ready for FedFM site bundles.
"""

from __future__ import annotations

import gzip
import io
import os
import struct
import subprocess
import time
import urllib.error
import urllib.request
import zlib
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

NYGC_BASE_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20201028_3202_phased"
)
BGZF_MAX_BLOCK_BYTES = 65536
_TABIX_MIN_SHIFT = 14
_TABIX_LINEAR_STEP = 1 << _TABIX_MIN_SHIFT
_USER_AGENT = "appfl-bio-suite-mage-fedfm/1.0"


def normalize_chrom(chrom: int | str) -> str:
    """Return canonical autosomal chromosome label ``'1'`` .. ``'22'``."""
    raw = str(chrom).strip()
    if raw.lower().startswith("chr"):
        raw = raw[3:]
    if not raw.isdigit() or not 1 <= int(raw) <= 22:
        raise ValueError(f"Expected autosomal chromosome 1..22, got {chrom!r}")
    return str(int(raw))


def nygc_vcf_url(chrom: int | str, base_url: str = NYGC_BASE_URL) -> str:
    """HTTPS URL for the NYGC 30x GRCh38 phased VCF on ``chrom``."""
    c = normalize_chrom(chrom)
    return (
        f"{base_url.rstrip('/')}/"
        f"CCDG_14151_B01_GRM_WGS_2020-08-05_chr{c}.filtered.shapeit2-duohmm-phased.vcf.gz"
    )


def nygc_tbi_url(chrom: int | str, base_url: str = NYGC_BASE_URL) -> str:
    """HTTPS URL for the matching Tabix index ``.tbi`` file."""
    return f"{nygc_vcf_url(chrom, base_url=base_url)}.tbi"


def reg2bins(beg: int, end: int) -> list[int]:
    """Compute Tabix bin IDs overlapping 0-based half-open ``[beg, end)``."""
    if beg < 0 or end <= beg:
        raise ValueError(f"Invalid 0-based half-open interval [{beg}, {end})")
    end_adj = end - 1
    bins = [0]
    for offset, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        lo = offset + (beg >> shift)
        hi = offset + (end_adj >> shift)
        bins.extend(range(lo, hi + 1))
    return bins


@dataclass(frozen=True)
class TabixIndex:
    """Parsed contents of a SAMtools/HTSlib Tabix ``.tbi`` index."""

    names: tuple[str, ...]
    bins_by_ref: tuple[dict[int, tuple[tuple[int, int], ...]], ...]
    linear_by_ref: tuple[tuple[int, ...], ...]

    def resolve_ref_index(self, chrom: int | str) -> int:
        """Locate the reference index matching ``chrN`` or ``N``."""
        c = normalize_chrom(chrom)
        for candidate in (f"chr{c}", c):
            if candidate in self.names:
                return self.names.index(candidate)
        raise KeyError(f"Chromosome {chrom!r} not present in Tabix index: {self.names[:5]}...")


def _unpack_from(buf: bytes, pos: int, fmt: str) -> tuple[tuple, int]:
    size = struct.calcsize(fmt)
    if pos + size > len(buf):
        raise ValueError("Truncated Tabix index stream")
    return struct.unpack_from(fmt, buf, pos), pos + size


def parse_tbi(raw_bytes: bytes) -> TabixIndex:
    """Parse a BGZF-compressed or already-decompressed Tabix ``.tbi`` payload."""
    data = gzip.decompress(raw_bytes) if raw_bytes[:2] == b"\x1f\x8b" else raw_bytes
    if len(data) < 36 or data[:4] != b"TBI\x01":
        raise ValueError(f"Not a valid Tabix index (header={data[:4]!r})")
    (n_ref, _fmt, _col_seq, _col_beg, _col_end, _meta, _skip, l_nm), pos = _unpack_from(
        data, 4, "<8i"
    )
    if n_ref < 0 or l_nm < 0 or pos + l_nm > len(data):
        raise ValueError("Corrupt Tabix header dimensions")
    names_raw = data[pos : pos + l_nm]
    pos += l_nm
    names = tuple(part.decode("ascii") for part in names_raw.split(b"\x00")[:-1])
    if len(names) != n_ref:
        names = tuple(
            part.decode("ascii") for part in names_raw.rstrip(b"\x00").split(b"\x00") if part
        )
        if len(names) != n_ref:
            raise ValueError(f"Expected {n_ref} reference names, parsed {len(names)}")

    bins_list: list[dict[int, tuple[tuple[int, int], ...]]] = []
    linear_list: list[tuple[int, ...]] = []
    for _ in range(n_ref):
        (n_bin,), pos = _unpack_from(data, pos, "<i")
        ref_bins: dict[int, tuple[tuple[int, int], ...]] = {}
        for _ in range(n_bin):
            (bin_id, n_chunk), pos = _unpack_from(data, pos, "<Ii")
            if n_chunk > 0:
                chunks, pos = _unpack_from(data, pos, f"<{2 * n_chunk}Q")
                ref_bins[bin_id] = tuple((chunks[2 * i], chunks[2 * i + 1]) for i in range(n_chunk))
            else:
                ref_bins[bin_id] = ()
        bins_list.append(ref_bins)

        (n_intv,), pos = _unpack_from(data, pos, "<i")
        if n_intv > 0:
            lin, pos = _unpack_from(data, pos, f"<{n_intv}Q")
            linear_list.append(tuple(lin))
        else:
            linear_list.append(())

    return TabixIndex(
        names=names,
        bins_by_ref=tuple(bins_list),
        linear_by_ref=tuple(linear_list),
    )


def _merge_ranges(ranges: Sequence[tuple[int, int]], merge_gap: int) -> list[tuple[int, int]]:
    if not ranges:
        return []
    ordered = sorted(ranges)
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for lo, hi in ordered[1:]:
        if lo <= merged[-1][1] + merge_gap + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [(lo, hi) for lo, hi in merged]


def query_byte_ranges(
    tbi: TabixIndex,
    chrom: int | str,
    start_bp: int,
    end_bp: int,
    file_size: int | None = None,
    merge_gap: int = 262_144,
) -> list[tuple[int, int]]:
    """Compute inclusive compressed byte ranges ``[(lo, hi), ...]`` for ``[start_bp, end_bp]``.

    Adds ``BGZF_MAX_BLOCK_BYTES`` after each ``chunk_end >> 16`` block offset so
    that the terminal BGZF block is always fetched intact.
    """
    ref_idx = tbi.resolve_ref_index(chrom)
    beg0 = max(0, int(start_bp) - 1)
    end0 = max(beg0 + 1, int(end_bp))
    ref_bins = tbi.bins_by_ref[ref_idx]
    linear = tbi.linear_by_ref[ref_idx]
    lin_slot = beg0 // _TABIX_LINEAR_STEP
    min_off = linear[lin_slot] if lin_slot < len(linear) else 0

    raw_ranges: list[tuple[int, int]] = []
    for bid in reg2bins(beg0, end0):
        for chunk_beg, chunk_end in ref_bins.get(bid, ()):
            if min_off > 0 and chunk_end <= min_off:
                continue
            lo = max(chunk_beg >> 16, min_off >> 16 if min_off > 0 else 0)
            hi = (chunk_end >> 16) + BGZF_MAX_BLOCK_BYTES
            if file_size is not None:
                hi = min(hi, file_size - 1)
            if hi >= lo:
                raw_ranges.append((lo, hi))
    return _merge_ranges(raw_ranges, merge_gap=merge_gap)


def iter_bgzf_blocks(buf: bytes) -> Iterator[tuple[int, bytes]]:
    """Yield ``(byte_offset_in_buf, uncompressed_payload)`` for complete BGZF blocks."""
    pos = 0
    total = len(buf)
    while pos + 18 <= total:
        if buf[pos : pos + 4] != b"\x1f\x8b\x08\x04":
            break
        xlen = int.from_bytes(buf[pos + 10 : pos + 12], "little")
        ext_end = pos + 12 + xlen
        if ext_end > total:
            break
        bsize = None
        sub = pos + 12
        while sub + 4 <= ext_end:
            si1, si2 = buf[sub], buf[sub + 1]
            slen = int.from_bytes(buf[sub + 2 : sub + 4], "little")
            if si1 == ord("B") and si2 == ord("C") and slen == 2 and sub + 6 <= ext_end:
                bsize = int.from_bytes(buf[sub + 4 : sub + 6], "little")
                break
            sub += 4 + slen
        if bsize is None:
            break
        block_end = pos + bsize + 1
        if block_end > total:
            break
        if bsize == 27 and block_end == pos + 28:
            # Standard empty BGZF EOF marker
            pos = block_end
            continue
        cdata_start = ext_end
        cdata_end = block_end - 8
        if cdata_end < cdata_start:
            break
        payload = zlib.decompress(buf[cdata_start:cdata_end], -zlib.MAX_WBITS)
        yield pos, payload
        pos = block_end


def http_content_length(url: str, timeout: float = 30.0, retries: int = 4) -> int:
    """Retrieve remote file length via HTTPS ``HEAD`` with retry."""
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": _USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                cl = resp.headers.get("Content-Length")
                if cl is None:
                    raise RuntimeError(f"Missing Content-Length header from {url}")
                return int(cl)
        except (urllib.error.URLError, TimeoutError, RuntimeError) as exc:
            last_exc = exc
            time.sleep(min(2.0**attempt, 8.0))
    raise RuntimeError(f"HEAD failed for {url}: {last_exc}") from last_exc


def http_get_bytes(
    url: str,
    start: int | None = None,
    end: int | None = None,
    timeout: float = 60.0,
    retries: int = 4,
) -> bytes:
    """Download a full URL or an inclusive HTTP ``Range: bytes=start-end`` slice."""
    headers = {"User-Agent": _USER_AGENT}
    if start is not None and end is not None:
        headers["Range"] = f"bytes={start}-{end}"
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, zlib.error) as exc:
            last_exc = exc
            time.sleep(min(2.0**attempt, 8.0))
    raise RuntimeError(f"GET failed for {url} range=({start},{end}): {last_exc}") from last_exc


def ensure_tbi(chrom: int | str, cache_dir: Path, base_url: str = NYGC_BASE_URL) -> Path:
    """Download and cache the chromosome ``.tbi`` index if missing."""
    c = normalize_chrom(chrom)
    tbi_dir = Path(cache_dir) / "tbi"
    tbi_dir.mkdir(parents=True, exist_ok=True)
    target = tbi_dir / f"chr{c}.vcf.gz.tbi"
    if not target.is_file() or target.stat().st_size == 0:
        tmp = target.with_suffix(".tbi.tmp")
        tmp.write_bytes(http_get_bytes(nygc_tbi_url(c, base_url=base_url)))
        tmp.replace(target)
    return target


def ensure_vcf_header(chrom: int | str, cache_dir: Path, base_url: str = NYGC_BASE_URL) -> Path:
    """Fetch and cache the VCF header lines (through ``#CHROM``) from compressed stream start."""
    c = normalize_chrom(chrom)
    hdr_dir = Path(cache_dir) / "headers"
    hdr_dir.mkdir(parents=True, exist_ok=True)
    target = hdr_dir / f"chr{c}.header.txt"
    if target.is_file() and target.stat().st_size > 0:
        return target

    url = nygc_vcf_url(c, base_url=base_url)
    lines: list[str] = []
    offset = 0
    chunk = 262_144
    carry = b""
    done = False
    while not done and offset < 4_000_000:
        raw = carry + http_get_bytes(url, offset, offset + chunk - 1)
        consumed_end = 0
        for rel_off, payload in iter_bgzf_blocks(raw):
            xlen = int.from_bytes(raw[rel_off + 10 : rel_off + 12], "little")
            sub = rel_off + 12
            while sub + 4 <= rel_off + 12 + xlen:
                if raw[sub : sub + 2] == b"BC":
                    bsize = int.from_bytes(raw[sub + 4 : sub + 6], "little")
                    consumed_end = rel_off + bsize + 1
                    break
                sub += 4 + int.from_bytes(raw[sub + 2 : sub + 4], "little")
            for line in payload.decode("utf-8", errors="replace").splitlines(keepends=True):
                if line.startswith("#"):
                    lines.append(line)
                    if line.startswith("#CHROM\t"):
                        done = True
                        break
                else:
                    done = True
                    break
            if done:
                break
        carry = raw[consumed_end:]
        offset += chunk

    if not any(line.startswith("#CHROM\t") for line in lines):
        raise RuntimeError(f"Could not locate #CHROM header in {url}")
    tmp = target.with_suffix(".txt.tmp")
    tmp.write_text("".join(lines), encoding="utf-8")
    tmp.replace(target)
    return target


def extract_vcf_window(
    chrom: int | str,
    start_bp: int,
    end_bp: int,
    out_vcf_gz: Path,
    cache_dir: Path,
    base_url: str = NYGC_BASE_URL,
) -> Path:
    """Slice ``[start_bp, end_bp]`` from NYGC 1000G phased VCF into ``out_vcf_gz``."""
    out_vcf_gz = Path(out_vcf_gz)
    if out_vcf_gz.is_file() and out_vcf_gz.stat().st_size > 0:
        return out_vcf_gz
    out_vcf_gz.parent.mkdir(parents=True, exist_ok=True)

    c = normalize_chrom(chrom)
    tbi = parse_tbi(ensure_tbi(c, cache_dir, base_url=base_url).read_bytes())
    header_path = ensure_vcf_header(c, cache_dir, base_url=base_url)
    url = nygc_vcf_url(c, base_url=base_url)
    total_size = http_content_length(url)
    ranges = query_byte_ranges(tbi, c, start_bp, end_bp, file_size=total_size)

    chrom_ok = {f"chr{c}", c}
    seen_keys: set[str] = set()
    tmp = out_vcf_gz.with_suffix(".gz.tmp")
    with gzip.open(tmp, "wt", encoding="utf-8") as out_fh:
        out_fh.write(header_path.read_text(encoding="utf-8"))
        for lo, hi in ranges:
            raw = http_get_bytes(url, lo, hi)
            # Align to first BGZF member in case `lo` landed on block start cleanly
            start_pos = raw.find(b"\x1f\x8b\x08\x04")
            if start_pos < 0:
                continue
            buf = raw[start_pos:]
            text_stream = io.StringIO()
            for _off, payload in iter_bgzf_blocks(buf):
                text_stream.write(payload.decode("utf-8"))
            for line in text_stream.getvalue().splitlines(keepends=True):
                if not line or line.startswith("#"):
                    continue
                tab1 = line.find("\t")
                tab2 = line.find("\t", tab1 + 1)
                tab3 = line.find("\t", tab2 + 1)
                if tab3 < 0:
                    continue
                rec_chr = line[:tab1]
                if rec_chr not in chrom_ok:
                    continue
                pos = int(line[tab1 + 1 : tab2])
                if pos < start_bp:
                    continue
                if pos > end_bp:
                    continue
                vid = line[tab2 + 1 : tab3]
                key = f"{pos}:{vid}"
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                out_fh.write(line)
    tmp.replace(out_vcf_gz)
    return out_vcf_gz


def extract_locus_plink(
    chrom: int | str,
    start_bp: int,
    end_bp: int,
    out_prefix: Path,
    keep_samples_file: Path,
    cache_dir: Path,
    plink2_bin: str | os.PathLike[str] = "plink2",
    maf: float = 0.01,
    snps_only: bool = True,
    rm_dup_mode: str | None = "exclude-mismatch",
    base_url: str = NYGC_BASE_URL,
    overwrite: bool = False,
) -> Path:
    """Materialize a clean biallelic PLINK1 fileset for a MAGE cis window."""
    out_prefix = Path(out_prefix)
    bed = Path(f"{out_prefix}.bed")
    bim = Path(f"{out_prefix}.bim")
    fam = Path(f"{out_prefix}.fam")
    if (
        not overwrite
        and bed.is_file()
        and bim.is_file()
        and fam.is_file()
        and bed.stat().st_size > 3
    ):
        return out_prefix

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    c = normalize_chrom(chrom)
    vcf_gz = Path(cache_dir) / "vcf" / f"chr{c}_{start_bp}_{end_bp}.vcf.gz"
    extract_vcf_window(c, start_bp, end_bp, vcf_gz, cache_dir=cache_dir, base_url=base_url)

    cmd = [
        str(plink2_bin),
        "--vcf",
        str(vcf_gz),
        "--keep",
        str(keep_samples_file),
        "--chr",
        c,
        "--from-bp",
        str(start_bp),
        "--to-bp",
        str(end_bp),
        "--max-alleles",
        "2",
    ]
    if snps_only:
        cmd.extend(["--snps-only", "just-acgt"])
    if rm_dup_mode:
        cmd.extend(["--rm-dup", rm_dup_mode])
    cmd.extend(
        [
            "--maf",
            str(maf),
            "--geno",
            "0.0",
            "--make-bed",
            "--threads",
            "1",
            "--out",
            str(out_prefix),
        ]
    )
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        raise RuntimeError(
            f"plink2 failed (rc={res.returncode}) for {out_prefix.name}:\n"
            f"{res.stderr or res.stdout}"
        )
    if not bed.is_file() or not bim.is_file() or not fam.is_file():
        raise RuntimeError(f"Expected PLINK1 output at {out_prefix}.{{bed,bim,fam}}")
    return out_prefix
