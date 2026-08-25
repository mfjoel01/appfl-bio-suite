#!/usr/bin/env bash
# Turnkey HAPNEST downloader for the FedFM data-simulation pipeline.
#
# Pulls EBI BioStudies S-BSST936 (HAPNEST pre-generated synthetic cohort).
# The published layout is already in the form we need: one combined PLINK
# binary per chromosome with all six superpopulations merged.
#
#   genotypes/synthetic_v1_chr-<N>.{bed,bim,fam}   for N=1..22
#   synthetic_v1.sample                            1,008,000 rows, superpop code per row
#
# This script:
#   1. Downloads the .sample population label file (~4 MB).
#   2. Downloads the requested per-chromosome PLINK binaries (default: all 22).
#   3. Builds population_manifest.tsv by joining .sample row-by-row with chr1.fam.
#   4. Symlinks the chrN.{bed,bim,fam} layout the pipeline expects on top of the
#      downloaded files (no copy, no extra disk).
#
# Resumable: wget --continue and aria2c --continue=true both pick up partials.
# Idempotent: re-running skips already-complete files.
#
# Total dataset size: ~1.5 TB for all 22 chromosomes (chr1.bed = 134 GB).
# Use CHROMS="1" to limit to chr1 (~135 GB) which matches the default
# `chromosome: 1` in simulation_config.yaml.
#
# Usage:
#   scripts/download_hapnest.sh [dest_dir]
#   CHROMS="1" scripts/download_hapnest.sh
#   FORCE=1 scripts/download_hapnest.sh

set -euo pipefail

# ---------------------------------------------------------------------------
ACCESSION="S-BSST936"
BASE_URL="https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/936/${ACCESSION}/Files"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${1:-${REPO_ROOT}/data/raw/hapnest}"
DEST="$(mkdir -p "${DEST}" && cd "${DEST}" && pwd)"
GENO_DIR="${DEST}/genotypes"
mkdir -p "${GENO_DIR}"

CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
FORCE="${FORCE:-0}"

LOG_DIR="${REPO_ROOT}/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/download_hapnest_$(date +%Y%m%dT%H%M%S).log"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "${LOG_FILE}" >&2; }
die() { log "FATAL: $*"; exit 1; }

# Approximate file sizes (bytes), used for disk-space pre-check.
# bed sizes scale ~linearly with chromosome length.
# Values derived from EBI Content-Length headers.
declare -A BED_SIZE_GB=(
  [1]=135  [2]=147  [3]=121  [4]=117  [5]=108  [6]=102  [7]=95
  [8]=91   [9]=72   [10]=82  [11]=82  [12]=80  [13]=58  [14]=54
  [15]=48  [16]=51  [17]=44  [18]=46  [19]=36  [20]=37  [21]=22
  [22]=29
)

log "HAPNEST download — dest=${DEST}"
log "Chromosomes requested: ${CHROMS}"
log "Log: ${LOG_FILE}"

# ---------------------------------------------------------------------------
# Disk space pre-check
# ---------------------------------------------------------------------------
required_gb=2  # baseline: .sample + manifest
for c in ${CHROMS}; do
  required_gb=$((required_gb + ${BED_SIZE_GB[${c}]:-50} + 1))
done
avail_kb=$(df -P "${DEST}" | awk 'NR==2 {print $4}')
avail_gb=$((avail_kb / 1024 / 1024))
log "Estimated download size: ${required_gb} GB. Available at ${DEST}: ${avail_gb} GB."
if [[ "${avail_gb}" -lt "${required_gb}" ]]; then
  die "Insufficient disk space: need ${required_gb} GB, have ${avail_gb} GB. \
Free space or set CHROMS to a smaller list."
fi

# ---------------------------------------------------------------------------
# Downloader selector
# ---------------------------------------------------------------------------
USE_ARIA2=0
if command -v aria2c >/dev/null 2>&1; then
  USE_ARIA2=1
  log "Using aria2c (parallel, resumable)."
else
  log "aria2c not found; using wget --continue (slower)."
  log "Hint: install aria2c for 5-10× speedup on Polaris (conda install -c conda-forge aria2)."
fi

# Single-file fetch with resume.
fetch_one() {
  local url="$1" out="$2"
  if [[ -s "${out}" && "${FORCE}" != "1" ]]; then
    # Verify size against remote Content-Length.
    local remote_size local_size
    remote_size=$(curl -sIL "${url}" | awk -F': ' 'tolower($1)=="content-length"{gsub(/\r/,"",$2); print $2}' | tail -n1)
    local_size=$(stat -c %s "${out}" 2>/dev/null || echo 0)
    if [[ -n "${remote_size}" && "${local_size}" == "${remote_size}" ]]; then
      log "  exists & size OK: $(basename "${out}") (${local_size} bytes)"
      return 0
    fi
    log "  partial: $(basename "${out}") local=${local_size} remote=${remote_size:-?} — resuming"
  fi
  if [[ "${USE_ARIA2}" == "1" ]]; then
    aria2c \
      --continue=true --auto-file-renaming=false \
      --max-connection-per-server=8 --split=8 \
      --retry-wait=30 --max-tries=20 --timeout=120 \
      --console-log-level=warn --summary-interval=30 \
      --dir="$(dirname "${out}")" --out="$(basename "${out}")" \
      "${url}" 2>&1 | tee -a "${LOG_FILE}"
  else
    wget --continue --tries=20 --waitretry=30 --timeout=120 \
         -O "${out}" "${url}" 2>&1 | tee -a "${LOG_FILE}"
  fi
}

# ---------------------------------------------------------------------------
# 1. Sample (population labels)
# ---------------------------------------------------------------------------
log "Step 1/4: synthetic_v1.sample"
fetch_one "${BASE_URL}/synthetic_v1.sample" "${DEST}/synthetic_v1.sample"
n_sample=$(wc -l < "${DEST}/synthetic_v1.sample")
log "  ${n_sample} rows in synthetic_v1.sample"
[[ "${n_sample}" -eq 1008000 ]] \
  || log "  WARNING: expected 1,008,000 rows, got ${n_sample}"

# ---------------------------------------------------------------------------
# 2. Per-chromosome PLINK trios
# ---------------------------------------------------------------------------
log "Step 2/4: per-chromosome PLINK binaries"
for c in ${CHROMS}; do
  for ext in fam bim bed; do
    url="${BASE_URL}/genotypes/synthetic_v1_chr-${c}.${ext}"
    out="${GENO_DIR}/synthetic_v1_chr-${c}.${ext}"
    log "  chr${c}.${ext}"
    fetch_one "${url}" "${out}"
  done
done

# ---------------------------------------------------------------------------
# 3. population_manifest.tsv (FID, IID, superpopulation)
# ---------------------------------------------------------------------------
MANIFEST="${DEST}/population_manifest.tsv"
log "Step 3/4: building ${MANIFEST}"
# Use chr1.fam if it was downloaded; otherwise pick any downloaded fam.
ref_fam=""
for c in ${CHROMS}; do
  if [[ -s "${GENO_DIR}/synthetic_v1_chr-${c}.fam" ]]; then
    ref_fam="${GENO_DIR}/synthetic_v1_chr-${c}.fam"
    break
  fi
done
[[ -n "${ref_fam}" ]] || die "No .fam file downloaded — cannot build manifest"

n_fam=$(wc -l < "${ref_fam}")
[[ "${n_fam}" -eq "${n_sample}" ]] \
  || die "Row-count mismatch: ${ref_fam}=${n_fam} vs synthetic_v1.sample=${n_sample}"

# Stream-join row by row: FID IID from .fam, superpop from .sample.
awk 'BEGIN{OFS="\t"; print "FID","IID","superpopulation"}
     NR==FNR { sp[NR]=$1; next }
     { print $1, $2, sp[FNR] }' \
  "${DEST}/synthetic_v1.sample" "${ref_fam}" > "${MANIFEST}"

n_manifest=$(($(wc -l < "${MANIFEST}") - 1))
log "  wrote ${MANIFEST} (${n_manifest} individuals)"

# Per-superpop sanity check.
log "  superpop counts:"
awk -F'\t' 'NR>1 {print $3}' "${MANIFEST}" | sort | uniq -c | tee -a "${LOG_FILE}"

# ---------------------------------------------------------------------------
# 4. Symlink the chr<N>.{bed,bim,fam} layout the pipeline expects.
# ---------------------------------------------------------------------------
log "Step 4/4: creating chr<N>.{bed,bim,fam} symlinks"
for c in ${CHROMS}; do
  for ext in bed bim fam; do
    src="${GENO_DIR}/synthetic_v1_chr-${c}.${ext}"
    dst="${DEST}/chr${c}.${ext}"
    if [[ ! -s "${src}" ]]; then
      log "  WARNING: ${src} missing — skipping chr${c}.${ext} symlink"
      continue
    fi
    # Relative symlink so the link still resolves if the dir is moved.
    ln -sfn "genotypes/synthetic_v1_chr-${c}.${ext}" "${dst}"
  done
done

log "Done."
log "Pipeline-ready inputs:"
ls -lh "${DEST}/chr"*.bed 2>/dev/null | tee -a "${LOG_FILE}" || true
ls -lh "${MANIFEST}" | tee -a "${LOG_FILE}"
log "Next: ./scripts/run_all.sh"
