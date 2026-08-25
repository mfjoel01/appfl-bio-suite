#!/usr/bin/env bash
# One-shot installer for PLINK 1.9 and PLINK 2.0 static Linux binaries.
# Drops them into <repo>/vendor/bin/, which scripts/_env.sh prepends to PATH.
#
# Idempotent: skips download if a working binary already lives in vendor/bin/.
# No conda, no root, no module — just a static binary fetch.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
VENDOR_BIN="${REPO_ROOT}/vendor/bin"
mkdir -p "${VENDOR_BIN}"

# Known-good build dates verified against s3.amazonaws.com on 2026-06-03.
# Override either of these via env if you want a different build.
PLINK1_DATE="${PLINK1_DATE:-20241022}"
PLINK2_DATE="${PLINK2_DATE:-20250129}"
PLINK2_TRACK="${PLINK2_TRACK:-alpha6}"

PLINK1_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_${PLINK1_DATE}.zip"
PLINK2_URL="https://s3.amazonaws.com/plink2-assets/${PLINK2_TRACK}/plink2_linux_x86_64_${PLINK2_DATE}.zip"

log() { printf '[install_plink] %s\n' "$*" >&2; }

install_one() {
  local name="$1" url="$2" binname="$3"
  local target="${VENDOR_BIN}/${binname}"
  if [[ -x "${target}" ]] && "${target}" --version >/dev/null 2>&1; then
    log "${name}: already installed ($("${target}" --version 2>&1 | head -n1))"
    return 0
  fi
  local tmpdir
  tmpdir="$(mktemp -d)"
  log "${name}: downloading ${url}"
  if ! curl -fL --retry 5 --retry-delay 10 -o "${tmpdir}/p.zip" "${url}"; then
    rm -rf "${tmpdir}"
    log "ERROR: download failed for ${url}"
    log "If the build date is stale, find a current one at"
    log "  https://www.cog-genomics.org/plink/   (1.9)"
    log "  https://www.cog-genomics.org/plink/2.0/  (2.0)"
    log "and re-run with e.g. PLINK1_DATE=YYYYMMDD ./scripts/install_plink.sh"
    return 1
  fi
  unzip -q -o "${tmpdir}/p.zip" -d "${tmpdir}/x"
  # PLINK ships the executable named ``plink`` or ``plink2`` in the zip root.
  local src
  src="$(find "${tmpdir}/x" -maxdepth 2 -type f -name "${binname}" -executable | head -n1)"
  if [[ -z "${src}" ]]; then
    # PLINK 2's zip may not set the executable bit on all filesystems.
    src="$(find "${tmpdir}/x" -maxdepth 2 -type f -name "${binname}" | head -n1)"
  fi
  [[ -n "${src}" ]] || { log "ERROR: ${binname} not found inside zip"; rm -rf "${tmpdir}"; return 1; }
  install -m 0755 "${src}" "${target}"
  rm -rf "${tmpdir}"
  log "${name}: installed to ${target} ($(${target} --version 2>&1 | head -n1))"
}

install_one "PLINK 1.9" "${PLINK1_URL}" "plink"
install_one "PLINK 2.0" "${PLINK2_URL}" "plink2"

log "Done. Run ./scripts/run_sampling.sh — the _env.sh helper will put vendor/bin on PATH."
