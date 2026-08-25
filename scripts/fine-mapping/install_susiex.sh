#!/usr/bin/env bash
# One-shot installer for SuSiEx (cross-ancestry fine-mapping).
# Clones the repo into <repo>/vendor/SuSiEx and drops the prebuilt STATIC
# binary into <repo>/vendor/bin/, which scripts/_env.sh prepends to PATH.
#
# SuSiEx is NOT a Python library — it is a C++ CLI we shell out to (like plink).
# The upstream repo ships a fully-static linux_x86_64 binary in bin_static/,
# so no compiler/OpenMP is needed. If that binary ever fails to run on a host,
# fall back to `cd vendor/SuSiEx/src && make all` (needs GCC + OpenMP), which
# writes vendor/SuSiEx/bin/SuSiEx.
#
# Idempotent: skips if a working vendor/bin/SuSiEx already exists.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
VENDOR="${REPO_ROOT}/vendor"
VENDOR_BIN="${VENDOR}/bin"
SUSIEX_DIR="${VENDOR}/SuSiEx"
SUSIEX_URL="${SUSIEX_URL:-https://github.com/getian107/SuSiEx.git}"
mkdir -p "${VENDOR_BIN}"

log() { printf '[install_susiex] %s\n' "$*" >&2; }

target="${VENDOR_BIN}/SuSiEx"
if [[ -x "${target}" ]] && "${target}" --help >/dev/null 2>&1; then
  log "already installed: ${target}"
  exit 0
fi

if [[ ! -d "${SUSIEX_DIR}/.git" ]]; then
  log "cloning ${SUSIEX_URL} -> ${SUSIEX_DIR}"
  git clone --depth 1 "${SUSIEX_URL}" "${SUSIEX_DIR}"
else
  log "repo present at ${SUSIEX_DIR}"
fi

static="${SUSIEX_DIR}/bin_static/SuSiEx"
built="${SUSIEX_DIR}/bin/SuSiEx"
if [[ -f "${static}" ]]; then
  install -m 0755 "${static}" "${target}"
  log "installed static binary -> ${target}"
elif [[ -x "${built}" ]]; then
  install -m 0755 "${built}" "${target}"
  log "installed compiled binary -> ${target}"
else
  log "no prebuilt binary found; attempting 'make all' (needs GCC + OpenMP)"
  ( cd "${SUSIEX_DIR}/src" && make all )
  install -m 0755 "${built}" "${target}"
  log "compiled + installed -> ${target}"
fi

if "${target}" --help >/dev/null 2>&1; then
  log "OK: $("${target}" --help 2>&1 | head -n1)"
else
  log "ERROR: ${target} did not run --help successfully"; exit 1
fi
