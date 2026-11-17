#!/usr/bin/env bash
# Called by --worker-command. The suite is mounted as source, never pip-installed
# into ALMA3's runtime (which would resolve incompatible Torch pins).
set -euo pipefail
: "ALMA3_SIF:?Set ALMA3_SIF to a pinned 3.1.2 Apptainer image}"
: "ALMA3_SUITE_SOURCE:?Set ALMA3_SUITE_SOURCE to this worktree's src directory}"
# Bind shared paths explicitly on Sophia when they are not automatically visible.
# The runtime overlay must provide sklearn for the suite's existing metrics.
exec apptainer exec --cleanenv \
    --env "PYTHONPATH=${ALMA3_SUITE_SOURCE}" \
    --bind "${ALMA3_SUITE_SOURCE}:${ALMA3_SUITE_SOURCE}:ro" \
    "${ALMA3_SIF}" python "$@"
