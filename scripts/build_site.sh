#!/usr/bin/env bash
# Assemble the project site published at https://mfjoel01.github.io/appfl-bio-suite/
#
#     scripts/build_site.sh [OUTDIR]          # default: _site
#
# The pages in website/ are hand-written. The network-map preview under OUTDIR/map is
# exported here, from website/demo-federation.yaml -- twenty invented institutions,
# never local/federation.yaml, which names real partners and carries their contacts.
# check_site_preview.py enforces that rather than trusting this comment.
#
# .github/workflows/pages.yml runs this script, so previewing locally exercises the same
# assembly that deploys:
#
#     pip install -e '.[watch]' -c constraints.txt
#     scripts/build_site.sh local/site-preview
#     python -m http.server -d local/site-preview 8000
#
# A plain file:// open will not work -- the viewer fetches its data from alongside
# itself, and browsers refuse that for local files.
#
# APPFL_BIO_SUITE overrides the CLI, for a checkout that is not pip-installed.

set -euo pipefail

out="${1:-_site}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cli="${APPFL_BIO_SUITE:-appfl-bio-suite}"
python="${PYTHON:-python3}"

cd "$root"

# Replaced rather than merged: a stale map.html from an older viewer is exactly the
# thing a preview is supposed to catch.
rm -rf "$out"
mkdir -p "$out"

cp -R website/. "$out"/
cp src/appfl_bio_suite/core/watch_assets/suite-logo.png "$out"/assets/suite-logo.png

"$cli" watch export \
  --federation website/demo-federation.yaml \
  --out "$out/map" \
  --title "APPFL Bio Suite federation network" >/dev/null

# Serve the tree as-is. Redundant under the Actions deployment, and the thing that keeps
# docs/ from being rewritten if Pages is ever pointed at a branch instead.
touch "$out/.nojekyll"

"$python" scripts/check_site_preview.py "$out/map/network.map.json"

echo "site: $root/$out"
