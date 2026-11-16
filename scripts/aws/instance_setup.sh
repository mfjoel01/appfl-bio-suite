#!/usr/bin/env bash
# Turn a fresh Ubuntu 24.04 EC2 instance into the live federation map, or refresh one.
#
# Run as root ON THE INSTANCE, from a checkout of this repository at the commit to serve.
# Nobody runs it by hand: scripts/aws/live_map.py puts it in the instance's first boot
# (`deploy`) and re-runs it over SSM (`update`). It is idempotent -- every step replaces
# what the previous run left.
#
# What it builds:
#
#   /opt/live-map/venv     the suite, installed against constraints.txt
#   /opt/live-map/site     scripts/build_site.sh output: the same pages GitHub Pages serves,
#                          with the map exported to poll for endpoint status
#   live-map-prober        systemd unit: `appfl-bio-suite watch live`, rewriting
#                          site/map/endpoints.status.json
#   caddy                  serves site/ on 80/443, with a certificate for $LIVE_MAP_HOST
#
# Configuration is /etc/live-map.env, written by the deployer:
#
#   LIVE_MAP_HOST       name to request a certificate for; empty serves plain HTTP only
#   LIVE_MAP_REGION     AWS region holding the parameters below
#   LIVE_MAP_PARAMS     SSM parameter prefix: <prefix>/targets, /globus-client-id,
#                       /globus-client-secret
#   LIVE_MAP_INTERVAL   seconds between probes
#
# The Globus secret and the endpoint UUIDs are never written by this script. The prober's
# launcher reads them from SSM at start, with the instance's role; see live_map_prober.py.

set -euo pipefail

[[ $EUID -eq 0 ]] || { echo "run as root" >&2; exit 1; }

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
home=/opt/live-map
# shellcheck disable=SC1091
source /etc/live-map.env

export DEBIAN_FRONTEND=noninteractive
apt-get update -q
apt-get install -y -q python3.12-venv python3-dev build-essential git caddy

id -u livemap >/dev/null 2>&1 || useradd --system --home-dir "$home" --shell /usr/sbin/nologin livemap
mkdir -p "$home"
# The Globus SDK keeps its token cache under $HOME; the rest of $home stays root's.
install -d -o livemap -g livemap -m 0700 "$home/.globus_compute"

# The pinned install every site in a federation runs, as on GitHub Pages' build. The CPU
# wheel index keeps torch -- a dependency of the package, unused here -- to ~200 MB
# instead of several GB of CUDA libraries this instance has no use for.
[[ -x $home/venv/bin/python ]] || python3.12 -m venv "$home/venv"
PIP_EXTRA_INDEX_URL=https://download.pytorch.org/whl/cpu \
  "$home/venv/bin/pip" install --no-cache-dir -q -e "$repo[watch]" boto3 -c "$repo/constraints.txt"

# Built beside the live directory and swapped in, so a failed build leaves the old site up.
rm -rf "$home/site.new"
(
  cd "$repo"
  SITE_LIVE_STATUS=1 APPFL_BIO_SUITE="$home/venv/bin/appfl-bio-suite" PYTHON="$home/venv/bin/python" \
    scripts/build_site.sh "$home/site.new"
)
chown -R livemap: "$home/site.new"
rm -rf "$home/site.old"
[[ -d $home/site ]] && mv "$home/site" "$home/site.old"
mv "$home/site.new" "$home/site"
rm -rf "$home/site.old"

cat > /etc/systemd/system/live-map-prober.service <<EOF
[Unit]
Description=APPFL Bio Suite live map: endpoint status prober
After=network-online.target
Wants=network-online.target

[Service]
User=livemap
EnvironmentFile=/etc/live-map.env
Environment=HOME=$home
ExecStart=$home/venv/bin/python $repo/scripts/aws/live_map_prober.py $home/site/map
Restart=always
RestartSec=15
RuntimeDirectory=live-map
RuntimeDirectoryMode=0700
NoNewPrivileges=true
ProtectSystem=strict
ReadWritePaths=$home

[Install]
WantedBy=multi-user.target
EOF

# The status file is the one thing that changes; everything else may be cached briefly.
# With a host name Caddy obtains and renews the certificate itself and redirects HTTP for
# that name; the bare :80 site keeps the IP address answering either way.
{
  if [[ -n ${LIVE_MAP_HOST:-} ]]; then echo "$LIVE_MAP_HOST, :80 {"; else echo ":80 {"; fi
  cat <<EOF
	root * $home/site
	encode gzip
	file_server
	@status path /map/endpoints.status.json
	header @status Cache-Control "no-store"
	header -Server
}
EOF
} > /etc/caddy/Caddyfile

systemctl daemon-reload
systemctl enable -q caddy live-map-prober
systemctl reload-or-restart caddy
systemctl restart live-map-prober

echo "live map: ${LIVE_MAP_HOST:-(no host name)}  commit $(git -C "$repo" rev-parse --short HEAD)"
