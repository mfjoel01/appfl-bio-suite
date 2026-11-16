#!/usr/bin/env python
"""Start `appfl-bio-suite watch live` on the EC2 instance, with its secrets from SSM.

The systemd unit written by instance_setup.sh runs this. It exists so that the sensitive
inputs -- the prober's Globus service-account secret, the real endpoint UUIDs, and the
secret of any site service account the probe list routes an endpoint to -- are fetched
with the instance's IAM role each time the service starts, and are never in user-data,
in the unit file, or in the repository checkout. Rotating any of them is
`live_map.py update`, which rewrites the parameters and restarts the service.

    live_map_prober.py <exported map directory>

Reads LIVE_MAP_REGION, LIVE_MAP_PARAMS and LIVE_MAP_INTERVAL from the environment
(/etc/live-map.env).
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2

    import boto3

    prefix = os.environ["LIVE_MAP_PARAMS"].rstrip("/")
    ssm = boto3.client("ssm", region_name=os.environ["LIVE_MAP_REGION"])

    def parameter(name: str) -> str:
        response = ssm.get_parameter(Name=f"{prefix}/{name}", WithDecryption=True)
        return response["Parameter"]["Value"]

    # systemd's RuntimeDirectory: tmpfs, mode 0700, removed when the service stops.
    runtime = Path(os.environ.get("RUNTIME_DIRECTORY", "/run/live-map"))
    targets = runtime / "targets.json"
    targets.touch(mode=0o600)
    targets.write_text(parameter("targets"), encoding="utf-8")

    env = dict(os.environ)
    env["GLOBUS_COMPUTE_CLIENT_ID"] = parameter("globus-client-id")
    env["GLOBUS_COMPUTE_CLIENT_SECRET"] = parameter("globus-client-secret")
    # The secrets of any sites the probe list asks about as their own service account,
    # under the variable names its routes give. Absent when no route is used.
    try:
        env.update(json.loads(parameter("site-secrets")))
    except ssm.exceptions.ParameterNotFound:
        pass

    cli = Path(sys.executable).with_name("appfl-bio-suite")
    command = [
        str(cli), "watch", "live",
        "--targets", str(targets),
        "--out", argv[1],
        "--interval", os.environ.get("LIVE_MAP_INTERVAL", "60"),
    ]  # fmt: skip
    os.execve(command[0], command, env)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
