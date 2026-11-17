#!/usr/bin/env python
"""Put the federation map on AWS, live, with one command -- and take it down with another.

    python scripts/aws/live_map.py deploy     # create everything, print the URL
    python scripts/aws/live_map.py status     # is it up, and what does it say
    python scripts/aws/live_map.py update     # serve a newer commit / new targets / new secret
    python scripts/aws/live_map.py destroy --yes

WHAT IT STANDS UP
-----------------
One small EC2 instance serving the project site -- the same pages GitHub Pages serves,
built by the same scripts/build_site.sh -- plus a prober that asks Globus Compute every
minute whether each partner's endpoint is connected, and a map that shows the answer.
docs/coordinator/live-map-aws.md has the design; core/live_status.py has the reasoning.

Everything is created through the API and tagged ``Project=<name>``, so nothing here
needs the AWS console, and ``destroy`` can find all of it again:

  * a CloudFront distribution in front of the instance, which is what gives the page an
    ``https://…cloudfront.net`` address with a valid certificate and no domain of ours
    (``--domain`` serves your own name from the instance instead; ``--no-cloudfront``
    serves plain HTTP on the address only). Caching is off: the status file changes
    every minute, and a cached copy is exactly the stale page this exists to avoid;
  * an Elastic IP, so the address survives the instance being replaced;
  * a security group open on 80 and 443 only -- there is no SSH and no key pair;
    ``update`` reaches the instance through Systems Manager;
  * an instance role that may read this deployment's three parameters and nothing else;
  * three SSM parameters: the probe targets (real endpoint UUIDs), and the Globus service
    account's id and secret (encrypted). They are the only sensitive inputs, and they go
    from this machine to Parameter Store to the prober's memory -- never into user-data,
    which anyone with ec2:DescribeInstanceAttribute can read.

WHAT YOU NEED ONCE
------------------
1. AWS credentials in the environment or ~/.aws (``--profile`` to pick one).
2. A Globus service account for the prober: app.globus.org -> Settings -> Developers ->
   "Register a service account". Save its id and secret as
   ``local/live-map/globus-client.json``: ``{"client_id": "...", "client_secret": "..."}``.
   Multi-user endpoints have answered status requests from identities that do not own
   them; an endpoint that refuses this one shows as "unknown" on the map, not offline.
3. The probe list, from your real federation:
   ``appfl-bio-suite watch targets --out local/live-map/targets.json``
   An endpoint owned by a site's own service account (single-user mode) answers status
   requests only from that account. Give the list a route for it with
   ``--site-credentials``, and export that client's secret before deploying -- a separate
   secret minted for the prober, since that account can submit tasks, not only ask.

The commit being served must be on GitHub: the instance clones the public repository.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

REPO_URL = "https://github.com/mfjoel01/appfl-bio-suite.git"
ROOT = Path(__file__).resolve().parents[2]

DEFAULT_NAME = "appfl-bio-suite-live-map"
DEFAULT_REGION = "us-east-1"
# x86, 2 GB. x86 because the instance installs the full pinned suite, and on x86 that
# install is proven on every push by the Pages build; on arm64 several pinned packages
# (psutil<6, among others) have no wheel and would compile on the instance. 2 GB plus
# the swap file instance_setup.sh adds is enough for the install, which is the peak.
DEFAULT_INSTANCE_TYPE = "t3.small"
DEFAULT_TARGETS = ROOT / "local" / "live-map" / "targets.json"
DEFAULT_GLOBUS = ROOT / "local" / "live-map" / "globus-client.json"

_NAME = re.compile(r"^[a-z][a-z0-9-]{2,40}$")
_HOST = re.compile(r"^(?=.{4,253}$)([a-z0-9]([a-z0-9-]{0,61}[a-z0-9])?\.)+[a-z]{2,}$")
_SHA = re.compile(r"^[0-9a-f]{40}$")
_UUID = re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$")


class DeployError(RuntimeError):
    """Something the operator has to fix. The message says what."""


def say(message: str) -> None:
    print(message, flush=True)


# ---------------------------------------------------------------------------
# inputs: pure, and tested without AWS
# ---------------------------------------------------------------------------


def resolve_ref(ref: str | None, root: Path = ROOT) -> str:
    """The full commit to serve, checked to be one the instance can fetch from GitHub."""

    def git(*args: str) -> str:
        return subprocess.run(
            ["git", "-C", str(root), *args], check=True, capture_output=True, text=True
        ).stdout.strip()

    try:
        sha = git("rev-parse", "--verify", f"{ref or 'HEAD'}^{{commit}}")
        published = git("branch", "-r", "--contains", sha)
    except subprocess.CalledProcessError as exc:
        raise DeployError(f"cannot resolve commit {ref or 'HEAD'}: {exc.stderr.strip()}") from exc
    if not published:
        raise DeployError(
            f"commit {sha[:10]} is not on GitHub yet. The instance clones the public "
            "repository, so push the branch first (git push -u origin HEAD)."
        )
    return sha


def public_site_ids(root: Path = ROOT) -> set[str]:
    """Every site id the published map draws."""
    import yaml

    partners = json.loads((root / "website" / "partners.json").read_text(encoding="utf-8"))
    federation = yaml.safe_load((root / "website" / "federation.yaml").read_text(encoding="utf-8"))
    return {p["id"] for p in partners.get("partners", [])} | {
        s["id"] for s in federation.get("sites", [])
    }


def read_targets(path: Path, known_ids: set[str] | None = None) -> tuple[str, list[str], list[str]]:
    """The targets file as compact JSON to store, warnings, and the secrets it names.

    The last is the environment variables its credential routes read, sorted. The
    values are read from this machine's environment at deploy time and stored
    encrypted beside the targets; the file itself never holds one.
    """
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise DeployError(
            f"no probe targets at {path}. Write them from your real federation:\n"
            f"    appfl-bio-suite watch targets --out {path}"
        ) from None
    except json.JSONDecodeError as exc:
        raise DeployError(f"{path} is not JSON: {exc}") from exc

    sites = raw.get("sites") if isinstance(raw, dict) else None
    if not isinstance(sites, dict) or not sites or set(raw) - {"sites", "routes"}:
        raise DeployError(f"{path}: expected a non-empty `sites` object, and optionally `routes`")
    for site, uuids in sites.items():
        if not isinstance(uuids, list) or not uuids:
            raise DeployError(f"{path}: site {site!r} needs a list of endpoint UUIDs")
        for uuid in uuids:
            if not isinstance(uuid, str) or not _UUID.match(uuid.lower()):
                raise DeployError(f"{path}: site {site!r}: not an endpoint UUID: {uuid!r}")

    probed = {uuid.lower() for uuids in sites.values() for uuid in uuids}
    routes = raw.get("routes") or {}
    if not isinstance(routes, dict):
        raise DeployError(f"{path}: `routes` must map endpoint UUIDs to routes")
    for uuid, route in routes.items():
        if (
            uuid.lower() not in probed
            or not isinstance(route, dict)
            or set(route) != {"client_id", "client_secret_env"}
            or not _UUID.match(str(route["client_id"]).lower())
            or not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", str(route["client_secret_env"]))
        ):
            raise DeployError(
                f"{path}: route for {uuid} needs exactly client_id and client_secret_env, "
                "for an endpoint listed under `sites` -- and never the secret itself"
            )
    secret_envs = sorted({route["client_secret_env"] for route in routes.values()})

    warnings = []
    if known_ids is not None:
        for site in sorted(set(sites) - known_ids):
            warnings.append(
                f"site {site!r} is not on the published map (website/partners.json), so its "
                "status will never be shown. Rename it with `watch targets --rename`."
            )
    document = {"sites": sites, **({"routes": routes} if routes else {})}
    body = json.dumps(document, separators=(",", ":"))
    if len(body) > 4096:
        raise DeployError(f"{path}: over the 4 KB a standard SSM parameter holds")
    return body, warnings, secret_envs


def read_site_secrets(names: list[str], environ: dict | None = None) -> dict[str, str]:
    """The secrets a probe list's routes name, from this machine's environment."""
    environ = os.environ if environ is None else environ
    missing = [name for name in names if not environ.get(name)]
    if missing:
        raise DeployError(
            "the probe list asks about some endpoints as a site's own service account. "
            f"Export the secret(s) first, in this shell only: {', '.join(missing)}\n"
            "Generate a separate secret on that client for the prober, so it can be revoked "
            "without touching the coordinator's."
        )
    return {name: environ[name] for name in names}


def read_globus_client(path: Path, environ: dict | None = None) -> tuple[str, str] | None:
    """The prober's Globus service account, from the environment or a local file."""
    environ = os.environ if environ is None else environ
    client_id = environ.get("LIVE_MAP_GLOBUS_CLIENT_ID")
    secret = environ.get("LIVE_MAP_GLOBUS_CLIENT_SECRET")
    if not (client_id and secret):
        if not path.is_file():
            return None
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
            client_id = raw["client_id"]
            secret = raw.get("client_secret") or raw["secret"]
        except (json.JSONDecodeError, KeyError, TypeError) as exc:
            raise DeployError(
                f'{path}: expected {{"client_id": "...", "client_secret": "..."}}'
            ) from exc
    if not _UUID.match(str(client_id).lower()) or not secret:
        raise DeployError("the Globus client id must be a UUID, with a non-empty secret")
    return str(client_id), str(secret)


def ec2_public_dns(ip: str, region: str) -> str:
    """EC2's public name for an address: the origin CloudFront fetches from.

    CloudFront needs a name, not an address. Not sslip.io and its kind: Argonne's resolvers
    answer those with a sinkhole, so a page on one does not open from the lab network.
    """
    dashed = ip.replace(".", "-")
    if region == "us-east-1":
        return f"ec2-{dashed}.compute-1.amazonaws.com"
    return f"ec2-{dashed}.{region}.compute.amazonaws.com"


def distribution_config(
    *, name: str, origin: str, cache_policy: str, request_policy: str, reference: str
) -> dict:
    """CloudFront in front of the instance: HTTPS to viewers, plain HTTP to the origin.

    The origin leg carries nothing but the public page. Caching is disabled and every
    viewer header but Host is forwarded, so the instance compresses as usual and answers
    to its own name.
    """
    methods = {"Quantity": 2, "Items": ["GET", "HEAD"]}
    return {
        "CallerReference": reference,
        "Comment": name,
        "Enabled": True,
        "Origins": {
            "Quantity": 1,
            "Items": [
                {
                    "Id": "instance",
                    "DomainName": origin,
                    "CustomOriginConfig": {
                        "HTTPPort": 80,
                        "HTTPSPort": 443,
                        "OriginProtocolPolicy": "http-only",
                    },
                }
            ],
        },
        "DefaultCacheBehavior": {
            "TargetOriginId": "instance",
            "ViewerProtocolPolicy": "redirect-to-https",
            "CachePolicyId": cache_policy,
            "OriginRequestPolicyId": request_policy,
            "AllowedMethods": {**methods, "CachedMethods": methods},
            "Compress": False,
        },
        "PriceClass": "PriceClass_All",
        "HttpVersion": "http2and3",
        "IsIPV6Enabled": True,
        "ViewerCertificate": {"CloudFrontDefaultCertificate": True},
    }


def render_user_data(*, sha: str, host: str, region: str, name: str, interval: int) -> str:
    """First-boot script: fetch the commit and hand over to scripts/aws/instance_setup.sh.

    Holds no secret -- user-data is readable through the EC2 API and from the instance's
    metadata service. Every interpolated value is validated here because this is a root
    shell script assembled from command-line arguments.
    """
    if not _SHA.match(sha):
        raise DeployError(f"not a full commit id: {sha!r}")
    if host and not _HOST.match(host):
        raise DeployError(f"not a host name: {host!r}")
    if not _NAME.match(name):
        raise DeployError(f"--name must match {_NAME.pattern}")
    if not re.match(r"^[a-z]{2}(-[a-z]+)+-\d$", region):
        raise DeployError(f"not an AWS region: {region!r}")
    if not 15 <= int(interval) <= 3600:
        raise DeployError("--interval must be between 15 and 3600 seconds")
    return f"""#!/bin/bash
set -euxo pipefail
exec > >(tee -a /var/log/live-map-setup.log) 2>&1

cat > /etc/live-map.env <<'EOF'
LIVE_MAP_HOST={host}
LIVE_MAP_REGION={region}
LIVE_MAP_PARAMS=/{name}
LIVE_MAP_INTERVAL={int(interval)}
EOF

export DEBIAN_FRONTEND=noninteractive
apt-get update -q
apt-get install -y -q git

git init -q /opt/live-map/repo
cd /opt/live-map/repo
git remote add origin {REPO_URL}
git fetch -q --depth 1 origin {sha}
git checkout -q FETCH_HEAD
bash scripts/aws/instance_setup.sh
"""


def with_alias(config: dict, alias: str, certificate_arn: str) -> dict:
    """A distribution config that also answers to ``alias``, with that name's certificate."""
    config = copy.deepcopy(config)
    config["Aliases"] = {"Quantity": 1, "Items": [alias]}
    config["ViewerCertificate"] = {
        "ACMCertificateArn": certificate_arn,
        "SSLSupportMethod": "sni-only",
        "MinimumProtocolVersion": "TLSv1.2_2021",
    }
    return config


def dns_request(alias: str, cloudfront: str, validation: dict | None) -> str:
    """The records to ask whoever runs the alias's DNS for -- both at once, one ticket."""
    zone = alias.split(".", 1)[1] if "." in alias else alias
    rows = []
    if validation:
        rows.append((validation["Name"], validation["Type"], validation["Value"],
                     "lets AWS issue the certificate for the name"))  # fmt: skip
    rows.append((alias, "CNAME", cloudfront, "points the name at the site"))
    width = max(len(r[0]) for r in rows)
    lines = [
        f"  To use {alias}, ask whoever runs the {zone} DNS to add:",
        "",
        *[f"    {name:<{width}}  {kind}  {value}    ({why})" for name, kind, value, why in rows],
        "",
        "  Then run this same command again; it attaches the name once the certificate is",
        f"  issued. Until then the site stays at https://{cloudfront}.",
    ]
    return "\n".join(lines)


def update_commands(sha: str, host: str | None = None) -> list[str]:
    """The update script. ``host`` None keeps the name the instance serves HTTPS for;
    a name sets it, and "" clears it (CloudFront or plain HTTP in front instead)."""
    if not _SHA.match(sha):
        raise DeployError(f"not a full commit id: {sha!r}")
    if host and not _HOST.match(host):
        raise DeployError(f"not a host name: {host!r}")
    # AWS-RunShellScript runs these with /bin/sh (dash on Ubuntu): no bash-only options.
    serve = (
        []
        if host is None
        else [f"sed -i 's/^LIVE_MAP_HOST=.*/LIVE_MAP_HOST={host}/' /etc/live-map.env"]
    )
    return [
        "set -eu",
        *serve,
        "cd /opt/live-map/repo",
        f"git fetch -q --depth 1 origin {sha}",
        "git checkout -q -f FETCH_HEAD",
        "bash scripts/aws/instance_setup.sh",
    ]


# ---------------------------------------------------------------------------
# AWS
# ---------------------------------------------------------------------------


class Cloud:
    """This deployment's resources, found by tag."""

    def __init__(self, name: str, region: str, profile: str | None):
        if not _NAME.match(name):
            raise DeployError(f"--name must match {_NAME.pattern}")
        try:
            import boto3
            from botocore.exceptions import BotoCoreError, ClientError
        except ImportError:
            raise DeployError("boto3 is not installed: pip install boto3") from None

        self.name, self.region = name, region
        self.ClientError = ClientError
        try:
            session = boto3.Session(region_name=region, profile_name=profile)
            self.ec2 = session.client("ec2")
            self.iam = session.client("iam")
            self.ssm = session.client("ssm")
            self.cloudfront = session.client("cloudfront")
            # CloudFront only takes certificates from us-east-1, whatever the instance's region.
            self.acm = session.client("acm", region_name="us-east-1")
            self.account = session.client("sts").get_caller_identity()["Account"]
        except (BotoCoreError, ClientError) as exc:
            raise DeployError(
                f"no usable AWS credentials ({exc}).\n"
                "Set AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY (and AWS_SESSION_TOKEN for "
                "temporary ones), or configure ~/.aws/credentials and pass --profile."
            ) from exc
        self.tags = [{"Key": "Project", "Value": name}, {"Key": "Name", "Value": name}]
        self.filter = [{"Name": "tag:Project", "Values": [name]}]

    def tagged(self, kind: str) -> list[dict]:
        return [{"ResourceType": kind, "Tags": self.tags}]

    # -- lookups ------------------------------------------------------------

    def instance(self) -> dict | None:
        states = ["pending", "running", "stopping", "stopped"]
        pages = self.ec2.describe_instances(
            Filters=[*self.filter, {"Name": "instance-state-name", "Values": states}]
        )
        found = [i for r in pages["Reservations"] for i in r["Instances"]]
        if len(found) > 1:
            ids = ", ".join(i["InstanceId"] for i in found)
            raise DeployError(f"more than one instance is tagged Project={self.name}: {ids}")
        return found[0] if found else None

    def distribution(self) -> dict | None:
        """This deployment's CloudFront distribution, found by its comment."""
        for page in self.cloudfront.get_paginator("list_distributions").paginate():
            for item in page["DistributionList"].get("Items", []):
                if item.get("Comment") == self.name:
                    return item
        return None

    def _managed_policy(self, kind: str, name: str) -> str:
        # Looked up by name rather than written in as ids, so nothing here goes stale.
        if kind == "cache":
            listing = self.cloudfront.list_cache_policies(Type="managed")["CachePolicyList"]
            pairs = [(i["CachePolicy"]["CachePolicyConfig"]["Name"], i["CachePolicy"]["Id"])
                     for i in listing.get("Items", [])]  # fmt: skip
        else:
            listing = self.cloudfront.list_origin_request_policies(Type="managed")[
                "OriginRequestPolicyList"
            ]
            pairs = [
                (i["OriginRequestPolicy"]["OriginRequestPolicyConfig"]["Name"],
                 i["OriginRequestPolicy"]["Id"])
                for i in listing.get("Items", [])
            ]  # fmt: skip
        for found, policy in pairs:
            if found == name:
                return policy
        raise DeployError(f"CloudFront has no managed policy named {name}")

    def ensure_distribution(self, origin: str) -> dict:
        found = self.distribution()
        if found is not None:
            if found["Origins"]["Items"][0]["DomainName"] != origin:
                current = self.cloudfront.get_distribution_config(Id=found["Id"])
                config = current["DistributionConfig"]
                config["Origins"]["Items"][0]["DomainName"] = origin
                self.cloudfront.update_distribution(
                    Id=found["Id"], IfMatch=current["ETag"], DistributionConfig=config
                )
            return found
        config = distribution_config(
            name=self.name,
            origin=origin,
            cache_policy=self._managed_policy("cache", "Managed-CachingDisabled"),
            request_policy=self._managed_policy("request", "Managed-AllViewerExceptHostHeader"),
            reference=f"{self.name}-{int(time.time())}",
        )
        return self.cloudfront.create_distribution_with_tags(
            DistributionConfigWithTags={
                "DistributionConfig": config,
                "Tags": {"Items": self.tags},
            }
        )["Distribution"]

    def project_certificates(self) -> list[str]:
        """Certificates this deployment requested, found by tag."""
        found = []
        for page in self.acm.get_paginator("list_certificates").paginate():
            for summary in page["CertificateSummaryList"]:
                arn = summary["CertificateArn"]
                tags = self.acm.list_tags_for_certificate(CertificateArn=arn).get("Tags", [])
                if {"Key": "Project", "Value": self.name} in tags:
                    found.append(arn)
        return found

    def ensure_certificate(self, alias: str) -> dict:
        """The certificate for ``alias``, requesting one if needed: arn, status, DNS record."""
        if not _HOST.match(alias):
            raise DeployError(f"not a host name: {alias!r}")
        arn = None
        for candidate in self.project_certificates():
            detail = self.acm.describe_certificate(CertificateArn=candidate)["Certificate"]
            if detail["DomainName"] == alias and detail["Status"] in (
                "PENDING_VALIDATION",
                "ISSUED",
            ):
                arn = candidate
        if arn is None:
            arn = self.acm.request_certificate(
                DomainName=alias,
                ValidationMethod="DNS",
                IdempotencyToken=re.sub(r"\W", "", alias)[:32],
                Tags=self.tags,
            )["CertificateArn"]
            say(f"  requested a certificate for {alias}")
        for _ in range(12):
            detail = self.acm.describe_certificate(CertificateArn=arn)["Certificate"]
            options = detail.get("DomainValidationOptions") or [{}]
            record = options[0].get("ResourceRecord")
            if record or detail["Status"] == "ISSUED":
                return {"arn": arn, "status": detail["Status"], "record": record}
            time.sleep(5)
        raise DeployError(
            "AWS has not produced the certificate's DNS record yet; run again shortly"
        )

    def attach_alias(self, distribution_id: str, alias: str, certificate_arn: str) -> bool:
        current = self.cloudfront.get_distribution_config(Id=distribution_id)
        config = current["DistributionConfig"]
        if alias in (config.get("Aliases") or {}).get("Items", []):
            return False
        self.cloudfront.update_distribution(
            Id=distribution_id,
            IfMatch=current["ETag"],
            DistributionConfig=with_alias(config, alias, certificate_arn),
        )
        return True

    def address(self) -> dict | None:
        found = self.ec2.describe_addresses(Filters=self.filter)["Addresses"]
        return found[0] if found else None

    def parameter(self, key: str) -> str:
        return f"/{self.name}/{key}"

    def has_parameters(self, *keys: str) -> bool:
        names = [self.parameter(k) for k in keys]
        return not self.ssm.get_parameters(Names=names)["InvalidParameters"]

    # -- parameters ---------------------------------------------------------

    def put_inputs(
        self,
        targets: str | None,
        globus: tuple[str, str] | None,
        site_secrets: dict[str, str] | None = None,
    ) -> None:
        def put(key: str, value: str, kind: str = "String") -> None:
            self.ssm.put_parameter(Name=self.parameter(key), Value=value, Type=kind, Overwrite=True)

        if targets is not None:
            put("targets", targets, "SecureString")
            say("  stored probe targets")
            if site_secrets:
                put("site-secrets", json.dumps(site_secrets), "SecureString")
                say(f"  stored {len(site_secrets)} site service-account secret(s)")
            else:
                # New targets with no routes: a secret left from an older list is one
                # nobody is using, and the least useful thing to keep around.
                self.ssm.delete_parameters(Names=[self.parameter("site-secrets")])
        if globus is not None:
            put("globus-client-id", globus[0])
            put("globus-client-secret", globus[1], "SecureString")
            say("  stored the Globus service account")

    # -- pieces -------------------------------------------------------------

    def ensure_role(self) -> str:
        role = f"{self.name}-instance"
        trust = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Principal": {"Service": "ec2.amazonaws.com"},
                    "Action": "sts:AssumeRole",
                }
            ],
        }
        try:
            self.iam.create_role(
                RoleName=role,
                AssumeRolePolicyDocument=json.dumps(trust),
                Description="APPFL Bio Suite live map: read its own parameters, accept SSM",
                Tags=self.tags,
            )
        except self.iam.exceptions.EntityAlreadyExistsException:
            pass
        # Lets `update` run a command on the instance, in place of SSH.
        self.iam.attach_role_policy(
            RoleName=role, PolicyArn="arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore"
        )
        resource = f"arn:aws:ssm:{self.region}:{self.account}:parameter/{self.name}/*"
        policy = {
            "Version": "2012-10-17",
            "Statement": [{"Effect": "Allow", "Action": "ssm:GetParameter", "Resource": resource}],
        }
        self.iam.put_role_policy(
            RoleName=role, PolicyName="read-own-parameters", PolicyDocument=json.dumps(policy)
        )
        try:
            self.iam.create_instance_profile(InstanceProfileName=role, Tags=self.tags)
        except self.iam.exceptions.EntityAlreadyExistsException:
            pass
        try:
            self.iam.add_role_to_instance_profile(InstanceProfileName=role, RoleName=role)
        except self.iam.exceptions.LimitExceededException:
            pass  # already holds the role
        return role

    def ensure_security_group(self) -> str:
        found = self.ec2.describe_security_groups(
            Filters=[{"Name": "group-name", "Values": [self.name]}]
        )["SecurityGroups"]
        if found:
            return found[0]["GroupId"]
        vpcs = self.ec2.describe_vpcs(Filters=[{"Name": "is-default", "Values": ["true"]}])["Vpcs"]
        if not vpcs:
            raise DeployError(
                f"{self.region} has no default VPC. Create one (aws ec2 create-default-vpc) "
                "or choose another --region."
            )
        group = self.ec2.create_security_group(
            GroupName=self.name,
            Description="APPFL Bio Suite live map: HTTP and HTTPS in, nothing else",
            VpcId=vpcs[0]["VpcId"],
            TagSpecifications=self.tagged("security-group"),
        )["GroupId"]
        self.ec2.authorize_security_group_ingress(
            GroupId=group,
            IpPermissions=[
                {
                    "IpProtocol": "tcp",
                    "FromPort": port,
                    "ToPort": port,
                    "IpRanges": [{"CidrIp": "0.0.0.0/0"}],
                    "Ipv6Ranges": [{"CidrIpv6": "::/0"}],
                }
                for port in (80, 443)
            ],
        )
        return group

    def ensure_address(self) -> dict:
        return self.address() or self.ec2.allocate_address(
            Domain="vpc", TagSpecifications=self.tagged("elastic-ip")
        )

    def image(self, instance_type: str) -> str:
        info = self.ec2.describe_instance_types(InstanceTypes=[instance_type])["InstanceTypes"][0]
        arch = "arm64" if "arm64" in info["ProcessorInfo"]["SupportedArchitectures"] else "amd64"
        # Canonical publishes the current image id here; no AMI id to go stale in this file.
        name = (
            f"/aws/service/canonical/ubuntu/server/24.04/stable/current/{arch}/hvm/ebs-gp3/ami-id"
        )
        return self.ssm.get_parameter(Name=name)["Parameter"]["Value"]

    def launch(self, *, instance_type: str, role: str, group: str, user_data: str) -> str:
        request = dict(
            ImageId=self.image(instance_type),
            InstanceType=instance_type,
            MinCount=1,
            MaxCount=1,
            IamInstanceProfile={"Name": role},
            SecurityGroupIds=[group],
            UserData=user_data,
            BlockDeviceMappings=[
                {
                    "DeviceName": "/dev/sda1",
                    "Ebs": {"VolumeSize": 20, "VolumeType": "gp3", "Encrypted": True},
                }
            ],
            MetadataOptions={"HttpTokens": "required", "HttpEndpoint": "enabled"},
            TagSpecifications=self.tagged("instance") + self.tagged("volume"),
        )
        # A just-created instance profile takes a few seconds to exist everywhere.
        for attempt in range(12):
            try:
                return self.ec2.run_instances(**request)["Instances"][0]["InstanceId"]
            except self.ClientError as exc:
                if "Invalid IAM Instance Profile" not in str(exc) or attempt == 11:
                    raise
                time.sleep(5)
        raise AssertionError("unreachable")


def fetch_status(base: str, timeout: float = 8.0) -> dict | None:
    try:
        with urllib.request.urlopen(f"{base}/map/endpoints.status.json", timeout=timeout) as reply:
            return json.load(reply)
    except (urllib.error.URLError, OSError, ValueError):
        return None


def print_status(status: dict) -> None:
    say(f"  endpoints checked {status.get('checked_at')}:")
    for site, entry in status.get("sites", {}).items():
        say(f"    {site:16} {entry['state']:9} {entry['online']}/{entry['endpoints']}")


def wait_for_site(ip: str, host: str, minutes: int) -> bool:
    """Poll until the prober has published. Returns whether HTTPS is answering too."""
    deadline = time.monotonic() + minutes * 60
    status = None
    while time.monotonic() < deadline:
        status = fetch_status(f"http://{ip}")
        if status:
            break
        say(f"  installing… ({int((deadline - time.monotonic()) / 60)} min left before giving up)")
        time.sleep(30)
    if not status:
        raise DeployError(
            f"the site did not come up within {minutes} minutes. The instance is still "
            "there; its setup log is the EC2 console output:\n"
            "    python scripts/aws/live_map.py status --log"
        )
    print_status(status)
    if not host:
        return False
    # The certificate is requested on first start and usually lands within a minute.
    for _ in range(12):
        if fetch_status(f"https://{host}"):
            return True
        time.sleep(10)
    return False


def wait_for_distribution(domain: str, minutes: int = 25) -> bool:
    """Poll the page through CloudFront until it answers. New distributions take minutes."""
    deadline = time.monotonic() + minutes * 60
    while time.monotonic() < deadline:
        if fetch_status(f"https://{domain}"):
            return True
        say(f"  CloudFront deploying… ({int((deadline - time.monotonic()) / 60)} min left)")
        time.sleep(30)
    return False


def serve_through_cloudfront(cloud: Cloud, ip: str, alias: str | None = None) -> str | None:
    """Create (or find) the distribution and wait for it; attach ``alias`` when it can.

    Returns the best address that answers: the alias once its certificate is issued and
    its DNS points here, else the cloudfront.net name, else None.
    """
    distribution = cloud.ensure_distribution(ec2_public_dns(ip, cloud.region))
    domain = distribution["DomainName"]
    say(f"  CloudFront {domain}")
    if not wait_for_distribution(domain):
        return None
    if not alias:
        return f"https://{domain}"
    certificate = cloud.ensure_certificate(alias)
    if certificate["status"] != "ISSUED":
        say(dns_request(alias, domain, certificate["record"]))
        return f"https://{domain}"
    if cloud.attach_alias(distribution["Id"], alias, certificate["arn"]):
        say(f"  attached {alias}; CloudFront takes a few minutes to spread it")
    for _ in range(20):
        if fetch_status(f"https://{alias}"):
            return f"https://{alias}"
        time.sleep(30)
    say(f"  {alias} is attached but not reaching the site yet. Its DNS record must read:")
    say(f"    {alias}  CNAME  {domain}")
    return f"https://{domain}"


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------


def gather_inputs(args, cloud: Cloud | None, *, required: bool):
    targets, site_secrets = None, None
    if args.targets.is_file() or (required and not (cloud and cloud.has_parameters("targets"))):
        targets, warnings, secret_envs = read_targets(args.targets, public_site_ids())
        for warning in warnings:
            say(f"  warning: {warning}")
        missing = [name for name in secret_envs if not os.environ.get(name)]
        if missing and not required and cloud is not None and cloud.has_parameters("targets"):
            # A content-only update should not need a partner's secret in the shell. The
            # stored list and its secrets stay as they are, consistent with each other.
            say(f"  keeping the stored probe list (export {', '.join(missing)} to replace it)")
            targets = None
        else:
            site_secrets = read_site_secrets(secret_envs)
    globus = read_globus_client(args.globus_client)
    if globus is None and required:
        if not (cloud and cloud.has_parameters("globus-client-id", "globus-client-secret")):
            raise DeployError(
                f"no Globus service account for the prober. Register one at app.globus.org "
                f"(Settings -> Developers -> Register a service account) and save it as\n"
                f'    {args.globus_client}:  {{"client_id": "...", "client_secret": "..."}}\n'
                "or set LIVE_MAP_GLOBUS_CLIENT_ID and LIVE_MAP_GLOBUS_CLIENT_SECRET."
            )
    return targets, globus, site_secrets


def cmd_deploy(args) -> int:
    sha = resolve_ref(args.ref)
    if args.dry_run:
        gather_inputs(args, None, required=True)
        host = args.domain or ""
        front = (
            "your domain" if args.domain
            else "plain HTTP" if args.no_cloudfront
            else f"CloudFront as {args.alias}" if args.alias
            else "CloudFront"
        )  # fmt: skip
        target = f"{args.instance_type} in {args.region}"
        say(f"would deploy commit {sha[:10]} to {target}, via {front}\n")
        say(render_user_data(sha=sha, host=host, region=args.region, name=args.name,
                             interval=args.interval))  # fmt: skip
        return 0

    cloud = Cloud(args.name, args.region, args.profile)
    say(f"AWS account {cloud.account}, {args.region}")
    if cloud.instance():
        say("already deployed -- refreshing it instead.")
        return cmd_update(args, cloud=cloud, sha=sha)

    targets, globus, site_secrets = gather_inputs(args, cloud, required=True)
    cloud.put_inputs(targets, globus, site_secrets)
    role = cloud.ensure_role()
    group = cloud.ensure_security_group()
    address = cloud.ensure_address()
    ip = address["PublicIp"]
    host = args.domain or ""
    say(f"  address {ip}" + (f"  ({host})" if host else ""))
    if args.domain:
        say(f"  point an A record for {args.domain} at {ip}; the certificate follows it")

    instance = cloud.launch(
        instance_type=args.instance_type,
        role=role,
        group=group,
        user_data=render_user_data(
            sha=sha, host=host, region=args.region, name=args.name, interval=args.interval
        ),
    )
    say(f"  instance {instance} launching (commit {sha[:10]})")
    cloud.ec2.get_waiter("instance_running").wait(InstanceIds=[instance])
    cloud.ec2.associate_address(InstanceId=instance, AllocationId=address["AllocationId"])

    say("  waiting for the install; the pinned suite takes several minutes")
    https = wait_for_site(ip, host, args.wait)
    return finish(cloud, args, ip, host, https)


def finish(cloud: Cloud, args, ip: str, host: str, https: bool) -> int:
    """Put CloudFront in front if wanted, and print the address to use."""
    url = f"https://{host}" if https else f"http://{ip}"
    if not (args.domain or args.no_cloudfront):
        served = serve_through_cloudfront(cloud, ip, args.alias)
        if served:
            url = served
        else:
            say("  CloudFront is not answering yet; run `status` in a few minutes")
    elif host and not https:
        say(f"  (https://{host} is not answering yet; the certificate may still be issuing)")
    say("")
    say(f"live map: {url}/")
    return 0


def run_remote(cloud: Cloud, instance_id: str, commands: list[str], comment: str) -> dict:
    """Run a shell script on the instance through Systems Manager, and wait for it."""
    try:
        command = cloud.ssm.send_command(
            InstanceIds=[instance_id],
            DocumentName="AWS-RunShellScript",
            Comment=comment[:100],
            Parameters={"commands": commands, "executionTimeout": ["1800"]},
        )["Command"]["CommandId"]
    except cloud.ssm.exceptions.InvalidInstanceId:
        raise DeployError(
            f"{instance_id} is not reachable through Systems Manager yet. Its agent registers "
            "a few minutes after boot; try again shortly, or use `status --log`."
        ) from None
    while True:
        time.sleep(5)
        try:
            result = cloud.ssm.get_command_invocation(CommandId=command, InstanceId=instance_id)
        except cloud.ssm.exceptions.InvocationDoesNotExist:
            continue
        if result["Status"] not in ("Pending", "InProgress", "Delayed"):
            return result


LOG_COMMANDS = [
    "echo '== services'; systemctl is-active caddy live-map-prober || true",
    "echo '== prober'; journalctl -u live-map-prober --no-pager -n 40 -o cat",
    "echo '== web server'; journalctl -u caddy --no-pager -n 15 -o cat",
    "echo '== status file'; ls -la /opt/live-map/site/map/ | grep -i status || echo none",
    "echo '== setup (last lines)'; tail -n 25 /var/log/live-map-setup.log",
]


def cmd_logs(args) -> int:
    cloud = Cloud(args.name, args.region, args.profile)
    instance = cloud.instance()
    if instance is None:
        say("nothing is deployed.")
        return 1
    result = run_remote(cloud, instance["InstanceId"], LOG_COMMANDS, "live map logs")
    say((result.get("StandardOutputContent") or "").rstrip())
    error = (result.get("StandardErrorContent") or "").strip()
    if error:
        say("== stderr\n" + error[-3000:])
    return 0


def cmd_update(args, cloud: Cloud | None = None, sha: str | None = None) -> int:
    sha = sha or resolve_ref(args.ref)
    cloud = cloud or Cloud(args.name, args.region, args.profile)
    instance = cloud.instance()
    if instance is None:
        raise DeployError("nothing is deployed. Run `deploy`.")
    if instance["State"]["Name"] != "running":
        raise DeployError(f"instance {instance['InstanceId']} is {instance['State']['Name']}")

    targets, globus, site_secrets = gather_inputs(args, cloud, required=False)
    cloud.put_inputs(targets, globus, site_secrets)

    # A plain update leaves how the site is served alone; only an explicit choice changes it.
    host = args.domain if args.domain else "" if (args.alias or args.no_cloudfront) else None
    say(f"  updating {instance['InstanceId']} to commit {sha[:10]}")
    result = run_remote(
        cloud, instance["InstanceId"], update_commands(sha, host), f"live map -> {sha[:10]}"
    )
    tail = (result.get("StandardOutputContent") or "").strip().splitlines()[-3:]
    for line in tail:
        say(f"    {line}")
    if result["Status"] != "Success":
        say((result.get("StandardErrorContent") or "").strip()[-2000:])
        raise DeployError(f"the update {result['Status'].lower()}; the previous site is still up")
    if not (args.domain or args.no_cloudfront):
        serve_through_cloudfront(cloud, instance["PublicIpAddress"], args.alias)
    return cmd_status(args, cloud=cloud)


def cmd_status(args, cloud: Cloud | None = None) -> int:
    cloud = cloud or Cloud(args.name, args.region, args.profile)
    instance = cloud.instance()
    if instance is None:
        say("nothing is deployed.")
        return 1
    ip = instance.get("PublicIpAddress")
    say(f"instance {instance['InstanceId']}  {instance['InstanceType']}  "
        f"{instance['State']['Name']}  {ip or 'no address'}")  # fmt: skip
    if getattr(args, "log", False):
        output = cloud.ec2.get_console_output(InstanceId=instance["InstanceId"], Latest=True)
        say(output.get("Output", "(no console output yet)")[-6000:])
    if not ip:
        return 1
    status = fetch_status(f"http://{ip}")
    if not status:
        say("  the site is not answering")
        return 1
    print_status(status)
    base = f"http://{ip}"
    distribution = cloud.distribution()
    if distribution is not None:
        domain = distribution["DomainName"]
        names = (distribution.get("Aliases") or {}).get("Items", []) + [domain]
        answering = [n for n in names if fetch_status(f"https://{n}", timeout=8)]
        if answering:
            base = f"https://{answering[0]}"
        else:
            say(f"  CloudFront {domain} is not answering yet ({distribution['Status']})")
    say(f"live map: {base}/")
    return 0


def cmd_destroy(args) -> int:
    cloud = Cloud(args.name, args.region, args.profile)
    if not args.yes:
        raise DeployError(
            f"this deletes the CloudFront distribution, the instance, its address, role, "
            f"security group and stored parameters tagged Project={args.name} in "
            f"{args.region}. Re-run with --yes."
        )
    # A distribution can only be deleted once disabled, and disabling takes AWS minutes;
    # start it first so it runs while everything else is torn down.
    distribution = cloud.distribution()
    if distribution is not None:
        current = cloud.cloudfront.get_distribution_config(Id=distribution["Id"])
        if current["DistributionConfig"]["Enabled"]:
            current["DistributionConfig"]["Enabled"] = False
            cloud.cloudfront.update_distribution(
                Id=distribution["Id"],
                IfMatch=current["ETag"],
                DistributionConfig=current["DistributionConfig"],
            )
        say(f"  disabling CloudFront {distribution['DomainName']}")
    instance = cloud.instance()
    if instance:
        cloud.ec2.terminate_instances(InstanceIds=[instance["InstanceId"]])
        say(f"  terminating {instance['InstanceId']}")
        cloud.ec2.get_waiter("instance_terminated").wait(InstanceIds=[instance["InstanceId"]])
    for address in cloud.ec2.describe_addresses(Filters=cloud.filter)["Addresses"]:
        cloud.ec2.release_address(AllocationId=address["AllocationId"])
        say(f"  released {address['PublicIp']}")
    for group in cloud.ec2.describe_security_groups(
        Filters=[{"Name": "group-name", "Values": [cloud.name]}]
    )["SecurityGroups"]:
        cloud.ec2.delete_security_group(GroupId=group["GroupId"])
        say(f"  deleted security group {group['GroupId']}")

    role = f"{cloud.name}-instance"
    missing = cloud.iam.exceptions.NoSuchEntityException
    steps = [
        lambda: cloud.iam.remove_role_from_instance_profile(
            InstanceProfileName=role, RoleName=role
        ),
        lambda: cloud.iam.delete_instance_profile(InstanceProfileName=role),
        lambda: cloud.iam.detach_role_policy(
            RoleName=role, PolicyArn="arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore"
        ),
        lambda: cloud.iam.delete_role_policy(RoleName=role, PolicyName="read-own-parameters"),
        lambda: cloud.iam.delete_role(RoleName=role),
    ]
    for step in steps:
        try:
            step()
        except missing:
            pass
    say(f"  deleted role {role}")

    keys = ("targets", "site-secrets", "globus-client-id", "globus-client-secret")
    names = [cloud.parameter(k) for k in keys]
    cloud.ssm.delete_parameters(Names=names)
    say("  deleted stored parameters")

    if distribution is not None:
        say("  waiting for CloudFront to finish disabling (often 5-15 minutes)")
        deadline = time.monotonic() + 40 * 60
        while True:
            reply = cloud.cloudfront.get_distribution(Id=distribution["Id"])
            if reply["Distribution"]["Status"] == "Deployed":
                cloud.cloudfront.delete_distribution(Id=distribution["Id"], IfMatch=reply["ETag"])
                say(f"  deleted CloudFront {distribution['DomainName']}")
                break
            if time.monotonic() > deadline:
                say("  CloudFront is still disabling; run `destroy --yes` again later to delete it")
                break
            time.sleep(30)
    for arn in cloud.project_certificates():
        try:
            cloud.acm.delete_certificate(CertificateArn=arn)
            say("  deleted certificate")
        except cloud.acm.exceptions.ResourceInUseException:
            say("  certificate still in use by CloudFront; run `destroy --yes` again later")
    say("gone.")
    return 0


def parser() -> argparse.ArgumentParser:
    top = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--region", default=os.environ.get("AWS_REGION", DEFAULT_REGION))
    common.add_argument("--profile", default=None, help="AWS profile in ~/.aws")
    common.add_argument("--name", default=DEFAULT_NAME, help="tags and names every resource")

    inputs = argparse.ArgumentParser(add_help=False)
    inputs.add_argument("--ref", default=None, help="commit to serve (default: HEAD)")
    inputs.add_argument("--targets", type=Path, default=DEFAULT_TARGETS)
    inputs.add_argument("--globus-client", type=Path, default=DEFAULT_GLOBUS)

    serving = argparse.ArgumentParser(add_help=False)
    front = serving.add_mutually_exclusive_group()
    front.add_argument(
        "--domain", default=None, help="serve HTTPS for your own host name from the instance"
    )
    front.add_argument(
        "--no-cloudfront", action="store_true", help="plain HTTP on the address only"
    )
    front.add_argument(
        "--alias",
        default=None,
        help="your own name for the CloudFront address, e.g. map.appfl.ai: prints the DNS "
        "records to ask for, then attaches the name on a later run",
    )

    commands = top.add_subparsers(dest="command", required=True)
    deploy = commands.add_parser("deploy", parents=[common, inputs, serving], help="create it")
    deploy.add_argument("--instance-type", default=DEFAULT_INSTANCE_TYPE)
    deploy.add_argument("--interval", type=int, default=60, help="seconds between probes")
    deploy.add_argument("--wait", type=int, default=25, help="minutes to wait for the install")
    deploy.add_argument("--dry-run", action="store_true", help="check inputs, call nothing")
    deploy.set_defaults(run=cmd_deploy)

    update = commands.add_parser(
        "update", parents=[common, inputs, serving], help="refresh it in place"
    )
    update.set_defaults(run=cmd_update)

    status = commands.add_parser("status", parents=[common], help="is it up")
    status.add_argument("--log", action="store_true", help="print the instance's console log")
    status.set_defaults(run=cmd_status)

    logs = commands.add_parser("logs", parents=[common], help="the server's recent logs")
    logs.set_defaults(run=cmd_logs)

    destroy = commands.add_parser("destroy", parents=[common], help="delete all of it")
    destroy.add_argument("--yes", action="store_true")
    destroy.set_defaults(run=cmd_destroy)
    return top


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        return args.run(args)
    except DeployError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
