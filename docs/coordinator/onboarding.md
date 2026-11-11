# Package and register partners

One federation YAML describes the experiment and participants. Use the existing
`federation.yaml.example` as the starting point. Endpoint UUIDs may be omitted until the
partner registers. The launch configuration builder refuses sites without a UUID.

## Generate bundles

```bash
appfl-bio-suite experiment package gwas --federation local/federation.yaml \
  --out local/onboarding
```

Without `--wheel`, partners install the pinned release from PyPI
(`appfl-bio-suite[<extras>]==<this version>`) during setup. To ship the code inside the
ZIP instead, so that its checksum is bound into the specification, add the published
wheel itself:

```bash
pip download "appfl-bio-suite==0.3.2" --no-deps -d dist/
appfl-bio-suite experiment package gwas --federation local/federation.yaml \
  --out local/onboarding --wheel dist/appfl_bio_suite-0.3.2-py3-none-any.whl
```

Use the PyPI wheel rather than one you built locally, so that every partner's checksum
matches the published release.

This emits a ZIP per participating site. Use repeated `--site SITE_ID` options for a
subset. `partner-bundle` also emits the machine-readable specification for compatibility
with the existing single-site workflow. Send each partner its own ZIP.

Each ZIP contains a schema-versioned specification, exact runtime constraints, the
partner's client configuration and data requirements, a local settings template, and
the existing endpoint/data guides. The specification records a digest of experiment
settings; receipts bind to the complete specification. ZIP creation uses the generated
manifest's allowlist, excluding unrelated files left in the output directory. Checksums
identify revisions; they do not authenticate a sender. Distribute bundles through your
normal trusted collaboration channel.

The partner runs `appfl-bio-suite join site.zip`, or prepares `site.yaml` first and runs
`appfl-bio-suite join site.zip --local site.yaml`. The guided command remembers cluster
settings and hands off endpoint activation explicitly, for whichever endpoint mode the
bundle declares (below).

## Endpoint mode and host-local workers

Each site entry declares how its endpoint runs, and the bundle carries that declaration.
The partner never chooses it, and `join` never infers it from who runs the command:

```yaml
experiments:
  fine-mapping:
    sites:
      - site: site-east
        client_id: east
        endpoint_mode: single-user   # default: multi-user
```

- **`multi-user`** (default): an administrator with root activates a MEP that maps your
  Globus identity onto the service account. You reach it with your own login.
- **`single-user`**: the endpoint runs as the service account itself and registers under
  a **Globus service-account client you issue to that site**. It needs no root and no
  identity mapping. It accepts tasks only from that client, so you reach it through a
  `--site-credentials` route for the same client.

For a host with no batch scheduler (a cloud VM, a workstation), set the site's
`scheduler: local` and omit `partition`, `queue` and `account`. Workers then run on the
endpoint host through parsl's `LocalProvider`, and `join` asks only for the data and
work directories.

For a single-user site, mint one service-account client **per site** in a Globus Auth
project, never one shared across the federation: sites holding the same client can
submit to each other's endpoints. Create two secrets for it: one the site places in a
0600 file for `site activate --credentials-file`, one you keep for the coordinator. Send
the client ID and the site's secret through your normal secure channel, never in the
bundle.

Build the coordinator's routing manifest from the [service-account
guide](service-account-trial.md#coordinator-routing-manifest). A run needs a route for
**every** endpoint in it. A multi-user site in the same run as single-user ones gets an
explicit route saying you reach it with your own login:

```json
{
  "version": 1,
  "endpoints": {
    "<site-east endpoint UUID>": {
      "client_id": "<site-east service-account client ID>",
      "client_secret_env": "BIOSIM_SITE_EAST_SECRET",
      "expected_user": "finemap_svc"
    },
    "<site-north endpoint UUID>": {"coordinator_login": true, "expected_user": "finemap_svc"}
  }
}
```

There is still no fallback: an endpoint without a route is refused, so a forgotten
entry cannot send a site's tasks under the wrong identity.

The partner can change scheduler account/queue, local paths and network interface.
Everything experiment-owned comes from your federation configuration. A missing local
value is reported as a named validation error before installation. See
[partner onboarding](../partner/onboarding.md) for the partner commands.

## Register and verify

When a partner returns `registration.json`:

```bash
appfl-bio-suite experiment register registration.json \
  --bundle local/onboarding/gwas_Site1.zip \
  --federation local/federation.yaml \
  --updated-federation local/federation.site1.yaml \
  --out local/site1-ready.json --timeout 600
```

For a single-user site, add `--site-credentials local/site-credentials.json` (or set
`APPFL_BIO_SUITE_GLOBUS_CREDENTIALS`). The probe then authenticates as that endpoint's
route, and registration refuses to start without one rather than waiting on a 403.

This rejects a stale/wrong-site receipt, verifies that experiment settings have not
changed, and submits a real worker probe using your Globus login (or the site's route). It checks the service
account, runtime versions and directory access. Only success produces a ready record and
a new federation YAML containing the partner's endpoint UUID and local data/output paths.
Use that new federation file as input when registering the next partner. Existing output
files are never overwritten by the CLI. YAML comments are not preserved in the generated
configuration; the source remains intact.

Use the resulting federation file for the existing experiment preflight and launch
commands. Readiness covers infrastructure, not scientific data correctness; the existing
data and governance validation still applies. A queue timeout is not readiness: rerun
with a longer timeout after checking the queue. The CLI does not wait indefinitely for
executor shutdown after a timeout.

## Distribution and development

Releases are published on [PyPI](https://pypi.org/project/appfl-bio-suite/). The wheel
and source distribution both carry the partner-only documents, templates and
constraints, so bundle generation works from a plain PyPI install; a source checkout is
not required. Generate bundles with the same release your partners will run: a bundle
records the suite version that generated it, and `join` refuses a different one.

To test unreleased changes, build a wheel from your checkout with `python -m build`,
install it into a Python 3.12 environment, and include it with
`experiment package --wheel`. Partners install the bootstrap tool from the ZIP and setup
uses that same wheel automatically. `site setup --install-source` remains available for
a separately supplied artifact. Publish a new version with the
[release procedure](releasing.md) before using such changes in a real federation.

The v1 flow uses venv runtimes. Container-backed runtimes can be added behind the same
bundle/setup interface later. It does not provision cluster accounts, alter institutional
endpoints, or claim to eliminate the first administrator onboarding step.

## Validation boundary

Automated tests cover bundle/schema validation, guided input, reusable profiles,
resumable installation, endpoint ownership protection, receipt matching, and coordinator
configuration updates. Packaging is checked from an installed wheel outside the checkout.
Endpoint activation and remote execution are mocked in those tests; they do not establish
that a real partner's scheduler or network is configured correctly.

Before a first live deployment, exercise setup, administrator activation, local receipt
creation, and coordinator registration against a dedicated test MEP. Confirm the remote
probe runs as the intended service account on a scheduled worker with access to the
shared installation, data, and output directories. No live endpoint was started as part
of implementing this feature.

Single-user activation is covered the same way: refusal as root or as any other account,
refusal without service-account credentials or with a credentials file others can read,
endpoint ownership, restart of a managed endpoint, and registration through a route.
`globus-compute-endpoint configure`/`start` are mocked there too. Its first live
deployment (a cloud VM) is the validation step; the
[service-account trial](service-account-trial.md) records the cross-site isolation checks
that remain open. Both endpoint modes keep the bundle revision and readiness checks.
