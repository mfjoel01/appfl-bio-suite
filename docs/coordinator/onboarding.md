# Package and register partners

One federation YAML describes the experiment and participants. Use the existing
`federation.yaml.example` as the starting point. Endpoint UUIDs may be omitted until the
partner registers. The launch configuration builder refuses sites without a UUID.

## Generate bundles

```bash
appfl-bio-suite experiment package gwas --federation local/federation.yaml \
  --out local/onboarding --wheel dist/appfl_bio_suite-0.2.1-py3-none-any.whl
```

Build the wheel first with `python -m build`. Including `--wheel` makes the ZIP carry
the setup tool and experiment code; its checksum is bound into the specification.
Omit `--wheel` to install from the pinned published release instead.

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
settings and hands off MEP activation explicitly. The single-user Globus pattern is
still under separate development and will be incorporated after validation.

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

This rejects a stale/wrong-site receipt, verifies that experiment settings have not
changed, and submits a real worker probe using your Globus login. It checks the service
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

Build a wheel and source distribution with `python -m build`. Both carry the partner-only
documents, templates and constraints. Bundle generation works from an installed wheel;
a source checkout is no longer required. GitHub release wheels can be used before PyPI
publication. Publish the matching version/tag using the existing release procedure.

Until this branch has a published release, install the built wheel into a Python 3.12
environment and include it with `experiment package --wheel`. Partners install the
bootstrap tool from the ZIP and setup uses that same wheel automatically.
`site setup --install-source` remains available for a separately supplied artifact.
Do not use the old published 0.2.1 tag to test these new commands. No release/tag is created by this change.

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

The [single-user service-account trial](service-account-trial.md) remains a separate
experimental path. Its eventual integration must retain bundle revision and readiness
checks while replacing the MEP activation/identity-mapping step.
