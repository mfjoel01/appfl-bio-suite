# Experimental: Hyperstack and Argonne single-user endpoint trial

> **Experimental implementation, not a production onboarding default.**
> Hyperstack–Argonne execution and scheduler validation are pending.
> Dashboard account minting and automated partner onboarding are not implemented.

## Pilot acceptance checks

- [ ] Run owner execution and both cross-site denial checks between Hyperstack and Argonne.
- [ ] Complete two APPFL training rounds using an Argonne scheduler allocation.
- [ ] Exercise endpoint restart, failed submissions, retry, cancellation, and worker cleanup.
- [ ] Test secret rotation/revocation with independent token caches and running work.
- [ ] Finish recoverable credential provisioning and teardown before dashboard integration.
- [ ] Add root-free partner bundles and operator guidance before changing onboarding defaults.

This implementation is for review and controlled testing; these unchecked items are not completed
by the local training tests below.

This branch supports an opt-in coordinator route for distinct service accounts. The
default multi-user flow remains available. The dashboard does not yet mint these
accounts automatically.

Use two Globus Auth service accounts: one for the Hyperstack site and one for the
Argonne site. Each endpoint and the coordinator authenticate as that site's identity.
Separate endpoint/coordinator secrets for the same identity help rotation, but have
the same permissions. Neither site receives the other site's credentials.

## What we need before connecting

Provide the Hyperstack SSH command (for example `ssh ubuntu@<public-IP>`) and the
Argonne hostname and username. If using existing endpoints instead, provide their
endpoint UUIDs, owning service-account client IDs, and expected Unix usernames.
Do not paste passwords, SSH private keys, or Globus secrets into chat.

The specific Argonne system determines its scheduler, allocation, queue, and worker
initialization. Those settings cannot be inferred from “an Argonne node.” Run the
training worker in an appropriate allocation; a login-node health check is not a
scheduler-backed training test.

## Prepare the installations and endpoints

The coordinator and both workers need Python 3.12 and the branch's pinned dependencies:
APPFL 1.10.0, Compute SDK/endpoint 4.9.0, Globus SDK 4.4.1, and Torch 2.8.0.
Install the trial branch from a checkout or its locally built wheel. This development
work is not published on PyPI. For a checkout (see the [tracked setup examples](../../examples/service-account-trial/)
for cloning, a virtual environment, and the CPU installation variant):

```bash
python -m pip install -e . -c constraints.txt
```

The synthetic trial below requires no cohort download or experiment data staging.

At each site, start a new endpoint under the local Unix account intended to run jobs.
Set that site's `GLOBUS_COMPUTE_CLIENT_ID` and `GLOBUS_COMPUTE_CLIENT_SECRET` through
the site's protected environment/secrets mechanism. These must be set on every restart.
Use a dedicated Compute directory so existing endpoints are unaffected:

```bash
export GLOBUS_COMPUTE_USER_DIR="$HOME/.globus_compute_biosim_trial"
globus-compute-endpoint configure biosim-trial
# Configure the engine/provider for this machine before starting.
globus-compute-endpoint start biosim-trial </dev/null
```

For Hyperstack, a `GlobusComputeEngine` with `LocalProvider` and one worker is sufficient.
For Argonne, use the appropriate scheduler provider and worker environment. Do not
apply a privileged manager configuration or identity mapping: this is a single-user
endpoint. Existing [partner templates](../../partners/templates/) provide provider
examples, but are not a turnkey single-user deployment bundle.

Keep the process supervised/running during the trial. Record the UUID printed by each
endpoint. A Globus service account does not create a Unix account; the worker runs with
the local account's filesystem permissions.

## Coordinator routing manifest

Copy [the routing example](../../examples/service-account-trial/site-credentials.example.json)
to `local/site-credentials.json`, then replace the placeholders with real values.
The two real endpoint UUIDs are the keys:

```json
{
  "version": 1,
  "endpoints": {
    "10000000-0000-0000-0000-000000000001": {
      "client_id": "20000000-0000-0000-0000-000000000001",
      "client_secret_env": "BIOSIM_HYPERSTACK_SECRET",
      "expected_user": "ubuntu"
    },
    "10000000-0000-0000-0000-000000000002": {
      "client_id": "20000000-0000-0000-0000-000000000002",
      "client_secret_env": "BIOSIM_ARGONNE_SECRET",
      "expected_user": "your-argonne-username"
    }
  }
}
```

The sample UUIDs and usernames must be replaced. Client IDs are public identifiers;
secrets are supplied through the two named coordinator environment variables. Inline
secrets are rejected. Keep these variables out of experiment YAML and partner bundles.

```bash
export APPFL_BIO_SUITE_GLOBUS_CREDENTIALS="$PWD/local/site-credentials.json"
```

All requested endpoints must have routes and populated secret variables. There is no
fallback to the coordinator's personal login when a route is missing. Each route must
have a distinct client identity. Credential references are never inserted into configs
sent to workers. Apps have independent in-memory tokens and executors have independent
task groups; credentials can be reauthenticated on the next coordinator run.

## Test execution and isolation first

```bash
python -m appfl_bio_suite.core.site_trial --check-endpoints --timeout 600
```

This submits harmless account/runtime probes to both endpoints, checks the expected
local users and Python 3.12, and attempts both cross-site submissions. A successful
check requires owner execution to succeed and cross-site requests to return the specific
`403 ENDPOINT_ACCESS_FORBIDDEN` error. A timeout, offline endpoint, or unrelated error
does not count as isolation evidence. An allowlist must permit the probe function for
this test; `FUNCTION_NOT_PERMITTED` is not an identity-isolation result.

## Run two synthetic training rounds

Prepare a new output directory (the command refuses to overwrite an existing trial):

```bash
python -m appfl_bio_suite.core.site_trial --out-dir local/hyperstack-argonne-trial
python -m appfl_bio_suite.core.drivers.globus_compute_driver \
  --server-config local/hyperstack-argonne-trial/server.yaml \
  --client-config local/hyperstack-argonne-trial/clients.yaml \
  --get-sample-size
```

Each worker generates 16 synthetic examples, trains a small linear model for two steps
per round, and returns its parameters. The coordinator uses the real APPFL scheduler
and sample-size-weighted FedAvg for two rounds. Workers write temporary logs below
`/tmp/appfl-site-trial-<run-id>/`; the coordinator writes logs inside the trial directory.
The model and dataset source are shipped by APPFL, so these generated files are only
needed on the coordinator. A successful exit and both sites returning in each of the
two rounds are the acceptance criteria.

SDK shutdown cancels local futures; it does not guarantee termination of remote work
already running. If a run is interrupted or stalls, inspect/stop trial workers and
scheduler jobs explicitly before retrying.

## Use the same routing for a real experiment

Once the trial succeeds, keep the manifest in the coordinator environment and use:

```bash
appfl-bio-suite endpoint smoke --experiment gwas --federation local/federation.yaml
appfl-bio-suite run gwas --federation local/federation.yaml
```

Alternatively pass `--site-credentials local/site-credentials.json` explicitly to
`endpoint status`, `endpoint smoke`, or `run`. The manifest remains outside generated
worker configs. `--dry-run` validates routes and secret availability without authenticating
or submitting. Manifest routing applies only to the Globus Compute driver.

After the trial, stop only the trial endpoints/jobs. Delete their endpoint registrations
and retire trial secrets when no longer needed; leave existing endpoints intact. The
experimental provisioning probe is not yet a production credential lifecycle manager.

## Software validation, 2026-11-05

For the focused PR against current `main`, 114 targeted tests passed, plus a
separate regression check for rejecting site credentials on serial simulations.
Lint, formatting, diff, and reusability checks also passed. These are local checks;
the pilot acceptance checks above remain outstanding.

- 62 focused checks passed against published APPFL 1.10.0 in the existing environment.
- A fresh virtual environment installed the actual suite wheel and published APPFL wheel.
  All 25 routing, synthetic training, endpoint-smoke, and probe tests passed against
  that installed package; `pip check` reported no broken requirements.
- The clean training check used Python 3.12, Torch 2.8.0+cpu, torchvision 0.23.0+cpu,
  NumPy 1.26.4, Compute SDK/endpoint 4.9.0, Globus SDK 4.4.1, and dill 0.3.9.
- Ruff, formatting, diff checks, and the repository reusability check passed.

Training validation used the real APPFL server, communicator, client training, and
aggregation, with an in-process test transport. It is not evidence of successful SSH
deployment, cross-machine authorization, message queues, or Argonne scheduler execution.
Those are the next pilot acceptance checks.

The [tracked CPU constraints](../../examples/service-account-trial/constraints-cpu.txt)
and [installation steps](../../examples/service-account-trial/README.md) reproduce the
clean test environment without a locally built bundle.
Its `+cpu` Torch build is intentional. The repository's standard pin check still requires
the exact production Torch version string; use the direct synthetic-trial driver above
for this CPU pilot, rather than treating it as a production experiment environment.
