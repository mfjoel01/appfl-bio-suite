# Experimental: Evaluating single-user service-account endpoints

> **Experimental implementation, not a production onboarding default.**
> Hyperstack–Argonne execution and scheduler validation are pending.
> Dashboard account minting and automated partner onboarding are not implemented.

## Development update: 2026-11-04

The branch now contains opt-in per-site routing in the training driver, endpoint status,
and endpoint smoke commands. The coordinator manifest references environment variables;
credentials are not added to worker configs. Each site uses an independent app, executor,
and task group while retaining APPFL's central scheduler and aggregation loop.

See the [Hyperstack–Argonne trial runbook](service-account-trial.md) for the isolated
two-site execution check and synthetic two-round training preparation commands. The
cross-machine trial has not run yet. The September findings below are historical;
the coordinator routing described there as missing is now implemented for testing.

## Recommendation

Use **one Globus Auth service account per site within each federation**, with a
single-user Globus Compute endpoint at that site. Treat the site/federation pair as
the default credential boundary; split further when experiments have different data
permissions or local execution accounts.

This keeps partner setup simpler: export two credentials and start an endpoint without
running a privileged multi-user manager or configuring identity mapping. A site receives
only its own service-account credentials. The coordinator needs credentials for each site,
or a separately validated delegation mechanism.

One federation-wide identity with separate secrets is appropriate only when every
participant intentionally shares the same execution and endpoint-ownership authority.
**Separate secrets allow separate rotation; they do not create separate permissions.**

| Property | One identity per federation | One identity per site/federation |
| --- | --- | --- |
| Partner endpoint needs root | No | No |
| Partner credential reaches | Endpoints owned by the shared identity | Endpoints owned by that site's identity |
| Attribution in Compute | Shared federation identity | Site identity; not individual people or secrets |
| Rotation | Separate secrets possible | Separate secrets possible per site |
| Coordinator integration | One authenticated client | Route each site to its authenticated client and task group |
| Credential compromise | Potential federation-wide access | Site/federation boundary, unless more permissions were granted |
| Local execution identity | Account running the endpoint | Account running the endpoint |

A Globus service account is an OAuth identity. It does **not** create or select a Unix
account. Running the endpoint as a dedicated, restricted local account remains desirable;
obtaining that account may require local IT even though running the endpoint does not
require root. The smoke-test user/UID still matters: compare it with the intended local
execution account in either mode.

## Online research, checked 2026-09-25

Globus explicitly supports endpoint authentication through
`GLOBUS_COMPUTE_CLIENT_ID` and `GLOBUS_COMPUTE_CLIENT_SECRET`. Its documentation also
states that an unprivileged single-user endpoint accepts its owner's tasks and runs
under the starting local account. See the pinned
[4.9 endpoint guide](https://globus-compute.readthedocs.io/en/4.9.0/endpoints/endpoints.html#client-identities)
and the [current guide](https://globus-compute.readthedocs.io/en/stable/endpoints/endpoints.html#starting-the-endpoint).
The current guide resolved to 4.17.1 during this review; the experiment below used the
repository's 4.9.0 pin, without upgrading it.

Service accounts have their own identities and permissions, independent of the human
who created them. Merely creating the accounts under one Auth project does not grant
the coordinator submission rights to all of their endpoints. See
[Understanding Service Accounts](https://globus-sdk-python.readthedocs.io/en/stable/user_guide/usage_patterns/service_accounts/understanding_service_accounts.html).

The [Auth API reference](https://docs.globus.org/api/auth/reference/#access-control)
documents a 30-minute administrator session requirement for editing developer resources.
It separately documents client-credential creation/deletion and token revocation.
Do not assume an ordinary stored human refresh token enables unattended provisioning,
or that deleting a client secret proves that every previously issued token has stopped
working. The exact revocation behavior must be tested with independent credential caches.

The [4.9 configuration reference](https://globus-compute.readthedocs.io/en/4.9.0/endpoints/config_reference.html)
documents function allowlists and subscription-gated additional endpoint administrators.
An administrative role should not be assumed to provide the desired single-user
submission delegation without a specific test.

## Live validation on 2026-09-25

Checkout: `globus-service-account-endpoints`, starting commit `3278522`, also verified
against the remote branch. Existing uncommitted changes to the probe and endpoint smoke
test were preserved.

Environment: Python 3.12, globus-compute-sdk/endpoint 4.9.0, globus-sdk 4.4.1.
The installed APPFL reports 1.10.0 and is an editable checkout outside this repository;
APPFL integration observations below describe that installed source and this branch's
driver, not a newly installed pristine wheel.

### Fresh authentication

Using a separate in-memory token cache for every credential:

| Existing probe | Credential | Result |
| --- | --- | --- |
| `probe-shared` | alpha | Authenticated as shared client `c8791dcf…` |
| `probe-shared` | beta, previously marked revoked | HTTP 401 `UNAUTHORIZED` |
| `probe-persite` | alpha | Authenticated as `2d126934…` |
| `probe-persite` | beta | Authenticated as `496ca6b3…` |

Reading the shared client's credential list using the existing human login returned
HTTP 403 with a session-policy requirement. No new secrets were minted, and none of the
existing probe identities or credentials were deleted.

### Execution and authorization matrix

Three temporary endpoints were started on this host, one under each usable identity.
They used `ThreadPoolEngine(max_workers=1)`, independent configuration/token directories,
and a harmless function returning 42 and the Unix UID. These tests exercise
real Globus Auth, Compute, and worker round trips; they do not exercise a partner HPC
scheduler, patient data, or an APPFL training round.

| Caller | Target owner | Result |
| --- | --- | --- |
| Shared alpha | Shared identity | Executed, returned 42 |
| Per-site alpha | Per-site alpha | Executed, returned 42 |
| Per-site beta | Per-site beta | Executed, returned 42 |
| Per-site beta | Per-site alpha | HTTP 403 `ENDPOINT_ACCESS_FORBIDDEN` |
| Per-site alpha | Per-site beta | HTTP 403 `ENDPOINT_ACCESS_FORBIDDEN` |
| Shared alpha | Per-site alpha | HTTP 403 `ENDPOINT_ACCESS_FORBIDDEN` |

Every successful task ran with UID 41746, the local user that started the endpoints.
All temporary processes were stopped and all six temporary endpoint registrations
(from the initial run and corrected rerun) were deleted. Existing endpoints were left
alone.

Local evidence, intentionally gitignored:

- Runner: `local/globus_service_account/run_matrix_20260925.py`
- Initial run: `local/globus_service_account/matrix-20260925-152433/results.json`
- Corrected matrix: `local/globus_service_account/matrix-20260925-152600/results.json`

### Integration failure found by the first run

The first identity's task succeeded, but subsequent identities' own-endpoint submissions
returned `404 TASK_GROUP_NOT_FOUND`. Inspection found that SDK 4.9
`Executor._default_task_group_id` is shared across the process. A new executor can inherit
a task group belonging to another identity even when given a different `Client`.

Setting the new executor's public `task_group_id` property to `None` **before submission**
allowed the service to create a group under the correct identity; the corrected matrix
above then passed. Production routing must retain separate task groups/result handling
per identity. Merely constructing multiple clients does not solve this issue.

### Limits of the historical evidence

The earlier 2026-09-16 notes reported that two shared secrets authenticated as the same
identity, that beta could execute on alpha's endpoint, and that an already-issued token
still worked after beta's secret was deleted. The saved state still associates both
shared credentials with the same client ID.

The two-live-secret comparison and immediate post-deletion token behavior were **not
repeated** today: beta is already revoked, and administrator credential access requires
a fresh interactive session. Furthermore, the old probe used the SDK's per-client token
cache. Without a transcript establishing cache isolation, a cached request cannot prove
which secret originally obtained that token. Treat that historical token result as
provisional, not a measured revocation guarantee.

The identity-level access model and today's per-site denial matrix support the
per-site recommendation independently of that historical result.

## Probe changes and local checks

`scripts/globus_service_account_probe.py` now:

- Uses explicit `ClientApp` instances for per-site authentication, without switching
  process-wide credential variables.
- Stores tokens in a protected probe database, with a namespace for each client and
  credential ID. Two secrets for one identity no longer reuse one another's tokens.
- Makes `--fresh` use an empty in-memory cache without deleting the operator's normal
  endpoint tokens or destroying the cached-token comparison.
- Uses the same app for Auth identity checks and Compute operations.
- Clears inherited executor task groups before submission.
- Labels a timeout inconclusive instead of reporting it as an authorization refusal.
- Prints hashed secret fingerprints and refuses to overwrite an existing provisioned set.

The probe is still an experimental tool, not a production credential lifecycle service.
Provisioning needs durable recovery from partial API failures, and teardown needs stronger
ownership checks and recovery-state preservation before automated onboarding adopts it.

Validation commands:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m pytest \
  tests/test_globus_service_account_probe.py \
  tests/test_endpoint_smoke.py tests/test_preflight.py tests/test_docs.py -q
python -m ruff check scripts/globus_service_account_probe.py \
  tests/test_globus_service_account_probe.py
```

Result: **42 tests passed** and Ruff passed. Fresh `whoami` and cached `endpoints`
commands also completed against the live service using the corrected probe.

## What is still needed for a federation run

The current driver constructs APPFL's `GlobusComputeServerCommunicator`. The installed
communicator creates one `Client` and one `Executor`, then passes that executor to every
site. It does not route credentials by endpoint. Possessing all the secrets therefore
does **not** make the per-site design a configuration-only change.

The next implementation should:

1. Resolve each endpoint to a coordinator-side credential reference, outside partner
   bundles and experiment configuration sent to workers.
2. Construct explicit per-site apps/clients/executors, with independent task groups,
   function registration, result collection, retries, cancellation, and shutdown.
3. Keep the experiment's model aggregation centralized; separate authentication must not
   create separate federated training runs.
4. Make partner setup and readiness checks distinguish single-user execution from
   multi-user identity mapping, while checking the expected Unix account in both.
5. Validate a complete two-site APPFL round and a scheduler-backed worker round trip.

A coordinator can use a different secret from the site's endpoint secret for each
site identity. This helps rotation, but both secrets still carry the same identity
permissions, and it does not distinguish coordinator actions from site actions in
Compute's identity-level audit trail. Individual-user attribution requires a different
delegation or multi-user design.

## Function restrictions and revocation

`allowed_functions` is useful, but it is not a sandbox or a sufficient boundary against
a credential that also owns the endpoint. Endpoint-owner control and function arguments
must be included in the threat model.

In this installation, APPFL's communicator ships `dataset_source`, and its client agent
can load dataset, model, and trainer source from configuration. An allowlisted training
wrapper that accepts executable source still permits caller-supplied code. A design
promising “reviewed functions only” must constrain arguments, code loading, serialization,
and endpoint configuration authority as well as function UUIDs.

Do not claim that deleting an Auth client automatically removes Compute endpoints or
kills queued/running jobs; that lifecycle was not tested. Similarly, deleting the
departing site's endpoint does not prevent a still-valid shared-identity token from
targeting other sites' endpoints.

Offboarding must distinguish new credential grants, existing tokens, endpoint/worker
processes, scheduler jobs, and accessible data. For immediate containment, the site must
be able to stop its endpoint and workers/jobs and remove local data access. Per-site
identities make that response local to the affected site/federation boundary; shared
identities make selective containment substantially harder.
