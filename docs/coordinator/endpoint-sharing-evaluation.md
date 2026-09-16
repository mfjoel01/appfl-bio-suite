# Evaluating the shared service-account endpoint design

A proposal has been circulated for sharing Globus Compute endpoints across a federation
without needing root on every partner's login node:

1. APPFL registers a Globus account of its own.
2. Per federation, a *service account* (a Globus Auth `client_identity` client) is created
   under it.
3. Its client ID and secret are distributed to the federation.
4. Every partner registers an ordinary single-user endpoint using those credentials.
5. All endpoints then belong to the one service account, so anything holding the
   credentials can submit to any of them.

Globus have confirmed the mechanism works, and suggested minting a **separate secret per
partner** so that any one partner can be revoked by deleting their secret.

This document records what was actually measured rather than assumed. The probe that
produced it is `scripts/globus_service_account_probe.py`.

> **This repository already tried this.** `architecture.md` records: *"An earlier design
> used a shared confidential client — one credential distributed to every site. It worked,
> but anyone holding it could submit arbitrary functions to any site, and everything ran
> under a single machine identity with no attribution."* The proposal is a return to that
> design plus per-partner secrets. The question worth answering is therefore narrow: **do
> per-partner secrets fix the two things that caused it to be dropped?**

---

## What the proposal gets right

**The mechanism is real and first-class, not a workaround.** In globus-compute 4.9.0,
`sdk/auth/globus_app.py` branches on the two environment variables: with both
`GLOBUS_COMPUTE_CLIENT_ID` and `GLOBUS_COMPUTE_CLIENT_SECRET` set it builds a
`globus_sdk.ClientApp` (client-credentials grant) instead of a `UserApp`. The endpoint CLI
reaches the same function through `globus_compute_endpoint/auth.py`, so `configure`,
`start`, `whoami` and `delete` all authenticate as the client identity. No browser login
happens on the partner's machine, which matters on a headless login node.

**It genuinely removes the root requirement.** Multi-user endpoints gate on
`is_privileged()` (`endpoint/utils/__init__.py:120`), which is true only for uid 0, the
literal user `root`, or a process holding the relevant Linux capabilities. Nothing about
the service-account design needs any of that: each partner runs a normal single-user
endpoint as themselves.

**It does not disturb an existing login.** Tokens are namespaced by
`_resolve_namespace()` (`sdk/auth/token_storage.py`) as `clientprofile/<env>/<client_id>`
when client credentials are set, versus `user/<env>` otherwise. Both live in the same
`~/.globus_compute/storage.db` without colliding, so a partner can run their own endpoints
and a federation endpoint on one machine.

**Per-secret revocation is supported.** `globus_sdk.AuthClient` exposes
`create_client_credential`, `get_client_credentials` and `delete_client_credential`. A
client identity may hold several credentials at once, so the suggested one-secret-per-
partner scheme is implementable exactly as described.

---

## What it costs

### 1. Every partner can reach every other partner's cluster

Authorization for a single-user endpoint is decided **server-side, by identity**. The
endpoint daemon performs no owner check of its own. So when every endpoint is owned by one
client identity, the credential that lets a partner register their own endpoint is the same
credential that lets them submit arbitrary functions to everyone else's — and delete or
reconfigure those endpoints.

Per-partner secrets do not change this. All of them authenticate as the *same identity*;
they differ only in which string opens the door. This is the first of the two reasons the
design was dropped here, and the suggestion does not address it.

For a biomedical federation this is the crux. The current model's guarantee to a partner is
"we can run only what you have mapped one named identity to run, and you revoke us by
editing one line". The replacement guarantee is "any institution in the federation, and
anyone who has obtained any copy of any federation secret, can execute code on your cluster
against your data".

### 2. There is no per-partner attribution

A client-credentials token's subject is the client identity, not the credential used to
obtain it. Compute therefore records one actor for the whole federation. You can revoke a
partner individually but you cannot tell, from the service's records, which partner did
anything — including which one submitted the job that touched patient data.

Revocation granularity and audit granularity are not the same property, and the proposal
supplies only the first. *(Whether Globus Auth retains credential-level issuance records
that a project admin can retrieve is worth asking directly — it is not exposed through any
SDK call found here.)*

### 3. Tasks stop running as a service account

Under the multi-user design, identity mapping enforces that a task runs as a named local
account (`gwas_svc`, `flamby_svc`) regardless of who started the endpoint. Under this
design, a task runs as **whoever started the endpoint** — an ordinary partner's Unix user,
with that user's full filesystem access.

This is a silent change to something this codebase treats as evidence: the `user` field
returned by `where_am_i()` (`core/endpoint.py`) exists specifically to prove the partner's
mapping took effect. With service-account endpoints that field still returns a value, but
the value no longer means anything, and the smoke test no longer tests what its docstring
says it tests.

### 4. Secrets are long-lived, plaintext, and copyable

The secret must persist on every partner's filesystem — an env var in a shell profile or a
systemd unit — because the endpoint needs it at every restart to refresh tokens. It does
not expire. Nothing prevents a partner from forwarding it. A single compromised login node
exposes credentials that reach every site in the federation.

### 5. Revocation is not immediate

Deleting a credential stops that partner obtaining *new* tokens. It does not invalidate
access tokens already issued, and it does not change who owns the endpoints they already
registered — those stay under the federation identity. Cutting off a partner promptly means
deleting the credential *and* deleting or disowning their endpoints.

### 6. Minting secrets is not as automatable as it sounds

The suggestion is to mint secrets "during some registration or automated process". Measured
here: `POST /v2/api/clients` returned

```
403 FORBIDDEN — To access this project you must have an identity with admin privileges
in session within the last 30 minutes.
session_required_policies: [...]
```

The call is programmatic, but Globus Auth enforces a **session-freshness policy** on
project administration that a stored refresh token does not satisfy. A human with admin
rights has to have authenticated interactively within the last half hour. Fully unattended
onboarding therefore needs a different arrangement (for example making a service account an
admin of the project) — which is worth confirming with Globus before designing a
registration flow around it.

Relatedly, `create_project` rejects a project with no administrator
(`Either admin_ids or admin_group_ids has to be defined`). Use `admin_group_ids` with a
Globus group, so "who can mint and revoke federation secrets" is not one person's account.

---

## The variant that costs nothing and fixes most of it

Issue one client identity **per site**, not per federation.

Every property the proposal wants is kept: no root anywhere, no Globus group and no support
ticket, credentials created and revoked by the APPFL project owner, partners running
ordinary single-user endpoints. What changes is that a secret unlocks exactly one site.

* **Blast radius** collapses from the whole federation to one site. A leaked secret reaches
  one cluster.
* **Attribution** becomes real: distinct identities, distinct owners, distinct records.
* **Revocation** becomes complete — delete the client and the endpoints under it go with
  it, rather than being orphaned under a still-live shared identity.

The thing it gives up is any-to-any submission, and APPFL does not use that. The Globus
Compute driver (`core/drivers/globus_compute_driver.py`) has a central server dispatching
to each client endpoint; sites never submit to one another. So the coordinator, which holds
every secret because it minted them, can still reach every site — and no site can reach
another.

If any-to-any really is needed for some future topology, that is the point at which the
shared identity earns its cost, and it should be argued on that basis rather than adopted
by default.

---

## Regardless of shape: pin `allowed_functions`

Endpoints accept an `allowed_functions` allowlist of registered function UUIDs
(`endpoint/config/config.py:93`), enforced by the service rather than by the daemon. Set it
on every federation endpoint and a stolen or misused credential can still only run the
experiment's reviewed functions, not arbitrary code. Under the shared-identity shape this
is the only thing standing between one partner's credential and another partner's cluster,
so it stops being optional.

Note that `admins` — the supported way to give another identity administrative access to an
endpoint you own — requires an active Globus subscription. That gating is very likely why
this design is being considered at all, and is worth confirming before building around a
workaround.
