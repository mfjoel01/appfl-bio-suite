# Standing up your own federation

**Start here.** This takes you from a clone to a running federated experiment across
institutions you recruit yourself. It assumes no prior contact with anyone involved in
building this suite.

You are the **coordinator**: the one site that dispatches work. Everyone else is a
**partner**: they run a compute endpoint, authorize your identity, and their cluster
executes tasks locally on data that never leaves it.

---

## What you are actually setting up

```
        YOU (coordinator)                  no data moves; only tasks and results
        ┌──────────────────┐
        │ driver process   │
        │ your Globus id   │
        │ aggregator       │
        └────────┬─────────┘
                 │  Globus Compute cloud relay (dispatch only)
     ┌───────────┼───────────┐
     ▼           ▼           ▼
  Partner A   Partner B   Partner C          each runs a MULTI-USER endpoint
  their data  their data  their data         that maps YOUR identity -> a local
  stays put   stays put   stays put          service account on THEIR cluster
```

Two things make this work, and both are worth understanding before you start:

**You run a single-user endpoint; partners run multi-user endpoints.** A single-user
endpoint accepts tasks only from the identity that registered it — exactly right for you,
since you are the only one dispatching. A multi-user endpoint can accept tasks from
*other* identities and map each to a local account. That is what each partner runs.

**Nothing secret is distributed.** The only thing you hand a partner is your Globus
identity string, which is not a credential. Knowing it grants nothing; the partner's own
endpoint configuration is what decides to honour it. There is no shared token to
distribute, rotate, or leak.

---

## Before you start

- **Python 3.12.** Not 3.11, not 3.13. The federation ships function bytecode between
  machines and it is only portable within a minor version.
- **A Globus account.** Free. <https://app.globus.org>
- **A machine that can stay up for the length of a run.** The driver process must outlive
  the work it dispatches. A login node in `tmux` is fine.
- **Partners who will cooperate.** Each needs an administrator willing to run a
  privileged process on their cluster. This is usually the long pole — start the
  conversation before you finish the technical setup.

You do **not** need: a shared filesystem, a VPN, matching operating systems, matching
schedulers, or any inbound network access to your machine.

---

## Step 1 — Install

```bash
git clone https://github.com/<your-org>/appfl-bio-suite.git
cd appfl-bio-suite

python -m venv .venv && source .venv/bin/activate
pip install -e ".[all]" -c constraints.txt
```

The `-c constraints.txt` is not optional decoration. It pins the versions every site in
the federation must share. Without it, pip resolves `globus-compute-sdk` to a release that
cannot start the endpoint you need.

Verify:

```bash
appfl-bio-suite preflight --check pins
```

Everything must pass. If it does not, fix it now — version skew does not fail at install
time, it fails as a deserialization error partway into a run on someone else's cluster.

**If you are running fine-mapping**, it also needs two C++ binaries that pip cannot
provide. Once, with network access:

```bash
scripts/fine-mapping/install_plink.sh      # PLINK 1.9 -> vendor/bin/
scripts/fine-mapping/install_susiex.sh     # SuSiEx    -> vendor/bin/
```

Both are coordinator-side only; no partner is asked to install either.
`appfl-bio-suite preflight --experiment fine-mapping` checks for them.

## Step 2 — Prove the install works, before involving anyone

Run a complete federation on your own machine. No partners, no endpoints, no network:

```bash
appfl-bio-suite simulate gwas --scenario ci-tiny --out /tmp/gwas-ci
appfl-bio-suite run gwas --config loopback --driver serial --data-root /tmp/gwas-ci
```

That generates a synthetic cohort, splits it across two simulated sites, runs the local
analysis at each, and combines the results, into `local/output/gwas-loopback/`. If it
completes, your install is correct end to end.

No federation config is needed for this, which is why it comes before Step 4: a loopback
run describes its own simulated sites. Once you do have a `local/federation.yaml`, the
same command uses it instead, with each site's paths repointed at this machine — they
otherwise name directories on partner clusters.

Every implemented experiment has the same pair of commands. For fine-mapping:

```bash
appfl-bio-suite simulate fine-mapping --scenario ci-tiny --out /tmp/fm-ci
appfl-bio-suite run fine-mapping --config loopback --driver serial --data-root /tmp/fm-ci
```

Do this first. Every hour spent debugging your own environment while a partner waits is an
hour of someone else's goodwill.

## Step 3 — Find your Globus identity

```bash
python -c "
from globus_compute_sdk import Client
c = Client()
print(c.login_manager.get_auth_client().oauth2_userinfo())
"
```

A browser flow opens the first time. You want the `preferred_username` — usually something
like `you@your-institution.edu`. Note the `sub` UUID too.

This string is what every partner will authorize. It is not a secret.

## Step 4 — Describe your federation

One file, and it is the only one you edit:

```bash
cp federation.yaml.example local/federation.yaml
$EDITOR local/federation.yaml
```

The example is fully worked and fictional. Replace:

- `coordinator.identity` — from Step 3
- `coordinator.identity_id` — the `sub` UUID
- `coordinator.organization` and `contact` — these appear in generated partner documents
- `sites` — one entry per participating institution, with their scheduler and account
- `experiments.<name>.sites` — who is running what, with their endpoint UUIDs

You will not have the endpoint UUIDs yet. Leave the placeholders; you fill them in when
partners send them back in Step 6.

`local/` is gitignored, permanently. Your real federation config, with live endpoint IDs,
never enters version control.

> **Everything in the suite reads from this file.** No generator, config, document
> template, or check contains a coordinator identity, an endpoint UUID, or a site name. If
> you ever find yourself editing code to describe your federation, that is a bug — please
> report it.

## Step 5 — Set up your own endpoint (optional)

Only if you are also a training site. If you only dispatch, skip this and delete the
`coordinator.endpoint` block.

```bash
globus-compute-endpoint configure my-endpoint
$EDITOR ~/.globus_compute/my-endpoint/config.yaml
```

Yours is **single-user** — no identity mapping, because it only ever accepts your own
tasks. See [reference-deployment.md](reference-deployment.md) for a complete working
example, and copy the `worker_init` lines from it: each one prevents a specific failure.

```bash
globus-compute-endpoint start my-endpoint
appfl-bio-suite endpoint smoke
```

## Step 6 — Onboard each partner

Generate their bundle:

```bash
appfl-bio-suite partner-bundle gwas --site site-north
```

That produces a directory containing:

- **`1-endpoint-setup.md`** — how to set up their endpoint. Same for every experiment.
- **`2-<experiment>.md`** — what differs for this experiment, with a verification command
  and its expected output.
- **`your-config-block.yaml`** — their configuration, every value already filled in.
- **`example_identity_mapping_config.json`** — their identity mapping, ready to install.
- **`user_config_template.yaml.j2`** — their endpoint template, with their scheduler,
  account and queue already set.
- **`validate_mapping.py`** — so they can check their own work without waiting on you.

Send them the directory. There are **no placeholders** for them to interpret — that is
deliberate, and it is the single biggest improvement over how this used to be done.
Essentially every question partners asked was about which placeholder applied to them.

They send back their endpoint UUID and one confirmation line. Put the UUID into
`local/federation.yaml`.

> Documents in a bundle are copied from `docs/partner/` and nowhere else. That is enforced
> in code, not by convention, so "did I accidentally include my internal notes" is not a
> judgment you have to make correctly under time pressure.

## Step 7 — Verify each partner before running

```bash
appfl-bio-suite endpoint smoke --experiment gwas
```

This submits a trivial task and reports **which POSIX account it ran as**. That is the
only real proof the identity mapping took effect — an endpoint being "online" tells you a
daemon is connected and nothing more.

If the user in the result is not the experiment's service account, the mapping is not
working. See [troubleshooting.md](troubleshooting.md#identity-mapping) — and note that a
403 and a 422 have completely different causes.

Do this for every site **before** a real run, not during one.

## Step 8 — Generate and distribute data

FLamby partners download a public dataset themselves; skip to Step 9.

For GWAS you generate the data and distribute per-site bundles:

```bash
appfl-bio-suite simulate gwas --scenario three-site-skewed --out local/data/gwas/runs/three-site-skewed
```

Each `local/data/gwas/runs/three-site-skewed/SiteN/data/` directory is one partner's bundle. Send each site
theirs, by whatever transfer mechanism you and they already trust.

Every run writes a manifest recording the scenario, every seed, the suite commit, input
and output checksums, and package versions. Keep it — it is what lets you prove later
which data produced which result.

Read [the GWAS data documentation](../experiments/gwas/DATA.md) before using any output
for a publication. It is explicit about what this pipeline does and does not reproduce.

Fine-mapping additionally writes a DUO profile into each bundle and registers every
bundle with DRS, so a site can verify it holds the bundle you cut for it and can enforce
its own data use terms. Wiring those into `federation.yaml` is optional and additive; see
[ga4gh.md](ga4gh.md).

## Step 9 — Run

```bash
appfl-bio-suite preflight --experiment gwas
appfl-bio-suite run gwas --dry-run     # inspect the resolved configs first
appfl-bio-suite run gwas
```

Run it under `tmux` or `screen`. The driver must outlive the scheduler jobs it triggers on
partner clusters, and those can sit in a queue.

For a first multi-site run, start small. Prove the transport completes before asking for
twenty rounds — with sites in different time zones and different queue policies, a
synchronous run proceeds at the speed of the slowest site, every round.

---

## When something goes wrong

[troubleshooting.md](troubleshooting.md) is organized by symptom and gives the mechanism
rather than just the fix. Most entries there cost someone days.

The three failures worth knowing about before you hit them:

1. **`^` or `$` in an identity mapping.** Looks correct, can never match. Catch it with
   `appfl-bio-suite identity validate`.
2. **An endpoint started by an unprivileged user.** Starts successfully, silently ignores
   its identity mapping. "It started fine" proves nothing.
3. **Version skew.** Does not fail at install time. Fails partway into a run, on someone
   else's machine.

---

## What to read next

| | |
| --- | --- |
| [architecture.md](architecture.md) | How the pieces fit, and why it is built this way |
| [reference-deployment.md](reference-deployment.md) | A complete working deployment, as an example |
| [troubleshooting.md](troubleshooting.md) | Symptom-indexed, with mechanisms |
| [releasing.md](releasing.md) | Upgrading pins without breaking the federation |
| [../experiments/](../experiments/) | Per-experiment design, data, and runbooks |
