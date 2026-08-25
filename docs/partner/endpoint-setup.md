# Part 1 — Setting up your endpoint

You are being asked to set up a **multi-user Globus Compute endpoint (MEP)** on your
cluster. When a task arrives from {{ coordinator_organization }}, your cluster runs it
**locally** and returns only results — nothing else leaves your machine. You authorize
one Globus identity; no secret is shared, and you never run the experiment yourself.

**This is an admin task.** You need privileged (root) access on a login node, a scheduler
account and {{ scheduler_slot }}, and `conda`.

Everything you need is in this document. Nothing is sent separately.

> **Already run another Globus Compute MEP on this cluster?** Your account, cluster and
> conda work carries over, but do **not** reuse that endpoint. Build a separate one
> (`{{ endpoint_name }}` with service account `{{ service_account }}`) so the two stay
> independent, with separate scheduler accounting.

This is Part 1 of 2 and is the same for every experiment. Part 2 covers what differs for
yours.

---

## Your details

| | |
| --- | --- |
| Your institution | {{ site_name }} |
| Your site ID | `{{ client_id }}` |
| Endpoint to create | `{{ endpoint_name }}` |
| Service account to create | `{{ service_account }}` |
| Scheduler | {{ scheduler }} |
| Identity to authorize | `{{ coordinator_identity }}` |

---

## Step 1 — Get the software

```bash
python -m venv ~/appfl-env       # or: conda create -n appfl_env python=3.12 -y
source ~/appfl-env/bin/activate  # or: conda activate appfl_env
python -V                        # must print 3.12.x
```

Python must be **{{ python_version }}**. We ship function bytecode between machines, and
bytecode is only portable within a minor version — a mismatch fails before any of your
code runs. The patch level does not need to match ours.

```bash
pip install "appfl-bio-suite[{{ partner_extras }}]"
```

Verify:

```bash
globus-compute-endpoint version
```

It **must** print `{{ globus_compute_version }}`. If it prints something else, stop and
tell us. Every site in the federation runs the same version; a mismatch surfaces as a
deserialization failure several rounds into a run, which is the most expensive kind of
bug for everyone involved.

## Step 2 — Check that multi-user mode works

```bash
globus-compute-endpoint configure --multi-user true _probe
globus-compute-endpoint delete --yes _probe
```

Both lines must finish with **no error**. If either errors, stop and send us the output
before continuing — we will align the whole federation on one version rather than have
you work around it.

## Step 3 — Create the service account

Create a local POSIX account named **`{{ service_account }}`**. It must be able to:

- activate the environment from Step 1
- submit {{ scheduler }} jobs
- read and write `{{ output_dir }}`

Every task we send runs as this account and nothing else.

## Step 4 — Configure the endpoint

> If your centre already operates a Globus Compute MEP, do not build your own — send them
> this whole section and ask them to add it.

Run as the **privileged (root)** user:

```bash
globus-compute-endpoint configure --multi-user true {{ endpoint_name }}
```

That creates `~/.globus_compute/{{ endpoint_name }}/` with several files already filled
in. You only need to change two of them.

### 4a — `config.yaml` (usually no change)

`configure` already wrote this and already pointed it at the identity-mapping file it
created alongside. **Leave that path exactly as generated** — do not retype it. If
`display_name` came out as `null`, set it to `{{ endpoint_name }}`. That is the only edit.

### 4b — the identity mapping

Replace the entire contents of the file that `identity_mapping_config_path` names —
`example_identity_mapping_config.json`, in the same folder — with the block below. We have
included it pre-filled as `example_identity_mapping_config.json` in this bundle.

```json
[
  {
    "DATA_TYPE": "expression_identity_mapping#1.0.0",
    "mappings": [
      {
        "source": "{username}",
        "match": "{{ coordinator_identity_match }}",
        "output": "{{ service_account }}"
      }
    ]
  }
]
```

This says: *if the submitter is `{{ coordinator_identity }}`, run the task as
`{{ service_account }}`.* That mapping is the entire gatekeeper. No other identity is
accepted.

> ### Do not add `^` or `$` anchors to `match`
>
> This is the single most expensive mistake on this project — it cost one site over a
> week, and a wrong diagnosis was sent to another before it was understood.
>
> `match` is **not** a normal regular expression. The Globus expression mapper *escapes*
> any `^` or `$` you write, turning them into literal characters, and then wraps your
> pattern in its own `^...$`. So writing `"^{{ coordinator_identity_match }}$"` compiles
> to a pattern that matches only a literal string beginning with a caret — it can never
> match anything, and every submission is rejected with
> `422 SEMANTICALLY_INVALID … Identity failed to map to a local user name`.
>
> The unanchored form above is still an exact whole-string match, because the mapper adds
> the anchors for you. Omitting them is not a loosening.
>
> Keep the backslash before each `.` — it escapes the dot so it matches only a literal
> dot. Do **not** escape anything else: the mapper accepts only a small set of escapes
> and rejects the whole file if it sees one it does not recognize.

The endpoint **polls this file** and picks up changes within about five seconds. A restart
is not required after editing it.

### 4c — the worker template

Edit `user_config_template.yaml.j2`. We have included a pre-filled copy in this bundle;
check the account and {{ scheduler_slot }} are right for you and use it as-is.

```yaml
engine:
  type: GlobusComputeEngine
  max_workers_per_node: 1
  strategy: simple
  provider:
    type: {{ provider_type }}
    account: {{ account }}
    {{ scheduler_slot }}: {{ scheduler_slot_value }}
    nodes_per_block: 1
    init_blocks: 0
    min_blocks: 0
    max_blocks: 1
    walltime: 00:30:00
    launcher:
      type: {{ launcher_type }}
    worker_init: |
      module load {{ conda_module }}
      source "$(conda info --base)/etc/profile.d/conda.sh"
      conda activate appfl_env
      export TMPDIR=/tmp
      export OPENBLAS_NUM_THREADS=1
      export OMP_NUM_THREADS=1
      export MKL_NUM_THREADS=1
idle_heartbeats_soft: 0
idle_heartbeats_hard: 5760
```

Every line in `worker_init` prevents a specific failure we have already hit. Please keep
all of them:

- **`export TMPDIR=/tmp`** — schedulers often set `$TMPDIR` to a long path. The worker
  pool opens a Unix socket underneath it, and the socket path limit is about 108
  characters. Overflow it and workers never register: you see
  *"workers failed to register"* or *"AF_UNIX path too long"*.
- **the three `*_NUM_THREADS=1` lines** — numpy's linear-algebra backend starts one
  thread per core on import. On a busy node that exhausts the per-user thread budget and
  you get `RuntimeError: can't start new thread`.
- **`idle_heartbeats_soft: 0` / `idle_heartbeats_hard: 5760`** — these keep your worker
  block warm between tasks. Without them the block is torn down and re-queued, which on a
  busy queue turns a run that should take minutes into one that takes hours.

> **Workers start but never connect back?** The engine could not auto-detect which network
> interface faces your compute nodes. Add this under `engine:`:
>
> ```yaml
>   address:
>     type: address_by_interface
>     ifname: <the interface your compute nodes use>
> ```

## Step 5 — Verify the mapping before you start

Reading the JSON by eye does not catch the anchor problem, because a wrong file looks
correct. Check it mechanically, in the same environment as the endpoint:

```bash
appfl-bio-suite identity validate \
    /root/.globus_compute/{{ endpoint_name }}/example_identity_mapping_config.json \
    --identity '{{ coordinator_identity }}' \
    --expect '{{ service_account }}'
```

It prints what your `match` actually compiles to and says PASS or FAIL in about a second,
with no round trip to us. Get **PASS** before continuing. If it fails, no amount of
restarting will help — the message names the fix.

## Step 6 — Start the endpoint

Run as the **same privileged user** that ran Step 4:

```bash
globus-compute-endpoint start {{ endpoint_name }}
globus-compute-endpoint list      # copy the Endpoint ID (UUID)
```

> **Two traps here, each of which has cost a site a week.**
>
> **1. Starting as an ordinary user silently disables identity mapping.** The endpoint
> still starts, with no error and nothing that looks wrong — but it ignores
> `identity_mapping_config_path` entirely and accepts only the identity that started it.
> Every submission from us is then rejected with `403`. "It started fine" is not evidence
> that Step 4b took effect.
>
> **2. Endpoint UUIDs are per-account.** The ID lives in
> `<that user's home>/.globus_compute/{{ endpoint_name }}/endpoint.json`. An endpoint you
> built under your personal account and one built under `root` are *different endpoints
> with different UUIDs*. If you configured as yourself and later restarted as root, the
> UUID you send us must be the one `globus-compute-endpoint list` prints **as root**.

## Step 7 — Send us

1. Your **endpoint UUID** from Step 6.
2. One line confirming `{{ coordinator_identity }}` is mapped to `{{ service_account }}`.
3. Confirmation that `globus-compute-endpoint version` still prints
   `{{ globus_compute_version }}`.

Then **leave the endpoint running**. Part 2 covers the data step for your experiment, and
after that we drive everything from our side.

---

## Troubleshooting

**`422 SEMANTICALLY_INVALID … Identity failed to map to a local user name`**
Good news: your endpoint is privileged and its mapping *is* being read. The expression
just did not match us. Almost always the `^`/`$` anchors — see the callout in Step 4b.
Run the Step 5 validator; it answers this in a second. Next most likely: you edited a
different file from the one `identity_mapping_config_path` names, which is easy to do if
that path still points into a personal home directory after you switched to root. No
restart is needed after a fix.

**`403 ENDPOINT_ACCESS_FORBIDDEN`**
This is *not* a mapping-expression problem — do not start editing the JSON. It almost
always means the endpoint was started by an ordinary user, so identity mapping was skipped
(see Step 6). Confirm who owns the running process, restart as the privileged user, then
**re-read the UUID** — it will have changed.

**`404`**
Wrong UUID, or the endpoint was deleted. Note that `delete` followed by `configure` issues
a **new** UUID; `stop`/`start`/`restart` preserve it.

**Version is not {{ globus_compute_version }}**
You installed without the pin. Reinstall and tell us — every site must match.

**"workers failed to register" / "AF_UNIX path too long"**
`export TMPDIR=/tmp` is missing from `worker_init`.

**`RuntimeError: can't start new thread`**
The three `*_NUM_THREADS=1` lines are missing from `worker_init`.

**First task works, later ones stall**
Your block is being torn down between tasks and re-queued. Check
`idle_heartbeats_soft: 0` and `idle_heartbeats_hard: 5760` are present, and prefer a
{{ scheduler_slot }} with short queue waits.

**Nothing happens at all**
The queue is busy and the job is waiting in line. That is normal; tell us if it persists.

Anything else — send us `~/.globus_compute/{{ endpoint_name }}/*.log` and we will take it
from there. Contact: {{ coordinator_contact }}
