# Bundle-based onboarding

Prepared for **{{ site_name }}** ({{ client_id }}), experiment **{{ experiment }}**.
Coordinator identity: **{{ coordinator_identity }}**.

The coordinator supplies one ZIP per partner. It contains `onboarding.json`, local
settings in `site.yaml`, and the experiment's data instructions. The bundle identifies
only your participation; keep your edited `site.yaml` separate from the original ZIP.
The onboarding commands require the release that generated the bundle.

**Your endpoint:** {% if single_user %}a **single-user** endpoint that runs as
`{{ service_account }}` itself and authenticates as a Globus service-account client your
coordinator issues to your site. No root access and no identity mapping are needed.
{% else %}a **multi-user** endpoint (MEP), which an administrator with root starts and
which maps the coordinator's identity onto `{{ service_account }}`.
{% endif %}**Workers run:** {% if local_scheduler %}on this host itself, with no batch
scheduler, so there is no scheduler account or queue to choose.{% else %}as
{{ scheduler }} jobs.{% endif %} Both choices are fixed by the bundle.

## Install the setup tool

Unpack the ZIP and install the tool in a Python 3.12 environment. The suite is published
on [PyPI](https://pypi.org/project/appfl-bio-suite/); install exactly the release that
generated your bundle:

```bash
python3.12 -m zipfile -e site.zip site-bundle
python3.12 -m venv bootstrap
bootstrap/bin/python -m pip install --extra-index-url https://download.pytorch.org/whl/cpu \
  -c site-bundle/runtime-constraints.txt "{{ install_spec }}"
source bootstrap/bin/activate
```

No `python3.12`? Ubuntu 22.04, for one, ships 3.10. Install 3.12 alongside it (for
example from the deadsnakes PPA, or with `uv python install 3.12`) and use that
interpreter for these commands and for `--python` below.

If your coordinator included a wheel in the ZIP, install that file instead by replacing
`"{{ install_spec }}"` with `site-bundle/*.whl`. The bootstrap environment supplies the
CLI; setup creates the runtime used by the workers.

## Join in one guided command

{% if single_user %}Run it as `{{ service_account }}`, the account the endpoint will run
as. Create that account first if it does not exist; it needs no special privileges.
{% endif %}After installing BIOSIM, run:

```bash
appfl-bio-suite join site.zip
```

The command asks for {% if not local_scheduler %}the scheduler account, {{ scheduler_slot }},
and {% endif %}absolute data/work paths, showing coordinator-provided values as
suggestions. Invalid answers are reported at the question. It prepares the runtime and
prints the exact activation and registration commands for this installation. It never
activates an endpoint as a side effect of join.

You can also prepare your settings **before** pointing BIOSIM at the bundle:

```yaml
{% if not local_scheduler %}account: my-allocation
{{ scheduler_slot }}: cpu
{% endif %}data_dir: /shared/project/data
work_dir: /shared/project/work
```

Then run:

```bash
appfl-bio-suite join site.zip --local site.yaml
```

Values supplied in the file are used without asking again. Missing values are prompted.
Use `--non-interactive` for automation: bundle defaults fill omitted fields and invalid
or missing required values fail before installation. Optional `interface` belongs in
the local file when the cluster needs it.

Scheduler settings are remembered in the application config directory's `cluster.yaml`.
Use `--profile /path/to/cluster.yaml` to choose a profile for a different cluster, or
`--no-profile` to disable reading/saving profiles. Profiles never retain experiment
identities, credentials, endpoint UUIDs, or data/work paths. Explicit local values take
precedence over a profile. A profile for another scheduler is refused.

The installation defaults to `WORK_DIR/.biosim/BUNDLE_REVISION`; choose `--state` to
place it elsewhere{% if not local_scheduler %} on a shared filesystem{% endif %}. Resolved
local settings are saved there as `site.yaml`. To resume, supply the same bundle and
`--state`; already recorded settings are not prompted again. Changing the bundle or
local settings requires a new directory.

The data can arrive before or after `join`; it must be in place before activation.
{% if single_user %}
Follow the activation command printed by join (see [Activate the endpoint](#activate-the-endpoint)
below). It prints the endpoint UUID; then run:
{% else %}
Follow the administrator activation command printed by join. Once the administrator
supplies the endpoint UUID, run as the service account:
{% endif %}
```bash
appfl-bio-suite join site.zip --state /shared/project/onboarding --endpoint-id UUID
```

Use your actual state directory from the first command. A prepared installation resumes
verification without reinstalling packages. Send its `registration.json` to the
coordinator for the remote probe. No registration is sent automatically.

## Individual commands and local settings

```bash
appfl-bio-suite site init site.zip --out site.yaml
```

Edit {% if not local_scheduler %}the scheduler account, {{ scheduler_slot }}, {% endif %}the
data directory and work directory. Use absolute paths visible from {% if local_scheduler %}this
host{% else %}compute nodes{% endif %}. An optional `interface` selects the endpoint's interface
facing the workers. Unknown keys and missing settings are errors. Experiment settings,
coordinator identity and software versions come from the bundle and cannot be
overridden in this file.

Stage your experiment data using the experiment guide in the ZIP before activation.
The onboarding ZIP contains software configuration, not the scientific datasets.

## Prepare the runtime

Run setup as the experiment's service account, or an administrator preparing an
installation accessible to that account. Choose a state directory {% if not local_scheduler %}shared with compute
nodes and {% endif %}writable by the service account. Python 3.12 must already be available.

```bash
appfl-bio-suite site setup site.zip --local site.yaml \
  --state /shared/project/onboarding --python /path/to/python3.12
```

Setup creates an isolated Python environment, installs the experiment's dependencies
under the bundle's constraints, checks compatibility and work-directory access, and
renders the endpoint configuration. It can be rerun after an interrupted installation.
Changing inputs requires a new state directory. Nothing edits an unrelated Python
environment.

Without a bundled wheel, setup installs `{{ install_spec }}` from PyPI. If the
coordinator included a wheel, setup verifies its checksum and installs it with the
experiment extras automatically. No source override is needed in either case.

To install a separately delivered wheel instead, supply
`--install-source '/absolute/path/appfl_bio_suite-VERSION-py3-none-any.whl[gwas]'`
(substitute the extra for your experiment). The coordinator must distribute the matching
wheel.

## Activate the endpoint
{% if single_user %}
Your coordinator gives you a Globus **service-account client ID and secret** for this
site. Never paste the secret into a chat, a ticket or a shared file. As
`{{ service_account }}`, put it in a file only that account can read:

```bash
mkdir -p ~/.config && install -m 600 /dev/null ~/.config/{{ endpoint_name }}.env
# then edit it to contain exactly these two lines:
#   GLOBUS_COMPUTE_CLIENT_ID=<client UUID from your coordinator>
#   GLOBUS_COMPUTE_CLIENT_SECRET=<secret from your coordinator>
```

Then, as `{{ service_account }}` and using the prepared runtime's executable:

```bash
/shared/project/onboarding/runtime/bin/appfl-bio-suite site activate \
  --state /shared/project/onboarding --credentials-file ~/.config/{{ endpoint_name }}.env
```

Activation refuses to run as root or as any account other than `{{ service_account }}`,
because the endpoint runs every task as whoever starts it. It checks the runtime, the
data and the work directory as that account, configures a dedicated single-user endpoint
named `{{ endpoint_name }}`, starts it in the background, and prints its **UUID**. It
refuses to alter an endpoint it did not create. The credentials are read from the file,
never echoed, and never placed in any configuration sent to anyone.

It also writes `{{ endpoint_name }}.service` in the state directory: a systemd unit that
keeps the endpoint running across logouts, crashes and reboots. Installing it needs root
once; the command prints the three lines to run. Without it, the endpoint stops when the
host reboots and must be activated again.

Then, as `{{ service_account }}`:

```bash
appfl-bio-suite site verify --state /shared/project/onboarding --endpoint-id UUID
```

A single-user endpoint accepts tasks only from its own service-account client, so you
cannot submit a test task to it from your personal Globus login. That is the isolation
working, not a fault. The coordinator holds a credential for the same client and
performs the remote probe.
{% else %}
An administrator must create the named service account and give it scheduler access.
The account must be able to traverse the state directory, activate its runtime, read the
data, and write the work directory. Setup does not create POSIX accounts or allocations.

As root, using the prepared runtime's executable:

```bash
/shared/project/onboarding/runtime/bin/appfl-bio-suite site activate \
  --state /shared/project/onboarding
```

Activation verifies runtime access as the service account, creates a dedicated multi-user
endpoint, installs the mapping and worker template, and starts it. It refuses to alter an
endpoint it does not own. Existing institutional endpoints should continue to use the
manual integration instructions in Part 1; this command manages dedicated endpoints.
Globus authentication may require an interactive login.

As the same administrator, obtain the UUID:

```bash
/shared/project/onboarding/runtime/bin/globus-compute-endpoint list
```

As the service account:

```bash
appfl-bio-suite site verify --state /shared/project/onboarding --endpoint-id UUID
```
{% endif %}
Send the resulting `registration.json` to the coordinator. It records local verification;
the coordinator then verifies a real remote worker before accepting the registration.
Use `--runtime-only` to check local runtime access without producing a receipt.

## What readiness means

The remote probe checks the actual POSIX account, Python and pinned runtime versions,
data-directory access, and the ability to create a temporary file in the work directory.
It does not read or transmit scientific records. A successful probe establishes that the
endpoint, scheduler and runtime work together. Experiment-specific data validation and
consent checks still run through the existing experiment/preflight workflow. FLamby's
external dataset/package preparation remains in its experiment guide.

Containers and fully offline dependency installation are not part of this initial
workflow. Python packages are fetched during setup; cluster modules and scheduler
permissions remain local administrator responsibilities.
