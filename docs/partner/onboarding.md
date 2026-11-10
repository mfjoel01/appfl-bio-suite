# Bundle-based onboarding

Prepared for **{{ site_name }}** ({{ client_id }}), experiment **{{ experiment }}**.
Coordinator identity: **{{ coordinator_identity }}**.

The coordinator supplies one ZIP per partner. It contains `onboarding.json`, local
settings in `site.yaml`, and the experiment's data instructions. The bundle identifies
only your participation; keep your edited `site.yaml` separate from the original ZIP.
The onboarding commands require the release that generated the bundle.

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

If your coordinator included a wheel in the ZIP, install that file instead by replacing
`"{{ install_spec }}"` with `site-bundle/*.whl`. The bootstrap environment supplies the
CLI; setup creates the runtime used by the workers.

## Join in one guided command

After installing BIOSIM, run:

```bash
appfl-bio-suite join site.zip
```

The command asks for the scheduler account, queue/partition, and absolute data/work
paths, showing coordinator-provided values as suggestions. Invalid answers are reported
at the question. It prepares the runtime and prints the exact activation and registration
commands for this installation. It never activates an endpoint as a side effect of join.

You can also prepare your settings **before** pointing BIOSIM at the bundle:

```yaml
account: my-allocation
partition: cpu
data_dir: /shared/project/data
work_dir: /shared/project/work
```

For PBS Pro, use `queue` instead of `partition`. Then run:

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
place it elsewhere on a shared filesystem. Resolved local settings are saved there as
`site.yaml`. To resume, supply the same bundle and `--state`; already recorded settings
are not prompted again. Changing the bundle or local settings requires a new directory.

**Current endpoint support:** this workflow supports the multi-user endpoint (MEP).
The single-user Globus pattern is actively being tested separately and is not yet
available through `join`. The endpoint handoff is kept separate so that pattern can be
incorporated once validated; it is never silently inferred from the current user.

For MEP, follow the administrator activation command printed by join. Once the
administrator supplies the endpoint UUID, run as the service account:

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

Edit the scheduler account, partition (Slurm) or queue (PBS Pro), data directory and
work directory. Use absolute paths visible from compute nodes. An optional `interface`
selects the endpoint's interface facing those nodes. Unknown keys and missing settings
are errors. Experiment settings, coordinator identity and software versions come from
the bundle and cannot be overridden in this file.

Stage your experiment data using the experiment guide in the ZIP before verification.
The onboarding ZIP contains software configuration, not the scientific datasets.

## Prepare the runtime

Run setup as the experiment's service account, or an administrator preparing an
installation accessible to that account. Choose a state directory shared with compute
nodes and writable by the service account. Python 3.12 must already be available.

```bash
appfl-bio-suite site setup site.zip --local site.yaml \
  --state /shared/project/onboarding --python /path/to/python3.12
```

Setup creates an isolated Python environment, installs the experiment's dependencies
under the bundle's constraints, checks compatibility and local path access, and renders
the endpoint configuration. It can be rerun after an interrupted installation. Changing
inputs requires a new state directory. Nothing edits an unrelated Python environment.

Without a bundled wheel, setup installs `{{ install_spec }}` from PyPI. If the
coordinator included a wheel, setup verifies its checksum and installs it with the
experiment extras automatically. No source override is needed in either case.

To install a separately delivered wheel instead, supply
`--install-source '/absolute/path/appfl_bio_suite-VERSION-py3-none-any.whl[gwas]'`
(substitute the extra for your experiment). The coordinator must distribute the matching
wheel.

## Activate the endpoint

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
