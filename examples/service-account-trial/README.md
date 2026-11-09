# Single-user service-account trial setup

These files support the [Hyperstack–Argonne trial guide](../../docs/coordinator/service-account-trial.md).
This opt-in implementation still needs live cross-site and scheduler validation.

Prepare Python 3.12 on the coordinator and each worker. Clone the release tag for these
example files and constraints, then install the same release from PyPI:

```bash
git clone --branch v0.3.0 https://github.com/mfjoel01/appfl-bio-suite.git
cd appfl-bio-suite
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install "appfl-bio-suite==0.3.0" -c constraints.txt
python -m pip check
```

For the synthetic CPU pilot, replace the install command above with:

```bash
python -m pip install "appfl-bio-suite==0.3.0" -c constraints.txt \
  -c examples/service-account-trial/constraints-cpu.txt \
  --extra-index-url https://download.pytorch.org/whl/cpu
```

The CPU constraints supplement the production constraints. The standard production
preflight requires the exact production Torch version string; use the trial guide's
direct synthetic driver for this CPU pilot. Dependencies require internet access.

On the coordinator only, copy the routing example:

```bash
mkdir -p local
cp examples/service-account-trial/site-credentials.example.json local/site-credentials.json
```

Replace both endpoint UUIDs, client IDs, and expected Unix usernames. Supply the two
named secret environment variables through your protected secrets mechanism. Never
put secret values in this file or distribute the coordinator manifest to partners.
Each site receives credentials only for its own service-account identity.

Follow the trial guide to configure providers, start endpoints, verify isolation, and
run two synthetic rounds. Argonne requires a suitable scheduler allocation and worker
environment; these examples do not choose those site-specific settings.
