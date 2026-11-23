# Unreleased

Fixes from the first live single-user site, a cloud VM running Ubuntu 24.04. Each of
these stopped or misled the onboarding as written in 0.3.2.

- **Activation no longer reports a dead endpoint as started.** Globus Compute launches
  its interchange as bare `globus-compute-endpoint`, so an endpoint started by absolute
  path, without the runtime on `PATH`, registered and then died. `site activate` printed
  its UUID anyway. Activation now puts the runtime first on `PATH`, waits to see the
  endpoint stay up, and fails with the log's location if it does not.
- **The generated systemd unit keeps the endpoint running.** It sets `PATH`, and the
  single-user configuration now sets `detach_endpoint: false`. Before, the endpoint
  daemonized under `Type=simple`, systemd saw a clean exit, and the service went
  inactive with nothing to restart. **Re-run `site activate` and reinstall the unit on a
  0.3.2 single-user site.**
- The printed systemd hand-off stops the background copy as the service account, not
  as root, which could not see the endpoint.
- **Credentials are checked before anything is configured.** Activation asks Globus
  Auth whether the client ID and secret belong together and says which field is the
  client ID (the Client UUID, not the project UUID or a secret's ID). A value still
  wrapped in the guide's `<...>` is refused by name. The guide shows bare values, and
  `join` names the credentials file it expects.
- **`appfl-bio-suite site check-data BUNDLE`** loads the staged data with the task's own
  arguments. The command the GWAS and fine-mapping guides printed passed only a
  directory, so a cohort that declares data-use terms was refused as if the data were
  wrong.
- `site setup` and `join` name the remedy when `python3.12` cannot create an environment
  (Ubuntu 24.04 without `python3.12-venv`), and the guide checks for it first.
- **`texttable` is pinned to 1.7.1.** `globus-compute-endpoint list` crashes with any
  later release. `torchvision`, `wandb` and `contourpy` are pinned so pip stops
  backtracking through them on every install.
- The partner guide's commands use the bundle's real filename, the fine-mapping guide
  covers a `.zip` delivery and states the space needed, and a link to a guide that is
  not in the bundle is gone.
- **The onboarding ZIP carries one setup path.** `experiment package` no longer ships
  Part 1's manual endpoint guide, worker template, identity mapping, validator or
  config block. Part 1's template activated a hand-made conda environment under the same
  filename `site setup` writes for its own runtime, so a partner copying it ran workers
  in an environment nothing had verified. The experiment guide is now `1-<experiment>.md`
  and describes the `join` setup; FLamby is installed into the prepared runtime.
  `partner-bundle` still writes the manual bundle.
- **Sites can be added to an experiment at any time.** A bundle's experiment digest
  covers the shared settings and its own site entry, not the other sites or
  `secure_aggregation` (which lists them). Adding, removing or registering another site
  no longer stops a partner who has not registered yet from registering.
- `bundle_manifest.json` records each file's SHA-256 and the bundle digest a receipt
  carries, so a ZIP can be matched to a fresh packaging and to a partner's
  `registration.json`.
- `experiment register` accepts a site that names its own `service_account`. It
  compared the receipt with the experiment's default account, which the bundle had
  never recorded for that site, so every such site was refused.

# appfl-bio-suite 0.3.2

Cross-site fine-mapping over Globus Compute, and onboarding for sites with no root and
no batch scheduler. Found while preparing a cloud-VM site; none of it had been exercised
end to end before.

- **Fine-mapping results now come back.** Globus Compute endpoints refuse task results
  over 10 MiB, and one fine-mapping locus already exceeds that, so 0.3.1 could not
  complete a cross-site run. A site whose payload is over 6 MiB now keeps it in
  `output_dir/spool/` and returns a small handle. The Globus Compute driver reads it back
  through the same endpoint and credentials in 6 MiB chunks, verifies the sha-256 the
  site reported, deletes the site's copy, and aggregates the payload unchanged. `run`
  enables this for the Globus Compute driver only; serial and TES runs keep the inline
  payload. **This changes the shipped trainer.**
- **Single-user endpoints through the bundle.** A site entry may declare
  `endpoint_mode: single-user`. Its bundle tells `join` so, and `site activate` then runs
  as the service account (never root). It takes the site's Globus service-account client
  from a 0600 file, starts a single-user endpoint, prints its UUID, and writes a systemd
  unit. `experiment register --site-credentials` probes it as that client.
- **`scheduler: local`** runs workers on the endpoint host itself, for cloud VMs and
  workstations. No account or queue is asked for.
- **Mixed credential routes.** A routing manifest entry may say
  `"coordinator_login": true` for one endpoint, so a site you reach with your own login
  can share a run with service-account sites. An endpoint with no route is still refused.
- Bundles whose specification exceeds 1 MB (any DRS-verified fine-mapping cohort of
  realistic size) are read from their ZIP; the limit is now 64 MB.
- `site setup` no longer requires the data to be staged first; activation and
  verification still check it.
- A missing scheduler setting is named, rather than reported as a pattern mismatch.
- `bundle_manifest.json` no longer records an unrelated repository's commit when the
  suite is installed from PyPI into a virtual environment inside one.

## Upgrading from 0.3.1

Every site in a federation runs the same release, and the shipped trainer changed.
Upgrade partners first and the coordinator last (`docs/coordinator/releasing.md`), then
regenerate partner bundles: a bundle cut with 0.3.1 refuses to set up under 0.3.2. The
experiment digest now includes each site's `endpoint_mode`, so registrations must use
bundles generated by 0.3.2. Federation files need no change: `endpoint_mode` defaults to
`multi-user` and `scheduler` to `slurm`, as before.

## Installation

Python 3.12 is required:

```bash
python -m pip install "appfl-bio-suite[all]==0.3.2" \
    -c https://raw.githubusercontent.com/mfjoel01/appfl-bio-suite/v0.3.2/constraints.txt
```

Coordinator scripts, example configurations and `federation.yaml.example` are not in the
wheel. Take them from the `v0.3.2` tag.

The MultiSuSiE backend is pinned to a validated upstream revision, which PyPI metadata
cannot carry, so it is installed alongside the extra:

```bash
python -m pip install "appfl-bio-suite[multisusie]==0.3.2" \
    -r https://raw.githubusercontent.com/mfjoel01/appfl-bio-suite/v0.3.2/requirements-multisusie.txt \
    -c https://raw.githubusercontent.com/mfjoel01/appfl-bio-suite/v0.3.2/constraints.txt
```

FLamby and the PLINK/SuSiEx binaries keep their documented separate installation
procedures.

# appfl-bio-suite 0.3.1

A packaging and documentation release. The modules shipped to partner workers and all
dependency pins are unchanged from 0.3.0.

- Partner bundles install the suite from PyPI at this exact version,
  `appfl-bio-suite[<extras>]==0.3.1`, instead of from the Git tag. A bundle built with
  `--wheel` still installs the wheel it carries.
- `preflight --check pins` works from any directory. With no `constraints.txt` in the
  working directory or a checkout, it checks against the copy packaged in the wheel.
- When the MultiSuSiE or network-map extra is missing, the error gives the PyPI install
  command with this release's constraints.
- The README, coordinator guides, experiment runbooks and website install from PyPI. A
  test fails if a documented install command names any other release.

## Upgrading from 0.3.0

Every site in a federation runs the same release. Upgrade partners first and the
coordinator last (`docs/coordinator/releasing.md`), then regenerate partner bundles: a
bundle cut with 0.3.0 refuses to set up under 0.3.1.

## Installation

Python 3.12 is required:

```bash
python -m pip install "appfl-bio-suite[all]==0.3.1" \
    -c https://raw.githubusercontent.com/mfjoel01/appfl-bio-suite/v0.3.1/constraints.txt
```

Coordinator scripts, example configurations and `federation.yaml.example` are not in the
wheel. Take them from the `v0.3.1` tag.

The MultiSuSiE backend is pinned to a validated upstream revision, which PyPI metadata
cannot carry, so it is installed alongside the extra:

```bash
python -m pip install "appfl-bio-suite[multisusie]==0.3.1" \
    -r https://raw.githubusercontent.com/mfjoel01/appfl-bio-suite/v0.3.1/requirements-multisusie.txt \
    -c https://raw.githubusercontent.com/mfjoel01/appfl-bio-suite/v0.3.1/constraints.txt
```

FLamby and the PLINK/SuSiEx binaries keep their documented separate installation
procedures.

# appfl-bio-suite 0.3.0

This release adds portable partner onboarding, a MultiSuSiE inference backend, and
fixed-roster secure aggregation for fine-mapping. It is also the suite's first release
on [PyPI](https://pypi.org/project/appfl-bio-suite/0.3.0/), with the same version as
the Git release.

- Generate portable partner bundles and use guided setup and readiness checks.
- Configure service-account endpoints for independently administered sites, with per-site service account overrides.
- Select SuSiEx or MultiSuSiE inference over the same APPFL fine-mapping payloads.
- Use fixed-roster secure pooling with both fine-mapping backends.
- Verify GA4GH DRS nested bundle checksums and handle resilient TES execution with concurrent site cancellation.
- Reproduce GWAS polygenic scores independently of the host's BLAS reduction order.
- Install versioned wheels with packaged partner guides, templates, and consent data.

DNA methylation remains a single-machine experimental workflow. Local and CI tests
do not establish live cross-site endpoint readiness; run the documented partner
readiness checks before a federation.
