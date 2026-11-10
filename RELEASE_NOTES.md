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
