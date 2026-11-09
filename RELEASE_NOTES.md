# appfl-bio-suite 0.3.0

This release adds portable partner onboarding, a MultiSuSiE inference backend, and
fixed-roster secure aggregation for fine-mapping. It also prepares the suite's first
PyPI distribution with the same version as the Git release.

- Generate portable partner bundles and use guided setup and readiness checks.
- Configure service-account endpoints for independently administered sites, with per-site service account overrides.
- Select SuSiEx or MultiSuSiE inference over the same APPFL fine-mapping payloads.
- Use fixed-roster secure pooling with both fine-mapping backends.
- Verify GA4GH DRS nested bundle checksums and handle resilient TES execution with concurrent site cancellation.
- Reproduce GWAS polygenic scores independently of the host's BLAS reduction order.
- Install versioned wheels with packaged partner guides, templates, and consent data.

## Installation

Python 3.12 is required:

```bash
python -m pip install "appfl-bio-suite[all]==0.3.0"
```

Use the `v0.3.0` checkout for coordinator scripts, example configurations, and
`constraints.txt`. Generated partner bundles continue to reference that immutable tag.
The Globus Compute and scientific dependency pins are unchanged.

MultiSuSiE users must install the backend separately from the matching checkout:

```bash
python -m pip install ".[multisusie]" -r requirements-multisusie.txt -c constraints.txt
```

The backend remains pinned to its validated upstream revision. PyPI distribution
metadata cannot contain a direct Git dependency. FLamby and the PLINK/SuSiEx binaries
also retain their documented separate installation procedures.

DNA methylation remains a single-machine experimental workflow. Local and CI tests
do not establish live cross-site endpoint readiness; run the documented partner
readiness checks before a federation.
