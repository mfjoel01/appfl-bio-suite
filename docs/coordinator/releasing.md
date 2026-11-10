# Upgrading and releasing

## The rule

**A version bump is a deliberate, tested change — never a routine one.**

Every site in a federation must run byte-identical versions of the Globus Compute stack.
Skew does not fail at install time on your laptop. It fails as a deserialization error
several rounds into a run on somebody else's cluster, in a different country, hours after
you launched. That is the most expensive class of bug this project has.

This is why the suite pins exactly and ships `constraints.txt`, and it is the honest cost
of depending on `appfl` rather than forking it. The benefit was never "automatically stays
current" — that would reintroduce the exact problem. The benefit is that an upgrade is a
one-line change you chose, instead of a merge you inherited.

## Upgrading a pin

1. Change the one line in `constraints.txt`, and in `pyproject.toml` if the pin appears
   there too.
2. Reinstall cleanly: `pip install -e ".[all]" -c constraints.txt`
3. `appfl-bio-suite preflight --check pins`
4. Run the full test suite. The GWAS simulation parity tests are the sensitive ones: they
   assert byte-identical output against the pre-migration pipeline, so a change in numpy,
   pandas or pandas-plink that alters a low-order bit will show up there.
5. Run a loopback federation for both experiments.
6. **Then** ask every partner to upgrade, and re-run `endpoint smoke` at each site before
   the next real run.

Step 6 is the one people skip. Do not skip it.

## Order of operations across a federation

Upgrade the coordinator **last**, not first. A partner on the new version and a
coordinator on the old one is a state you can fix by finishing the rollout; a coordinator
on the new version and partners on the old one means nothing runs until every partner has
acted, and partners act on their own schedule.

## Packages whose version can change results

Recorded in every simulation run manifest, because a number that stops reproducing is
usually one of these having moved:

`numpy`, `pandas`, `scipy`, `scikit-learn`, `pandas-plink`, `xarray`

If a parity test fails after an upgrade, compare the manifest's `packages` block against
the current environment before assuming the code is wrong.

### And one thing that is not a package

**The machine.** A BLAS matrix-vector product adds its terms in an order the library
chooses, and OpenBLAS chooses by CPU model and thread count, each rounding the low-order
bits differently. The GWAS polygenic score used to be one, so the scores and every
phenotype and summary file computed from them followed the hardware. CI passed or failed
depending on which host it landed on, even after BLAS was pinned to one thread.

`_compute_pgs` now computes the product itself, in the reference-BLAS order, without
calling BLAS. The golden checksums are the legacy scripts' output with OpenBLAS pinned to
the kernel that adds in that order (`OPENBLAS_CORETYPE=Prescott`, set by
`tests/legacy_harness.py`), so they hold on every machine.
`test_pgs_is_the_reference_blas_product` fails if the product goes back through BLAS. A
parity failure that comes and goes between re-runs means something follows the machine
again: find it, and do not regenerate the golden file.

## Removing the compat shim

`src/appfl_bio_suite/core/compat.py` works around `appfl==1.10.0` importing a
`globus-compute-sdk` module that 4.9.0 removed. It is meant to be temporary.

`tests/test_upstream_shim.py::test_shim_is_still_necessary` **fails once the workaround is
no longer needed**. That failure is the signal, not a regression. When it fires:

1. Delete `src/appfl_bio_suite/core/compat.py`
2. Delete its call in `src/appfl_bio_suite/__init__.py`
3. Delete `tests/test_upstream_shim.py`
4. Note the removal here

Until then, be aware that a clean install of the pinned set is *functional* only because
of that shim. It is a real constraint, not a tidiness issue.

## Cutting a release

Releases are published on [PyPI](https://pypi.org/project/appfl-bio-suite/), and
generated partner bundles install one exact version from there,
`appfl-bio-suite[<extras>]==X.Y.Z` (`docs/partner/endpoint-setup.md`, Step 1). PyPI never
lets a published version's files be replaced, so every site that installs `==X.Y.Z` runs
the same wheel. `main` moves; a released version does not. Two sites that installed from
`main` a week apart are running different code with no way to say which, which is the
same failure mode `constraints.txt` exists to prevent — just one level up.

Each release also has a matching `vX.Y.Z` git tag. The publish workflow builds the PyPI
distributions from it, and coordinators check it out for the scripts and example files
that are not in the wheel.

Three things carry the version and must move together:

| | |
| --- | --- |
| `pyproject.toml` | `version` |
| `src/appfl_bio_suite/__init__.py` | `__version__` — `install_spec()` and `INSTALL_TAG` are derived from it |
| `CITATION.cff` | `version` |

Install commands in `README.md`, `docs/`, `examples/` and `website/` name the version
too. `tests/test_version_is_consistent.py` fails until every one of them is updated.
`RELEASE_NOTES.md` describes the release being cut.

Then tag and push:

```bash
git tag -a vX.Y.Z -m "appfl-bio-suite X.Y.Z"
git push origin vX.Y.Z
```

**Publish to PyPI before you send a bundle that names the version.** A partner whose
bundle pins a version that is not on PyPI yet gets
`No matching distribution found for appfl-bio-suite==X.Y.Z` on step one, with nothing in
it to act on.

Do not move a tag once a release has been published from it. `git tag -f` gives the
tag and the PyPI wheel different code under one version string — worse than skew,
because nothing detects it. Cut a new version instead.

## Release checklist

- [ ] `pytest` passes, including the slow and legacy-tree marked tests where the old
      source trees are still available
- [ ] `python scripts/check_reusability.py` passes
- [ ] A loopback federation runs for both implemented experiments
- [ ] `pip freeze` in a clean 3.12 environment matches `constraints.txt`
- [ ] Every experiment still has all four documents, non-empty
- [ ] `CITATION.cff` version matches `pyproject.toml` and `__version__`, and every
      documented install command names the new version
      (`tests/test_version_is_consistent.py`)
- [ ] `RELEASE_NOTES.md` describes this release
- [ ] After publishing (below): `pip install "appfl-bio-suite[gwas]==X.Y.Z"` from PyPI
      succeeds in a clean 3.12 environment — this is literally a partner's step one, so
      run it as one
- [ ] No live endpoint UUID, coordinator identity, or partner name appears in the tracked
      tree or in git history

## Publishing on PyPI

Use Python 3.12. The release version must agree in all three files above. Build
and validate both distributions before creating the GitHub release:

```bash
python -m pip install build twine packaging
python -m build
python -m twine check --strict dist/*
python scripts/check_distribution.py
```

`python -m build` builds the wheel from the source archive, checking that the archive
contains the inputs needed by the custom partner-resource packaging step. CI runs
these packaging checks for every pull request.

MultiSuSiE is installed separately with `pip install -r requirements-multisusie.txt`;
PyPI rejects direct Git URL requirements in distribution metadata, even in extras.
The `multisusie` extra installs the coordinator dependencies for that backend.

Publishing uses PyPI Trusted Publishing. The PyPI project `appfl-bio-suite` trusts
workflow `publish.yml` in GitHub repository `mfjoel01/appfl-bio-suite`, environment
`pypi`. No long-lived PyPI token is needed. A fork that publishes under another name
needs its own Trusted Publisher and GitHub Actions environment.

After release validation passes, merge the release branch, push the immutable
`vX.Y.Z` tag, and publish a GitHub release for that tag. `.github/workflows/publish.yml`
checks the tag against the package version, builds and validates the distributions,
and publishes them through PyPI Trusted Publishing. A manual workflow run against
that tag can retry a failure that occurred before upload.

Verify the published version with a clean Python 3.12 install:

```bash
python -m pip install "appfl-bio-suite[gwas]==X.Y.Z" \
    -c https://raw.githubusercontent.com/mfjoel01/appfl-bio-suite/vX.Y.Z/constraints.txt
appfl-bio-suite --version
appfl-bio-suite preflight --check pins   # checks against the wheel's constraints.txt
```

Do not overwrite tags or reuse a published PyPI version. Correct a released defect
with a new version. Keep the Git tag and PyPI package version aligned across sites.
