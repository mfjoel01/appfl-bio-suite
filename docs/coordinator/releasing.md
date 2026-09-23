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

**The machine's core count.** The polygenic-score accumulation in `phenotypes.py` is a
BLAS matrix-vector product, and given eight or more threads OpenBLAS splits its reduction
and sums the parts in a different order. The scores move in their low-order bits, every
phenotype and summary file inherits the difference, and nothing raises. It reproduces on
a laptop and fails on a 64-core node with the source unchanged, which is exactly what
happened once: a CI runner grew and the committed golden file stopped matching.

`_compute_pgs` therefore pins that accumulation to a single BLAS thread, via
`threadpoolctl`, which is why that package is a declared dependency of the `gwas` extra
rather than an incidental one from scikit-learn.
`test_pgs_accumulation_does_not_follow_the_machine` fails if the pin is removed. Do not
"fix" a thread-count parity failure by regenerating the golden file: the golden is the
single-threaded result, and it is the published one.

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

Partners install from a git tag, not from PyPI and not from `main`
(`docs/partner/endpoint-setup.md`, Step 1). `main` moves; a tag does not. Two sites that
installed from `main` a week apart are running different code with no way to say which,
which is the same failure mode `constraints.txt` exists to prevent — just one level up.

Three things carry the version and must move together:

| | |
| --- | --- |
| `pyproject.toml` | `version` |
| `src/appfl_bio_suite/__init__.py` | `__version__` — `INSTALL_TAG` is derived from it |
| `CITATION.cff` | `version` |

Then tag and push:

```bash
git tag -a v0.1.0 -m "appfl-bio-suite 0.1.0"
git push origin v0.1.0
```

**Push the tag before you send a bundle that references it.** A partner whose bundle names
`v0.1.0` while the tag exists only on your machine gets a bare
`fatal: Remote branch v0.1.0 not found`, on step one, with nothing in it to act on.

Do not move a tag once a partner has installed from it. `git tag -f` gives two sites the
same version string over different code — worse than skew, because nothing detects it.
Cut a new version instead.

## Release checklist

- [ ] `pytest` passes, including the slow and legacy-tree marked tests where the old
      source trees are still available
- [ ] `python scripts/check_reusability.py` passes
- [ ] A loopback federation runs for both implemented experiments
- [ ] `pip freeze` in a clean 3.12 environment matches `constraints.txt`
- [ ] Every experiment still has all four documents, non-empty
- [ ] `CITATION.cff` version matches `pyproject.toml` and `__version__`
- [ ] The `vX.Y.Z` tag is pushed to `origin`, and
      `pip install "appfl-bio-suite[gwas] @ git+<repo>.git@vX.Y.Z"` succeeds in a clean
      3.12 environment — this is literally a partner's step one, so run it as one
- [ ] No live endpoint UUID, coordinator identity, or partner name appears in the tracked
      tree or in git history
