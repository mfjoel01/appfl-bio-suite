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

**The machine.** A BLAS matrix-vector product adds its terms in whatever order the
library chooses. OpenBLAS chooses its kernel by CPU model, and its split by thread count,
and each choice rounds the low-order bits differently. Nothing raises; the numbers just
follow the hardware.

The GWAS simulation's polygenic score used to be such a product, and the GWAS parity
tests failed intermittently on CI because of it. The first golden checksums were
captured with the Haswell-family kernel of an AMD node. CI passed on hosts that select
that kernel and failed on the rest, so re-running an unchanged commit could turn red
green. Pinning BLAS to one thread removed the thread-count half of the dependence and
left the CPU half in place.

`_compute_pgs` now computes the product itself, in the order the reference BLAS defines,
with no BLAS call at all. The golden checksums are the legacy scripts' own output with
OpenBLAS pinned to the kernel that adds in that same order (`OPENBLAS_CORETYPE=Prescott`,
set by `tests/legacy_harness.py`), so they hold on every machine. Two tests keep it that
way, and unlike the golden comparison on its own, neither depends on which host CI lands
on:

- `test_golden_holds_on_an_older_cpu` re-runs the whole fixture pipeline in a fresh
  interpreter with a different OpenBLAS kernel and numpy's vectorized loops switched off.
  A score routed back through BLAS fails it on any x86-64 machine.
- `test_pgs_is_the_reference_blas_product` holds `_compute_pgs` to the same arithmetic
  written out term by term.

If a parity failure comes and goes between re-runs, the code is following the machine
again. Find the operation that does it. Do not regenerate the golden file, because a
checksum captured on one machine just fails on another. A failure on every machine after
a version bump is the case above: compare the environment against `constraints.txt`.

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
