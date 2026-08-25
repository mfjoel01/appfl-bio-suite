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

## Release checklist

- [ ] `pytest` passes, including the slow and legacy-tree marked tests where the old
      source trees are still available
- [ ] `python scripts/check_reusability.py` passes
- [ ] A loopback federation runs for both implemented experiments
- [ ] `pip freeze` in a clean 3.12 environment matches `constraints.txt`
- [ ] Every experiment still has all four documents, non-empty
- [ ] `CITATION.cff` version matches `pyproject.toml`
- [ ] No live endpoint UUID, coordinator identity, or partner name appears in the tracked
      tree or in git history
