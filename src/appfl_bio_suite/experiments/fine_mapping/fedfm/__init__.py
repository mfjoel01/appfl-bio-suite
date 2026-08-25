"""The FedFM science, vendored verbatim from the standalone simulation repository.

WHAT THIS PACKAGE IS
--------------------
Eight modules copied byte-for-byte from ``fedfm-simulation/src/``. They are the
statistical core of the experiment: the cohort simulator, the centralized SuSiEx
baseline, and the federated site->coordinator protocol whose derivation and proofs live
in ``docs/experiments/fine-mapping/reference/algo.md``.

They arrived already correct and already tested. The federated path is not an
approximation of the centralized one -- by Corollary 1 of algo.md it *equals* it, and
``tests/test_fine_mapping_federated.py`` holds it to that on a real fit. Code carrying a
proof and a test that pins its numbers is code to wrap, not to rewrite.

WHY VERBATIM, AND WHAT THAT COSTS
---------------------------------
Nothing in here has been reformatted, renamed, or restructured to match the suite's
conventions. Two consequences worth knowing before editing:

* Their docstrings still say ``src.fine_mapping`` where the module is now
  ``appfl_bio_suite.experiments.fine_mapping.fedfm.fine_mapping``. That is stale prose,
  not a stale import -- there are no cross-package imports to fix, because these modules
  only ever imported each other, and they still only import each other.
* They use their own config object (:class:`utils.SimulationConfig`, a pydantic model
  loaded from a YAML that the standalone repository owns), not the suite's
  ``federation.yaml``. The adapter layer above translates; see ``simulation/schema.py``.

Keeping the diff against upstream empty is what makes it possible to answer "is this the
code the paper's numbers came from?" with ``cmp`` rather than with a reading. If a fix is
needed here, make it upstream and re-copy, or the answer stops being checkable.

WHERE THE SUITE PLUGS IN
------------------------
    fine_mapping/trainer.py     re-implements the SITE stage, self-contained, because it
                                is shipped to partner workers and may not import this
                                package. The duplication is deliberate and is tested for
                                agreement against ``fed_fine_mapping.site_aggregates``.
    fine_mapping/aggregator.py  calls THIS code for the coordinator stage -- pool_geno,
                                build_columns, pool_pheno, fed_finemap_instance -- so the
                                APPFL path and the standalone path run the same
                                statistics rather than two implementations of them.
    fine_mapping/simulation/    orchestrates sampling -> locus_selection -> phenotype_sim
                                under the suite's manifest and provenance contract.

MODULE ROLES
------------
    utils.py              config model, deterministic seeding, PLINK wrapper, .bed reader
    sampling.py           stage 1: disjoint per-site cohorts out of HAPNEST
    locus_selection.py    stage 2A: LD-divergence scoring, stratified locus choice
    phenotype_sim.py      stage 2B: causal architectures, effect sizes, phenotypes
    qc.py                 stage 3: relatedness, MAF, LD decay; one HTML report per site
    validation.py         stage 4: the five-check self-validation harness
    fine_mapping.py       stage 5: centralized SuSiEx baseline (the reference fit)
    fed_fine_mapping.py   stage 6: the federated path -- site aggregates, then pooling

None of it ever runs on a partner cluster. It is coordinator-side in full, which is why
it may depend on joblib, matplotlib, jinja2 and a PLINK/SuSiEx toolchain that no partner
is asked to install.
"""

__all__ = [
    "utils",
    "sampling",
    "locus_selection",
    "phenotype_sim",
    "qc",
    "validation",
    "fine_mapping",
    "fed_fine_mapping",
]
