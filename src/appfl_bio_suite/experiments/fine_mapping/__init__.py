"""Federated cross-ancestry statistical fine-mapping.

Single-round aggregate federated learning. Each site computes raw second-moment
aggregates over its own individuals at a locus and returns only those; the coordinator
sums them, standardizes once against the pooled moments, and runs SuSiEx. Genotypes and
phenotypes never leave the site they originate from.

The claim this experiment makes is stronger than the GWAS experiment's. Inverse-variance
meta-analysis *approximates* a pooled analysis under an assumption about the sites. This
does not approximate anything: by Corollary 1 of the derivation in
``docs/experiments/fine-mapping/reference/algo.md``, the federated fit **equals** the
centralized fit -- same credible sets, same posterior inclusion probabilities, up to
floating-point summation order. A discrepancy is a bug, never a cost of federating, and
the test suite is written to that standard.

WHAT LEAVES A SITE
------------------
Per (site, ancestry, locus), on the raw 0/1/2 dosage scale, over the locus window only::

    G = X'X   (M x M)      u = 1'X   (M)      n

and per instance (locus x architecture x replicate) additionally::

    c = X'y   (M)          q = 1'y            w = y'y

Nothing else. No genotype row, no phenotype value, and -- importantly -- no locally
standardized LD matrix. Standardizing at the sites and combining the results is biased,
because each site would centre against its own column means and delete the between-site
frequency variation before the coordinator can see it. That is the one real trap of the
design; algo.md A.10 gives the worked example where it produces 0.866 against a truth of
0.816, and ``tests/test_fine_mapping_federated.py`` pins that number.

Privacy is out of scope: the genotypes are synthetic, so aggregates ship in the clear.
No secure aggregation, no differential privacy. See DATA.md.

MODULE ROLES, which differ in a way that matters
------------------------------------------------
    dataset.py     SHIPPED to workers -- self-contained, no suite imports
    trainer.py     SHIPPED to workers -- self-contained, no suite imports
    aggregator.py  coordinator-side   -- may import freely
    plotting.py    coordinator-side   -- may import freely
    fedfm/         coordinator-side   -- the vendored science, verbatim from upstream
    simulation/    coordinator-side   -- never reaches a partner at all

See ``dataset.py``'s module docstring for what "shipped" means and why it constrains
imports, and ``fedfm/__init__.py`` for why the site stage exists twice.

THE COST OF FEDERATING THIS, IN BYTES
-------------------------------------
``G`` is M x M. At a realistic locus (M ~ 2,000 variants after the MAF filter) one
ancestry block is ~32 MB, and a site holding five ancestries uplinks ~160 MB per locus.
That is the honest price of shipping LD structure rather than marginal statistics, and it
is why a run is sharded by locus rather than sweeping all 100 in one exchange. The
per-locus payload is reported by ``preflight`` and configured by ``locus_n_shards`` /
``locus_shard_index``; see RUNBOOK.md.
"""

__all__ = ["dataset", "trainer", "aggregator", "plotting", "fedfm", "simulation"]
