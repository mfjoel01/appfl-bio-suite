"""Federated GWAS with polygenic-score evaluation.

Single-round summary-statistic federated learning. Each site runs a complete local GWAS
and returns per-variant effect sizes, standard errors and allele frequencies; the server
combines them by inverse-variance fixed-effect meta-analysis. Genotypes never leave the
site they originate from.

Module roles, which differ in a way that matters:

    dataset.py     SHIPPED to workers -- self-contained, no suite imports
    trainer.py     SHIPPED to workers -- self-contained, no suite imports
    aggregator.py  coordinator-side   -- may import freely
    plotting.py    coordinator-side   -- may import freely
    simulation/    coordinator-side   -- never reaches a partner at all

See dataset.py's module docstring for what "shipped" means and why it constrains imports.
"""
