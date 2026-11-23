# MAGE held-out federation comparison

Select one fine-mapped MAGE eGene per autosome from the v1.0 catalog. Exclude the extended MHC and genes within 2 Mb of an original panel TSS on the same chromosome. Within each chromosome take the smallest hex SHA-256 digest of the UTF-8 string formed by appending the versioned Ensembl gene ID to `MAGE-FedFM-holdout-v1/`. Eligibility and ranking are independent of credible-set location, size, PIP, and federation results. The committed manifest fixes all 22 loci before fitting.

Use the clean 730-person cohort, the same three population-defined sites, released expression and covariates, ±250 kb windows, and the original FWL SuSiE and five-ancestry MultiSuSiE settings. Compare each federated fit to its pooled fit. Report variant coverage, maximum absolute PIP difference, and exact credible-set membership across all 22 selected loci for each method. Retain failed or discordant loci in the denominator.
