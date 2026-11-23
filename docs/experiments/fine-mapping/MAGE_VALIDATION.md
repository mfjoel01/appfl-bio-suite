# Real-data validation: MAGE v1.0 + 1000 Genomes NYGC 30x

This document records the real-data validation of FedFM ([ABOUT.md](ABOUT.md)) on
**MAGE v1.0** (*Multi-ancestry Analysis of Gene Expression*; Taylor et al., *Nature* 2024),
LCL RNA-seq for 731 1000 Genomes individuals, paired with their **1000 Genomes NYGC 30x
GRCh38** genotypes. The design has three steps:

1. **Replicate the published analysis** with a pooled run of the same model on the same data.
2. **Show each site alone underperforms** that pooled run (same method, same settings).
3. **Show federated equals pooled**, with the covariate adjustment itself federated.

Multi-ancestry methods (MultiSuSiE, SuSiEx) are then reported against the same references.

## Cohort, sample swap, and federation

- **Expression and covariates:** the released inverse-normal TMM expression and the eQTL
  covariate matrix (`sex`, `PC1`–`PC5`, `PEER1`–`PEER60`) from Zenodo record `10535719`.
- **Sample swap:** library `SRR19762653` is labelled `NA11919` (CEU, XY) but is `HG00237`
  (GBR, XX), which also has its own library. The data confirm it: `NA11919`'s expression
  profile correlates best with `HG00237` (r = 0.58; median across samples ≈ 0) and it
  expresses XIST like a female sample. The published analysis used all 731 libraries, so the
  replication runs on `N = 731`; every other comparison uses the clean `N = 730`
  (`196 AFR`, `141 EUR`, `141 EAS`, `139 SAS`, `113 AMR`).
- **Genotypes:** ±250 kb GRCh38 windows sliced from the NYGC 30x phased VCFs over HTTP
  (`vcf_slice.py`), biallelic, no missing calls, MAF ≥ 0.01 in the pooled cohort. All 318
  published credible-set variants in the panel windows are present after extraction.
- **Three sites**, split by 1000 Genomes population: `site_atlantic` (`N = 252`: 85 EUR,
  84 AFR, 56 AMR, 27 SAS), `site_africa_eurasia` (`N = 253`: 112 AFR, 56 EUR, 56 SAS, 29 AMR),
  `site_asiapacific` (`N = 225`: 141 EAS, 56 SAS, 28 AMR).
- **What each site holds:** its genotypes, the **raw** expression for its individuals, and its
  rows of the covariate table (`covariates.tsv`). No phenotype preprocessing happens
  centrally (`build_mage_federation(phenotype="raw")`).

## Matching the published analysis

The published settings below are from the paper's Supplementary Methods (sections 10.1–10.5)
and agree with the authors' pipeline code
([MAGE analysis pipeline](https://github.com/mccoy-lab/MAGE/tree/main/analysis_pipeline/05_eQTL_mapping)).

| Setting | Published MAGE | Here |
| --- | --- | --- |
| Genotypes | NYGC 30x high-coverage calls, MAF > 0.01 in the 731 samples | same; MAF ≥ 0.01 |
| Samples | 731 (includes the swapped library; the paper does not mention the swap) | 731 for the replication, 730 otherwise |
| Expression | TMM-normalized, inverse-normal transformed over all samples | same (released values) |
| Covariates | sex, top 5 genotype PCs, 60 PEER factors, all computed on the 731 samples, removed from **genotypes and expression** | same (FWL) |
| Model | `susieR` 0.12.16 `susie_rss` on FastQTL z-scores, in-sample LD from the adjusted genotypes, one pooled model, L = 10, coverage 0.95, `min_abs_corr` 0.5 | SuSiE on the same sufficient statistics, via MultiSuSiE with one population |
| cis window | ±1 Mb of the TSS | **±250 kb** (all 46 published sets lie inside); ±1 Mb as a sensitivity run |
| Prior variance | `susieR` default (`optim`) | MultiSuSiE `early_EM` |

**Covariate adjustment must include the genotypes.** Regressing covariates out of expression
alone, as the first version of this benchmark did, removes the part of a variant's effect
that lies in the span of the genotype PCs and leaves it in the genotype as noise. At the
high-F<sub>ST</sub> loci this is large: the published lead variant at ERN1 has p = 1.3e-3
with expression-only adjustment and p = 8.7e-10 with both adjusted; NPNT goes from 5.1e-5 to
1.4e-13. The federated pipeline therefore adjusts both (next section).

### Federated covariate adjustment (FWL)

With `covariate_adjustment="fwl"`, each site also sends `Z'X`, `Z'Z`, `1'Z` per (locus,
ancestry) and `Z'y` per phenotype, alongside its usual `X'X`, `1'X`, `X'y`, `1'y`, `y'y`, `n`.
The coordinator sums them, fits `B = (Z1'Z1)⁻¹Z1'X` and `g = (Z1'Z1)⁻¹Z1'y` with
`Z1 = [1, Z]`, and rewrites each block's moments as those of `(X − Z1B, y − Z1g)`
(`covariates.py`). Allele frequencies (MAF filter, SuSiEx `.frq`) stay on the raw dosages.
The covariate moments are off by default, are masked like the others under secure
aggregation, and a one-site run is that site adjusting on its own.

## Stratified 16-locus panel

`select_mage_panel()` chooses 16 non-overlapping autosomal MAGE eGenes outside the extended
MHC whose published credible sets all lie within ±250 kb of the TSS: 46 published sets in all.

1. **PAGE GWAS colocalized:** ILRUN, CETP, NFATC1, TPM4
2. **High cross-ancestry F<sub>ST</sub>:** NPNT (0.78), PCMTD1 (0.72), ARIH1 (0.73), ERN1 (0.76)
3. **Many independent signals (8–9 published sets):** THOC3, DHRS4, ACOT1, RRP7A
4. **One broad published set:** ZNRF2P1, GS1-124K5.4, POT1, NUDT13

## Results

A published credible set is **matched** when a reported set has Jaccard ≥ 0.5 with it and
**identical** at Jaccard 1. "Lead is top" counts loci whose top-PIP variant is the published
lead. PIP r correlates our PIPs with the published PIPs over the published credible-set
variants.

### 1. Replication of the published analysis

| Pooled analysis | N | Sets | Singletons | Matched | Identical | Lead is top | PIP r |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Published MAGE (SuSiE) | 731 | 46 | 25 | — | — | — | — |
| **SuSiE, FWL (published setup)** | 731 | 58 | 28 | **41/46** | 30/46 | 12/16 | **0.91** |
| same, published ±1 Mb window | 731 | 55 | 27 | 43/46 | 31/46 | 12/16 | 0.94 |
| SuSiE, FWL | 730 | 58 | 26 | 37/46 | 26/46 | 12/16 | 0.91 |
| SuSiE, expression-only adjustment | 730 | 48 | 22 | 24/46 | 18/46 | 6/16 | 0.70 |

The pooled SuSiE run reproduces 41 of the 46 published sets (30 exactly) and their PIPs
(r = 0.91); at the published ±1 Mb window, 43 (31 exactly, r = 0.94). Of the four loci whose top variant differs, two have diffuse published leads
(PIP 0.10 at GS1-124K5.4, 0.13 at POT1) and at the other two (DHRS4, NFATC1) we give the
published lead PIP ≥ 0.97 and another variant an equal or higher PIP. It reports 12 more
sets than published at ±250 kb (POT1 +3; THOC3, DHRS4 +2; five loci +1) and 9 more at
±1 Mb; the prior-variance estimator and individual-level versus z-score inputs are the known
remaining differences. Dropping the swapped library (731 → 730) costs four matches.

### 2. Each site alone against the pooled analysis

Pooled SuSiE on all 730 is the reference (58 sets). Each site runs the same model and the
same covariate adjustment on its own data.

| SuSiE (FWL) | N | Sets | Pooled sets recovered | Median set size | Published matched |
| --- | --- | --- | --- | --- | --- |
| site_atlantic alone | 252 | 30 | 15/58 | 4 | 13/46 |
| site_africa_eurasia alone | 253 | 28 | 17/58 | 2 | 18/46 |
| site_asiapacific alone | 225 | 24 | 5/58 | 6.5 | 4/46 |
| **Federated, 3 sites** | 730 | 58 | **58/58** | 3 | 37/46 |

The same comparison for MultiSuSiE (five ancestry columns; pooled reference 57 sets):

| MultiSuSiE (FWL) | N | Sets | Pooled sets recovered | Median set size | Published matched |
| --- | --- | --- | --- | --- | --- |
| site_atlantic alone | 252 | 29 | 12/57 | 7 | 10/46 |
| site_africa_eurasia alone | 253 | 29 | 14/57 | 2 | 14/46 |
| site_asiapacific alone | 225 | 24 | 5/57 | 5.5 | 4/46 |
| **Federated, 3 sites** | 730 | 57 | **57/57** | 3 | 30/46 |

No site recovers more than a third of the pooled result. The Asia-Pacific site, which holds
all 141 EAS samples and little else, recovers the least; two of the three sites also report
wider sets than the pooled analysis (median up to 7 variants against 3).

### 3. Federated equals pooled

| Method (FWL) | Loci | Worst max \|ΔPIP\| | Identical credible sets |
| --- | --- | --- | --- |
| SuSiE, N = 730 | 16 | 1.1e-13 | 16/16 |
| SuSiE, N = 731 | 16 | 1.6e-13 | 16/16 |
| MultiSuSiE, 5 ancestries | 16 | 1.5e-12 | 16/16 |
| SuSiEx, 5 ancestries | 14 | 0 | 14/14 |

The covariate adjustment is computed from the sites' cross-moments, and the result equals
regressing the covariates out of the pooled data directly. SuSiEx does not converge at
ZNRF2P1 and DHRS4 in either arm (1,000 iterations), so its comparison covers 14 loci.

### Covariates a site can compute itself

The published covariates are pooled quantities: expression was rank-normalized across all 731
samples and the 60 PEER factors were estimated on all of them. A federation adjusts for them
exactly (above) but could not compute them. With `covariate_source="site_local"` each site
instead re-applies the inverse normal transform to its own samples and computes 15
expression PCs from all of its genes (MAGE's own choice for n = 142 subsets); sex and the
genotype PCs stay shared, since a site gets those by projecting onto the public 1000 Genomes
reference. Steps 2 and 3 hold under these covariates as well:

| SuSiE, site-derived covariates (N = 730) | Sets | Pooled sets recovered |
| --- | --- | --- |
| site_atlantic alone | 19 | 9/48 |
| site_africa_eurasia alone | 23 | 14/48 |
| site_asiapacific alone | 16 | 5/48 |
| **Federated, 3 sites** | 48 | **48/48** (max \|ΔPIP\| 7e-14) |

Fifteen of sixteen loci have identical sets; at NUDT13 two variants in perfect LD tie at the
95% coverage boundary (component PIP 0.023633 each) and each run includes a different one.

Agreement with the published analysis drops, though, and the drop comes from doing the
preprocessing per site, not from using PCs instead of PEER (pooled SuSiE, N = 731):

| Normalization | Hidden-factor covariates | Published sets matched |
| --- | --- | --- |
| pooled | 60 PEER, pooled (published) | 41/46 |
| pooled | 60 expression PCs, pooled | 42/46 |
| per site | 60 PEER, pooled | 31/46 |
| pooled | 15 expression PCs per site | 27/46 |
| per site | 15 expression PCs per site | 25/46 |
| per site | 30 expression PCs per site | 29/46 |
| per site | 45 expression PCs per site | 26/46 |

Pooled PCs reproduce the published sets as well as PEER does, so a federated PCA (sites sum
their statistics, then project their own samples onto shared loadings) would recover the
pooled-factor result; per-site normalization has a separate cost of about ten sets. The site
covariates also yield fewer credible sets (46–48 against 58), consistent with less of the
technical variation being removed. As elsewhere, matching the published sets measures
agreement with an analysis made on pooled preprocessing, not accuracy.

### Multi-ancestry methods

| Pooled analysis (FWL, N = 730 unless noted) | Sets | Matched | Identical | Lead is top | PIP r |
| --- | --- | --- | --- | --- | --- |
| SuSiE, one population (reference) | 58 | 37/46 | 26/46 | 12/16 | 0.91 |
| MultiSuSiE, 5 ancestries | 57 | 30/46 | 19/46 | 11/16 | 0.68 |
| SuSiEx, 5 ancestries | 22 | 14/46 | 10/46 | 7/16 | 0.47 |

Splitting 730 samples into five ancestries of 113–196 costs MultiSuSiE about a fifth of the
published matches: effects shared across ancestries are estimated from fewer samples each.
Agreement with a SuSiE analysis is not a measure of accuracy, so this does not say the
multi-ancestry sets are wrong; that needs known causal variants (the simulation benchmark).

#### Why SuSiEx falls short, and the fix

With one population SuSiEx should reduce to SuSiE, but as shipped it matches 22/46 published
sets where SuSiE matches 41. Three differences in SuSiEx explain it (pooled, FWL, N = 731):

| Run | Sets | Matched | Identical | PIP r |
| --- | --- | --- | --- | --- |
| SuSiE (reference) | 58 | 41/46 | 30/46 | 0.91 |
| SuSiE with SuSiEx's fixed prior | 59 | 38/46 | 27/46 | 0.91 |
| SuSiEx as shipped | 39 | 22/46 | 15/46 | 0.64 |
| SuSiEx, marginal-p filter off | 63 | 26/46 | 17/46 | 0.70 |
| SuSiEx, finite-sample adjustment | 45 | 32/46 | 24/46 | 0.82 |
| **SuSiEx, adjustment and filter off** | 60 | **39/46** | 29/46 | **0.91** |

1. **Effect standardization (largest).** SuSiEx turns each variant into
   `(BETA/SE)/sqrt(n)` (`data.cpp:311`), which is the marginal correlation only while
   `t² ≪ n`. eQTL leads break that: at t = 25 and n = 731 it reports 0.92 for a correlation
   of 0.68, the statistics no longer fit the LD matrix, and secondary signals are misplaced.
   `susieR::susie_rss` corrects for this when given `n`; the published analysis used it.
   `susiex_finite_sample_adjust` applies the same correction to SuSiEx's inputs
   (`fed_fine_mapping.finite_sample_adjust`); it also removes SuSiEx's two non-converging
   loci.
2. **Marginal-p filter.** `--pval_thresh` (1e-5) drops any credible set with no variant
   below that marginal p, which removes signals that only appear after conditioning.
3. **Fixed prior.** SuSiEx fixes the prior effect variance at the largest squared marginal
   effect (`data.cpp:482`) for every signal; SuSiE estimates one per signal. Forcing SuSiE
   onto SuSiEx's prior costs three sets, the remaining gap.

With five ancestries (N = 730):

| Run | Sets | Matched | Identical | PIP r |
| --- | --- | --- | --- | --- |
| SuSiEx as shipped | 22 | 14/46 | 10/46 | 0.47 |
| SuSiEx, finite-sample adjustment | 24 | 11/46 | 6/46 | 0.47 |
| SuSiEx, adjustment and filter off | 27 | 11/46 | 6/46 | 0.50 |
| MultiSuSiE (for comparison) | 57 | 30/46 | 19/46 | 0.68 |

The three differences account for SuSiEx's single-population gap but not its multi-ancestry
one: with 113–196 samples per ancestry it stays near 11–14/46 whatever these settings, against
30/46 for MultiSuSiE. The likely cause is SuSiEx's per-ancestry model at these sample sizes
(a fixed prior per ancestry and per-ancestry effects estimated separately); that is not
established here, and SuSiEx's multi-ancestry results on MAGE should not be read as a method
comparison.

The adjustment is opt-in (`susiex_finite_sample_adjust: true` in `federation.yaml`), so
existing SuSiEx results are unchanged; it is recommended for molecular QTLs. The federated
path applies the same function as the pooled one, and the offline tests check that federated
and pooled SuSiEx still agree with it on.

### Controls carried over from the first version

The first version of this benchmark (`run_validation_suite`, `local/output/mage-validation`)
used expression-only adjustment; its comparisons still hold as relative controls:

- **Sample swap (731 vs 730), SuSiEx:** worst max |ΔPIP| 0.0135, top variant unchanged at
  14/14 converged loci, identical sets at 11/14.
- **Site-local standardization (`naive_site_standardized`):** each site centres its own
  dosages and the per-site LD matrices and z-scores are combined. Median max |ΔPIP| against
  pooled is 0.004 and it overlaps as many published sets (15/46) as the correct federated
  arm; it diverges strongly at a few loci (worst 0.9985). The reported max |ΔLD| of 0.60–0.81 is
  mostly an artefact of that arm's implementation, which sets a variant's correlation to 0 at
  a site where it does not vary: the values are exact sample fractions (e.g. 112/139 SAS).
  Mean |ΔLD| is 0.014.

## Limitations

- **Pooled preprocessing.** The main analysis uses MAGE's covariates and normalization, which
  were computed on the pooled cohort. The federation adjusts for them exactly but does not
  compute them; computing them per site costs agreement with the published sets (above).
- **Window.** ±250 kb for the main design; the published ±1 Mb window is run as a sensitivity
  analysis for pooled SuSiE only.
- **Sixteen loci, one dataset.** There are no per-locus replicates, so differences between
  arms are descriptive, not tested.
- **The published sets are a reference, not ground truth.** They come from the same samples.

## Reproduction

```bash
# Stage MAGE, slice genotypes, build the federations (raw expression; ±1 Mb; site-derived
# covariates) and run the replication design, then draw its figures
PYTHONNOUSERSITE=1 PYTHONPATH=src python scripts/fine-mapping/run_mage_validation.py \
    --replication --maf-analysis 0.005 --max-workers 8
PYTHONPATH=src python -c "from pathlib import Path; \
from appfl_bio_suite.experiments.fine_mapping.real_data.figures_replication import \
render_replication_figures as r; \
r(Path('local/data/mage/federation_clean_fwl'), Path('local/output/mage-replication'))"

# Diagnostics behind the covariate and SuSiEx tables (see the script's docstring)
PYTHONPATH=src python scripts/fine-mapping/mage_replication_diagnostics.py covariates pooledINT_pooledPC60
PYTHONPATH=src python scripts/fine-mapping/mage_replication_diagnostics.py susiex 1pop_adjusted

# Offline tests (no network): synthetic MAGE fixture, federated FWL exactness, the full design
pytest -v tests/test_fine_mapping_mage.py tests/test_fine_mapping_covariates.py
```

Reports are written to `local/output/mage-replication/reports/` (`replication_vs_mage.tsv`,
`sites_vs_pooled.tsv`, `federated_vs_pooled.tsv`, per-locus `replication_loci_*.tsv`) and
figures to `local/output/mage-replication/figures/`: `figR1_replication` (published PIPs,
covariate adjustment of the leads, which analyses reproduce MAGE), `figR2_sites_and_federation`
(each site alone, federated against pooled) and `figR3_regional` (PIP tracks at ERN1 and NPNT).
The ±1 Mb federations need about 4x the genotype data and LD matrices up to ~12,000 variants
per locus; on a shared node, keep `--max-workers` low. `--all` additionally runs the
first-version suite and its figures (`local/output/mage-validation`), which use
expression-only adjustment.

### Reproduce the original-model and held-out checks

The original-model check uses the released **all-associations** nominal FastQTL results.
Install `pysam` in the benchmark Python environment and provide an R environment with
`susieR` 0.12.16. The following commands use the 731-sample ±1 Mb bundle built above.
The download command extracts the indexed nominal table from Zenodo record 10535719 and
requires about 4.7 GB on disk. Set `RSCRIPT` to the pinned environment's Rscript.

```bash
export PYTHONPATH=src
nominal=local/data/mage/nominal/eQTL_FastQTL_results.nominal_pass.allAssociations.MAGE.v1.0.txt.gz
python scripts/fine-mapping/mage_original_susie.py --nominal "$nominal" --download
python scripts/fine-mapping/mage_original_susie.py --nominal "$nominal" --rscript "${RSCRIPT:-Rscript}"
python scripts/fine-mapping/mage_original_susie.py --nominal "$nominal" --report

python scripts/fine-mapping/mage_heldout.py --stage build
python scripts/fine-mapping/mage_heldout.py --stage fit --workers 6
```

The R check uses the authors' slope/SE, in-sample covariate-adjusted LD, seed, model options,
and per-variant credible-set export convention. File loading uses the staged PLINK dosages
aligned to ALT and the released nominal variant IDs. The runner checks the tested variant
count for every gene against the release. The report compares set overlap in both
directions for the R fit and the earlier ±1 Mb pooled fit. The held-out comparison uses
the clean 730-person cohort and a [prespecified 22-locus manifest](mage/heldout_design.md).
The authors' [model invocation and export code](https://github.com/mccoy-lab/MAGE/blob/f929fb5106055a124c0d18e6431829d428717452/analysis_pipeline/05_eQTL_mapping/run_susie/scripts/run_susie.R)
record the upstream implementation.
