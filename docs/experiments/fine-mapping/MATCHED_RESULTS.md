# Completed matched-filter MultiSuSiE–SuSiEx comparison

All 931 cohort/locus tasks completed: 139,650 paired instances and 558,600 method/policy rows. The final acceptance gate passed on October 30, 2026. This replaces the earlier pilot-only evidence; the full matched-filter comparison is now available.

## Takeaway

Under the common purity ≥0.5 policy in the full 150,000-person held-out cohort, MultiSuSiE has slightly higher causal-variant recall, while SuSiEx has coverage closer to the nominal 95%, lower high-confidence FDR, and smaller credible sets. Any-causal capture and high-confidence causal yield show no clear difference. These results support a recall–calibration tradeoff, not universal estimator superiority. For calibrated, compact credible sets in this configuration, the evidence favors SuSiEx.

## Full-cohort held-out results

54 loci × 150 instances = 8,100 paired instances. Differences are MultiSuSiE minus SuSiEx; intervals are 95% paired locus-bootstrap intervals (2,000 replicates). Percentage-metric differences are percentage points. Mean set size averages all retained sets, unlike the earlier poster statistic for the smallest causal-containing set.

### Common purity ≥0.5

| Metric | MultiSuSiE | SuSiEx | Difference [95% CI] |
|---|---:|---:|---:|
| Any-causal capture (%) | 98.27 | 98.23 | +0.04 [-0.16, +0.25] |
| Causal-variant recall (%) | 69.85 | 68.58 | +1.27 [+0.96, +1.57] |
| Credible-set coverage (%) | 92.88 | 95.09 | -2.22 [-2.68, -1.77] |
| Mean credible-set size (variants) | 7.84 | 7.38 | +0.46 [+0.31, +0.61] |
| FDR at PIP >0.95 (%) | 4.23 | 3.19 | +1.05 [+0.56, +1.56] |
| Any true causal with PIP >0.95 (%) | 48.65 | 48.54 | +0.11 [-0.21, +0.44] |

### Unfiltered native sets

| Metric | MultiSuSiE | SuSiEx | Difference [95% CI] |
|---|---:|---:|---:|
| Any-causal capture (%) | 98.84 | 99.28 | -0.44 [-0.62, -0.28] |
| Causal-variant recall (%) | 78.75 | 91.65 | -12.90 [-13.48, -12.32] |
| Credible-set coverage (%) | 86.66 | 67.65 | +19.01 [+18.45, +19.58] |
| Mean credible-set size (variants) | 141.03 | 337.93 | -196.90 [-211.56, -182.79] |
| FDR at PIP >0.95 (%) | 4.23 | 3.19 | +1.05 [+0.56, +1.56] |
| Any true causal with PIP >0.95 (%) | 48.65 | 48.54 | +0.11 [-0.21, +0.44] |

Without the external purity filter, SuSiEx captures more causal variants but produces much larger sets and substantially lower set coverage. MultiSuSiE also falls below nominal coverage. Thus the output policy changes the recall comparison; unfiltered recall alone does not establish a better estimator. PIP metrics are identical across policies by design because post-hoc set filtering does not alter native PIPs.

## Other cohorts and benchmark

The full-cohort benchmark repeats the common-filter tradeoff: MultiSuSiE recall is 1.52 percentage points higher, coverage 1.77 points lower, FDR 0.59 points higher, and sets 0.46 variants larger. These four paired intervals exclude zero.

The selected 50,000-person design does not reproduce the full-cohort recall advantage: SuSiEx recall is slightly higher in each of the three fixed sampling draws (by 0.23–0.41 percentage points). Coverage and FDR differences in those draws have intervals spanning zero. Cohort-specific results should not be collapsed into a universal winner. All seven benchmark and seven evaluation cohort configurations, plus architecture/causal-mode strata, are retained in the TSV files.

## Interpretation limits and provenance

The people, phenotypes, eligible variant union and external set filters are matched. Native priors, active-component/PIP definitions and precision remain method-specific (MultiSuSiE rho=0.75; SuSiEx float32 LD). The development-selected cohort was fixed before this comparison and selected using MultiSuSiE. These are results for the frozen synthetic design and configured methods, not every possible tuning or population. Intervals are marginal, without multiplicity adjustment; overlapping cohort draws are not independent replications. Archived timings do not support speed claims.

`design.json` and `snapshot.json` record the frozen design. `accepted.json` hashes the complete numerical outputs. `comparison.tsv` contains the full paired grid; `summary.tsv` and `paired_differences.tsv` contain overall estimates and paired intervals; the corresponding `_stratified.tsv` files contain architecture/causal-mode results. The original MultiSuSiE poster remains the original pipeline report; this report supplies the completed matched-estimator comparison.

Artifact paths, poster scope and relocation verification are documented in
[the archive guide](MULTISUSIE_ARTIFACTS.md).

The [comparison figure guide](COMPARISON_FIGURES.md) lists the dedicated offline
report, ten comparison figures, complete tables and reproduction command in the
primary checkout.
