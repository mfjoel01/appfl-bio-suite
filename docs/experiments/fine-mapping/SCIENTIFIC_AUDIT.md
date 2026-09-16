**Scientific audit — 14 September 2026**

Implementation follow-up: [corrected protocol and rerun](SCIENTIFIC_RERUN.md). The findings below describe the audited pre-fix revision and remain as historical evidence.

**Verdict: the moment-pooling mathematics is sound under its stated conditions, but the current experiment does not support an unqualified claim of exact centralized/federated equivalence or calibrated fine-mapping. Scientific sign-off is withheld.** There are reproducible differences in the clean production results, stale deployment bundles, and errors in outcome definitions. This audit changes no production data, estimator, or existing results.

Audited repository: `appfl-bio-suite`, commit `765d11ecd8c5819c9272d2d5af0f0d956dbb92c4`. Scope includes the simulator, centralized and federated estimators, shipped trainer/aggregator, participation arms, validation, reporting, current local configuration, and available production artifacts. The separate `hive-watch-ui` checkout and remote partner installations were not independently audited.

**Evidence and verification**

- All **130 tests** in the nine `tests/test_fine_mapping_*.py` modules passed, including the real-binary centralized/federated acceptance tests. Log: [pytest.log](../../../local/output/fine-mapping-audit/pytest.log). Interpreter: `/soft/applications/miniconda3/3.12/bin/python`; Python 3.12.2, NumPy 1.26.4, pandas 3.0.2, SciPy 1.17.1. pandas/SciPy differ from the project pins; this is not a pinned-environment certification.
- Inspected **11,850 clean centralized rows**, **3,555 clean federated rows**, and five additional arms with **3,555 rows each**. All have empty execution-error fields. This does not establish statistical convergence.
- Compared all 533,532 BIM variant orientations at every site against the HAPNEST reference. The clean cohort has zero discrepancies; the legacy cohort does not.
- Re-ran a production failure, with all 150,000 enrolled individuals, using the checked-out implementation and installed PLINK/SuSiEx binaries. Then exchanged the two paths' summary-statistic and LD inputs to isolate the discrepancy.
- Reproduction scripts, paired differences, captured SuSiEx output, and machine-readable counts are retained in [the audit directory](../../../local/output/fine-mapping-audit/). Start with [inspection.json](../../../local/output/fine-mapping-audit/inspection.json), [reproduce_parity.py](../../../local/output/fine-mapping-audit/reproduce_parity.py), and [reproduction/results.json](../../../local/output/fine-mapping-audit/reproduction/results.json).

**1. High — the actual inputs violate the equality conditions, with material downstream differences.**

Across the 3,555 shared instance keys:

| Outcome | Disagreeing instances |
| --- | ---: |
| Number of credible sets | 138 |
| Sorted credible-set sizes | 158 |
| Number of causal variants captured | 139 |
| Any causal variant captured | 137 |
| Reported convergence flag | 137 |

Among 3,279 pairs with nonmissing `causal_pip_max`, 78 differ by more than 0.0001; the maximum absolute difference is 0.999766. These are not solely last-digit changes.

At `L0002_ncsl1_h2-0.001_rg1_rep2`, a fresh centralized run returns a two-SNP credible set containing the causal variant `chr1:62133393:C:T`, with causal score 0.272666. The federated run returns no set. Federated summary statistics with centralized LD still return no set; centralized summary statistics restricted to the federated SNP lists also return no set. Conversely, centralized summary statistics with federated LD retain the two-SNP set, with causal score 0.272663. Thus the summary-statistic input selection is sufficient to reproduce this failure; a general diagnosis of all 138 count discrepancies would require additional instances.

A further controlled split isolates **MAF-based summary-statistic restriction** in this example: restricting centralized summary statistics to each centralized MAF-filtered LD panel, while retaining its incomplete variants, loses the set. Removing only the three additional incomplete variants retains it (causal score 0.272667). This is not explained solely by missing genotypes or floating-point summation. See [isolated_filters.json](../../../local/output/fine-mapping-audit/reproduction/isolated_filters.json). On SNPs present in both input tables, absolute beta and SE agree exactly at the emitted precision in all six ancestries.

The coordinator writes statistics only for `FedColumn.variants`, after completeness and MAF filtering; centralized GWAS writes all estimable window variants. Additionally, the federated site stage discards any variant with a missing genotype, whereas PLINK uses variant-specific complete observations for GWAS. At this locus the centralized LD lists include one additional EUR and two additional CSA variants with incomplete genotypes. The production logs explicitly warn that equality does not hold. See [site filtering and column construction](../../../src/appfl_bio_suite/experiments/fine_mapping/fedfm/fed_fine_mapping.py), especially `site_aggregates`, `build_columns`, and `fed_finemap_instance`, and [centralized GWAS](../../../src/appfl_bio_suite/experiments/fine_mapping/fedfm/fine_mapping.py), `run_gwas`.

**Required correction:** define identical ancestry-specific sample, variant, missingness, MAF, allele, and input-serialization policies for both arms. For an exactness experiment, run the centralized comparator under the federated completeness rule, or implement matching sufficient statistics for the intended missingness treatment. Treat a different PLINK missingness policy as a separate sensitivity comparator. Make violated exactness conditions a failed acceptance gate, not only a log warning. Re-run paired production comparisons after correction, including full CS membership and variant-level probabilities.

**2. High — the configured federation still serves the legacy phenotype-generation defect.**

`local/federation.yaml` points at `local/data/fine-mapping-bundles`, including its causal manifest. All three bundles' BED/BIM/FAM files are the **same inodes** as the legacy `local/data/fine-mapping/processed` files; sampled phenotype files are also hard-linked to the legacy phenotype tree. They are not the clean artifacts.

| Site | Legacy variants reversed against reference | Legacy instances containing a reversed causal variant |
| --- | ---: | ---: |
| ANL | 9,398 | 421 |
| Covenant | 31,031 | 1,293 |
| MBZUAI | 9,777 | 459 |

The simulator computes `g = X beta` from each site's allele coding. Reversing a causal dosage changes the realized effect to `-beta`, plus an intercept shift. Harmonizing dosages later in the trainer cannot repair an already-generated phenotype. The corrected sampler now preserves allele order, and the clean cohort's BIM files are aligned, but that does not update existing bundles. See [phenotype generation](../../../src/appfl_bio_suite/experiments/fine_mapping/fedfm/phenotype_sim.py), `simulate_phenotype_for_site`, and [sampling](../../../src/appfl_bio_suite/experiments/fine_mapping/fedfm/sampling.py), `_extract_plink_for_site`.

**Required correction:** generate new bundles and provenance from the clean data, update the local deployment configuration, and validate their phenotype/variant identities before interpreting another federation run. Remote bundles require their own verification. Preserve the legacy artifacts as explicitly labeled legacy results.

**3. High — reported PIPs are component maxima, not overall inclusion probabilities.**

`_pip_by_snp` takes the maximum across `PIP(CSk)` columns. These are component probabilities. The installed SuSiEx source separately computes overall PIP as `1 - product(1 - alpha_k)` over retained components and writes `OVRL_PIP` for CS members. Two component probabilities of 0.4 yield overall inclusion 0.64, whereas this parser reports 0.4. This affects causal-PIP summaries and the downstream 0.5/0.95 threshold flags in [paper_plots.py](../../../src/appfl_bio_suite/experiments/fine_mapping/figures/paper_plots.py).

The parser also averages only causal variants present in the PIP map, silently changing the denominator when causal variants were filtered. `_write_rollup` labels the mean of `causal_pip_max` as `mean_causal_pip`.

**Required correction:** distinguish component probabilities, retained-component overall PIP, and any full-model PIP. Match SuSiEx's output convention explicitly, track causal-variant eligibility, and report unconditional causal recovery separately from conditional performance among eligible variants. Recompute affected summaries and figures. Evidence: [parser](../../../src/appfl_bio_suite/experiments/fine_mapping/fedfm/fine_mapping.py), `_pip_by_snp`/`parse_susiex`; installed `local/vendor/fine-mapping/SuSiEx/src/model.cpp`, `pip_overall`/`write_cs`.

**4. High — neither convergence nor 95% calibration is measured reliably by the summary table.**

`parse_susiex` sets `converged` from the presence of a CS or PIP column. SuSiEx can converge and then filter out all CSs; its writer emits `NULL` in that case. The summary flag consequently conflates a valid no-discovery result with optimizer failure. `run_susiex` discards the captured stdout/stderr needed for diagnosis, and working files are deleted by default.

`pap4_power_coverage` substitutes the instance-level probability of capturing any causal variant when per-CS details are unavailable, but still draws a nominal 95% coverage line and titles the plot as calibration. Power and credible-set coverage have different denominators. The footnote does not make their comparison valid. The optional detail harvester contains the correct per-CS membership information, but no such retained detail table was found under the current clean results tree. The SuSiEx study evaluates power, calibration, and resolution separately. [Primary study](https://www.nature.com/articles/s41588-024-01870-z).

**Required correction:** retain diagnostic status and per-CS/per-variant records during every run. Separate execution failure, nonconvergence, converged/no-CS, and converged/CS. Disable calibration plots without actual CS membership and truth. Report empirical coverage, false discoveries, power, and resolution with explicit denominators and uncertainty; do not infer 95% coverage from the requested CS level.

**5. High — phenotype “re-derivation” can pass when phenotypes are assigned to the wrong people.**

`validation._rederive_row` recomputes `Var(g) / Var(y)` and compares it with manifest heritability. It does not reconstruct the phenotype using the recorded noise seed or check the genotype–phenotype association. Any permutation of `y` preserves its variance, so this check can pass after destroying the biological relationship being tested. This contradicts the stronger re-derivation assurance in [DATA.md](DATA.md).

**Required correction:** independently reconstruct the expected phenotype vector, including the deterministic noise stream, and compare by FID/IID. Add a negative control that permutes phenotype assignment while preserving values and variance. Validate residual behavior as a separate model check. Existing manifest heritabilities in the clean cohort are within 2.86% of targets, which is reassuring about variance calibration but does not resolve this gap.

**6. Medium — the “50k, diversity only” control is not actually matched on analyzed sample size or information.**

The downsampled enrollment file contains 49,999 individuals. `min_gwas_n=1000` excludes the 833 EAS individuals; the five emitted ancestry columns total **49,166**, versus 50,000 in each solo arm. The one-person rounding difference is negligible; the excluded ancestry is a substantive mismatch in the declared control. See [arms.py](../../../src/appfl_bio_suite/experiments/fine_mapping/arms.py), `build_downsampled_site_dir`/`build_configs`, and `plan_columns` in the centralized implementation.

Even after fixing N, ancestry composition also changes allele frequencies, effect sizes at `rg<1`, and the distribution of site-specific phenotype noise. Therefore a remaining arm difference cannot be attributed uniquely to LD diversity. The borrowed-LD arm uses another site's **same-ancestry** reference panels, drawn from the same synthetic population pool; it tests external-panel sampling differences, not an EUR panel substituted for AFR.

**Required correction:** match N after ancestry filtering, report actual per-column N and genetic signal/noise, use several downsampling seeds, and describe the estimand as a composition/participation effect. To isolate LD diversity, add controlled effect/noise and ancestry-composition comparisons. Report borrowed-LD results without presupposing that a large mismatch must occur.

**7. Medium — causal eligibility and model capacity are not aligned with all simulated architectures.**

The wrapper leaves SuSiEx's defaults for ambiguous-SNP handling and number of signals. Its documented defaults remove A/T and C/G variants and fit at most five signals. [Official interface](https://github.com/getian107/SuSiEx). In the clean manifest, **3,333 of 11,850 instances** include at least one such ambiguous causal variant. All **30 divergent instances** contain **seven union causal variants**, despite their architecture IDs specifying two per ancestry. Thus some truths are excluded by default and the divergent stress case exceeds the default signal capacity. These may be intentional stressors, but are not currently separated adequately in outcome summaries.

**Required correction:** explicitly configure and record these options. In this controlled reference-aligned simulation, either retain resolvable ambiguous variants or report their exclusion as an eligibility stratum. Set and assess signal capacity against the union causal set; report shared/divergent architectures and ancestry-specific truth separately. Thirty divergent instances concentrated in a few locus/architecture assignments are insufficient for broad claims about ancestry-private causal architectures.

**8. Medium — the simulation parameters support narrower scientific claims than the prose implies.**

`draw_effect_sizes` sets a correlation among raw per-allele effect draws. With one to three causal variants, ancestry-specific frequencies/LD, and structural zeros in divergent cases, this is not a guarantee of the same realized genetic correlation for every phenotype. `run_effect_size_correlation` pools fully nonzero beta rows and therefore validates the shared-effect drawing mechanism, not all realized cross-ancestry genetic signal.

Noise is calibrated from the variance of raw genetic values across each ancestry-mixed site. That variance includes between-ancestry genetic-mean differences. Pooled ancestry regressions consequently need not have the nominal site-level heritability or homogeneous residual variance. The claim in `fine_mapping.py` that heteroscedasticity costs nothing in correctness is stronger than the validation supports: equality of two implementations does not establish calibrated standard errors or posterior uncertainty.

**Required correction:** label `rg` as the effect-draw correlation parameter and distinguish marginal site heritability from within-ancestry signal. Measure per-site/per-ancestry residual variances and association calibration; include homogeneous-noise and heteroscedastic sensitivity scenarios. Bound conclusions to the synthetic HAPNEST chr1 design rather than treating arbitrary site compositions as representative institutional cohorts. HAPNEST is a synthetic-data resource whose realism is evaluated through selected genetic properties. [HAPNEST study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10493177/).

**9. Medium — uncertainty and locus-selection provenance need stronger treatment.**

Most rate intervals use an ordinary Wilson interval over rows. Replicates share loci, LD, cohorts, and often architecture parameters; multiple credible sets share an instance. For claims generalized across loci, rows are not independent locus draws. Use a paired, locus-clustered bootstrap, with a clearly stated target population and sensitivity to the 79-locus sample. Do not present 11,850 simulation rows as 11,850 independent genomic regions.

The clean configuration reuses the legacy selected loci. Although a within-site r-squared matrix is invariant to an allele flip, `_pick_tag_snps` first concatenates unharmonized site genotypes and prunes that pooled matrix. That operation is not invariant to site-specific flips. Thus the archived assertion that locus selection was entirely unaffected is too broad. Reuse can define a fixed benchmark panel, but it does not prove that a fully corrected pipeline would select the same panel. Cached window scores are also keyed by file existence rather than a validated input/config fingerprint.

**Required correction:** either document these 79 loci as a frozen legacy-selected panel, with harmonized re-scoring, or regenerate selection under the corrected coding. Separate LD magnitude, frequency differences, and between-site divergence in interpretation; the score alone does not identify the effect of diversity. Record input/config/tool hashes and invalidate stale caches.

**What is supported now**

The additive raw-moment identities and subsequent pooled centering are correct for aligned, identically handled data. Tests exercise these identities and shipped/vendored parity. The clean BIM files fix the observed allele-coding defect, architecture/result counts agree with 79 loci and the specified replicate subsets, and the participation merger pairs shared instance keys. These are useful foundations, but none overrides the production discrepancies above.

**Acceptance before scientific sign-off**

1. Replace stale deployment bundles and verify clean phenotype identities independently.
2. Align analysis inputs and demonstrate production-scale paired equality under the declared exactness conditions, with no silent completeness mismatch.
3. Correct probability/status/coverage definitions and retain enough raw output to recompute them.
4. Explicitly handle excluded causal variants, divergent union capacity, analyzed sample-size matching, and noise/composition controls.
5. Re-run affected analyses and publish paired, locus-aware uncertainty with bounded claims.

This audit did not regenerate the entire cohort, rerun every production instance, perform a new full-cohort kinship scan, or inspect remote partner files. Installed binary hashes are recorded below to identify the fresh reproduction:

```
SuSiEx  50d26474ea71dfc9ad75e46355a3d4a188364ed0301ea7831ccefb159ba78acd
plink   013114acc9db095b78588ee25bae523419db4965b04c111c1c0e2fc2b43409d5
plink2  28e45a65689bacf2b6675898abc6fc2eb2a5de19ae21a60ccec0dd37e6b21ce2
```
