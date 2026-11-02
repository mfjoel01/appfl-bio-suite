# Corrected fine-mapping rerun — 14 September 2026

The [scientific audit](SCIENTIFIC_AUDIT.md) led to protocol
`complete-variants-pooled-maf-v3`. Previous result tables and legacy bundles remain
historical artifacts. Scientific sign-off is conditional on full phenotype validation,
complete-grid parity, and the observed calibration results.

Both paths use the same complete-genotype, pooled-ancestry MAF-eligible variants for
association statistics and LD. Centralized GWAS uses PLINK2 independently; its effect
alleles are explicitly aligned to the genotype reference, reversing beta and the test
statistic on an allele swap while retaining SE and p. Centralized in-sample LD is computed
directly from pooled genotypes in float64 and serialized to float32. PLINK's own LD
arithmetic had produced larger rounding differences at production scale. Per-ancestry
sample counts are checked, and shared cohort keep files are replaced atomically so a
parallel shard cannot read a partially written list.

Before inference, a hard gate requires each ancestry's GWAS and LD lists to agree in
SNP identifier, position, allele coding and order. Missing or incompatible entries stop
the analysis batch after diagnostics are saved. The final paired gate also compares the
ancestry-specific input-list fingerprints, sample sizes and options **between** centralized
and federated paths, then checks every fit status, credible-set membership and retained
PIP. Matching pooled SNP unions alone is insufficient. Missing instances or any mismatch
fail acceptance; there is no "99.x% is close enough" path. PIP comparisons retain the
pre-existing output-precision tolerance (`rtol=2e-5`, `atol=1e-12`). Converged and jointly
nonconverged pairs are counted separately: matching failures do not establish posterior
accuracy. This establishes observed equivalence on the 79-locus benchmark if it passes,
not a universal numerical guarantee or a posterior-calibration theorem.

The borrowed external-LD sensitivity arm deliberately uses a separate reference panel.
Its input-list differences are recorded explicitly; it is excluded from the centralized
versus federated exactness claim. This exception does not relax either paired path.

SuSiEx settings are explicit: ten signals, 1,000 iterations, tolerance 1e-6, 95% sets,
p-value filter 1e-5, and retained ambiguous SNPs. Ambiguous alleles have a verified common
reference coding. Ten signals exceed the largest simulated causal union (seven).
Probability summaries use `1 - product(1 - alpha)` over **retained components**, matching
the emitted SuSiEx output; they do not reconstruct components the binary omits. Successful
fits without retained sets (`NULL`) and nonconvergence (`FAIL`) are separate outcomes.
Excluded truths count as zero unconditional recovery; eligible-truth summaries are
separate. Unknown eligibility after a failed fit is missing. Commands, sumstats,
credible-set membership and probabilities are archived per instance. See the
[SuSiEx paper](https://www.nature.com/articles/s41588-024-01870-z) and
[SuSiE component-union implementation](https://github.com/stephenslab/susieR/blob/master/R/susie_get_functions.R).

The expanded rerun regenerates effects, phenotypes, causal manifests and seed records
**after** harmonized locus re-scoring. Effects are saved in float64. MAF metadata now
preserves every causal variant per site instead of overwriting earlier variants in a
multi-causal instance. Phenotype generation refuses inconsistent site/reference coding.
Validation independently reconstructs every phenotype in FID/IID order from current
genotypes, saved effects and deterministic noise, and checks every instance/site exactly
once. It also reports per-site/per-ancestry signal variance, noise variance and genetic
means, structural completeness, effect-draw correlation and variance-ratio calibration.
The established clean genotype cohorts are reused after reference-coding, identity,
composition and disjointness checks. Their underlying HAPNEST pool and sampling design
are unchanged inputs to this fixed-panel experiment.

The design retains 79 historically selected loci, re-scores them using harmonized
genotypes and assigns new within-panel divergence strata. It is a frozen benchmark panel,
not a new genome-wide random sample. All 15 architectures and ten replicates are included
in each of nine runs: pooled centralized, pooled federated, three sites separately,
Covenant statistics with same-ancestry ANL LD, and three 50,000-person federation samples.
Seeds 20260601, 20260602 and 20260603 match **analyzed** N after the small EAS column is
excluded. Remaining arm differences include composition, frequencies, effects and noise;
they do not isolate an LD-only effect. Borrowing ANL panels tests independent
same-ancestry references, not European LD substituted for African LD.

Rates and paired-arm differences use a bootstrap over whole loci, retaining replicates
and credible sets within each locus. Coverage means the fraction of returned credible
sets containing a true causal variant, distinct from per-instance power. Shared and
divergent causal modes had separate summaries and figure groups; the divergent mode has
since been withdrawn (see
[Ancestry-divergent architectures](#ancestry-divergent-architectures)). The effect-draw `rg`
parameter is distinguished from realized genetic-value correlation, and site-marginal
heritability from within-ancestry heritability. Homogeneous-noise, null-association and
LD-isolation studies would be needed for stronger general calibration or causal claims.
Implementation equality and variance-ratio checks do not establish those claims.

The PBS workflow is implemented by
[`submit_scientific_rerun.py`](../../../scripts/fine-mapping/submit_scientific_rerun.py)
and [`scientific_rerun.pbs`](../../../scripts/fine-mapping/scientific_rerun.pbs):

1. Harmonized locus re-scoring array and complete merge.
2. Fresh phenotype/effect simulation array and manifest/seed merge.
3. Full phenotype-validation array and acceptance merge.
4. Fresh bundles and DRS content identifiers; publish the corrected local configuration.
5. A 54-task analysis array, capped at three concurrent tasks.
6. Complete-grid, input-list and output parity gates, summaries and paired-arm estimates.
7. Full available EDA, centralized, federated, participation and exemplar plots. Figure
   exceptions fail the job. Missing optional inputs are documented in the plot manifest;
   genome-wide candidate-window plots do not apply to a frozen selected panel.

The figure stage harvests saved outputs, uses corrected union PIPs in locus plots,
separates causal modes and shows paired locus-bootstrap differences rather than a
post-hoc "best solo site" decomposition. It records `plots.accepted.json` with every
written figure. The stage does not silently accept a rendering exception.

Prepare and submit **after committing** the corrections, from the repository root:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/soft/applications/miniconda3/3.12/bin/python \
  scripts/fine-mapping/submit_scientific_rerun.py \
  --base-config local/configs/fine-mapping/pipeline_config.clean.yaml \
  --root local/output/fine-mapping-full-rerun-20260914 --account "$PBS_ACCOUNT" --submit
```

The source, executables, task configurations, subset manifests and package versions are
frozen and checked at job start. Source is exported using `git archive` from the committed
revision, so concurrent uncommitted work cannot enter the PBS snapshot.
`environment.json` records the committed revision;
`jobs.json` records the actual PBS IDs and commands. A repeated `--submit-only` invocation
reuses recorded job IDs instead of duplicating submissions. Set `PBS_ACCOUNT` to your
Sophia allocation; submission requires an explicit `--account`. These are centralized and standalone federated HPC computations. The two
placeholder partner endpoints still require actual deployment and corrected bundle
transfer before a remote federation can run.

The earlier chain `185778`–`185784` was cancelled when the user expanded the scope to fresh
simulation and the complete figure suite. Its diagnostics and cancellation record remain
under `local/output/fine-mapping-rerun-20260914/`; the replacement uses the separate
`local/output/fine-mapping-full-rerun-20260914/` directory. `analysis.accepted.json` and
`plots.accepted.json` mean their respective execution gates passed, not automatic
scientific sign-off. Inspect coverage, eligibility, nonconvergence and noise diagnostics
before making inferential claims.

Verification before submission: 216 tests pass, including regenerated phenotype-vector
reconstruction, causal-MAF completeness, signed allele harmonization, missing-variant
rejection, per-mode credible-set joins and the PBS dependency graph. Lint, shell syntax
and whitespace checks pass. A synthetic rendering smoke test produced 24 available
figures (11 inference, 5 federation and 8 EDA); those smoke figures are test fixtures.
The audited production case `L0002_ncsl1_h2-0.001_rg1_rep2` passes the explicit input-list
gate for all six ancestries and retains identical PIPs (maximum absolute difference 0)
and the same causal-containing credible set after explicit GWAS allele harmonization.
These checks do not replace the scheduled full-grid parity and calibration assessment.

Scheduler incident and recovery: on 14 September, replacement array `185829[]`
finished with two successful subjobs (1 and 3) and two terminated subjobs (0 and 2,
exit 143 / SIGTERM). PBS cancelled dependent jobs `185830`–`185838` at 20:35 UTC;
phenotype simulation, analysis and plotting had not started. The failed jobs used
about 2.2 GB and 25 minutes against requests of 120 GB and four hours. Their logs
contain no Python traceback. Both terminations coincided with successful completion
of another subjob on the same physical node (`sophia-gpu-07`). This is consistent
with interference during shared-node cleanup; the available logs do not identify
who sent SIGTERM, so that cause remains unconfirmed.

Submission now requests `place=exclhost` for every stage. It prevents separate jobs
from sharing a physical host, as defined by the [ALCF PBS placement documentation](https://docs.alcf.anl.gov/running-jobs/).
Sophia rejected an exclusive one-GPU request with `invalid ngpus= value`; the
accepted request reserves a full node (`ngpus=8:ncpus=256:mem=960gb`). The site hook
normalizes placement to `scatter:excl`. Jobs `185867`–`185876` were accepted at
20:46 UTC: re-scoring is queued for resources and subsequent stages are held by
their success dependencies. Whole-node reservation can increase queue time and
resource charges. This is a scheduling mitigation; it does not relax any scientific
acceptance gate.
The original frozen computation revision and data remain unchanged. To recover the
terminal dependency graph while retaining output and validated caches:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/soft/applications/miniconda3/3.12/bin/python \
  scripts/fine-mapping/submit_scientific_rerun.py \
  --root local/output/fine-mapping-full-rerun-20260914 --account "$PBS_ACCOUNT" --retry-failed
```

Retry verifies that every recorded job is terminal, archives its job IDs and PBS
history under `attempts/`, then submits a new graph. It refuses to duplicate queued,
held or running work. Re-scoring reuses its validated per-window caches. New job IDs
are recorded in `jobs.json`; prior job histories are retained for diagnosis.

## Outcome — 17 September 2026

The graph completed. Every execution gate passed: `bundles.accepted.json`,
`analysis.accepted.json` (nine runs), `parity.json` and `plots.accepted.json`. Gate
passage is not scientific sign-off; the findings below are the review those gates
were built to make possible.

### Exactness

All 11,850 paired instances agree. Zero mismatches, maximum absolute PIP difference
**0.0**, and every pair converged, so the count is not inflated by jointly-failed fits.
Descending to the credible-set members strengthens this: across **102,150** harvested
members the two paths agree exactly on membership, `cs_pip` and `ovrl_pip`. Twelve of the
fourteen shared-mode inference figures are byte-identical between the centralized and
federated trees, which is the same result read off the rendered output.

The claim's scope is PIPs, credible-set membership and fit status — not every emitted
number. Reviewing the per-ancestry statistics attached to those members found one
residual difference, confined to **AFR**: 2,628 of 55,275 variants (4.8%) carry a
different `beta_AFR`, `se_AFR` and `frq_AFR` between the paths. The bound is
2.44 × 10⁻⁴ absolute on beta (4.8 × 10⁻⁴ relative), 1.38 × 10⁻⁵ on the standard error,
and exactly one unit in the last emitted decimal on the frequency. It is essentially a fixed
per-variant property rather than scattered noise: 2,624 of the 2,628 differ at every
instance in which they appear and only four are inconsistent. The five other ancestries
agree exactly. AFR is the largest pooled column — 60,000 people across all three sites —
which is where accumulation order in the moment path has the most terms to disagree
about. None of it reaches inference: the PIPs and sets built from these
statistics are identical. The mechanism is consistent with summation order and is not
established by this review; what is established is the bound and the absence of
propagation.

### Simulation validity

Five of five validation checks pass. All **35,550** phenotype instances reconstruct from
current genotypes, saved effects and the deterministic noise stream, in FID/IID order;
worst relative heritability error 1.1 × 10⁻¹⁰ and worst phenotype-vector error
4.2 × 10⁻⁹ SD. Site heritability lands within 0.43% of target at the median and 2.9% at
the worst against a 20% tolerance. Effect-draw correlation recovers its `rg` parameter to
within 0.0087, 0.0006 and 1.0 × 10⁻¹⁰ at 0.5, 0.7 and 1.0. Fourteen structural invariants
hold. All 237 LD matrices are free of NaN and within tolerance, though every one of them
overshoots |r| = 1 slightly on rounding — maximum 1.0043 against a 1.01 floor, and a
maximum diagonal deviation of 4.3 × 10⁻³. The negative control fires as intended: the
wrong site-level model is materially off for 99.9% of instances at `rg` 0.5 and 0.7, and
for 2.4% at `rg` 1.0, where the two models coincide by construction.

### Credible-set coverage

Coverage is the fraction of **returned** credible sets containing a true causal variant,
against a nominal 95%. Pooled over both causal modes, with a locus-clustered bootstrap:

| arm | analyzed N | coverage | 95% CI | credible sets |
| --- | ---: | ---: | --- | ---: |
| all three, full | 150,000 | 0.9631 | 0.9594–0.9668 | 15,093 |
| MBZUAI alone | 50,000 | 0.9570 | 0.9531–0.9608 | 8,828 |
| Covenant alone | 50,000 | 0.9560 | 0.9516–0.9604 | 10,907 |
| Covenant stats + ANL LD | 50,000 | 0.9553 | 0.9511–0.9595 | 10,903 |
| ANL alone | 50,000 | 0.9472 | 0.9414–0.9525 | 8,876 |
| all three, n matched (3 draws) | 50,000 | 0.9461 | 0.9422–0.9498 | 21,543 |

Every arm is within 1.3 points of nominal, which is the headline. The deviations are
nonetheless resolvable at this number of sets, and they are not symmetric: four arms sit
**above** nominal and the pooled matched-N control sits **below** it, by 0.4 points. The
sign tracks power rather than participation — better-powered arms are conservative, the
weakest arm is slightly liberal — so the aggregate figure should not be read as evidence
that coverage is unconditionally calibrated. Individual matched-N draws span
0.9400–0.9515, a spread comparable to each draw's own interval, which is why the three
are pooled here rather than quoted separately.

The aggregate hides the one place coverage genuinely fails. Conditioning on `rg` = 1.0,
the only level of the star design that spans all three heritabilities:

| arm | h² = 0.0005 | h² = 0.001 | h² = 0.005 |
| --- | --- | --- | --- |
| all three, full | **0.949** (0.939–0.959) | 0.966 | 0.969 |
| Covenant alone | 0.938 (0.923–0.952) | 0.961 | 0.954 |
| MBZUAI alone | 0.901 (0.871–0.928) | 0.954 | 0.974 |
| ANL alone | 0.878 (0.853–0.903) | 0.950 | 0.967 |
| all three, n matched (draw 1) | **0.799** (0.755–0.842) | 0.946 | 0.973 |

At the weakest signal the full federation is the only arm whose 95% sets still cover at
95%. Every 50,000-person arm under-covers there, down to 80% for the matched-N control,
whose 293 sets in that cell come from a 9.9% discovery rate — the few sets a weak arm
does return are selected on having got lucky. This is a selection effect on which
instances produce a set at all, and it is the sharpest limitation in the run: **a 95%
credible set from a single site at h² = 0.0005 is not a 95% credible set.** Coverage
also declines from the low to the high LD-divergence stratum in every arm, but
`pap7_ld_divergence` establishes that this stratum tracks LD *abundance* rather than
cross-site disagreement, so the honest reading is that more LD gives wider sets and
slightly lower coverage, not that ancestral divergence breaks calibration.

### False positives

PIP thresholds are conservative everywhere. A PIP > *t* threshold should hold false
discovery below 1 − *t*; both are met with room to spare, in every arm:

| threshold | bound | observed range | full federation |
| --- | ---: | --- | ---: |
| PIP > 0.5 | 0.50 | 0.128–0.170 | 0.128 (0.118–0.138) |
| PIP > 0.95 | 0.05 | 0.018–0.035 | 0.019 (0.015–0.024) |

The full federation has the lowest false-discovery rate at both thresholds while
returning by far the most high-confidence variants (6,225 at PIP > 0.95 against
2,336–3,372 elsewhere), so the gain is not bought by committing less often. Note that
this denominator counts any non-causal variant as false, including one in complete LD
with the causal variant, which is unresolvable rather than wrong; the rates are therefore
upper bounds.

### Comparative performance

Power is the probability that a returned set captures a causal variant. Paired
locus-clustered bootstrap, full federation minus the arm, in percentage points:

| arm | power | full federation minus arm |
| --- | ---: | --- |
| all three, full | 0.9433 | — |
| Covenant alone | 0.7773 | +16.6 (15.8–17.4) |
| Covenant stats + ANL LD | 0.7769 | +16.6 (15.9–17.4) |
| MBZUAI alone | 0.6257 | +31.8 (30.8–32.7) |
| ANL alone | 0.6224 | +32.1 (31.2–32.9) |
| all three, n matched (3 draws) | 0.489–0.497 | +44.7 to +45.4 |

Federating helps on every axis at once: power 0.94 against 0.62–0.78, causal recall 0.61
against 0.35–0.44, mean credible-set size 6.8 variants against 9.5–12.7, mean purity
0.981 against 0.937–0.955, and 1.27 sets per instance against 0.75–0.92. It finds more
signals and localizes each one more tightly.

**The power gain is sample size, not diversity — but diversity is not idle.** The two act
at different stages and the marginal figures conflate them.

At the discovery stage the matched-N control points against diversity, decisively. Hold
analyzed N at 50,000 and the three-site composition is the *worst* arm in the run, below
every solo site, reproducibly across three independent draws. SuSiEx fits one effect per
ancestry column, so power tracks the largest column rather than the total, and the solo
arms win by concentration — Covenant puts 47,500 of its 50,000 people into AFR, where the
matched-N federation's largest column is 20,339. The controlled version of that comparison
holds both N and column count fixed, ANL (5 columns, largest 30,000) against matched-N
draw 1 (5 columns, largest 20,339): **13.3 pp apart (12.3–14.3) in ANL's favour.**

At the localization stage it reverses. Restricting to the **4,987** instances where both
arms returned a set — which removes the differential-discovery selection that makes a
marginal resolution comparison unreadable — the diverse arm's best credible set is
**2.16 variants smaller (95% CI 1.75–2.58)**, 8.50 against 11.05, with higher mean purity
(0.963 against 0.954) and higher causal PIP (0.583 against 0.530). Capture rate within
that subset is essentially equal (0.965 against 0.970), so this is sharper localization of
the same signals, not a different set of them. That is the canonical cross-ancestry
fine-mapping benefit — differing LD breaking ties between correlated variants — and it is
present here.

So at fixed N, ancestral composition trades discovery for resolution: fewer signals found,
each one pinned down better. The full federation escapes the trade because it is not at
fixed N — 150,000 people buy both, and it leads every arm on power *and* on set size
(6.8 variants). Because matching N does not match allele frequencies, realized effects or
site noise, the 2.16-variant gap is consistent with the LD-diversity mechanism rather than
an isolated measurement of it. Labelling either bar simply "diversity" would misstate the
finding in one direction or the other.

Borrowing a same-ancestry external LD panel is, on this design, close to free: Covenant
statistics against ANL's AFR panel reach power 0.7769 against Covenant's own 0.7773 and
coverage 0.9553 against 0.9560. This is the cheap alternative to an O(M²) uplink and it
does not visibly break here. The panels are drawn from the same synthetic HAPNEST pool
for the same ancestries, so this bounds the cost of *external-panel sampling*, not the
cost of a genuinely mismatched reference. It is excluded from the exactness comparison by
construction.

### Limitations

- **Ancestry-divergent architectures are not evaluated, and the mode is now off.** See
  [Ancestry-divergent architectures](#ancestry-divergent-architectures) below. Nothing in
  this run speaks to ancestry-private causal variants.
- **Coverage is conditional on power.** See the h² table above. Nothing here licenses a
  95% claim for a single site at the weakest signal.
- **`rg` is the effect-draw correlation**, not realized genetic-value correlation, and the
  grid is a star rather than a full factorial: `rg` 0.5 and 0.7 exist only at h² = 0.001.
  Any marginal over `rg` is confounded with h², and pooling them reverses the sign of the
  effect — within h² = 0.001 power rises with `rg` in every arm (full federation
  0.969 → 0.974 → 0.978; ANL 0.581 → 0.597 → 0.719). The published rollups now carry
  h² inside the `rg` grouping for this reason.
- **Arm differences are composition and participation effects.** Matching analyzed N does
  not match allele frequencies, realized effects or site noise, so no arm difference
  isolates an LD-only cause.
- **79 frozen loci on synthetic HAPNEST chr1.** Locus-clustered intervals treat these as
  the locus sample; they are a fixed benchmark panel, not a genome-wide draw, and the
  11,850 rows are not 11,850 independent regions.
- **No homogeneous-noise or null-association arm**, so calibration of the association
  model itself is untested. Implementation equality and variance-ratio checks do not
  supply it.
- **Still centralized and standalone-federated HPC computation.** The two placeholder
  partner endpoints remain undeployed.

### Ancestry-divergent architectures

Withdrawn. The mode simulated a causal set that differs by ancestry — a union of shared
variants plus one private to each superpopulation — as a stress test of the assumption the
whole method rests on, that the causal variant is the same in every ancestry. The question
is the right one and the compute was negligible: 2.7 core-hours, 0.23% of the run. The
design could not answer it, for two independent reasons.

It was far too small. `min_per_stratum: 1` selected three (locus, architecture) pairs, so
the mode was 30 of 11,850 instances per arm, yielding 3 to 31 credible sets. Per-arm
coverage ranged from 0.75 to 1.00 and recall from 0.014 to 0.138; at those denominators
neither is an estimate.

More instances would not have rescued it, because the comparison is confounded. A
divergent instance draws a union of `ncsl - n_private` shared variants plus
`n_private_per_pop` for each of six superpopulations, while keeping the architecture's
target h². An architecture **labelled** `ncsl2` therefore carried **seven** causal
variants where its shared counterpart carried two, making each effect roughly 3.5 times
fainter, and the shared grid spans `ncsl` 1–3 only, so there is no seven-causal cell to
compare against. The observed recall drop — 0.615 shared against 0.138 divergent for the
full federation, with every arm falling 27 to 48 points — mixes ancestry-private
architecture with much weaker per-variant effects and cannot separate them.

`ancestry_divergent_causal` is consequently `enabled: false` in the shipped configuration,
and no future run emits a divergent figure group. The implementation is retained, and
`_assign_ancestry_divergent_flags` now warns once at simulation setup whenever the grid
lacks a shared cell at the union size — which is exactly the condition that made this run's
contrast unreadable. Re-enabling it responsibly needs both a larger `min_per_stratum` and
a matched `ncsl` cell equal to `n_shared + 6 x n_private_per_pop`.

The completed 14 September run is left intact as a frozen artifact: its divergent figures
and `plots.accepted.json` still list what that run produced, with a `NOT-EVALUATED.md`
marker placed in each divergent figure directory. Nothing in this repository cites a
divergent number as a result.

### Corrections applied during this review

The estimates above are the run's own; nothing was recomputed. Four reporting defects
found while checking figure labels, denominators and uncertainty against the results were
fixed, and the figures and stratum rollups re-rendered from the frozen results at the
revision that carries the fix:

1. **Four captions claimed Wilson intervals.** Audit finding 9 replaced row-independent
   intervals with a locus-clustered bootstrap; `clustered_ratio` was adopted but the
   prose in `pap1`, `pap2`, `pap4` and `fed6` was not updated. The drawn intervals were
   always the clustered ones, so no number changes. The unused `wilson`/`binom_summary`
   helpers are deleted so the assumption cannot return.
2. **`pap1` reported the wrong denominator.** Each panel plots one arm of the star design
   but the footnote printed the whole frame: 11,820 instances against the 7,100 and 7,090
   actually drawn, and 30 against 10 and 20 in divergent mode.
3. **The stratum rollup pooled h² inside the `rg` marginal**, putting 2,340 instances in
   the `rg` = 1.0 row against 780 in the others and reversing the apparent direction of
   the `rg` effect. `h2_target` is now part of that grouping in all three writers.
4. **Two matched-N arms reached a published figure as raw identifiers.** The draws were
   built inline by the submit script and labelled in the figure module, with no shared
   source of truth, so `federation_50k_seed2` and `federation_50k_seed3` printed beside
   the other arms' prose labels. `MATCHED_N_SEEDS` is now the single definition both read.
   The same panel's cohort-size and column-count rows were positioned for four upright
   tick labels and were overrun by eight rotated ones; past four arms those numbers move
   to the footnote, where they are legible and carry the prose labels too. The stale
   "14.6 pp" in `fed6`'s docstring is corrected to the measured 13.3 pp.

`figures.before-caption-fix/` and `*.before-h2-split` retain the published versions.
Three regression tests pin all four.
