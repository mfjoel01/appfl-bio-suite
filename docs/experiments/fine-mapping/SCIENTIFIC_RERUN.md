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
divergent causal modes have separate summaries and figure groups. The effect-draw `rg`
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
