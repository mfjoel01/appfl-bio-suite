# MultiSuSiE coordinator pipeline

## Current status — 2 November 2026

The full pipeline is complete: 240,900 fits across 2,281 tasks, plus 11,850
independently formed centralized comparison fits. All comparisons passed numerical
PIP tolerance; one exact credible-set membership mismatch is retained and reported.
The editable poster, print PDFs and visual audit are complete. No jobs remain active.

The subsequent matched-filter comparison is also complete: 139,650 paired instances
across 931 tasks. See [the comparison design](MATCHED_COMPARISON.md) and
[results with paired confidence intervals](MATCHED_RESULTS.md). With common purity
filtering, the full-cohort result is modestly higher recall for MultiSuSiE versus
better-calibrated, smaller sets and lower FDR for SuSiEx. Rankings depend on cohort
and filter policy; these results do not establish universal superiority.

All local artifacts now reside in the primary checkout. See
[artifact locations and archive verification](MULTISUSIE_ARTIFACTS.md).
The original poster reports the MultiSuSiE pipeline; the matched comparison report
is a separate deliverable. MultiSuSiE is now selectable in the remote APPFL
aggregator as well as the standalone coordinator pipeline. Default SuSiEx behavior
is unchanged. The completed large-grid studies used locally staged synthetic sites;
they are not evidence of a live remote MultiSuSiE deployment.

Before remote-backend integration, targeted validation was 46 passed, 2 skipped; lint, formatting and
reusability checks passed. The skipped tests require optional executables. Real-backend
paired comparison and checkpoint reuse passed. Historical test counts below refer
to their dated runs.

## Model and inputs

The estimator is [MultiSuSiE](https://github.com/jordanero/MultiSuSiE), pinned to
`10351ebbbc0ed4c891c455000eb6ec7800d0b717`. Its `susie_multi_ss` sufficient-statistic API
accepts the centered second moments available from the existing site protocol. The
adapter calls this API directly, avoiding rounded GWAS text and LD serialization:

```
XtX = sum(G) - outer(sum(u), sum(u)) / sum(n)
Xty = sum(c) - sum(u) * sum(q) / sum(n)
yty = sum(w) - sum(q)**2 / sum(n)
```

All sums and centering are within ancestry. The fit uses raw dosage effect sizes,
float64 arithmetic, estimated residual and prior variances, ten signals, 1,000 maximum
iterations, convergence tolerance 1e-6, and 95% credible sets. Genotypes are not scaled
to unit variance. This is the existing linear additive phenotype experiment, without
additional covariates.

`rho=0.75` is an explicit working prior, matching the upstream RSS default. It is never
read from the simulation's `rg`, effect sizes or causal manifest. Sensitivity runs should
use independently specified values such as 0, 0.5, 0.75 and 0.95, with separate output
directories. Exactly 1 is rejected because the correlation matrix would be singular.

The adapter applies pooled per-ancestry MAF >= 0.005, as in our SuSiEx protocol, and
excludes monomorphic variants. MAC defaults to zero rather than upstream RSS's 20 to
avoid silently changing eligibility; `--min-mac 20` is available as a separate sensitivity
analysis. Variants retained in any ancestry enter the union. Ineligible ancestry/variant
combinations have zero rows/columns and zero cross-products, rather than discarding the
variant from every ancestry. Within an ancestry, the existing complete-variant policy
and allele harmonization are retained. Mismatched or duplicate site contributions fail.

## Switch methods without changing the workflow

Both methods use the same `run fine-mapping` command, site payloads, sharding,
driver choices and reporting paths. Override just the backend; use distinct results
directories to preserve each fit:

```bash
appfl-bio-suite run fine-mapping --config loopback --driver serial \
  --data-root /path/to/simulation --backend susiex --results-dir local/output/run-susiex
appfl-bio-suite run fine-mapping --config loopback --driver serial \
  --data-root /path/to/simulation --backend multisusie --results-dir local/output/run-multisusie
```

For remote runs, keep your federation and driver unchanged:

```bash
appfl-bio-suite run fine-mapping --federation local/federation.yaml \
  --driver globus_compute --backend multisusie --results-dir local/output/remote-multisusie
```

`--backend` overrides `experiments.fine-mapping.fine_mapping_backend` for that run.
`--results-dir` overrides `fine_mapping_output_dir`; `--out-dir` still selects where
resolved launch configurations are written, not inference results. These two new flags
apply only to fine-mapping. The SuSiEx default remains unchanged when neither is supplied.

Both backends write the same `data/fed_fm_results.tsv` column schema, rollup paths,
`data/inference_backend.json`, and per-instance `data/instances/<instance>/fit.json`,
`variants.tsv`, and `credible_sets.tsv`. The common posterior tables expose `snp_id`
and `pip`, and component/set membership with component PIP and coverage. Native
backend details are retained additionally: MultiSuSiE posterior arrays in the instance
directory, SuSiEx compressed outputs in `artifacts/<instance>/`. Existing instance
outputs cannot be silently overwritten. Both backends reject unsuccessful fits before
reporting the exchange ready. With unknown truth, both leave truth-based metrics missing.

This is an interface swap, not identical statistical behavior. PIP definitions, priors,
precision and filters remain method-specific and are recorded. MultiSuSiE's per-ancestry
minimum-p-value and CS-purity summaries are currently missing in the shared table;
SuSiEx populates its own native diagnostics. Neither missing value means zero.
The full matched-filter study remains the appropriate comparison for estimator claims.

## Select MultiSuSiE for an APPFL federation

Install the optional backend on the coordinator from this checkout:

```bash
python -m pip install '.[multisusie]'
```

Set these fields in your existing federation YAML (keep your site/endpoint settings):

```yaml
experiments:
  fine-mapping:
    fine_mapping_backend: multisusie
    multisusie_rho: 0.75
    multisusie_min_abs_corr: 0.0
    multisusie_min_mac: 0.0
    n_signals: 10
    max_iter: 1000
    tol: 1.0e-6
    credible_set_level: 0.95
    maf: 0.005
```

```bash
appfl-bio-suite preflight --experiment fine-mapping --check env --federation local/federation.yaml
appfl-bio-suite run fine-mapping --driver globus_compute --federation local/federation.yaml
```

The launcher forwards backend selection to `FineMappingAggregator`. Sites run the
same shipped trainer and send the same raw-moment tensors; partners do not install
MultiSuSiE. Backend-aware preflight requires the optional Python backend, not SuSiEx
or PLINK, for this inference path. Simulation preparation may still need PLINK.
Select `fine_mapping_backend: susiex` (the default) to retain the original path.

Shared settings are `n_signals`, `max_iter`, `tol`, `credible_set_level`, and `maf`.
`pval_thresh` and `keep_ambiguous` are SuSiEx-specific; MultiSuSiE does not apply a
significance filter. Its purity/MAC/rho settings are the prefixed fields above.
The remote path does not silently retry with a higher iteration ceiling: a
nonconverged fit retains diagnostic files and raises before reporting a ready result.
Use a fresh output directory for a rerun; existing native fit directories are never
overwritten. `n_workers` in the server configuration controls instance parallelism.

Results retain the standard `data/fed_fm_results.tsv` and rollups. Backend/options
are recorded in `data/inference_backend.json`; native MultiSuSiE PIPs, credible sets,
posterior arrays and backend source/version provenance are retained under
`data/instances/<instance>/`, including when `keep_work` is false. Native PIPs are
not redefined by CS filtering. Without a causal manifest, truth-based metrics are
missing rather than false or zero. CS purity summaries are currently unreported
(missing) in the remote result table; the native CS membership remains available.

Regression tests use actual shipped trainer payloads and compare decoded-aggregate
fits with independently formed standalone moments, PIPs, sets and posterior arrays,
including serial/parallel execution. They also check configuration propagation,
backend-aware preflight and failure on nonconvergence. These are local protocol
integration tests; no live partner endpoint was contacted for this validation.
All 12 backend integration tests passed, including CLI selection for both drivers,
identical result/artifact schemas, no-truth behavior and overwrite protection.
The original SuSiEx APPFL-versus-standalone parity, configuration, preflight and
documentation regression checks also passed.

## Run a pilot

Install into a coordinator environment with Python 3.12:

```bash
python -m pip install '.[multisusie]'
```

Use an existing **synthetic** simulation's `pipeline_config.yaml`. The standalone pilot
reads each site's staged data using the existing `site_aggregates` implementation;
inference receives only those aggregates. This command does not dispatch remote tasks.
Use compute-node resources for production-sized genotype windows.

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python -m appfl_bio_suite.experiments.fine_mapping.multisusie_experiment \
  --config /absolute/path/to/pipeline_config.yaml \
  --out local/output/multisusie-pilot-rho075 \
  --limit 2 --rho 0.75
```

`--locus L0000` optionally restricts the selected locus. The default limit is one
instance. Selection is sorted by locus, architecture and replicate: this is a plumbing
pilot, not a representative scientific benchmark. The output directory must be new.
Existing full-federation, solo-site or matched-N arm configurations can be supplied;
the planned ancestry counts are checked before inference.

Each instance writes:

- `variants.tsv`: native MultiSuSiE PIPs and ancestry-specific eligibility flags.
- `credible_sets.tsv`: included sets, component probabilities and claimed coverage.
- `posterior.npz`: all component probabilities, fitted effects and variance estimates.
- `fit.json`: options, population sizes, moment fingerprint, backend source hash/version,
  convergence and elapsed time, including diagnostic failures.

The pilot's `results.tsv` scores truth only after fitting. Nonconverged fits retain
diagnostic output but receive no discovery/coverage scores and cause a nonzero exit.
Native PIPs include active components before CS filtering; they must not be relabeled
as SuSiEx's retained-component PIPs. Default `--min-abs-corr 0` follows the MultiSuSiE
RSS interface. An optional purity cutoff filters sets, not the native PIPs.

## Validation

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python -m pytest tests/test_fine_mapping_multisusie.py -q
```

Tests use the real backend when installed: federated moments versus independently
centered centralized data, identical PIPs/CSs within numerical tolerance, ancestry-only
censoring, malformed inputs, nonconvergence, and the existing miniature site filesets.
Backend-dependent tests skip when the optional dependency is absent.

## Historical execution record

The following dated notes preserve how the run and acceptance criteria evolved.
Job IDs, old worktree paths, pending states and test counts are historical; use the
current status and artifact guide above for the completed run.

Initial validation on 18 September 2026: all nine new tests passed. The combined
MultiSuSiE, federated and scientific regression files passed 41 tests, with 12 skipped
because optional legacy/runtime dependencies were unavailable. A saved CLI smoke run
at `local/output/multisusie-smoke/rho075` used three synthetic sites, three ancestries,
60 variants, and two phenotype instances, with the default ten signals. Both fits
converged and their returned sets captured all simulated causal variants. This validates
execution only; two high-signal instances do not establish calibration or superiority.

## Paired comparison on the existing cohort

`multisusie_compare` runs current SuSiEx, SuSiEx-unfiltered and MultiSuSiE on an explicit
TSV panel with `locus_id`, `architecture_id`, and `replicate` columns. The SuSiEx wrapper
now accepts an optional `min_purity` setting; calls omitting it keep their old behavior.
Both SuSiEx fits use ten signals and disable multistep fitting. The unfiltered fit sets
`min_purity=0` and `pval_thresh=1`; the current fit uses 0.5 and 1e-5 respectively.

The independent centralized check reads genotypes directly from the original HAPNEST
fileset, selects the enrolled individuals using FID/IID, aligns pooled phenotypes, and
forms ancestry moments over those rows. It does not add site Gram matrices. The gate
requires matching variant metadata, ancestry eligibility, sample sizes, fit settings,
convergence and CS membership, with PIPs agreeing at `atol=1e-7`, `rtol=1e-5`.

```bash
PYTHONPATH=src python -m appfl_bio_suite.experiments.fine_mapping.multisusie_compare \
  --config /path/to/pipeline_config.yaml --panel /path/to/panel.tsv \
  --out local/output/paired-pilot --susiex /path/to/SuSiEx
```

`multisusie_batch prepare` freezes source, the binary, configurations, and a prescribed
panel from the completed September run. The panel takes the median-divergence locus
within each of low/medium/high strata, four specified causal architectures, and replicate
zero. It never selects loci based on earlier inference outcomes. The three configurations
are full federation, the existing 50,000-person federation draw, and Covenant alone.
There are 36 cases, 108 method fits, and 12 additional centralized checks for the full
federation. The other arms use the validated federated path without an independent
centralized rerun. This small, deterministic panel cannot establish general superiority.

```bash
PYTHONPATH=src python -m appfl_bio_suite.experiments.fine_mapping.multisusie_batch prepare \
  --base /path/to/fine-mapping-full-rerun-20260914 \
  --root local/output/multisusie-comparison --checkout "$PWD" --susiex /path/to/SuSiEx
```

Submit the generated `job.pbs` with an allocation/resource request and **absolute**
`FM_ROOT` and `FM_PYTHON` environment variables. It runs the three arms in sequence and
then validates and merges them. The script sets eight BLAS threads; request at least
eight CPUs and enough memory for the full source genotype window. No GPU computation is
used. Source, binary, configuration, selected phenotype, sample metadata, and package
checks are repeated before each arm. Large BED files are checked by size and modification
time rather than content hashes; this is a drift check, not a cryptographic data identity.

The final gate rejects missing/duplicated instance-method pairs, failed or nonconverged
fits, missing arms and incomplete centralized parity evidence. Successful runs write
`comparison.tsv`, descriptive `summary.tsv`, and `accepted.json`. CS coverage and
variant-level false discoveries are distinct metrics. `site_stage_s` is a locus-wide
cost repeated on its instance rows and must not be summed as independent per-fit cost.

On 18 September 2026 the frozen pilot was submitted as PBS job
`186940`, under
`local/output/multisusie-comparison-20260918`. Submission is not a successful result;
inspect `accepted.json`, `summary.tsv`, `federation/parity.json`, and the job log before
interpreting it. Its selected loci are L0037 (low), L0021 (medium), and L0043 (high).

After adding the comparison and batch gates, 71 targeted tests passed with the existing
SuSiEx/PLINK binaries on PATH, including the original centralized/federated tests and
the real three-method comparison. There were no skips in that run. At that September check,
the production pilot remained queued for production-tagged nodes. The completed
full-grid and matched results linked above supersede that preliminary status.

In this historical pilot, MultiSuSiE used native moments while SuSiEx retained six-digit
summary statistics and float32 LD protocol. The variant union must agree between fits;
the numerical representations differ. An RSS-versus-moment sensitivity check remains
necessary before attributing small differences solely to the model. The completed
matched comparison uses 17-digit summary statistics and reports full-grid calibration;
broader rho sensitivity remains future work. Remote APPFL backend selection was
added on 2 November; see the current configuration instructions below.

## Full-grid continuation (23 October 2026)

The September pilot failed in upstream `susie_multi_ss` for a zero-credible-set fit:
that branch returns a tuple, and optional variant-label handling called `.append`.
The adapter now omits the backend's optional labels and maps returned integer indices
through its own validated variant metadata. The pinned backend is unchanged. A real
backend null-signal regression test covers this failure.

The fresh pilot is `local/output/multisusie-comparison-20261023`, PBS `198219`.
The original failed run is preserved. Submission does not establish acceptance.

`multisusie_pipeline` adds a checkpointed full-grid run at
`local/output/multisusie-full-20261023`:

- Benchmark: seven arms, 79 loci, 150 instances/locus: 82,950 fits, plus 11,850
  independently formed centralized full-federation fits for parity.
- Development: all 20 prespecified candidates, with three sampling seeds for mixed
  cohorts and one for each solo cohort: 101,250 fits on the 25 development loci.
- Selection: rerank only the MultiSuSiE development results using the original
  coverage/FDR screen and power/tie-breaking rule. No prior-estimator winner is assumed.
- Evaluation: 54 held-out loci and 8,100 instances per cohort; the selected cohort
  (three draws if mixed, one if solo) and the four baseline cohorts.
- Finish: validate exact expected instance keys, convergence and benchmark parity;
  summarize locus-bootstrap intervals and build the editable/print poster.

Each locus saves input phenotype hashes, fit artifacts, cumulative `results.tsv`,
parity evidence, and a checksummed acceptance record. Partial or failed fits are
archived before a retry. Converged fits and accepted loci can be reused on restart;
changed input phenotypes and tampered accepted result tables are rejected. Jobs
continue through independent loci and report all task failures. Failed or incomplete
phases cannot enter selection or the final poster. Code, configs and backend source
are frozen in the run; large genotype BEDs use size/mtime identity checks.

Historical PBS chain (full IDs and submission commands are in `submission.*.json`):

| Stage | Job |
| --- | --- |
| Benchmark, eight shards | `198235[]` |
| Development, eight shards | `198236[]` |
| Development selection | `198238` |
| Held-out evaluation, eight shards | `198239[]` |
| Validation, summaries, poster build/export | `198243` |

Historical worktree status command (do not use against the relocated archive; see
[MULTISUSIE_ARTIFACTS.md](MULTISUSIE_ARTIFACTS.md)):

```bash
local/venv-multisusie/bin/python scripts/fine-mapping/multisusie_status.py \
  --root local/output/multisusie-full-20261023
```

The copied poster was built in the worktree at `local/scratch/poster`
(now `local/scratch/poster-multisusie` in the primary checkout). Its original SuSiEx outputs and audits
are preserved in `reference-native`; these are not MultiSuSiE evidence. New outputs
are `native/multisusie-fine-mapping-poster.pptx`, the half-size PDF, a full-size
`*-print-72x36.pdf`, and a PNG preview. `multisusie_data.py` requires accepted complete
MultiSuSiE grids. `finalize_multisusie.py` verifies export paragraph preservation,
finite chart caches, embedded fonts, and print dimensions. At this stage, a final visual review was
still required; it passed on 29 October. `layout-check/` is an explicitly
watermarked mechanical test using the original data, never a result poster.

Validation on this continuation: 48 targeted tests passed (47 in the combined run
and one additional held-out-selection test), 12 optional-dependency skips; this includes actual-backend no-CS handling, checkpoint restart/tampering,
and the three-method miniature comparison. Full-grid jobs were still running/queued
at the time this section was written; inspect acceptance files before citing findings.
The execution remains local synthetic site stages and aggregate-only coordinator
inference. That historical run did not exercise a remote APPFL MultiSuSiE backend.

### Recovery on 24 October

The three-method comparison completed and passed all 108 fits. The full-grid arrays
stopped before completion: several subjobs reached 12 hours; one worker pool died;
two recorded fits reached the 1,000-iteration ceiling without convergence. At recovery
planning there were 189 accepted benchmark arm/loci and 459 accepted development
arm/loci. All successful fits and immutable source/input fingerprints are preserved.

The recovery overlay is `scripts/fine-mapping/recover_multisusie.py`, frozen separately
under `recovery-20261024`. It imports the original frozen implementation, leaving that
snapshot unchanged. Remaining work is split into 31 benchmark batches (up to 12 loci,
four workers) and 28 development batches (up to 32 loci, twelve workers). Benchmark
concurrency is halved to reduce memory pressure; the original worker crash is not
assumed to be proven out-of-memory. Evaluation uses 64 small shards.

The two capped fits first run in a diagnostic job with a 10,000-iteration ceiling.
Their moment fingerprints must equal the originals. The model, rho, precision,
initialization, variant filters, and tolerance remain unchanged. Recovery retries a
new nonconverged 1,000-iteration fit at 10,000, preserving the failed attempt and the
ELBO trace. A fit still unconverged at 10,000 continues to fail acceptance. This policy
is based on convergence only, never truth or discovery metrics.

Pooled/federated parity permits the iteration ceiling alone to differ when both fits
converged. Model options, sample sizes, eligibility, exact CS membership, and existing
PIP tolerances are still required. This is the same stopping rule with more time for
fits that need it, not a relaxed numerical convergence criterion.

Jobs: diagnostic `198503`, benchmark `198504[]`, development `198505[]`, selection
`198506`, evaluation `198507[]`, final validation/export `198508`. The recovery arrays
are conditional on successful diagnostics. See the latest recovery records rather
than the superseded October 23 scheduler chain. Two targeted recovery tests passed:
only the ceiling changes on retry; parity still rejects changed priors and lack of
convergence.

A separate lightweight manager (`watch_multisusie_recovery.py`) records state in
`recovery-20261024/manager.status.json`. It retries interrupted/incomplete batches
up to three rounds with smaller batches and fewer workers, replacing downstream
scheduler-canceled jobs after prerequisites pass. It stops for inspection on persistent
10,000-iteration nonconvergence or executed selection/export failures. This is scheduler
recovery only; acceptance rules are not relaxed. A third targeted test verifies reduced
batch sizing and removal of failed dependencies on retry.


The initial diagnostic converged at 2,139 iterations but correctly failed its exact
moment fingerprint check: single-phenotype staging differed from the original full
batch. Those outputs and submissions are preserved under `recovery-20261024`.
The corrected overlay under `recovery-20261024-r2` reconstructs the complete original
batch and verifies exact moments before fitting. Its current jobs are diagnostic
`198511`, benchmark `198512[]`, development `198513[]`, selection `198514`, evaluation
`198515[]`, and export `198516`. The first diagnostic input fingerprint has passed.
Five recovery tests pass, covering retries, strict parity, smaller scheduler retries,
no-eligible selection, and actual backend reconstruction of batched moments.

If no development candidate passes the prespecified screen, the selection records
`name: null` and the held-out evaluation runs the four baseline cohorts only. It does
not force an ineligible winner or change the thresholds. The poster reports “None
eligible” and omits the selected-cohort arm in this case. The manager status at that stage was
`recovery-20261024-r2/manager.status.json`.


## 28 October continuation

Development is complete and accepted: 101,250 fits across 1,350 tasks.
Selection was rerun after a premature scheduler dependency released it before retry
completion. The prespecified winner is `csa_to_anl_afr_1000` (49,000 Covenant + 1,000
ANL; 48,500 AFR + 1,500 EUR overall). Development coverage is 0.945057 and FDR is
0.012837. Evaluation did not influence selection; selection.lock.json freezes it.

Jobs at the 28 October checkpoint: benchmark `199674[]` (19 remaining loci, one per subjob), held-out
evaluation `199675[]` (378 cohort/loci, 56,700 fits), final export `199676`.
The monitor at that checkpoint was `recovery-20261028/manager-v2.py`; its status is in that
recovery directory. Eight recovery tests passed. The monitor can now retry a
premature downstream job only when its exact missing task acceptance file exists;
other executed failures still require inspection.

### Scientific finding and revised completion criterion

The original strict parity gate failed at `L0000_ncsl2_h2-0.005_rg1_rep4`: both
fits converged in eight iterations with identical options, aligned variants and
sample sizes, and maximum PIP difference 5.3013e-15, but CS boundary membership
differed. Original parity evidence is archived under recovery-20261028/original-parity;
known-cs-mismatch.json records both sets. Raw fits and CS exports are unchanged.

Completion now requires converged numerical PIP parity at the original tolerances,
while exact CS membership is a separately reported endpoint. This is an explicit
change from the original joint exact-CS acceptance criterion, not a claim that the
original strict test passed. New comparisons record strict_pass and
credible_sets_exact separately. The poster removes “Lossless” and universal CS
agreement claims and reports the number of CS-mismatched pairs. Prior accepted
comparisons already passed the stricter criterion. No credible sets were edited to
force agreement. Final data and poster acceptance were still pending at this checkpoint; both
passed on 29 October.


### Final-locus continuation on 29 October

All benchmark tasks passed collection (82,950 fits; 11,850 pooled comparisons).
The evaluation retry loop exhausted scheduler retries with one locus left,
`federation/L0063`, after an exit-143 termination. Its 115 converged fit checkpoints
were preserved. `parallel_multisusie_locus.py` stages the original complete phenotype
batch, forks eight workers that each own a distinct fit directory, and delegates
final scoring and acceptance to the original frozen `work_item`. Phenotype hashes,
selection hashes, model options, and convergence checks remain in force. A regression
test compares serial/parallel moment hashes and every saved posterior array exactly.
The recovery and parallel tests pass (nine tests). Job 199917 completed all remaining
fits; job 199918 follows its validation to collect results and export the poster.


### Completed 29 October

All 2,281 main-grid tasks passed: 82,950 benchmark, 101,250 development and 56,700
held-out fits (240,900 total), plus 11,850 centralized benchmark comparison fits.
All pooled comparisons passed numerical PIP tolerance; one exact CS membership
mismatch is reported. The separate comparison pilot also remains accepted.
Final jobs 199917, 199918 and 199925 exited zero; no MultiSuSiE jobs remain active.
The final PPTX, standard and 72x36-inch PDFs, preview, integrity audit and completed
visual-review record are now accessible from the primary checkout at
`local/scratch/poster-multisusie/native`. The final results panel
was inspected at higher resolution and its selected-cohort annotation corrected.
completion.json and RUN_STATUS.md in the run root are the final status records.

PBS submission requires your allocation: pass `--account` to
`submit_multisusie_pipeline.py`, or set `FM_ACCOUNT` for recovery submission.
