# Full matched-filter MultiSuSiE–SuSiEx comparison

This run addresses the output-filter confounding in the original comparison.
It does not presuppose which estimator will win.

Run root: `local/output/matched-multisusie-susiex-20261029`.

## Fixed design

Reuse every accepted instance from the 79-locus benchmark and the held-out evaluation:
82,950 + 56,700 = **139,650 paired instances**, across 931 cohort/locus tasks.
The development-selected cohort remains fixed for both estimators; no new cohort
selection uses these comparison outcomes. The existing MultiSuSiE posterior files
and original poster remain unchanged.

SuSiEx is run on the same people, phenotypes, ancestry partitions, eligible variants
and pooled moments, with its significance filter disabled (`pval_thresh=1`) and
internal CS-purity filter disabled (`min_purity=0`). Both methods use ten signals,
95% CS coverage, tolerance 1e-6 and the same convergence-retry ceiling policy
(1,000 then 10,000 iterations). Nonconverged fits cannot pass acceptance.

Two output policies are applied to both estimators:

1. Unfiltered native credible sets.
2. The same post-hoc minimum absolute LD-correlation threshold of 0.5. For a set,
   purity is the maximum, over ancestries in which all set members are eligible,
   of the minimum absolute within-set correlation. Both methods use the same
   float64 centered Gram matrices to calculate this value.

Native PIPs are unchanged by the post-hoc filter. This avoids silently giving one
method a PIP definition based on filtered components while treating the other
method's PIPs as unfiltered. SuSiEx still emits its own native active/retained
components; no attempt is made to invent unreported posterior components.

The native models and priors remain distinct: MultiSuSiE rho is fixed at 0.75.
SuSiEx receives 17-significant-digit summary statistics and its required float32 LD;
MultiSuSiE uses native float64 raw moments. These implementation differences are
recorded in design.json. This is a comparison of the configured estimators under
matched external filters, not a claim of identical mathematical models.

## Analysis and acceptance

Report power (any causal capture), causal recall, CS coverage, mean CS size,
FDR at PIP >0.95, and high-confidence causal yield. Report each method and paired
MultiSuSiE-minus-SuSiEx differences with 2,000 jointly resampled locus-bootstrap
replicates (fixed seed 20261029). Also stratify by architecture and causal mode.
Undefined ratios/intervals are left missing, never replaced by zero.

Recall must be interpreted together with calibration, FDR and set size. Nominal
95% coverage is a calibration target, not merely a quantity to maximize. No winner
is inferred from power alone, and no hyperparameter/filter is tuned on evaluation
results. Archived MultiSuSiE runtime is retained as provenance but is explicitly
marked incomparable; it cannot establish speed superiority.

Each accepted task contains every expected instance × two methods × two policies,
with convergence checks and checksums. Per-instance checkpoints resume interrupted
jobs. SuSiEx outputs and MultiSuSiE source files used for scoring are hashed.

## Execution

The immutable source/input plan is in snapshot.json. A real full-cohort preflight
runs three prespecified signal architectures at the first locus before releasing
the full array. The array runs up to 16 locus jobs concurrently, eight independent
instances per job. The monitor retries only unfinished tasks and stops for
persistent nonconvergence or unresolved failures.

Inspect monitor.status.json and submission.*.json for current scheduler IDs.
Final outputs, created only after complete acceptance:

- comparison.tsv: complete paired grid (558,600 method/policy rows).
- summary.tsv and paired_differences.tsv: overall estimates and paired intervals.
- summary_stratified.tsv and paired_differences_stratified.tsv: architecture/causal-mode results.
- RESULTS.md and accepted.json: report and completion gate.

Four tests passed, including actual execution of both backends, the shared purity
rule, unchanged PIPs after CS filtering, paired-bootstrap identity, and checkpoint
reuse. The cluster preflight remains a separate production-data gate.

## Completed results (November 2, 2026 verification)

See [the complete interpretation and paired intervals](MATCHED_RESULTS.md).

All 931 tasks passed, yielding 139,650 paired instances / 558,600 rows.
Collector 200217 exited successfully on October 30; no jobs remain active.
Frozen inputs, final output hashes, every task acceptance/result hash, convergence,
and unique paired keys were rechecked successfully. See the run root's
`final-verification.json`, `accepted.json`, and `RESULTS.md`.

With common purity ≥0.5 in the full-cohort held-out evaluation, MultiSuSiE versus
SuSiEx gives causal recall 69.85% versus 68.58%, CS coverage 92.88% versus 95.09%,
FDR at PIP >0.95 of 4.23% versus 3.19%, and mean CS size 7.84 versus 7.38.
All four paired differences have marginal 95% locus-bootstrap intervals excluding
zero. Any-causal capture and high-confidence yield differences have intervals
spanning zero. This supports a recall–calibration tradeoff; SuSiEx is better for
calibrated, compact sets in this configuration. Unfiltered results and the selected
50k cohort give different recall rankings, so no universal winner is claimed.
The former pilot-only evidence gap is resolved.

For a new PBS run, set `account`, `python`, `workers`, and `concurrency` in
`scheduler.json` before starting the monitor. Use your own scheduler allocation.
