# Five-site class-mix sweep for the St. Jude discussion

Status: complete. All 21 seed/convergence points and 168 model records are
committed. The three-seed specification was locked before production training.
This is a synthetic mechanism demonstration using ALMA3 CV:e11, not clinical
validation, real hospital data, or an institutional deployment.

## Question

As five sites move from sharply different class proportions to the same class
proportions, do separately trained local models approach a shared federated
model? Show every local model, their mean, and their best performance. Do not
claim federation beats every local unless the measured best-local comparison
supports it. The previous five-site experiment had one nearly perfect local.

## Controlled design

| Site | Extreme counts: classes 5 / 8 / 58 / 76 | Identical endpoint |
|---|---|---|
| S0 | 97 / 1 / 1 / 1 | 30 / 20 / 20 / 30 |
| S1 | 2 / 96 / 1 / 1 | 30 / 20 / 20 / 30 |
| S2 | 1 / 1 / 96 / 2 | 30 / 20 / 20 / 30 |
| S3 | 1 / 1 / 1 / 97 | 30 / 20 / 20 / 30 |
| S4 | 49 / 1 / 1 / 49 | 30 / 20 / 20 / 30 |

Convergence values: 0, 0.05, 0.10, 0.25, 0.5, 0.75, 1. Linear interpolation is rounded while
preserving BOTH 100 patients per site and pooled class totals 150/100/100/150.
Extra resolution near zero captures changes when rare classes gain their first few examples.
These extra points were specified before development training results.
All four classes remain present at every site. Class-wise patient prefixes from
fixed independently generated reservoirs pair training samples across points.
The same class-balanced validation patients (25/class/site) are used at every
point and for every model within a data seed. The task does not become easier
because the evaluation distribution changes.

Site-specific 70-degree signature rotations, additive shifts and probe coverage
are held fixed. The endpoint has identical **class ratios**, not identical
patient measurements or identical site measurement mechanisms. Existing model
weights are restored before each arm and point; there is no warm-start path
that gives later points extra training.

Compare zero-shot, known-support zero-shot, pooled, FedAvg and local S0–S4.
All trained arms use the same released initialization, heads-only updates,
learning rate 0.3, clipping norm 1, ten rounds and two steps per round. The sweep
is nonprivate; FedAvg uses the existing secure-aggregation protocol simulation.
Patient-level gradients are still clipped individually when computed in batches.
Pooled/federated training uses 500 patients versus 100 per local model; this
measures combining data, not an equal-data-budget algorithm advantage.

The metric is mean balanced accuracy over the five evaluation sites. The best
local is chosen separately within each seed on these same evaluation cohorts:
an intentionally strong oracle comparator, not a deployable selection policy.
Three seeds (42–44) supply descriptive paired bootstrap intervals. Seed 41 is
reserved for development checks. No test cohorts or calibrated decision policies
are used in this sweep; synthetic validation results must be labeled accordingly.

## Reproducibility and recovery

Specification: `src/appfl_bio_suite/experiments/methylation/configs/alma3-skew-sweep.json`.
Run/resume instructions: [ALMA3-HANDOFF.md](ALMA3-HANDOFF.md).
Each round has an atomic progress file and checkpoint hash. Each completed point
publishes compact scores and makes a git commit when `--commit-progress` is used.
Weights and patient representations remain outside git. The source release is
immutable, and release licenses and attribution accompany all derivative models.

Presentation artifacts: `sweep.png`, `sweep.svg`, `sweep.pdf`,
`all-locals.png`, `site-mixes.png`, `summary.csv`, and `paired-gaps.csv`.


## Measured results

Balanced accuracy averaged across three data seeds on the same balanced,
shifted validation cohorts at all five sites:

| Class-mix convergence | FedAvg | Mean local | Best local |
|---|---:|---:|---:|
| 0 (extreme skew) | 100.00% | 72.72% | 80.53% |
| 0.05 | 100.00% | 91.08% | 99.40% |
| 0.10 | 100.00% | 97.45% | 100.00% |
| 0.25 | 100.00% | 99.79% | 100.00% |
| 0.50 | 100.00% | 99.96% | 100.00% |
| 0.75 | 100.00% | 99.99% | 100.00% |
| 1 (identical ratios) | 100.00% | 100.00% | 100.00% |

Federation exceeds every measured local model at the extreme point in all three
seeds. Its average advantage over the per-seed best local is 19.47 percentage
points. Local catch-up happens early: the best local matches federation at 0.10,
and mean local accuracy reaches 99.79% at 0.25, while sites still have very
different class ratios. Pooled training also scores 100% throughout.

Use the [full sweep PDF](figures/alma3/skew-sweep/sweep.pdf) together with the
[early-convergence PDF](figures/alma3/skew-sweep/sweep-early.pdf). Supporting
figures show [all five locals](figures/alma3/skew-sweep/all-locals.png),
[site class counts](figures/alma3/skew-sweep/site-mixes.png), and
[each local at its own site](figures/alma3/skew-sweep/own-site.pdf).
At extreme skew, federation also beats all five locals on their own balanced
site cohorts in each seed. This rules out cross-site evaluation as the sole
explanation for that endpoint gap.

[Aggregate scores and audit](results/alma3-skew-sweep/REPORT.md) include paired
seed intervals, checkpoint hashes, and a supplementary
[four-class unshifted retention view](figures/alma3/skew-sweep/retention.png).
Retention in the other 105 classes was not evaluated.

## Interpretation limits

This demonstrates a mechanism on synthetic source-template data with deliberately
constructed measurement shifts. It is not evidence of clinical accuracy or an
expected gain at St. Jude. Local training uses ordinary unweighted samples;
class-balanced local resampling and reweighting were not tested. Therefore this
compares the implemented local baseline, not every possible local strategy.
Federation combines 500 patients while each local has 100. Three data seeds
provide descriptive uncertainty only. The endpoint equalizes class ratios while
preserving site measurement differences. Real test cohorts remain untouched.

Integrity audit passed: 168 model records, 1,680 site/evaluation records, and
150 unique final checkpoint hashes verified; site and pooled counts preserved;
25 validation patients per class per site; zero-shot scores identical across
levels; no test access. Sixteen implementation tests passed. A real-worker
serial/batched gradient comparison differed by at most 7.45e-9; resuming the
development run preserved all 141 checkpoint hashes and modification times.
