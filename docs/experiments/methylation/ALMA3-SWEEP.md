# Five-site class-mix sweep for the St. Jude discussion

Status: development endpoint check running; the three-seed specification is
committed before that check finishes. Results will be added when complete.
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

Expected presentation artifacts: `sweep.png`, `sweep.svg`, `sweep.pdf`,
`all-locals.png`, `site-mixes.png`, `summary.csv`, and `paired-gaps.csv`.
