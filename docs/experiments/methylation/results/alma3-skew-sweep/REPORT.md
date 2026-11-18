# Five-site synthetic class-mix sweep

Analysis performed using ALMA3.

Validation-only mechanism study. This is not clinical evidence or a DP release.

| Convergence | FedAvg | Mean local | Best local | FedAvg minus best (pp) |
|---|---:|---:|---:|---:|
| 0 | 100.00% | 72.72% | 80.53% | 19.47 |
| 0.05 | 100.00% | 91.08% | 99.40% | 0.60 |
| 0.1 | 100.00% | 97.45% | 100.00% | 0.00 |
| 0.25 | 100.00% | 99.79% | 100.00% | 0.00 |
| 0.5 | 100.00% | 99.96% | 100.00% | 0.00 |
| 0.75 | 100.00% | 99.99% | 100.00% | 0.00 |
| 1 | 100.00% | 100.00% | 100.00% | 0.00 |

Best local is the maximum per seed over all five local models, evaluated on the same balanced cohorts at all sites. It is an oracle comparator, not a selected deployable model.

Interpolation preserves all site counts and pooled class totals. Patient reservoirs and evaluation cohorts are paired across sweep points. Site measurement shifts remain fixed; only class ratios converge.

Federated/pooled training uses 500 patients; each local model uses 100. Learning rate, clipping, number of rounds and local steps are fixed. No early stopping or per-point retuning. Zero-shot known-class support is an explicit no-training control.

Intervals resample data seeds jointly across convergence points and compare models within seed. With few seeds these intervals are descriptive, not a clinical-population uncertainty estimate. No calibrated decision policy is fitted in this sweep.

See sweep.png/svg/pdf, all-locals.png, site-mixes.png, summary.csv and paired-gaps.csv.

Local training uses ordinary unweighted examples; class-balanced local resampling or reweighting is not evaluated. This does not establish superiority over all possible local training strategies.

Own-site comparisons and four-class unshifted retention are supplementary views of the same recorded predictions, not additional training or test data.
