# SuSiEx versus MultiSuSiE: results and figures

The completed matched-filter study has a dedicated report in the primary
`appfl-bio-suite` checkout:

- `local/output/susiex-vs-multisusie/index.html`: offline report and figure gallery.
- `local/output/susiex-vs-multisusie/RESULTS.md`: detailed findings and paired intervals.
- `local/output/susiex-vs-multisusie/graphs/`: ten figures, each as PNG, PDF and SVG.
- `local/output/susiex-vs-multisusie/tables/`: complete estimates and paired differences,
  both overall and stratified by architecture and causal mode.
- `local/output/susiex-vs-multisusie/report-manifest.json`: input/output checksums,
  generator hash, figure titles and captions.

These are local scientific artifacts, excluded from Git. The numerical source is
`local/output/matched-multisusie-susiex-20261029/`, containing all 139,650 paired
instances / 558,600 method-policy rows. No fits were rerun or pooled across cohorts
for these figures. See [the design](MATCHED_COMPARISON.md),
[findings](MATCHED_RESULTS.md), and [archive guide](MULTISUSIE_ARTIFACTS.md).

## Figure inventory

| Filename pattern | Contents |
| --- | --- |
| `benchmark_full_cohort` / `evaluation_full_cohort` | Both methods' six metric estimates, with 95% intervals, for the full 150k cohort under both policies |
| `{phase}_{policy}_cohorts` | MultiSuSiE minus SuSiEx for all seven cohort configurations, six metrics, with paired 95% intervals |
| `{phase}_{policy}_architectures` | Full-cohort recall, coverage and FDR differences by architecture and causal mode, with paired 95% intervals |

`phase` is `benchmark` or `evaluation`; `policy` is `common_purity_0.5` or
`unfiltered`. The benchmark has 79 loci; held-out evaluation has 54. Architecture
labels use `ncsl` for causal count, `h2` for simulated heritability and `rg` for
simulated effect correlation. MultiSuSiE's prior rho remains fixed at 0.75.

Intervals use 2,000 locus-bootstrap replicates. Difference plots use joint resampling
of both methods. They are marginal intervals without multiplicity adjustment.
Positive differences mean a larger metric, not necessarily a better result. Coverage
should be judged relative to nominal 95%; lower FDR and smaller credible sets are
preferable only in conjunction with useful recall and calibration. Mean set size
includes all retained sets, unlike the older poster's smallest causal-containing set
statistic. Cohort draws overlap and are not independent replicates.

Native PIP metrics are identical across the two policies because the external filter
changes credible sets, not PIPs. Missing estimates/intervals are omitted rather than
shown as zero. Native model/precision differences remain, and archived runtimes do
not support comparative speed plots.

## Reproduce the report

From a checkout containing the PR's plotting script, with coordinator plotting
packages installed:

```bash
python scripts/fine-mapping/plot_matched_comparison.py \
  --root /path/to/appfl-bio-suite/local/output/matched-multisusie-susiex-20261029 \
  --out /path/to/new-comparison-report
```

Use a new destination. The script verifies all five numerical outputs against the
original acceptance record before rendering; relocated input tables are resolved by
basename under `--root`, preserving historical checksums. It copies the four summary
TSVs into the report and leaves the full per-instance table in the source run.

The complete source is also exported under
`local/archive/fine-mapping-multisusie-20261102/` in the primary checkout. Before the
PR is merged, the renderer can be run from that archive's `scripts/fine-mapping/`
with its `local/venv-multisusie/bin/python`. This does not require the deleted worktree.

The original MultiSuSiE PPTX/PDF poster is preserved separately at
`local/scratch/poster-multisusie/`; the new report specifically compares the two methods.
