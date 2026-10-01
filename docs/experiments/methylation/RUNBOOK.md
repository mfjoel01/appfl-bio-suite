# Federated methylation — runbook

Everything runs on the coordinator's machine, in one process. There is no partner step,
no endpoint and no federation file.

## Before you start

Python 3.12, from a clone of this repository:

```bash
pip install -e ".[methylation]" -c constraints.txt
```

A CPU is enough. Each run writes into a new directory and refuses to overwrite one, so
give every command a fresh `--out-dir`.

## Prove the install

```bash
appfl-bio-suite run methylation --driver serial --config ci-tiny \
  --out-dir local/output/methylation-ci
```

This generates 240 synthetic patients with 512 CpGs, assigns four sites and trains six
models for three rounds. It then evaluates five coverage levels. It takes seconds and
needs no downloaded data.

## What a healthy run looks like

One `Round k/N complete` line per round, then a summary such as
`{'output': 'local/output/methylation-ci', 'models': 6, 'evaluations': 30}`.
`status.json` ends as `{"complete": true, "stage": "complete"}`, and `REPORT.md` states
the federated-minus-centralized gap and the federated-minus-local gain for each site.

## Choosing a scenario

`--config` takes a packaged scenario name — `ci-tiny`, `poc`, `stress`, `tcga-laml`,
`uneven` or `uneven-curves`; see [ABOUT.md](ABOUT.md#configurations) — or a path to a YAML file of
`RunConfig` fields. Copy a packaged file from
`src/appfl_bio_suite/experiments/methylation/configs/` to start one. An
unknown or misspelled field is rejected by name. `--dry-run` writes the resolved
`config.json` without training.

`--driver` must be `serial`. Remote drivers, `--watch` and `--federation` are refused,
because this experiment has no partner execution yet.

For GPU execution, set `device: cuda`. Runs enforce deterministic algorithms, so also
export `CUBLAS_WORKSPACE_CONFIG=:4096:8`. GPU execution has not been validated. Use your
allocation's own scheduler for anything larger than `poc`.

## Real data: TCGA-LAML

Build the cohort once — see [DATA.md](DATA.md#tcga-laml-from-the-gdc) for what it
selects and the data-use terms — then pass it explicitly:

```bash
python -m appfl_bio_suite.experiments.methylation.gdc --out local/tcga-laml-fab
appfl-bio-suite run methylation --driver serial --config tcga-laml \
  --data-root local/tcga-laml-fab/cohort.npz --out-dir local/output/methylation-tcga-laml
```

The test split from this configuration has already been inspected once. Do any further
tuning with `--validation-only`, and disclose any reuse of the split.

## Developing without touching the test set

Select settings on validation patients only, then evaluate the test set once:

```bash
python -m appfl_bio_suite.experiments.methylation --run --validation-only \
  --config candidate.yaml --out local/output/methylation-candidate-01
```

This mode never predicts test patients and leaves calibration temperatures at 1. Its
metrics report `n_validation` and `"evaluation_split": "validation"`, and it draws no
figures. Add `--input cohort.npz` for supplied data. A config can instead declare
`validation_only: true`, which `appfl-bio-suite run` then honours too. Add
`track_validation: true` to score every model on validation patients after each round,
at every coverage level, into `validation_history.json`. Scoring never changes training.
The protocol behind the `uneven` scenario:

1. Declare the candidates and the selection metric before running any of them.
2. Run every candidate with `--validation-only` and keep all of them, including the losers.
3. Freeze the chosen YAML, then run it once without `--validation-only`.
4. Report the test result as measured. Do not rerun seeds to pick a test result.

## Paired stress study

```bash
python -m appfl_bio_suite.experiments.methylation.stress \
  --out local/output/methylation-stress
```

For each seed (default 42, 43 and 44; shorten with `--seeds 42`), this runs `stress`
twice: once with unconstrained sites, once with two patients of every class reserved
at each site. Both conditions share data, holdouts, panel and site capacities. The study
fails if the pooled arm differs between them. The parent directory gets `summary.json`
and a `REPORT.md` of mean and population standard deviation across seeds.

## Figures

Every completed run draws `site_comparison.png`, `coverage_curve.png` and the federated
`confusion_matrix.png`. For a synthetic scenario, the analysis figures go further. For
`uneven`, the tested run plus its validation-only curves:

```bash
appfl-bio-suite run methylation --driver serial --config uneven \
  --out-dir local/output/methylation-uneven
appfl-bio-suite run methylation --driver serial --config uneven-curves \
  --out-dir local/output/methylation-uneven-curves
python -m appfl_bio_suite.experiments.methylation.figures \
  --run local/output/methylation-uneven --curves local/output/methylation-uneven-curves \
  --out local/output/methylation-uneven-figures
```

| Figure | What it answers | Scored on |
| --- | --- | --- |
| `01_who_gains` | How much each site gains, or loses, by joining | Test |
| `02_training_rounds` | How every model improves round by round | Validation |
| `03_sparse_input` | How accuracy holds up as fewer CpGs are observed | Validation |
| `04_class_recall` | Which classes each model gets right | Test |
| `05_site_composition` | How skewed each site's classes are | Training data |
| `A1_data_overview`, `A2_cpg_methylation` | Appendix: what the synthetic data looks like | Training data |

`uneven-curves` repeats `uneven`'s training for 30 rounds and scores it at seven coverage
levels. It never predicts test patients. The module refuses a curves run that does not
share the tested run's cohort, splits and sites, or that did not follow its training
round for round, so the curves describe the tested models. Without `--curves`, figures
02 and 03 are left out.

Every input's checksums are checked first. Each figure is written as PNG and SVG, with
`figure_data.json` recording every plotted value and file hash. Intervals are 95%
bootstrap intervals computed from the confusion matrices. With the same fonts and package
versions, rebuilding produces identical files. Plots use the suite's shared Argonne
style; see `core/plot_style.py` for fonts.

## Reading the output

| File | Meaning |
| --- | --- |
| `status.json` | Last completed stage and round; `complete: true` only at the very end |
| `config.json` | Every resolved setting |
| `cohort.npz` | The input cohort, before any synthetic site shift |
| `preprocessing.json` | Split and site membership, class mapping, probe panel, shifts |
| `checkpoint.pt` | All six model state dictionaries, replaced atomically each round |
| `history.json` | Per-round training loss for every arm and client |
| `validation_history.json` | With `track_validation`: per-round validation scores and confusion matrices |
| `metrics.json` | E1–E4 scores per model and coverage level, with temperatures |
| `site_class_counts.csv` | Training composition of each simulated site |
| `REPORT.md` | Federated gap to centralized, and gain over each local model |
| `manifest.json` | Versions, commit and dirty state, source and output checksums |

The pooled and federated bars in `site_comparison.png` repeat one shared-test score per
site; they are not site-specific test sets. A `null` confident accuracy means no
prediction reached 0.90 confidence. Calibration picks each model's temperature from a
fixed grid on [0.05, 20], which includes 1, using validation patients across all
coverage levels.

## When it is not healthy

| Message | Meaning and fix |
| --- | --- |
| `output directory already exists` | Runs never overwrite; choose a new `--out-dir` |
| `unknown setting(s) ...` | A misspelled field in the YAML |
| `coverage count exceeds selected panel size` | Lower `coverage_counts` or raise `panel_size` |
| `needs at least three patients for holdouts` | A subtype too small to split; relabel or drop it |
| `training class threshold leaves fewer than two classes` | Lower `min_class` |
| `site capacity is too small ...`, `not enough training patients of class ...` | Lower `min_site_class_patients` |
| `CUDA requested but unavailable` | Set `device: cpu`, or run on a GPU node |
| `training loss became nonfinite` | Lower `learning_rate`; check the input for degenerate rows |

A run whose `status.json` still says `complete: false` stopped early. Its checkpoint
holds the completed rounds but not optimizer state, and nothing resumes it. Rerun from
the saved `config.json` into a fresh directory.

## Tests

```bash
pytest tests/test_methylation.py tests/test_methylation_gdc.py
```

Both files are offline: the GDC tests stub every network call. They cover deterministic
reruns, leakage controls, missingness, NPZ validation, sample-weighted aggregation,
validation-only isolation, the CLI paths and the loader's selection rules.
