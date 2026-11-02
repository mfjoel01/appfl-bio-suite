# Federated methylation — design and record of runs

> **Work in progress.** This is a single-node research simulation. It does not run on
> partner endpoints: `partner-bundle` and the remote drivers refuse it. The results below
> are recorded as measured. None of them yet supports near-pooled federated performance
> on real data.

## What this experiment is

A **DNA-methylation subtype classifier trained across simulated hospitals**, built toward
leukemia classification from sparse, nanopore-like methylation calls. A hospital that
sees few patients of a subtype cannot train a good classifier alone. A shared model could
help it without anyone pooling patient data. The open questions are how close federated
training comes to pooling everything, and whether the smallest site gains.

Every arm trains in one process, on one machine, from identical initial weights:

| | Arm | What it trains on |
| --- | --- | --- |
| E1 | Centralized | All training patients pooled — the upper reference |
| E2 | Local | Each of four sites alone — the per-site lower reference |
| E3 | Federated | The four sites under APPFL sample-weighted FedAvg |
| E4 | Coverage | E1–E3 evaluated at several observed-CpG counts |

The serial driver calls APPFL's `FedAvgAggregator` directly. It isolates the learning
question from transport. Remote agents, Globus Compute dispatch and privacy mechanisms are
not implemented here.

## Design

| | |
| --- | --- |
| Input | Beta values, samples × CpGs, NaN = unmeasured ([DATA.md](DATA.md)) |
| Encoding | Per-probe calls: +1 methylated, −1 unmethylated, 0 unobserved |
| Model | Residual MLP: LayerNorm, SiLU, dropout; no batch statistics |
| Optimizer | AdamW, reset every communication round for every arm |
| Loss | Class-weighted cross-entropy |
| Federation | Four simulated sites, synchronous FedAvg, weighted by sample count |
| Calibration | One temperature per model, fitted on validation patients only |
| Metrics | Balanced accuracy, macro-F1, per-class recall, ECE, confident-call rate and accuracy at ≥ 0.90 |

### Why patients, not samples, are the unit of every split

Repeated measurements of one patient stay together. About 20% of each subtype's patients
are held out for test, and 10% of the remainder for validation, before anything is
fitted. The variable-CpG panel, the class mapping and the class counts are all derived
from training patients only. Metrics are computed per sample on the shared test set.

### Simulated sites

Training patients are assigned to four sites with fixed capacities: the smallest site
receives `small_fraction` of them and the rest are split evenly. Each class draws site
preferences from a Dirichlet(`alpha`) distribution, which produces label skew. Setting
`min_site_class_patients` reserves that many patients of every class at every site
first, so a local model's deficit cannot be explained only by an absent class. Synthetic
cohorts may add a small per-site batch shift to training values. Supplied cohorts never
receive one.

### Sparse input

Each covered probe yields one biased binary call:
P(+1 | β) = σ(−0.155 + 1.517 · logit β). Training draws coverage log-uniformly between
`min_coverage` and the full panel. Evaluation fixes the observed-CpG count per level, and
every model sees identical calls. Missing betas are never observed.

### Equal data passes, not equal optimizer steps

Per round, the pooled model, each local model and each federated client make one pass
over their own data. A pooled epoch therefore takes more optimizer steps than a client
epoch. Pooled and federated training use class weights from pooled training counts.
Local models use their own site's counts. The comparison fixes communication rounds, not
compute.

### Relationship to Lamprey

The caller-bias coefficients and the residual-MLP design follow the published methods of
Lamprey (Achterberg et al., medRxiv 2026). This is an independent, smaller
implementation. No Lamprey code or weights are used. It omits Lamprey's fitted read-depth
distribution, multiple reads per CpG, ensemble, per-coverage calibration and broader
taxonomy. It must not be presented as Lamprey or as a reproduction of its results.

## Configurations

Packaged scenarios, selected with `--config NAME`. Any `RunConfig` field can be set in a
YAML file passed by path instead.

| Config | Data | Purpose |
| --- | --- | --- |
| `ci-tiny` | 240 synthetic patients, 512 CpGs, 3 classes | Install check; trains in seconds |
| `poc` | 768 synthetic patients, 50,000 CpGs, 6 classes | Full-size execution check |
| `stress` | 768 patients, 2,048 CpGs, 2% informative | Hard synthetic case for the paired study |
| `tcga-laml` | The GDC cohort in [DATA.md](DATA.md), supplied explicitly | First real-data run |
| `uneven` | 2,400 patients, unequal classes, label noise, missing values | Illustrative uneven federation |
| `uneven-curves` | `uneven`, trained 30 rounds, scored at seven coverage levels | Validation-only learning and sparsity curves for `uneven` |

## Record of runs

All runs below were made on 2026-10-01 with APPFL 1.10.0 on CPU. Each run directory
holds its configuration, preprocessing record, checkpoints, per-round history, metrics,
figures and a manifest of versions, commit and checksums. Run artifacts and data stay
outside version control. The synthetic runs are regenerable from their configuration
and seed.

### Full-size synthetic execution check (`poc`, seed 42)

Splits of 552 training, 60 validation and 156 test patients; site sizes 166, 166, 165
and 55. Pooled and federated models scored a balanced accuracy of 1.0 at every coverage
level. Local scores were 1.0, 0.83, 0.83 and 0.67. Each local shortfall corresponded
exactly to the classes absent from that site's training data.

**What it shows:** the pipeline runs end to end at full panel size. The generator is
too easy to compare arms, because every probe is informative.

### Paired stress study (`stress`, seeds 42–44)

Only 2% of probes carry signal. Each seed runs two conditions that share data, holdouts,
panel and site capacities. In one, sites are unconstrained. In the other, every site
holds at least two patients of each class. Batch shifts are off, so the pooled arm is
identical across conditions, and this was checked. Balanced accuracy at 2,048 observed
CpGs with every class present, as mean ± population SD across seeds:

| Pooled | Federated | Smallest site alone |
| ---: | ---: | ---: |
| 0.479 ± 0.073 | 0.297 ± 0.070 | 0.192 ± 0.032 |

**What it shows:** federation helps the smallest site, by 10.5 points, but stays 18.2
points below pooled. At 32 observed CpGs, federated scores are near the six-class chance
level of 0.167.

### TCGA-LAML, adult AML FAB morphology (`tcga-laml`, seed 42)

The first real-data run, on the 148-patient open GDC cohort described in
[DATA.md](DATA.md). Splits of 109 training, 11 validation and 28 test patients; site
sizes 33, 33, 32 and 11, each holding every class. The analysis was fixed before
training. Balanced accuracy at 50,000 requested CpGs (48,932 observed on average):

| Pooled | Federated | Local S1 | Local S2 | Local S3 | Local S4 (11 patients) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 56.25% | 43.75% | 50.00% | 25.00% | 43.75% | 43.75% |

The federated model predicts only M2 and M5: its recall for M1/M2/M4/M5 is
0/100/0/75%. The pooled model also misses M1. Across 1k–50k CpGs, pooled scores
43.75–56.25% and federated 37.50–43.75%. With 28 test patients, one patient moves a
class recall by 12.5–25 points.

**What it shows:** no evidence yet of near-pooled performance, or of a gain for the
smallest site, on real data. **This test split has now been inspected.** Further tuning
must use validation only, and any reuse of this split must be disclosed.

### Uneven synthetic illustration (`uneven`, seed 2026)

A deliberately engineered scenario meant to illustrate one pattern: strong centralized
performance, federation somewhat lower, and a small skewed site far lower. The setup:

- 2,400 patients in six classes in proportions 30/23/17/13/10/7.
- 5% informative probes, 4% profile/label discordance and 3% missing values.
- Four sites of 542, 541, 541 and 104 training patients. The smallest holds 83 patients
  of one class and 3–7 of each other.

Four candidate configurations were scored on validation patients only, with
`--validation-only`, which never predicts test patients. A regression test checks that
test values cannot change its results. The development targets were 94–97% centralized,
90–94% federated and below 75% for the smallest site.
The frozen configuration was then evaluated once on 480 test patients:

| | Centralized | Federated | S1 | S2 | S3 | S4 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Balanced accuracy | 97.7% | 90.1% | 92.2% | 48.7% | 64.0% | 19.1% |
| Macro-F1 | 0.966 | 0.906 | 0.911 | 0.360 | 0.582 | 0.075 |

**What it shows:** within this one chosen scenario, federation recovers most of the
pooled performance for sites that hold too little of most classes. Each local model
scores near zero on exactly the classes its site barely holds. Federation does not beat
the best local site, S1, and their 95% intervals overlap: 86.9–93.0% and 89.8–94.4%.
Treat this as an illustration of a mechanism, not a performance estimate.

**Learning and sparsity curves** (`uneven-curves`, validation patients only). The same
models were trained for 30 rounds and scored after every round. They reproduce the
selection runs' validation scores exactly at rounds 10, 15 and 30, and their first 15
rounds are identical to the tested run.

- **Training length.** From round 21 on, federated stays within 2 points of centralized,
  and both reach 96.8% at round 30. The tested gap reflects the 15-round budget, not a
  converged difference.
- **The smallest site.** Training alone, it swings between 17% and 68% from one round to
  the next. Its tested 19.1% is one draw from an unstable model.
- **Sparse input.** Accuracy falls steeply with fewer observed CpGs. After 30 rounds,
  centralized and federated both score about 80% at 1,024 CpGs, 50% at 256 and 31% at 64.
  The scenario trained on at least 819 CpGs per patient, so it was never built for very
  sparse input.

## What is not done

- **Partner execution.** No shipped trainer, transport, partner bundle or remote driver.
  A partner protocol must first fix the data, the CpG panel, local preprocessing and the
  outputs returned.
- **Real nanopore input.** Sparse calls are simulated from array betas.
- **A pediatric cohort.** The planned pediatric ALL cohort, GSE49031, has been withdrawn
  ([DATA.md](DATA.md)). The real-data run is adult AML morphology, not pediatric subtype.
- **Untouched evaluation.** The TCGA-LAML and uneven test splits are both inspected. A
  confirmatory claim needs a new, prespecified evaluation.
- **GPU execution and privacy mechanisms** are unvalidated and unimplemented,
  respectively.

## References

- Achterberg et al. *Lamprey*. medRxiv, 2026. doi:10.64898/2026.07.02.26356825.
  Methods only; article licence CC BY-NC-ND. Code and weight terms are separate.
- [APPFL](https://github.com/APPFL/APPFL) — `FedAvgAggregator`, release 1.10.0.
- The Cancer Genome Atlas, TCGA-LAML, via the
  [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/projects/TCGA-LAML).
- Nordlund et al. *Genome Biology* 14:r105, 2013 (GSE49031; data since withdrawn).
