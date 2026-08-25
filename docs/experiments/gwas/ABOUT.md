# Federated GWAS — design and record of runs

## What this experiment is

A **federated genome-wide association study** with **polygenic-score evaluation**, run
across institutions on different continents without any of them sharing raw genomes.

It reproduces the simulation design of *"Simulation Benchmarking of Differential Privacy
and Federated Learning for Secure Genetic Analysis and Polygenic Risk Scoring"* (Joel et
al.) on all-chromosome synthetic European-ancestry data.

Two traits, chosen because they exercise different analyses:

- **BMI** — continuous; per-variant association by covariate-adjusted linear regression.
- **T2D** — binary; per-variant association by a logistic score test.

## Single-round summary-statistic federated learning

This is **not** iterative model training, and that distinction drives most of the design:

1. Each site runs a complete local GWAS on its own cohort and computes, per variant, an
   effect size (`BETA`), a standard error (`SE`), and a minor-allele frequency (`MAF`),
   plus a local polygenic-score metric.
2. Each site sends **only** those summary statistics. Never genotypes, never phenotypes.
3. The server combines them by **inverse-variance fixed-effect meta-analysis** and reports
   pooled per-variant statistics, genome-wide hits, and pooled score metrics.

Because only per-variant `BETA`/`SE`/`MAF` cross the wire, the payload is a few megabytes
regardless of cohort size, and raw genomic data never leaves the site it originates from.

There is exactly one exchange per site. `num_global_epochs: 1` is correct and not a
placeholder — asking for more would recompute identical results at full cost.

### Why inverse-variance weighting

Each site's estimate is weighted by its precision:

```
w_i  = 1 / SE_i²
beta = Σ(w_i · beta_i) / Σw_i
SE   = √(1 / Σw_i)
```

Larger and less noisy sites count for more, automatically, with no weight to configure.
Under the assumption that all sites estimate the same underlying effect, this recovers
what a pooled analysis of the combined cohort would have found — which is precisely the
claim the experiment exists to test, and why a central baseline over the pooled data is
worth computing alongside it.

## Design

| | |
| --- | --- |
| Federation pattern | Single round, summary statistics |
| Server | Custom `MetaAnalysisAggregator` (coordinator-side, not shipped) |
| Client | Custom `SiteGWASTrainer` (shipped to workers, self-contained) |
| Covariates | Intercept, standardized age, sex |
| BMI test | Residualize phenotype and genotypes on covariates; two-sided t test |
| T2D test | Logistic score test against a covariate-only null; 1-df chi-square |
| Significance | P < 5×10⁻⁸ |
| Data | Synthetic; simulated and distributed by the coordinator — see [DATA.md](DATA.md) |

### Why a score test for T2D

Fitting a separate logistic regression per variant is not tractable at genome scale — that
is 240,000 model fits per site. The score test evaluates every variant against a single
covariate-only null fit. It is the standard approach and near-equivalent at the effect
sizes that matter.

### Why residualization for BMI

Residualizing both the phenotype and the genotypes on the covariates gives the same
per-variant coefficient as fitting the full model, at a fraction of the cost: the covariate
fit happens once rather than once per variant.

### Polygenic scores are evaluated on an independent replicate

Each site scores its own cohort using its own GWAS effect estimates, then evaluates that
score against a **separate phenotype replicate** drawn from a different seed — BMI by R²,
T2D by AUROC. Evaluating against the phenotypes the weights were estimated from would
inflate both badly.

Site-level metrics are pooled by evaluation sample size.

### Why not the hosted web UI

The hosted APPFLx form accepts a fixed algorithm plus model, loss and metric files, and a
dataloader returning a `torch.utils.data.Dataset`. There is **no field for a custom
trainer or aggregator**. This experiment is defined by both, so the hosted form
structurally cannot launch it. Self-driving also keeps the run scriptable and reproducible,
including timing instrumentation.

---

## Record of runs

### 5-site full-scale run

The reference result from the pre-migration pipeline, over the original 100,000-sample
cohort at full variant scale.

> **This run used the original input cohort, which is not regenerable from this
> repository** (see [DATA.md](DATA.md#reproducibility-what-is-and-is-not-regenerable)).
> Running `published-five-site` with the default cohort provider reproduces the design and
> every parameter, but uses different genotypes and therefore produces different numbers.

| | |
| --- | --- |
| Sites | 5 |
| Total N | 100,000 |
| Variants tested | 239,801 (all chromosomes, GRCh37) |
| Scaling | `data_sim_scaling` 1.0, `variant_scaling` 1.0 |
| Split | 18,032 / 9,237 / 25,028 / 40,752 / 6,951 |

**Per-site and pooled metrics**

| Site | GWAS N | Eval N | Local BMI R² | Local T2D AUROC |
| --- | --- | --- | --- | --- |
| Site1 | 18,032 | 18,032 | 0.2520 | 0.6932 |
| Site2 | 9,237 | 9,237 | 0.2699 | 0.7025 |
| Site3 | 25,028 | 25,028 | 0.2505 | 0.6843 |
| Site4 | 40,752 | 40,752 | 0.2537 | 0.6905 |
| Site5 | 6,951 | 6,951 | 0.2599 | 0.6790 |
| **Meta (N-weighted)** | **100,000** | **100,000** | **0.2545** | **0.6897** |

**Genome-wide significant hits** (P < 5×10⁻⁸)

| Trait | Variants | P < 5×10⁻⁸ | P < 1×10⁻⁵ | Max −log₁₀P |
| --- | --- | --- | --- | --- |
| BMI | 239,801 | 33 | 79 | 55.1 |
| T2D | 239,801 | 14 | 33 | 33.8 |

### Earlier validation runs

- **Single-site over Globus Compute** — validated end to end, 18,032 samples,
  BMI R² ≈ 0.26, T2D AUROC ≈ 0.70, with genome-wide BMI hits produced. Re-confirmed on a
  later date with an independent run.
- **3-site functional run**, all three clients dispatching to one endpoint with disjoint
  cohorts: `NUM_CLIENTS=3`, `TOTAL_GWAS_N=52,297`, weighted BMI R² 0.2627, T2D AUROC
  0.6950, 33 genome-wide hits. This also confirmed an import cleanup was
  behaviour-preserving: Site1's local metrics were identical to the earlier single-site
  run.
- **Cross-institution execution confirmed** — a round-trip probe against a partner's
  multi-user endpoint returned a result executing as the expected service account on a
  real compute node at that institution, proving the full authorization chain.

### Migration note

The simulation pipeline was migrated into this repository with an output-identity
guarantee: `tests/test_gwas_simulation_parity.py` asserts all 35 outputs are byte-identical
to the pre-migration scripts for a fixed seed, against checksums captured by running those
scripts unmodified.

The training path was refactored for self-containment — the trainer no longer imports
sibling modules, which removed the `PYTHONPATH` requirement from every partner's endpoint
configuration. Per-site plotting moved to the coordinator, which has all the information
needed to produce it and which keeps matplotlib off partner clusters. Numerical outputs are
unchanged.

### Template for a run of record

```
### YYYY-MM-DD — <what this run was>

| | |
| --- | --- |
| Sites | N |
| Total N | |
| Variants | |
| Scenario | |
| Sim manifest | <sha256 of run_manifest.json> |
| Suite commit | <sha> |

| Site | GWAS N | Local BMI R² | Local T2D AUROC |

Meta: ...
Hits: BMI n, T2D n

Reproduce:
    appfl-bio-suite simulate gwas --scenario <name> --out <dir>
    appfl-bio-suite run gwas
```

Record the **simulation manifest hash** alongside the commit. The commit says which code
ran; the manifest says which data it ran on, and without both the run is not reproducible.

---

See [DATA.md](DATA.md) for the data pipeline and its reproducibility limits, and
[RUNBOOK.md](RUNBOOK.md) for how to run it.
