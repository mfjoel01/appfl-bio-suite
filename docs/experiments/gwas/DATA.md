# Federated GWAS — how the data comes into existence

This experiment generates its own data. The coordinator runs the simulation and
distributes per-site bundles; partners receive data and never simulate.

**Read the reproducibility section before using any output in a publication.** It is
explicit about what this pipeline does and does not reproduce.

---

## Reproducibility: what is and is not regenerable

### The original cohort cannot be regenerated

The GWAS results previously reported from this work were computed from a specific input
cohort: a PLINK1 triple of **100,000 synthetic European-ancestry samples over ~240,000
variants**, roughly 6 GB.

That cohort was generated in a Google Colab notebook. **The notebook was never
version-controlled and is not recoverable.** Nothing in the original code or documentation
recorded where the file came from; it appeared only as a path constant, and the pipeline
began at the step *after* it.

This is stated plainly rather than worked around, because it is a reproducibility claim in
a paper-facing repository:

> **The original input cohort is not regenerable. Results computed from it can be
> recomputed by anyone holding a copy of the file, and cannot be regenerated from this
> repository alone.**

### What this repository does provide

The generation step is now a first-class, versioned part of the pipeline — **step 0**,
ahead of phenotype simulation — with two providers:

| Provider | What it does | Use it when |
| --- | --- | --- |
| `synthetic_ld` (default) | Generates a cohort from a documented model. Deterministic from a seed, no download. | You do not have the original cohort — which is most people. |
| `plink_prefix` | Uses an existing PLINK1 triple. | You hold the original, or your own real cohort. |

So the whole chain — from nothing to per-site data — is runnable from this repository
alone. **But the substitute cohort is a different cohort from a different generative
process. It reproduces the pipeline, not the published numbers.** Do not describe output
from it as reproducing any previously reported figure.

Every run's manifest records which provider was used, so no output is ever ambiguous about
its own provenance.

### If you hold the original cohort

Point a scenario at it:

```yaml
cohort:
  provider: plink_prefix
  plink_prefix: /path/to/EUR.synthetic.100k.ld.maf    # no file extension
```

Everything downstream — phenotype simulation, the site split, every seed — is byte-identical
to the pre-migration pipeline, verified by `tests/test_gwas_simulation_parity.py`.

---

## The pipeline

Three stages, run by one command.

```
step 0   cohort.py       genotypes                 -> PLINK1 triple
step 1   phenotypes.py   + PGS Catalog weights     -> phenotypes, covariates, scores
step 2   bundler.py      split by scenario         -> per-site datasets
```

```bash
appfl-bio-suite simulate gwas --scenario three-site-skewed --out local/data/gwas/runs/three-site-skewed
```

(`--out` is a per-run directory, so repeated simulations never mix with staged inputs or
with each other.)

### Step 0 — the cohort

With the default provider, genotypes are generated from a documented model:

1. **Allele frequencies** are drawn from a Beta distribution and clipped. Beta(0.8, 2.4)
   is skewed toward low frequencies, which is the qualitative shape of a real
   site-frequency spectrum — most variants are rare.
2. **Linkage disequilibrium** is induced by blocks: within a block, each variant's latent
   value mixes a shared block factor with its own noise, weighted so that variants in a
   block correlate at approximately `ld_correlation`. Independent variants would make
   every association test cleanly independent and the simulation unrealistically easy.
3. **Genotypes** come from thresholding two independent haplotype draws against each
   variant's frequency, giving an allele count of 0, 1 or 2.
4. **Missingness** is introduced at a small rate, so the mean-imputation path downstream
   is genuinely exercised rather than dead code.

This is a coarse model of population structure. It is adequate for exercising and
demonstrating the pipeline. It is **not** a substitute for real or well-validated
synthetic genotypes if the goal is a claim about genetics rather than about federation.

### Step 1 — phenotypes, covariates, polygenic scores

Two traits, chosen because they exercise different analyses:

- **BMI** — continuous, tested by covariate-adjusted linear regression
- **T2D** — binary, tested by a logistic score test

Real PGS Catalog scoring files are used when present beside the cohort. When they are
absent, matched substitutes are generated against the cohort's own variants so the
pipeline remains runnable end to end — with a warning, and with their checksums in the
manifest, so a run using substitutes is never mistaken for one using published weights.

Two independent phenotype replicates are drawn, from different seeds:

- `phenotypes_gwas.csv` — the replicate the GWAS is run on
- `phenotypes_pgs_eval.csv` — an independent replicate the polygenic score is evaluated on

**The separation is load-bearing.** Evaluating a score against the phenotypes its weights
were estimated from would inflate every reported R² and AUROC. Covariates are drawn once
and shared by both.

Key parameters, all defaulting to the published design:

| | Default | |
| --- | --- | --- |
| `h2_t2d` / `h2_bmi` | 0.20 / 0.25 | Heritability; sets environmental noise variance |
| `t2d_prevalence` | 0.119 | Liability threshold for case assignment |
| `effect_scale_*` | 0.50 | Scale on the genetic effect |
| `seed_covariates` | 42 | |
| `seed_gwas_replicate` | 42 | |
| `seed_pgs_replicate` | 43 | Distinct — this is what keeps the replicates independent |

### Step 2 — the site split

The cohort is shuffled once and cut into consecutive, non-overlapping slices. Each site
receives a complete self-contained dataset: its own genotype subset plus the matching
phenotype, covariate and score rows.

All sites share the full variant set — meta-analysis requires it — and no samples. The
pipeline verifies disjointness and fails if it is violated: overlapping sites would make
the federation train twice on the same people while still reporting the full sample size.

Site sizes come from the scenario, not from code. The published design uses a deliberately
unequal split (18032 / 9237 / 25028 / 40752 / 6951) because cross-site heterogeneity is
the thing federated analysis has to cope with, and a scenario where every site is the same
size tests the easy case.

---

## Scenarios

```bash
appfl-bio-suite simulate gwas --list-scenarios
```

| | |
| --- | --- |
| `ci-tiny` | 2 sites, 400 samples. Seconds. Proves the chain works; no meaningful genetics. |
| `two-site-balanced` | 2 equal sites, 20k samples. A clean baseline. |
| `three-site-skewed` | 3 unequal sites, matching the first three of the published split. |
| `published-five-site` | The published design. Reproduces the pipeline and parameters, **not** the original genotypes. |

### Defining your own

Copy any scenario in `src/appfl_bio_suite/experiments/gwas/configs/simulation/`, edit, and
pass the path:

```bash
appfl-bio-suite simulate gwas --scenario /path/to/my-scenario.yaml
```

The knobs worth knowing: `sites` (the split), `cohort.n_samples` and `cohort.n_variants`
(scale), and `phenotypes.data_sim_scaling` (a variant fraction, for fast smoke runs).

---

## Provenance

Every run writes `run_manifest.json` recording the scenario, every RNG seed, the suite git
commit, input and output checksums, package versions, platform, and timestamps.

The pre-migration pipeline emitted a summary of diagnostics but no record of *how* the
outputs were produced. For work supporting a publication, that means the numbers could not
be regenerated with confidence — not because the code was wrong, but because nobody could
prove which code and which inputs produced them.

Confirm a rerun reproduced a previous one:

```bash
appfl-bio-suite simulate gwas --verify local/data/gwas/runs/three-site-skewed/run_manifest.json
```

That re-checksums every output and reports any difference — a mechanical answer rather
than eyeballing summary statistics.

**If verification fails after an environment change,** compare the manifest's `packages`
block against the current environment before assuming the code is wrong. A moved `numpy`,
`pandas` or `pandas-plink` is the usual cause.

---

## Numerical stability

The simulation is **frozen**: `tests/test_gwas_simulation_parity.py` asserts that all 35
outputs are byte-identical to those of the pre-migration pipeline for a fixed seed,
against checksums captured by running the original scripts unmodified.

Four ways this breaks silently, all guarded by that test:

1. **RNG algorithm.** The site split uses numpy's *legacy* global generator
   (`np.random.seed` + `np.random.shuffle`). Modernizing to `default_rng` is a different
   algorithm and produces a completely different assignment of individuals to sites, with
   no error.
2. **Population vs. sample standard deviation.** `ndarray.std()` is ddof=0;
   `pandas.Series.std()` is ddof=1. Substituting one rescales every score.
3. **Chunked accumulation.** Scoring is chunked at 2000 variants. Floating-point addition
   is not associative, so the chunk size affects low-order bits.
4. **Operation order in the PGS overlap.** The groupby-join-filter sequence determines
   which duplicate variants survive.

Regenerating the golden checksums is a decision to change published numerical behaviour,
not a routine fix.

---

## What partners receive, and what they are told

Each site receives only its own bundle. The six files the loader requires are listed in
their setup guide.

**A partner's own setup guide never mentions that the data is simulated.** That is
deliberate, and it is why this document lives outside `docs/partner/` and is never
bundled: a partner should not learn how their data came into existence by reading their
own instructions. If they ask, tell them — but it is not their setup step.

## Dependencies

Simulation is coordinator-side only. Its dependencies live in the `gwas-sim` extra and
never land on a partner, who installs `appfl-bio-suite[gwas]`.

**No `plink2` binary is required.** The pipeline writes `score_bmi.txt` and
`score_t2d.txt` in plink2's `--score` input format for optional external use, but nothing
here invokes plink2 — scoring is done in numpy. (Earlier documentation implied otherwise;
it was describing an intention rather than the code.)
