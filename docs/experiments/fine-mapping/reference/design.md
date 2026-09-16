> Current protocol: [September scientific corrections and rerun](../SCIENTIFIC_RERUN.md). Historical numerical results and earlier exactness/calibration claims below require the new validation gates.

> **ARCHIVED RECORD (migrated 2026-08-21).** This document describes the standalone
> standalone `fedfm-simulation` tree, which is retired. The body
> below is preserved exactly as written; read its paths through this mapping: the source
> tree → this repository (`appfl-bio-suite`); its `data/` → `local/data/fine-mapping/`;
> its `scripts/` → `scripts/fine-mapping/` (archived originals in [pbs/](pbs/)). The
> runnable procedure is [docs/experiments/fine-mapping/RUNBOOK.md](../RUNBOOK.md) — do not
> run commands quoted here.

# FedFM: A Federated Fine-Mapping Simulation Pipeline for Multi-Ancestry Cohorts

**Working draft — design and first end-to-end realisation**
**Author:** M. Joel (Argonne National Laboratory)
**Date:** 2026-06-17
**Code:** the standalone `fedfm-simulation` tree (retired; migrated into this repository)

---

## Abstract

We describe the design and first end-to-end realisation of a controlled-
truth simulation pipeline supporting federated, multi-ancestry statistical
fine-mapping experiments. The pipeline partitions the HAPNEST synthetic
cohort (1,008,000 individuals, six superpopulations) into three
non-overlapping site cohorts of 50,000 individuals each, designed to
mirror the ancestry compositions of three real federated research partners:
Argonne National Laboratory (US), Covenant University (Nigeria), and
Mohamed bin Zayed University of Artificial Intelligence (UAE). Across these
sites we select up to 100 loci on chromosome 1, stratified by inter-site
linkage-disequilibrium (LD) divergence, generate causal genetic
architectures parameterised by (number of causal SNPs, heritability,
cross-ancestry genetic correlation), and simulate phenotypes whose effect
sizes are indexed by each individual's ancestral superpopulation rather
than by the cohort to which the individual was assigned. The first end-to-
end run on Polaris produced 79 non-overlapping loci, 15 causal architectures
per locus, 10 replicates per architecture, and 35,550 phenotype files
across the three sites. The realised data package matches the design
invariants to within sampling noise (median heritability calibration error
0.4%, empirical cross-ancestry effect-size correlations within 5% of
target). The package is intended as input to downstream multi-ancestry
Bayesian fine-mapping methods (e.g., SuSiEx) under a federated execution
model, and the repository now carries two such consumers side by side: a
centralized SuSiEx baseline that pools each ancestry across sites outright,
and a federated path in which each site emits only raw second-moment
aggregates over its own individuals and a coordinator sums and standardizes
them (§11.3). The two are held to *equality*, not similarity — the theory
says the federated fit reconstructs the centralized one exactly, so the
comparison between them is a correctness test rather than a benchmark.

---

## 1 Background and motivation

Statistical fine-mapping aims to identify the small set of variants causally
responsible for a genome-wide-association-study (GWAS) signal at a locus,
typically by computing posterior inclusion probabilities (PIPs) over candidate
variants conditional on observed marginal summary statistics and an LD
reference panel. Multi-ancestry fine-mapping methods exploit between-population
differences in LD structure to narrow credible sets — variants that are in
strong LD with a causal variant in one ancestry may not be in another, so
the intersection of credible sets across ancestries is often substantially
smaller than any individual set. SuSiEx (Yuan et al., 2024) is the current
reference implementation for this paradigm.

Real multi-ancestry collaborations rarely have the luxury of pooling
individual-level genotypes across institutions. Data-residency,
patient-privacy, and institutional-review constraints typically require that
genotypes and phenotypes remain within site boundaries, with only summary
statistics or model parameters communicated between sites. *Federated*
fine-mapping is the natural extension of multi-ancestry methods to this
setting. Empirical evaluation of federated fine-mapping requires a
ground-truth simulation: a setting in which (a) the causal variants are known
by construction, (b) the per-population effect sizes are known, and (c) the
per-site cohorts can be assembled and re-assembled flexibly to probe how
ancestry composition, sample-size imbalance, and inter-site LD divergence
affect fine-mapping quality. This document describes the simulator we are
building to that purpose.

## 2 Synthetic cohort

We use the **HAPNEST** synthetic-genotype resource (EBI BioStudies S-BSST936),
which provides 1,008,000 simulated individuals partitioned uniformly across
six 1000-Genomes-style superpopulations (EUR, AFR, AMR, EAS, CSA, MID;
168,000 individuals per superpop). HAPNEST genotypes are generated by a
haplotype-conditional resampling procedure that preserves both local LD
structure and superpopulation-specific allele-frequency spectra, so the
resulting cohort is realistic at the locus scale while avoiding the
re-identification risks of real biobanks.

The on-disk layout we use is `synthetic_v1_chr-N.{bed,bim,fam}` for
N ∈ {1, …, 22} (PLINK-1.9 binary). A companion file
`synthetic_v1.sample` provides the superpopulation label for each of the
1,008,000 individuals, ordered identically to the rows of the `.fam` files.
We materialise a **population manifest** by row-wise joining the two files
into a TSV with columns `(FID, IID, superpopulation)`. The manifest is the
canonical mapping from individual ID to ancestral label and is used at every
downstream stage to look up the per-individual ancestry without re-reading
HAPNEST.

For the work described here we restrict to **chromosome 1** (`chr1`,
533,532 variants). Chromosome 1 is the largest human chromosome, contains
well-studied fine-mapping benchmark loci, and gives a tractable per-job
memory footprint (≈8 GB raw, ≈1.5 GB after MAF and per-site filtering).
Extending to additional chromosomes is a configuration change rather than a
code change.

## 3 Site cohorts

We define three federated *sites*, each enrolling 50,000 unique individuals
from the HAPNEST pool. The compositions are configured to mirror the
real-world ancestry mix at each partner institution:

| Site       | Locale            | n      | Ancestry composition (n) |
| ---------- | ----------------- | -----: | ------------------------- |
| **ANL**    | USA               | 50,000 | EUR 30 000, AFR 7 500, AMR 7 500, EAS 2 500, CSA 2 500 |
| **Covenant** | Nigeria          | 50,000 | AFR 47 500, EUR 1 500, CSA 1 000 |
| **MBZUAI** | UAE               | 50,000 | MID 25 000, CSA 15 000, AFR 5 000, EUR 5 000 |

**Disjointness invariant.** No individual is assigned to more than one site.
The total demand of 150,000 individuals is well below the HAPNEST cap of
1,008,000, and below the per-superpop cap of 168,000 for every
superpopulation (the largest per-superpop demand is AFR, at 7,500 + 47,500
+ 5,000 = 60,000, ≈ 36% of available AFR). Disjointness is enforced by a
deterministic, per-superpop shuffle of the HAPNEST individual list, followed
by a greedy partition in a fixed site order — `("anl", "covenant",
"mbzuai")`. The shuffle is seeded with
`SeedSequence(master_seed, "shuffle", <superpop>)`, so the partition is
reproducible across machines and across re-runs. After partitioning we verify
globally that the three site ID sets are pairwise disjoint and that the per-
site superpopulation tallies match the configured targets exactly. From the
2026-06-03 run:

```
Site anl:      50000 individuals (AFR=7500, AMR=7500, CSA=2500, EAS=2500, EUR=30000)
Site covenant: 50000 individuals (AFR=47500, CSA=1000, EUR=1500)
Site mbzuai:   50000 individuals (AFR=5000, CSA=15000, EUR=5000, MID=25000)
```

We materialise each site as a PLINK-1.9 binary fileset
(`{site}_chr1.bed/.bim/.fam`) via `plink --keep`, retaining the full chr1
variant list at this stage (variant pruning happens locus-by-locus during
locus selection). We also write a per-site manifest TSV `(FID, IID,
superpopulation)` so that downstream stages can recover each individual's
ancestry without re-joining against HAPNEST.

**Known defect: the site filesets are not allele-harmonized.** That `plink
--keep --make-bed` was issued *without* `--keep-allele-order`, so PLINK 1.9
re-assigned A1 to the minor allele **within each site's own cohort**. Because
the sites have different ancestry mixes, the same variant can end up coded the
opposite way round at different sites: against the HAPNEST coding, 9,398 /
31,031 / 9,777 of 533,532 chr1 variants are reversed at ANL / Covenant /
MBZUAI, and ~6% of variants differ in orientation between ANL and Covenant.

**Locus selection (§4) is unaffected**: it reads one site at a time and scores
windows on r² matrices, and r² is invariant to a flip applied to a whole
column.

**Phenotype simulation (§5) is affected, and this is the serious one.**
`phenotype_sim` computes each individual's genetic value as `β · x` from that
individual's *own site's* `.bed` (`_compute_genetic_value`), with no allele
handling anywhere in the module. Where a site codes a causal variant the other
way round, its individuals get `β·(2−x) = −β·x + 2β` — the opposite effect
direction. The core invariant of §5.2, *β is indexed by superpopulation and not
by site*, therefore does not hold for those variants: the realised effect
direction is indexed by site after all.

Measured on the current data package: **1,582 of 22,305 distinct causal
variants (7.1%) are coded inconsistently across the three sites**, touching
**1,595 of 11,850 manifest rows (13.5%)**. Verified directly on
`L0002 / ncsl1_h2-0.005_rg1 / rep 3` (causal `chr1:62560271:G:T`, flipped at
Covenant only), correlating each site's phenotype against the HAPNEST-coded
dosage over EUR individuals: ANL −0.054, MBZUAI −0.059, **Covenant +0.035**.

The consequence for any downstream per-ancestry analysis is *attenuation*, not
noise: a pooled ancestry column mixes `+β` and `−β` contributions in proportion
to its per-site sizes. Writing the retained fraction as
`Σ_k ±n_{p,k} / Σ_k n_{p,k}`, the affected variants retain a median of 0.918
(EUR), **0.583 (AFR)** and 0.892 (CSA) of their marginal signal. AFR is worst
because Covenant holds 47,500 of the 60,000 AFR individuals and carries the
most flips (1,286 of the causal variants), so a Covenant-only flip inverts the
sign of the pooled AFR marginal relative to the stored β and leaves 58% of its
magnitude. Single-site ancestries (AMR, EAS, MID) are untouched: a uniform flip
within one site is just `β → −β` consistently, which changes no PIP and no
credible set. Nothing is fully cancelled, and because SuSiEx scores variants
through `z²` a per-population sign inversion costs nothing by itself — the cost
is the lost magnitude, i.e. power, on ~13.5% of instances. It is paid equally
by both fine-mapping paths, so it does not threaten the §11.3 equality result;
it depresses the absolute power numbers both of them report.

**For the federated path it is also a correctness issue**, because that path is
the one place two sites' genotype matrices are *combined*: a flipped site
contributes `2 − x` to a summed Gram where the others contribute `x`, and that
sums in silently ([`algo.md` §7](algo.md) condition (i); [§A.10](algo.md)
trap 2). `src/fed_fine_mapping.py` therefore recodes each site's dosages
against a published canonical variant list before forming any moment, and the
coordinator re-checks. The centralized baseline never notices, because it cuts
its window straight from the HAPNEST fileset.

**The fix is to re-cut the site filesets with `--keep-allele-order` and re-run
phenotype simulation.** Locus selection can resume from its cache (its inputs
do not change). Until that happens: any new code reading more than one site's
`.bed` must harmonize first, and power numbers on the affected 13.5% of
instances are understated.

## 4 Locus selection

Across the three sites we select 100 loci on chr1, **stratified by
inter-site LD divergence**. The intuition is that fine-mapping resolution
is bounded by LD: at a locus where all three sites have nearly identical LD
structure, multi-ancestry methods reduce to single-ancestry methods, while
at a locus where LD differs substantially across sites we expect
multi-ancestry methods to materially shrink credible sets. By stratifying
the candidate locus list along this axis we obtain a benchmark that
characterises fine-mapping performance as a function of the very property
that multi-ancestry methods are designed to exploit.

### 4.1 Candidate windows

We tile chr1 with **1 Mb sliding windows in 500 kb steps** (50% overlap),
producing a candidate list of ≈500 windows. Window size is the smallest
unit for which LD divergence is meaningful at modern variant density, and is
matched to the typical SuSiEx working-window scale.

### 4.2 Per-window LD computation

For each candidate window, and for each of the three sites independently, we
compute a tag-SNP LD matrix as follows:

1. **MAF filter** at the *site-specific* threshold of 0.01.
2. **Greedy LD prune** to a target of ≈500 tag SNPs per window, using a
   variant-count window of 50 with step 5 and an r² pruning threshold of
   0.995. The pruner walks the window left-to-right and retains a candidate
   variant if its r² with every previously retained variant in a sliding
   neighbourhood is below threshold.
3. **r² matrix** over the surviving tag SNPs, computed as
   `R = (Xᵀ X) / n_samples` after column-wise z-scoring with mean imputation
   of missing dosages. We use the **population variance convention**
   (divisor `n`, not `n − 1`), so the diagonal would equal exactly 1.0 in
   the no-missing-data case. In practice HAPNEST chr1 carries a small
   missing-genotype rate (≈0.4%) so the empirical diagonals sit just
   below 1 (median deviation ≈ 4 × 10⁻³ in the realised data). This is
   the same convention as `plink --r2 --matrix` and is what subsequent
   Frobenius-norm comparisons assume.

The r² matrix size is bounded by the tag-SNP target (≈500 × 500 = 250k
float64 entries ≈ 2 MB per site per window), so the full per-window memory
footprint across three sites is on the order of 6–10 MB, comfortably below
any node-RAM budget. The expensive operation is the underlying genotype
read; we use a minimal pure-NumPy `.bed` reader that decodes the four
PLINK bit codes (00 = hom-A1 → dosage 2, 01 = missing → NaN,
10 = het → 1, 11 = hom-A2 → 0) and accumulates dosages column-wise.

### 4.3 Pairwise divergence and stratification

For each window we compute the three pairwise Frobenius distances
`‖R_i − R_j‖_F` between site LD matrices, and summarise the window with
the mean of those three distances. The resulting per-window divergence score
is binned into terciles (low / medium / high), and 33 / 33 / 34 windows are
sampled without replacement from the three terciles respectively. A final
overlap-removal pass drops any selected window whose physical interval
overlaps a previously selected one. Window-level computation is
embarrassingly parallel and is dispatched with `joblib.Parallel` (loky
backend); the configured `n_workers` is 8 on Sophia and 32 on a Polaris
compute node.

The overlap-removal pass operates in two passes. The first pass, which is
called inline by `locus_selection.py`, deduplicates *within* each stratum:
it sorts each tercile by within-stratum extremeness (highest score in the
high tercile, lowest in the low tercile, closest-to-median in the medium
tercile) and admits windows greedily. The second pass, applied as a
post-processing step, removes residual *cross-stratum* overlaps: windows
selected for different terciles can still share physical intervals (e.g.
a high-divergence window adjacent to a medium-divergence window), so we
do a global greedy sweep along the chromosome and drop the window of lower
within-stratum extremeness from each remaining overlap pair. On the chr1
realised data this two-pass dedup yields 79 non-overlapping
loci, with strata 25 high / 26 medium / 28 low. Per-stage caching of the
LD-matrix `.npz` files means re-running locus selection after configuration
changes does not re-pay the genotype-read cost.

The 1 Mb window with 500 kb step deliberately produces 50%-overlapping
candidate windows; the inline + post-processing dedup is the price of that
overlap. An alternative tile would be non-overlapping 1 Mb windows, which
would eliminate the need for either pass at the cost of halving the
candidate pool.

## 5 Causal architectures

At each locus we instantiate a grid of causal architectures, parameterised
by three factors:

* **ncsl** ∈ {1, 2, 3} — number of causal SNPs per locus.
* **h²** ∈ {5 × 10⁻⁴, 10⁻³, 5 × 10⁻³} — per-locus narrow-sense heritability
  (so per-genome heritability scales with the number of selected loci).
* **rg** ∈ {0.5, 0.7, 1.0} — cross-ancestry genetic correlation of effect
  sizes at the causal SNPs.

Three factorial modes are supported:

* **minimal** (9 architectures): ncsl × h², with rg fixed at 1.0.
* **extended** (15): minimal plus ncsl × rg at h² = 10⁻³, deduplicated.
* **full** (27): the full ncsl × h² × rg Cartesian product.

The default is **extended**, which retains the full ncsl-and-rg sweep
relevant for fine-mapping while collapsing the heritability axis to a single
mid-range value when rg is varied. Each architecture is realised in **10
independent replicates** that re-draw the causal-SNP positions and the
per-population effect-size vector under different seed components. The total
realised architecture count is therefore 79 loci × 15 architectures ×
10 replicates = 11,850 phenotype simulations per site.

### 5.1 Causal-SNP selection

Within each locus, causal SNPs are drawn from the locus's tag-SNP list with
two constraints intended to make the simulated GWAS signal realistic and
non-degenerate:

* **Common-vs-rare split.** Candidate variants are partitioned into a
  common stratum (MAF ≥ 0.05) and a rare stratum
  (0.01 ≤ MAF < 0.05). At least one causal variant is drawn from each
  stratum whenever `ncsl ≥ 2`, so that no replicate concentrates all
  signal in the common-variant regime where single-ancestry fine-mapping is
  already strong.
* **Ancestry-specific availability.** Each causal SNP must have
  MAF ≥ 0.01 in every superpopulation present at the cohort sites. This
  prevents the pathological case in which a variant is monomorphic in a
  site and thereby trivially absent from that site's GWAS signal.

### 5.2 Per-ancestry effect-size draw

The core correctness invariant of this simulator is that **effect sizes
are indexed by the individual's superpopulation, not by the site**.

> **Caveat (found 2026-08-20).** In the *realised* data package this invariant
> is broken for 13.5% of instances, not by this stage's logic but by the
> allele-order defect in the site filesets it reads (§3): where a site codes a
> causal variant the other way round, its individuals receive `−β` rather than
> `+β`. The draw described below is correct; the dosages it is multiplied
> against are not consistently coded. §3 has the measured blast radius and the
> fix.

For
each causal SNP we draw a 6-vector β ∈ ℝ⁶ — one entry per HAPNEST
superpopulation — from a multivariate Gaussian with covariance matrix Σ
where Σ_ii = 1, Σ_ij = rg (i ≠ j). The cross-population
correlation rg is the simulator's principal degree of freedom for
modelling cross-ancestry effect-size heterogeneity. We sample via Cholesky
decomposition with a small ridge (`+ 1e-8 I`) at rg = 1.0 to keep the
factorisation numerically stable at the boundary of the PSD cone.

The genetic value of individual *i* with genotype vector
g_i ∈ ℝ^{ncsl} and superpopulation label
pop(*i*) ∈ {1, …, 6} is

> y_i^(g) = Σ_k g_{i,k} · β_{k, pop(i)}

i.e. the same genotype contributes differently to the phenotype of an EUR
individual and an AFR individual at the *same* site. This per-individual
ancestry lookup is what makes the simulation a meaningful test bed for
multi-ancestry fine-mapping: any method that assumes a single site-level β
will produce systematically inflated or attenuated effect-size estimates
relative to ground truth, in proportion to the within-site ancestry
diversity. Crucially, the lookup is by **ancestry**, not by **site**: an
AFR individual at ANL and an AFR individual at Covenant share the same β
entry by construction.

### 5.3 Phenotype noise

Phenotypes are generated under a linear-additive model
y_i = y_i^(g) + ε_i, with ε_i ∼ 𝒩(0, σ²_ε). The noise variance is
calibrated *per site* so that the *empirical* heritability on that site
matches the target h². Concretely, we measure
σ²_g = Var(y^(g)) over the site's individuals and set

> σ²_ε = σ²_g · (1 − h²) / h²

so that h² = σ²_g / (σ²_g + σ²_ε) by construction. We then assert that
the realised empirical heritability is within a relative tolerance
(`h2_tolerance_relative = 0.20`) of the target; this slack absorbs
finite-sample variance in the genetic-variance estimator without admitting
gross mis-specifications.

### 5.4 Ancestry-divergent causal sets

The baseline model of §5.2 shares a single causal SNP set across all
superpopulations and lets only the *effect sizes* diverge (via `rg`). An
optional **ancestry-divergent** mode instead lets the causal *set itself*
differ across superpopulations — the regime in which multi-ancestry
fine-mapping has its largest theoretical advantage, because a variant that is
causal in one ancestry is simply absent from the signal in another.

A divergent instance draws a **union** causal set of `n_shared` variants
common to all superpopulations plus `n_private_per_pop` variants private to
each, where `n_shared = ncsl − n_private_per_pop` (so each superpopulation
still has exactly `ncsl` causal variants). We represent this with the *same*
`(n_union, n_pops)` effect-size matrix as the baseline, but with **structural
zeros**: a shared variant carries the usual MVN(0, Σ_rg) draw across all pop
columns, while a private variant of pop *p* carries a single 𝒩(0, 1) entry in
column *p* and zero elsewhere. Because the genetic value is the per-individual
lookup `g_i = Σ_k g_{i,k} · β_{k, pop(i)}` of §5.2, the zeros automatically
restrict each individual to their own superpopulation's causal set — no change
to the genetic-value or phenotype machinery is required, and the §10.5
bit-exact re-derivation check applies unchanged.

The mode is **distinct from** the ancestry-specific-MAF override also present
in the code (`ancestry_specific_causal`), which diverges a single variant's
*allele frequency* across sites while keeping the causal set shared; the two
are orthogonal and independently toggled. Following the same design as that
override, divergent instances are selected by a global seeded draw of
`min_per_stratum` (locus, architecture) pairs per divergence stratum, computed
on the full locus set before any sharding so single-node and multi-node runs
agree. Only architectures with `ncsl > n_private_per_pop` are eligible (a
shared variant must remain); otherwise the instance falls back to a shared
set, recorded as `causal_mode = shared` in the manifest. The mode was
**enabled on 2026-07-15** (`ancestry_divergent_causal.enabled = true`) after an
end-to-end sanity check on real genotypes (3 flagged loci; union sets and
structural-zero β verified). A prod re-run to fold the regime into the data
package is **in progress**; until it lands, the realised
on-disk package remains the shared-only 2026-06-17 data, so the §10 numbers
below predate divergence and will be refreshed after the re-run.

## 6 Quality control

Three QC stages run after sampling and phenotype simulation, with HTML
reports written to `reports/`:

* **Relatedness** via `plink2 --make-king-table`, flagging pairs with
  kinship coefficient above 0.05 (3rd-degree or closer). HAPNEST individuals
  are simulated to be unrelated, so this is expected to be empty; we run it
  as a sanity check on our partitioning code.
* **Per-superpop MAF concordance** — for each site and each superpopulation
  present at that site, we tabulate the MAF distribution and check that the
  empirical MAFs of common SNPs deviate from the corresponding HAPNEST
  reference by no more than 3 standard deviations.
* **LD decay** — empirical r² vs. physical distance, binned into 1 kb
  intervals up to 1 Mb, computed on a random sample of five loci per site.
  This is a fingerprint of the local LD architecture and is expected to
  agree closely with the HAPNEST source LD; gross deviations would indicate
  a bug in the genotype-read or LD-computation paths.

## 7 Reproducibility infrastructure

### 7.1 Deterministic seeding

A single 64-bit `master_seed` (20260601 in the default configuration) is
threaded through every stochastic decision. Per-decision sub-seeds are
derived as

> `np.random.SeedSequence(master_seed, *components).spawn(1)[0]`

where the variadic `components` argument is hashed using a stable
FNV-1a-based string hash applied to each non-integer element. This makes
the sub-seeds independent of Python's per-process `hash()` randomisation
and reproducible across machines. Concretely, the sub-seed for the
per-superpop shuffle in §3 is
`derive_seed(master_seed, "shuffle", "EUR")` and so on. We have unit tests
asserting bit-exact reproducibility of every sub-seed.

### 7.2 Configuration

All parameters described above live in a single YAML file
(`config/simulation_config.yaml`) validated by a Pydantic v2 schema at
load time. Source code contains no magic numbers; every threshold, count,
and path is sourced from the validated config object. This is what makes
the pipeline a configuration-bound experiment harness rather than a
hand-tuned script.

### 7.3 Software environment

The pipeline depends on Python ≥ 3.10 with numpy, pandas, scipy, pydantic,
PyYAML, joblib, matplotlib, and pytest; and on external binaries
`plink` (1.9, build 2024-10-22) and `plink2` (2.0, build 2025-01-29).
A `scripts/install_plink.sh` helper downloads static Linux binaries into
`vendor/bin/` so the pipeline runs without root or conda-environment
mutation. `scripts/_env.sh` resolves a Python interpreter with the
required scientific stack at job start, preferring the ALCF conda module
(absolute path `/soft/applications/conda/2024-08-08/mconda3/bin/python`)
and falling back to `$PYTHON` or `$PATH` if that is unavailable. The same
`_env.sh` works on both Sophia (interactive / development) and Polaris
(batch) because `/soft` is mounted identically on both systems.

### 7.4 Tests

The simulator ships with 82 unit tests across six modules:

* `test_sampling.py` (6) — disjointness, per-superpop tallies, seed
  determinism, manifest schema.
* `test_locus_selection.py` (12) — candidate-window tiling, MAF filter,
  greedy LD prune, r² matrix diagonal = 1, Frobenius norm symmetry,
  stratification quantile boundaries, overlap drop.
* `test_phenotype_sim.py` (17) — including a critical
  `test_beta_by_superpopulation_lookup` that verifies a mixed-ancestry
  cohort produces phenotypes consistent with per-individual β lookup
  (and inconsistent with a single site-level β).
* `test_utils.py` (5) — config-loading errors, `derive_seed`
  reproducibility, PLINK `.bed` bit-code decoding.
* `test_fine_mapping.py` (14) — SuSiEx `.cs` / `.snp` / `.summary` parsing
  (perfect recovery, multi-CS, the "no credible set" case), the ancestry
  column planner, LD-file readiness, binary resolution.
* `test_fed_fine_mapping.py` (28) — the federated path (§11.3). Unit layer:
  additivity of the site aggregates over an arbitrary row partition
  ([`algo.md` §7](algo.md) Theorem 1), the pooled-vs-local standardization
  worked example of [`algo.md` §A.10](algo.md) (0.816 against 0.866), the
  closed-form pooled GWAS against `numpy.linalg.lstsq`, the SuSiEx LD-panel
  byte layout, and the allele-flip and missing-genotype harmonization rules.
  Acceptance layer: both fine-mapping paths run end-to-end over a miniature
  three-site package and their credible sets and PIPs compared — see §11.3.

The full suite passes in ~7 s on a single Sophia core. Of the 28 federated
tests, 17 are pure Python and 11 shell out to PLINK 1.9, PLINK 2 and SuSiEx to
build the centralized comparator; those 11 skip if the binaries are absent.

## 8 Execution model

The **simulator** is structured as four stages, each invokable as
`scripts/run_<stage>.sh [config_path]`:

1. **sampling** — read population manifest, partition individuals across
   sites, run PLINK `--keep` to materialise per-site `.bed/.bim/.fam`.
2. **locus_selection** — tile chr1, compute per-site LD matrices, score
   divergence, stratify and sample 100 loci, write per-locus tag-SNP
   lists and a summary TSV.
3. **phenotype_sim** — for each (site, locus, architecture, replicate),
   draw per-ancestry effect sizes, compute genetic values via the
   per-individual ancestry lookup, calibrate noise to target h², write
   phenotype CSVs and a ground-truth manifest.
4. **qc** — relatedness, MAF concordance, LD decay; emit HTML reports.

Two **analysis** stages consume the finished package and are run separately,
after it exists (§11.3). Both shard by locus with the same
`--shard-index / --n-shards / --merge` contract as `phenotype_sim`, and both
must be launched through a thin script rather than `python -m`, because loky
workers import the instance function by qualified module name:

5. **fine_mapping** (`scripts/run_fine_mapping.py`) — the centralized SuSiEx
   baseline: six ancestry columns pooled across sites, one GWAS and one
   in-sample LD panel each, cut from a single per-locus window.
6. **fed_fine_mapping** (`scripts/run_fed_fine_mapping.py`) — the federated
   path: per-(site, ancestry) raw moments, summed and standardized at a
   coordinator that never sees a genotype.

On Sophia (the interactive ALCF system used for development) we run all four
simulator stages back-to-back as a single shell pipeline. On Polaris (the production
batch system) we submit a single PBS-Pro job
(`scripts/submit_polaris.pbs`) that skips the sampling stage (whose output is
already on the shared parallel filesystem) and runs locus selection,
phenotype simulation, and QC inside one `preemptable`-queue allocation of
one node, 32 cores, 5-hour walltime. A bumped configuration
`config/simulation_config_polaris.yaml` is generated at job start to raise
`locus_selection.n_workers` from 8 to 32. We initially attempted the
`debug` queue (1 h cap, ≤2 nodes); the first run reached 386 of 496
candidate windows in 47 min and was wall-killed at the 1 h mark, motivating
the move to `preemptable`. We measured `n_workers=64` to be roughly 2×
*slower* than `n_workers=32` on this workload, consistent with the LD
computation being memory-bandwidth-bound on the single-NUMA-domain EPYC
node, and reverted to 32 workers.

A second design decision motivated by the failed `debug` run was to add
resume-from-cache logic at the top of `score_window`: each per-site r²
matrix is written to `data/loci/ld_matrices/{window_id}_{site}.npz` as
soon as it is computed, and on re-entry the function reads the cache and
skips straight to divergence scoring if all three site files are present.
This makes the locus-selection stage restart-safe under preemption.

## 9 Status as of 2026-06-17

The pipeline has been run end-to-end on Polaris (a clean full recompute on
the sharded `prod` queue, superseding the earlier 2026-06-09 single-node run)
and the resulting data package validated against the design invariants below.
Specifically:

* **Sampling: complete.** The three per-site PLINK filesets exist at
  `data/processed/{anl,covenant,mbzuai}/{site}_chr1.{bed,bim,fam}` with
  exactly 50,000 individuals and 533,532 variants each, matching the
  configured compositions to the individual.
* **HAPNEST: complete.** The full chr1 fileset and population manifest
  (1,008,000 rows) are materialised under `data/raw/hapnest/`.
* **Locus selection: complete.** Of 496 candidate 1 Mb / 500 kb-step
  windows, 458 were scorable and were binned into terciles by the
  stratified draw; the two-pass overlap-removal (§4.3) then yielded
  **79 non-overlapping loci** (25 high / 26 medium / 28 low). The cached
  `.npz` r² matrices for all 458 scored windows (1,374 files across the
  three sites) are preserved under `data/loci/ld_matrices/` for downstream
  LD-aware analyses.
* **Phenotype simulation: complete.** 79 loci × 15 architectures × 10
  replicates × 3 sites = **35,550 phenotype files** (11,850 per site) under
  `data/ground_truth/phenotypes/{site}/L*_*.pheno`, plus matching β
  arrays under `data/ground_truth/effect_sizes/` and a unified
  `data/ground_truth/causal_manifest.tsv` with one row per (locus,
  architecture, replicate).
* **QC: complete.** Three HTML reports under `reports/qc_{site}.html`
  covering relatedness, per-superpop MAF concordance, and LD decay.
* **Wall-clock budget.** The clean recompute ran in ~23 min end-to-end on
  the `prod` queue (10 nodes, `place=scatter`, each embarrassingly-parallel
  stage sharded across nodes with 32 joblib workers each), versus ~3 h 40 min
  for the earlier single-node run. QC (dominated by
  `plink2 --make-king-table` on 50k-individual cohorts) runs concurrently on
  three dedicated nodes and is the long pole.

## 10 Empirical validation

We validated the resulting data package by checking each design invariant
against the realised outputs.

### 10.1 Heritability calibration

For every (locus, architecture, replicate, site) triple, the simulator
reports `empirical_h2_<site>` alongside the target `h2_target`. The relative
error |empirical / target − 1| should lie within the configured
`h2_tolerance_relative = 0.20`. In practice it is well inside that bound
at every site:

| Site | Median rel. err | 95th pctile | Max |
|---|---:|---:|---:|
| ANL | 0.0042 | 0.0126 | 0.0246 |
| Covenant | 0.0042 | 0.0123 | 0.0286 |
| MBZUAI | 0.0043 | 0.0124 | 0.0271 |

Median error is ~0.4% across all 11,850 manifest rows; the largest
single-replicate deviation is under 3%. The configured 20% tolerance is
therefore conservative by more than an order of magnitude.

### 10.2 Cross-ancestry effect-size correlation

For each (locus, architecture, replicate) we draw the 6-vector of
per-superpop effect sizes from 𝒩(0, Σ) with Σ_ii = 1, Σ_ij = rg.
Aggregating the β rows of each `rg` architecture across the realised data
and computing the 6 × 6 sample correlation matrix yields mean off-diagonal
correlations that track the target closely: **0.491** at rg = 0.5 (range
[0.475, 0.503], |err| 0.009) and **0.699** at rg = 0.7 (range [0.686, 0.708],
|err| 8 × 10⁻⁴). At `rg = 1.0` (the singular boundary, sampled via Cholesky
with a `+10⁻⁸ · I` ridge) the realised β rows are identical across
superpopulations (mean off-diagonal 1.0000 to reported precision), i.e.
floor noise from the ridge.

### 10.3 Structural completeness

We verified the file inventories and manifest shapes match the
combinatorics of the configured grid:

| Invariant | Observed | Expected |
|---|---:|---:|
| Selected loci | 79 | (≤ 100 by config, after overlap dedup) |
| Phenotype files per site | 11,850 | 79 × 15 × 10 |
| Total phenotype files | 35,550 | 11,850 × 3 |
| Effect-size `.npy` files | 11,850 | 79 × 15 × 10 |
| Manifest rows (causal + seeds) | 11,850 each | 79 × 15 × 10 |
| ncsl matches `\| causal_snp_ids \|` | 11,850 / 11,850 | all match |
| Causal SNP bp within locus window | 11,850 / 11,850 | all in window |
| Phenotype rows per file | 50,000 | per-site n |
| Phenotype y has no NaN | passes | finite |

All structural invariants hold exactly.

### 10.4 Tag-SNP r² matrices

Cached `.npz` r² matrices were inspected for numerical sanity. With
HAPNEST's ≈ 0.4% missing-genotype rate and mean imputation of missing
dosages, the diagonal sits at 1 − (n_missing / n_total). The realised
median diagonal deviation across all 1,374 cached files is 4.6 × 10⁻⁴
(max 5.4 × 10⁻³), consistent with the missing-data rate, and no NaN
entries are present. The same mean-imputation floor means the imputed
columns are not exactly unit-variance, so a small number of off-diagonal
entries marginally exceed 1 (realised max |R| = 1.0043, i.e. 4 × 10⁻³
above unity — the same floor as the diagonal deviation). This corrects an
earlier statement that no entry exceeds 1; the overshoot is bounded by the
missing-data floor and well within the 1 % tolerance applied by the
validation harness. Per-window tag-SNP counts agree across the three sites
by construction.

### 10.5 Bit-exact phenotype re-derivation

The per-individual β lookup has now been verified end-to-end against the
realised production data, closing the gap noted in earlier revisions. For
every (locus, architecture, replicate, site) instance we independently
recompute the genetic value

> g_i = Σ_k g_{i,k} · β_{k, pop(i)}

directly from the raw `data/processed/{site}/{site}_chr1.bed`, the per-
replicate `data/ground_truth/effect_sizes/*.npy`, and the per-individual
ancestry column of `data/processed/{site}/{site}_manifest.tsv` (aligned to
`.fam` order), using a reimplementation of the lookup that does **not**
import the production helper. Because the realised phenotype `y = g + ε`
carries a noise term whose RNG state is not persisted, `y` itself is not
reproducible byte-for-byte; the deterministic, checkable quantity is the
genetic component. The simulator recorded
`empirical_h2_<site> = Var(g)/Var(y)`, so the bit-exact tie is

> Var(g_rederived) / Var(y_stored)  ≟  manifest.empirical_h2_<site>

Across all **35,550** (site × locus × arch × replicate) instances the
re-derived ratio reproduces the recorded empirical heritability to a median
relative error of ≈ 2.4 × 10⁻⁸ and a worst case of 1.6 × 10⁻⁷, far inside
the 10⁻⁵ tolerance.

To confirm the check actually exercises the *per-ancestry* lookup and not
merely some additive model, we recompute `g` under a deliberately wrong
**site-level** model that assigns every individual the site's dominant-
superpopulation β column, and verify it fails to reproduce the manifest
heritability wherever it genuinely differs. As expected from the
construction, at rg = 1.0 the per-population β columns are identical (to the
cholesky-ridge floor) so the two models coincide and the site-level model is
*not* rejected (rejected on only 2.3 % of rg = 1.0 instances, i.e. ridge
floor); at rg ∈ {0.5, 0.7} the site-level model is materially wrong on
99.9 % of the ancestry-mixed (ANL, MBZUAI) instances. The per-individual
lookup is therefore only *discriminated* by the rg < 1 architectures — a
useful caveat for downstream experiment design.

This check is one of five in the exhaustive validation harness
(`src/validation.py`, run via `scripts/submit_polaris_validation.pbs`),
which re-runs §10.1–§10.4 alongside it and writes per-instance CSVs plus a
`reports/validation/SUMMARY.txt` rollup.

## 11 Next steps

In rough order of priority:

1. ~~Add the bit-exact phenotype re-derivation described in §10.5 to the
   validation harness.~~ **Done** — implemented as `src/validation.py`
   (check 1 of 5) and run on Polaris via
   `scripts/submit_polaris_validation.pbs`; see §10.5.
2. ~~Add an explicit "ancestry-divergent causal" mode in which the causal
   variant set is allowed to differ across superpopulations (in addition
   to the effect-size correlation already captured by rg).~~ **Built** —
   see §5.4 (`ancestry_divergent_causal`, union-set + masked-β
   representation; validated by two new structural checks in
   `src/validation.py` and unit tests in `tests/test_phenotype_sim.py`).
   **Enabled on 2026-07-15** and sanity-checked on real genotypes; the prod
   re-run to fold the regime into the data package is in progress (see §5.4).
3. ~~Wire SuSiEx into a downstream `src/fine_mapping.py` that consumes the
   per-site phenotypes and emits credible sets.~~ **Built (2026-07-15)** —
   `src/fine_mapping.py` runs, per (locus, architecture, replicate), a GWAS per
   population column (`plink2 --glm` on the locus window, which also serves as
   the LD reference panel), reformats to SuSiEx summary-statistics format, runs
   SuSiEx across those columns, and scores the returned credible sets against
   the ground-truth causal SNPs in `causal_manifest.tsv` — emitting standard
   fine-mapping metrics (credible-set size, PIP of the true causal, capture
   rate, purity) rolled up by architecture and by (LD-divergence stratum, rg).
   SuSiEx is a C++ CLI we shell out to (installed by `scripts/install_susiex.sh`
   into `vendor/bin/`), not a Python import; the federated formulation it
   realises is derived in [`algo.md`](algo.md). Proven end-to-end (perfect
   single-causal recovery on the strong-signal class); the full-package run is
   pending the divergence re-run of step 2. The `data/loci/tag_snps/` artifact
   referenced in earlier drafts is unused — windows are taken directly from the
   per-site filesets by base-pair range.

   **Updated (2026-07-27): the SuSiEx population columns are now one per
   superpopulation, pooled across sites.** Earlier drafts fed SuSiEx one column
   per *site*, but a site cohort is ancestry-mixed, so its marginal effect is a
   frequency-weighted blend of several true per-ancestry effects and its LD is
   ancestry-averaged — misspecified against the per-superpopulation effect model
   of §5.2 and blunting the very cross-ancestry LD contrast the locus strata are
   built on. The stage now pools each ancestry across every site holding it
   (AFR = 7,500 ANL + 47,500 Covenant + 5,000 MBZUAI = 60,000) into **six
   columns covering all 150,000 individuals exactly once**: EUR 36,500,
   AFR 60,000, AMR 7,500, EAS 2,500, CSA 18,500, MID 25,000. Each column gets
   its own GWAS and in-sample LD panel over its full pooled cohort, both cut
   from a single per-locus window with `plink --keep`; `min_gwas_n` drops any
   ancestry whose pooled cohort is too small to inform a GWAS.

   This is the **centralized baseline** — pooling by ancestry with no regard to
   site boundaries is exactly the fit obtained if all three sites' data lived in
   one place, and by [`algo.md` §8](algo.md) Corollary 1 it is what a correct
   federated implementation must reproduce. It is therefore the upper bound
   against which federated variants are reported. Two caveats on that ceiling:
   pooling only helps ancestries that span sites (AMR and EAS are ANL-only and
   gain nothing), and it presumes the shared-causal-variant assumption that the
   §5.4 ancestry-divergent instances deliberately violate.

   **Built (2026-08-20): the federated path.** `src/fed_fine_mapping.py`
   implements [`algo.md` §7 and §9](algo.md) as a second, independent module
   alongside the centralized one — both exist and run separately, since the
   centralized module is the comparator the federated one is validated against
   and changing it would void the comparison. Each site reads only its own
   `data/processed/{site}/{site}_chr1.{bed,bim,fam}`, manifest and `.pheno`
   files and emits, per (site, ancestry), the raw unstandardized moments
   `G = X'X`, `c = X'y`, `u = 1'X`, `n`, `q = 1'y`, `w = y'y` — twelve blocks
   (ANL holds five ancestries, Covenant three, MBZUAI four). Nothing else
   crosses the site boundary. The coordinator sums each ancestry's blocks and
   standardizes **once** against the pooled moments (eq 9.1–9.2), including the
   MAF filter; standardizing at the sites instead is biased, which is the one
   real trap of the design ([`algo.md` §9](algo.md), worked in §A.10). Having
   no genotypes, it cannot call `plink --r`, so it writes SuSiEx's LD panel
   files itself in the layout the C++ reader expects, and derives each
   ancestry's `BETA`/`SE`/`P` in closed form from the moments
   ([`algo.md` §9.1](algo.md)). SuSiEx is then invoked exactly as the
   centralized stage invokes it. Entry point
   `scripts/run_fed_fine_mapping.py`, sharded by locus like its twin; outputs
   `reports/fed_fine_mapping/fed_fm_results.tsv` with the same schema as
   `fm_results.tsv`, so the two are directly comparable.

   **Agreement measured.** `tests/test_fed_fine_mapping.py` runs both paths
   over a miniature three-site package and asserts identical credible sets and
   equal PIPs; the residual decomposes into a summary-statistic half that is
   *bit-exact* (asserted with zero tolerance) and an LD half that agrees to ≤2
   units in the last place of the float32 both paths store `R_p` in. On the
   production package, one instance (L0000, `ncsl1_h2-0.0005_rg1`, rep 0) gave
   summary statistics byte-identical to `plink2 --glm` across all 9,848
   variants of the six ancestry columns, the same single-variant credible set,
   and max |ΔPIP| = 6.9×10⁻¹⁶ over 2,033 variants. See
   [`algo.md` §11.1](algo.md). **Not yet run at scale** — one locus, one
   instance; the 11,850-instance federated run is the next step, and the
   deliverable there is federated-vs-centralized agreement over the whole grid,
   not a performance gap.

   Two harmonization conditions from [`algo.md` §7](algo.md) had to be handled
   explicitly, and both are properties of the data package rather than of the
   estimator: the site filesets are not allele-harmonized (§3), and the
   comparator's per-variant / pairwise missing-data handling cannot be
   expressed in one `(G, u, n)` triple, so a site drops any window variant it
   cannot observe completely and reports the count (at L0000: 6 / 1 / 3 / 1
   window variants from the EUR / AMR / CSA / MID columns). Where that count is
   nonzero the two paths fine-map slightly different variant sets and
   Corollary 1's equality is no longer exact.

   Cost, measured on Sophia for L0000 (M = 2,046 variants, 150,000
   individuals): the site stage is ~8.5 min per locus for all twelve blocks and
   is *amortized over every instance at that locus*, since `G`, `u` and `n` do
   not depend on the phenotype; the coordinator then spends ~9 s per instance.
   The centralized comparator spends ~21 s of per-locus setup and ~12 s per
   instance, so at the configured 150 instances per locus the two are within a
   few percent of each other overall.
4. **Deferred (decision 2026-07-15): stay on chr1 until the very end.** The
   extension to chr1+chr11 and then the full autosome is a configuration change
   that scales linearly with chromosome length (~30× chr1 chr-wide), held off
   deliberately so the causal-architecture and fine-mapping methodology is fully
   settled on the single-chromosome package first.
5. **Not needed for the synthetic package; required before any real-data
   deployment: harden the aggregate uplink.** Everything in this pipeline runs
   on HAPNEST synthetic genotypes, so no disclosure control is warranted today
   and the aggregates of [`algo.md` §9](algo.md) may be shipped in the clear —
   which, as of 2026-08-20, is literally what `src/fed_fine_mapping.py` does:
   no secure aggregation, no DP, the twelve blocks move as plain arrays inside
   one process. Everything below is therefore unimplemented by design, and is
   the gap between this simulator and a deployment.
   The moment the same protocol is pointed at real participants it needs the
   following, because the exactness fix of §9 — ship *raw, unstandardized*
   moments and standardize once at the coordinator — is also a **privacy
   regression** relative to shipping $(R_p, \hat\beta_p, n_p)$ directly.
   Standardization is lossy, and what it discards is the column means, i.e. the
   allele frequencies; sending the pre-image sends those too. In rough order of
   value:

   a. **Single-site ancestries should never form a Gram.** AMR and EAS are
      ANL-only and MID is MBZUAI-only (§3), so for those columns the holding
      site already has the complete pooled cohort and can compute the minimal
      sufficient statistic $T_p = (R_p, \hat\beta_p, n_p, \tau_p^2)$
      ([`algo.md` §6](algo.md)) itself. Shipping $T_p$ instead of
      $(G_{p,k}, c_{p,k}, u_{p,k}, \dots)$ is identical in result, strictly
      lower in disclosure, and puts no allele frequencies on the wire. Three of
      the twelve (site, ancestry) blocks are covered by this alone.

   b. **Secure aggregation for the nine cross-site blocks.** EUR, AFR and CSA
      span all three sites. The coordinator only ever needs $\sum_k G_{p,k}$,
      never the per-site terms, so pairwise additive masks that cancel in the
      sum reduce its view to the pooled per-ancestry total. Because the
      protected operation is literally addition, this preserves Theorem 1
      ([`algo.md` §7](algo.md)) bit-for-bit — the federated fit stays *exactly*
      equal to the centralized baseline of §11.3. It also hides the small
      blocks inside large ones: Covenant's 1,000-individual CSA block
      disappears into an 18,500 pooled total. This is the main item.

   c. **Freeze and version cohorts; never re-release a block whose membership
      changed.** If the same (site, ancestry) block is emitted twice with one
      individual added, $G_2 - G_1 = xx^\top$, and since the shipped dosages are
      raw and non-negative, $x_j = \sqrt{(G_2-G_1)_{jj}}$ recovers that person's
      genotype across the whole locus exactly. The $O(M)$ vectors are worse:
      $u_2 - u_1 = x$ directly. This is a release policy, not a code change.

   d. **Minimum block size.** Blocks with $n$ below roughly the variant count
      are the reconstruction-risk regime — $G$ is rank-deficient and genotypes
      lie on the $\{0,1,2\}$ lattice. With $M$ averaging 2,290 per locus
      (range 1,155–4,125; §4.1) that flags Covenant CSA ($n=1{,}000$) and
      Covenant EUR ($n=1{,}500$); merge or drop them.

   **Explicitly not recommended: differential privacy on $G_{p,k}$.** A Gaussian
   mechanism at $\varepsilon=1$ needs per-entry noise on the order of the Gram's
   $\ell_2$ sensitivity $\lVert x\rVert^2 \approx 0.48M$, which propagates to
   roughly $\pm0.3$ absolute error on *every* LD correlation even at AFR's
   pooled $n=60{,}000$, and far worse below that. Fine-mapping resolution is LD
   precision — the cross-ancestry contrast of §4.3 works by distinguishing
   $r\approx0.98$ from $r\approx0.9$ — so that noise removes the mechanism the
   benchmark exists to measure, and a noised $G$ is no longer PSD, breaking the
   ELBO monotonicity SuSiEx converges on. If a formal guarantee is required,
   apply DP to the *pooled* $c_p$ only: it is $O(M)$ rather than $O(M^2)$, its
   sensitivity does not grow with $n$ while the signal does, and it leaves the
   LD panels untouched. Noising a subset of the uplink is not a guarantee at all
   — releasing $u_{p,k}$ in the clear leaves $\varepsilon=\infty$ by (c).

---

## Appendix A — Configured parameters (chr1 default)

| Symbol | Name in YAML | Value |
|---|---|---|
| `n_sites` | `sites` | 3 (ANL, Covenant, MBZUAI) |
| `n_per_site` | `sites.<site>.n` | 50,000 |
| `superpopulations` | `superpopulations` | EUR, AFR, AMR, EAS, CSA, MID |
| `master_seed` | `master_seed` | 20260601 |
| `chromosome` | `chromosome` | 1 |
| Window size | `locus_selection.window_size_bp` | 1,000,000 |
| Step size | `locus_selection.step_size_bp` | 500,000 |
| Loci selected | `locus_selection.n_loci` | 100 (configured target) |
| Tercile counts | `locus_selection.strata` | 33 / 33 / 34 (target); 25 / 26 / 28 realised after overlap dedup → 79 loci |
| Tag-SNP target | `locus_selection.ld_tag_snp_target` | 500 |
| MAF filter | `locus_selection.maf_filter` | 0.01 |
| LD prune r² | `locus_selection.ld_r2_prune_threshold` | 0.995 |
| ncsl grid | `architecture.ncsl` | {1, 2, 3} |
| h² grid | `architecture.h2` | {5e-4, 1e-3, 5e-3} |
| rg grid | `architecture.rg` | {0.5, 0.7, 1.0} |
| Factorial mode | `architecture.factorial_mode` | `extended` (15 arch) |
| Replicates | `architecture.replicates` | 10 |
| Ancestry-divergent causal | `architecture.ancestry_divergent_causal` | `enabled: true` (as of 2026-07-15), `n_private_per_pop: 1`, `min_per_stratum: 1` (§5.4) |
| h² tolerance | `phenotype.h2_tolerance_relative` | 0.20 |
| SuSiEx columns | (fixed) | one per superpopulation, pooled across sites — 6 columns, 150,000 individuals |
| Min pooled GWAS n | `fine_mapping.min_gwas_n` | 1000 (keeps all 6; smallest pooled cohort is EAS at 2,500) |
| LD-panel MAF filter | `--maf` (CLI, both FM stages) | 0.005 (SuSiEx's own default; applied to the *pooled* per-ancestry frequency) |
| Credible-set level / p filter | `--level`, `--pval-thresh` (CLI) | 0.95, 1e-5 |
| Federated blocks per locus | (fixed) | 12 = (site, ancestry) pairs: ANL 5, Covenant 3, MBZUAI 4 |
| Kinship threshold | `qc.kinship_threshold` | 0.05 |

## Appendix B — Repository layout

```
fedfm-simulation/
├── config/
│   └── simulation_config.yaml      # single source of truth (Pydantic-validated)
├── data/
│   ├── raw/hapnest/                # HAPNEST chr1 + population manifest
│   ├── processed/{anl,covenant,mbzuai}/  # per-site PLINK filesets + manifest
│   ├── loci/ld_matrices/           # per-locus tag SNPs + r² (locus_selection)
│   └── ground_truth/
│       ├── effect_sizes/           # per-arch per-replicate β
│       └── phenotypes/{site}/      # per-(site,locus,arch,rep) y
├── src/
│   ├── utils.py                    # config, seeding, .bed reader, PLINK shell-out
│   ├── sampling.py                 # §3
│   ├── locus_selection.py          # §4
│   ├── phenotype_sim.py            # §5
│   ├── qc.py                       # §6
│   ├── validation.py               # §10.5 (5-check harness)
│   ├── fine_mapping.py             # §11.3 centralized SuSiEx baseline
│   └── fed_fine_mapping.py         # §11.3 federated path (algo.md §7, §9)
├── tests/                          # 82 tests, all passing (§7.4)
├── scripts/
│   ├── _env.sh                     # Python + vendor/bin resolver
│   ├── install_plink.sh            # static PLINK 1.9 + 2.0 to vendor/bin/
│   ├── install_susiex.sh           # static SuSiEx to vendor/bin/
│   ├── download_hapnest.sh         # turnkey HAPNEST fetch
│   ├── run_<stage>.sh              # per-stage wrapper
│   ├── run_validation.py           # thin entry (loky needs an importable module)
│   ├── run_fine_mapping.py         # thin entry — centralized fine-mapping
│   ├── run_fed_fine_mapping.py     # thin entry — federated fine-mapping
│   ├── plot_fine_mapping.py        # 4 paper figures from fm_results.tsv
│   └── submit_polaris.pbs          # ALCF Polaris job (debug queue)
├── vendor/
│   ├── bin/                        # plink, plink2, SuSiEx (static)
│   └── SuSiEx/                     # upstream source, cited by file:line in algo.md
├── reports/                        # QC HTML, validation, fine_mapping/, fed_fine_mapping/
├── logs/                           # per-run structured logs
└── paper/
    ├── design.md                   # this document
    ├── algo.md                     # SuSiEx + the federated derivation
    └── progress.md                 # running log / TODO
```

## Appendix C — References

* Sophia Wharrie *et al.*, "HAPNEST: efficient, large-scale generation
  and evaluation of synthetic datasets for genotypes and phenotypes."
  *Bioinformatics* 39 (2023). EBI BioStudies accession S-BSST936.
* Kai Yuan *et al.*, "Fine-mapping across diverse ancestries drives the
  discovery of putative causal variants underlying human complex traits and
  diseases." *Nature Genetics* 56 (2024). SuSiEx reference implementation.
* Christopher C. Chang *et al.*, "Second-generation PLINK: rising to the
  challenge of larger and richer datasets." *GigaScience* 4 (2015).
* Argonne Leadership Computing Facility, "Polaris user guide,"
  https://docs.alcf.anl.gov/polaris/.
