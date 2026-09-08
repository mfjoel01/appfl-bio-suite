# Federated cross-ancestry fine-mapping — design and record of runs

## What this experiment is

**Statistical fine-mapping across institutions without pooling anybody's genome.** Given a
locus that a GWAS has implicated, fine-mapping asks the harder question: *which* variant in
it is causal? The answer is a **credible set** — the smallest set of variants that probably
contains the causal one — and it depends on the correlation structure between variants,
which is a property of the cohort rather than of any single variant.

That dependence is what makes fine-mapping harder to federate than association testing. A
meta-analysis needs one number per variant per site. A credible set needs the whole
linkage-disequilibrium matrix, and LD is exactly the thing that does not decompose into
per-variant summaries.

Cross-ancestry, because it is the lever that makes fine-mapping work. Two ancestries that
disagree about which variants are correlated will disagree about which variant in a
credible set could be the causal one — and the intersection of their disagreements is much
smaller than either set alone. The federation here is three sites with sharply different
ancestry compositions (EUR-dominant, AFR-dominant, MID-dominant), which is the realistic
version of that lever and the reason a federation is worth building rather than simulating.

The estimator is **SuSiEx**, a C++ command-line tool the coordinator shells out to.

## The claim is exactness, not approximation

This is the distinction that shapes everything below.

The [GWAS experiment](../gwas/ABOUT.md) combines sites by inverse-variance meta-analysis,
which *approximates* what a pooled analysis would have found, under an assumption that all
sites estimate the same effect. That is the standard and it is defensible.

This experiment does not approximate anything. By **Corollary 1** of the derivation in
[reference/algo.md](reference/algo.md), the federated fit **equals** the centralized fit —
same credible sets, same posterior inclusion probabilities, up to floating-point summation
order — provided three harmonization conditions hold. A discrepancy is a bug in the
federated path, never a cost of federating, and the test suite is written to that standard:
`tests/test_fine_mapping_federated.py::test_federated_equals_centralized` runs both paths
over one package and asserts they agree.

### Why it can be exact

Everything SuSiEx needs is a function of second moments, and a matrix product over rows
decomposes over **any** partition of those rows (Theorem 1). So each site can compute its
own block and the coordinator can add them:

```
site k emits, per (site, ancestry, locus), on the raw 0/1/2 dosage scale:

    G = X'X   (M x M)      u = 1'X   (M)      n
    c = X'y   (M)          q = 1'y            w = y'y

coordinator:   X'X = Σ_k G_k     1'X = Σ_k u_k     n = Σ_k n_k     ...
```

Nothing else crosses the site boundary. No genotype row, no phenotype value, and — the
part that matters — **no locally standardized LD matrix**.

### The one real trap: sum first, standardize second

Standardizing at the sites and combining the results looks equivalent and is not. Each site
would centre against its own column means, deleting the between-site allele-frequency
variation before the coordinator could see it.

[algo.md §9 / A.10](reference/algo.md) works the smallest example: six people, two SNPs,
split three and three, where one SNP's frequency differs between the sites. Pooling raw
moments and standardizing once gives a correlation of **0.816**. Standardizing at each site
first gives **0.866** — from either site, and therefore from any weighted average of them.
No reweighting recovers the right answer, because the information was destroyed before the
averaging.

`tests/test_fine_mapping_federated.py::test_standardizing_at_the_sites_is_biased` pins both
numbers. It is the single most important test in the experiment.

## Design

| | |
| --- | --- |
| Federation pattern | Single round, raw second moments |
| Server | `FineMappingAggregator` (coordinator-side, not shipped) |
| Client | `SiteFineMappingTrainer` (shipped to workers, self-contained) |
| Estimator | SuSiEx, one population column per superpopulation |
| Columns | 6 — EUR, AFR, AMR, EAS, CSA, MID — each pooling that ancestry across every site |
| Credible sets | 95% coverage; marginal p-value filter 1×10⁻⁵ |
| Uplink | O(M²) per (locus, ancestry); ~32 MB per block at a realistic locus |
| Data | Synthetic; simulated and distributed by the coordinator — see [DATA.md](DATA.md) |

One exchange per site. `num_global_epochs: 1` is correct and not a placeholder — the site
moments do not change, so a second round recomputes them at full cost.

### Why one column per ancestry rather than one per site

SuSiEx takes a list of population columns and fine-maps them jointly. The obvious mapping
is one column per site, and the first complete run used it. It is wrong for this design:
the simulator indexes effect sizes by **superpopulation**, not by site, so an
ancestry-mixed site column is not one homogeneous effect and violates the assumption
SuSiEx's model makes about a column.

Pooling by ancestry across sites gives six columns covering all 150,000 individuals exactly
once, each a single homogeneous effect. The earlier site-column results are superseded and
**not** comparable column-for-column, so they are kept in the coordinator's `local/`
directory rather than published here.

### The coordinator writes SuSiEx's LD panel by hand

The centralized path gets its LD reference panel from `plink --r square bin4` on genotypes
it holds. The coordinator holds none, so it writes the three panel files itself, in the
byte layout the SuSiEx C++ reader demands — a float32 matrix that must be exactly `M² × 4`
bytes against the `_ref.bim` line count, and a `.frq` whose alleles must repeat each bim
line's exactly. SuSiEx aborts on a mismatch rather than degrading, which is why "SuSiEx ran
at all" in the loopback run is a real signal.

### Allele harmonization is not free

Theorem 1 condition (i) wants one ordered, allele-harmonized variant list. The per-site
filesets were cut with `plink --keep --make-bed` and *without* `--keep-allele-order`, so
PLINK 1.9 set A1 to each site's own minor allele: about 6% of chr1 variants are coded
oppositely at two of the three sites. A flipped site contributes `2 − x` where the others
contribute `x`, which sums in silently and negates every off-diagonal that variant touches.

Each site therefore recodes against an agreed reference variant list before forming a
single moment, and the coordinator re-checks rather than trusting the report. This is why
the bundle carries `reference_variants.tsv` and why the run log reports a per-site flip
count.

### Governance and provenance are part of the design, not documentation around it

The experiment is assembled out of four GA4GH standards, each closing a hole the design
otherwise had:

| | | what it answers |
| --- | --- | --- |
| **DUO** | Data Use Ontology | may this site's data be used for this study? |
| **DRS** | Data Repository Service | which bytes did this site compute over? |
| **TRS** | Tool Registry Service | which tool version computed them? |
| **TES** | Task Execution Service | dispatch that tool, to a site that speaks TES |

Two of them change what a run can silently get wrong. A site's consent code lives in its
own bundle and is enforced by its own worker before a genotype is opened — so a study a
site's terms do not permit is refused where the data is, not where the coordinator is. And
each bundle is content-addressed, so *this site ran on the bundle I cut for it* is a
string comparison rather than an assumption; a site running last month's bundle otherwise
produces well-formed aggregates over the wrong individuals, pools without complaint, and
is wrong in no visible way.

The third makes results attributable: the pin covers the source of the two modules APPFL
ships to workers, so "which code produced this credible set" is a lookup. The fourth is
optional and adds a second transport.

All of it is opt-in per experiment, and a federation declaring none of it runs exactly as
it did before. **→ [../../coordinator/ga4gh.md](../../coordinator/ga4gh.md)**

### What is deliberately not addressed

**Privacy.** The genotypes are synthetic, so aggregates ship in the clear. There is no
secure aggregation and no differential privacy. Note what DUO does and does not do here:
it decides whether a study may run against a site's data, and it does not make the uplink
private. A site whose terms permit a study is still shipping an `X'X`. A real deployment would need both, and
[reference/design.md §11.5](reference/design.md) says what that would involve; none of it
is implemented here. An `X'X` over a small cohort is not a harmless object.

## Record of runs

### Centralized baseline at full scale — 2026-08-02

The reference fit: six ancestry columns, each pooling that ancestry across all three sites.
This is what SuSiEx returns with all the data in one place.

| | |
| --- | --- |
| Loci | 79 (of 100 requested; the rest dropped by the non-overlap constraint) |
| Architectures | 15 (`extended` grid) × 10 replicates |
| Instances | 11,850 |
| Failures | 0 |
| Overall power (any causal captured) | **0.881** |
| Median instance runtime | 17.4 s |

Power rises with heritability and with cross-ancestry genetic correlation, and falls as the
number of causal variants per locus rises — all as expected:

| ncsl | h² | rg | power | mean credible sets | median best CS size | mean causal PIP |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.0005 | 1.0 | 0.819 | 0.86 | 3 | 0.567 |
| 1 | 0.001 | 1.0 | 0.935 | 0.96 | 2 | 0.683 |
| 1 | 0.005 | 1.0 | 0.930 | 0.98 | **1** | 0.839 |
| 3 | 0.0005 | 1.0 | 0.657 | 0.75 | 4 | 0.481 |
| 3 | 0.005 | 1.0 | 0.939 | 2.09 | **1** | 0.909 |

Full tables: [results/fine_mapping/](results/fine_mapping/), figures in
[results/fine_mapping/figures/](results/fine_mapping/figures/).

### Data-package validation — 5/5 PASS

The five-check harness over the realised package: phenotype re-derivation (35,550/35,550
instances match to a worst relative error of 1.6×10⁻⁷), heritability calibration,
cross-ancestry effect-size correlation (recovering rg = 0.5/0.7/1.0 to within 0.01),
structural completeness (14 checks), and LD matrix sanity over 1,374 files.
[results/validation/SUMMARY.txt](results/validation/SUMMARY.txt).

### Federated path — proven exact, not yet run at full scale

The federated implementation is complete and its equality with the centralized path is
demonstrated, but **the production sweep has not been run**. What exists:

- **Unit and acceptance tests**, 28 of them, including the end-to-end
  federated-equals-centralized assertion over a miniature three-site package with real
  SuSiEx. The summary-statistics half is bit-exact; the LD half agrees with PLINK's own
  panel to within 2 ulp in float32, which is the precision PLINK itself writes.
- **A production spot-check** (locus L0000, `ncsl1_h2-0.0005_rg1`, rep 0): the
  coordinator's closed-form summary statistics were byte-identical to `plink2 --glm` on all
  9,848 variants, the credible set was the same, and the maximum PIP difference was
  6.9×10⁻¹⁶.
- **The APPFL path**, which reproduces the standalone federated driver's result table
  exactly — `tests/test_fine_mapping_appfl_parity.py`, zero tolerance on every statistic.

What does not exist is a 79-locus federated sweep to put beside the centralized table
above. [RUNBOOK.md](RUNBOOK.md) has the command; it is a scheduler allocation away.

### Known defect in the ground truth

Found 2026-08-20, and it affects the numbers above. The allele-order problem described
earlier reached the *phenotypes*: `phenotype_sim` multiplies β by that site's dosage with no
allele handling, so a site coding a causal variant the other way round gave its people
`−β`. **1,582 of 22,305 causal variants (7.1%) are coded inconsistently across sites,
affecting 1,595 of 11,850 manifest rows (13.5%).**

The consequence is attenuated pooled signal at those variants — a median 92% of magnitude
retained for EUR, 89% for CSA, and **58% for AFR**, which is worst because Covenant holds
47,500 of the 60,000 AFR individuals and has the most flips. Single-site ancestries
(AMR, EAS, MID) are untouched: a uniform flip is just β → −β and changes no PIP.

Two things follow. The reported power **understates** the design's true power on the
affected instances. And it does **not** threaten the exactness result, because both paths
read the same phenotypes and are affected identically.

The fix is to re-cut the site filesets with `--keep-allele-order` and re-run the phenotype
stage; locus selection is unaffected and resumes from its cache. Until then, treat the
power figures as a lower bound.

### Migration note

This experiment was ported from a standalone repository rather than written here. Its eight
science modules are vendored **byte-for-byte** under
`src/appfl_bio_suite/experiments/fine_mapping/fedfm/`, and
`tests/test_fine_mapping_configs.py` asserts that with `cmp` against the upstream checkout
when one is available. Its 82 upstream tests were ported with nothing changed but their
import paths.

What this repository added is the APPFL layer — a shipped trainer, an aggregator that
delegates every statistic to the vendored code, per-site bundling, and a synthetic cohort
provider so the pipeline runs without a 135 GB download. The site stage exists twice as a
result, because shipped source may not import the suite; the two implementations are
compared at zero tolerance by `tests/test_fine_mapping_shipped_parity.py`.

### Template for a run of record

```
### YYYY-MM-DD — <what this run was>

| | |
| --- | --- |
| Path | centralized | federated (standalone) | federated (APPFL) |
| Sites | 3 |
| Loci | |
| Instances | |
| Scenario | |
| Sim manifest | <sha256 of run_manifest.json> |
| Suite commit | <sha> |

Power (any causal captured): ...
Median best CS size: ...

Reproduce:
    appfl-bio-suite simulate fine-mapping --scenario <name> --out <dir>
    appfl-bio-suite run fine-mapping
```

Record the **simulation manifest hash** alongside the commit. The commit says which code
ran; the manifest says which data it ran on, and without both the run is not reproducible.

---

## Reference material

Carried over from the standalone repository, under [reference/](reference/):

- **[algo.md](reference/algo.md)** — SuSiEx, and the federated derivation with its proofs.
  Theorem 1, Corollary 1 and the §9 standardization order are what the aggregator
  implements; read this before changing anything in `aggregator.py`.
- **[design.md](reference/design.md)** — the simulator's design and its realised numbers.
- **[pbs/](reference/pbs/)** — the PBS scripts the published run was actually submitted
  with, kept verbatim as a record. The portable equivalents are in
  `scripts/fine-mapping/`.

The standalone repository's running log and TODO list are not carried over. They were a
working record of one operator on one allocation — dated entries, queue notes, and the
paths of a retired tree — so they live in the coordinator's `local/` directory under the
rule in [local/README.md](../../../local/README.md) rather than in the published docs.
What that log established about the method is in `algo.md` and `design.md`; what it
established about running the thing is in [RUNBOOK.md](RUNBOOK.md).

See [DATA.md](DATA.md) for the data pipeline and its limits, and [RUNBOOK.md](RUNBOOK.md)
for how to run it.
