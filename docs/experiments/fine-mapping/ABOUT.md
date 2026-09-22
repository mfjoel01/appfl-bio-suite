> Current protocol: [September scientific corrections and rerun](SCIENTIFIC_RERUN.md). Its validation gates have passed and the run is complete; the findings and limitations are in that document's [outcome section](SCIENTIFIC_RERUN.md#outcome--17-september-2026). Numerical results below that predate the rerun are historical.

# Federated cross-ancestry fine-mapping

**Statistical fine-mapping across institutions without pooling anybody's genome.**

A GWAS says *this locus matters*. Fine-mapping asks the harder question: *which variant
in it is causal?* The answer is a **credible set** — the smallest set of variants that
probably contains the causal one — and it depends on the correlation structure between
variants, which is a property of the cohort rather than of any single variant.

That dependence is what makes fine-mapping harder to federate than association testing.
A meta-analysis needs one number per variant per site. A credible set needs the whole
linkage-disequilibrium matrix, and LD is exactly the thing that does not decompose into
per-variant summaries.

The estimator is **SuSiEx**, a C++ command-line tool the coordinator shells out to.

## The claim is exactness, not approximation

The [GWAS experiment](../gwas/ABOUT.md) combines sites by inverse-variance meta-analysis,
which *approximates* what a pooled analysis would have found. That is the standard and it
is defensible.

This experiment does not approximate anything. By **Corollary 1** of the derivation in
[reference/algo.md](reference/algo.md), the federated fit **equals** the centralized fit —
same credible sets, same posterior inclusion probabilities, up to floating-point summation
order — provided three harmonization conditions hold. A discrepancy is a bug in the
federated path, never a cost of federating, and the test suite is written to that
standard: `tests/test_fine_mapping_federated.py::test_federated_equals_centralized` runs
both paths over one package and asserts they agree.

### Why it can be exact

Everything SuSiEx needs is a function of second moments, and a matrix product over rows
decomposes over **any** partition of those rows (Theorem 1). So each site computes its own
block and the coordinator adds them:

```
site k emits, per (site, ancestry, locus), on the raw 0/1/2 dosage scale:

    G = X'X   (M x M)      u = 1'X   (M)      n
    c = X'y   (M)          q = 1'y            w = y'y

coordinator:   X'X = Σ_k G_k     1'X = Σ_k u_k     n = Σ_k n_k     ...
```

Nothing else crosses the site boundary. No genotype row, no phenotype value, and — the
part that matters — **no locally standardized LD matrix**.

### The one real trap: sum first, standardize second

Standardizing at the sites and combining the results looks equivalent and is not. Each
site would centre against its own column means, deleting the between-site allele-frequency
variation before the coordinator could see it.

[algo.md §9 / A.10](reference/algo.md) works the smallest example: six people, two SNPs,
split three and three. Pooling raw moments and standardizing once gives a correlation of
**0.816**. Standardizing at each site first gives **0.866** — from either site, and
therefore from any weighted average of them. No reweighting recovers the right answer,
because the information was destroyed before the averaging.

`tests/test_fine_mapping_federated.py::test_standardizing_at_the_sites_is_biased` pins
both numbers. It is the single most important test in the experiment.

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
| Data | Synthetic — see [DATA.md](DATA.md) |

One exchange per site. `num_global_epochs: 1` is correct and not a placeholder — the site
moments do not change, so a second round recomputes them at full cost.

### Why one column per ancestry rather than one per site

SuSiEx takes a list of population columns and fine-maps them jointly. The obvious mapping
is one column per site, and it is wrong here: effect sizes are indexed by
**superpopulation**, not by site, so an ancestry-mixed site column is not one homogeneous
effect and violates the assumption SuSiEx's model makes about a column.

Pooling by ancestry across sites gives six columns covering every individual exactly once,
each a single homogeneous effect.

### The coordinator writes SuSiEx's LD panel by hand

The centralized path gets its LD reference panel from `plink --r square bin4` on genotypes
it holds. The coordinator holds none, so it writes the panel files itself, in the byte
layout the SuSiEx C++ reader demands — a float32 matrix that must be exactly `M² × 4`
bytes against the `_ref.bim` line count, and a `.frq` whose alleles must repeat each bim
line's exactly. SuSiEx aborts on a mismatch rather than degrading.

### Allele harmonization is not free

Theorem 1's first condition wants one ordered, allele-harmonized variant list. PLINK sets
A1 to each cohort's own minor allele, so a site can code a variant the other way round —
contributing `2 − x` where the others contribute `x`, which sums in silently and negates
every off-diagonal that variant touches. There is no error and no warning.

Each site therefore recodes against an agreed reference variant list before forming a
single moment, and the coordinator re-checks rather than trusting the report. This is why
a bundle carries `reference_variants.tsv` and why the run log reports a per-site flip
count.

### Governance and provenance are part of the design

Four GA4GH standards, each closing a hole the design otherwise had:

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
string comparison rather than an assumption.

All of it is opt-in per experiment. **→ [../../coordinator/ga4gh.md](../../coordinator/ga4gh.md)**

### What is deliberately not addressed

**Privacy.** The genotypes are synthetic, so aggregates ship in the clear. There is no
secure aggregation and no differential privacy. Note what DUO does and does not do: it
decides whether a study may run against a site's data; it does not make the uplink
private. A site whose terms permit a study is still shipping an `X'X`, and an `X'X` over a
small cohort is not a harmless object. [reference/design.md §11.5](reference/design.md)
says what a real deployment would involve; none of it is implemented here.

## Does federating actually help?

Correctness and value are different questions. The exactness result says a federated fit
*equals* the pooled one; it says nothing about whether pooling was worth the trouble. That
needs an arm in which a site works alone, so the experiment carries a set of
**participation arms** — the same loci, the same architectures, the same ground truth,
differing only in who took part:

| arm | who | ancestry columns |
| --- | --- | --- |
| one site alone | that site's cohort only | 3–5, depending on the site |
| all three, *n* matched | every site, down-sampled to one site's worth of people | 5 |
| all three | the whole federation | 6 |
| borrowed LD | one site's summary statistics against another's LD panel | 3 |

The matched arm holds analyzed sample size fixed. Remaining differences reflect
composition and participation, including allele frequencies, effects and phenotype noise.
They cannot be attributed uniquely to LD diversity.

The fourth is the objection rather than the design: the cheap alternative to shipping
O(M²) is "send me your summary statistics and I'll use my own LD panel." That is the
fine-mapping analogue of applying a model fitted on one cohort to another, and it is the
LD-mismatch failure mode. Whether it produces well-calibrated credible sets is the test of
whether the uplink is justified.

Site-alone arms need no code — `plan_columns`, `materialize_column_keeps` and
`pooled_phenotype` are all scoped by `cfg.sites`, so a config listing one site *is* a
site-alone fit. See `experiments/fine_mapping/arms.py`.

### What the arms found

From the completed corrected run; full numbers, intervals and limitations in
[SCIENTIFIC_RERUN.md](SCIENTIFIC_RERUN.md#outcome--17-september-2026).

| arm | analyzed N | power | mean set size | coverage |
| --- | ---: | ---: | ---: | ---: |
| all three | 150,000 | **0.943** | **6.8** | 0.963 |
| Covenant alone | 50,000 | 0.777 | 9.6 | 0.956 |
| borrowed LD | 50,000 | 0.777 | 9.6 | 0.955 |
| MBZUAI alone | 50,000 | 0.626 | 9.5 | 0.957 |
| ANL alone | 50,000 | 0.622 | 12.7 | 0.947 |
| all three, *n* matched | 50,000 | 0.489–0.497 | 9.2 | 0.946 |

Federating wins on every axis at once — more power, more sets per locus, tighter sets,
higher purity, and the lowest false-discovery rate at both PIP thresholds. So the uplink
is justified, and the O(M²) cost buys something real.

**But the *power* gain is sample size rather than diversity, and the matched arm shows
it.** Hold analyzed N at 50,000 and the three-site composition becomes the *worst* arm in
the run, below every single site, reproducibly across three independent draws. The reason
is in the design above: SuSiEx fits one effect per ancestry column, so power tracks the
**largest column** rather than the total. Covenant alone puts 47,500 of its 50,000 people
into AFR; the matched federation's largest column is 20,339. Splitting a fixed cohort
across ancestries costs 13.3 pp (12.3–14.3) in the controlled comparison that holds both N
and column count fixed. Adding cohorts buys it back several times over.

**Diversity earns its keep at the other end, though.** Among the 4,987 instances where both
arms returned a set — comparing like with like, rather than letting the weaker arm look
sharp because it only reports its easiest finds — the diverse arm's credible sets are
**2.16 variants smaller (95% CI 1.75–2.58)**, with higher purity and higher causal PIP, at
an essentially equal capture rate. That is the textbook cross-ancestry benefit: different
LD across ancestries breaks ties between correlated variants.

So at fixed N the trade is fewer signals found, each pinned down better. The full
federation sidesteps it by not being at fixed N — 150,000 people buy both. Calling either
matched-N bar simply "diversity" misstates the finding in one direction or the other.

The borrowed-LD objection survives its test on this design: Covenant's statistics against
ANL's AFR panel land within 0.0004 of Covenant's own power and 0.0007 of its coverage.
That bounds the cost of *external-panel sampling* between same-ancestry panels drawn from
one synthetic pool. It is not a mismatched reference, which is the failure mode the
objection is really about.

The one place credible sets stop being trustworthy is weak signal. At h² = 0.0005 the
full federation still covers at 0.949, but every 50,000-person arm under-covers — down to
**0.799** for the matched control, whose sets in that cell survive a 9.9% discovery rate
and are selected on having got lucky. A 95% credible set from a single site at the weakest
simulated signal is not a 95% credible set.

## Figures

The coordinator draws its diagnostics in-process at the end of a run. A fuller set —
exploratory, paper-form, and federation-cost figures — comes from a batch renderer; see
[the figures package](../../../src/appfl_bio_suite/experiments/fine_mapping/figures/README.md).

## Where the code lives

```
src/appfl_bio_suite/experiments/fine_mapping/
├── dataset.py     SHIPPED to workers — self-contained, no suite imports
├── trainer.py     SHIPPED to workers — self-contained, no suite imports
├── aggregator.py  coordinator-side   — may import freely
├── plotting.py    coordinator-side   — the aggregator's figure entry point
├── figures/       coordinator-side   — the full figure set
├── fedfm/         coordinator-side   — the maintained science fork; see UPSTREAM.json
└── simulation/    coordinator-side   — never reaches a partner at all
```

`fedfm/` is a maintained fork with upstream hashes and changed-module provenance in
`fedfm/UPSTREAM.json`. The shipped site stage is self-contained; its numerical results
are compared against the coordinator implementation at zero tolerance by
`tests/test_fine_mapping_shipped_parity.py`.

## Running it

**→ [RUNBOOK.md](RUNBOOK.md)**
