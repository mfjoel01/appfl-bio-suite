> Current protocol: [September scientific corrections and rerun](SCIENTIFIC_RERUN.md). Its validation gates have passed and the run is complete; the findings and limitations are in that document's [outcome section](SCIENTIFIC_RERUN.md#outcome--17-september-2026). Numerical results below that predate the rerun are historical.

# Fine-mapping — the data

Synthetic throughout. No real genotypes are used, which is why aggregates ship in the
clear and why the ground truth is known exactly — the point of the package is that every
credible set can be scored against the causal variants that actually generated the
phenotype.

## Two providers

| Provider | What it is | When |
| --- | --- | --- |
| `hapnest` | the real pool — HAPNEST synthetic genotypes, six superpopulations, ~1M individuals | the published design; needs a staged download (~135 GB for chr1) |
| `synthetic_hapnest` | a generated substitute with the same shape | CI, install checks, and anywhere the download is impractical |

Both produce the same directory layout, so everything downstream is provider-agnostic.

## The pipeline

**Step 1 — per-site cohorts.** Sites are stratified out of the pool to declared ancestry
compositions, with **no individual reused across sites**. Disjointness is asserted before
any extraction and again after bundling: an individual appearing twice would be counted
twice by the pooled moments while `n` reports them once, biasing every standardized
quantity downstream with nothing to indicate it.

The three site identities are ancestry-contrasting by design — one EUR-dominant, one
AFR-dominant, one MID-dominant. That contrast is the reason a federation is worth building
rather than simulating.

> The site ids are not free: the vendored sampler enumerates a fixed list and silently
> skips any other name, so the scenario schema refuses unknown ids rather than quietly
> producing a smaller cohort.

**Step 2 — locus selection.** Candidate windows are tiled across the chromosome, scored
for cross-site LD divergence (the Frobenius norm of pairwise r² differences), and
stratified into low/medium/high bands. Overlapping selections are dropped greedily, so
asking for 100 loci yields fewer.

> **Read the score carefully.** It is bounded by the LD magnitudes it differences — a
> window where every site has r² ≈ 0 cannot score high — so it correlates strongly with
> how much LD is in the window and does *not* cleanly isolate "the sites disagree about
> LD". Treat it as an LD-structure covariate, not as a diversity variable.

**Step 3 — causal architectures and phenotypes.** For each (locus × architecture ×
replicate): pick causal variants, draw per-superpopulation effect sizes from MVN(0, Σ_rg),
and generate a Gaussian phenotype at each site. The grid is a star design over three
parameters — `ncsl` (causal variants per locus), `h2` (per-locus heritability), and `rg`
(cross-ancestry genetic correlation).

Effect sizes are indexed by **superpopulation, not by site**. An individual's genetic
value uses the β for their own ancestry, wherever they happen to be enrolled. That is what
makes an ancestry column a single homogeneous effect and what makes pooling by ancestry
the right column definition.

**Step 4 — per-site bundles.** The pipeline's own layout assumes one machine can see
everything, which is exactly the assumption a federation removes. This step rearranges it
into one self-contained directory per site.

**Step 5 — DRS objects and data use terms.** Each bundle is content-addressed and carries
its site's consent profile, so a site can refuse a study its terms do not permit, at the
site, before a genotype is opened.

### What is deliberately not in a bundle

No ground truth. A site holds its own genotypes and phenotypes and nothing else — it
cannot score its own credible sets, and it cannot see another site's individuals. The
causal manifest stays with the coordinator.

## Scenarios

```bash
appfl-bio-suite simulate fine-mapping --list-scenarios
```

A scenario carries the upstream pipeline configuration **verbatim** under a `pipeline:`
key; the suite supplies only `paths:` and `cohort:`. So a scenario is readable side by
side with the standalone repository's config, and diffing them answers "did the migration
change any parameter?" mechanically — `tests/test_fine_mapping_configs.py` does exactly
that.

`ci-tiny` is deliberately unrealistic in two ways, both called out in the file, because a
fixture that quietly stops testing the hard case is worse than no fixture: its
heritability is roughly a hundred times a real per-locus value (so a smoke test can tell
"the pipeline works" from "the fit is broken"), and its missingness is above the pool's
(so the complete-case path is actually exercised).

## Per-site QC

Relatedness, per-superpopulation MAF concordance, and LD decay, as one HTML report per
site. It deliberately skips PCA: the ancestry labels are known by construction here, so a
PCA would be confirming the simulator rather than checking the data.

## Validation

A five-check harness over the realised package: phenotype re-derivation from the recorded
β and genotypes, heritability calibration against target, cross-ancestry effect-size
correlation against target `rg`, structural completeness, and LD matrix sanity.

```bash
scripts/fine-mapping/run_stage.py validation --config <pipeline_config.yaml>
```

Re-derivation is the load-bearing one: it recomputes each phenotype from scratch and
compares, so a silent change in the simulator shows up as a mismatch rather than as a
plausible-looking number.
