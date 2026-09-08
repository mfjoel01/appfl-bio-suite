# Federated fine-mapping — how the data comes into existence

Every experiment in this suite has this document, and its content differs sharply between
them. For FLamby it is "download this public dataset". For this one it is a five-stage
simulation that produces a three-site cohort, a set of loci chosen for how much their LD
structure diverges between sites, phenotypes under a grid of causal architectures, and a
ground-truth answer key that **stays with the coordinator**.

```
appfl-bio-suite simulate fine-mapping --scenario ci-tiny --out /tmp/fm-ci
```

## The two providers, and which one produced a given result

The pipeline subsamples **HAPNEST**, a pre-generated synthetic multi-ancestry cohort
published on EBI BioStudies as accession `S-BSST936`. Chromosome 1 alone is about 135 GB;
the full 22-chromosome set is 1.5 TB.

That is a reasonable price for the production run and an unreasonable one for finding out
whether your install works. So there are two cohort providers, and the run manifest records
which was used, so no output is ever ambiguous about its own provenance.

### `hapnest` — the real pool

What the published run used, and the only provider whose output means anything about
genetics. Stage it once:

```bash
scripts/fine-mapping/download_hapnest.sh              # all 22 chromosomes (~1.5 TB)
CHROMS="1" scripts/fine-mapping/download_hapnest.sh   # chr1 only (~135 GB) — the default
export FEDFM_HAPNEST_DIR=/path/to/data/raw/hapnest
```

The script is resumable and idempotent — a partial file is detected by comparing local size
against the EBI `Content-Length` header — and it aborts with a clear error if the target
volume has insufficient space. It fetches the population labels, the per-chromosome PLINK
filesets, joins them into `population_manifest.tsv`, and symlinks the names the pipeline
expects.

The staged directory is named by `$FEDFM_HAPNEST_DIR` rather than written into the
scenario. A staged pool is a fact about a machine, not about a scenario, and a committed
scenario carrying one coordinator's path means every other coordinator edits a tracked file
and then carries that edit forever.

### `synthetic_hapnest` — a runnable substitute

Generates a small pool with HAPNEST's *shape*: the same file layout, the same
superpopulation labelling, the same variant-major PLINK1 encoding, a fraction of the size.
Deterministic from a seed, no download, seconds to run. This is what makes the whole chain
runnable from this repository alone, which is most of what makes it checkable.

**It is not HAPNEST.** It is a different pool from a different generative process, it does
not reproduce the published figures, and it must not be described as doing so. It
reproduces the *pipeline*, not the *numbers*.

Its model is documented in full in `simulation/cohort.py`, because a substitute whose own
provenance is unclear repeats the problem it exists to solve. In outline: allele
frequencies from a Beta spectrum, then per-superpopulation frequencies from the
Balding–Nichols model with `ld_divergence` playing the role of Fst; LD blocks whose
within-block correlation is itself perturbed per superpopulation, so ancestries differ in
*how correlated* a block is and not only in its frequencies; two thresholded haplotype
draws per individual; a little missingness so the complete-case path is exercised.

The per-ancestry LD divergence is the part that matters here. Fine-mapping is entirely a
statement about correlation structure, and a pool with one shared LD pattern would make the
low/medium/high divergence strata indistinguishable and the locus-selection stage vacuous.

## The pipeline

Five stages. `simulate` runs all of them in one process.

### Step 0 — the genotype pool

Stage HAPNEST, or generate the substitute. See above.

### Step 1 — per-site cohorts

Three sites, **no individual reused across sites**, stratified by superpopulation to the
composition each site declares:

| Site | n | Composition | Dominant |
| --- | --- | --- | --- |
| `anl` | 50,000 | EUR 30,000 · AFR 7,500 · AMR 7,500 · EAS 2,500 · CSA 2,500 | EUR |
| `covenant` | 50,000 | AFR 47,500 · EUR 1,500 · CSA 1,000 | AFR |
| `mbzuai` | 50,000 | MID 25,000 · CSA 15,000 · AFR 5,000 · EUR 5,000 | MID |

Disjointness is asserted **before** any PLINK extraction, and again after bundling. It is
the premise of the whole federation: an individual appearing twice would be counted twice
by the pooled moments while `n` reports them once, which biases every standardized quantity
downstream with nothing to indicate it.

> **The three site ids are not free.** The vendored sampler enumerates exactly `anl`,
> `covenant`, `mbzuai` in a fixed order and silently skips any other name, producing a
> cohort smaller than the scenario asked for. Rather than work around that — which would
> mean editing code whose value is being unedited — the scenario schema refuses any other
> site id with a message saying where the list lives. A fourth site is a change to make
> upstream and re-vendor.

### Step 2 — locus selection

Candidate windows are tiled across the chromosome, each scored for **cross-site LD
divergence** — the Frobenius norm of pairwise r² differences between sites — and the
selection is stratified into low, medium and high divergence bands. Overlapping selections
are dropped greedily by score, which is why asking for 100 loci yields 79.

Divergence is the independent variable the whole experiment turns on: it is what makes
cross-ancestry fine-mapping sharper than single-ancestry fine-mapping, so the loci are
chosen to span it rather than sampled at random.

### Step 3 — causal architectures and phenotypes

For each (locus × architecture × replicate) instance: pick causal variants, draw
per-superpopulation effect sizes from MVN(0, Σ_rg), and generate a Gaussian phenotype at
each site.

The architecture grid is a star design over three parameters — `ncsl` (causal variants per
locus), `h2` (per-locus heritability), and `rg` (cross-ancestry genetic correlation). The
`extended` mode is 15 points: the full ncsl × h2 grid at rg = 1.0, plus the ncsl × rg grid
at h2 = 0.001, deduplicated. Ten replicates each.

Effect sizes are indexed by **superpopulation, not by site**. An individual's genetic value
uses the β for their own ancestry, wherever they happen to be enrolled. That is what makes
an ancestry column a single homogeneous effect and what makes pooling by ancestry the right
column definition — see [ABOUT.md](ABOUT.md).

> See ABOUT.md's *Known defect in the ground truth*: in the realised package this invariant
> is violated for 13.5% of instances by an allele-order problem, and those instances carry
> site-dependent effect directions.

### Step 4 — per-site bundles

The pipeline's own layout assumes one machine can see everything, which is exactly the
assumption a federation removes. This step rearranges it into one self-contained directory
per site:

```
<out>/<site>/data/
├── site_genotypes.{bed,bim,fam}   that site's cohort, nobody else's
├── site_manifest.tsv              FID, IID, superpopulation
├── reference_variants.tsv         the agreed canonical allele coding
├── selected_loci.tsv              which windows this study fine-maps
├── phenotypes/<instance>.pheno    one per locus × architecture × replicate
└── DATA_USE.json                  this site's DUO terms (see step 5)
```

`reference_variants.tsv` is reference metadata — the same annotation table a public panel
or a consortium protocol distributes before any data is touched. It says nothing about any
individual, and without it the sites' Grams are not addable.

### What is deliberately NOT in a bundle

**`ground_truth/causal_manifest.tsv`** — the answer key. It stays at the coordinator. A site
holding it could score its own credible sets, and a benchmark whose answers travel with its
inputs is not measuring what it claims to. The aggregator reads it from the coordinator's
copy via `aggregator_kwargs.causal_manifest`; a real federation leaves that unset and gets
credible sets it cannot score, which is correct.

Also not in a bundle: any other site's anything.

### Step 5 — DRS objects and data use terms

Every bundle is registered as a GA4GH DRS object and written a `DATA_USE.json`. Both land
in `<out>/drs_registry.json` and `<out>/<site>/data/` respectively.

**The consent code travels with the data.** It is written *into* the bundle, not kept at
the coordinator, because the copy that decides is the copy the site holds — its worker
reads it and refuses a study the terms do not permit before opening a genotype file. The
coordinator's copy is for the offline preview and has no authority.

**The simulated terms are illustrative, and every generated profile says so.** The
individuals do not exist and no consent was given, so any terms here are invented. They
are declared anyway because the alternative is a federation whose governance machinery is
exercised for the first time against a partner's real dataset. The three shipped sites
declare *different* terms —

| site | permission | modifiers |
| --- | --- | --- |
| `anl` | `DUO:0000042` general research use | — |
| `covenant` | `DUO:0000006` health/medical/biomedical research | `DUO:0000018` not for profit, non-commercial only |
| `mbzuai` | `DUO:0000042` general research use | `DUO:0000019` publication required |

— so the loopback run matches a real consent code rather than a trivially permissive one.
Adding population-origins research to the run's declared purposes makes `covenant` refuse,
which is the behaviour worth having proven before a partner's terms are in play. Change
them in the scenario's `data_use:` block.

Ids are content-addressed: a blob's id is its sha-256, a bundle's is a Merkle hash over
its members. The run manifest already checksums everything, so what DRS adds is a *name*
for those checksums that travels — into a client config, a task document, a results table,
and a site that has never seen the manifest.

`DATA_USE.json` is written *after* the object is registered, deliberately: the profile
carries the object's `drs_uri`, and a file cannot be inside its own checksum. The site-side
verifier excludes exactly that one filename; anything else in a bundle the object does not
list is still a hard error.

**→ [../../coordinator/ga4gh.md](../../coordinator/ga4gh.md)**

## Per-site QC, and the PCA it deliberately skips

`scripts/fine-mapping/run_stage.py qc` writes one HTML report per site — KING relatedness
against `qc.kinship_threshold`, per-superpopulation MAF summaries, LD decay curves, and the
cross-site divergence summary. The published reports are
[results/qc_anl.html](results/qc_anl.html) and siblings.

What the report does **not** contain is an ancestry PCA, and that is a scope decision, not
an oversight. QC operates on chr1 only because the rest of the simulation does, and a
chr1-only PCA is not a valid ancestry sanity check; `fedfm/qc.py` skips it with a logged
warning rather than plotting something misleading. The config already carries the settings
for the real thing (`qc.pca_n_components: 20`, `qc.pca_reference: "1kg_phase3"`), and
enabling it is a known follow-on with three steps:

1. Cut **full-genome** per-site binaries using the same `<site>_ids.txt` keep files the
   sampling stage already writes — the cohort assignment does not change, only its extent.
2. Stage the 1KG Phase 3 PLINK reference binaries.
3. Project each site's individuals onto the 1KG PC space (e.g. `plink2 --score`).

## Scenarios

```bash
appfl-bio-suite simulate fine-mapping --list-scenarios
```

| Scenario | Provider | Sites | Enrolled | What it is for |
| --- | --- | --- | --- | --- |
| `ci-tiny` | `synthetic_hapnest` | 3 | 1,500 | CI and install validation. Seconds. Produces no meaningful genetics. |
| `three-site-hapnest` | `hapnest` | 3 | 150,000 | The published design. Needs the staged pool; hours on a many-core node. |

A scenario carries the upstream pipeline configuration **verbatim** under a `pipeline:`
key, and the suite supplies only `paths:` (from `--out`) and `cohort:` (where the genotypes
come from). So a scenario is readable side by side with the standalone repository's
`config/simulation_config.yaml`, and diffing them answers "did the migration change any
parameter?" mechanically. `tests/test_fine_mapping_configs.py` does exactly that.

### `ci-tiny` is deliberately unrealistic in two ways

Both are called out in the file, because a fixture that quietly stops testing the hard case
is worse than no fixture:

- **h² = 0.3**, roughly a hundred times a real per-locus heritability. At 1,500 individuals
  spread over six ancestry columns, a realistic effect produces no detectable signal,
  SuSiEx correctly returns no credible sets, and the smoke test cannot distinguish "the
  pipeline works" from "the fit is broken".
- **missing rate 5×10⁻⁶**, above HAPNEST's 2×10⁻⁷. Scaled naively it would drop nothing and
  leave the complete-case path untested; the obvious overcorrection (5×10⁻⁴) drops a
  quarter of every window and makes the exactness warning fire on every locus, which trains
  people to ignore it.

### Defining your own

Copy a shipped scenario, edit, and pass the path:

```bash
cp src/appfl_bio_suite/experiments/fine_mapping/configs/simulation/ci-tiny.yaml my.yaml
appfl-bio-suite simulate fine-mapping --scenario ./my.yaml --out /scratch/mine
```

Validation refuses, with the reason, a scenario that names an unknown site, enrols more
individuals than the pool holds, demands more of a superpopulation than the pool would
supply, or names an unknown provider. Each of those otherwise fails much later and much
less legibly.

## Provenance

Every run writes `run_manifest.json`: the scenario in full, every seed, the suite commit
(with a `-dirty` suffix if the tree had uncommitted changes), the Python and platform
versions, the versions of every package that can move a number, checksums of the inputs,
and checksums of the outputs.

```bash
appfl-bio-suite simulate fine-mapping --verify /scratch/fm/run_manifest.json
```

re-checksums the outputs and reports drift. That answers "did this rerun produce the same
data?" mechanically instead of by eyeballing summary statistics.

Only the bundles and the answer key are checksummed, not every intermediate. The pipeline's
scratch is large, uninteresting, and in places not bit-reproducible across PLINK builds;
including it would make `--verify` report drift that says nothing about the data.

A run also writes `pipeline_config.yaml` — the resolved upstream config — so the vendored
stages that this suite does not wrap can be pointed at the package without hand-writing
one.

## Running it at scale

`simulate` is one process. The published scenario is 100 loci × 15 architectures × 10
replicates over 150,000 individuals, and two of its stages shard over loci; a single
process would take days.

The stages are therefore also reachable one at a time, and
`scripts/fine-mapping/submit_polaris_simulate.pbs` runs the sharded sequence across nodes.
See [RUNBOOK.md](RUNBOOK.md) and `simulation/stages.py`.

For sizing, the published run's measured wall times — one Polaris node, 32 cores, chr1
only, 50,000 individuals per site:

| Stage | Wall time | Dominated by |
| --- | --- | --- |
| Sampling | 5–10 min | PLINK `--keep` I/O |
| Locus selection | 1–3 h | LD scoring; `locus_selection.n_workers` controls parallelism |
| Phenotype simulation | 30–60 min | the instance grid (100 loci × 15 architectures × 10 replicates requested) |
| QC | 5–10 min per site | |

Memory is the other axis: one 1 Mb window is roughly 5k variants × 50k samples × 4 bytes
≈ 1 GB peak, and locus selection fans out over joblib subprocesses, so its peak is about
`n_workers × 1 GB`. More workers is not automatically faster either — 64 measured about 2×
*slower* than 32 on this memory-bandwidth-bound workload; see
[reference/design.md](reference/design.md) for that measurement and the queue history.

What you give up is manifest strength: a sharded run's manifest records the scenario, the
seeds and the output checksums, but nothing witnessed each shard, so it cannot attest that
every shard ran the code it names. The manifest says so in its notes.

## Dependencies

Simulation is **coordinator-side only** and its dependencies live in the separate
`finemapping-sim` extra, so they never land on a partner:

```bash
pip install -e '.[finemapping-sim]' -c constraints.txt
scripts/fine-mapping/install_plink.sh      # PLINK 1.9, into vendor/bin/
scripts/fine-mapping/install_susiex.sh     # SuSiEx, into vendor/bin/
```

PLINK and SuSiEx are C++ binaries and not pip-installable. Both are coordinator-side: a
partner computes second moments in numpy and never fits anything, which is why the partner
extra is two packages. `appfl-bio-suite preflight --experiment fine-mapping` checks for
both before a run rather than after a queue wait.

## What partners receive, and what they are told

The bundle above, generated by `appfl-bio-suite partner-bundle fine-mapping --site <id>`,
with their config already filled in. They are told the data is simulated, what the files
are, and one verification command with its expected output. They are not told the causal
variants, and they are not told anything about another site.
See [the partner guide](../../partner/experiments/fine-mapping.md).
