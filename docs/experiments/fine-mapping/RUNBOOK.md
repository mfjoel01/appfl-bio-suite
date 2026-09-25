> Current protocol: [September scientific corrections and rerun](SCIENTIFIC_RERUN.md). Its validation gates have passed and the run is complete; the findings and limitations are in that document's [outcome section](SCIENTIFIC_RERUN.md#outcome--17-september-2026). Numerical results below that predate the rerun are historical.

# Fine-mapping — how to run it

Three things you can run, in increasing order of what they need.

| | needs | gives |
| --- | --- | --- |
| **Loopback** | nothing but the repo | the whole federation in one process, no network |
| **Centralized baseline** | a data package | the reference fit, all data in one place |
| **Federated** | a data package, or partners | the same fit, computed without pooling |

## 1. Make a data package

```bash
appfl-bio-suite simulate fine-mapping --list-scenarios
appfl-bio-suite simulate fine-mapping --scenario ci-tiny --out /scratch/fm
```

`ci-tiny` runs in seconds and produces no meaningful genetics — it is for checking the
install. The published scenario needs a staged genotype pool and is not a single-process
job; see [DATA.md](DATA.md).

The package it writes is self-contained: per-site bundles, ground truth, selected loci,
and a run manifest. Record the **manifest hash alongside the commit** — the commit says
which code ran, the manifest says which data it ran on, and without both a run is not
reproducible.

## 2. Check before you launch

```bash
appfl-bio-suite preflight fine-mapping --data-root /scratch/fm
```

Reports the per-locus uplink, the ancestry columns that survive `min_gwas_n`, and the
GA4GH checks (consent terms, bundle digests, tool pin). A DRS mismatch is fatal rather
than a warning: if a site computed over an object other than the one the run's provenance
claims, every number in the results table is attributed to the wrong data.

## 3. Run it

```bash
# everything in one process — no network, no scheduler
appfl-bio-suite run fine-mapping --config loopback --driver serial --data-root /scratch/fm
```

For a real federation, point the config at your partners' endpoints and drop
`--driver serial`. Sharding is by locus; a run larger than one exchange can carry is split
with `locus_n_shards` / `locus_shard_index`.

The coordinator-side stages that are **not** the federation — the centralized baseline,
the standalone federated driver, validation, QC, figures — go through:

```bash
scripts/fine-mapping/run_stage.py                       # list them
scripts/fine-mapping/run_stage.py centralized --config /scratch/fm/pipeline_config.yaml
scripts/fine-mapping/run_stage.py federated   --config /scratch/fm/pipeline_config.yaml \
    --n-shards 4 --shard-index 0
scripts/fine-mapping/run_stage.py federated   --config /scratch/fm/pipeline_config.yaml \
    --merge --n-shards 4
```

None of them contacts a partner. `scripts/fine-mapping/submit_polaris_*.pbs` wrap the
sharded versions for a scheduler; adapt the directive block to your own site.

## Participation arms

To ask what federating buys, run the same loci with different cohorts taking part:

```bash
python -m appfl_bio_suite.experiments.fine_mapping.arms \
    --base-config <pipeline_config.yaml> \
    --out-dir <arms/> \
    --arms anl,covenant,mbzuai,federation_50k,ld_borrowed \
    --federation-results <the all-sites fm_results.tsv>
```

It writes one config per arm, runs each through the **same** vendored driver, and merges
the results into `fm_results_by_arm.tsv` with an `arm` column. `--build-only` writes the
configs and stops, which is the cheap way to check what each arm's columns will be before
committing to a sweep.

Two arms are not plain site subsets and are worth knowing about:

- `federation_50k` down-samples every site to a third, stratified within site and
  ancestry, into a shadow `processed/` tree whose PLINK filesets are symlinks. It holds
  analyzed sample size at exactly 50,000 after ancestry exclusions. The corrected rerun
  uses three sampling seeds and five ancestry columns. Remaining differences are
  composition/participation effects, including differences in genetic signal and noise.
- `ld_borrowed` pairs one site's summary statistics with another's LD panel, restricted to
  the ancestries both hold. It is the shortcut the O(M²) uplink exists to avoid.

Then draw the comparison:

```bash
python -m appfl_bio_suite.experiments.fine_mapping.figures.render \
    --results <arms/>/fm_results_by_arm.tsv --by-arm <arms/>/fm_results_by_arm.tsv \
    --data-root <package> --out-dir <figures/>
```

## Reading the output

```
<out>/
├── data/
│   ├── fed_fm_results.tsv              one row per (locus, architecture, replicate)
│   ├── fed_fm_rollup_by_*.tsv
│   ├── fed_fm_site_summary.csv         what each site sent
│   ├── ga4gh_provenance.json           what was dispatched, what each site attests to
│   └── drs_outputs.json
├── graphs/                             figures, drawn in-process by the aggregator
└── logs/
```

Per instance, the columns that matter:

| Column | Meaning |
| --- | --- |
| `n_credible_sets` | how many credible sets SuSiEx returned |
| `cs_sizes_json` | their sizes — smaller is sharper resolution |
| `any_causal_captured` | did any credible set contain a true causal variant |
| `best_cs_size` | size of the smallest credible set that did |
| `causal_pip_max` | posterior inclusion probability of the true causal variant |
| `min_p_<ANCESTRY>` | best marginal p-value in that column — diagnose "no sets" here |
| `converged` | SuSiEx produced sets or PIPs at all |

`fed_fm_site_summary.csv` is the only place the per-site uplink and harmonization counts
survive a run.

**An empty result is not a failure.** Where no signal is resolvable, returning no credible
set is the correct answer; it should concentrate at low heritability rather than spread
evenly across the design.

## Figures

The aggregator draws the paper-analogue set into `graphs/`. For everything —
exploratory figures over the package, and the federation-cost figures — use the batch
renderer:

```bash
python -m appfl_bio_suite.experiments.fine_mapping.figures.render \
    --results <out>/data/fed_fm_results.tsv \
    --data-root /scratch/fm \
    --out-dir <out>/figures
```

Only `--results` and `--out-dir` are required; each further input unlocks specific figures
and is skipped with a note when absent. Two are worth knowing about:

- `--detail-dir` — output of `figures.detail`, a sampled re-run that keeps SuSiEx's
  working directories. The results table summarises each fit to one row and deletes them,
  which discards the per-credible-set records **coverage needs**, the population-specific
  causal probabilities, and the per-variant PIPs.
- `--federated` + `--run-log` — a second results table over the same package, for the
  federated-vs-centralized parity figure.

Full index and the mapping onto the SuSiEx paper's own figures:
[figures/README.md](../../../src/appfl_bio_suite/experiments/fine_mapping/figures/README.md).

## When it is not healthy

| Symptom | Where to look |
| --- | --- |
| no credible sets anywhere | `min_p_*` — if nothing is near 5×10⁻⁸ there is no signal to fine-map |
| SuSiEx aborts on the LD panel | the panel is byte-exact or it is rejected; check `M² × 4` against the `_ref.bim` line count |
| federated ≠ centralized | first check the run log for variants dropped as *not fully observed* — Corollary 1 assumes both paths see the same variant list |
| a site's uplink is huge | it is O(M²) per ancestry; shard by locus |
| DRS mismatch | the site ran on a different bundle than the one cut for it — re-cut and re-dispatch |
