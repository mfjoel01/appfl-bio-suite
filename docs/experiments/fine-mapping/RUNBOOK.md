# Federated fine-mapping — runbook

## Before you start

You need the coordinator install, both C++ binaries, and a data package.

```bash
pip install -e '.[finemapping-sim]' -c constraints.txt
scripts/fine-mapping/install_plink.sh      # PLINK 1.9  -> vendor/bin/
scripts/fine-mapping/install_susiex.sh     # SuSiEx     -> vendor/bin/
```

Both installers need outbound network access once and are idempotent. Neither is
pip-installable, and neither is a partner's problem — the site stage is numpy.

**Prove the install before recruiting anyone.** Skip to
[Validating without partners](#validating-without-partners); it takes about a minute and
exercises the entire chain including the byte layout SuSiEx's reader demands.

## 1. Generate and distribute data

```bash
appfl-bio-suite simulate fine-mapping --list-scenarios
appfl-bio-suite simulate fine-mapping --scenario three-site-hapnest --out /scratch/fm-run
```

For the published scenario this wants the HAPNEST pool staged and
`$FEDFM_HAPNEST_DIR` set — see [DATA.md](DATA.md) — and it is hours of work in one
process. At that scale use the sharded job instead:

```bash
qsub -A <allocation> -v FM_RUN_DIR=/scratch/fm-run \
     scripts/fine-mapping/submit_polaris_simulate.pbs
```

Then generate one bundle per partner and send it:

```bash
appfl-bio-suite partner-bundle fine-mapping --site site-north
```

**Keep `/scratch/fm-run/ground_truth/causal_manifest.tsv`.** It is the answer key, it is
deliberately absent from every bundle, and the aggregator needs it to score the run.

Record in `local/federation.yaml`, per site: their `data_dir`, their `output_dir`, their
`expected_samples`, and the path to the answer key under `causal_manifest`.

## 2. Preflight

```bash
appfl-bio-suite preflight --experiment fine-mapping
```

Fine-mapping adds one check the other experiments do not have: **external binaries**. It
is a hard failure when an experiment is named, because a missing SuSiEx is not discovered
until after every site has computed and transferred its second moments — the expensive
part of the run.

A bare `appfl-bio-suite preflight` demotes it to a warning, because a coordinator running
only the GWAS experiment is entitled to a clean report.

The `data` group reports each site honestly: validated where the path resolves here,
skipped where it does not, because `data_dir` is a path on a partner's cluster and from
here it is not merely absent but unknowable.

### The GA4GH checks

```bash
appfl-bio-suite preflight --check ga4gh --experiment fine-mapping
```

Offline, in about a second, and each prevents a failure that otherwise costs a scheduler
queue wait:

| check | fails when |
| --- | --- |
| `duo ontology` | the vendored DUO release will not load |
| `data use` | a site's declared terms do not permit this study |
| `drs` | a site's `drs_uri` does not resolve, or names a file rather than a bundle |
| `trs pin` | the installed site stage is not the version this federation pinned |

`data use` is the one that stops a launch. It is the same evaluation each site's worker
performs on its own copy of its terms — a site that refuses here will refuse there — so
this is how you find out before spending three sites' allocations to be told.

`appfl-bio-suite ga4gh duo check` runs just that one, and prints the study alongside the
per-site verdict. A site whose terms you have no copy of is reported as unchecked rather
than as permitted; its own worker still enforces them.

Skipped entirely when the experiment declares no `ga4gh` block.
**→ [../../coordinator/ga4gh.md](../../coordinator/ga4gh.md)**

## 3. Size the run before you launch it

**This is the step that decides whether the run completes.** Read it even if you skip the
rest.

A site's uplink is O(M²) per (locus, ancestry). At a realistic locus M is around 2,000
variants after the MAF filter, so:

| | |
| --- | --- |
| One ancestry block, float64 | ~32 MB |
| A site holding five ancestries, one locus | ~160 MB |
| Three sites, one locus | ~480 MB |
| Three sites, 79 loci | **~38 GB** |

A full sweep in one exchange will not complete. Shard it: launch the run several times,
bumping `locus_shard_index` each time, and concatenate the result tables.

```yaml
# local/federation.yaml
fine-mapping:
  locus_n_shards: 8
  locus_shard_index: 0      # 1, 2, ... 7 on subsequent launches
  uplink_gram_dtype: float32
```

`float32` halves the dominant term and is **exact**: Gram entries are integers bounded by
4n, and float32 holds every integer below 2²⁴, so at any realistic cohort size the round
trip is lossless. The trainer verifies the bound before downcasting and falls back to
float64 if it does not hold.

Both sharding fields are federation-wide rather than per-site, deliberately. Two sites on
different shards would each contribute to loci the other skipped, and every affected locus
would be fine-mapped on a smaller cohort than its reported `n` — with nothing downstream
able to detect it.

To find out what a real locus actually costs before committing, run one:

```yaml
  locus_limit: 1
```

The trainer logs each site's payload size, and the aggregator logs every site's block
count. Multiply.

## 4. Dry run

```bash
appfl-bio-suite run fine-mapping --dry-run --out-dir local/configs/generated
```

Writes the exact server and client configs APPFL will be handed, and stops. This is the
fastest way to see what a partner's endpoint will be sent. Check that every `*_path` is an
absolute path into the installed package (those are read on **your** machine and shipped as
source) and that every `data_dir` and `output_dir` is an absolute path on the **partner's**
cluster.

## 5. Launch

```bash
appfl-bio-suite run fine-mapping
```

Preflight runs first; a hard failure stops the launch.

## What a healthy run looks like

Each site logs its own progress, then the coordinator logs what arrived:

```
Site1: reading site index from /home/finemap_svc/data/anl
Site1: 636 of 533532 variants recoded to the reference allele order
Site1: locus L0000 -- 5 ancestry block(s), 150 instance(s)
Site1: 5 genotype block(s), 750 phenotype block(s), 163.4 MB uplink

fine-mapping across 3 site(s)
  Site1: n=50,000, ancestries ['AFR','AMR','CSA','EAS','EUR'], 5 genotype block(s), ...
  L0000: 12 block(s) from 3 site(s) -> AFR(n=60000, M=1987), EUR(n=36500, M=1994), ...
wrote 150 result row(s) -> local/output/fine-mapping/data/fed_fm_results.tsv
```

Four things to check:

1. **Recoded variant counts are nonzero and differ between sites.** Around 6% is expected
   on real chr1 data. All-zero means the reference list was cut from the same filesets the
   sites hold, which defeats the check.
2. **Every ancestry column's `n` matches the composition you recorded.** The aggregator
   pools whoever contributed; a site whose bundle went missing shows up as a smaller `n`,
   not as an error.
3. **Uplink sizes match your estimate.** If they are 10× larger, the MAF filter is not
   biting and M is larger than you sized for.
4. **The dropped-variant warning is rare.** See below.

## When it is not healthy

**`SuSiEx: ... aborting`, immediately, on every instance.** The LD panel the coordinator
wrote does not match what the reader expects. The reader checks that `.ld.bin` is exactly
`M² × 4` bytes against the `_ref.bim` line count, that every bim line's chromosome equals
`--chr`, and that the `.frq` repeats each bim line's alleles. It aborts rather than
degrading. Run the loopback below — it exercises exactly this and nothing else has to be
working for it to fail informatively.

**`N window variant(s) dropped as not fully observed at some site`.** Expected in small
numbers; HAPNEST's missing rate costs about four variants in two thousand. Where it fires,
the federated and centralized paths fine-map slightly different variant sets and Corollary
1's equality no longer bites for that locus — which is why it is a warning rather than
silence. A *large* count means a site's fileset has real missingness, and the run's
comparability to the centralized baseline is gone.

**`site X disagrees with site Y on the coding of N variant(s)`.** Harmonization failed.
Either the sites were cut from different filesets, or their `reference_variants.tsv` files
differ. Compare the bundles' `run_manifest.json` checksums.

**`no credible sets` on every instance, with no error.** SuSiEx ran and legitimately found
nothing. Check the `min_p_*` columns in the results: if the best marginal p-value is above
`pval_thresh` in every column, the signal is genuinely too weak — too few individuals per
ancestry, or too low a per-locus heritability. This is a real result, not a failure.

**A worker dies partway with no traceback.** Almost always memory. One `G` at a large locus
is tens of megabytes and the trainer holds every block for a locus at once. Lower
`locus_limit`, or raise the MAF filter to cut M.

**Everything else** — endpoints, identity mapping, version skew, `AF_UNIX path too long` —
is in [../../coordinator/troubleshooting.md](../../coordinator/troubleshooting.md), which
is shared across experiments.

## Reading the output

```
local/output/fine-mapping/
├── data/
│   ├── fed_fm_results.tsv              one row per (locus, architecture, replicate)
│   ├── fed_fm_rollup_by_architecture.tsv
│   ├── fed_fm_rollup_by_stratum_rg.tsv
│   ├── fed_fm_site_summary.csv         what each site sent
│   ├── ga4gh_provenance.json           what was dispatched, and what each site attests to
│   └── drs_outputs.json                a DRS object per results file
├── graphs/                             the four figures
└── logs/
```

The columns that matter, per instance:

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
survive after a run.

`ga4gh_provenance.json` carries the two halves of the record and keeps them apart, because
they are attested by different parties. `dispatched` is what this coordinator sent — the
data use request, the tool pin, the DRS object it expected each site to hold. `attested` is
what each site says it actually did: its DUO decision term by term, the object it verified
its bundle against, the pin it was handed. The coordinator can only record the first; the
site is the only witness to the second.

They are compared on arrival. A site that computed over a **different DRS object** than the
run names is a hard failure, because every number it contributed would otherwise be
attributed to data it did not read.

**The comparison that matters** is this table against the centralized baseline's, column
for column. They have identical schemas precisely so that they can be diffed. Produce the
baseline with:

```bash
scripts/fine-mapping/run_stage.py centralized --config /scratch/fm-run/pipeline_config.yaml
```

or, at scale:

```bash
qsub -A <allocation> -v FM_RUN_DIR=/scratch/fm-run,FM_STAGE=centralized \
     scripts/fine-mapping/submit_polaris_finemap.pbs
```

Any disagreement beyond floating-point summation order is a bug in the federated path. See
[ABOUT.md](ABOUT.md).

Two things about that diff are worth knowing before you read one as a failure:

**The `min_p_*` columns agree only to about 1e-6, and that is expected.** The centralized
path reads its marginal p-values out of plink2's `.glm.linear` text, which carries roughly
six significant figures; the federated path computes them in-process at full double
precision. Every other column is exact. Compare `min_p_*` with a tolerance, not with
`diff` — a relative difference around 1e-6 in those six columns and nowhere else is the
signature of agreement, not of drift.

**Both sides must use the same MAF filter, and by default they do.** It comes from the
scenario's `fine_mapping.maf`, which `simulate` writes into `pipeline_config.yaml`; the
centralized stage, the standalone federated stage, and a `--driver serial` loopback run
all read it from there. Passing `--maf` to one side and not the other fine-maps two
different variant sets, and the resulting credible-set disagreement looks exactly like the
bug this comparison exists to find.

## Validating without partners

The loopback run: the full federation on one machine, no Globus, no scheduler.

```bash
appfl-bio-suite simulate fine-mapping --scenario ci-tiny --out /tmp/fm-ci
appfl-bio-suite run fine-mapping --config loopback --driver serial --data-root /tmp/fm-ci
```

It needs no `federation.yaml` — one is synthesized describing the simulated sites on this
machine — and it wires up the answer key automatically, because a loopback run is the one
case where the coordinator legitimately holds it: it simulated the data seconds ago.
Without that the run would complete and report credible sets it could not score, which
looks identical to success whether or not the statistics are right.

If it completes, the install is correct end to end: the simulation produced bundles in the
layout the loader expects, each site computed its second moments, the coordinator pooled
and standardized them, SuSiEx read the panel the coordinator wrote by hand, and credible
sets came back scored against the truth.

Expect three loci, each capturing its causal variant with PIP 1.0 — the fixture uses a
deliberately enormous effect size so that a credible set actually comes back. Zero credible
sets means something is wrong; see above.

## The other stages

`scripts/fine-mapping/run_stage.py` reaches the coordinator-side stages that are not part
of the federation — the centralized baseline, the standalone federated path, the validation
harness, per-site QC, and the figures.

```bash
scripts/fine-mapping/run_stage.py                      # list them
scripts/fine-mapping/run_stage.py validation --config /scratch/fm-run/pipeline_config.yaml
```

The shardable ones take `--n-shards N --shard-index I`, then `--merge --n-shards N` to
reduce. None of them contacts a partner.
