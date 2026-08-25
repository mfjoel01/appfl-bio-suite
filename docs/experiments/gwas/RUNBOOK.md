# Federated GWAS — runbook

## Before you start

- Every partner has completed setup and sent an endpoint UUID.
- Every partner has their data bundle in place and has confirmed the sample count.
- `local/federation.yaml` has each site's UUID, `data_dir`, and `expected_samples`.
- You are on a machine that can stay up for the run, in `tmux` or `screen`.

## 1. Generate and distribute data

```bash
appfl-bio-suite simulate gwas --scenario three-site-skewed --out local/data/gwas/runs/three-site-skewed
```

Point `--out` at a per-run directory like this one, so repeated simulations never mix
with staged inputs or with each other.

Each `local/data/gwas/runs/three-site-skewed/SiteN/data/` is one partner's bundle. Send
each site theirs.

**Keep `run_manifest.json`.** It records the scenario, every seed, the suite commit, and
input and output checksums — it is what lets you prove later which data produced which
result. Record its hash in the run of record.

Read [DATA.md](DATA.md) before using any output for a publication.

## 2. Preflight

```bash
appfl-bio-suite preflight --experiment gwas
```

Checks the environment, version pins, that every shipped config parses, and that every
endpoint reports online. Hard failures block a launch; warnings do not.

Then verify each partner for real:

```bash
appfl-bio-suite endpoint smoke --experiment gwas
```

The `user` in each result must be the experiment's service account. If it is the account
that started the endpoint, identity mapping is not taking effect. `endpoint status` alone
proves only that a daemon is connected.

**Every site must agree on `variant_scaling`.** The aggregator requires all sites to be on
the same variant set and will refuse to combine mismatched payloads — but it is much
cheaper to catch that here than after every site has run a full genome-wide scan.

You do not keep them in step by hand: `variant_scaling` and `hit_p_threshold` are set once
under `experiments.gwas` in federation.yaml, and every site's generated config takes them
from there.

## 3. Dry run

```bash
appfl-bio-suite run gwas --dry-run
```

Check each site's `data_dir` is the path that site actually confirmed, and that
`variant_scaling` is identical across all of them.

## 4. Launch

```bash
appfl-bio-suite run gwas
```

One round. Each site receives exactly one task, runs a complete local GWAS, and returns
summary statistics.

**Expect it to be slow at each site and quiet in between.** A genome-wide scan over tens
of thousands of individuals takes real time — this is not a per-round exchange where
silence means something is wrong.

## What a healthy run looks like

```
[setup] 3 client(s), 1 global round(s)
[train] dispatching round 1 to all clients
[train] Site2 returned round 1/1
[train] Site1 returned round 1/1
[train] Site3 returned round 1/1
meta-analysis across 3 site(s)
Federated Learning Training Completed!
```

Sites return in whatever order they finish; smaller cohorts usually come back first.

## When it is not healthy

**A site never returns.** Either its scheduler has not started the block, or the analysis
exceeded its walltime. A full scan can take a while — if a site's walltime is tight, lower
`variant_scaling` in federation.yaml for a first pass, which lowers it for every site
at once.

**`site payloads carry N variants but the metadata describes M`.** Sites analyzed
different variant sets. Almost always mismatched `variant_scaling`, or one site holding a
bundle from a different simulation run. Compare their `run_manifest.json` checksums.

**`missing required input files`.** That site's bundle was unpacked one level too deep.
The six required files must sit directly in `data_dir`. The error names exactly which are
missing and lists what it found.

**`a polygenic score has zero variance`.** Every effect estimate at that site was NaN,
which usually means genotypes and phenotypes do not actually correspond — a bundle
assembled from mismatched parts.

**`422` or `403` on dispatch.** Identity mapping. See
[troubleshooting](../../coordinator/troubleshooting.md#identity-mapping) — the two have
different causes.

## Reading the output

Under the aggregator's `output_dir`:

```
data/appfl_meta_gwas_bmi.csv.gz     pooled per-variant statistics, BMI
data/appfl_meta_gwas_t2d.csv.gz     pooled per-variant statistics, T2D
data/appfl_meta_gwas_hits.csv       significant variants (or the top 100 if none)
data/appfl_site_pgs_metrics.csv     per-site score metrics
data/appfl_meta_summary.csv         one-line summary of the whole run
graphs/                             Manhattan and QQ plots, meta and per-site
```

Start with `appfl_meta_summary.csv`: client count, total N, and the weighted score
metrics. Then the QQ plots — deviation from the diagonal at the tail is signal; deviation
along its whole length suggests uncontrolled structure or a systematic problem rather than
a discovery.

Per-site plots are rendered by the coordinator from the returned summary statistics, so
you can compare sites directly without asking anyone to send you anything.

For a run of record, capture the site count, total N, variant count, scenario name, the
**simulation manifest hash**, the suite commit, per-site and pooled metrics, and hit
counts. Add it to [ABOUT.md](ABOUT.md#record-of-runs).

## Validating without partners

```bash
appfl-bio-suite simulate gwas --scenario ci-tiny --out /tmp/gwas-ci
appfl-bio-suite run gwas --config loopback --driver serial --data-root /tmp/gwas-ci
```

Runs the whole chain on one machine in seconds: simulation, per-site data, local analysis
at two simulated sites, and the meta-analysis. Proves the install end to end.

It does **not** exercise the transport — serialization, identity mapping and scheduler
behaviour are exactly what a loopback run cannot test. Use `endpoint smoke` for those.
