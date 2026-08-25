# Part 2 — Federated GWAS: placing your data

Assumes [Part 1](../endpoint-setup.md) is done and your endpoint is running.

---

## About the data

You will receive a **data bundle** from {{ coordinator_organization }} containing your
site's cohort. It is synthetic data prepared for this study — not patient data, yours or
anyone's.

Your bundle contains **{{ expected_samples }} individuals**. Every site in the federation
gets a different, non-overlapping set of individuals, and all sites share the same variant
set — that is what makes the results combinable.

**What leaves your cluster:** per-variant effect sizes, standard errors and allele
frequencies, plus two summary numbers. Genotypes and phenotypes never leave. The analysis
runs entirely on your machine and returns a payload of a few megabytes.

---

## Step 1 — Place your bundle

Unpack the archive we sent you somewhere `{{ service_account }}` can read:

```bash
sudo -iu {{ service_account }}
mkdir -p {{ data_dir }}
cd {{ data_dir }}
tar xzf /path/to/the/bundle/you/received.tar.gz --strip-components=1
```

Your data directory is:

```
{{ data_dir }}
```

That exact path is already in the configuration we hold for you, so if you put it
somewhere else, tell us.

## Step 2 — Check the contents

The directory must contain these six files, **directly inside it** — not in a nested
subdirectory:

```
EUR.synthetic.100k.ld.maf.bed
EUR.synthetic.100k.ld.maf.bim
EUR.synthetic.100k.ld.maf.fam
phenotypes_gwas.csv
phenotypes_pgs_eval.csv
covariates.csv
```

A few additional files may be present. They are harmless and unused.

```bash
ls {{ data_dir }}
```

## Step 3 — Verify, and send us the number

Still as `{{ service_account }}`, in the Part 1 environment:

```bash
python -c "
from appfl_bio_suite.experiments.gwas.dataset import get_dataset
d, _ = get_dataset(data_dir='{{ data_dir }}', site_id='{{ client_id }}')
print('OK, samples:', len(d))
"
```

**Expected output: `OK, samples: {{ expected_samples }}`**

That command validates every required file is present and readable **as the account that
will actually run the analysis**, which is the check that matters. If it prints the
expected number, you are done.

If it reports missing files, the archive was probably unpacked one level too deep — the
six files must sit directly in `{{ data_dir }}`.

If it prints a *different* sample count, tell us rather than adjusting anything. It means
you have a different bundle from the one we recorded for you.

---

## Step 4 — Send us

Just one line:

> {{ client_id }} data in place at {{ data_dir }}, {{ expected_samples }} samples.

**Leave the endpoint running.** We drive the run from our side. It is a single round —
your cluster receives one task, runs the local analysis, and returns summary statistics —
so there is nothing to monitor and nothing further for you to do.

---

## Troubleshooting

**`FileNotFoundError` naming missing files**
The archive was unpacked into a subdirectory. The six files must be directly inside
`{{ data_dir }}`. The error message lists exactly which ones are missing and what it found
instead.

**`data_dir does not exist`**
The path is resolved on **your** cluster by the worker, not on ours. It must be an
absolute path readable by `{{ service_account }}`.

**Permission denied**
`{{ service_account }}` cannot read `{{ data_dir }}`. Check ownership and permissions on
the directory and every parent of it.

**The sample count does not match**
Do not adjust anything. Send us the number you got.

**The task runs for a long time**
A full analysis over all variants takes a while — it is a genome-wide scan over
{{ expected_samples }} individuals. It is one task, not one per round. If it exceeds your
walltime, tell us and we will lower the variant count for a first pass.

Contact: {{ coordinator_contact }}
