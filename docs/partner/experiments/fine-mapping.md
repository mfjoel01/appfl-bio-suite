# Part 2 — Federated fine-mapping: placing your data

Assumes [Part 1](../endpoint-setup.md) is done and your endpoint is running.

---

## About the data

You will receive a **data bundle** from {{ coordinator_organization }} containing your
site's cohort. It is synthetic data prepared for this study — not patient data, yours or
anyone's.

Your bundle contains **{{ expected_samples }} individuals**. Every site in the federation
gets a different, non-overlapping set of individuals, and all sites are given the same
variant list and the same set of genomic regions to analyse — that is what makes the
results combinable.

**What leaves your cluster:** sums and sums of products over your own individuals — a
matrix of variant-by-variant co-occurrence counts, per-variant totals, and a handful of
scalars, for each region and each ancestry group you hold. Genotypes and phenotypes never
leave, and neither does anything about any individual. The analysis runs entirely on your
machine.

**The payload is larger than you may expect** — hundreds of megabytes rather than a few,
because the statistic being combined is a variant-by-variant matrix rather than one number
per variant. That is the nature of this analysis, not a fault. If your endpoint has an
outbound transfer limit that this would run into, tell us and we will split the run into
smaller pieces; that is a change on our side only.

**Nothing to install beyond Part 1.** The analysis is numpy and pandas. There is no
statistical package, no genetics toolchain, and no compiler involved on your side.

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

**Check the free space first.** The genotype file is the bulk of it and is not small.

## Step 2 — Check the contents

The directory must contain these six files and one subdirectory, **directly inside it** —
not in a nested subdirectory:

```
site_genotypes.bed
site_genotypes.bim
site_genotypes.fam
site_manifest.tsv
reference_variants.tsv
selected_loci.tsv
phenotypes/          (one .pheno file per analysis instance — there are many)
```

```bash
ls {{ data_dir }}
ls {{ data_dir }}/phenotypes | wc -l
```

## Step 3 — Verify, and send us the number

Still as `{{ service_account }}`, in the Part 1 environment:

```bash
python -c "
from appfl_bio_suite.experiments.fine_mapping.dataset import get_dataset
d, _ = get_dataset(data_dir='{{ data_dir }}', site_id='{{ client_id }}')
print('OK, samples:', len(d))
print('    ancestries:', d.composition)
print('    regions:', len(d.loci), ' instances:', len(d.available_instances))
"
```

**Expected output starts with: `OK, samples: {{ expected_samples }}`**

That command validates every required file is present and readable **as the account that
will actually run the analysis**, which is the check that matters. It also checks that the
files agree with each other — that every individual in the genotype file has an ancestry
label, and that the phenotype files cover them — because a bundle whose parts came from
different runs is the failure that is hardest to spot later.

If it reports missing files, the archive was probably unpacked one level too deep — the
six files must sit directly in `{{ data_dir }}`.

If it prints a *different* sample count, tell us rather than adjusting anything. It means
you have a different bundle from the one we recorded for you.

---

## Step 4 — Send us

Just one line:

> {{ client_id }} data in place at {{ data_dir }}, {{ expected_samples }} samples.

**Leave the endpoint running.** We drive the run from our side. It is a single round —
your cluster receives one task, computes its summaries, and returns them — so there is
nothing to monitor and nothing further for you to do.

We may run it more than once, covering different genomic regions each time. Each is again
a single round and needs nothing new from you.

---

## Troubleshooting

**`FileNotFoundError` naming missing files**
The archive was unpacked into a subdirectory. The six files must be directly inside
`{{ data_dir }}`. The error lists exactly which are missing and what it found instead.

**`data_dir does not exist`**
The path is resolved on **your** cluster by the worker, not on ours. It must be an
absolute path readable by `{{ service_account }}`.

**Permission denied**
`{{ service_account }}` cannot read `{{ data_dir }}`. Check ownership and permissions on
the directory and every parent of it.

**`individual(s) ... have no ancestry label`**, or **`names individuals that are not in`**
The genotype file and `site_manifest.tsv` describe different cohorts, which normally means
two bundles' files were mixed. Do not try to reconcile them — tell us, and we will re-cut
your bundle.

**The sample count does not match**
Do not adjust anything. Send us the number you got.

**The task runs out of memory**
The analysis holds one region's genotypes and their cross-products at a time. If a worker
is killed, tell us the region count in the error and how much memory your workers get; we
lower the number of regions per task from our side.

**The task runs for a long time**
It reads a genomic region out of your genotype file and forms matrix products over your
cohort, once per region. Long is normal; hours is worth telling us about.

Contact: {{ coordinator_contact }}
