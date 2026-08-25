# Part 2 — Fed-Heart-Disease: getting your data

Assumes [Part 1](../endpoint-setup.md) is done and your endpoint is running.

This part is short: you download a small public dataset, verify it loads, and tell us one
number.

---

## About the data

The dataset is **Fed-Heart-Disease** from the public
[FLamby](https://github.com/owkin/FLamby) benchmark — about **40 KB** of tabular
cardiology records from the classic UCI Heart Disease study, released under **CC BY 4.0**
with no registration or approval step.

It involves **none of your own patient data**. The benchmark has four natural centres —
the four hospitals of the original study — and you have been assigned **centre
{{ center }}**. The centres genuinely differ from one another, which is the point of the
benchmark.

Expected size of your training split: **{{ expected_train_samples }} samples**. You will
confirm that number below.

---

## Step 1 — Install FLamby

Run this **as `{{ service_account }}`**, not as yourself. FLamby records where it put the
data *inside its own package directory*, so whoever runs the download determines who can
find the data afterwards. Running it as the service account avoids a whole class of
"file not found" failures that only appear mid-run.

```bash
sudo -iu {{ service_account }}     # or however you switch to the service account
conda activate appfl_env           # the same environment from Part 1
```

Get FLamby from GitHub — it is not published on PyPI:

```bash
git clone https://github.com/owkin/FLamby.git
cd FLamby
```

Now install it **without running its `setup.py`**:

```bash
python -c "import site, pathlib; \
  pathlib.Path(site.getsitepackages()[0], 'flamby.pth').write_text(str(pathlib.Path.cwd()))"
python -c "import flamby; print('flamby importable from', flamby.__file__)"
```

That writes a `.pth` file pointing at your checkout, which puts FLamby on the import path
and nothing else.

> ### Do NOT run `pip install -e .` or `make install` here
>
> This is deliberate and it matters.
>
> FLamby's `setup.py` hooks `install`, `develop`, **and `egg_info`** to unconditionally
> shell out to an external package index. The `egg_info` hook means it fires on mere
> metadata generation — so even a resolver *looking at* the package triggers it. It
> reaches a host many clusters cannot get to, and it can move packages that our version
> pins depend on.
>
> `make install` is worse: it builds an entirely separate conda environment from FLamby's
> own unpinned `environment.yml` — an environment your endpoint never activates, so the
> workers would not see it anyway.
>
> The `.pth` approach above sidesteps all of that and leaves every pinned version
> untouched.

Confirm nothing moved:

```bash
globus-compute-endpoint version    # must STILL print {{ globus_compute_version }}
```

If that number changed, stop and send it to us before going further.

## Step 2 — Download your centre's data

Still as `{{ service_account }}`:

```bash
mkdir -p {{ data_dir }}
cd flamby/datasets/fed_heart_disease/dataset_creation_scripts
python download.py --output-folder {{ data_dir }}
```

A few seconds, about 40 KB.

That path is a suggestion, not a requirement — FLamby records where the data went inside
its own package directory and resolves it internally, so we never need to know. Any
location `{{ service_account }}` can read and write will do. If you use a different one,
you do not need to tell us.

If you ever move the data afterwards, tell FLamby where it went rather than just moving
the folder:

```bash
python update_config.py --new-path <new absolute path>
```

## Step 3 — Verify, and send us the number

Still as `{{ service_account }}`, in the same environment:

```bash
python -c "
from flamby.datasets.fed_heart_disease import FedHeartDisease
d = FedHeartDisease(train=True, center={{ center }}, pooled=False)
print('OK, samples:', len(d))
"
```

**Expected output: `OK, samples: {{ expected_train_samples }}`**

If you get that number, your data step is complete. There is nothing to transfer and no
path for us to know — FLamby resolves it internally on your machine.

If you get a *different* number, tell us rather than adjusting anything. It means the
centre assignment or the download differs from what we expect, and we would rather find
that now than discover mid-run that two sites trained on the same data.

---

## Step 4 — Send us

Just one line:

> Fed-Heart-Disease centre {{ center }} loads, {{ expected_train_samples }} training samples.

That is everything. Combined with the endpoint UUID from Part 1, we have what we need.

**Leave the endpoint running.** We drive the run from our side; there is nothing further
for you to do, and nothing to monitor.

---

## Troubleshooting

**`ModuleNotFoundError: flamby`**
Step 1 ran in the wrong environment or as the wrong user. Re-run the verification command
as `{{ service_account }}` with the Part 1 environment active.

**FLamby cannot find the dataset / `dataset_location.yaml` error**
The download ran as a different user, or into a path `{{ service_account }}` cannot read,
or the folder was moved afterwards. Fix the path with `update_config.py --new-path <abs
path>` in `dataset_creation_scripts/`, or re-run the download as `{{ service_account }}`.

**The sample count does not match**
Do not adjust it to fit. Send us what you got — the discrepancy is information.

**`globus-compute-endpoint version` changed after Step 1**
Something moved a pinned package. Send us the output of
`pip list | grep -E '^(click|globus-sdk|torch|numpy) '` and we will send the exact restore
command.

Contact: {{ coordinator_contact }}
