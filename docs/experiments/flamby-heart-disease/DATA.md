# Fed-Heart-Disease — obtaining the data

**Each site downloads its own data. Nothing is distributed by the coordinator, and there
is no simulation stage for this experiment.**

That is the main reason this experiment is cheap to onboard a partner onto: there is no
data transfer, no data-use agreement, and no per-site bundle to build.

---

## What the data is

The **UCI Heart Disease** dataset, packaged by FLamby with its natural four-centre split.
The centres are the four collection sites of the original study — Cleveland, Hungary,
Switzerland, and VA Long Beach — so the heterogeneity between them is real and documented
rather than synthetic.

| | |
| --- | --- |
| Total size | ~40 KB |
| Centres | 4 |
| Features | 13 |
| Task | Binary classification |
| Licence | CC BY 4.0 — no registration, no approval step |

Approximate training-split sizes per centre: 199, 172, 30, 85. Confirm the exact number
for your assigned centre with the verification command below and record it in
`federation.yaml` as `expected_train_samples`, so preflight can catch a site that has
silently loaded the wrong centre.

## Licence and attribution

FLamby ships **no data** — it ships downloaders, and every dataset requires accepting the
upstream licence. Fed-Heart-Disease is the exception that makes it practical: CC BY 4.0
with no gate.

Attribute the original UCI Heart Disease authors in any write-up.

## How a site obtains it

The full procedure, with the exact commands, is in the partner-facing guide:
[docs/partner/experiments/flamby-heart-disease.md](../../partner/experiments/flamby-heart-disease.md).

In outline:

1. Clone FLamby from GitHub. It is not published on PyPI. The coordinator's reference
   checkout is `https://github.com/owkin/FLamby` at commit `edacf54d` ("Update README.md
   (#306)", 2024-06-19) — pin to it if upstream ever changes in a way that matters; the
   `fed_heart_disease` loader and downloader have been stable for years.
2. Install it **without running its `setup.py`** — via a `.pth` file pointing at the
   checkout.
3. Run FLamby's downloader for the heart-disease dataset.
4. Verify by loading the assigned centre and checking the sample count.

### Why not `pip install -e .`

FLamby's `setup.py` hooks `install`, `develop`, **and `egg_info`** to unconditionally
shell out to an external package index. The `egg_info` hook means it fires on mere
metadata generation, so even a resolver inspecting the package triggers it. It reaches a
host many clusters cannot get to, and it can move packages our version pins depend on.

`make install` is worse: it builds a separate conda environment from FLamby's own unpinned
`environment.yml` — an environment the endpoint never activates.

The `.pth` approach bypasses `setup.py` entirely and leaves every pin byte-identical.

### Run the download as the service account

FLamby records the downloaded data path *inside its own installed package directory*, so
whoever runs the download determines who can find the data afterwards.

If the download runs as a person and tasks run as a service account, the worker may read
a pointer to a path it cannot reach — and fail late, mid-round, with an error that looks
nothing like a permissions problem.

## Verifying

```bash
python -c "
from flamby.datasets.fed_heart_disease import FedHeartDisease
d = FedHeartDisease(train=True, center=<N>, pooled=False)
print('samples:', len(d))
"
```

Run **as the service account**, in the environment the endpoint's `worker_init` activates.
A count that differs from the expected one means the centre assignment or the download
differs from what the coordinator recorded — worth resolving before a run rather than
discovering that two sites trained on the same centre.
