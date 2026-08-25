# Fed-Heart-Disease — design and record of runs

## What this experiment is

A **cross-silo federated learning benchmark** run across independently administered HPC
centres, using [FLamby](https://github.com/owkin/FLamby) — the healthcare federated
learning benchmark from Owkin (NeurIPS 2022 Datasets & Benchmarks).

**The scientific point is not the model.** FLamby datasets are small, public, and already
have published centralized and federated baselines. The point is the **infrastructure
claim**: that APPFL over Globus Compute can drive a genuine multi-round federated training
job across separately administered clusters on different continents — with no shared
secret, no data movement, and no VPN — and reach the published baseline.

Fed-Heart-Disease is the right vehicle precisely because it is boring and reproducible:

- **Natural splits.** The client partition is not synthetic IID slicing. Each centre is a
  real hospital collection site, so the heterogeneity is genuine and documented.
- **A published number to hit.** The FLamby paper reports pooled, local-only, and
  federated baselines per dataset.
- **Small and CPU-only.** 40 KB of tabular data and a logistic regression. Every failure
  is unambiguously an infrastructure failure, never a compute failure — which is exactly
  the property you want when the infrastructure is what you are testing.
- **No registration gate.** The upstream data is CC BY 4.0 with no approval step, so a
  partner can be running the same afternoon rather than waiting weeks on a data-use
  agreement.

## Design

| | |
| --- | --- |
| Task | Binary classification (presence of heart disease) |
| Data | UCI Heart Disease, via FLamby. ~40 KB, 4 natural centres |
| Model | Logistic regression: `Linear(13, 1)` + sigmoid |
| Loss | Binary cross-entropy |
| Metric | Accuracy |
| Federation | Multi-round, synchronous FedAvg |
| Local training | 100 steps of Adam per round, lr 0.001 |
| What crosses the wire | Model weights, every round, both directions |

### Why accuracy and not AUC

The centres are small and unbalanced — one has a single class present in most batches, and
`roc_auc_score` raises on such a batch. That turns a metric into a crash partway through
validation. Accuracy degrades gracefully. This mirrors FLamby's own baseline code.

### Why every site validates on the pooled test set

`get_dataset` returns each site its own centre for training (`pooled=False`) and the
**pooled** test set for validation. Every site therefore scores against identical held-out
data, which is what makes per-round numbers comparable across sites. This is the
benchmark's own protocol.

It does mean every site holds the full test set. For 40 KB that costs nothing, but it is
worth stating explicitly in any write-up so it is not mistaken for a leak.

### Site-to-centre assignment

Each participating site trains on one FLamby centre, assigned in `federation.yaml`. The
dataset has four; `num_clients: 4` is the dataset's natural split size, not the number of
participating sites. Three sites drawing from a four-centre split is fine.

Past four sites, either two share a centre — defensible, but it must be disclosed — or the
experiment moves to a dataset with more centres. Fed-TCGA-BRCA has six and is still tiny
and tabular, which makes it the natural scale-up.

### Weighting

`client_weights_mode` selects how sites are weighted in the average:

- `equal` — every site counts the same
- `sample_size` — sites weighted by cohort size

The centres are quite unbalanced (199 vs 30 training samples), so `sample_size` is usually
the more defensible choice for a write-up. It costs one extra round trip to collect the
counts. **Record which mode a run used** — it is not recoverable from the results.

### Synchronous versus asynchronous

Synchronous FedAvg proceeds at the speed of the slowest site, every round. With sites in
different time zones and different queue policies, that is the likely limiting factor of
a first multi-site run rather than anything about the model.

APPFL ships asynchronous aggregation for exactly this situation. Switching is a config
change; it is arguably also the more interesting result, since it characterizes what
intercontinental federated learning actually costs.

---

## Record of runs

Runs of record go here: date, configuration, git commit, results, and how to reproduce.
This is what makes the repository citable.

### No completed baseline yet

**As of this repository's creation, no full multi-site run of this experiment has
completed.** The pre-migration tree contains a single output artifact from a two-site
attempt, and it holds one logging header line and nothing else — no rounds, no weights.

That is recorded here rather than omitted, because it bears directly on what the migration
could verify. The FLamby migration was checked for **structural** parity — the new code
emits equivalent client and server configurations and task payloads for the same inputs —
and **not** for numerical parity, because there is no baseline result to compare against.

Do not treat a first green run from this repository as evidence of parity with the
pre-migration tree. It would be evidence that this repository works, which is a different
and weaker claim.

### What was verified before the first run

- The partner-side path was proven on the pre-migration tree: a round-trip probe against a
  partner's multi-user endpoint returned a result executing as the expected service
  account, confirming the whole authorization chain.
- Structural parity of generated configurations.
- A loopback federation completes end to end on a single machine.

### Template for a run of record

```
### YYYY-MM-DD — <what this run was>

| | |
| --- | --- |
| Sites | N |
| Rounds | N |
| Weighting | equal / sample_size |
| Suite commit | <sha> |
| Config | <variant> |

Results:
| Site | Local-only accuracy | Federated accuracy |
| --- | --- | --- |

Pooled baseline: ...
FLamby published baseline: ...

Reproduce with:
    appfl-bio-suite run flamby-heart-disease --config <variant>
```

---

## References

- FLamby: <https://github.com/owkin/FLamby> · <https://arxiv.org/abs/2210.04620>
- UCI Heart Disease: <https://archive.ics.uci.edu/dataset/45/heart+disease> (CC BY 4.0)

See also [DATA.md](DATA.md) for how to obtain the data and [RUNBOOK.md](RUNBOOK.md) for
how to run it.
