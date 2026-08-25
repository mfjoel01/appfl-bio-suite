# Architecture

## Why this is a downstream package, not a fork of APPFL

Three options were considered.

**Use APPFL `main` as-is.** Rejected. It ships examples that demonstrate the API but none
of the scaffolding a partner needs to actually run an experiment — no onboarding
documents, no bundle generation, no identity validation, no preflight. Not something you
can hand to a collaborator.

**Fork APPFL.** Rejected. A perpetual merge burden to track upstream, and it bloats a
project that is not ours with experiment-specific material that does not belong there.

**A standalone package depending on `appfl` from PyPI.** Chosen. No fork maintenance, no
upstream pull requests needed to ship an experiment, no bloat added to APPFL.

### Why that is possible

APPFL already exposes the extension points a fork would be needed for. A trainer resolves
three ways — `trainer_path` (a filesystem path), `trainer_source` (source text shipped
over the wire), or a bare name looked up in `appfl.algorithm.trainer`. Datasets, models,
losses, metrics and aggregators use the same mechanism.

So custom trainers, models, loaders and aggregators live entirely in this package and are
referenced by path. Nothing needs to be merged upstream for an experiment here to run.

### One honest caveat

The value of depending on APPFL rather than forking it is **not** "automatically stays up
to date". That would reintroduce exactly the version-skew problem this suite exists to
prevent. The value is that upgrading becomes a one-line, deliberately tested change
instead of a merge.

And it is not entirely free today: `appfl==1.10.0` on PyPI imports a `globus-compute-sdk`
module that the pinned 4.9.0 removed, so a clean install of the pinned set cannot import
the production dispatch path without a narrow workaround. See
`src/appfl_bio_suite/core/compat.py`, which is written to delete itself when upstream
fixes it.

---

## The trust model

```
   COORDINATOR                                          PARTNER
   single-user endpoint                                 multi-user endpoint (root)
   accepts only its own tasks                           maps ONE identity -> ONE account

   identity: you@example.org  ──── dispatch ────►  match: you@example\.org
                                                   output: gwas_svc
                                                        │
                                                        ▼
                                                   task runs as gwas_svc
                                                   on their compute node,
                                                   against their local data
```

**What crosses the network:** the task (a function and its arguments), and the result.
For GWAS the result is per-variant summary statistics. For FLamby it is model weights.
Raw data never moves.

**What is shared in advance:** the coordinator's Globus identity string. That is all. It
is not a credential — knowing it grants nothing, because the partner's own endpoint
configuration is what decides to honour it.

**What each side controls:** the partner controls whether to accept you at all, which
local account your tasks run as, and what that account can read. Revoking access is a
one-line edit to a file their endpoint polls. The coordinator controls what gets
dispatched.

An earlier design used a shared confidential client — one credential distributed to every
site. It worked, but anyone holding it could submit arbitrary functions to any site, and
everything ran under a single machine identity with no attribution. The multi-user
endpoint model needs no shared secret and gives per-identity attribution.

---

## Where code runs, and why it matters

This is the distinction that causes the most confusion, so it is worth stating plainly.

| | Read on | Executed on | May import |
| --- | --- | --- | --- |
| `dataset.py` | coordinator | **partner worker** | stdlib + the experiment's extra |
| `trainer.py` | coordinator | **partner worker** | stdlib + the experiment's extra |
| `model.py` | coordinator | **partner worker** | stdlib + the experiment's extra |
| `aggregator.py` | coordinator | coordinator | anything |
| `plotting.py` | coordinator | coordinator | anything |
| `simulation/` | coordinator | coordinator | anything |

APPFL reads the first three on the driver, inlines their **source text**, and ships the
text. The worker executes it standalone, in a temporary directory, with whatever the
partner installed.

Two consequences:

1. **A shipped module may not import from `appfl_bio_suite`.** A partner installs the
   experiment's extra, not this package. Enforced by `tests/test_shipped_modules.py`,
   which walks the AST rather than trusting a comment.

2. **`dataset_path` and `trainer_path` are resolved on *your* machine**, for every client,
   including partners on other continents. A partner is never asked for one. The opposite
   holds for `data_dir` and `output_dir`, which are used on the worker and must be
   absolute paths on their cluster.

The original GWAS trainer violated the first rule by importing two sibling modules. That
is the entire reason every partner needed a `PYTHONPATH` line in their endpoint
configuration, and a missing one was the largest single source of partner-side breakage.
Making shipped modules self-contained removed that line from the setup guide.

---

## Simulation is coordinator-side only

Everything else in this suite eventually reaches a partner cluster. Simulation never does.
It runs on your hardware and produces per-site bundles that are then distributed.

Two consequences, both enforced by packaging rather than by documentation:

- The shipped-module import rules do not apply to simulation code. It may import freely
  and use heavy dependencies.
- Its dependencies must never land on a partner. They live in a separate `gwas-sim`
  extra; a partner installs `appfl-bio-suite[gwas]` and gets none of them.

Not every experiment simulates. FLamby partners download a public dataset, so that
experiment has no simulation stage at all — and the pipeline tolerates that rather than
requiring an empty stub to satisfy an interface.

---

## One config file

`federation.yaml` declares the coordinator identity, the participating sites, their
endpoint UUIDs and schedulers, and per-site data assignments. Every generator, config,
document template and check reads from it.

Nothing in `src/` contains a coordinator identity, an endpoint UUID, or a site name.
That is checked mechanically by `scripts/check_reusability.py`, which runs in CI.

The consequence is the property the whole suite is organized around: standing up a new
federation is filling in one file. If describing your federation ever requires editing
code, that is a bug.

---

## Layout

```
src/appfl_bio_suite/
├── cli.py                  one entry point for every operational task
├── core/                   experiment-agnostic
│   ├── config.py           federation.yaml: load, validate, cross-reference
│   ├── experiments.py      the registry -- the only place that knows what exists
│   ├── identity.py         identity-mapping validation
│   ├── endpoint.py         status and round-trip smoke tests
│   ├── preflight.py        environment, pins, configs, endpoints
│   ├── launch.py           federation.yaml -> APPFL configs; run
│   ├── partner.py          bundle generation
│   ├── simulation.py       pipeline contract, run manifests, provenance
│   ├── compat.py           one upstream workaround, written to be deleted
│   └── drivers/            globus_compute (production), serial (loopback)
└── experiments/
    ├── flamby_heart_disease/
    ├── gwas/
    └── fine_mapping/       fine-mapping; also fedfm/ (vendored verbatim upstream)
```

Adding an experiment means writing a registry entry in `core/experiments.py`, creating the
package it names, and writing four documents. The CLI, preflight, bundle generator and
test suite pick it up with no further changes. `fine_mapping` exists, unimplemented,
specifically to keep that claim honest.
