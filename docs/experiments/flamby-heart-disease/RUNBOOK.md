# Fed-Heart-Disease — runbook

## Before you start

- Every partner has completed their setup and sent an endpoint UUID.
- `local/federation.yaml` has each site's UUID, centre, and output directory.
- You are on a machine that can stay up for the run, in `tmux` or `screen`.

## Preflight

```bash
appfl-bio-suite preflight --experiment flamby-heart-disease
```

Checks the environment, the version pins, that every shipped config parses, and that every
endpoint reports online. **Hard failures block a launch; warnings do not** — a warning is
something that is wrong on some clusters and fine on others.

Then verify each partner for real:

```bash
appfl-bio-suite endpoint smoke --experiment flamby-heart-disease
```

The `user` in each result must be the experiment's service account. If it is the account
that started the endpoint, identity mapping is not taking effect — that endpoint would
accept your tasks for the wrong reason and run them with the wrong permissions. Fix it
before running.

`endpoint status` alone is not sufficient. It tells you a daemon is connected and nothing
more.

## Dry run first

```bash
appfl-bio-suite run flamby-heart-disease --dry-run
```

Writes the resolved server and client configs without launching. Worth reading once per
federation change — it is the fastest way to see exactly what each partner's endpoint will
be sent, with every path resolved.

Check: each site's `client_id` maps to the intended `centre`, and no two sites share one.

## Launch

```bash
appfl-bio-suite run flamby-heart-disease
```

Start with a small number of rounds. Multi-round training exercises the transport
`rounds x sites` times, so anything flaky about an intercontinental link shows up
repeatedly rather than never. Prove a run completes before asking for twenty.

## What a healthy run looks like

```
[setup] 3 client(s), 5 global round(s)
[train] dispatching round 1 to all clients
[train] Site1 returned round 1/5
[train] sent updated global model to Site1
...
Federated Learning Training Completed!
```

Signs it is healthy:

- Every site returns each round. A site that returns round 1 and then goes quiet is
  usually a worker block being torn down and re-queued.
- Round times are roughly consistent. A round that takes much longer than its predecessor
  usually means a queue wait, not a compute problem.
- Validation accuracy moves in the first few rounds.

## When it is not healthy

**A site never returns round 1.** Its scheduler has not started the block. Check
`endpoint status`, then ask them to check their queue. This is normal for a few minutes.

**A site returns round 1, then rounds slow dramatically.** The block is being torn down
between rounds. Check `idle_heartbeats_soft: 0` and `idle_heartbeats_hard: 5760` in their
endpoint template.

**Everything stalls at the same round.** Synchronous FedAvg proceeds at the speed of the
slowest site. Check whether one site is queueing.

**A `422` or `403` on dispatch.** Identity mapping. See
[troubleshooting](../../coordinator/troubleshooting.md#identity-mapping) — the two have
completely different causes, and a 403 is not a mapping-expression problem.

**A deserialization error mid-run.** Version skew. Compare the Globus Compute stack and
Python minor version at every site. Note that a patch-level warning on every submit is
benign and expected.

## Reading the output

Results land in the server's `logging_output_dirname`. Per-site logs land in each
partner's `output_dir`, on their cluster.

For a run of record, capture: the number of sites, rounds, weighting mode, the suite
commit, per-site local-only accuracy, the federated accuracy, and the pooled baseline.
Add it to [ABOUT.md](ABOUT.md#record-of-runs) — `client_weights_mode` in particular is not
recoverable from the results afterwards.

## Validating without partners

```bash
appfl-bio-suite run flamby-heart-disease --config loopback --driver serial
```

Runs the whole federation on one machine with simulated sites. Proves the install, the
configs, the shipped modules, the trainer and the aggregator. It does **not** exercise the
transport — serialization, identity mapping and scheduler behaviour are exactly what a
loopback run cannot test.
