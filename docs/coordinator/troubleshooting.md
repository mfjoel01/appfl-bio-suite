# Troubleshooting

Every entry here cost real debugging time, some of it days. They are written with the
**mechanism**, not just the symptom, because in almost every case the symptom pointed
somewhere other than the cause.

Where a problem can be caught by a check instead of by reading, it is — the relevant
command is named. Documentation that has to be read is the weakest control available.

---

## Quick index by symptom

| You see | Go to |
| --- | --- |
| `422 … Identity failed to map to a local user name` | [Identity mapping](#identity-mapping) |
| `403 ENDPOINT_ACCESS_FORBIDDEN` | [Unprivileged endpoint start](#an-unprivileged-start-silently-disables-identity-mapping) |
| `404` on submit | [UUIDs are per-account](#endpoint-uuids-are-per-account) |
| `Invalid \ expression` when loading a mapping | [Escapes in `match`](#only-a-few-escapes-are-legal-in-match) |
| "workers failed to register" / "AF_UNIX path too long" | [TMPDIR](#tmpdir-overflows-the-unix-socket-path-limit) |
| `RuntimeError: can't start new thread` | [Thread fan-out](#numpys-thread-fan-out-exhausts-the-process-limit) |
| Workers start but never connect back | [NIC detection](#workers-start-but-never-connect-back) |
| First round fine, later rounds crawl | [Block teardown](#worker-blocks-are-torn-down-between-rounds) |
| `Configuration … contains an 'engine' field` | [Version 4.12](#globus-compute-endpoint-412-rejects-the-single-user-config-layout) |
| `SystemError: unknown opcode` | [Python minor mismatch](#python-minor-version-mismatch) |
| `ModuleNotFoundError: globus_compute_sdk.sdk.login_manager` | [The appfl/SDK conflict](#appfl-1100-and-globus-compute-sdk-490-conflict) |
| A pinned package moved on its own | [Utilities move pins](#installing-utilities-moves-pinned-packages) |
| Same version reported, different behaviour | [`~/.local` shadowing](#user-site-packages-shadow-the-environment) |
| `ModuleNotFoundError` for your own helper module | [Shipped module imports](#shipped-modules-cannot-import-siblings) |
| "Another instance of this endpoint is running" | [Stale pidfile](#stale-pidfile-on-a-shared-filesystem) |
| Endpoint reports online but tasks never run | [Status is not liveness](#online-is-not-the-same-as-working) |

---

## Identity mapping

### `^` and `$` anchors silently break every submission

**Symptom.** Every task is rejected with:

```
422 SEMANTICALLY_INVALID — Request payload failed validation:
Identity failed to map to a local user name. (LookupError)
```

**Mechanism.** The `match` field in a multi-user endpoint's identity-mapping file is
**not** a plain regular expression. `ExpressionIdentityMapping._compile_match` escapes any
`^ $ + { } [ ]` you supply into *literal characters*, then wraps the whole pattern in its
own `"^" + match + "$"`.

So `"^me@example\\.org$"` compiles to `^\^me@example\.org\$$` — a pattern that matches
only a string literally beginning with a caret. It can never match, with either `source`
field.

```
match written          compiles to                    matches?
^me@example\.org$      ^\^me@example\.org\$$          never
me@example\.org        ^me@example\.org$              yes
```

**Fix.** Drop the anchors. Keep the `\.`. Dropping them does not loosen the match — the
mapper anchors it for you.

**Catch it in a second, offline:**

```bash
appfl-bio-suite identity validate <mapping.json> --experiment <name>
```

**Cost when this was live:** over a week across two partner sites, plus one confident
wrong diagnosis (`{email}` vs `{username}`) that survived a full partner round trip and
changed nothing, because the anchors were the real cause with either field.

### Only a few escapes are legal in `match`

**Symptom.** `InvalidMappingError: … Invalid \ expression: \-` and the endpoint refuses
the entire mapping document.

**Mechanism.** The mapper accepts only `\.` `\?` `\*` `\|` `\(` `\)` `\\` and `\0`–`\9`.
Anything else after a backslash is a hard error.

This bites anyone who builds the pattern with Python's `re.escape()`, which escapes `-` —
so an identity like `a.researcher@example-university.edu` becomes
`a\.researcher@example\-university\.edu` and the file becomes invalid.

**Fix.** Escape only the dots. The suite's bundle generator does this correctly; if you
are hand-writing a mapping, do not use `re.escape`.

### `{username}` vs `{email}` is not the problem

Both work once the anchors are gone. `{username}` is what Globus Compute's own shipped
example uses and what has been verified end to end, so the suite generates that — but if
someone tells you the field name is the cause, they are repeating a diagnosis that was
investigated and found wrong.

### Records without a `status` are silently skipped

`map_identities` ignores any identity record whose `status` is not `used` or `private`.
A hand-built test record without one matches nothing, which looks exactly like a bad
regex and has produced a false negative mid-debug. The suite's validator sets it.

### The mapping file is polled, not uploaded

`PosixIdentityMapper` runs a thread that stats the file every few seconds and reloads on
change. **A restart is not needed** to apply a mapping edit.

An earlier version of our own guide said the config was uploaded to Globus at endpoint
start and required a restart. That is wrong at the mechanism level: the file's *contents*
are never sent anywhere. Only the path is stored, locally. A restart *is* still needed to
change **who runs** the endpoint, which is unrelated.

Consequence worth internalizing: **a 422 does not prove the fault is service-side.**
Mapping is evaluated on the partner's endpoint.

### An unprivileged start silently disables identity mapping

**Symptom.** `403 ENDPOINT_ACCESS_FORBIDDEN`. The endpoint looks completely healthy.

**Mechanism.** Identity mapping is only honoured when the endpoint process is privileged.
Started as an ordinary user, the endpoint logs a warning, **starts successfully anyway**,
ignores `identity_mapping_config_path` entirely, and accepts only the identity that
started it.

The silent success is what makes this expensive. "It started fine" is not evidence that
the mapping took effect.

**Fix.** Restart as root — then re-read the UUID, because it will have changed.

### Endpoint UUIDs are per-account

The ID lives in `<that user's home>/.globus_compute/<name>/endpoint.json`. An endpoint
configured under a personal account and one configured under `root` are **different
endpoints with different UUIDs**. After switching to root, any UUID sent earlier is dead.

Also: `stop`/`start`/`restart` preserve the UUID; `delete` followed by `configure` issues
a new one.

### Why the endpoint must be on the partner's side

A multi-user endpoint maps submitters to local accounts **on the machine it runs on**. A
coordinator-hosted one would map partners to accounts on the coordinator's cluster, which
cannot grant access to anything on theirs — the wrong direction entirely.

It also needs root, which a coordinator typically does not have on a shared HPC system,
while each partner's administrator does have it on their own.

---

## Version alignment

### `globus-compute-endpoint` 4.12 rejects the single-user config layout

**Symptom.** `Configuration … contains an 'engine' field; endpoint will not start.`

**Mechanism.** 4.12 redesigned `start` to always launch the multi-user manager, which
rejects the single-user `engine`-in-config layout that a coordinator endpoint uses.

**Fix.** Pin `==4.9.0`. Never install with `>=4.7.0`, which is what APPFL's own `setup.py`
declares and which resolves to 4.12 today.

```bash
appfl-bio-suite preflight --check pins
```

### Python minor-version mismatch

**Symptom.** `SystemError: unknown opcode`, before any of your code runs.

**Mechanism.** Dill ships function bytecode. Bytecode is portable within a minor version
and not across one, so a 3.12 driver against 3.10 workers fails at decode.

**Fix.** Standardize the federation on one minor version.

**Patch levels do NOT need to match, and you should not make partners chase them.** A
3.12.11 driver against 3.12.13 workers prints:

```
    SDK: Python 3.12.11/Dill 0.3.9
Workers: Python 3.12.13/Dill 0.3.9
This may cause serialization issues.
```

That warning is benign here: the SDK compares full version strings, but CPython's bytecode
magic number is frozen across a minor series (verified: `MAGIC_NUMBER` = 3531 on 3.12.11,
the 3.12 value) and dill matches. It has been observed live against a working partner
endpoint with no ill effect.

An earlier version of our guide overstated this as "fails before your code ever runs".
That is true of a *minor* mismatch and false of a patch one. Keep patch skew as the first
suspect for a real deserialization error involving model payloads — but do not preemptively
force alignment.

### appfl 1.10.0 and globus-compute-sdk 4.9.0 conflict

**Symptom.** `ModuleNotFoundError: No module named 'globus_compute_sdk.sdk.login_manager'`
when importing APPFL's Globus Compute communicator.

**Mechanism.** `appfl==1.10.0` as published on PyPI imports that module at module scope
(`globus_compute_server_communicator.py`, line 18). The 4.9 SDK line removed it. Both pins
are individually required — 4.9.0 for the reasons above, 1.10.0 because it is current — so
a clean install of the pinned set cannot import the production dispatch path.

The import is dead weight on our path: `AuthorizerLoginManager` is used only in the hosted
APPFLx token-auth branch, which a self-driving coordinator never reaches.

**Fix.** `appfl_bio_suite.core.compat` registers a placeholder before APPFL is imported.
It is a no-op wherever the problem is absent, and `tests/test_upstream_shim.py` **fails
once it is no longer needed**, so it gets deleted rather than rotting.

If you see this error, something imported `appfl` before `appfl_bio_suite`. Import the
suite first.

### Installing utilities moves pinned packages

**Symptom.** A working endpoint breaks on its next restart, having changed nothing.

**Mechanism.** The globus-compute stack is tightly pinned. Installing an unrelated tool
into the same environment lets pip resolve its dependencies against yours. Observed live:
`pip install globus-cli` upgraded `click` 8.1.8 → 8.4.2 and `globus-sdk` 4.4.1 → 4.8.0,
both violating `globus-compute-endpoint 4.9.0`'s requirements. The *running* daemon kept
working, because its modules were already imported. The next restart would not have.

**Fix.** Install one-off tools in a throwaway virtualenv, never in the environment that
runs the endpoint or the driver.

### User-site packages shadow the environment

**Symptom.** Two machines report the same version and behave differently. Or the driver
works and the workers do not, with no version difference visible anywhere.

**Mechanism.** A package in `~/.local/lib/pythonX.Y/site-packages` shadows the
environment's copy. Workers usually run with `PYTHONNOUSERSITE` set (activating a conda
environment sets it), so they use the environment's copy while an interactive driver
session uses `~/.local`'s — a silent version skew.

**`pip list` reports the same version either way.** It cannot reveal this.

**Fix.** Check the actual file:

```bash
python -c "import torch; print(torch.__version__, torch.__file__)"
```

Install with `PYTHONNOUSERSITE=1 python -m pip install --no-deps --ignore-installed <pkg>`.
`appfl-bio-suite preflight --check env` flags it.

---

## Workers and scheduling

### `TMPDIR` overflows the Unix socket path limit

**Symptom.** The block starts and logs normally, then the task fails with
`BadStateException: workers failed to register`. Block stderr shows
`OSError: AF_UNIX path too long`.

**Mechanism.** The parsl worker pool opens a `multiprocessing` SyncManager Unix socket
under `$TMPDIR`. Schedulers often set `$TMPDIR` to a long per-job path, and the socket
path limit is about 108 characters.

**Fix.** `export TMPDIR=/tmp` in `worker_init`. The suite's endpoint template includes it.

### numpy's thread fan-out exhausts the process limit

**Symptom.** `RuntimeError: can't start new thread`, or `Interchange failed to start`.
Sometimes the login node refuses even `echo` with
`fork: retry: Resource temporarily unavailable`.

**Mechanism.** OpenBLAS and OpenMP each start one thread per logical core on every numpy
import. On a 64-core shared login node with other users active, that exceeds the per-user
thread budget.

**Fix.** Export all three in `worker_init` **and** in the driver environment:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

The driver needs them too — it imports numpy as well, and it runs on the same busy login
node. `appfl-bio-suite run` sets them for the driver process automatically.

### Workers start but never connect back

**Symptom.** The block runs, workers appear to start, the interchange times out. A
config copied between two clusters can also fail outright with `Errno 99 Cannot assign
requested address`, when the interface is named differently on each.

**Mechanism.** The engine could not auto-detect which network interface faces the compute
nodes.

**Fix.** Add to the endpoint config, under `engine:`:

```yaml
address:
  type: address_by_interface
  ifname: <the interface the compute nodes use>
```

Set `interface:` on the site in `federation.yaml` and the generated template includes it.

### Worker blocks are torn down between rounds

**Symptom.** The first round completes promptly; later rounds take minutes each. A run
that should take minutes takes hours.

**Mechanism.** By default an idle block is released. A multi-round federated run sends one
task per site per round, so the block is torn down and re-queued every round — and each
re-queue costs a full scheduler wait.

**Fix.**

```yaml
idle_heartbeats_soft: 0
idle_heartbeats_hard: 5760
```

This matters far more for FLamby (20 rounds x N sites) than for GWAS (one task per site).

### Stale pidfile on a shared filesystem

**Symptom.** *"Another instance of this endpoint is running (perhaps on another login
node?)"* but the PID does not exist on this host.

**Mechanism.** `daemon.pid` lives on the shared filesystem and records a host-local PID.
Another login node's PID means nothing here.

**Fix.** Confirm the cloud status is `offline`, confirm `ps -p <PID>` is empty, then remove
the pidfile. Avoid it by always starting the daemon on the same host.

### "online" is not the same as "working"

`get_endpoint_status` reports that a daemon is connected. It does **not** tell you that
your identity is accepted, that mapping resolves, that the scheduler accepts the job, that
the environment activates, or that the worker can return a result. Every one of those has
failed independently while the endpoint reported `online`.

Use a real round trip:

```bash
appfl-bio-suite endpoint smoke <site> --experiment <name>
```

The `user` in the result is the load-bearing field: on a partner's multi-user endpoint it
must be the experiment's service account. If it is the account that started the endpoint,
mapping is not taking effect.

Also note the reverse: the **cloud** status is authoritative, not local
`globus-compute-endpoint list`, which reads a possibly-stale pidfile.

---

## Code shipped to workers

### Shipped modules cannot import siblings

**Symptom.** `ModuleNotFoundError` on the worker for a module that sits right next to the
one being shipped.

**Mechanism.** APPFL reads dataset and trainer files on the coordinating driver, inlines
their **source text**, and ships the text. The receiving worker executes it standalone in
a temporary working directory. `Path(__file__).parent` is meaningless there and siblings
do not exist.

The original GWAS trainer did `from gwas_config import ...`, which is the entire reason
every partner needed a `PYTHONPATH` entry in `worker_init` — and a missing one was the
largest single source of partner-side breakage on that project.

**Fix.** Shipped modules are self-contained. Enforced by
`tests/test_shipped_modules.py`, which walks the AST and rejects any import outside the
allowlist. That test is why the `PYTHONPATH` line no longer appears in the partner guide.

### `dataset_path` is resolved on the driver, never on the partner

`dataset_path` and `trainer_path` are read on **your** machine. Their contents are shipped.
A partner is never asked for either, and any document that asks for one is wrong — the
value would be silently unused.

The opposite holds for `data_dir`, `output_dir` and every logging path: those are used on
the worker and must be absolute paths on the partner's own cluster. A relative path there
resolves inside the endpoint's task working directory, which is not where anyone expects.

---

## FLamby specifically

### FLamby's `setup.py` is hostile — do not run it

`FLamby/setup.py` hooks `install`, `develop`, **and `egg_info`** to unconditionally run

```
pip install large-image-source-openslide --find-links https://girder.github.io/large_image_wheels
```

with no error checking. The `egg_info` hook means it fires on **mere metadata
generation** — so even a resolver inspecting the package triggers it. It reaches an
external index many clusters cannot access, and it can move pinned packages.

`make install` is worse: it builds a separate conda environment from FLamby's own unpinned
`environment.yml`, which the endpoint never activates.

**Fix.** Install via a `.pth` file pointing at the checkout, which bypasses `setup.py`
entirely and leaves every pin byte-identical. That is what the partner guide now says.

### `dataset_location.yaml` lives inside the installed package

FLamby records the downloaded data path inside its own package directory. Whoever runs the
download determines who can find the data afterwards.

**Fix.** Run the download **as the service account** the endpoint maps tasks to, and
verify with one import as that account before reporting the endpoint ready. Otherwise the
worker may read a pointer to a path it cannot reach, and fail late, mid-round.

---

## Traps specific to one deployment

The machines this suite was built on hit a handful of problems that are shapes of trap
rather than symptoms with a fix: two clusters sharing a filesystem but not a scheduler, a
second Python shadowing the first, and a queue that looked generous but never scheduled.
They are written up once, with the config that provoked them, in
[reference-deployment.md](reference-deployment.md#traps-this-deployment-hit).

If your deployment has a constraint of that kind, set `coordinator.host_check` in
`federation.yaml`. It warns by default; `host_check_enforce: true` makes it fatal.
