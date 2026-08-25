# Reference deployment

A complete, working deployment — **as an example, not a requirement**.

Everything here is specific to the machines this suite was built and first run on. Your
coordinator node is a different machine, possibly with a different scheduler, possibly a
cloud VM with no scheduler at all. Nothing in this document is encoded as a preflight
failure, and none of it is assumed anywhere in the code.

Read it for the *shape* of the decisions, and for the class of trap to look for on your
own systems.

---

## The deployment

| | |
| --- | --- |
| Coordinator node | An HPC login node, PBS Pro scheduler |
| Coordinator endpoint mode | Single-user |
| Python | 3.12.11 (partners on 3.12.13 — patch differences are fine) |
| Environment | conda, activated identically by the driver and by `worker_init` |
| Partner sites | Three continents, mixed SLURM and PBS, all multi-user endpoints |

---

## A working single-user endpoint config

This is a coordinator's own endpoint — no identity mapping, because it only ever accepts
its own tasks.

```yaml
display_name: my-endpoint
engine:
  type: GlobusComputeEngine
  max_workers_per_node: 1

  # Intra-cluster traffic only, so the curvezmq auth thread is skipped. Do not copy this
  # if your workers are reachable from outside your cluster.
  encrypted: false

  # The engine could not reliably auto-detect which NIC faces the compute nodes.
  # The interface name is cluster-specific -- on a sibling cluster sharing the same
  # filesystem, the equivalent interface was named differently and a copied config
  # failed with "Errno 99 Cannot assign requested address".
  address:
    type: address_by_interface
    ifname: bond0

  # The 4.9 schema wants the string form. A dict here fails validation.
  strategy: simple

  provider:
    type: PBSProProvider
    account: MY-ALLOCATION
    queue: debug
    walltime: 00:30:00
    cpus_per_node: 32
    nodes_per_block: 1
    init_blocks: 0
    min_blocks: 0

    # The queue has a per-user limit of one queued job, and parsl allocates climbing
    # block IDs even at max_blocks: 1. Holding any other job in that queue makes every
    # dispatch fail.
    max_blocks: 1

    select_options: ngpus=4
    scheduler_options: "#PBS -l filesystems=home:grand:eagle"
    launcher:
      type: MpiExecLauncher
      bind_cmd: --cpu-bind
      overrides: --depth=64 --ppn 1

    worker_init: |
      module use /soft/modulefiles
      module load conda
      source $(conda info --base)/etc/profile.d/conda.sh
      conda activate appfl_env
      # Scheduler $TMPDIR overflows the ~108-char Unix socket limit; workers then
      # never register.
      export TMPDIR=/tmp
      # numpy's BLAS starts one thread per core on import; on a busy login node that
      # exhausts the per-user thread budget.
      export OPENBLAS_NUM_THREADS=1
      export OMP_NUM_THREADS=1
      export MKL_NUM_THREADS=1
      export PYTHONUNBUFFERED=1
      echo "worker up: $(python --version) on $(hostname)"

# Keep the block warm between rounds. Without these it is torn down and re-queued each
# round, which on a busy queue turns minutes into hours.
idle_heartbeats_soft: 0
idle_heartbeats_hard: 5760
```

The four `worker_init` exports and the two `idle_heartbeats` settings are the parts worth
copying verbatim to any cluster. Everything above them is local detail.

## Bringing it up

```bash
module use /soft/modulefiles && module load conda
conda activate appfl_env
globus-compute-endpoint version            # must match the federation pin

# A pidfile on a shared filesystem records a host-local PID, so one left by another
# login node blocks startup here for no real reason.
rm -f ~/.globus_compute/my-endpoint/daemon.pid

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  globus-compute-endpoint start my-endpoint

appfl-bio-suite endpoint smoke
```

Healthy log lines: `Awaiting messages from queue`, `Scaling out by 1 blocks`,
`Allocated block ID 0`. Climbing block IDs or repeated `qsub: would exceed` means the
queue problem above.

---

## Traps this deployment hit

Listed as **shapes of problem**, because yours will differ in detail.

### Two clusters sharing a filesystem are not interchangeable

Two clusters shared `/home` and a project filesystem but had separate schedulers, separate
process tables, and different hardware. Consequences:

- The environment's torch was a Cray-MPI build that could not import on the sibling
  cluster at all (`OSError: libmpi_gtl_cuda.so.0`). The driver had to run on one specific
  cluster, and a shared `$HOME` made the wrong one look correct until the import failed.
- `qstat` and `ps` on one said nothing about the other, so an endpoint daemon started on
  the wrong one submitted jobs to a scheduler nobody was watching.
- The network interface had a different name on each.
- There was no non-interactive ssh between them, so nothing could be scripted across.

**The general shape:** a shared filesystem creates the illusion of one machine. If your
site has anything similar, set `coordinator.host_check` in `federation.yaml`.

### A second Python installation shadowing the first

`~/.local/lib/python3.12/site-packages` held a generic torch build that imported anywhere.
Invoking the environment's Python directly, without activating it properly, silently
picked up that copy — while workers, whose `worker_init` activated the environment (which
sets `PYTHONNOUSERSITE`), used the environment's build.

**`pip list` reported the same version in both cases.** The only way to see it was
`torch.__file__`.

**The general shape:** version reporting is not version truth. `preflight --check env`
looks for this.

### The queue that looked more generous was worse

One queue had a per-user limit of a single queued job but a dedicated node pool and
scheduled in seconds. Another had a higher job cap but left jobs pending indefinitely on a
resource tag. A third could not host the endpoint at all because it required a minimum
node count far above the one node an endpoint needs.

**The general shape:** pick the queue that schedules quickly at one node, not the one with
the most generous limits.

---

## What was NOT needed

Worth recording, because each was considered and would have added real complexity:

- **A VPN or any inbound network access.** Endpoints connect outward to the Globus Compute
  relay. Nothing listens on the coordinator.
- **A shared filesystem between institutions.** Data never moves.
- **A shared credential.** The coordinator's identity string is not a secret, and it is
  the only thing distributed.
- **Matching operating systems, schedulers, or hardware.** Only the Python minor version
  and the Globus Compute stack need to match.
- **A coordinator-hosted multi-user endpoint.** It would map submitters to accounts on the
  coordinator's cluster, which cannot grant access to anything on a partner's — the wrong
  direction entirely.
