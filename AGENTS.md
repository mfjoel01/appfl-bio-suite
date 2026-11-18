# Repository working preferences

The user explicitly requests regular commits so another agent can continue if a
session or credit budget ends. Commit meaningful, tested checkpoints, including
specifications and recovery instructions before long experiments. Commit only
task-owned changes; do not sweep unrelated staged changes into a commit.

For the ongoing ALMA3 work, read:
- `docs/experiments/methylation/ALMA3-HANDOFF.md` for current status and exact resume commands.
- `docs/experiments/methylation/ALMA3-SWEEP.md` for the presentation study and controls.

Before starting/resuming an experiment, check for an existing live process and
its progress files; do not launch duplicate runs. Keep weights/checkpoints and
patient-level data out of git. Commit compact aggregate results and figures.
Preserve failed/negative results and distinguish synthetic validation from
clinical evidence. Do not alter an active run's hashed training code or spec.
