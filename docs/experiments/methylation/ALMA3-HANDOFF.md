# ALMA3 continuation log

## User preferences

Commit completed work and meaningful checkpoints regularly. Keep this log current
so another agent can resume if the session ends. Do not leave long experiments
without a committed configuration, run instructions, and output locations.

## Current state

Branch: `feat/alma3-federated-finetune` in `appfl-bio-suite-alma3`.
The original implementation and five-site validation are complete. See ALMA3.md.
Member `CV:e11`. Commands below use three machine-specific settings:
- `ALMA3_RELEASE`: the extracted, unmodified v3.1.2 release (see `scripts/methylation/fetch_alma3_release.py`).
- `ALMA3_PYTHON`: the isolated worker environment (Python 3.12, Torch 2.14.1, ALMA3 3.1.2; see ALMA3.md).
- `COORDINATOR_PYTHON`: the suite's own environment with APPFL installed.
Paths and hosts in the committed run records use the same placeholders.
Run from this worktree with `PYTHONPATH=src NUMPY_MADVISE_HUGEPAGE=0 OPENBLAS_NUM_THREADS=1`.
Weights/checkpoints stay outside git; licenses accompany derivatives.
Latest completed output: `local/output/alma3-skew-sweep`.
Sixteen ALMA3 tests pass; lint/format checks pass. Real test cohorts untouched.

## Active request

Create a presentation-ready sweep for St. Jude: five sites interpolate from
extreme class-ratio imbalance to identical class ratios. Show whether local
models catch up with federation. The existing extreme point does NOT show
federation beating every local: local-S2 is 99.98% versus FedAvg 100%. Treat the
proposed curve as a hypothesis; show all locals and the best local, not only
an average that hides this result. Keep patient counts, optimization budgets,
and evaluation distributions fixed across the sweep. Preserve negative results.

## Completion (2026-11-18)

The requested sweep and presentation bundle are complete. No training remains.
All 21 points / 168 model records completed; integrity audit passed on 1,680
evaluation records and 150 unique final checkpoints. See ALMA3-SWEEP.md for the
measured results and presentation links. Full output: `local/output/alma3-skew-sweep`.
Compact results, figures, specifications and checkpoint hashes are committed.
Weights and patient-level data remain outside git.

FedAvg scores 100% throughout; extreme mean/best local scores are 72.72/80.53%.
Best local catches up at convergence 0.10; all locals reach 100% at the identical
ratio endpoint. Do not imply the advantage lasts until site ratios are identical.
Scoped methylation lint/format checks pass. Repository-wide lint also reports
four pre-existing errors in unrelated fine-mapping scripts.
Sixteen tests passed; real-worker gradient parity and no-retraining resume checks
passed. Local class-balanced reweighting/resampling was not evaluated.

To audit and rerender existing results without training:

```bash
python scripts/methylation/audit_alma3_sweep.py local/output/alma3-skew-sweep
NUMPY_MADVISE_HUGEPAGE=0 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src \
$COORDINATOR_PYTHON -m appfl_bio_suite.experiments.methylation.alma3.sweep_report local/output/alma3-skew-sweep
```

The following sections preserve historical progress and recovery commands;
references to running jobs or pending checks below are superseded by this status.

## Sweep protocol (implementation in progress)

Committed specification: `src/appfl_bio_suite/experiments/methylation/configs/alma3-skew-sweep.json`.
Seven convergence levels 0, .05, .10, .25, .5, .75, 1; three data seeds 42–44.
Extreme rows: [97,1,1,1], [2,96,1,1], [1,1,96,2], [1,1,1,97], [49,1,1,49].
Identical endpoint: [30,20,20,30] at every site. Every row has 100 patients;
pool totals remain [150,100,100,150] at every point. Integer rounding must
preserve both row and column totals and at least one patient per class.
Evaluation: fixed, class-balanced validation (25/class/site), with common
patient reservoirs paired across sweep points. Preserve site measurement shifts
throughout; only class ratios converge. No test access, DP noise, or clinical
calibration claim in this sweep. Compare zero-shot, pooled, FedAvg, and all
five locals. Report the best local as well as their mean. Use development seed
41 for a mechanics check, then locked seeds 42–44; do not tune toward a desired
winner. Checkpoint/commit compact per-point scores as experiments progress.

## Implementation checkpoint

Sweep runner: `python -m appfl_bio_suite.experiments.methylation.alma3.sweep`.
Resume requires an identical specification, worker identity and training-code
hashes. Each arm saves atomic round progress plus a hashed checkpoint. Finished
arms are reused on resume; initial checkpoints are not rewritten. Completed
points publish compact scores and optionally create git commits automatically.

Sixteen tests passed before final resume/format cleanup. Real CV:e11 batched
per-patient clipping versus serial: maximum parameter error 7.45e-9 over three
steps on 12 patients. Batched gradients preserve per-patient clipping and only
sum afterward. No DP noise is used in this sweep.

Development command (seed 41, endpoints) and full run use the same module.
Full run/resume command (add `--resume` only when output already exists):

```bash
NUMPY_MADVISE_HUGEPAGE=0 OPENBLAS_NUM_THREADS=1 PYTHONPATH=src \
$COORDINATOR_PYTHON -m appfl_bio_suite.experiments.methylation.alma3.sweep \
 --config src/appfl_bio_suite/experiments/methylation/configs/alma3-skew-sweep.json \
 --out local/output/alma3-skew-sweep \
 --release $ALMA3_RELEASE \
 --worker-command $ALMA3_PYTHON \
 --publish docs/experiments/methylation/results/alma3-skew-sweep --commit-progress
```

Sweep checkpoint: seed 41, convergence 0 complete. Scores in `docs/experiments/methylation/results/alma3-skew-sweep-development`; resume command below remains valid.

Sweep checkpoint: seed 41, convergence 1 complete. Scores in `docs/experiments/methylation/results/alma3-skew-sweep-development`; resume command below remains valid.

Sweep checkpoint: seed 41, convergence 0 complete. Scores in `docs/experiments/methylation/results/alma3-skew-sweep-development`; resume command below remains valid.

Sweep checkpoint: seed 41, convergence 1 complete. Scores in `docs/experiments/methylation/results/alma3-skew-sweep-development`; resume command below remains valid.

## Development result and parallel production runner

Development seed 41 is complete at both endpoints: FedAvg 100% throughout;
mean local 72.52% -> 99.96%; best local 79.80% -> 100%. No training settings were
changed after inspecting these results. Scores were automatically committed.
A resume integrity check is in progress; compare hashes AND mtimes against
`local/alma3-sweep-before-resume.json` when it finishes.

The full run will use three independent CPU seed workers under one controller,
which alone writes progress commits (avoiding concurrent git writes). Each
worker caps Torch/BLAS thread counts. There are no GPU or multi-member jobs.
Use this command for the full parallel run; add `--resume` if its output exists:

```bash
NUMPY_MADVISE_HUGEPAGE=0 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src \
$COORDINATOR_PYTHON -m appfl_bio_suite.experiments.methylation.alma3.sweep_parallel \
 --config src/appfl_bio_suite/experiments/methylation/configs/alma3-skew-sweep.json \
 --out local/output/alma3-skew-sweep \
 --release $ALMA3_RELEASE \
 --worker-command $ALMA3_PYTHON \
 --publish docs/experiments/methylation/results/alma3-skew-sweep
```

Check `execution.json`, `progress.json`, and `seed-42.log` / `seed-43.log` /
`seed-44.log` under the output directory. `COMPLETE.json` is written only after
all workers and combined figures succeed. A failed worker leaves resumable
round checkpoints; the controller records `FAILED.json`. Verify old processes
have stopped before resuming. The sequential command above remains available,
but do NOT use it against a parallel-run output directory.

Resume verification passed: all 141 checkpoint hashes and mtimes unchanged; no retraining. Full parallel sweep is ready to launch.

## Full sweep running

Launched with implementation commit `18272de` on a cluster login node.
Initial controller PID 2954105; coordinator seed PIDs 2954454/2954455/2954456.
These PIDs may become stale: inspect matching command lines and progress before
resuming. Output `local/output/alma3-skew-sweep`; 21 seed/level points expected.
The controller commits compact results automatically. All code and specifications
were committed before launch. Do not modify hashed training modules mid-run.

When complete: inspect the combined report and figures; verify 168 model records
(3 seeds × 7 levels × 8 arms), fixed class-balanced evaluation, pooled training
counts and same zero-shot scores across levels. Copy final figures/summary CSVs
into docs, update ALMA3-SWEEP.md with measured results, and commit. Preserve the
best-local comparator and explicitly label this as synthetic validation.

Parallel sweep checkpoint: 2 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 3 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 4 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 6 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 7 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 8 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 9 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 10 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 12 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 13 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 15 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 16 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 17 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 18 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 19 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 20 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

Parallel sweep checkpoint: 21 completed seed/level points. Durable scores: `docs/experiments/methylation/results/alma3-skew-sweep`.

## Platform extension (2026-11-18)

User clarified the intended application: hospital-specific nanopore adaptation
plus differing class ratios. Added optional read-count/calling-error/purity
simulation, a separate rotation-free platform sweep config and ALMA3-PLATFORM.md.
The original completed sweep is preserved. New platform training has NOT run;
no accuracy or adaptation gain is claimed for this extension. Run to a new
output directory; old results cannot resume against modified training hashes.
