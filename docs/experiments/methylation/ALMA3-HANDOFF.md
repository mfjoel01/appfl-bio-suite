# ALMA3 continuation log

## User preferences

Commit completed work and meaningful checkpoints regularly. Keep this log current
so another agent can resume if the session ends. Do not leave long experiments
without a committed configuration, run instructions, and output locations.

## Current state

Branch: `feat/alma3-federated-finetune` in `appfl-bio-suite-alma3`.
The original implementation and five-site validation are complete. See ALMA3.md.
Real release: `../scratch/alma3-release-v3.1.2`, member `CV:e11`.
Worker Python: `/lus/grand/projects/GeomicVar/mjoel/BIOSIM_SUITE/scratch/alma3-runtime/bin/python`.
Coordinator Python: `/home/mfjoel01/.conda/envs/appfl_env/bin/python`.
Run from this worktree with `PYTHONPATH=src NUMPY_MADVISE_HUGEPAGE=0 OPENBLAS_NUM_THREADS=1`.
Weights/checkpoints stay outside git; licenses accompany derivatives.
Latest completed output: `local/output/alma3-five-sites-skewed-validation`.
Thirteen ALMA3 tests pass; lint/format checks pass. Real test cohorts untouched.

## Active request

Create a presentation-ready sweep for St. Jude: five sites interpolate from
extreme class-ratio imbalance to identical class ratios. Show whether local
models catch up with federation. The existing extreme point does NOT show
federation beating every local: local-S2 is 99.98% versus FedAvg 100%. Treat the
proposed curve as a hypothesis; show all locals and the best local, not only
an average that hides this result. Keep patient counts, optimization budgets,
and evaluation distributions fixed across the sweep. Preserve negative results.

## Next steps

1. Commit baseline implementation and durable compact result summaries.
2. Implement a resumable sweep, paired patient sampling, fixed evaluation cohorts,
   explicit interpolation and figures with per-local and paired uncertainty.
3. Run a small development check, then a locked multi-seed validation sweep.
4. Commit scripts/configs/results/figures and update this log at each milestone.

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
/home/mfjoel01/.conda/envs/appfl_env/bin/python -m appfl_bio_suite.experiments.methylation.alma3.sweep \
 --config src/appfl_bio_suite/experiments/methylation/configs/alma3-skew-sweep.json \
 --out local/output/alma3-skew-sweep \
 --release /lus/grand/projects/GeomicVar/mjoel/BIOSIM_SUITE/scratch/alma3-release-v3.1.2 \
 --worker-command /lus/grand/projects/GeomicVar/mjoel/BIOSIM_SUITE/scratch/alma3-runtime/bin/python \
 --publish docs/experiments/methylation/results/alma3-skew-sweep --commit-progress
```

Sweep checkpoint: seed 41, convergence 0 complete. Scores in `docs/experiments/methylation/results/alma3-skew-sweep-development`; resume command below remains valid.
