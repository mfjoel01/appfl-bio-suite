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
