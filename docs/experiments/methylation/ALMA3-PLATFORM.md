# Platform adaptation extension

Implemented, pending model-training evaluation. Configuration:
`src/appfl_bio_suite/experiments/methylation/configs/alma3-platform-sweep.json`.
Use the existing `alma3.sweep_parallel` command in ALMA3-HANDOFF.md with this
configuration and a NEW output directory `local/output/alma3-platform-sweep`.
Never resume the previous sweep with changed training code.

This separates platform observation from class composition. It retains five
sites, seven class-mix levels, three seeds and the same local/pooled/federated
comparators. Signature rotation, additive batch shift and subtype drift are
zero. The new optional measurement model applies only to shifted cohorts:

- Patient purity drawn uniformly from 0.55–0.95, mixed with a neutral beta=0.5
  background. This is a placeholder, not a measured healthy-cell methylome.
- Per-CpG Poisson read depths with site means 8, 6.5, 5, 3.5 and 2.
- False-positive methylation calls 2–8% and false negatives 8–14% across sites.
- Binomial methylated-read counts divided by depth; zero-depth CpGs are missing.
- Existing sparse CpG observation masks remain; all four classes remain at all
  sites with the original sweep's fixed patient counts and pooled class totals.

Parameters are illustrative stress conditions, not estimates of Oxford Nanopore
performance. This simulates methylation features after calling/alignment, not raw
nanopore signals or reads. Site IDs represent technical profiles, not geography,
ancestry, income or care quality. Unshifted evaluation removes these effects and
retains the original four-class retention comparator. RNG streams are disjoint
between patients/splits; effects never use patient labels.

Do not assume the resulting federation curve will match the earlier rotation
experiment. In particular, a stronger pretrained baseline or irrecoverable
measurement noise may reduce adaptation gains. Report those outcomes unchanged.
Class-balanced local baselines and held-out-site evaluation remain future
controls; they are not implemented by this extension.

Before making claims about a hospital, replace these assumptions with aggregate
summaries from its pipeline: observed CpG fraction and depth distribution,
caller/chemistry versions, specimen type and blast fraction, and diagnostic
class counts. A normal-cell reference would improve the purity model. Avoid
using geography as a surrogate for these measurements.

Published ALMA work already includes nanopore evaluation; do not claim ALMA has
never encountered this modality or was trained only in Europe without checking
the exact release's provenance. Also distinguish the publication's model from
the specific ALMA3 CV:e11 neural release used here.
Sources: https://www.nature.com/articles/s41467-025-62005-4 and
https://github.com/f-marchi/ALMA-classifier .
