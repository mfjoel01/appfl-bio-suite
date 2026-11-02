# Completed MultiSuSiE artifacts

The pipeline completed on 29 October 2026 and the matched comparison completed on
30 October. On 2 November, all worktree-local artifacts were moved into the primary
`appfl-bio-suite` checkout using same-filesystem renames, with inode identity checked.
The source worktree is no longer required for artifact storage. It was left in place
for deletion; its local artifact entries are convenience symlinks into the archive.

## Locations in the primary checkout

These are local coordinator artifacts, excluded from Git. They are not downloads
provided by a fresh clone. Paths below are relative to the primary checkout.

| Artifact | Path |
| --- | --- |
| Matched comparison report and numerical outputs | `local/output/matched-multisusie-susiex-20261029/` |
| Full MultiSuSiE grid, selection and numerical parity evidence | `local/output/multisusie-full-20261023/` |
| Three-method pilot | `local/output/multisusie-comparison-20261023/` |
| Editable poster | `local/scratch/poster-multisusie/native/multisusie-fine-mapping-poster.pptx` |
| 72×36-inch print PDF | `local/scratch/poster-multisusie/native/multisusie-fine-mapping-poster-print-72x36.pdf` |
| Poster source, standard PDF, preview and audits | `local/scratch/poster-multisusie/` |
| Complete archive and exported source checkout | `local/archive/fine-mapping-multisusie-20261102/` |

The output and poster paths above link into the archive. The existing primary-checkout
`local/scratch/poster` was not overwritten. The MultiSuSiE poster reports the original
pipeline; the subsequent matched-estimator findings are in
[MATCHED_RESULTS.md](MATCHED_RESULTS.md), not inserted into that poster.

The archive retains every output/checkpoint, unsuccessful attempt, scheduler log,
frozen source/configuration snapshot, poster source/export, backend vendor checkout,
Python environment and poster dependency directory. Its `README-HANDOFF.md` and
`RELOCATION.json` record the local path map and source commit.

## Completion and verification

The full pipeline has 240,900 accepted main-grid fits plus 11,850 pooled comparisons.
All comparisons passed numerical PIP tolerances; one exact CS membership mismatch
remains explicitly reported. Its `completion.json` records completion and poster hashes.
The matched comparison has 139,650 accepted paired instances / 558,600 method-policy
rows; `accepted.json` hashes its final numerical outputs.

Relocation verification checked 2,457 input/result hashes, all 931 matched task
acceptance records, and BED size/mtime identities using the new paths. Seven poster
export/source hashes were checked separately. Records in the archive are
`relocation-verification.json` and `poster-relocation-verification.json`.

To repeat the main relocation check from the primary checkout:

```bash
python local/archive/fine-mapping-multisusie-20261102/verify_relocation.py
```

Historical manifests and scripts retain their original absolute paths. Rewriting
those files would invalidate frozen provenance. The verifier maps `old_checkout`
to `archive_checkout` using `RELOCATION.json`; use the same mapping when inspecting
historical paths. Shared input data already stored in the primary checkout remain
external dependencies and were not duplicated.

This is a completed-run archive. Prepare a new run with current paths to rerun the
experiment; do not restart the historical scheduler or monitor scripts. Archived
virtual-environment activation scripts and executable shebangs retain the old prefix.
Use `local/venv-multisusie/bin/python` within the archive directly with `PYTHONPATH`
set to the desired archived `src` or frozen `code/src`, or recreate the environment
with the optional `multisusie` extra. Raw results and reports need no backend install.
