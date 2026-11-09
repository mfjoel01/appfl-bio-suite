# Federated methylation — the data

Two sources feed the same input contract: a seeded synthetic generator, and the open
TCGA-LAML methylation cohort from the NCI Genomic Data Commons. Holdouts, the CpG panel
and site assignment are made at run time from the unsplit cohort, so every arm of a run
shares them.

## Input contract

`MethylationData` in `dataset.py` validates these arrays, stored as an NPZ file:

| Field | Shape | Meaning |
| --- | --- | --- |
| `X` | samples × CpGs | Floating beta values in [0, 1]; NaN denotes missing data |
| `y` | samples | Nonempty Unicode subtype labels |
| `cpg_ids` | CpGs | Unique ordered Unicode probe identifiers |
| `patient_ids` | samples | Unicode patient identifiers; required for training |
| `sample_ids` | samples | Unique Unicode measurement identifiers; required for training |
| `site_id` | samples, optional | Reserved; the simulation assigns sites and rejects it |

Loading disables pickle. Invalid dimensions, infinities, out-of-range betas, duplicate
probe or sample IDs and empty identifiers fail validation. Repeated measurements must
share a patient ID — never fabricate distinct IDs for them. A patient whose samples
carry conflicting labels is rejected. Each subtype needs at least three patients so it
can appear in the training, validation and test splits.

## Synthetic cohorts

`synthetic_data()` draws background betas from Beta(2, 5). A fraction
`signal_fraction` of probes is informative. Those probes are divided into disjoint
per-class panels, and a patient's own panel is drawn from Beta(5, 2) instead. CpGs and
classes are named `synthetic_*`. No score on this data says anything about leukemia.

The generator's settings are `RunConfig` fields, so a scenario YAML controls them:

| Field | Effect | Default |
| --- | --- | --- |
| `samples`, `cpgs`, `classes`, `seed` | Cohort size and reproducibility | 240, 512, 3, 42 |
| `signal_fraction` | Share of probes carrying class signal | 1.0 |
| `class_weights` | Relative class sizes, one weight per class | Balanced |
| `profile_noise` | Share of patients given another class's profile, label unchanged | 0 |
| `missing_fraction` | Share of beta values set to NaN | 0 |

The defaults reproduce the original balanced generator exactly. `profile_noise` is
declared label/profile discordance. It sets an irreducible error rate, applied
identically for every arm, and does not model a biological mechanism.

```bash
appfl-bio-suite simulate methylation --scenario ci-tiny --out local/methylation-data
appfl-bio-suite simulate methylation --verify local/methylation-data/manifest.json
```

`simulate` writes an unsplit `cohort.npz` and a manifest recording every setting, the
generator's checksum, the suite commit and the cohort's SHA-256. `--verify` re-checksums
the cohort against it. A run given that cohort through `--data-root` recognizes it as
synthetic only when the manifest's checksum matches the file.

## TCGA-LAML from the GDC

The real-data cohort is **adult acute myeloid leukemia**, labeled by **FAB morphology**:
M1, M2, M4 and M5. It is not pediatric, and it is not molecular subtype.

```bash
python -m appfl_bio_suite.experiments.methylation.gdc --out local/tcga-laml-fab
```

The loader queries the GDC API for open Illumina 450K *Methylation Beta Value* files in
project TCGA-LAML. It joins each file to its case by exact case UUID and cross-checks the
patient submitter ID. On 2026-10-01 the query returned 194 files, each the single primary
peripheral-blood sample of a distinct patient:

| FAB label | Patients | Policy |
| --- | ---: | --- |
| M1, M2, M4 | 42 each | Include |
| M5 | 22 | Include |
| M0, M3 | 19 each | Exclude: below 20 patients |
| M6, M7 | 3 each | Exclude: below 20 patients |
| Not classified | 2 | Exclude: no FAB label |

That selects **148 patients**. Class counts were fixed from the cohort before any model
was fitted. The loader refuses to build if a future snapshot drops a selected class below
20 patients. It also refuses ambiguous joins, conflicting labels and multiple beta files
per patient.

Downloads use four workers. Each file is checked against its GDC size and MD5, retried on
failure, and cached under `raw/` by its GDC file UUID. An interrupted download restarts
with the same command and reuses only verified files. A completed cohort is never
overwritten. Run one loader per output directory.

| Output | Contents |
| --- | --- |
| `metadata.json` | Cached API responses, with retrieval time |
| `raw/` | Verified source files |
| `cohort.npz` | 148 × 482,421 `cg` probes in source order, float32, NaNs preserved |
| `manifest.json` | Class counts, selected rows, and hashes of the cohort, metadata and loader |

The cohort is about 13.8% missing, and its 225 MiB file is not normalized, imputed or
variance-filtered. The 50,000-probe panel is chosen at run time from training patients
only. The `tcga-laml` config applies the class threshold (`min_class: 10`) to training
patients, because the 20-patient gate applies to the whole cohort. It reserves one
patient of each class per site and disables batch shifts. The four sites are simulated:
TCGA is a single collection, not a federation.

**Terms.** The [GDC data analysis policy](https://gdc.cancer.gov/analyze-data/data-analysis-policies)
permits analysis and publication of open data and prohibits attempts to reidentify
participants. It is a data-use policy, not a Creative Commons licence. Acknowledge TCGA
and the GDC in any publication. Downloaded data stays in `local/` and out of version
control.

## Leakage controls

- Patients, not samples, are split: repeated measurements never straddle a holdout.
- Class counts come from training patients. A category below `min_class` maps to
  `other`, and a run left with fewer than two classes is rejected.
- Variance ranks CpGs on training rows only. All-missing, singly observed and constant
  training features are excluded. The panel is identical whatever the held-out values
  are, and a test checks this.
- `preprocessing.json` records the probe order, original column indices, label mapping,
  patient and sample IDs, disjoint split indices and site membership of every run.
- Sparse-call simulation never observes a missing beta. Requested coverage is capped by
  each sample's observed probes, and the actual mean coverage is reported.

## Other public cohorts considered

- **GSE49031** (NOPHO pediatric ALL, 450K). This was the original plan. The GEO record
  reports that raw and processed data were removed on 2024-04-26 over patient privacy
  concerns. No patient-level data was downloaded, and none should be sought elsewhere.
- **TARGET-AML** (pediatric AML). 113 primary patients have open betas, but no labeling
  gives four classes of at least 20 patients. By primary cytogenetic code only "Other"
  (41) and "Normal" (31) qualify. By FAB, only M2 (23) and M4 (21) do.
- **ALMA** (Marchi et al., 2025) was not used.
