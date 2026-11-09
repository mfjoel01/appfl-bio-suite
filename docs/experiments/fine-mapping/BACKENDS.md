# Fine-mapping inference backends

The APPFL coordinator supports `susiex` (default) and `multisusie`. Both consume pooled
moments from the same shipped trainer, including optional [secure aggregation](SECURE_AGGREGATION.md).
Workers never install an inference backend. This selector applies to `appfl-bio-suite run
fine-mapping`; the standalone `fedfm` driver and centralized baseline remain SuSiEx.

## Install on the coordinator

For SuSiEx, use the existing `finemapping-sim` extra and
`scripts/fine-mapping/install_susiex.sh` / `install_plink.sh`.

For MultiSuSiE, from the repository root:

```bash
python -m pip install '.[finemapping-multisusie]' -c constraints.txt
```

The optional extra pins the upstream [MultiSuSiE package](https://github.com/jordanero/MultiSuSiE)
to commit `a89e091d45297d983118438f116c3eb5bc1219ed`. Its GPL-3.0 license applies to
that separately installed package. MultiSuSiE inference needs no SuSiEx or PLINK binary.
Simulation/data preparation can still require PLINK. It is not installed by the default
partner extra or the existing `all` extra.

## Configure and run

In `federation.yaml`, under `experiments.fine-mapping`:

```yaml
fine_mapping_backend: multisusie  # switch to susiex for the original engine
multisusie_rho: 0.75              # scalar cross-ancestry effect correlation: 0 <= rho < 1
multisusie_min_abs_corr: 0.5      # credible-set purity threshold
credible_set_level: 0.95
n_signals: 10
max_iter: 1000
tol: 1.0e-6
pval_thresh: 1.0e-5
maf: 0.005
```

Use the existing `appfl-bio-suite run fine-mapping` command with that federation config.
For direct APPFL server YAML, the selector is `server_configs.aggregator_kwargs.backend`.
The federation-aware preflight checks the selected backend before dispatching site work.
Missing MultiSuSiE raises an installation error, with no fallback to another engine.

To enable masked transport, additionally set `secure_aggregation.run_id` and
`secure_aggregation.participants`, and provision worker-local keys as documented in
[SECURE_AGGREGATION.md](SECURE_AGGREGATION.md). Every retry/shard needs fresh authorization.
Selecting an inference backend does not enable secure aggregation automatically.

## Scientific behavior and artifacts

The adapter calls upstream `MultiSuSiE.multisusie_rss` with pooled effect estimates,
standard errors, phenotype sample variances, sample sizes, and float64 arrays. It retains
the existing LD panel precision (float32 on disk). Summary statistics for this backend
are written with 17 significant digits. Upstream's remaining model defaults apply.

Within each ancestry the LD and summary-statistic lists must match. Across ancestries,
only variants passing the pooled MAF/complete-variant filters in **every** ancestry are
fitted. Both axes of each LD matrix are reordered to the common list, and position or
allele disagreement fails. Excluded IDs and fitted IDs are retained in `cs.execution.json`.
This explicit common-variant policy differs from SuSiEx's handling of ancestry-specific
variants; equivalent results between the two models are not expected.

MultiSuSiE's retained credible sets must also contain a marginal p-value below
`pval_thresh`, matching the purpose of that setting in SuSiEx. This filter affects sets,
not full-model PIPs. `keep_ambiguous: false` excludes A/T and C/G variants from the fit.

Result tables identify `fine_mapping_backend: multisusie`, the analysis protocol, and
`pip_definition: multisusie_full_model`. SuSiEx keeps its existing
`overall_retained_components` PIP definition. These quantities must not be compared as
though they were the same posterior. Convergence status and excluded causal variants
are explicit; nonconverged runs do not receive scored posterior probabilities.

Archived `cs.execution.json.gz`, `cs.posterior.npz.gz`, `cs.snp.gz`, `cs.cs.gz`, and
`cs.summary.json.gz` record model settings, full posterior arrays, PIPs, and retained sets,
even when scratch output is removed. PIP readers distinguish the backend definitions.

## Validate

With both backends installed and SuSiEx/PLINK on PATH:

```bash
PYTHONPATH=src python -m pytest tests/test_fine_mapping_multisusie.py \
  tests/test_fine_mapping_secure_aggregation.py tests/test_fine_mapping_shipped_parity.py \
  tests/test_launch_and_loopback.py tests/test_fine_mapping_appfl_parity.py
```

Tests cover reordered/subset ancestry panels, rejected allele mismatches, full-model PIP
semantics, a real MultiSuSiE fit against a direct upstream call, and secure/plaintext
trainer-to-inference comparisons for each backend. These are local checks, not a
cryptographic audit or a distributed-endpoint certification.
