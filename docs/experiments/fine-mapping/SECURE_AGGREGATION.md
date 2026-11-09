# Secure aggregation for fine-mapping

The APPFL fine-mapping trainer and aggregator support opt-in pairwise masked sums.
The default remains plaintext summary-statistic transport. The standalone `fedfm` driver
and centralized baseline do not use this protocol.

## APPFL native compatibility

This uses APPFL's trainer/aggregator lifecycle, but not its native `SecureAggregator`.
The pinned APPFL 1.10.0 implementation's `mask_update` unconditionally converts inputs
to float32 and adds unit-scale Gaussian masks. For example, `2**40 + 1` becomes
`2**40`, and at that magnitude the mask itself is rounded away. It therefore cannot
preserve the moment protocol's precision or hide large moments. A regression probe in
`test_fine_mapping_secure_aggregation.py` checks the actual installed implementation.
An `appfl`/`appfl-native` protocol request fails explicitly rather than silently
falling back or claiming native protection. The integer pairwise protocol remains in use.

Both [SuSiEx and MultiSuSiE](MULTISUSIE.md) run downstream of the same secure pooling
step; the key setup and transport requirements are identical for both backends.

## Scope and threat model

Each site masks `G = X'X`, `u = 1'X`, `c = X'y`, and the `(q, w, n)` phenotype tensor
before returning its payload. The coordinator sums the masked tensors before decoding
any moments. Opposite pairwise masks cancel in integer modular arithmetic; the existing
pooling and standardization code then operates on one pooled block per ancestry.

This protects individual sites' moments from an honest-but-curious coordinator that
does not possess the pairwise keys. At least two noncolluding sites must contribute to
every block. With two sites, either site can infer the other's moments if it also
obtains the aggregate. This is not differential privacy, malicious-party-secure MPC,
or a defense against compromised workers or a coordinator shipping malicious trainer
code. Participants must trust/review the executed trainer and authorize the run locally.
Transport authentication and endpoint access controls are still required. The protocol
has not undergone an independent cryptographic audit.

**Public metadata remains public:** site identities, sample sizes, ancestry composition,
variant identities/order/alleles, loci, phenotype instance names, missing-variant counts,
and GA4GH provenance. Aggregate statistics and downstream results are also disclosed.
The implementation does not hide cohort membership metadata or promise that aggregate
statistics cannot reveal information about individuals.

## Requirements

- At least two sites, with a fixed roster for the whole exchange.
- Identical ordered genotype and phenotype blocks, including the same ancestry groups,
  locus windows, and ordered, allele-harmonized variant lists at every site. Select common
  ancestries with `pops`, and prepare a common complete-variant list before launch.
  Heterogeneous ancestry layouts and site-specific missing-variant filtering are rejected;
  there is no zero-filling or implicit intersection after masking.
- All sites must return. Dropout, mixed plaintext/masked payloads, a changed roster/run,
  layout disagreement, and inconsistent pairwise keys abort aggregation.
- Disable lossy compression. Masked tensors use all 64 bits of `torch.int64` and must
  arrive unchanged. The default configuration disables compression.

## Configure the run

In `federation.yaml`, under `experiments.fine-mapping`, set:

```yaml
pops: [EUR]  # an example; every selected site must hold the same requested blocks
secure_aggregation:
  run_id: "study-unique-random-id-shard-0-attempt-1"
  participants: [site_a, site_b, site_c]  # exact client_id values in this experiment
```

Use a fresh, unpredictable run ID for **every launch, retry, and shard** (for example,
`python -c 'import secrets; print(secrets.token_hex(32))'`). The launcher propagates the
public configuration to both roles and checks the roster against the configured sites.
Omit `secure_aggregation` to use the original plaintext path.

For each unordered site pair, provision an independent random 32-byte secret through a
trusted site-to-site channel. Generate it with
`python -c 'import secrets; print(secrets.token_hex(32))'` and give the same value to the
two peers only. Never put secrets in federation YAML, coordinator storage, shared output,
logs, or submitted configuration. Fresh pairwise secrets for every attempt are recommended.

Each worker has a private JSON authorization file. For site_a:

```json
{
  "client_id": "site_a",
  "run_id": "study-unique-random-id-shard-0-attempt-1",
  "participants": ["site_a", "site_b", "site_c"],
  "keys": {
    "site_b": "<64 hex characters shared only with site_b>",
    "site_c": "<64 hex characters shared only with site_c>"
  }
}
```

Restrict the file and its directory to the worker account (`chmod 600` for the file), and
set `APPFL_FINEMAP_SECAGG_KEYS` to its absolute path in the **worker process environment**.
The worker must be able to create a sibling `.used` file. This exclusive marker prevents
a second emission with the same authorization file, including concurrent retries. A failed
attempt may consume it: provision a fresh file, run ID and keys rather than deleting the
marker. Do not copy an authorization file to another path to bypass this guard.

Launch with the normal `appfl-bio-suite run fine-mapping` workflow. Serial/loopback runs
can test the math but do not provide a privacy boundary: one process controls all sites.
The protocol intentionally provides no automatic key distribution or dropout recovery.

## Numerical behavior and outputs

Moments are rounded to fixed point with scale `2**20`, masked and summed modulo `2**64`,
then decoded as signed integers and divided by the scale. A site rejects nonfinite values
or any scaled magnitude at least `2**62 / number_of_sites`, preventing overflow of the
honest aggregate. Each pooled element's quantization error is at most
`number_of_sites / (2 * 2**20)`, apart from final float64 rounding. Integer genotype
moments and counts remain exact within the supported range. Phenotype moments can change
slightly; bit-identical SuSiEx or MultiSuSiE results or identical threshold decisions are not guaranteed.
`uplink_gram_dtype` does not reduce the masked wire size: every masked element is 8 bytes.

HMAC-SHA256 derives domain-separated pair keys from the run, roster and layout; SHAKE256
expands per-tensor seeds into masks. Pairwise confirmations catch accidental mismatched
keys/contexts before decoding. These checks do not authenticate malicious contributions.

Successful secure runs write `data/secure_aggregation.json` containing protocol version,
run ID, roster and scale. Existing site summaries and GA4GH provenance retain the original
per-site metadata. No pairwise keys or decoded per-site moments are written by this path.

Tests cover modular recovery, negative phenotype moments, exact genotype sums, the actual
shipped trainer and decoder, launch wiring, retry prevention, and failure cases. Run:

```bash
PYTHONPATH=src python -m pytest tests/test_fine_mapping_secure_aggregation.py \
  tests/test_fine_mapping_shipped_parity.py tests/test_launch_and_loopback.py \
  tests/test_fine_mapping_secure_backends.py
```
