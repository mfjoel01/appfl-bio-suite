# appfl-bio-suite

Cross-silo **federated learning experiments in computational biology**, built on
[APPFL](https://github.com/APPFL/APPFL) and Globus Compute.

Institutions train on their own data, on their own hardware, and return only results.
No data movement, no shared credential, no VPN, no inbound network access.

```
    YOU (coordinator)                          PARTNER SITES
    ┌──────────────┐
    │   driver     │ ──── task ────►    their cluster, their data
    │  aggregator  │ ◄─── result ───    (weights, or summary statistics)
    └──────────────┘
```

---

## The experiments

| | What it is | Data |
| --- | --- | --- |
| **[FLamby Fed-Heart-Disease](docs/experiments/flamby-heart-disease/ABOUT.md)** | Multi-round FedAvg on a public tabular benchmark with natural hospital splits | Partners download it; ~40 KB, CC BY 4.0 |
| **[Federated GWAS](docs/experiments/gwas/ABOUT.md)** | Single-round summary-statistic meta-analysis with polygenic-score evaluation | Simulated and distributed by the coordinator |
| **[Federated fine-mapping](docs/experiments/fine-mapping/ABOUT.md)** | Single-round cross-ancestry credible-set construction from per-site second moments — provably equal to a pooled fit, not an approximation of one | Simulated and distributed by the coordinator |

## Quickstart

```bash
git clone https://github.com/mfjoel01/appfl-bio-suite.git
cd appfl-bio-suite
python -m venv .venv && source .venv/bin/activate
pip install -e ".[all]" -c constraints.txt
```

Prove the install works — a complete federation on one machine, no partners, no network,
and nothing to fill in first:

```bash
appfl-bio-suite simulate gwas --scenario ci-tiny --out /tmp/gwas-ci
appfl-bio-suite run gwas --config loopback --driver serial --data-root /tmp/gwas-ci
```

Results land in `local/output/gwas-loopback/`. This is the same pair of commands CI runs.

Then stand up a real federation:

```bash
cp federation.yaml.example local/federation.yaml
$EDITOR local/federation.yaml
appfl-bio-suite preflight
```

**→ [docs/coordinator/new-federation.md](docs/coordinator/new-federation.md)** takes you
from here to a running experiment across institutions you recruit yourself.

## The CLI

```
appfl-bio-suite preflight [--experiment X]          environment, pins, configs, endpoints
appfl-bio-suite endpoint status [SITE]              is it online
appfl-bio-suite endpoint smoke [SITE]               round-trip: does it actually work
appfl-bio-suite identity validate <mapping.json>    catch mapping bugs offline, in a second
appfl-bio-suite simulate EXPERIMENT --scenario NAME  generate and split synthetic data
appfl-bio-suite run EXPERIMENT                      run it
appfl-bio-suite partner-bundle EXPERIMENT --site X  a partner's complete setup package
```

---

## How it works

**One config file.** `federation.yaml` declares who you are, who your partners are, and
what each is doing. Every generator, config, document template and check reads from it.
Nothing in `src/` contains a coordinator identity, an endpoint UUID, or a site name, and
CI verifies that.

**Generated partner bundles, not templates.** `partner-bundle` emits a directory with a
partner's configuration fully resolved — their site ID, their data assignment, *your*
identity, their scheduler's provider block — plus their two setup documents with no
placeholders left to interpret.

**A closed set for partner-facing content.** Bundles are built from `docs/partner/` and
nowhere else, enforced in code, so nothing outside that tree can reach a partner.

**Exact version pinning.** Every site runs byte-identical versions of the Globus Compute
stack. Skew does not fail at install time; it fails as a deserialization error several
rounds into a run on somebody else's cluster. `constraints.txt` pins it and
`preflight --check pins` verifies it.

**Failures that name the fix.** The
[troubleshooting registry](docs/coordinator/troubleshooting.md) is indexed by symptom and
gives the mechanism alongside the fix.

## Reproducibility

- **Every simulation run writes a manifest** recording the scenario, every seed, the suite
  commit, input and output checksums, and package versions.
  `simulate gwas --verify <manifest>` confirms a rerun mechanically.
- **The GWAS simulation is deterministic.** For a fixed seed, all 35 outputs are
  byte-identical, verified against recorded checksums.
- **Known limits are documented.**
  [docs/experiments/gwas/DATA.md](docs/experiments/gwas/DATA.md) states which input cohort
  is not regenerable, and what the shipped substitute does and does not reproduce.

## Documentation

Organized by audience.

| | |
| --- | --- |
| [docs/README.md](docs/README.md) | Index: which document, for whom, in what order |
| [docs/coordinator/](docs/coordinator/) | Running a federation |
| [docs/partner/](docs/partner/) | What a partner receives — a closed set |
| [docs/experiments/](docs/experiments/) | Per-experiment design, data, and runbooks |

## Requirements

Python 3.12. A Globus account. Partners willing to run a privileged process on their own
clusters.

Not required: a shared filesystem, a VPN, matching operating systems or schedulers, or
inbound network access.

## Licence

MIT. See [LICENSE](LICENSE) — and note the third-party terms recorded there for APPFL,
FLamby, and PGS Catalog data.
