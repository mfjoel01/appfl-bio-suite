# `local/` — the coordinator's working directory

Everything in this directory except this file is gitignored, permanently, from the repo's
first commit. That is deliberate and it is not a convenience: sensitive material in git
history is not really removable, so the rule has to predate any content.

## Why it exists

The repo is public. Running a federation is not. A coordinator accumulates material that is
genuinely useful to them and genuinely must not be published:

- live endpoint UUIDs, which are access-relevant even though they are not secrets
- partner configs mid-negotiation, before a site has agreed to be named
- correspondence with partner administrators
- run outputs, scratch scripts, and half-finished analyses

None of that belongs in a public tree, and judging it case by case at commit time is how
leaks happen. Making it a directory boundary means the judgment is made once.

## What goes here — the actual layout (post-migration, 2026-08-21)

```
local/
├── federation.yaml            # YOUR real federation: identity, sites, endpoint UUIDs (not yet authored)
├── configs/
│   ├── flamby-heart-disease/  # migrated driver configs; run from the .suite-adapted siblings, originals are history
│   ├── gwas/                  # migrated Globus Compute driver/endpoint configs + scripts; same .suite-adapted convention
│   ├── endpoint-updates/      # staged ~/.globus_compute config replacements (see its README before applying)
│   └── env-updates/           # staged .pth files repointing appfl_env's editable installs at local/vendor/
├── data/                      # staged experiment inputs: fine-mapping/, flamby-heart-disease/, gwas/ (provenance in data/README.md)
├── output/                    # run results, per experiment (incl. migrated gwas/ and flamby-heart-disease/ run history)
├── logs/fedfm/                # migrated fine-mapping PBS/driver logs
├── papers/                    # migrated write-ups: fine-mapping/, flamby-heart-disease/, gwas/ (verbatim bodies + migration headers)
├── runs/gwas-5site-fullscale/ # the migrated 5-site full-scale GWAS run tree (6 GB, results of record)
├── vendor/
│   ├── fine-mapping/          # vendored binaries (plink, plink2, SuSiEx) — repo-root `vendor` symlinks here
│   ├── APPFL/                 # APPFL clone @ 21bd9436 — the live appfl_env's editable-install target (import appfl resolves here once env-updates is applied)
│   └── FLamby/                # FLamby clone @ edacf54d — the live appfl_env's .pth install target (import flamby resolves here once env-updates is applied)
├── worker_modules/gwas/       # live endpoint PYTHONPATH target — what Globus Compute workers import; keep byte-identical to what partners validated
├── partner_bundles/gwas/      # site data bundles as sent to partners (+ checksums)
├── archive/                   # pristine pre-migration mirrors and git provenance (patches, bundles) — do not edit
└── MIGRATION/                 # migration manifests and audits — start here for "where did X go?"
```

Start by copying the committed example and filling it in:

```bash
cp federation.yaml.example local/federation.yaml
$EDITOR local/federation.yaml
```

Every generator, config, doc template, and preflight check in the suite reads from that one
file. See [docs/coordinator/new-federation.md](../docs/coordinator/new-federation.md).

## What does NOT go here

**Partner-facing documentation.** This is the one rule people get wrong, and it has already
cost a partner real time.

In an earlier version of this project the partner setup guide lived in a git-ignored
directory and travelled only by email. A partner ended up working from a stale copy whose
step numbers were off by one and asked questions about a step that no longer existed.

Every partner-facing document is therefore a tracked file under `docs/partner/`. A partner
runs `git pull` and is current, and you can say "you are on commit X" and have that mean
something. If you find yourself writing a partner-facing doc in `local/`, it belongs in
`docs/partner/` instead — and if it cannot go there because it names another partner or
carries internal detail, then that content should not be in a partner's hands at all.
