# Documentation

Organized by **audience first, experiment second**. Every document has exactly one
audience, and the audience is a directory boundary rather than a convention — which is
what makes "did I accidentally send a partner my internal notes" structurally impossible
rather than a judgment call.

---

## I am setting up a federation

Read in this order:

1. **[coordinator/new-federation.md](coordinator/new-federation.md)** — start here. Clone
   to running experiment, assuming no prior contact with anyone who built this.
2. [coordinator/architecture.md](coordinator/architecture.md) — how the pieces fit and why
   it is built this way.
3. [coordinator/reference-deployment.md](coordinator/reference-deployment.md) — a complete
   working deployment, as an example. Not a requirement.
4. [coordinator/troubleshooting.md](coordinator/troubleshooting.md) — symptom-indexed,
   with mechanisms. Skim it once before you need it.
5. [coordinator/releasing.md](coordinator/releasing.md) — upgrading pins without breaking
   the federation.

## I am a partner site

You should have received a **generated bundle** containing your two documents with every
value already filled in. Use those, not the templates here — the templates contain
unrendered placeholders.

If you did not receive a bundle, ask the coordinator to run
`appfl-bio-suite partner-bundle <experiment> --site <you>`.

The sources live in [partner/](partner/): Part 1 is endpoint setup, identical for every
experiment; Part 2 is what differs for yours.

## I want to understand or reproduce an experiment

Each experiment has the same four documents, always:

| | |
| --- | --- |
| `experiments/<name>/ABOUT.md` | What it is, how it is designed, and the record of runs |
| `experiments/<name>/DATA.md` | How the data comes into existence |
| `experiments/<name>/RUNBOOK.md` | How to run it, and what a healthy run looks like |
| `partner/experiments/<name>.md` | What a partner does (Part 2 of their setup) |

- [experiments/flamby-heart-disease/](experiments/flamby-heart-disease/)
- [experiments/gwas/](experiments/gwas/)
- [experiments/fine-mapping/](experiments/fine-mapping/) — federated cross-ancestry
  fine-mapping. Also carries [reference/](experiments/fine-mapping/reference/), the
  derivation and design notes it was ported with, and
  [results/](experiments/fine-mapping/results/), the completed runs.

Uniform on purpose: no experiment gets a special extra document, and none is omitted. A
missing one is a test failure.

---

## Why one experiment's documentation lives in two trees

`docs/partner/` is a **closed set**. The bundle generator copies from that tree and
nowhere else — enforced in code, with a test — so nothing outside it can reach a partner.

The cost is that an experiment's documentation is split between `docs/partner/experiments/`
and `docs/experiments/`. That is the right trade: the failure it prevents is a
coordinator's internal notes, another partner's name, or the fact that a dataset is
simulated reaching someone who should not receive it. The inconvenience it creates is
following a cross-link.

A concrete example of why the boundary matters: a GWAS partner's own setup guide never
mentions that their data was simulated. That belongs in `experiments/gwas/DATA.md`, which
is coordinator- and reviewer-facing and is never bundled.
