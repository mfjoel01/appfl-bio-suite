# Fine-mapping figures

Twenty-eight figures in three groups, from one command. Everything here is
**coordinator-side**: none of it is shipped to a partner, which is why matplotlib stays
off every partner's dependency surface.

```bash
python -m appfl_bio_suite.experiments.fine_mapping.figures.render \
    --results <fm_results.tsv> \
    --data-root local/data/fine-mapping \
    --detail-dir local/output/fine-mapping/detail-run \
    --federated <fed_fm_results.tsv> --run-log <a federated run log> \
    --out-dir local/output/fine-mapping/figures
```

Only `--results` and `--out-dir` are required. Every other input unlocks specific
figures and is skipped with a note when absent, so a partial run still produces
everything it can.

## The groups

| Module | Prefix | What it is about |
| --- | --- | --- |
| `eda.py` | `eda*` | the simulated package, before any fine-mapping |
| `paper_plots.py` | `pap*` | the SuSiEx paper's own figure forms, over these results |
| `federation.py` | `fed*` | the cost, the correctness, and the *value* of federating |

### `pap*` — what maps from the paper, and what does not

Reference: Yuan et al., medRxiv 2023.01.07.23284293v4.

The paper varies the **discovery configuration** — which populations were combined, at
what sample size. This experiment holds that fixed (six ancestry columns over 150,000
people) and varies the **locus** instead. So its x-axis becomes ours:

| paper | here |
| --- | --- |
| which populations were combined | cross-site LD-divergence stratum |
| discovery sample size | per-locus h² |
| cross-population r_g | cross-ancestry r_g (unchanged) |
| method (SuSiEx / PAINTOR / …) | analysis path (federated / centralized) |

| Figure | Paper analogue |
| --- | --- |
| `pap1_recall_grid` | Supp. Figs 10–24 — power along the two design arms |
| `pap2_diversity_panels` | **Figure 2 a–d** — power / high-PIP / coverage / resolution |
| `pap3_convergence` | Figure 3a — did the fit succeed, and how often |
| `pap4_power_coverage` | Figure 3b — the yield/calibration trade-off |
| `pap5_resolution` | Figure 4a, 4b — PIP and credible-set distributions |
| `pap6_high_confidence` | Figure 4c, 4d — counts at PIP > 0.95 |
| `pap7_ld_divergence` | Supp. Fig 39 — but see the caveat below |
| `pap8_scalability` | Supp. Figs 8–9 — runtime, and what drives it |
| `pap9_coding_quality` | Ext. Data Fig 6 — quality issues against PIP |
| `pap10_population_probability` | **Figure 4 / Discussion** — per-population causal probability |
| `pap11_locuszoom_*` | **Figure 5**, Ext. Data Figs 4–5 — stacked association + PIP tracks |
| `pap12_effect_concordance` | Figure 4f, 4g — effect sizes across ancestries |

**Extended Data Figure 7 (VEP functional impact) has no analogue and is not attempted**:
the genotypes are synthetic, so no variant carries a functional annotation to bin by.
Figure 1 and Extended Data Figure 1 are method schematics rather than plots of data.

## Five-arm fed6 comparison

`fed6_what_federation_buys` shows exactly five arms: ANL alone, Covenant alone,
MBZUAI alone, a selected 50,000-person composition, and the full federation.
The selected arm is named `federation_smart_50k`; legacy proportional subsamples,
individual search candidates, repeat seeds and borrowed-LD controls stay in the
supporting tables. A missing selected arm is an error, never a relabeled old arm.

Use `scripts/fine-mapping/submit_composition.py` to prepare and submit the search
against a completed scientific run. It preserves the source run, partitions genomic
regions before selecting a cohort, tests the five AFR/EUR allocations with two site
recruitment policies and three sampling seeds, then freezes the development winner.
All five displayed arms are evaluated on the same held-out instances. Repeated
sampling seeds are combined into one bar, with locus-clustered uncertainty; copies
of the solo baselines are pairing controls, not additional independent observations.

The objective is unconditional causal-capture power. Coverage, false discovery rate
at PIP > 0.95, and high-confidence yield are reported alongside it. A winner is not
guaranteed to outperform a solo site. These comparisons vary ancestry composition,
site recruitment and noise as well as sample size; they do not isolate LD diversity.

## Argonne colors and Helvetica

All scientific plotting entry points use `core/plot_style.py`, including EDA,
paper panels, federation panels, QC report images, and GWAS Manhattan/QQ plots.
The packaged `core/plot_assets/argonne-palette.css` is the supplied palette;
light-mode categorical, sequential and diverging tokens are read directly from it.

Install licensed Helvetica fonts on the rendering machine, or set
`BIOSIM_FONT_PATH=/path/to/Helvetica-font-directory` (a single font file also works).
Use `BIOSIM_REQUIRE_HELVETICA=1` to reject rendering if Helvetica is unavailable.
Without the font, previews explicitly warn and use DejaVu Sans. The composition
run's `plots.json` records the actual font and the palette checksum. Font binaries
are not distributed with this repository.

```bash
/your/python scripts/fine-mapping/submit_composition.py \
  --source /path/to/completed-run --root /path/to/new-composition-run --submit
```

The dependency graph is development → selection → evaluation → finish. A parallel
plot job recolors the existing supported figures; finish adds the five-arm fed6.
`development_ranking.tsv`, `selection.json`, `evaluation_summary.tsv` and
`paired_vs_smart.tsv` document the choice and its held-out performance. Candidate
configs record exact ancestry and site counts. Covenant-first recruitment may use
fewer than three sites when their capacity suffices; the full federation always
uses all sites.

## The divergence stratum is confounded — read `pap7`'s docstring

The locus stratification is scored by the Frobenius norm of pairwise r² differences
between sites. That quantity is bounded by the magnitudes it differences — a window where
every site has r² ≈ 0 cannot score high — so it does not separate "the sites disagree
about LD" from "there is a lot of LD here". On the published package the two correlate at
**r = +0.945**, and a magnitude-free version of the same question shows no relationship
with resolution at all. `pap7` is drawn to make that visible rather than to hide it. Do
not read panels a–c as an LD-diversity result.

## The two figures that needed new plumbing

`fm_results.tsv` summarises each SuSiEx run to one row and deletes the working
directory. That is the right shape for an 11,850-instance sweep, but it discards three
things several of these figures need — all of which SuSiEx already computed:

* **per-credible-set records**, so coverage can be computed at all. The summary table
  reports how many causal variants were captured and how many sets were returned, but
  not *which* set held the causal variant, and coverage is a per-set quantity.
* **`POST-HOC_PROB_POP*`**, the population-specific causal probability — the quantity
  the paper thresholds at 0.8, and the only per-population statement the method makes.
* **per-variant PIP**, the PIP track of a LocusZoom-style panel.

The obvious fix — widening `parse_susiex` — is not available: `fedfm/` is vendored
byte-for-byte from upstream and `tests/test_fine_mapping_configs.py` asserts that with
`cmp`. So `detail.py` re-runs a stratified sample with `--keep-work` and reads the same
files from outside:

```bash
python -m appfl_bio_suite.experiments.fine_mapping.figures.detail \
    --config <pipeline_config.yaml> --out-dir local/output/fine-mapping/detail-run \
    --loci-per-stratum 3 --reps 3 --n-workers 16
```

A sample rather than the full sweep, on purpose: these are population quantities, a few
hundred instances estimate them far more finely than a figure can show, and the full
sweep is ~56 CPU-hours for records that would not move a line. **The headline metrics
stay the full 11,850 instances**; only the per-set panels use the sample, and each says
so in its own caption.

## The visual system

`figstyle.py` holds one palette, one set of rcParams, and the shared helpers. Every
palette in it was run through a colour validator (lightness band, chroma floor,
colour-vision-deficiency separation, normal-vision separation, surface contrast) rather
than chosen by eye. The recorded results, on the light surface these render on:

```
ANCESTRY (6)   adjacent-pair gates PASS  (worst CVD ΔE 9.1, normal-vision ΔE 19.6)
               ALL-PAIR gates FAIL       (green↔orange ΔE 3.2 protan)
SITE (3)       all-pair gates PASS       (worst CVD ΔE 15.3, normal-vision ΔE 20.8)
PATH (2)       all-pair gates PASS       (worst CVD ΔE 24.7, normal-vision ΔE 33.6)
ordinal ramps  monotone lightness, adjacent ΔL ≥ 0.06, light end clears 2:1
```

The ancestry all-pairs failure is load-bearing: six ancestries can be told apart in a
bar or a stack but **not in one scatter**, which is why `pap12` facets one ancestry pair
per panel exactly as the paper's Figure 4f/4g does. Do not add a seventh hue.

## Where the aggregator fits

`../plotting.py` is the aggregator's in-process entry point and still exports
`write_figures`, so a federated run produces its diagnostics without a second command.
It now delegates here. The old figure names map on as:

```
fig1_power_grid.png     -> pap1_recall_grid.png
fig2_cs_size.png        -> pap5_resolution.png     (panel b)
fig3_causal_pip.png     -> pap5_resolution.png     (panel a)
fig4_stratum_power.png  -> pap7_ld_divergence.png  (panel a)
```
