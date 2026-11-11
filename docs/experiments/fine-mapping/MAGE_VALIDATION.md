# Real-data validation: MAGE v1.0 + 1000 Genomes NYGC 30x

This document records the real-data federated fine-mapping benchmark validating FedFM
([ABOUT.md](ABOUT.md)) on **MAGE v1.0** (*Multi-ancestry Analysis of Gene Expression*;
Taylor et al., *Nature* 2024) paired with matched **1000 Genomes NYGC 30x GRCh38**
individual-level genotypes. Existing synthetic protocols ([DATA.md](DATA.md),
[SCIENTIFIC_RERUN.md](SCIENTIFIC_RERUN.md)) and shipped worker modules
(`dataset.py`, `trainer.py`) remain unchanged.

## Cohort, sample-swap remediation, and federation design

- **Expression & covariates:** MAGE v1.0 gene expression BED files and released covariate
  matrices (`sex`, `PC1`–`PC5`, `PEER1`–`PEER60`) staged from Zenodo record `10535719`
  (raw `N = 731` expression columns across 5 continental groups and 26 populations).
- **Documented sample-swap audit:** `audit_and_remediate_sample_swap()` programmatically
  flags `NA11919` (`SRR19762653`, batch 3, annotated `CEU`/`XY`), whose duplicate RNA-seq
  data actually belong to `HG00237` (retained separately as `SRR19762247`, batch 7, while the
  dropped batch-3 donor record is absent from eQTL inputs). Excluding `NA11919` yields the
  clean cohort (`N = 730` unique donors: `196 AFR`, `141 EUR`, `141 EAS`, `139 SAS`, `113 AMR`
  across 26 populations), while retaining a raw `N = 731` federation for paired sensitivity tests.
- **Genotypes:** Regional ±250 kb GRCh38 NYGC 30x VCF slices extracted directly via HTTP
  `Range` queries over `.vcf.gz.tbi` bin/linear indices (`vcf_slice.py`) and converted to
  PLINK1 `.bed/.bim/.fam` files.
- **Three-site geographic partition:** Clean donors are partitioned disjointly across
  `site_atlantic` (`N = 252`; 9 pops: `85 EUR, 84 AFR, 56 AMR, 27 SAS`), `site_africa_eurasia`
  (`N = 253`; 9 pops: `112 AFR, 56 EUR, 56 SAS, 29 AMR`), and `site_asiapacific`
  (`N = 225`; 8 pops: `141 EAS, 56 SAS, 28 AMR`), preserving multi-site representation for
  four ancestries while concentrating `EAS` at Asia-Pacific.

## Stratified 16-locus validation panel

`select_mage_panel()` chooses 16 non-overlapping autosomal MAGE eGenes (excluding extended
MHC, ±250 kb TSS cis-windows, totaling 46 published MAGE credible sets across 4 strata):

1. **PAGE GWAS colocalized (`page_coloc`):** `M02-ILRUN`, `M11-CETP`, `M13-NFATC1`, `M14-TPM4`
2. **High cross-ancestry $F_{\text{ST}}$ (`high_fst`):** `M00-NPNT` ($F_{\text{ST}}=0.78$),
   `M06-PCMTD1` ($0.72$), `M10-ARIH1` ($0.73$), `M12-ERN1` ($0.76$)
3. **High allelic heterogeneity (`multi_cs`, 8–9 CS):** `M01-THOC3`, `M08-DHRS4`,
   `M09-ACOT1`, `M15-RRP7A`
4. **Moderate single-CS resolution (`moderate_cs`):** `M03-ZNRF2P1`, `M04-GS1-124K5.4`,
   `M05-POT1`, `M07-NUDT13`

## Summary of validation outcomes (MAF ≥ 0.005, 16 loci)

| Comparison / Arm | Result |
| --- | --- |
| **SuSiEx: Federated vs Pooled** (`14/16` converged; same 2 loci `nonconverged`) | `worst_max_pip_abs_diff = 0.0`, `mean_pip_pearson_r = 1.0`, `frac_exact_cs_identical = 1.0`, `mean_cs_union_jaccard = 1.0` |
| **MultiSuSiE: Federated vs Pooled** (`16/16` converged) | `worst_max_pip_abs_diff = 4.64e-12` (median `1.46e-14`), `mean_pip_pearson_r = 1.0`, `frac_top_snp_identical = 1.0`, `frac_exact_cs_identical = 0.9375` (`mean_cs_union_jaccard = 0.999962`) |
| **Failure mode: site-local standardization** (`naive_site_standardized`) | Distorts within-ancestry LD (`max \|ΔLD\| = 0.597–0.806`) and Z-scores (`max \|ΔZ\| = 1.14–1.89`), driving `worst_max_pip_abs_diff = 0.9985` |
| **Covariate orthogonalization check** (`FWL (X_perp, Y_perp)` vs `Y_resid + raw X`) | Diverges up to `worst_max_pip_abs_diff = 0.9796` (`r = 0.9135`, `exact CS match = 0.357`) when global PCs/PEER factors correlate with ancestry structure |
| **Sample-swap sensitivity** (`N = 731` vs `730`) | `worst_max_pip_abs_diff = 0.01351`, `mean_pip_pearson_r = 0.9999`, top-SNP match `100%`, exact CS match `78.6%` (`11/14`) |
| **Federated MultiSuSiE resolution & MAGE concordance** | `16/16` converged, `74` CS (`19` singletons, `35` CS ≤ 5 variants, median min CS size `2.0`), **`38/46` (83%)** published CS captured (`100%` of published leads inside CS, mean best CS Jaccard `0.8190`) |
| **Federated SuSiEx resolution & single-cohort baselines** | Fed SuSiEx captures `15/46` (33%); solo sites capture `8–10/46` (17–22%); single-ancestry pooled arms capture `1–8/46` (2–17%) |

## Reproduction

```bash
# 1. Stage MAGE files, slice 1000G NYGC 30x regions, build federations, verify loopback, run benchmark & figures
PYTHONPATH=src python scripts/fine-mapping/run_mage_validation.py \
    --all --maf-analysis 0.005 --max-workers 16

# 2. Run offline MAGE real-data regression & shipped-pipeline unit tests (no network needed)
pytest -v tests/test_fine_mapping_mage.py
```

Publication figures (PDF + PNG) are written under `<out-root>/figures/` (`fig1_cohort_and_design` through `fig5_regional_locus_exemplars`), and machine-readable tables under `<out-root>/reports/`.
