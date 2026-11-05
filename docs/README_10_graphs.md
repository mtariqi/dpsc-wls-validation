# Scripts for the 10 graphs: PTCH1, PTCH2 and SMO in Axin2⁺ DPSCs

**Author:** Md Tariqul Islam
**Affiliation:** Laboratory of Professor Wei Hsu, The ADA Forsyth Institute, Somerville, MA, USA
**Date:** 4 November 2026
**Repository:** https://github.com/mtariqi/dpsc-wls-validation

## The 10 graphs

| Gene | Graphs | Files in `results/figures_composite/individual/` |
|---|---|---|
| PTCH1 | 2 violins | `PTCH1_GSE164157_violin`, `PTCH1_GSE202476_violin` |
| PTCH2 | 2 UMAPs + 2 violins | `PTCH2_GSE164157_UMAP`, `PTCH2_GSE164157_violin`, `PTCH2_GSE202476_UMAP`, `PTCH2_GSE202476_violin` |
| SMO | 2 UMAPs + 2 violins | `SMO_GSE164157_UMAP`, `SMO_GSE164157_violin`, `SMO_GSE202476_UMAP`, `SMO_GSE202476_violin` |

Each graph is saved as PNG (300 dpi) and PDF.

## Scripts

| Script | Role |
|---|---|
| `scripts/07_composite_figures.R` | Makes the 10 graphs (this folder) |
| `R/wls_validation_functions.R` | Shared functions sourced by script 07 (in the repository) |
| `scripts/02_preprocess.R` | Earlier step that built the GSE202476 object (QC, clustering, labels) |
| `dpsc-wls-human-scRNAseq/scripts/06_annotation.R` | Earlier step that built the annotated GSE164157 object |

## Command used

Run from the root of the `dpsc-wls-validation` repository:

```bash
Rscript scripts/07_composite_figures.R \
  ~/dpsc-wls-human-scRNAseq/data/processed/GSE164157_annotated.rds \
  "PTCH1:violin,PTCH2:umap_violin,SMO:umap_violin"
```

## Method

- **UMAPs:** Seurat `FeaturePlot` of log-normalised expression on each dataset's UMAP (Harmony-integrated for GSE164157); grey = not detected; axes UMAP_1 / UMAP_2.
- **Violins:** Seurat `VlnPlot` by group. **AXIN2⁺ pulp cells** (Axin2⁺ DPSCs, following Yoshida et al. 2025) = cells with at least one AXIN2 count among non-perivascular pulp cells (pulp fibroblast, MSC progenitor, odontoblast labels); boxed in blue. Groups with fewer than 50 cells are not shown (AXIN2⁺ pulp cells are kept if at least 10).
- **Datasets:** GSE164157 (fresh pulp, 5 donors, 23,174 cells; 56 AXIN2⁺ pulp cells) and GSE202476 (young tooth, 11,298 cells; 301 AXIN2⁺ pulp cells).
- **Table:** `results/figures_composite/composite_pct_expressing_by_Yoshida_group_PTCH1_PTCH2_SMO.csv` gives the number of cells, % expressing and mean expression per gene, dataset and group.
- **Environment:** R 4.6, Seurat 5.6.0; exact package versions in `renv.lock`. Random seed 2026.
