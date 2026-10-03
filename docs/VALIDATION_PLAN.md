# Cross-dataset validation of WLS in DPSCs

## Scope (Professor Hsu)

**One question:** Is Gpr177/Wntless (*WLS* in human) expressed in, and enriched in, the DPSC population?

Nothing else is in scope. DEG, GSEA, ligand–receptor analysis, WNT co-expression and disease comparisons are on hold unless he asks for them.

## Datasets

| Dataset | Data | Design | Comparison | Statistics |
|---|---|---|---|---|
| GSE164157 | 5 donors, fresh pulp | tissue | DPSC vs other pulp cells | done in [dpsc-wls-human-scRNAseq](https://github.com/mtariqi/dpsc-wls-human-scRNAseq) |
| GSE185222 | 4 teeth: DTP_01 sound, DTP_04 enamel caries, DTP_02 + DTP_05 deep caries | tissue | DPSC vs other pulp cells, **sound tooth (DTP_01) only** | descriptive (1 donor) |
| GSE202476 | 1 donor, 13 y, immature third molar | tissue | DPSC vs other pulp cells | descriptive (1 donor) |
| GSE227731 (= scRNA-seq of PRJNA946721) | 1 DPSC + 1 PDLSC library, processed count matrix | stem cells | DPSC vs PDLSC; WLS per DPSC sub-cluster | descriptive (1 library each) |

## Rules kept identical to GSE164157

- QC: `nFeature_RNA > 200`, `nCount_RNA > 500`, `percent.mt < 20`; scDblFinder; Harmony by sample when >1 sample.
- Same marker panels and the same four DPSC definitions (A–D). Definition A is the primary result.
- ≥3 donors: edgeR QL pseudobulk. Fewer donors: % WLS⁺, pooled log2FC and cell-level Wilcoxon, labelled **descriptive**.
- Caries samples (GSE185222): one-line check only (`secondary_conditions_check.csv`).

## Run order

```bash
Rscript scripts/00_install_packages.R
Rscript scripts/01_download_data.R            # GSE185222, GSE202476, GSE227731
Rscript scripts/02_preprocess.R GSE185222
# review results/GSE185222/cluster_labels_auto.csv
# save corrected labels as data/annotation/GSE185222_cluster_labels.csv, re-run 12
Rscript scripts/03_WLS_analysis.R GSE185222
# repeat 12-13 for GSE202476
Rscript scripts/02_preprocess.R GSE227731
Rscript scripts/03_WLS_analysis.R GSE227731
Rscript scripts/04_summary_table.R              # final table
```

## Deliverable

`results/WLS_cross_dataset_summary.csv`: one row per dataset, plus a WLS UMAP and dot plot per dataset. The script prints a one-line answer (e.g. "WLS is enriched in DPSCs in X of 4 datasets").

## Before running: things to confirm

- [ ] Definition A thresholds in `R/wls_validation_functions.R` (`stemcell > 0`, `perivasc > 0`) match `06_annotation.R`.
- [x] GSE185222 sample conditions: from the authors' code (github.com/vclabsysbio/scRNAseq_Dentalpulp, 02_Integration.R).
- [ ] GSE185222: confirm on GEO that GSM5608430 is titled DTP_05.
- [ ] GSE227731: check that `cell_regex` splits DPSC and PDLSC cells (script 12 prints example cell names and stops if not).
- [ ] GSE227731: confirm from the paper's Methods whether the cells were cultured, and how many donors were pooled.
