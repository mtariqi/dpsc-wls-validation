# WLS (Gpr177/Wntless) in dental pulp stem cells: cross-dataset validation

**Question (Prof. Wei Hsu, ADA Forsyth Institute):** Is Gpr177/Wntless (*WLS* in human) expressed in, and enriched in, the dental pulp stem cell (DPSC) population?

This repository tests that single question in independent human scRNA-seq datasets. The discovery analysis (GSE164157, 5 donors) is in [mtariqi/dpsc-wls-human-scRNAseq](https://github.com/mtariqi/dpsc-wls-human-scRNAseq); this repository uses the same QC, marker panels and DPSC definitions so the results are directly comparable.

## Datasets

| Dataset | Tissue | Comparison | Donors | Status |
|---|---|---|---|---|
| [GSE164157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164157) | Fresh pulp, 5 donors | DPSC vs other pulp cells | 5 | ✅ discovery (other repo) |
| [GSE185222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185222) | Fresh pulp; sound tooth (DTP_01) used | DPSC vs other pulp cells | 1 | ⬜ |
| [GSE202476](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202476) | Fresh pulp, 13 y, immature third molar | DPSC vs other pulp cells | 1 | ⬜ |
| [GSE227731](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227731) ([PRJNA946721](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA946721)) | DPSCs and PDLSCs | DPSC vs PDLSC | 1 library each | ⬜ |

## Methods in brief

- **QC:** `nFeature_RNA > 200`, `nCount_RNA > 500`, `percent.mt < 20`; doublets removed with scDblFinder; Harmony by sample when there is more than one.
- **DPSCs:** four definitions, as in the discovery repo. A (primary): mesenchymal sub-clusters with mean stem score above background (as in 06_annotation.R). B: top 20% stem score. C: ≥3 stem markers. D: MCAM⁺.
- **Statistics:** ≥3 donors: edgeR quasi-likelihood pseudobulk. Fewer donors: % WLS⁺, pooled log2FC and cell-level Wilcoxon, reported as **descriptive**.

## Run

```bash
Rscript scripts/00_install_packages.R
Rscript scripts/01_download_data.R
Rscript scripts/02_preprocess.R GSE202476     # then confirm cluster labels, re-run
Rscript scripts/03_WLS_analysis.R GSE202476
# repeat 02-03 for GSE185222 and GSE227731
Rscript scripts/04_summary_table.R
```

Step-by-step details and checkpoints: [`docs/VALIDATION_PLAN.md`](docs/VALIDATION_PLAN.md).

## Repository structure

```text
dpsc-wls-validation/
├── config/validation_datasets.yaml   # datasets, QC thresholds, discovery reference values
├── R/wls_validation_functions.R      # shared functions (QC, annotation, DPSC definitions, tests)
├── data/
│   ├── metadata/                     # sample sheets (tracked)
│   ├── annotation/                   # confirmed cluster labels (tracked)
│   ├── raw/                          # GEO downloads (not tracked)
│   └── processed/                    # Seurat objects (not tracked)
├── docs/                             # plan and setup notes
├── results/<accession>/              # tables and figures per dataset
├── results/WLS_cross_dataset_summary.csv   # final answer, one row per dataset
└── scripts/00-04
```

## Scope

Only WLS in DPSCs. DEG, GSEA, ligand–receptor and pathway analyses are deliberately out of scope.

## License

Code: MIT.
