# ==============================================================================
# Title:        add_headers_and_style.R - standard headers and code style
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Created:      2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      Adds a standard header (title, author, affiliation, purpose,
#               inputs, outputs, usage, dependencies) to every R script in a
#               repository and re-formats the code to the tidyverse style guide
#               with the styler package. Code logic is not changed.
# Usage:        Rscript tools/add_headers_and_style.R <repo_path> [--no-style]
#               e.g. Rscript tools/add_headers_and_style.R .
#                    Rscript tools/add_headers_and_style.R ~/dpsc-wls-human-scRNAseq
# Notes:        Existing top-of-file comments are kept under "Details".
#               Running the tool twice does not add a second header.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath(".")
do_style <- !("--no-style" %in% args)

AUTHOR <- "Md Tariqul Islam"
AFFILIATION <- c(
  "Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,",
  "Somerville, MA, USA"
)
UPDATED <- format(Sys.Date(), "%Y-%m-%d")
# One line on how the code was developed. Edit to follow your lab or programme policy.
DEVELOPMENT <- c(
  "Written with AI assistance (Claude, Anthropic); reviewed, run and",
  "validated by the author."
)
MARKER <- "# Author:       Md Tariqul Islam"

# ---- Script descriptions -------------------------------------------------------
# Scripts not listed here still receive a header, with a generic description.
META <- list(
  "00_install_packages.R" = list(
    title = "Install R packages and record the environment",
    purpose = "Installs the CRAN and Bioconductor packages used by the pipeline and records exact versions in renv.lock.",
    inputs = "None.", outputs = "renv.lock; project package library (renv/).",
    usage = "Rscript scripts/00_install_packages.R"
  ),
  "01_download_data.R" = list(
    title = "Download validation datasets from GEO",
    purpose = "Downloads the processed count matrices of GSE185222, GSE202476 and GSE227731 with GEOquery.",
    inputs = "None (GEO accession numbers are set in the script).",
    outputs = "data/raw/<GEO>/ with the supplementary files and FILE_LIST.txt.",
    usage = "Rscript scripts/01_download_data.R"
  ),
  "02_preprocess.R" = list(
    title = "Quality control, doublets, clustering, annotation and DPSC definitions",
    purpose = "Loads one dataset, filters cells (genes > 200, UMIs > 500, mitochondrial < 20%), removes doublets (scDblFinder), integrates samples (Harmony), clusters, labels clusters by marker-panel module scores and defines candidate DPSCs (definitions A-D).",
    inputs = "data/raw/<GEO>/; data/metadata/<GEO>_samples.tsv; config/validation_datasets.yaml; optional data/annotation/<GEO>_cluster_labels.csv.",
    outputs = "data/processed/<GEO>_preprocessed.rds; results/<GEO>/qc_summary.csv, cluster_labels_auto.csv, mesenchyme_subcluster_scores.csv, umap_cell_types.png.",
    usage = "Rscript scripts/02_preprocess.R <GEO accession>"
  ),
  "03_WLS_analysis.R" = list(
    title = "WLS expression and enrichment in candidate DPSCs",
    purpose = "Measures WLS (% expressing cells, mean expression) by cell type and tests DPSCs against other cells: edgeR pseudobulk with >= 3 donors, descriptive Wilcoxon otherwise.",
    inputs = "data/processed/<GEO>_preprocessed.rds; config/validation_datasets.yaml.",
    outputs = "results/<GEO>/WLS_tests.csv, WLS_by_cell_type.csv, WLS_umap.png, WLS_dotplot.png.",
    usage = "Rscript scripts/03_WLS_analysis.R <GEO accession>"
  ),
  "04_summary_table.R" = list(
    title = "Cross-dataset summary table",
    purpose = "Collects the main DPSC test from each dataset into one table, with GSE164157 as the reference row.",
    inputs = "results/<GEO>/WLS_tests.csv; config/validation_datasets.yaml.",
    outputs = "results/WLS_cross_dataset_summary.csv.",
    usage = "Rscript scripts/04_summary_table.R"
  ),
  "05_figures_for_Hsu.R" = list(
    title = "Cell-type and WLS figures; Yoshida et al. 2025 DPSC definition",
    purpose = "Makes labelled UMAPs, marker dot plots, WLS UMAPs and violins for all datasets; applies the Yoshida et al. (2025) marker panels and AXIN2+ pulp-cell DPSC definition; tests WLS in AXIN2+ vs AXIN2- pulp cells with depth-adjusted logistic regression.",
    inputs = "data/processed/<GEO>_preprocessed.rds; annotated GSE164157 object from the discovery repository.",
    outputs = "results/figures_for_Hsu/ (PNG, PDF, *_WLS_by_Yoshida_groups.csv, *_Yoshida_depth_check.csv).",
    usage = "Rscript scripts/05_figures_for_Hsu.R <path to GSE164157_annotated.rds>"
  ),
  "06_gene_panel_figures.R" = list(
    title = "UMAP and violin plots for any gene list",
    purpose = "Makes one UMAP and one violin plot per gene and dataset with the same UMAPs and cell-type groups as the WLS figures (used for GLI1, PTCH1, HHIP, AXIN2, PTCH2, SMO).",
    inputs = "data/processed/<GEO>_preprocessed.rds; annotated GSE164157 object.",
    outputs = "results/figures_gene_panel/<GENE>/ (PNG, PDF); gene_panel_pct_expressing_*.csv.",
    usage = "Rscript scripts/06_gene_panel_figures.R <GSE164157 rds> [GENES] [DATASETS]"
  ),
  "07_composite_figures.R" = list(
    title = "Composite figures in Prof. Hsu's layout",
    purpose = "For GSE164157 and GSE202476: labelled UMAP, gene UMAP and violin by Yoshida group with AXIN2+ pulp cells boxed; also saves every panel separately.",
    inputs = "data/processed/GSE202476_preprocessed.rds; annotated GSE164157 object.",
    outputs = "results/figures_composite/ (composites, individual/ panels, composite_pct_expressing_by_Yoshida_group_*.csv).",
    usage = "Rscript scripts/07_composite_figures.R <GSE164157 rds> \"GENE[:layout],...\""
  ),
  "wls_validation_functions.R" = list(
    title = "Shared functions for the validation pipeline",
    purpose = "Marker panels, data loading, QC, doublet removal, clustering, cluster annotation, DPSC definitions (A-D) and WLS test functions shared by scripts 02-05.",
    inputs = "Sourced by other scripts.", outputs = "None (function definitions).",
    usage = "source(\"R/wls_validation_functions.R\")"
  ),
  # discovery repository (dpsc-wls-human-scRNAseq)
  "02_download_GEO.R" = list(
    title = "Download GSE164157 from GEO",
    purpose = "Downloads the GSE164157 count matrices for five donors with GEOquery.",
    inputs = "None.", outputs = "data/raw/GSE164157/.", usage = "Rscript scripts/02_download_GEO.R"
  ),
  "04_QC.R" = list(
    title = "Quality control of GSE164157",
    purpose = "Filters cells with > 200 genes, > 500 UMIs and < 20% mitochondrial reads.",
    inputs = "data/processed/GSE164157_raw.rds.", outputs = "data/processed/GSE164157_qc.rds; QC tables and figures.",
    usage = "Rscript scripts/04_QC.R"
  ),
  "05_clustering.R" = list(
    title = "Normalisation, PCA and first clustering",
    purpose = "Log-normalises counts, selects variable genes, runs PCA and clusters cells before integration.",
    inputs = "data/processed/GSE164157_qc.rds.", outputs = "data/processed/GSE164157_clustered.rds.",
    usage = "Rscript scripts/05_clustering.R"
  ),
  "06a_batch_doublets.R" = list(
    title = "Doublet removal and Harmony integration",
    purpose = "Removes doublets per donor with scDblFinder and integrates donors with Harmony (20 PCs); clusters at resolution 0.5.",
    inputs = "data/processed/GSE164157_clustered.rds.", outputs = "data/processed/GSE164157_integrated.rds; results/06a_batch/.",
    usage = "Rscript scripts/06a_batch_doublets.R"
  ),
  "06_annotation.R" = list(
    title = "Cell-type annotation and candidate DPSCs",
    purpose = "Labels clusters with marker-panel module scores (confirmed manually), re-clusters the mesenchyme (Harmony dims 1-15, resolution 0.4) and calls candidate DPSCs (stem score > 0).",
    inputs = "data/processed/GSE164157_integrated.rds; data/annotation/cluster_labels.csv.",
    outputs = "data/processed/GSE164157_annotated.rds, GSE164157_mesenchyme.rds; results/06_annotation/.",
    usage = "Rscript scripts/06_annotation.R"
  ),
  "07_WLS_expression.R" = list(
    title = "WLS expression by donor and pseudobulk tests",
    purpose = "Measures % WLS-positive cells per donor and group; tests DPSCs against other groups with edgeR quasi-likelihood pseudobulk (~ donor + group) and an lme4 mixed model.",
    inputs = "data/processed/GSE164157_annotated.rds.", outputs = "results/07_WLS/.",
    usage = "Rscript scripts/07_WLS_expression.R"
  ),
  "07b_WLS_sensitivity.R" = list(
    title = "Sensitivity analyses for the WLS result",
    purpose = "Repeats the WLS comparison under four DPSC definitions, tests WLS-WNT co-expression with depth adjustment and summarises per-donor QC.",
    inputs = "data/processed/GSE164157_annotated.rds, GSE164157_mesenchyme.rds.", outputs = "results/07b_sensitivity/.",
    usage = "Rscript scripts/07b_WLS_sensitivity.R"
  )
)

wrap <- function(text, prefix, width = 80) {
  lines <- strwrap(text, width = width - nchar(prefix))
  paste0(c(prefix, rep(strrep(" ", nchar(prefix)), length(lines) - 1)), lines)
}

make_header <- function(file, details) {
  m <- META[[basename(file)]]
  if (is.null(m)) {
    m <- list(
      title = "Analysis script", purpose = "See Details below.", inputs = "See Details below.",
      outputs = "See Details below.", usage = paste("Rscript", file)
    )
  }
  bar <- paste0("# ", strrep("=", 78))
  sub <- paste0("# ", strrep("-", 78))
  h <- c(
    bar,
    paste0("# Title:        ", basename(file), " - ", m$title),
    "# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs",
    MARKER,
    paste0("# Affiliation:  ", AFFILIATION[1]),
    paste0("#               ", AFFILIATION[-1]),
    paste0("# Last updated: ", UPDATED),
    sub,
    wrap(m$purpose, "# Purpose:      "),
    wrap(m$inputs, "# Inputs:       "),
    wrap(m$outputs, "# Outputs:      "),
    wrap(m$usage, "# Usage:        "),
    "# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,",
    "#               patchwork, dplyr (exact versions in renv.lock)",
    paste0("# Development:  ", DEVELOPMENT[1]),
    paste0("#               ", DEVELOPMENT[-1])
  )
  if (length(details)) h <- c(h, sub, "# Details:", details)
  c(h, bar, "")
}

process_file <- function(file) {
  x <- readLines(file, warn = FALSE)
  if (any(x == MARKER)) {
    message("  already has a header: ", file)
    return(invisible(FALSE))
  }
  # existing top-of-file comment block: keep its text, drop banner and file-name lines
  n <- 0
  while (n < length(x) && grepl("^#", x[n + 1])) n <- n + 1
  top <- if (n > 0) x[seq_len(n)] else character(0)
  keep <- top[!grepl("^#\\s*[=-]{5,}\\s*$", top) &
    !grepl(paste0("^#\\s*(R/|scripts/)?", basename(file), "\\s*$"), top) &
    !grepl("^#\\s*$", top)]
  body <- if (n > 0) x[-seq_len(n)] else x
  while (length(body) && body[1] == "") body <- body[-1]
  writeLines(c(make_header(file, keep), body), file)
  message("  header added: ", file)
  invisible(TRUE)
}

files <- list.files(repo, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("/renv/|/tools/", files)]
message("Repository: ", repo, " (", length(files), " R scripts)")
invisible(lapply(files, process_file))

if (do_style) {
  if (!requireNamespace("styler", quietly = TRUE)) {
    stop("Install styler first: Rscript -e 'renv::install(\"styler\")'")
  }
  message("Applying tidyverse style (styler) ...")
  res <- styler::style_file(files, style = styler::tidyverse_style)
  message(sum(res$changed), " file(s) re-formatted")
}
message("Done. Check that every script still parses (see SCRIPTS.md), then commit.")
