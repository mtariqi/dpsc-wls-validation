# ==============================================================================
# Title:        01_download_data.R - Download validation datasets from GEO
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      Downloads the processed count matrices of GSE185222, GSE202476
                and GSE227731 with GEOquery.
# Inputs:       None (GEO accession numbers are set in the script).
# Outputs:      data/raw/<GEO>/ with the supplementary files and FILE_LIST.txt.
# Usage:        Rscript scripts/01_download_data.R
# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,
#               patchwork, dplyr (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ------------------------------------------------------------------------------
# Details:
# Authotr: Md Tariqul Islam | Northeastern University | Prof. Hsu Lab | ADA Forsyth Inst.
# Downloads processed count files for the GEO validation datasets.
# GSE227731 is the processed scRNA-seq record of PRJNA946721 (no alignment needed).
# Run from the repository root:  Rscript scripts/01_download_data.R
# ==============================================================================

suppressPackageStartupMessages({ library(GEOquery) })
options(timeout = 3600)

for (acc in c("GSE185222", "GSE202476", "GSE227731")) {
  out <- file.path("data/raw", acc)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  message("== ", acc)

  # 1. Series-level supplementary files (usually a _RAW.tar with per-sample files)
  try(getGEOSuppFiles(acc, baseDir = "data/raw", makeDirectory = TRUE))
  for (t in list.files(out, "\\.tar$", full.names = TRUE)) untar(t, exdir = out)

  # 2. If the series had no matrices, fall back to sample-level files
  have_counts <- length(list.files(out, "matrix\\.mtx|\\.h5$|count|umi|\\.csv",
                                   recursive = TRUE, ignore.case = TRUE)) > 0
  if (!have_counts) {
    gse <- getGEO(acc, GSEMatrix = FALSE)
    for (gsm in names(GSMList(gse))) try(getGEOSuppFiles(gsm, baseDir = out))
  }

  files <- list.files(out, recursive = TRUE)
  writeLines(files, file.path(out, "FILE_LIST.txt"))
  message(length(files), " files; see ", file.path(out, "FILE_LIST.txt"))
}
message("\nNext: check the FILE_LIST.txt files, ",
        "then run scripts/02_preprocess.R <accession>")
