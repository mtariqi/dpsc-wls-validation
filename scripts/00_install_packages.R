# ==============================================================================
# Title:        00_install_packages.R - Install R packages and record the environment
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      Installs the CRAN and Bioconductor packages used by the
                pipeline and records exact versions in renv.lock.
# Inputs:       None.
# Outputs:      renv.lock; project package library (renv/).
# Usage:        Rscript scripts/00_install_packages.R
# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,
#               patchwork, dplyr (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ------------------------------------------------------------------------------
# Details:
# Authotr: Md Tariqul Islam | Northeastern University | Prof. Hsu Lab | ADA Forsyth Inst.
# 00_install_packages.R  -  everything this repository needs (fresh R library)
# Run from the repository root:  Rscript scripts/00_install_packages.R
# ==============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c("Seurat", "harmony", "Matrix", "dplyr", "ggplot2", "yaml",
          "data.table", "hdf5r", "R.utils", "renv")
bioc <- c("GEOquery", "edgeR", "scDblFinder", "SingleCellExperiment")

install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)

# Record exact versions for reproducibility (creates renv.lock).
# renv::init() switches to an empty project library, so hydrate() then copies
# the packages installed above from the user library into it.
if (!file.exists("renv.lock")) renv::init(bare = TRUE, restart = FALSE)
renv::hydrate(packages = c(cran, bioc))
renv::snapshot(prompt = FALSE)

for (p in c(cran, bioc)) cat(sprintf("%-22s %s\n", p, as.character(packageVersion(p))))
