# =============================================================================
# 00_install_packages.R  -  everything this repository needs (fresh R library)
# Run from the repository root:  Rscript scripts/00_install_packages.R
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c("Seurat", "harmony", "Matrix", "dplyr", "ggplot2", "yaml",
          "data.table", "hdf5r", "R.utils", "renv")
bioc <- c("GEOquery", "edgeR", "scDblFinder", "SingleCellExperiment")

install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)

# Record exact versions for reproducibility (creates renv.lock)
if (!file.exists("renv.lock")) renv::init(bare = TRUE, restart = FALSE)
renv::snapshot(prompt = FALSE)

for (p in c(cran, bioc)) cat(sprintf("%-22s %s\n", p, as.character(packageVersion(p))))
