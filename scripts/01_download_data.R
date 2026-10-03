# =============================================================================
# 01_download_data.R
# Downloads processed count files for the GEO validation datasets.
# GSE227731 is the processed scRNA-seq record of PRJNA946721 (no alignment needed).
# Run from the repository root:  Rscript scripts/01_download_data.R
# =============================================================================
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
