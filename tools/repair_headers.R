# Repairs header continuation lines that are missing the leading "#".
args <- commandArgs(trailingOnly = TRUE)
for (repo in args) {
  files <- list.files(repo, "\\.R$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("/renv/", files)]
  for (f in files) {
    x <- readLines(f, warn = FALSE)
    bars <- which(grepl("^# =+$", x))
    if (length(bars) < 2 || !any(grepl("^# Author: ", x[bars[1]:bars[2]]))) next
    idx <- (bars[1] + 1):(bars[2] - 1)
    bad <- idx[!grepl("^#", x[idx]) & nzchar(trimws(x[idx]))]
    if (length(bad)) {
      x[bad] <- sub("^ ", "#", x[bad])
      writeLines(x, f)
      cat("repaired", length(bad), "line(s):", f, "\n")
    }
  }
}
