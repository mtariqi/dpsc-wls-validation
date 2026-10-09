# =============================================================================
# Authotr: Md Tariqul Islam | Northeastern University | Prof. Hsu Lab | ADA Forsyth Inst.
# 04_summary_table.R
# Final deliverable: one row per dataset answering the single question.
# Usage: Rscript scripts/04_summary_table.R
# =============================================================================
source("R/wls_validation_functions.R")
cfg <- load_config()
ref <- cfg$reference

rows <- list(data.frame(
  dataset = ref$dataset, design = ref$design, comparison = "DPSC vs all other cells",
  n_DPSC = 2847, pct_WLS_DPSC = ref$pct_WLS_DPSC, pct_WLS_ref = NA,
  log2FC = ref$log2FC_vs_all_other, p_value = ref$p_value, test = ref$test))

for (acc in names(cfg$datasets)) {
  f <- file.path("results", acc, "WLS_tests.csv")
  if (!file.exists(f)) { message("Not run yet: ", acc); next }
  t <- read.csv(f)
  main <- t[t$comparison %in% c("DPSC vs all other cells", "DPSC vs PDLSC") &
            t$definition %in% c("A_stem_subclusters", "sample label"), ]
  rows[[acc]] <- data.frame(dataset = acc, design = cfg$datasets[[acc]]$design,
    comparison = main$comparison, n_DPSC = main$n_DPSC,
    pct_WLS_DPSC = round(main$pct_WLS_DPSC, 1), pct_WLS_ref = round(main$pct_WLS_ref, 1),
    log2FC = round(main$log2FC, 2), p_value = signif(main$p_value, 2), test = main$test)
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
tab$WLS_expressed_in_DPSC <- ifelse(tab$pct_WLS_DPSC >= 5, "yes", "low")
tab$WLS_enriched_in_DPSC  <- ifelse(tab$log2FC > 0.5 & tab$p_value < 0.05, "yes",
                             ifelse(tab$log2FC < -0.5 & tab$p_value < 0.05, "depleted", "no"))
tab$WLS_enriched_in_DPSC[grepl("descriptive", tab$test) & tab$WLS_enriched_in_DPSC != "no"] <-
  paste(tab$WLS_enriched_in_DPSC[grepl("descriptive", tab$test) & tab$WLS_enriched_in_DPSC != "no"], "(1 donor)")

write.csv(tab, "results/WLS_cross_dataset_summary.csv", row.names = FALSE)
print(tab)
n_enr <- sum(grepl("^yes", tab$WLS_enriched_in_DPSC))
cat(sprintf("\nOne-line answer: WLS is enriched in DPSCs in %d of %d datasets.\n", n_enr, nrow(tab)))
