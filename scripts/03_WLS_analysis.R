# =============================================================================
# 03_WLS_analysis.R
# The single question: is WLS (Gpr177/Wntless) expressed in, and enriched in, DPSCs?
# Usage:  Rscript scripts/03_WLS_analysis.R GSE185222
# =============================================================================
source("R/wls_validation_functions.R")
set.seed(2026)

acc <- commandArgs(trailingOnly = TRUE)[1]
cfg <- load_config(); ds <- cfg$datasets[[acc]]
out <- file.path("results", acc)
obj <- readRDS(file.path("data/processed", paste0(acc, "_preprocessed.rds")))
gene <- target_symbol(obj, cfg)

if (ds$design == "tissue") {
  # Secondary conditions (e.g. caries): one-line check only, no disease analysis
  if (ds$primary_condition %in% obj$condition && length(unique(obj$condition)) > 1) {
    other <- subset(obj, condition != ds$primary_condition)
    chk <- wls_summary(other, gene, "A_stem_subclusters")
    write.csv(chk, file.path(out, "secondary_conditions_check.csv"), row.names = FALSE)
    obj <- subset(obj, condition == ds$primary_condition)
  } else if (!ds$primary_condition %in% obj$condition) {
    warning("primary_condition '", ds$primary_condition, "' not in sample sheet; using all samples")
  }

  obj$group <- ifelse(obj$A_stem_subclusters, "Candidate_DPSC", obj$cell_type)
  write.csv(wls_summary(obj, gene, "group"), file.path(out, "WLS_by_cell_type.csv"), row.names = FALSE)

  is_mes <- obj$cell_type %in% MESENCHYME
  tests <- list()
  for (d in c("A_stem_subclusters", "B_top20pct_stem_score",
              "C_3plus_stem_markers", "D_MCAM_positive")) {
    dp <- obj[[d]][, 1]
    tests[[paste(d, "all")]] <- cbind(definition = d, wls_test(obj, gene, dp, !dp, "DPSC vs all other cells"))
    tests[[paste(d, "mes")]] <- cbind(definition = d, wls_test(obj, gene, dp, is_mes & !dp, "DPSC vs other mesenchyme"))
  }
  # Primary definition A vs each other cell type
  dpA <- obj$A_stem_subclusters
  for (ct in unique(obj$cell_type)) {
    ref <- obj$cell_type == ct & !dpA
    if (sum(ref) >= 20)
      tests[[ct]] <- cbind(definition = "A_stem_subclusters",
                           wls_test(obj, gene, dpA, ref, paste("DPSC vs", ct)))
  }
  plot_group <- "group"
} else {
  # Stem-cell design: DPSC vs PDLSC, plus WLS per DPSC sub-cluster
  is_dpsc <- obj$cell_type == "DPSC"
  tests <- list(cbind(definition = "sample label",
                      wls_test(obj, gene, is_dpsc, !is_dpsc, "DPSC vs PDLSC")))
  dp <- subset(obj, cell_type == "DPSC")
  dp$group <- paste0("DPSC_sub", dp$seurat_clusters)
  write.csv(wls_summary(dp, gene, "group"), file.path(out, "WLS_by_DPSC_subcluster.csv"), row.names = FALSE)
  obj$group <- obj$cell_type
  write.csv(wls_summary(obj, gene, "group"), file.path(out, "WLS_by_cell_type.csv"), row.names = FALSE)
  plot_group <- "group"
}

tests <- do.call(rbind, tests); rownames(tests) <- NULL
tests$dataset <- acc
write.csv(tests, file.path(out, "WLS_tests.csv"), row.names = FALSE)
print(tests)

# ---- Figures: one UMAP + one dot plot -------------------------------------------
p1 <- FeaturePlot(obj, gene, order = TRUE) + ggtitle(paste0(gene, " (Gpr177) - ", acc))
p2 <- DotPlot(obj, features = gene, group.by = plot_group) + coord_flip() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(acc)
ggsave(file.path(out, "WLS_umap.png"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(out, "WLS_dotplot.png"), p2, width = 7, height = 3.5, dpi = 300)
message("Results in ", out)
