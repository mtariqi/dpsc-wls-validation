# ==============================================================================
# Title:        06_gene_panel_figures.R - UMAP and violin plots for any gene list
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      Makes one UMAP and one violin plot per gene and dataset with
                the same UMAPs and cell-type groups as the WLS figures (used
                for GLI1, PTCH1, HHIP, AXIN2, PTCH2, SMO).
# Inputs:       data/processed/<GEO>_preprocessed.rds; annotated GSE164157
                object.
# Outputs:      results/figures_gene_panel/<GENE>/ (PNG, PDF);
                gene_panel_pct_expressing_*.csv.
# Usage:        Rscript scripts/06_gene_panel_figures.R <GSE164157 rds> [GENES]
                [DATASETS]
# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,
#               patchwork, dplyr (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ------------------------------------------------------------------------------
# Details:
# Prof. Hsu, 2 Nov 2026: "the same analysis as Gpr177/Wls" for more genes:
# 4 UMAP plots and 4 violin plots per gene (one per dataset), nothing else.
# Genes: GLI1, PTCH1, HHIP (priority) and AXIN2 (validation).
# Groups and UMAPs are identical to 05_figures_for_Hsu.R, so the plots match the WLS figures.
# Usage (from the repository root):
#   Rscript scripts/06_gene_panel_figures.R ~/dpsc-wls-human-scRNAseq/data/processed/GSE164157_annotated.rds
# Optional: a comma-separated gene list as the 2nd argument, e.g. "GLI1,PTCH1,HHIP,AXIN2"
# ==============================================================================

source("R/wls_validation_functions.R")
suppressPackageStartupMessages(library(patchwork))
set.seed(2026)

args <- commandArgs(trailingOnly = TRUE)
discovery_rds <- args[1]
GENES <- if (length(args) >= 2) toupper(strsplit(args[2], ",")[[1]]) else c("GLI1", "PTCH1", "HHIP", "AXIN2")
# Optional 3rd argument: which datasets to plot, e.g. "GSE164157,GSE202476" (default: all four)
DATASETS <- if (length(args) >= 3) toupper(strsplit(args[3], ",")[[1]]) else
  c("GSE164157", "GSE202476", "GSE185222", "GSE227731")
message("Genes: ", paste(GENES, collapse = ", "), "   Datasets: ", paste(DATASETS, collapse = ", "))
out <- "results/figures_gene_panel"
ORDER <- c("Pulp_fibroblast", "Candidate_DPSC", "Perivascular_MCAM", "Odontoblast", "Endothelial",
           "Schwann_glia", "T_NK", "Myeloid", "B_plasma")
cfg <- load_config()

umap_name <- function(obj) {
  r <- grep("umap", Reductions(obj), value = TRUE, ignore.case = TRUE)
  r <- r[!grepl("mes", r)]
  pref <- r[grepl("harmony", r, ignore.case = TRUE)]
  if (length(pref)) pref[1] else r[1]
}
save_both <- function(p, file, w, h) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(file, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(file, ".pdf"), p, width = w, height = h, bg = "white")
}
missing_plot <- function(g, acc) {
  ggplot() + annotate("text", x = 0, y = 0, label = paste(g, "is not in the", acc, "count matrix"), size = 4) +
    theme_void() + ggtitle(paste0(g, " - ", acc))
}
prepare_tissue <- function(obj, dpsc_col) {
  grp <- as.character(obj$cell_type)
  grp[obj@meta.data[[dpsc_col]] %in% TRUE] <- "Candidate_DPSC"
  obj$plot_group <- factor(grp, levels = c(intersect(ORDER, unique(grp)), setdiff(unique(grp), ORDER)))
  obj
}

umaps <- setNames(lapply(GENES, function(x) list()), GENES)
vlns  <- setNames(lapply(GENES, function(x) list()), GENES)
tabs  <- list()

panel <- function(obj, acc, subtitle) {
  red <- umap_name(obj)
  min_cells <- if (ncol(obj) > 5000) 50 else 20
  big <- names(which(table(obj$plot_group) >= min_cells))
  ob <- subset(obj, plot_group %in% big)
  ob$plot_group <- droplevels(ob$plot_group)
  for (g in GENES) {
    if (!g %in% rownames(obj)) {
      umaps[[g]][[acc]] <<- missing_plot(g, acc); vlns[[g]][[acc]] <<- missing_plot(g, acc); next
    }
    x <- FetchData(obj, g)[, 1]
    if (max(x) == 0) {
      # Seurat colours an all-zero gene in one colour, which looks like expression: draw it grey instead
      p_u <- FeaturePlot(obj, g, reduction = red, raster = FALSE, cols = c("lightgrey", "lightgrey")) + NoLegend() +
        ggtitle(paste0(g, " - ", acc), subtitle = paste0(subtitle, " - not detected in any cell")) +
        labs(x = "UMAP_1", y = "UMAP_2")
    } else {
      p_u <- FeaturePlot(obj, g, reduction = red, order = TRUE, raster = FALSE) +
        ggtitle(paste0(g, " - ", acc), subtitle = subtitle) + labs(x = "UMAP_1", y = "UMAP_2")
    }
    p_v <- VlnPlot(ob, g, group.by = "plot_group", pt.size = 0.05) + NoLegend() +
      ggtitle(paste0(g, " - ", acc)) + labs(x = NULL)
    save_both(p_u, file.path(out, g, paste0(acc, "_UMAP")), 6, 5)
    save_both(p_v, file.path(out, g, paste0(acc, "_violin")), 7, 4.5)
    umaps[[g]][[acc]] <<- p_u; vlns[[g]][[acc]] <<- p_v
    tabs[[length(tabs) + 1]] <<- data.frame(gene = g, dataset = acc, group = obj$plot_group, x = x) %>%
      group_by(gene, dataset, group) %>%
      summarise(n_cells = n(), pct_expressing = round(100 * mean(x > 0), 1),
                mean_lognorm = round(mean(x), 3), .groups = "drop")
  }
}

# ---- GSE164157 (discovery) -------------------------------------------------------
if ("GSE164157" %in% DATASETS) {
  d <- readRDS(path.expand(discovery_rds))
  d <- subset(d, cell_type %in% c("Epithelial", "Erythrocyte", "Unassigned"), invert = TRUE)
  d <- prepare_tissue(d, "candidate_DPSC")
  panel(d, "GSE164157", paste0("Fresh pulp, 5 donors, ", ncol(d), " cells")); rm(d); gc()
}

# ---- GSE202476 and GSE185222 (sound tooth) --------------------------------------------
for (acc in intersect(c("GSE202476", "GSE185222"), DATASETS)) {
  o <- readRDS(file.path("data/processed", paste0(acc, "_preprocessed.rds")))
  pc <- cfg$datasets[[acc]]$primary_condition
  if (pc %in% o$condition) o <- subset(o, condition == pc)
  o <- prepare_tissue(o, "A_stem_subclusters")
  sub <- if (acc == "GSE202476") "Young tooth" else "Sound tooth (DTP_01)"
  panel(o, acc, paste0(sub, ", ", ncol(o), " cells")); rm(o); gc()
}

# ---- GSE227731 (DPSC sub-clusters and PDLSCs) -------------------------------------------
if ("GSE227731" %in% DATASETS) {
  o <- readRDS("data/processed/GSE227731_preprocessed.rds")
  o$plot_group <- ifelse(o$cell_type == "DPSC", paste0("DPSC_sub", o$seurat_clusters), "PDLSC")
  o$plot_group <- factor(o$plot_group, levels = c(sort(grep("^DPSC", unique(o$plot_group), value = TRUE)), "PDLSC"))
  panel(o, "GSE227731", paste0("DPSCs + PDLSCs, ", ncol(o), " cells")); rm(o); gc()
}

# ---- One page per gene: 4 UMAPs in a row, 4 violins in a 2 x 2 grid -----------------------
for (g in GENES) {
  save_both(wrap_plots(umaps[[g]], nrow = 1), file.path(out, g, paste0(g, "_all_datasets_UMAP_row")), 5.5 * length(umaps[[g]]), 5)
  save_both(wrap_plots(vlns[[g]], ncol = 2), file.path(out, g, paste0(g, "_all_datasets_violins")), 14, 9)
}
tab <- do.call(rbind, tabs)
# default run keeps the original file name; a custom gene list gets its own file, so earlier tables are not overwritten
csv_name <- if (length(args) >= 2) paste0("gene_panel_pct_expressing_", paste(GENES, collapse = "_"),
                                          if (length(args) >= 3) paste0("_", paste(DATASETS, collapse = "_")) else "", ".csv") else
  "gene_panel_pct_expressing_by_group.csv"
write.csv(tab, file.path(out, csv_name), row.names = FALSE)
print(tab %>% filter(group %in% c("Candidate_DPSC", "Pulp_fibroblast", "Perivascular_MCAM") | grepl("DPSC|PDLSC", group)) %>%
        group_by(gene, dataset, group) %>% summarise(cells = sum(n_cells), pct = median(pct_expressing), .groups = "drop"),
      n = 100)
message("Figures written to ", out)
