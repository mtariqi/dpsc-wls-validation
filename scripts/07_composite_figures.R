# ==============================================================================
# Title:        07_composite_figures.R - Composite figures in Prof. Hsu's layout
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      For GSE164157 and GSE202476: labelled UMAP, gene UMAP and
#               violin by Yoshida group with AXIN2+ pulp cells boxed; also
#               saves every panel separately.
# Inputs:       data/processed/GSE202476_preprocessed.rds; annotated GSE164157
#               object.
# Outputs:      results/figures_composite/ (composites, individual/ panels,
#               composite_pct_expressing_by_Yoshida_group_*.csv).
# Usage:        Rscript scripts/07_composite_figures.R <GSE164157 rds>
#               "GENE[:layout],..."
# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,
#               patchwork, dplyr (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ------------------------------------------------------------------------------
# Details:
# Composite figure in the layout Prof. Hsu assembled for GPR177/WLS (Nov 2026):
#   one row per dataset (GSE164157, GSE202476):
#   [labelled cell-type UMAP] [gene UMAP] [gene violin by Yoshida groups, AXIN2+ pulp boxed]
#   panel tags: A / B = UMAP pairs, C / D = violins
# Same objects, UMAPs, labels and AXIN2+ pulp (Yoshida DPSC) definition as 05_figures_for_Hsu.R.
# Usage (from the repository root):
#   Rscript scripts/07_composite_figures.R ~/dpsc-wls-human-scRNAseq/data/processed/GSE164157_annotated.rds
#   optional 2nd argument: genes, e.g. "PTCH2,SMO" (default "PTCH1,PTCH2,SMO,WLS")
# ==============================================================================

source("R/wls_validation_functions.R")
suppressPackageStartupMessages({
  library(patchwork)
  library(ggrepel)
})
set.seed(2026)

args <- commandArgs(trailingOnly = TRUE)
discovery_rds <- args[1]
# 2nd argument: genes, each optionally with a layout after ":"
#   full        = labelled UMAP + gene UMAP + violin per dataset (Prof. Hsu's WLS layout; default)
#   umap_violin = gene UMAP + violin per dataset
#   violin      = violins only
# e.g. "PTCH1:violin,PTCH2:umap_violin,SMO:umap_violin"
spec <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("PTCH1", "PTCH2", "SMO", "WLS")
GENES <- toupper(sub(":.*", "", spec))
LAYOUT <- setNames(ifelse(grepl(":", spec), tolower(sub(".*:", "", spec)), "full"), GENES)
stopifnot(all(LAYOUT %in% c("full", "umap_violin", "violin")))
message(paste(GENES, LAYOUT, sep = " -> ", collapse = "; "))
out <- "results/figures_composite"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cfg <- load_config()

ORDER <- c(
  "Pulp_fibroblast", "Candidate_DPSC", "Perivascular_MCAM", "Odontoblast", "Endothelial",
  "Schwann_glia", "T_NK", "Myeloid", "B_plasma"
)
PULP_TYPES <- c("Pulp_fibroblast", "MSC_progenitor", "Odontoblast")

umap_name <- function(obj) {
  r <- grep("umap", Reductions(obj), value = TRUE, ignore.case = TRUE)
  r <- r[!grepl("mes", r)]
  pref <- r[grepl("harmony", r, ignore.case = TRUE)]
  if (length(pref)) pref[1] else r[1]
}

prepare <- function(obj, dpsc_col) {
  # cell-type groups for the labelled UMAP (candidate DPSCs labelled), as in 05_figures_for_Hsu.R
  grp <- as.character(obj$cell_type)
  grp[obj@meta.data[[dpsc_col]] %in% TRUE] <- "Candidate_DPSC"
  obj$plot_group <- factor(grp, levels = c(intersect(ORDER, unique(grp)), setdiff(unique(grp), ORDER)))
  # Yoshida groups for the violin: AXIN2+ / AXIN2- among non-perivascular pulp cells
  ax <- if ("AXIN2" %in% rownames(obj)) FetchData(obj, "AXIN2")[, 1] > 0 else rep(FALSE, ncol(obj))
  pulp <- as.character(obj$cell_type) %in% PULP_TYPES
  yg <- as.character(obj$cell_type)
  yg[pulp & ax] <- "AXIN2+ pulp"
  yg[pulp & !ax] <- "AXIN2- pulp"
  obj$yosh_group <- factor(yg, levels = c(
    "AXIN2+ pulp", "AXIN2- pulp",
    sort(setdiff(unique(yg), c("AXIN2+ pulp", "AXIN2- pulp")))
  ))
  obj
}

cluster_labels <- function(obj, red) {
  clcol <- intersect(c("harmony_clusters", "seurat_clusters"), colnames(obj@meta.data))[1]
  emb <- Embeddings(obj, red)[, 1:2]
  df <- data.frame(
    x = emb[, 1], y = emb[, 2], group = as.character(obj$plot_group),
    cl = as.character(obj@meta.data[[clcol]])
  )
  min_n <- max(15, round(0.005 * nrow(df)))
  df %>%
    group_by(cl, group) %>%
    summarise(n = n(), x = median(x), y = median(y), .groups = "drop") %>%
    filter(n >= min_n) %>%
    group_by(group) %>%
    arrange(desc(n), .by_group = TRUE) %>%
    filter(!duplicated(paste(round(x / 2), round(y / 2)))) %>%
    ungroup()
}

row_panels <- function(obj, acc, g) {
  red <- umap_name(obj)
  lab <- cluster_labels(obj, red)
  p_annot <- DimPlot(obj, reduction = red, group.by = "plot_group", label = FALSE, raster = FALSE) +
    geom_text_repel(
      data = lab, aes(x = x, y = y, label = group), inherit.aes = FALSE, size = 2.6,
      fontface = "bold", box.padding = 0.3, min.segment.length = 0, max.overlaps = Inf, seed = 2026
    ) +
    ggtitle(acc) + labs(x = "UMAP_1", y = "UMAP_2") +
    theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))
  if (!g %in% rownames(obj)) {
    blank <- ggplot() +
      annotate("text", x = 0, y = 0, label = paste(g, "is not in the", acc, "count matrix")) +
      theme_void()
    return(list(p_annot, blank, blank))
  }
  x <- FetchData(obj, g)[, 1]
  p_feat <- if (max(x) == 0) {
    FeaturePlot(obj, g, reduction = red, raster = FALSE, cols = c("lightgrey", "lightgrey")) + NoLegend() +
      ggtitle(paste0(acc, ": ", g), subtitle = "not detected in any cell") + labs(x = "UMAP_1", y = "UMAP_2")
  } else {
    FeaturePlot(obj, g, reduction = red, order = TRUE, raster = FALSE) + ggtitle(paste0(acc, ": ", g)) +
      labs(x = "UMAP_1", y = "UMAP_2")
  }
  min_cells <- if (ncol(obj) > 5000) 50 else 20
  keep <- names(which(table(obj$yosh_group) >= min_cells))
  if (!"AXIN2+ pulp" %in% keep && sum(obj$yosh_group == "AXIN2+ pulp") >= 10) keep <- c("AXIN2+ pulp", keep)
  ob <- subset(obj, yosh_group %in% keep)
  ob$yosh_group <- droplevels(ob$yosh_group)
  p_vln <- VlnPlot(ob, g, group.by = "yosh_group", pt.size = 0.05) + NoLegend() + labs(x = NULL) +
    ggtitle(paste0(acc, ": ", g, " in Axin2+ DPSCs")) +
    theme(plot.title = element_text(size = 10), axis.text.x = element_text(size = 7))
  if ("AXIN2+ pulp" %in% levels(ob$yosh_group)) {
    p_vln <- p_vln + annotate("rect",
      xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
      fill = NA, colour = "#4C7BE5", linewidth = 0.9
    )
  }
  list(p_annot, p_feat, p_vln)
}

# ---- load the two datasets once ---------------------------------------------------
d <- readRDS(path.expand(discovery_rds))
d <- subset(d, cell_type %in% c("Epithelial", "Erythrocyte", "Unassigned"), invert = TRUE)
d <- prepare(d, "candidate_DPSC")
y <- readRDS("data/processed/GSE202476_preprocessed.rds")
y <- prepare(y, "A_stem_subclusters")

tab <- list()
for (g in GENES) {
  r1 <- row_panels(d, "GSE164157", g)
  r2 <- row_panels(y, "GSE202476", g)
  lay <- LAYOUT[[g]]
  if (lay == "full") {
    fig <- wrap_plots(c(r1, r2), ncol = 3, widths = c(1.15, 1, 1)) +
      plot_annotation(tag_levels = list(c("A", "", "C", "B", "", "D")))
    w <- 16
    h <- 9
  } else if (lay == "umap_violin") {
    fig <- wrap_plots(list(r1[[2]], r1[[3]], r2[[2]], r2[[3]]), ncol = 2, widths = c(1, 1.1)) +
      plot_annotation(tag_levels = list(c("A", "C", "B", "D")))
    w <- 11
    h <- 9
  } else {
    fig <- wrap_plots(list(r1[[3]], r2[[3]]), ncol = 2) + plot_annotation(tag_levels = "A")
    w <- 11
    h <- 4.5
  }
  # also save every panel on its own, so each graph can be copied separately
  ind <- file.path(out, "individual")
  dir.create(ind, showWarnings = FALSE)
  for (rr in list(list(r1, "GSE164157"), list(r2, "GSE202476"))) {
    if (lay %in% c("full", "umap_violin")) {
      ggsave(file.path(ind, paste0(g, "_", rr[[2]], "_UMAP.png")), rr[[1]][[2]], width = 6, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(ind, paste0(g, "_", rr[[2]], "_UMAP.pdf")), rr[[1]][[2]], width = 6, height = 5, bg = "white")
    }
    ggsave(file.path(ind, paste0(g, "_", rr[[2]], "_violin.png")), rr[[1]][[3]], width = 7, height = 4.5, dpi = 300, bg = "white")
    ggsave(file.path(ind, paste0(g, "_", rr[[2]], "_violin.pdf")), rr[[1]][[3]], width = 7, height = 4.5, bg = "white")
  }
  fig <- fig & theme(plot.tag = element_text(face = "bold", size = 14))
  stem <- paste0(g, "_", lay, "_GSE164157_GSE202476")
  ggsave(file.path(out, paste0(stem, ".png")), fig, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out, paste0(stem, ".pdf")), fig, width = w, height = h, bg = "white")
  for (o in list(list(d, "GSE164157"), list(y, "GSE202476"))) {
    obj <- o[[1]]
    if (!g %in% rownames(obj)) next
    x <- FetchData(obj, g)[, 1]
    tab[[length(tab) + 1]] <- data.frame(gene = g, dataset = o[[2]], group = obj$yosh_group, x = x) %>%
      group_by(gene, dataset, group) %>%
      summarise(
        n_cells = n(), pct_expressing = round(100 * mean(x > 0), 1), mean_lognorm = round(mean(x), 3),
        .groups = "drop"
      )
  }
  message("Saved ", stem)
}
tab <- do.call(rbind, tab)
write.csv(tab, file.path(out, paste0("composite_pct_expressing_by_Yoshida_group_", paste(GENES, collapse = "_"), ".csv")),
  row.names = FALSE
)
print(tab %>% filter(group %in% c("AXIN2+ pulp", "AXIN2- pulp")), n = 50)
