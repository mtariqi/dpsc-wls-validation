# ==============================================================================
# Title:        07_composite_figures.R - Composite figures in Prof. Hsu's layout
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-05
# ------------------------------------------------------------------------------
# Purpose:      For GSE164157 and GSE202476: cell-group UMAP (legend only, no
#               labels inside the plot), gene UMAP and violin by the same groups
#               (Axin2+ pulp, Axin2- pulp, other cell types); also saves every
#               panel separately.
# Inputs:       data/processed/GSE202476_preprocessed.rds; annotated GSE164157
#               object.
# Outputs:      results/figures_composite/ (composites, individual/ panels,
#               composite_pct_expressing_by_Yoshida_group_*.csv).
# Usage:        Rscript scripts/07_composite_figures.R <GSE164157 rds>
#               "GENE[:layout],..."
# Dependencies: R 4.6; Seurat 5, ggplot2, patchwork, dplyr
#               (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ------------------------------------------------------------------------------
# Details:
# Composite figure in the layout Prof. Hsu assembled for GPR177/WLS (Nov 2026):
#   one row per dataset (GSE164157, GSE202476):
#   [cell-group UMAP] [gene UMAP] [gene violin by the same groups]
#   panel tags: A / B = UMAP pairs, C / D = violins
# Changes on 2026-11-05 (requested by Prof. Hsu):
#   - no blue box in the violin plots
#   - cell-group UMAP has no labels inside the plot; groups are named in the legend only
#   - UMAP and violins use the same groups: Axin2+ pulp, Axin2- pulp, then other cell types
#   - each cell-group UMAP is also saved on its own (individual/<GEO>_cell_groups_UMAP)
# Same objects, UMAPs and Axin2+ pulp (Yoshida DPSC) definition as 05_figures_for_Hsu.R.
# Usage (from the repository root):
#   Rscript scripts/07_composite_figures.R ~/dpsc-wls-human-scRNAseq/data/processed/GSE164157_annotated.rds
#   2nd argument: genes with layouts, e.g. "PTCH1:violin,PTCH2:umap_violin,SMO:umap_violin"
# ==============================================================================

source("R/wls_validation_functions.R")
suppressPackageStartupMessages({
  library(patchwork)
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

ORDER <- c(
  "Pulp_fibroblast", "Candidate_DPSC", "Perivascular_MCAM", "Odontoblast", "Endothelial",
  "Schwann_glia", "T_NK", "Myeloid", "B_plasma"
)
PULP_TYPES <- c("Pulp_fibroblast", "MSC_progenitor", "Odontoblast")
# group names used in both the UMAP legend and the violin plots
AX_POS <- "Axin2+ pulp"
AX_NEG <- "Axin2- pulp"

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
  yg[pulp & ax] <- AX_POS
  yg[pulp & !ax] <- AX_NEG
  obj$yosh_group <- factor(yg, levels = c(
    AX_POS, AX_NEG,
    sort(setdiff(unique(yg), c(AX_POS, AX_NEG)))
  ))
  obj
}

row_panels <- function(obj, acc, g) {
  red <- umap_name(obj)
  # one fixed colour per group, so each group has the same colour in the UMAP and the violin
  pal <- setNames(scales::hue_pal()(nlevels(obj$yosh_group)), levels(obj$yosh_group))
  # cell-group UMAP: same groups as the violins, names in the legend only;
  # Axin2+ pulp cells are drawn on top so the few scattered cells stay visible
  p_annot <- DimPlot(obj,
    reduction = red, group.by = "yosh_group", label = FALSE, raster = FALSE,
    order = AX_POS, cols = pal
  ) +
    ggtitle(acc) + labs(x = "UMAP_1", y = "UMAP_2") +
    theme(legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))
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
  if (!AX_POS %in% keep && sum(obj$yosh_group == AX_POS) >= 10) keep <- c(AX_POS, keep)
  ob <- subset(obj, yosh_group %in% keep)
  ob$yosh_group <- droplevels(ob$yosh_group)
  p_vln <- VlnPlot(ob, g,
    group.by = "yosh_group", pt.size = 0.05,
    cols = unname(pal[levels(ob$yosh_group)])
  ) + NoLegend() + labs(x = NULL) +
    ggtitle(paste0(acc, ": ", g, " in Axin2+ DPSCs")) +
    theme(plot.title = element_text(size = 10), axis.text.x = element_text(size = 7))
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
    ggsave(file.path(ind, paste0(rr[[2]], "_cell_groups_UMAP.png")), rr[[1]][[1]],
      width = 7, height = 5, dpi = 300, bg = "white"
    )
    ggsave(file.path(ind, paste0(rr[[2]], "_cell_groups_UMAP.pdf")), rr[[1]][[1]],
      width = 7, height = 5, bg = "white"
    )
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
print(tab %>% filter(group %in% c(AX_POS, AX_NEG)), n = 50)
