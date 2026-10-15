# =============================================================================
# 05_figures_for_Hsu.R
# Figures requested by Prof. Hsu (email, 14 Oct 2026), for every dataset:
#   1. UMAP with each cluster labelled by cell type (candidate DPSCs labelled)
#   2. Marker dot plot that justifies each label (WLS added as the last column)
#   3. WLS FeaturePlot next to the labelled UMAP
#   4. WLS violin plot by cell type
#   5. One row of WLS UMAPs for all four datasets (incl. GSE164157)
# Usage (from the repository root):
#   Rscript scripts/05_figures_for_Hsu.R ~/dpsc-wls-human-scRNAseq/data/processed/<GSE164157 annotated>.rds
# The argument is the annotated GSE164157 Seurat object from the discovery repo.
# =============================================================================
source("R/wls_validation_functions.R")
suppressPackageStartupMessages(library(patchwork))
set.seed(2026)

args <- commandArgs(trailingOnly = TRUE)
discovery_rds <- if (length(args) >= 1) args[1] else NA
out <- "results/figures_for_Hsu"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

MARKERS <- list(
  Pulp_fibroblast   = c("COL1A1", "DCN", "LUM", "PDGFRA"),
  # AXIN2: Axin2+ DPSCs (Yoshida et al. 2025, J Dent Res). Shown in the dot plot only;
  # the DPSC definition (stem score of the 8 genes) is unchanged.
  Candidate_DPSC    = c("AXIN2", "THY1", "NT5E", "ENG", "NGFR", "FRZB", "LEPR", "PRRX1", "CXCL12"),
  Perivascular_MCAM = c("MCAM", "PDGFRB", "RGS5", "ACTA2"),
  Odontoblast       = c("DSPP", "DMP1", "PHEX"),
  Endothelial       = c("PECAM1", "CDH5", "VWF"),
  Schwann_glia      = c("PLP1", "MPZ", "SOX10"),
  T_NK              = c("CD3E", "CD3D", "NKG7"),
  Myeloid           = c("CD68", "LYZ", "C1QA"),
  B_plasma          = c("CD79A", "MZB1", "JCHAIN")
)
ORDER <- names(MARKERS)

# Marker panels from Yoshida et al. 2025 (J Dent Res; Fig. 1, Appendix Figs 1-2), human orthologs.
# Mouse-only Slfn3 is omitted.
YOSHIDA <- list(
  Macrophage    = c("AIF1", "C1QA", "C1QC", "CSF1R"),
  NK_T          = c("DPP4", "ITGB7", "CD37", "CCR2"),
  Glial         = c("PLP1", "ANK3", "MAL", "GJC3", "CADM4"),
  Perivascular  = c("LMOD1", "ACTA2", "MYH11", "TAGLN", "MYL9"),
  Pulp_Nestin   = c("NES"),
  Odontoblast   = c("DCN", "COL1A1", "MMP2", "TNC", "COCH", "DMP1"),
  PC2           = c("STC1", "CPM", "FAM167A", "TIAM2", "DES"),
  PC1           = c("ADGRL4", "ANKRD33B", "RBP7", "SH2D3C", "ALOX12"),
  DPSC1         = c("AXIN2", "BMP4", "BMP6", "BMX", "FBLN2"),
  DPSC2         = c("CDH4", "ACKR1", "SELP")
)
PULP_TYPES <- c("Pulp_fibroblast", "MSC_progenitor", "Odontoblast")   # non-perivascular mesenchyme

# Yoshida definition: DPSCs are Axin2+ pulp (Nestin+) cells that are NOT perivascular.
# Here: AXIN2 count > 0 among non-perivascular mesenchymal cells.
yoshida_figures <- function(obj, acc, subtitle, gene, min_cells) {
  ax <- if ("AXIN2" %in% rownames(obj)) FetchData(obj, "AXIN2")[, 1] > 0 else rep(FALSE, ncol(obj))
  pulp <- as.character(obj$cell_type) %in% PULP_TYPES
  grp <- as.character(obj$cell_type)
  grp[pulp & ax]  <- "AXIN2+ pulp (Yoshida DPSC)"
  grp[pulp & !ax] <- "AXIN2- pulp cells"
  obj$yosh_group <- grp
  lv <- c("AXIN2+ pulp (Yoshida DPSC)", "AXIN2- pulp cells",
          sort(setdiff(unique(grp), c("AXIN2+ pulp (Yoshida DPSC)", "AXIN2- pulp cells"))))
  obj$yosh_group <- factor(grp, levels = lv)

  expr <- FetchData(obj, gene)[, 1]
  tab <- data.frame(sample = obj$sample, group = obj$yosh_group, wls = expr) %>%
    group_by(sample, group) %>%
    summarise(n_cells = n(), pct_WLS = round(100 * mean(wls > 0), 1),
              mean_lognorm = round(mean(wls), 3), .groups = "drop")
  write.csv(tab, file.path(out, paste0(acc, "_WLS_by_Yoshida_groups.csv")), row.names = FALSE)

  keep <- names(which(table(obj$yosh_group) >= min_cells))
  oy <- subset(obj, yosh_group %in% keep)
  oy$yosh_group <- droplevels(oy$yosh_group)
  feats <- unique(unlist(YOSHIDA)); feats <- c(feats[feats %in% rownames(obj)], gene)
  p_dot <- DotPlot(oy, features = unique(feats), group.by = "yosh_group") + RotatedAxis() +
    ggtitle(paste0(acc, ": Yoshida et al. 2025 markers (human orthologs); WLS last")) +
    theme(axis.text.x = element_text(size = 7)) + labs(x = NULL, y = NULL)
  p_vln <- VlnPlot(oy, gene, group.by = "yosh_group", pt.size = 0.05) + NoLegend() +
    ggtitle(paste0(gene, " (Gpr177) - ", acc, ": AXIN2+ pulp cells (Yoshida DPSC) and other cells")) +
    labs(x = NULL) + theme(plot.title = element_text(size = 11))
  save_both(p_dot, paste0(acc, "_5_Yoshida_marker_dotplot"), 14, 4.8)
  save_both(p_vln, paste0(acc, "_6_WLS_violin_Yoshida_groups"), 8, 4.8)
  print(tab %>% group_by(group) %>% summarise(cells = sum(n_cells), median_pct_WLS = median(pct_WLS)))
}

umap_name <- function(obj) {
  r <- grep("umap", Reductions(obj), value = TRUE, ignore.case = TRUE)
  r <- r[!grepl("mes", r)]
  if (length(r) == 0) stop("No UMAP reduction found")
  r[1]
}

save_both <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}

# Build the plotting group: cell type, with candidate DPSCs as their own label
prepare <- function(obj, dpsc_col) {
  grp <- as.character(obj$cell_type)
  grp[obj@meta.data[[dpsc_col]] %in% TRUE] <- "Candidate_DPSC"
  keep_levels <- c(intersect(ORDER, unique(grp)), setdiff(unique(grp), ORDER))
  obj$plot_group <- factor(grp, levels = keep_levels)
  # groups with < 50 cells are kept on the UMAP but left out of dot/violin plots
  obj
}

# One text label per cluster (at the cluster's median position), named by the cell type
# most cells in that cluster carry, so labels sit on the clusters as in published UMAPs
cluster_labels <- function(obj, red) {
  clcol <- intersect(c("harmony_clusters", "seurat_clusters"), colnames(obj@meta.data))[1]
  emb <- Embeddings(obj, red)[, 1:2]
  df <- data.frame(x = emb[, 1], y = emb[, 2], group = as.character(obj$plot_group),
                   cl = as.character(obj@meta.data[[clcol]]))
  min_n <- max(15, round(0.005 * nrow(df)))
  lab <- df %>% group_by(cl, group) %>%
    summarise(n = n(), x = median(x), y = median(y), .groups = "drop") %>%
    filter(n >= min_n)
  # one label per cell type per cluster; drop repeats that sit on top of each other
  lab %>% group_by(group) %>% arrange(desc(n), .by_group = TRUE) %>%
    filter(!duplicated(paste(round(x / 2), round(y / 2)))) %>% ungroup()
}

tissue_figures <- function(obj, acc, subtitle, gene) {
  red <- umap_name(obj)
  min_cells <- if (ncol(obj) > 5000) 50 else 20   # small sound-tooth sample keeps 20
  big <- names(which(table(obj$plot_group) >= min_cells))
  obj_big <- subset(obj, plot_group %in% big)
  obj_big$plot_group <- droplevels(obj_big$plot_group)

  lab <- cluster_labels(obj, red)
  p_annot <- DimPlot(obj, reduction = red, group.by = "plot_group", label = FALSE, raster = FALSE) +
    ggrepel::geom_text_repel(data = lab, aes(x = x, y = y, label = group), inherit.aes = FALSE,
                             size = 3.3, fontface = "bold", box.padding = 0.4, min.segment.length = 0,
                             max.overlaps = Inf, seed = 2026) +
    ggtitle(acc, subtitle = subtitle) + theme(legend.text = element_text(size = 9))
  # perivascular markers are always shown when candidate DPSCs are present, because
  # candidate DPSCs sit in the perivascular population
  show <- intersect(ORDER, c(levels(obj_big$plot_group),
                             if ("Candidate_DPSC" %in% levels(obj_big$plot_group)) "Perivascular_MCAM"))
  feats <- unique(c(unlist(MARKERS[show]), gene))
  feats <- feats[feats %in% rownames(obj)]
  p_dot <- DotPlot(obj_big, features = feats, group.by = "plot_group") + RotatedAxis() +
    ggtitle(paste0(acc, ": markers used to assign cell types (", gene, " last)")) +
    theme(axis.text.x = element_text(size = 8)) + labs(x = NULL, y = NULL)
  p_wls <- FeaturePlot(obj, gene, reduction = red, order = TRUE, raster = FALSE) +
    ggtitle(paste0(gene, " (Gpr177) - ", acc), subtitle = subtitle)
  p_vln <- VlnPlot(obj_big, gene, group.by = "plot_group", pt.size = 0.05) + NoLegend() +
    ggtitle(paste0(gene, " (Gpr177) - ", acc)) + labs(x = NULL)

  save_both(p_annot, paste0(acc, "_1_UMAP_cell_types"), 7.5, 6)
  save_both(p_dot, paste0(acc, "_2_marker_dotplot"), 12, 4.5)
  save_both(p_annot + p_wls, paste0(acc, "_3_UMAP_cell_types_and_WLS"), 14, 6)
  save_both(p_vln, paste0(acc, "_4_WLS_violin"), 7, 4.5)
  yoshida_figures(obj, acc, subtitle, gene, min_cells)
  p_wls
}

wls_row <- list()

# ---- 1. Discovery dataset GSE164157 ---------------------------------------------
if (!is.na(discovery_rds)) {
  d <- readRDS(path.expand(discovery_rds))
  if (!"cell_type" %in% colnames(d@meta.data) || !"candidate_DPSC" %in% colnames(d@meta.data))
    stop("The GSE164157 object needs 'cell_type' and 'candidate_DPSC' columns (saved by 06_annotation.R)")
  d <- subset(d, cell_type %in% c("Epithelial", "Erythrocyte", "Unassigned"), invert = TRUE)
  d <- prepare(d, "candidate_DPSC")
  p_donor <- DimPlot(d, reduction = umap_name(d), group.by = "sample", shuffle = TRUE, raster = FALSE) +
    ggtitle("GSE164157: cells coloured by donor", subtitle = "Clusters contain cells from all donors after Harmony")
  save_both(p_donor, "GSE164157_0_UMAP_by_donor", 7.5, 6)
  wls_row[["GSE164157"]] <- tissue_figures(d, "GSE164157", paste0("Fresh pulp, 5 donors, ", ncol(d), " cells"),
                                           target_symbol(d, load_config()))
  rm(d); gc()
} else {
  message("No GSE164157 object given: skipping the discovery dataset")
}

# ---- 2. Validation tissue datasets -------------------------------------------------
cfg <- load_config()
for (acc in c("GSE202476", "GSE185222")) {
  obj <- readRDS(file.path("data/processed", paste0(acc, "_preprocessed.rds")))
  ds <- cfg$datasets[[acc]]
  if (ds$primary_condition %in% obj$condition) obj <- subset(obj, condition == ds$primary_condition)
  obj <- prepare(obj, "A_stem_subclusters")
  sub <- if (acc == "GSE202476") "Young tooth" else "Sound tooth (DTP_01)"
  wls_row[[acc]] <- tissue_figures(obj, acc, paste0(sub, ", ", ncol(obj), " cells"), target_symbol(obj, cfg))
  rm(obj); gc()
}

# ---- 3. GSE227731: DPSCs vs PDLSCs, with DPSC sub-clusters ---------------------------
obj <- readRDS("data/processed/GSE227731_preprocessed.rds")
gene <- target_symbol(obj, cfg)
red <- umap_name(obj)
obj$plot_group <- ifelse(obj$cell_type == "DPSC", paste0("DPSC_sub", obj$seurat_clusters), "PDLSC")
big <- names(which(table(obj$plot_group) >= 20))
lv <- c(sort(grep("^DPSC", big, value = TRUE)), "PDLSC")
obj_big <- subset(obj, plot_group %in% big)
obj_big$plot_group <- factor(obj_big$plot_group, levels = lv)
p_annot <- DimPlot(obj, reduction = red, group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("GSE227731", subtitle = paste0("DPSCs + PDLSCs, ", ncol(obj), " cells"))
p_wls <- FeaturePlot(obj, gene, reduction = red, order = TRUE) +
  ggtitle(paste0(gene, " (Gpr177) - GSE227731"), subtitle = paste0("DPSCs + PDLSCs, ", ncol(obj), " cells"))
p_vln <- VlnPlot(obj_big, gene, group.by = "plot_group", pt.size = 0.05) + NoLegend() +
  ggtitle(paste0(gene, " (Gpr177) - GSE227731")) + labs(x = NULL)
save_both(p_annot, "GSE227731_1_UMAP_cell_types", 7.5, 6)
p_sub <- DimPlot(obj, reduction = red, group.by = "plot_group", label = TRUE, repel = TRUE) +
  ggtitle("GSE227731: DPSC sub-clusters and PDLSCs", subtitle = paste0(ncol(obj), " cells"))
save_both(p_sub, "GSE227731_1b_UMAP_DPSC_subclusters", 7.5, 6)
feats227 <- unique(c(MARKERS$Candidate_DPSC, MARKERS$Perivascular_MCAM, MARKERS$Pulp_fibroblast,
                     "POSTN", "ASPN", "MKI67", "TOP2A", gene))
feats227 <- feats227[feats227 %in% rownames(obj_big)]
p_dot227 <- DotPlot(obj_big, features = feats227, group.by = "plot_group") + RotatedAxis() +
  ggtitle("GSE227731: stem, perivascular, fibroblast, PDL (POSTN, ASPN) and cycling (MKI67, TOP2A) markers; WLS last") +
  theme(axis.text.x = element_text(size = 8), plot.title = element_text(size = 10)) + labs(x = NULL, y = NULL)
save_both(p_dot227, "GSE227731_2_marker_dotplot", 12, 4)
ax227 <- if ("AXIN2" %in% rownames(obj)) FetchData(obj, "AXIN2")[, 1] > 0 else rep(FALSE, ncol(obj))
obj$yosh_group <- factor(ifelse(obj$cell_type == "DPSC", ifelse(ax227, "AXIN2+ DPSC", "AXIN2- DPSC"), "PDLSC"),
                         levels = c("AXIN2+ DPSC", "AXIN2- DPSC", "PDLSC"))
tab227 <- data.frame(group = obj$yosh_group, wls = FetchData(obj, gene)[, 1]) %>% group_by(group) %>%
  summarise(n_cells = n(), pct_WLS = round(100 * mean(wls > 0), 1), mean_lognorm = round(mean(wls), 3))
write.csv(tab227, file.path(out, "GSE227731_WLS_by_Yoshida_groups.csv"), row.names = FALSE)
print(tab227)
p_vy <- VlnPlot(obj, gene, group.by = "yosh_group", pt.size = 0.05) + NoLegend() +
  ggtitle(paste0(gene, " (Gpr177) - GSE227731: AXIN2+ vs AXIN2- DPSCs and PDLSCs")) + labs(x = NULL)
save_both(p_vy, "GSE227731_6_WLS_violin_Yoshida_groups", 6.5, 4.5)
save_both(p_annot + p_wls, "GSE227731_3_UMAP_cell_types_and_WLS", 14, 6)
save_both(p_vln, "GSE227731_4_WLS_violin", 7, 4.5)
wls_row[["GSE227731"]] <- p_wls
rm(obj, obj_big); gc()

# ---- 4. One row of WLS UMAPs for all datasets --------------------------------------
row <- wrap_plots(wls_row, nrow = 1)
save_both(row, "All_datasets_WLS_UMAP_row", 5.5 * length(wls_row), 5)
message("Figures written to ", out)
print(list.files(out))
