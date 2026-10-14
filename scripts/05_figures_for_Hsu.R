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
  Candidate_DPSC    = c("THY1", "NT5E", "ENG", "NGFR", "FRZB", "LEPR", "PRRX1", "CXCL12"),
  Perivascular_MCAM = c("MCAM", "PDGFRB", "RGS5", "ACTA2"),
  Odontoblast       = c("DSPP", "DMP1", "PHEX"),
  Endothelial       = c("PECAM1", "CDH5", "VWF"),
  Schwann_glia      = c("PLP1", "MPZ", "SOX10"),
  T_NK              = c("CD3E", "CD3D", "NKG7"),
  Myeloid           = c("CD68", "LYZ", "C1QA"),
  B_plasma          = c("CD79A", "MZB1", "JCHAIN")
)
ORDER <- names(MARKERS)

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
  # groups with < 20 cells are kept on the UMAP but left out of dot/violin plots
  obj
}

tissue_figures <- function(obj, acc, subtitle, gene) {
  red <- umap_name(obj)
  big <- names(which(table(obj$plot_group) >= 20))
  obj_big <- subset(obj, plot_group %in% big)
  obj_big$plot_group <- droplevels(obj_big$plot_group)

  p_annot <- DimPlot(obj, reduction = red, group.by = "plot_group", label = TRUE, repel = TRUE,
                     label.size = 3.5, raster = FALSE) +
    ggtitle(acc, subtitle = subtitle) + theme(legend.text = element_text(size = 9))
  feats <- unique(c(unlist(MARKERS[intersect(ORDER, levels(obj_big$plot_group))]), gene))
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
save_both(p_annot + p_wls, "GSE227731_3_UMAP_cell_types_and_WLS", 14, 6)
save_both(p_vln, "GSE227731_4_WLS_violin", 7, 4.5)
wls_row[["GSE227731"]] <- p_wls
rm(obj, obj_big); gc()

# ---- 4. One row of WLS UMAPs for all datasets --------------------------------------
row <- wrap_plots(wls_row, nrow = 1)
save_both(row, "All_datasets_WLS_UMAP_row", 5.5 * length(wls_row), 5)
message("Figures written to ", out)
print(list.files(out))
