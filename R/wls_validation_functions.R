# =============================================================================
# R/wls_validation_functions.R
# Shared helpers for the cross-dataset WLS validation.
# Marker panels and DPSC definitions are copied from 06_annotation.R and
# 07b_WLS_sensitivity.R in mtariqi/dpsc-wls-human-scRNAseq, so every dataset
# is analysed the same way as the discovery dataset GSE164157.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(yaml); library(dplyr)
  library(ggplot2); library(edgeR)
})

# ---- Marker panels (identical to GSE164157) ---------------------------------
STEM_CORE    <- c("THY1", "NT5E", "ENG", "NGFR", "FRZB", "LEPR", "PRRX1", "CXCL12")
PERIVASCULAR <- c("MCAM", "PDGFRB", "RGS5", "KCNJ8", "NOTCH3", "ACTA2")
FIBROBLAST   <- c("COL1A1", "DCN", "LUM", "PDGFRA")

CELLTYPE_PANELS <- list(
  Pulp_fibroblast = FIBROBLAST,
  Perivascular    = PERIVASCULAR,
  Endothelial     = c("PECAM1", "CDH5", "VWF", "CLDN5"),
  Schwann_glia    = c("PLP1", "MPZ", "S100B", "SOX10"),
  T_NK            = c("CD3E", "CD3D", "NKG7", "GZMA"),
  Myeloid         = c("CD68", "LYZ", "CD14", "C1QA"),
  B_plasma        = c("CD79A", "MS4A1", "JCHAIN", "MZB1"),
  Epithelial      = c("KRT14", "KRT5", "EPCAM"),
  Erythrocyte     = c("HBB", "HBA1")
)
MESENCHYME <- c("Pulp_fibroblast", "Perivascular")
EXCLUDED   <- c("Epithelial", "Erythrocyte")

# ---- Config and sample sheets -----------------------------------------------
load_config <- function(path = "config/validation_datasets.yaml") read_yaml(path)

read_sample_sheet <- function(path) {
  ss <- read.delim(path, stringsAsFactors = FALSE)
  if (any(grepl("TODO", unlist(ss[ss$include == "yes", ]))))
    stop("Sample sheet ", path, " still contains TODO entries. Fill it in first.")
  ss[ss$include == "yes", ]
}

target_symbol <- function(obj, cfg) {
  hit <- intersect(cfg$target_aliases, rownames(obj))
  if (length(hit) == 0) stop("Target gene not found (tried: ",
                             paste(cfg$target_aliases, collapse = ", "), ")")
  hit[1]
}

# ---- Loading count matrices of unknown GEO format ----------------------------
# Finds the files for one sample (by GSM / sample id) and reads them.
read_counts_any <- function(dir, id) {
  files <- list.files(dir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl(id, basename(files)) | grepl(id, dirname(files))]
  if (length(files) == 0) stop("No files for ", id, " under ", dir)

  mtx <- grep("matrix\\.mtx(\\.gz)?$", files, value = TRUE)
  h5  <- grep("\\.h5$", files, value = TRUE)
  tab <- grep("\\.(csv|tsv|txt)(\\.gz)?$", files, value = TRUE)
  tab <- tab[!grepl("barcodes|features|genes", basename(tab))]

  if (length(mtx) == 1) {
    bc   <- grep("barcodes\\.tsv(\\.gz)?$", files, value = TRUE)
    feat <- grep("(features|genes)\\.tsv(\\.gz)?$", files, value = TRUE)
    m <- ReadMtx(mtx = mtx, cells = bc[1], features = feat[1])
  } else if (length(h5) == 1) {
    m <- Read10X_h5(h5)
    if (is.list(m)) m <- m[["Gene Expression"]]
  } else if (length(tab) == 1) {
    df <- data.table::fread(tab, data.table = FALSE)
    genes <- df[[1]]; df <- df[, -1, drop = FALSE]
    m <- as(as.matrix(df), "dgCMatrix"); rownames(m) <- make.unique(genes)
  } else {
    stop("Could not identify a single count matrix for ", id,
         ". Files found:\n", paste(basename(files), collapse = "\n"))
  }
  m
}

# ---- QC, doublets, clustering -----------------------------------------------
qc_filter <- function(obj, cfg) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  subset(obj, nFeature_RNA > cfg$qc$min_features &
              nCount_RNA   > cfg$qc$min_counts &
              percent.mt   < cfg$qc$max_percent_mt)
}

remove_doublets <- function(obj) {
  sce <- scDblFinder::scDblFinder(as.SingleCellExperiment(obj), samples = "sample")
  obj$scDblFinder.class <- sce$scDblFinder.class
  subset(obj, scDblFinder.class == "singlet")
}

cluster_cells <- function(obj, cfg) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  red <- "pca"
  if (length(unique(obj$sample)) > 1) {
    obj <- harmony::RunHarmony(obj, group.by.vars = "sample", verbose = FALSE)
    red <- "harmony"
  }
  obj <- FindNeighbors(obj, reduction = red, dims = 1:cfg$clustering$n_pcs, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cfg$clustering$resolution, verbose = FALSE)
  RunUMAP(obj, reduction = red, dims = 1:cfg$clustering$n_pcs, verbose = FALSE)
}

# Cluster-level annotation by highest mean panel score.
# If data/annotation/<acc>_cluster_labels.csv exists (columns: cluster, cell_type),
# those confirmed labels are used instead, as in 06_annotation.R.
annotate_clusters <- function(obj, acc, out_dir) {
  for (nm in names(CELLTYPE_PANELS)) {
    g <- intersect(CELLTYPE_PANELS[[nm]], rownames(obj))
    obj <- AddModuleScore(obj, features = list(g), name = paste0("ct_", nm), seed = 2026)
  }
  sc <- grep("^ct_", colnames(obj@meta.data), value = TRUE)
  means <- aggregate(obj@meta.data[, sc], list(cluster = obj$seurat_clusters), mean)
  auto <- data.frame(cluster  = means$cluster,
                     cell_type = sub("^ct_(.*)1$", "\\1", sc[apply(means[, sc], 1, which.max)]))
  write.csv(cbind(auto, round(means[, sc], 3)),
            file.path(out_dir, "cluster_labels_auto.csv"), row.names = FALSE)

  confirmed <- file.path("data/annotation", paste0(acc, "_cluster_labels.csv"))
  lab <- if (file.exists(confirmed)) {
    message("Using confirmed labels: ", confirmed); read.csv(confirmed)
  } else {
    message("No confirmed labels yet; using automatic labels. Review ",
            file.path(out_dir, "cluster_labels_auto.csv"),
            " and save corrections to ", confirmed); auto
  }
  obj$cell_type <- lab$cell_type[match(as.character(obj$seurat_clusters), as.character(lab$cluster))]
  subset(obj, cell_type %in% EXCLUDED, invert = TRUE)
}

# ---- DPSC definitions A-D (from 07b_WLS_sensitivity.R) -----------------------
# A uses the stage-06 rule: mesenchymal sub-clusters with stem score above
# background (> 0) that are perivascular (perivascular score > 0).
# Check that these thresholds match 06_annotation.R before running.
define_dpsc <- function(obj, cfg) {
  mes <- subset(obj, cell_type %in% MESENCHYME)
  mes <- FindNeighbors(mes, reduction = if ("harmony" %in% Reductions(mes)) "harmony" else "pca",
                       dims = 1:cfg$clustering$n_pcs, verbose = FALSE)
  mes <- FindClusters(mes, resolution = cfg$clustering$mes_resolution,
                      cluster.name = "mes_subcluster", verbose = FALSE)
  stem <- intersect(STEM_CORE, rownames(mes)); peri <- intersect(PERIVASCULAR, rownames(mes))
  mes <- AddModuleScore(mes, features = list(stem), name = "stemcell", seed = 2026)
  mes <- AddModuleScore(mes, features = list(peri), name = "perivasc", seed = 2026)
  mes$stemcell <- mes$stemcell1; mes$perivasc <- mes$perivasc1

  sub_means <- aggregate(cbind(stemcell, perivasc) ~ mes_subcluster, mes@meta.data, mean)
  dpsc_subs <- sub_means$mes_subcluster[sub_means$stemcell > 0 & sub_means$perivasc > 0]

  X <- FetchData(mes, c(stem, "MCAM"))
  n_markers <- rowSums(X[, stem, drop = FALSE] > 0)
  defs <- data.frame(
    A_perivascular_subclusters = mes$mes_subcluster %in% dpsc_subs,
    B_top20pct_stem_score      = mes$stemcell >= quantile(mes$stemcell, 0.80),
    C_3plus_stem_markers       = n_markers >= 3,
    D_MCAM_positive            = X$MCAM > 0,
    row.names = colnames(mes))

  for (d in names(defs)) {
    obj[[d]] <- FALSE
    obj@meta.data[rownames(defs), d] <- defs[[d]]
  }
  list(obj = obj, subcluster_scores = sub_means, dpsc_subclusters = dpsc_subs)
}

# ---- WLS statistics ------------------------------------------------------------
# Detection and expression summaries per sample x group
wls_summary <- function(obj, gene, group_col) {
  df <- data.frame(sample = obj$sample, group = obj@meta.data[[group_col]],
                   expr = FetchData(obj, gene)[, 1])
  df %>% group_by(sample, group) %>%
    summarise(n_cells = n(), pct_WLS = 100 * mean(expr > 0),
              mean_lognorm = mean(expr), .groups = "drop")
}

# One test: is_dpsc (logical) vs is_ref (logical), within the given cells.
# >= 3 donors: edgeR QL pseudobulk ~ donor + group (as in 07_WLS_expression.R).
# < 3 donors : cell-level Wilcoxon (descriptive only; no donor replication).
wls_test <- function(obj, gene, is_dpsc, is_ref, label) {
  keep <- is_dpsc | is_ref
  grp  <- ifelse(is_dpsc[keep], "DPSC", "Ref")
  donor <- obj$donor[keep]
  counts <- LayerData(obj, layer = "counts")[, keep]
  expr   <- FetchData(obj, gene)[keep, 1]

  # pooled CPM log2FC (defined for any number of donors)
  pool <- sapply(c("DPSC", "Ref"), function(g) Matrix::rowSums(counts[, grp == g, drop = FALSE]))
  cpm  <- t(t(pool) / colSums(pool)) * 1e6
  lfc  <- log2((cpm[gene, "DPSC"] + 0.5) / (cpm[gene, "Ref"] + 0.5))

  n_donor <- length(unique(donor))
  if (n_donor >= 3) {
    key <- factor(paste(donor, grp, sep = "|"))
    pb  <- counts %*% sparse.model.matrix(~ 0 + key); colnames(pb) <- levels(key)
    meta <- data.frame(donor = sub("\\|.*", "", colnames(pb)), group = sub(".*\\|", "", colnames(pb)))
    meta$group <- factor(meta$group, levels = c("Ref", "DPSC"))
    y <- DGEList(as.matrix(pb)); y <- y[filterByExpr(y, group = meta$group), , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    # Paired (same donors in both groups) -> ~ donor + group; unpaired -> ~ group
    paired <- any(table(meta$donor) == 2)
    design <- if (paired) model.matrix(~ donor + group, meta) else model.matrix(~ group, meta)
    fit <- glmQLFit(estimateDisp(y, design), design)
    tt  <- topTags(glmQLFTest(fit, coef = "groupDPSC"), n = Inf)$table
    p   <- if (gene %in% rownames(tt)) tt[gene, "FDR"] else NA
    lfc <- if (gene %in% rownames(tt)) tt[gene, "logFC"] else lfc
    test <- sprintf("pseudobulk edgeR (%d donors, %s)", n_donor, if (paired) "paired" else "unpaired")
  } else {
    p <- wilcox.test(expr[grp == "DPSC"], expr[grp == "Ref"])$p.value
    test <- sprintf("cell-level Wilcoxon, descriptive (%d donor)", n_donor)
  }
  data.frame(comparison = label, n_DPSC = sum(grp == "DPSC"), n_ref = sum(grp == "Ref"),
             pct_WLS_DPSC = 100 * mean(expr[grp == "DPSC"] > 0),
             pct_WLS_ref  = 100 * mean(expr[grp == "Ref"] > 0),
             log2FC = unname(lfc), p_value = p, test = test)
}
