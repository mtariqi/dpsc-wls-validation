# =============================================================================
# 02_preprocess.R
# Load -> QC -> doublets -> cluster -> annotate -> DPSC definitions A-D,
# with the same settings as GSE164157 (scripts 04-07b in mtariqi/dpsc-wls-human-scRNAseq).
# Usage:  Rscript scripts/02_preprocess.R GSE185222
#         Rscript scripts/02_preprocess.R GSE202476
#         Rscript scripts/02_preprocess.R GSE227731
# =============================================================================
source("R/wls_validation_functions.R")
set.seed(2026)

acc <- commandArgs(trailingOnly = TRUE)[1]
cfg <- load_config()
if (is.na(acc) || !acc %in% names(cfg$datasets)) stop("Give one of: ", paste(names(cfg$datasets), collapse = ", "))
ds  <- cfg$datasets[[acc]]
ss  <- read_sample_sheet(ds$sample_sheet)
out <- file.path("results", acc); dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE)

# ---- Load ----------------------------------------------------------------------
if (isTRUE(ds$combined_matrix)) {
  # One count matrix for all samples (GSE227731): split cells by cell_regex
  f <- list.files(file.path("data/raw", acc), "count_matrix.*\\.(csv|tsv|txt)(\\.gz)?$",
                  recursive = TRUE, full.names = TRUE)[1]
  df <- data.table::fread(f, data.table = FALSE)
  genes <- df[[1]]; df <- df[, -1, drop = FALSE]
  m_all <- as(as.matrix(df), "dgCMatrix"); rownames(m_all) <- make.unique(genes)
  message("Example cell names: ", paste(head(colnames(m_all), 5), collapse = ", "))
  hits <- sapply(ss$cell_regex, function(r) grepl(r, colnames(m_all)))
  if (any(colSums(hits) == 0) || any(rowSums(hits) > 1))
    stop("cell_regex does not split the cells cleanly. Check the example cell names ",
         "above and edit cell_regex in ", ds$sample_sheet)
  objs <- lapply(seq_len(nrow(ss)), function(i) {
    o <- CreateSeuratObject(m_all[, hits[, i]], project = ss$sample_id[i], min.cells = 3)
    o$sample <- ss$sample_id[i]; o$donor <- ss$donor[i]; o$condition <- ss$condition[i]; o
  })
} else {
  objs <- lapply(seq_len(nrow(ss)), function(i) {
    m <- read_counts_any(file.path("data/raw", acc), ss$accession[i])
    o <- CreateSeuratObject(m, project = ss$sample_id[i], min.cells = 3)
    o$sample <- ss$sample_id[i]; o$donor <- ss$donor[i]; o$condition <- ss$condition[i]
    o
  })
}
obj <- if (length(objs) > 1) merge(objs[[1]], objs[-1], add.cell.ids = ss$sample_id) else objs[[1]]
obj <- JoinLayers(obj)
n0 <- ncol(obj)

# ---- QC and doublets -----------------------------------------------------------
obj <- qc_filter(obj, cfg); n1 <- ncol(obj)
obj <- remove_doublets(obj); n2 <- ncol(obj)
qc_tab <- data.frame(dataset = acc, cells_raw = n0, cells_after_QC = n1, singlets = n2)
write.csv(qc_tab, file.path(out, "qc_summary.csv"), row.names = FALSE); print(qc_tab)
write.csv(as.data.frame(table(sample = obj$sample)), file.path(out, "cells_per_sample.csv"), row.names = FALSE)

# ---- Clustering ----------------------------------------------------------------
obj <- cluster_cells(obj, cfg)

if (ds$design == "tissue") {
  # Tissue: annotate cell types, then define DPSCs within the mesenchyme
  obj <- annotate_clusters(obj, acc, out)
  res <- define_dpsc(obj, cfg)
  obj <- res$obj
  write.csv(res$subcluster_scores, file.path(out, "mesenchyme_subcluster_scores.csv"), row.names = FALSE)
  print(table(obj$cell_type))
  print(colSums(obj@meta.data[, grep("^[A-D]_", colnames(obj@meta.data))]))
} else {
  # Stem-cell dataset: every cell is a cultured stem cell; groups come from the
  # sample sheet (condition = DPSC / PDLSC). Sub-clusters are kept to check for
  # a WLS-high DPSC subset.
  obj$cell_type <- obj$condition
}

saveRDS(obj, file.path("data/processed", paste0(acc, "_preprocessed.rds")))

p <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE) + ggtitle(acc)
ggsave(file.path(out, "umap_cell_types.png"), p, width = 7, height = 5.5, dpi = 300)
message("Saved data/processed/", acc, "_preprocessed.rds. Review the cluster labels, then run 03_WLS_analysis.R")
