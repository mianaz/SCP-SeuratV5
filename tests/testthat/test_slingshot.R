context("Slingshot basic workflow")

library(testthat)
library(SCP)
library(Seurat)

make_small_processed <- function(ncells = 150, ngenes = 300) {
  set.seed(42)
  counts <- matrix(rpois(ncells * ngenes, lambda = 3), nrow = ngenes, ncol = ncells)
  rownames(counts) <- paste0("Gene", seq_len(ngenes))
  colnames(counts) <- paste0("Cell", seq_len(ncells))

  srt <- CreateSeuratObject(counts = counts)
  srt <- NormalizeData(srt, verbose = FALSE)
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 200, verbose = FALSE)
  srt <- ScaleData(srt, verbose = FALSE)
  srt <- RunPCA(srt, npcs = 10, verbose = FALSE)
  srt <- FindNeighbors(srt, dims = 1:10, verbose = FALSE)
  srt <- FindClusters(srt, resolution = 0.4, verbose = FALSE)
  srt <- RunUMAP(srt, dims = 1:10, verbose = FALSE)
  srt
}

test_that("RunSlingshot adds pseudotime and branch IDs", {
  skip_on_cran()
  skip_if_not_installed("slingshot")

  srt <- make_small_processed()
  expect_true("seurat_clusters" %in% colnames(srt@meta.data))

  srt2 <- RunSlingshot(srt, group.by = "seurat_clusters", reduction = "UMAP", dims = 1:2, show_plot = FALSE)
  # Expect pseudotime columns added
  pt_cols <- grep("^Pseudotime|^Lineage", colnames(srt2@meta.data), value = TRUE)
  expect_gt(length(pt_cols), 0)
  expect_true(any(grepl("BranchID$", colnames(srt2@meta.data))))
})

