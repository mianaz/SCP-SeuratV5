context("Seurat V5 basic workflows")

test_that("Standard_SCP runs on small synthetic data", {
  skip_on_cran()
  set.seed(1)
  mat <- matrix(rpois(2000*100, lambda = 2), nrow = 2000, ncol = 100)
  rownames(mat) <- paste0("g", seq_len(nrow(mat)))
  colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  srt <- Seurat::CreateSeuratObject(counts = mat)
  srt <- Standard_SCP(
    srt,
    normalization_method = 'LogNormalize',
    linear_reduction = 'pca',
    linear_reduction_dims = 20,
    linear_reduction_dims_use = 1:10,
    nonlinear_reduction = 'umap',
    neighbor_metric = 'euclidean',
    cluster_algorithm = 'louvain',
    cluster_resolution = 0.6
  )
  expect_true("seurat_clusters" %in% colnames(srt@meta.data))
  expect_true("pca" %in% Seurat::Reductions(srt))
})

test_that("External Chromium h5 path test (optional)", {
  h5 <- Sys.getenv("SCP_TEST_10X_H5", unset = NA)
  if (is.na(h5) || !file.exists(h5)) skip("SCP_TEST_10X_H5 not set")
  mat <- Seurat::Read10X_h5(h5)
  srt <- Seurat::CreateSeuratObject(counts = mat)
  if (ncol(srt) > 1000) srt <- subset(srt, cells = sample(colnames(srt), 1000))
  srt <- Standard_SCP(srt, normalization_method = 'LogNormalize', linear_reduction = 'pca', nonlinear_reduction = 'umap')
  expect_true("seurat_clusters" %in% colnames(srt@meta.data))
})

