context("PAGA and scVelo smoke tests (optional)")

library(testthat)
library(SCP)
library(Seurat)

py_mod_available <- function(mod) {
  if (!reticulate::py_available(initialize = FALSE)) return(FALSE)
  reticulate::py_module_available(mod)
}

make_clustered <- function(ncells = 300, ngenes = 500) {
  set.seed(7)
  counts <- matrix(rpois(ncells * ngenes, lambda = 3), nrow = ngenes, ncol = ncells)
  rownames(counts) <- paste0("Gene", seq_len(ngenes))
  colnames(counts) <- paste0("Cell", seq_len(ncells))
  srt <- CreateSeuratObject(counts = counts)
  srt <- Standard_SCP(srt,
                      normalization_method = "LogNormalize",
                      linear_reduction = "pca",
                      linear_reduction_dims = 20,
                      nonlinear_reduction = "umap",
                      nonlinear_reduction_dims = 2,
                      cluster_resolution = 0.6,
                      seed = 11)
  srt
}

test_that("RunPAGA executes when scanpy/anndata are available", {
  skip_on_cran()
  if (!py_mod_available("scanpy") || !py_mod_available("anndata")) {
    skip("scanpy/anndata not available; skipping PAGA smoke test")
  }
  srt <- make_clustered(ncells = 200, ngenes = 300)
  expect_no_error({
    srt2 <- RunPAGA(srt = srt, assay_X = "RNA", group_by = "seurat_clusters",
                    linear_reduction = "PCA", nonlinear_reduction = "UMAP",
                    return_seurat = TRUE, show_plot = FALSE)
    expect_s4_class(srt2, "Seurat")
    expect_true("paga" %in% names(srt2@misc))
  })
})

test_that("EnsureEnv works for scvelo when available", {
  skip_on_cran()
  if (!py_mod_available("scvelo") || !py_mod_available("anndata")) {
    skip("scvelo/anndata not available; skipping velocity env test")
  }
  expect_no_error({ EnsureEnv(required = c("scvelo", "anndata")) })
})
