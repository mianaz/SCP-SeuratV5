context("Spatial plotting fallbacks (metadata coords)")

library(testthat)
library(SCP)
library(Seurat)

create_spatial_like_srt <- function(ncells = 200, ngenes = 300) {
  set.seed(123)
  counts <- matrix(rpois(ncells * ngenes, lambda = 2), nrow = ngenes, ncol = ncells)
  rownames(counts) <- paste0("Gene", seq_len(ngenes))
  colnames(counts) <- paste0("Cell", seq_len(ncells))

  srt <- CreateSeuratObject(counts = counts)
  # add simple coords (no Images)
  srt$spatial_x <- runif(ncells, min = 0, max = 100)
  srt$spatial_y <- runif(ncells, min = 0, max = 100)
  # Add cluster labels
  srt$seurat_clusters <- sample(letters[1:4], size = ncells, replace = TRUE)

  srt <- NormalizeData(srt, verbose = FALSE)
  srt
}

test_that("SpatialFeaturePlot falls back to coords-only plotting", {
  skip_on_cran()
  srt <- create_spatial_like_srt()
  # choose a gene that exists post-normalization
  feat <- rownames(srt)[1]
  expect_no_error({
    p <- SpatialFeaturePlot(srt, features = feat, combine = TRUE)
    expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
  })
})

test_that("SpatialClusterPlot falls back to coords-only plotting", {
  skip_on_cran()
  srt <- create_spatial_like_srt()
  expect_no_error({
    p <- SpatialClusterPlot(srt, group.by = "seurat_clusters")
    expect_true(inherits(p, "ggplot"))
  })
})

test_that("SpatialDomainPlot falls back to coords-only plotting", {
  skip_on_cran()
  srt <- create_spatial_like_srt()
  # create a domain column
  srt$spatial_domain <- sample(LETTERS[1:3], size = ncol(srt), replace = TRUE)
  expect_no_error({
    p <- SpatialDomainPlot(srt, domain_key = "spatial_domain")
    expect_true(inherits(p, "ggplot"))
  })
})

