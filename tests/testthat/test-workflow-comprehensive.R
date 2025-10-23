# Basic workflow tests for SCPNext package
library(testthat)
library(SCPNext)
library(Seurat)

# Helper to create simple test data
create_test_data <- function(ncells = 100, ngenes = 200) {
  set.seed(42)
  counts <- matrix(
    rpois(ncells * ngenes, lambda = 5),
    nrow = ngenes,
    ncol = ncells,
    dimnames = list(
      paste0("Gene", 1:ngenes),
      paste0("Cell", 1:ncells)
    )
  )

  metadata <- data.frame(
    CellType = sample(c("TypeA", "TypeB", "TypeC"), ncells, replace = TRUE),
    Batch = sample(c("Batch1", "Batch2"), ncells, replace = TRUE),
    row.names = colnames(counts)
  )

  CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    project = "TestProject"
  )
}

# Test Standard_SCP basic functionality
test_that("Standard_SCP runs without errors", {
  skip_on_cran()

  srt <- create_test_data(ncells = 100, ngenes = 200)

  # Run with minimal parameters
  srt_processed <- Standard_SCP(
    srt = srt,
    normalization_method = "LogNormalize",
    nHVF = 50,
    linear_reduction = "pca",
    linear_reduction_dims = 10,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = 2,
    cluster_resolution = 0.5,
    seed = 42
  )

  # Basic checks
  expect_s4_class(srt_processed, "Seurat")
  expect_gt(length(VariableFeatures(srt_processed)), 0)
  expect_true(length(Reductions(srt_processed)) > 0)
})

# Test RunCellQC basic functionality
test_that("RunCellQC adds quality metrics", {
  skip_on_cran()

  srt <- create_test_data(ncells = 100, ngenes = 200)

  # Add mock mitochondrial genes
  rownames(srt)[1:5] <- paste0("MT-", rownames(srt)[1:5])

  srt_qc <- RunCellQC(
    srt = srt,
    return_filtered = FALSE,
    qc_metrics = c("umi", "gene", "mito")
  )

  # Check QC columns were added
  expect_true("CellQC" %in% colnames(srt_qc@meta.data))
  expect_true(all(srt_qc$CellQC %in% c("Pass", "Fail")))
})

# Test RunDEtest basic functionality
test_that("RunDEtest performs differential expression", {
  skip_on_cran()

  srt <- create_test_data(ncells = 100, ngenes = 200)
  srt <- NormalizeData(srt, verbose = FALSE)

  srt_de <- RunDEtest(
    srt = srt,
    group_by = "CellType",
    fc.threshold = 1,
    only.pos = FALSE
  )

  # Check DE results exist
  expect_true("DEtest_CellType" %in% names(srt_de@tools))
  expect_true("AllMarkers_wilcox" %in% names(srt_de@tools$DEtest_CellType))
})

# Test visualization functions
test_that("Basic visualization functions work", {
  skip_on_cran()

  srt <- create_test_data(ncells = 100, ngenes = 200)
  srt <- NormalizeData(srt, verbose = FALSE)
  srt <- FindVariableFeatures(srt, verbose = FALSE)
  srt <- ScaleData(srt, verbose = FALSE)
  srt <- RunPCA(srt, verbose = FALSE)
  srt <- RunUMAP(srt, dims = 1:10, verbose = FALSE)

  # Test CellDimPlot
  p1 <- CellDimPlot(
    srt = srt,
    group.by = "CellType",
    reduction = "umap"
  )
  expect_true(inherits(p1, "ggplot") || inherits(p1, "patchwork"))

  # Test CellStatPlot
  p2 <- CellStatPlot(
    srt = srt,
    stat.by = "CellType",
    plot_type = "bar"
  )
  expect_true(inherits(p2, "ggplot"))
})

# Test data type checking
test_that("check_DataType identifies data correctly", {
  srt <- create_test_data(ncells = 50, ngenes = 100)

  # Check raw counts
  status <- check_DataType(srt, layer = "counts", assay = "RNA")
  expect_equal(status, "raw_counts")

  # Check after normalization
  srt <- NormalizeData(srt, verbose = FALSE)
  status <- check_DataType(srt, layer = "data", assay = "RNA")
  expect_true(status %in% c("log_normalized_counts", "raw_normalized_counts"))
})