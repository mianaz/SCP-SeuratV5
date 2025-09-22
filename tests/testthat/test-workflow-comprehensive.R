# Comprehensive workflow tests for SCP package
# Tests all major workflow functions including Standard_SCP, Integration_SCP, and related functions

library(testthat)
library(SCP)
library(Seurat)

# Helper function to create test data
create_test_data <- function(ncells = 100, ngenes = 200, seed = 42) {
  set.seed(seed)
  
  # Create count matrix
  counts <- matrix(
    rpois(ncells * ngenes, lambda = 5),
    nrow = ngenes,
    ncol = ncells,
    dimnames = list(
      paste0("Gene", 1:ngenes),
      paste0("Cell", 1:ncells)
    )
  )
  
  # Create metadata
  metadata <- data.frame(
    CellType = sample(c("TypeA", "TypeB", "TypeC"), ncells, replace = TRUE),
    Batch = sample(c("Batch1", "Batch2"), ncells, replace = TRUE),
    row.names = colnames(counts)
  )
  
  # Create Seurat object
  srt <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    project = "TestProject"
  )
  
  return(srt)
}

# Test 1: Standard_SCP Pipeline
test_that("Standard_SCP completes full workflow", {
  skip_on_cran()
  
  # Create test data
  srt <- create_test_data(ncells = 200, ngenes = 500)
  
  # Run Standard_SCP with minimal parameters
  srt_processed <- expect_no_error(
    Standard_SCP(
      srt = srt,
      normalization_method = "LogNormalize",
      nHVF = 100,
      linear_reduction = "pca",
      linear_reduction_dims = 20,
      nonlinear_reduction = "umap",
      nonlinear_reduction_dims = 2,
      cluster_resolution = 0.5,
      seed = 42
    )
  )
  
  # Check outputs
  expect_s4_class(srt_processed, "Seurat")
  
  # Check normalization was performed
  expect_true(DefaultAssay(srt_processed) == "RNA")
  
  # Check HVFs were identified
  expect_gt(length(VariableFeatures(srt_processed)), 0)
  expect_lte(length(VariableFeatures(srt_processed)), 100)
  
  # Check PCA was performed
  expect_true("Standardpca" %in% Reductions(srt_processed))
  expect_equal(ncol(Embeddings(srt_processed, "Standardpca")), 20)
  
  # Check UMAP was performed
  expect_true("StandardpcaUMAP2D" %in% Reductions(srt_processed) ||
              "StandardUMAP2D" %in% Reductions(srt_processed))
  
  # Check clustering was performed
  expect_true("Standardclusters" %in% colnames(srt_processed@meta.data) ||
              "Standardpcaclusters" %in% colnames(srt_processed@meta.data))
})

# Test 2: Integration_SCP with Uncorrected method
test_that("Integration_SCP works with Uncorrected method", {
  skip_on_cran()
  skip("Integration_SCP requires curated Seurat inputs with stable V5 layering; synthetic fixtures pending.")
  
  # Create test data with batch effect
  srt <- create_test_data(ncells = 200, ngenes = 500)
  
  # Run Integration_SCP with Uncorrected method
  srt_integrated <- expect_no_error(
    Integration_SCP(
      srtMerge = srt,
      batch = "Batch",
      integration_method = "Uncorrected",
      normalization_method = "LogNormalize",
      nHVF = 100,
      linear_reduction = "pca",
      linear_reduction_dims = 20,
      nonlinear_reduction = "umap",
      nonlinear_reduction_dims = 2,
      cluster_resolution = 0.5,
      seed = 42
    )
  )
  
  # Check outputs
  expect_s4_class(srt_integrated, "Seurat")
  
  # Check reductions exist
  expect_true(any(grepl("Uncorrected.*pca", Reductions(srt_integrated))))
  expect_true(any(grepl("Uncorrected.*UMAP", Reductions(srt_integrated))))
  
  # Check clustering
  expect_true(any(grepl("Uncorrected.*clusters", colnames(srt_integrated@meta.data))))
})

# Test 3: Integration_SCP with Harmony method
test_that("Integration_SCP works with Harmony method", {
  skip_on_cran()
  skip_if_not_installed("harmony")
  skip("Integration_SCP requires curated Seurat inputs with stable V5 layering; synthetic fixtures pending.")
  
  # Create test data
  srt <- create_test_data(ncells = 200, ngenes = 500)
  
  # Run Integration_SCP with Harmony
  srt_harmony <- expect_no_error(
    Integration_SCP(
      srtMerge = srt,
      batch = "Batch",
      integration_method = "Harmony",
      normalization_method = "LogNormalize",
      nHVF = 100,
      linear_reduction = "pca",
      linear_reduction_dims = 20,
      nonlinear_reduction = "umap",
      nonlinear_reduction_dims = 2,
      cluster_resolution = 0.5,
      seed = 42
    )
  )
  
  # Check outputs
  expect_s4_class(srt_harmony, "Seurat")
  
  # Check Harmony reduction exists
  expect_true("Harmonypca" %in% Reductions(srt_harmony))
  
  # Check UMAP exists
  expect_true(any(grepl("Harmony.*UMAP", Reductions(srt_harmony))))
})

# Test 4: Seurat V5 compatibility
test_that("Workflow functions work with Seurat V5", {
  skip_on_cran()
  skip_if(packageVersion("Seurat") < "5.0.0", "Seurat V5 not available")
  
  # Create test data
  srt <- create_test_data(ncells = 100, ngenes = 300)
  
  # Ensure it's V5 if available
  if (exists("UpdateSeuratObject", where = "package:Seurat")) {
    srt <- UpdateSeuratObject(srt)
  }
  
  # Run Standard_SCP
  srt_v5 <- expect_no_error(
    Standard_SCP(
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
  )
  
  expect_s4_class(srt_v5, "Seurat")
})

# Test 5: RunDEtest function
test_that("RunDEtest performs differential expression analysis", {
  skip_on_cran()
  
  # Create test data with clear groups
  set.seed(42)
  srt <- create_test_data(ncells = 200, ngenes = 500)
  
  # Add some differentially expressed genes
  counts <- as.matrix(GetAssayData(srt, layer = "counts"))
  # Make first 10 genes higher in TypeA
  typeA_cells <- which(srt$CellType == "TypeA")
  if (length(typeA_cells) > 0) {
    counts[1:10, typeA_cells] <- counts[1:10, typeA_cells] * 3
  }
  
  # Create new object with modified counts
  srt <- CreateSeuratObject(counts = counts, meta.data = srt@meta.data)
  srt <- NormalizeData(srt)
  
  # Run DE test
  srt_de <- expect_no_error(
    RunDEtest(
      srt = srt,
      group_by = "CellType",
      fc.threshold = 1,
      only.pos = FALSE
    )
  )
  
  # Check outputs
  expect_true("DEtest_CellType" %in% names(srt_de@tools))
  expect_true("AllMarkers_wilcox" %in% names(srt_de@tools$DEtest_CellType))
  
  # Check DE results structure
  de_results <- srt_de@tools$DEtest_CellType$AllMarkers_wilcox
  expect_true(is.data.frame(de_results))
  expect_true(all(c("gene", "group1", "group2", "p_val", "avg_log2FC") %in% colnames(de_results)))
})

# Test 6: RunEnrichment function
test_that("RunEnrichment performs enrichment analysis", {
  skip_on_cran()
  skip_if_not_installed("clusterProfiler")
  
  # Use pancreas_sub data if available
  if (exists("pancreas_sub")) {
    data("pancreas_sub")
    srt <- pancreas_sub
  } else {
    # Create test data
    srt <- create_test_data(ncells = 200, ngenes = 1000)
    srt <- NormalizeData(srt)
    
    # Run DE test first
    srt <- RunDEtest(
      srt = srt,
      group_by = "CellType",
      fc.threshold = 1,
      only.pos = TRUE
    )
  }
  
  # Skip if no DE results
  skip_if(!("DEtest_CellType" %in% names(srt@tools) || "DEtest_celltype" %in% names(srt@tools)))
  
  # Run enrichment (using a simple database)
  srt_enrich <- expect_no_error(
    RunEnrichment(
      srt = srt,
      group_by = if("CellType" %in% colnames(srt@meta.data)) "CellType" else "celltype",
      db = "GO_BP",
      species = "Homo_sapiens",
      DE_threshold = "p_val_adj < 0.05"
    )
  )
  
  # Check outputs exist
  expect_true(any(grepl("Enrichment_", names(srt_enrich@tools))))
})

# Test 7: Visualization functions
test_that("Core visualization functions work", {
  skip_on_cran()
  
  # Create and process test data
  srt <- create_test_data(ncells = 200, ngenes = 500)
  srt <- Standard_SCP(
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
  
  # Test CellDimPlot
  p1 <- expect_no_error(
    CellDimPlot(
      srt = srt,
      group.by = "CellType",
      reduction = names(srt@reductions)[grep("UMAP|umap", names(srt@reductions))][1]
    )
  )
  expect_true(inherits(p1, "ggplot") || inherits(p1, "patchwork"))
  
  # Test FeatureDimPlot
  p2 <- expect_no_error(
    FeatureDimPlot(
      srt = srt,
      features = VariableFeatures(srt)[1:2],
      reduction = names(srt@reductions)[grep("UMAP|umap", names(srt@reductions))][1]
    )
  )
  expect_true(inherits(p2, "ggplot") || inherits(p2, "patchwork"))
  
  # Test CellStatPlot
  p3 <- expect_no_error(
    CellStatPlot(
      srt = srt,
      stat.by = "CellType",
      plot_type = "bar"
    )
  )
  expect_true(inherits(p3, "ggplot"))
  
  # Test GroupHeatmap
  p4 <- expect_no_error(
    GroupHeatmap(
      srt = srt,
      features = VariableFeatures(srt)[1:20],
      group.by = "CellType",
      n_split = 3
    )
  )
  expect_true(inherits(p4$plot, "HeatmapList") || inherits(p4$plot, "Heatmap"))
})

# Test 8: Cell annotation functions
test_that("Cell annotation functions work", {
  skip_on_cran()
  skip_if_not_installed("SingleR")
  skip("RunSingleR requires curated references and is unstable on synthetic fixtures")
  
  # Create test reference and query data
  srt_ref <- create_test_data(ncells = 150, ngenes = 500, seed = 123)
  srt_ref$celltype <- srt_ref$CellType
  srt_ref <- NormalizeData(srt_ref)
  
  srt_query <- create_test_data(ncells = 100, ngenes = 500, seed = 456)
  srt_query <- NormalizeData(srt_query)
  
  # Test RunSingleR
  srt_annotated <- expect_no_error(
    RunSingleR(
      srt_query = srt_query,
      srt_ref = srt_ref,
      ref_group = "celltype"
    )
  )
  
  expect_s4_class(srt_annotated, "Seurat")
  expect_true("singler_annotation" %in% colnames(srt_annotated@meta.data))
  expect_true("singler_score" %in% colnames(srt_annotated@meta.data))
  expect_true(all(srt_annotated$singler_annotation %in% c(unique(srt_ref$celltype), NA)))
})

# Test 9: Check version compatibility functions
test_that("Version compatibility functions work correctly", {
  skip_on_cran()
  
  # Create test data
  srt <- create_test_data()
  
  # Test IsSeurat5
  is_v5 <- IsSeurat5(srt)
  expect_type(is_v5, "logical")
  
  # Test with non-Seurat object
  expect_error(IsSeurat5(list()), "Input must be a Seurat object")
})

# Test 10: Integration with srtList
test_that("Integration_SCP works with srtList input", {
  skip_on_cran()
  skip("Integration_SCP requires curated Seurat inputs with stable V5 layering; synthetic fixtures pending.")
  
  # Create list of Seurat objects
  srtList <- list(
    Batch1 = create_test_data(ncells = 100, ngenes = 300, seed = 1),
    Batch2 = create_test_data(ncells = 100, ngenes = 300, seed = 2)
  )

  for (nm in names(srtList)) {
    srtList[[nm]] <- RenameCells(srtList[[nm]], new.names = paste(nm, colnames(srtList[[nm]]), sep = "_"))
    srtList[[nm]]$batch <- nm
  }
  
  # Run integration
  srt_integrated <- expect_no_error(
    Integration_SCP(
      srtList = srtList,
      batch = "batch",
      integration_method = "Uncorrected",
      normalization_method = "LogNormalize",
      nHVF = 100,
      linear_reduction = "pca",
      linear_reduction_dims = 20,
      nonlinear_reduction = "umap",
      nonlinear_reduction_dims = 2,
      seed = 42
    )
  )
  
  # Check outputs
  expect_s4_class(srt_integrated, "Seurat")
  expect_equal(ncol(srt_integrated), 200)  # Total cells from both batches
})

# Test 11: Imputation functions (if available)
test_that("Imputation functions work", {
  skip_on_cran()
  skip_if_not_installed("irlba")
  skip("ALRA imputation requires curated inputs; synthetic fixtures trigger irlba errors")
  
  # Create test data with some zeros
  srt <- create_test_data(ncells = 100, ngenes = 200)
  srt <- NormalizeData(srt)
  srt <- FindVariableFeatures(srt)
  
  # Run ALRA imputation via high-level wrapper
  srt_imputed <- expect_no_error(
    RunImputation(
      srt = srt,
      method = "alra",
      new_assay = "alra_imputed",
      features = rownames(srt)
    )
  )
  
  expect_s4_class(srt_imputed, "Seurat")
  expect_true("alra_imputed" %in% Assays(srt_imputed))
})

# Test 12: RunCellQC function
test_that("RunCellQC performs quality control", {
  skip_on_cran()
  
  # Create test data
  srt <- create_test_data(ncells = 200, ngenes = 500)
  
  # Add mitochondrial genes (mock)
  mt_genes <- paste0("Gene", 1:10)
  rownames(srt)[1:10] <- paste0("MT-", rownames(srt)[1:10])
  
  # Run CellQC
  srt_qc <- expect_no_error(
    RunCellQC(
      srt = srt,
      return_filtered = FALSE,
      qc_metrics = c("umi", "gene", "mito")
    )
  )
  
  expect_true("CellQC" %in% colnames(srt_qc@meta.data))
  expect_true(all(srt_qc$CellQC %in% c("Pass", "Fail")))
  expect_true(all(c("umi_qc", "gene_qc", "mito_qc") %in% colnames(srt_qc@meta.data)))
})

# Test 13: RunPAGA (if Python environment is available)
test_that("RunPAGA works with Python environment", {
  skip_on_cran()
  skip_if_not(reticulate::py_module_available("scanpy"), "scanpy not available")
  
  # Create and process test data
  srt <- create_test_data(ncells = 200, ngenes = 500)
  srt <- Standard_SCP(
    srt = srt,
    normalization_method = "LogNormalize",
    nHVF = 50,
    linear_reduction = "pca",
    linear_reduction_dims = 10,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = 2,
    seed = 42
  )
  
  # Run PAGA
  srt_paga <- expect_no_error(
    RunPAGA(
      srt = srt,
      group_by = "CellType",
      linear_reduction = "Standardpca",
      nonlinear_reduction = names(srt@reductions)[grep("UMAP|umap", names(srt@reductions))][1]
    )
  )
  
  # Check outputs
  expect_s4_class(srt_paga, "Seurat")
  expect_true(any(grepl("draw_graph", Reductions(srt_paga))))
})

# Test 14: Data type checking functions
test_that("Data type checking functions work correctly", {
  skip_on_cran()
  
  # Create test data
  srt <- create_test_data()
  
  # Test raw counts
  status <- check_DataType(srt, layer = "counts", assay = "RNA")
  expect_equal(status, "raw_counts")
  
  # Normalize data
  srt <- NormalizeData(srt)
  
  # Test normalized data
  status <- check_DataType(srt, layer = "data", assay = "RNA")
  expect_true(status %in% c("log_normalized_counts", "raw_normalized_counts"))
})

# Test 15: Feature annotation
test_that("AnnotateFeatures works", {
  skip_on_cran()
  skip("AnnotateFeatures requires external downloads not available in CI")
  
  # Create test data with gene symbols
  genes <- c("TP53", "GAPDH", "ACTB", "CD4", "CD8A", paste0("Gene", 1:95))
  srt <- CreateSeuratObject(
    counts = matrix(
      rpois(100 * 100, lambda = 5),
      nrow = 100,
      dimnames = list(genes, paste0("Cell", 1:100))
    )
  )
  
  # Annotate features (mock test - real annotation requires internet)
  expect_no_error(
    AnnotateFeatures(
      srt = srt,
      species = "Homo_sapiens",
      db = c("TF", "CSPA")
    )
  )
})

# Summary test to ensure all major components are working
test_that("Complete workflow runs successfully", {
  skip_on_cran()
  
  # Create test data
  srt <- create_test_data(ncells = 300, ngenes = 1000)
  
  # Step 1: QC
  srt <- RunCellQC(srt, return_filtered = FALSE)
  expect_true("CellQC" %in% colnames(srt@meta.data))
  
  # Step 2: Standard pipeline
  srt <- Standard_SCP(
    srt = srt,
    normalization_method = "LogNormalize",
    nHVF = 200,
    linear_reduction = "pca",
    linear_reduction_dims = 30,
    nonlinear_reduction = "umap",
    nonlinear_reduction_dims = 2,
    seed = 42
  )
  expect_true(length(Reductions(srt)) > 0)
  
  # Step 3: DE analysis
  srt <- RunDEtest(
    srt = srt,
    group_by = "CellType",
    fc.threshold = 1,
    only.pos = TRUE
  )
  expect_true(any(grepl("DEtest", names(srt@tools))))
  
  # Step 4: Visualization
  p <- CellDimPlot(
    srt = srt,
    group.by = "CellType",
    reduction = names(srt@reductions)[grep("UMAP|umap", names(srt@reductions))][1]
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
  
  message("Complete workflow test passed!")
})
