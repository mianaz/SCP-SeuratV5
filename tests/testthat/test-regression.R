# Regression tests for known bug fixes
# Tests to prevent re-introduction of previously fixed bugs
library(testthat)
library(SCPNext)
library(Seurat)

# ========== Bug Fix #2: HVF variable not updated after FindVariableFeatures ==========
# Fixed in R/SCP-workflow.R - Integration_SCP now updates HVF variable
# Previously: HVF stayed empty, ScaleData scaled 0 features, PCA failed

test_that("Integration updates HVF after FindVariableFeatures (Bug #2)", {
  skip_on_cran()

  data("panc8_sub")

  result <- Integration_SCP(
    srtMerge = panc8_sub,
    batch = "tech",
    integration_method = "Seurat",
    use_v5_workflow = TRUE,
    seed = 11
  )

  # Verify HVF variable is populated
  hvf <- VariableFeatures(result)
  expect_gt(
    length(hvf),
    1000,
    info = "HVF should contain ~2000 features after FindVariableFeatures"
  )

  # Verify ScaleData used correct features
  expect_true(
    "SeuratUMAP2D" %in% Reductions(result),
    info = "UMAP should be created (depends on valid PCA from scaled features)"
  )
})

# ========== Bug Fix #11: BiocParallel TERM2NAME subscript error ==========
# Fixed in R/SCP-analysis.R - Added tryCatch for subsetting in parallel workers
# Previously: Parallel workers couldn't access TERM2NAME causing subscript errors

test_that("Enrichment handles BiocParallel safely (Bug #11)", {
  skip_on_cran()

  data("pancreas_sub")

  # Run DE first
  pancreas_de <- RunDEtest(
    srt = pancreas_sub,
    group_by = "CellType",
    fc.threshold = 1,
    only.pos = FALSE
  )

  # Run enrichment - should not error with BiocParallel
  result <- expect_no_error(
    RunEnrichment(
      srt = pancreas_de,
      group_by = "CellType",
      db = "GO_BP",
      species = "Mus_musculus",
      DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
    ),
    info = "RunEnrichment should handle BiocParallel without subscript errors"
  )

  expect_s4_class(result, "Seurat")
})

# ========== Bug Fix #12: packageVersion("SCP") in RunSCExplorer ==========
# Fixed in R/SCP-app.R - Changed to packageVersion("SCPNext")
# Previously: Error "there is no package called 'SCP'"

test_that("SCExplorer uses correct package name (Bug #12)", {
  skip_on_cran()

  # This test verifies the fix is in place by checking the function code
  # We don't actually run SCExplorer (requires Shiny)

  code <- deparse(body(PrepareSCExplorer))
  code_text <- paste(code, collapse = "\n")

  # Should reference SCPNext, not SCP
  expect_true(
    grepl("SCPNext", code_text, fixed = TRUE),
    info = "PrepareSCExplorer should reference 'SCPNext' not 'SCP'"
  )

  expect_false(
    grepl('packageVersion\\("SCP"\\)', code_text),
    info = "Should not call packageVersion('SCP')"
  )
})

# ========== Seurat V5 Compatibility ==========
# Verify package works with both Seurat v4 and v5

test_that("Package maintains backward compatibility with Seurat v4/v5", {
  skip_on_cran()

  # Create simple Seurat object
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell", 1:10)
  srt <- CreateSeuratObject(counts = counts)

  # IsSeurat5 should work regardless of version
  result <- expect_no_error(
    IsSeurat5(srt),
    info = "IsSeurat5 should work with both v4 and v5"
  )

  expect_type(result, "logical")
  expect_length(result, 1)

  # Basic workflow should work
  srt_processed <- expect_no_error(
    {
      srt <- NormalizeData(srt, verbose = FALSE)
      srt <- FindVariableFeatures(srt, verbose = FALSE)
      srt <- ScaleData(srt, verbose = FALSE)
      srt
    },
    info = "Basic Seurat workflow should work"
  )

  expect_s4_class(srt_processed, "Seurat")
})

# ========== Data Type Detection ==========
# Verify check_DataType correctly identifies different data states

test_that("check_DataType correctly identifies data types", {
  skip_on_cran()

  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell", 1:10)
  srt <- CreateSeuratObject(counts = counts)

  # Test raw counts
  status <- check_DataType(srt, layer = "counts", assay = "RNA")
  expect_equal(status, "raw_counts")

  # Test normalized data
  srt <- NormalizeData(srt, verbose = FALSE)
  status <- check_DataType(srt, layer = "data", assay = "RNA")
  expect_true(status %in% c("log_normalized_counts", "raw_normalized_counts"))
})

# ========== UV Environment Basics ==========
# Verify UV environment functions work

test_that("UV environment checks work correctly", {
  skip_on_cran()

  # check_uv should not error
  result <- expect_no_error(
    check_uv(),
    info = "check_uv should run without errors"
  )

  expect_type(result, "logical")

  # uv_env_exists should not error
  result2 <- expect_no_error(
    uv_env_exists(),
    info = "uv_env_exists should run without errors"
  )

  expect_type(result2, "logical")
})

# ========== Palette Functions ==========
# Verify consolidated palette functions work

test_that("Consolidated palette functions work", {
  skip_on_cran()

  # Test GetPalette (consolidated function)
  colors <- GetPalette(n = 5, palette = "Paired")
  expect_type(colors, "character")
  expect_length(colors, 5)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))

  # Test with different palettes
  colors2 <- GetPalette(n = 3, palette = "Set1")
  expect_length(colors2, 3)

  # Test automatic palette selection (n > palette length)
  colors3 <- GetPalette(n = 50)
  expect_length(colors3, 50)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors3)))
})

# ========== P-value Combination ==========
# Verify consolidated p-value combination function

test_that("Consolidated CombinePvalues function works", {
  skip_on_cran()

  # Test Fisher's method
  pvals <- c(0.01, 0.05, 0.10, 0.20)
  result_fisher <- expect_no_error(
    CombinePvalues(pvals, method = "fisher"),
    info = "Fisher's method should work"
  )
  expect_type(result_fisher, "double")
  expect_true(result_fisher >= 0 && result_fisher <= 1)

  # Test Stouffer's method
  result_stouffer <- expect_no_error(
    CombinePvalues(pvals, method = "stouffer"),
    info = "Stouffer's method should work"
  )
  expect_type(result_stouffer, "double")
  expect_true(result_stouffer >= 0 && result_stouffer <= 1)

  # Test edge cases
  result_all_sig <- CombinePvalues(rep(0.001, 5), method = "fisher")
  expect_lt(result_all_sig, 0.05, info = "All significant p-values should combine to significant")

  result_all_nonsig <- CombinePvalues(rep(0.9, 5), method = "fisher")
  expect_gt(result_all_nonsig, 0.05, info = "All non-significant p-values should combine to non-significant")
})
