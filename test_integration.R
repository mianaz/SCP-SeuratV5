# Test refactored Uncorrected_integrate with v5 layer workflow
devtools::load_all()

# Load test data
data("panc8_sub")

cat("=== Testing Uncorrected_integrate with v5 layer workflow ===\n")
cat("Input: v5 object with", ncol(panc8_sub), "cells\n\n")

# Run integration
panc8_uncorrected <- Integration_SCP(
  srtMerge = panc8_sub,
  batch = "tech",
  integration_method = "Uncorrected",
  linear_reduction = "pca",
  linear_reduction_dims = 30,
  nonlinear_reduction = "umap",
  cluster_resolution = 0.6,
  seed = 11
)

cat("\n=== Integration Complete ===\n")
cat("Reductions:", paste(names(panc8_uncorrected@reductions), collapse=", "), "\n")

# Verify UMAP was created
if ("Uncorrectedumap2D" %in% names(panc8_uncorrected@reductions)) {
  cat("✓ SUCCESS: UMAP reduction created\n")
} else {
  cat("✗ FAILED: UMAP not found\n")
  cat("Available reductions:", names(panc8_uncorrected@reductions), "\n")
}
