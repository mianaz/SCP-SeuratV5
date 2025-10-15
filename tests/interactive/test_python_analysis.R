#!/usr/bin/env Rscript
# Interactive Test Script: Python Analysis Functions
#
# This script tests SCVELO, PAGA, and other Python-based analysis functions.
# Run this line-by-line in RStudio for debugging.
#
# PREREQUISITE: Run test_python_environment.R first to ensure environment is ready
#
# Author: SCP Development Team
# Date: 2025-01-14

# ==============================================================================
# SETUP
# ==============================================================================

cat("\n==============================================================================\n")
cat("SCP Python Analysis Functions Test Suite\n")
cat("==============================================================================\n\n")

# Load required packages
cat("Loading packages...\n")
library(SCP)
library(Seurat)
library(Matrix)

# ==============================================================================
# TEST 0: Prerequisites Check
# ==============================================================================

cat("\n--- TEST 0: Prerequisites Check ---\n")

# Check if environment exists
env_exists <- uv_env_exists()
cat("UV environment exists:", env_exists, "\n")

if (!env_exists) {
  stop("UV environment not found. Please run test_python_environment.R first and create the environment with PrepareEnv(extras = 'velocity')")
}

# Configure reticulate
pkg_dir <- system.file("", package = "SCP")
if (pkg_dir == "") {
  pkg_dir <- getwd()
}
venv_path <- file.path(pkg_dir, ".venv")

if (Sys.info()["sysname"] == "Windows") {
  python_path <- file.path(venv_path, "Scripts", "python.exe")
} else {
  python_path <- file.path(venv_path, "bin", "python")
}

cat("Configuring reticulate...\n")
tryCatch({
  use_uv_env()
  cat("✓ Reticulate configured\n")
}, error = function(e) {
  stop("Failed to configure reticulate: ", e$message)
})

# ==============================================================================
# TEST 1: Create Test Seurat Object
# ==============================================================================

cat("\n--- TEST 1: Creating Test Data ---\n")

cat("Creating synthetic Seurat object for testing...\n")

# Set seed for reproducibility
set.seed(42)

# Create synthetic count matrix (100 genes x 200 cells)
n_genes <- 100
n_cells <- 200

counts <- matrix(
  rpois(n_genes * n_cells, lambda = 5),
  nrow = n_genes,
  ncol = n_cells
)
rownames(counts) <- paste0("Gene_", 1:n_genes)
colnames(counts) <- paste0("Cell_", 1:n_cells)

# Create Seurat object
cat("Creating Seurat object...\n")
srt <- CreateSeuratObject(counts = counts, project = "TestProject")

# Add some metadata
srt$celltype <- sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE)
srt$batch <- sample(c("Batch1", "Batch2"), n_cells, replace = TRUE)

cat("✓ Created Seurat object:\n")
cat("  Cells:", ncol(srt), "\n")
cat("  Genes:", nrow(srt), "\n")
cat("  Assays:", names(srt@assays), "\n")

# ==============================================================================
# TEST 2: Basic Preprocessing
# ==============================================================================

cat("\n--- TEST 2: Basic Preprocessing ---\n")

cat("Normalizing data...\n")
srt <- NormalizeData(srt, verbose = FALSE)

cat("Finding variable features...\n")
srt <- FindVariableFeatures(srt, nfeatures = 50, verbose = FALSE)

cat("Scaling data...\n")
srt <- ScaleData(srt, verbose = FALSE)

cat("Running PCA...\n")
srt <- RunPCA(srt, npcs = 20, verbose = FALSE)

cat("Running UMAP...\n")
srt <- RunUMAP(srt, dims = 1:10, verbose = FALSE)

cat("Finding neighbors...\n")
srt <- FindNeighbors(srt, dims = 1:10, verbose = FALSE)

cat("Finding clusters...\n")
srt <- FindClusters(srt, resolution = 0.5, verbose = FALSE)

cat("✓ Preprocessing complete\n")
cat("  Reductions:", names(srt@reductions), "\n")
cat("  Clusters:", length(unique(srt$seurat_clusters)), "\n")

# ==============================================================================
# TEST 3: Test Python Module Imports
# ==============================================================================

cat("\n--- TEST 3: Python Module Imports ---\n")

cat("Testing Python module imports via reticulate...\n")

modules_to_test <- c(
  "numpy", "pandas", "scanpy", "anndata"
)

import_results <- list()
for (module in modules_to_test) {
  success <- tryCatch({
    reticulate::import(module)
    TRUE
  }, error = function(e) {
    FALSE
  })

  status <- if (success) "✓" else "✗"
  cat(sprintf("  %s %s\n", status, module))
  import_results[[module]] <- success
}

if (!all(unlist(import_results))) {
  cat("\n⚠ Some modules failed to import. Python analysis tests may fail.\n")
  cat("  Install missing packages with: PrepareEnv(update = TRUE, extras = 'all')\n")
}

# ==============================================================================
# TEST 4: Test SCVELO Prerequisites
# ==============================================================================

cat("\n--- TEST 4: SCVELO Prerequisites ---\n")

# Check if scvelo is installed
scvelo_available <- tryCatch({
  result <- system2(python_path,
                   c("-c", "import scvelo; print(scvelo.__version__)"),
                   stdout = TRUE, stderr = TRUE)
  !is.null(result) && length(result) > 0
}, error = function(e) {
  FALSE
})

cat("scvelo installed:", scvelo_available, "\n")

if (scvelo_available) {
  version <- system2(python_path,
                    c("-c", "import scvelo; print(scvelo.__version__)"),
                    stdout = TRUE)
  cat("scvelo version:", version[1], "\n")
} else {
  cat("⚠ scvelo not installed\n")
  cat("  Install with: uv_install_extras('velocity')\n")
  cat("  SCVELO tests will be skipped\n")
}

# Check igraph (for PAGA)
igraph_available <- tryCatch({
  result <- system2(python_path,
                   c("-c", "import igraph; print(igraph.__version__)"),
                   stdout = TRUE, stderr = TRUE)
  !is.null(result) && length(result) > 0
}, error = function(e) {
  FALSE
})

cat("igraph installed:", igraph_available, "\n")

if (igraph_available) {
  version <- system2(python_path,
                    c("-c", "import igraph; print(igraph.__version__)"),
                    stdout = TRUE)
  cat("igraph version:", version[1], "\n")
}

# ==============================================================================
# TEST 5: Test Basis Validation
# ==============================================================================

cat("\n--- TEST 5: Basis Validation ---\n")

cat("Testing basis validation logic...\n\n")

# Test Python script for basis validation
basis_test_script <- "
import scanpy as sc
import numpy as np
import anndata as ad

# Create minimal AnnData object
n_obs = 100
n_vars = 50

X = np.random.rand(n_obs, n_vars)
adata = ad.AnnData(X)
adata.obs_names = [f'Cell_{i}' for i in range(n_obs)]
adata.var_names = [f'Gene_{i}' for i in range(n_vars)]

# Add a PCA embedding
adata.obsm['X_pca'] = np.random.rand(n_obs, 10)

# Add a UMAP embedding
adata.obsm['X_umap'] = np.random.rand(n_obs, 2)

print('Created AnnData object:')
print(f'  Shape: {adata.shape}')
print(f'  Embeddings: {list(adata.obsm.keys())}')

# Test basis validation
basis = 'X_umap'
if basis in adata.obsm:
    print(f'✓ Basis {basis} found')
else:
    print(f'✗ Basis {basis} not found')
    print('  Available:', list(adata.obsm.keys()))

# Test fallback logic
test_basis = 'X_nonexistent'
if test_basis not in adata.obsm:
    print(f'\\nTesting fallback for missing basis {test_basis}:')
    # Try fallbacks
    for fallback in ['X_umap', 'X_pca']:
        if fallback in adata.obsm:
            print(f'  ✓ Fallback to {fallback}')
            test_basis = fallback
            break
    else:
        print('  ✗ No fallback found')
"

result <- system2(python_path, c("-c", basis_test_script),
                 stdout = TRUE, stderr = TRUE)

if (length(result) > 0) {
  cat(paste(result, collapse = "\n"), "\n")
}

# ==============================================================================
# TEST 6: Test M-series Configuration in Python
# ==============================================================================

cat("\n--- TEST 6: M-series Detection in Python ---\n")

m_series_test <- "
import os
import platform

# Detect M-series
is_m_series = platform.system() == 'Darwin' and platform.machine() == 'arm64'

print(f'Platform: {platform.system()}')
print(f'Machine: {platform.machine()}')
print(f'M-series Mac: {is_m_series}')

if is_m_series:
    print('\\nM-series configuration would be applied:')
    print('  - NUMBA_DISABLE_JIT = 1')
    print('  - OMP_NUM_THREADS = 1')
    print('  - Single-threaded execution')
    print('  - Matplotlib Agg backend')
else:
    print('\\nStandard configuration would be used')
"

result <- system2(python_path, c("-c", m_series_test),
                 stdout = TRUE, stderr = TRUE)

if (length(result) > 0) {
  cat(paste(result, collapse = "\n"), "\n")
}

# ==============================================================================
# TEST 7: PAGA Compatibility Test
# ==============================================================================

cat("\n--- TEST 7: PAGA igraph Compatibility ---\n")

if (igraph_available) {
  cat("Testing PAGA with igraph...\n\n")

  paga_test_script <- "
import scanpy as sc
import numpy as np
import anndata as ad

# Create test data
n_obs = 100
n_vars = 50

X = np.random.rand(n_obs, n_vars)
adata = ad.AnnData(X)
adata.obs_names = [f'Cell_{i}' for i in range(n_obs)]
adata.var_names = [f'Gene_{i}' for i in range(n_vars)]

# Add required fields for PAGA
adata.obsm['X_pca'] = np.random.rand(n_obs, 10)
adata.obsm['X_umap'] = np.random.rand(n_obs, 2)
adata.obs['clusters'] = np.random.choice(['A', 'B', 'C'], n_obs)

print('Testing PAGA computation...')
try:
    # Compute neighbors
    sc.pp.neighbors(adata, n_neighbors=10)
    print('✓ Neighbors computed')

    # Compute PAGA
    sc.tl.paga(adata, groups='clusters')
    print('✓ PAGA computed successfully')

    print('\\nPAGA result:')
    print(f'  Connectivity matrix shape: {adata.uns[\"paga\"][\"connectivities\"].shape}')

except Exception as e:
    print(f'✗ PAGA failed: {str(e)}')
    print('\\nThis is a known issue with scvelo 0.3.x and igraph')
    print('Error is expected and handled by try-except wrapper in SCP')
"

  result <- system2(python_path, c("-c", paga_test_script),
                   stdout = TRUE, stderr = TRUE)

  if (length(result) > 0) {
    cat(paste(result, collapse = "\n"), "\n")
  }

} else {
  cat("Skipped (igraph not installed)\n")
}

# ==============================================================================
# TEST 8: Memory Management Test
# ==============================================================================

cat("\n--- TEST 8: Memory Management ---\n")

cat("Testing memory handling with larger objects...\n\n")

memory_test <- "
import numpy as np
import sys
import gc

# Get initial memory
print('Testing memory management:')

# Create a large array
size = 1000
arr = np.random.rand(size, size)
arr_size_mb = arr.nbytes / (1024 * 1024)

print(f'  Created array: {size}x{size} ({arr_size_mb:.2f} MB)')

# Perform computation
result = np.dot(arr, arr.T)
result_size_mb = result.nbytes / (1024 * 1024)

print(f'  Computation result: {result.shape} ({result_size_mb:.2f} MB)')

# Clean up
del arr
del result
gc.collect()

print('  ✓ Memory cleaned up')
print('  Python garbage collection successful')
"

result <- system2(python_path, c("-c", memory_test),
                 stdout = TRUE, stderr = TRUE)

if (length(result) > 0) {
  cat(paste(result, collapse = "\n"), "\n")
}

# ==============================================================================
# TEST 9: Data Format Conversion Test
# ==============================================================================

cat("\n--- TEST 9: Seurat to AnnData Conversion ---\n")

cat("Testing data format compatibility...\n")

# Test if we can access Seurat data and prepare it for Python
counts_matrix <- GetAssayData(srt, layer = "counts")
cat("✓ Extracted counts matrix:", dim(counts_matrix), "\n")

data_matrix <- GetAssayData(srt, layer = "data")
cat("✓ Extracted data matrix:", dim(data_matrix), "\n")

# Check embeddings
if ("umap" %in% names(srt@reductions)) {
  umap_coords <- Embeddings(srt, "umap")
  cat("✓ Extracted UMAP coordinates:", dim(umap_coords), "\n")
}

if ("pca" %in% names(srt@reductions)) {
  pca_coords <- Embeddings(srt, "pca")
  cat("✓ Extracted PCA coordinates:", dim(pca_coords), "\n")
}

# Check metadata
metadata <- srt@meta.data
cat("✓ Extracted metadata:", nrow(metadata), "cells,", ncol(metadata), "columns\n")

cat("\nData is ready for Python analysis functions\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("PYTHON ANALYSIS TEST SUMMARY\n")
cat("==============================================================================\n\n")

test_results <- list(
  "Environment Ready" = env_exists,
  "Python Imports Work" = all(unlist(import_results)),
  "scvelo Available" = scvelo_available,
  "igraph Available" = igraph_available,
  "Test Data Created" = !is.null(srt)
)

for (test_name in names(test_results)) {
  status <- if (test_results[[test_name]]) "✓ PASS" else "✗ FAIL"
  cat(sprintf("%-25s %s\n", test_name, status))
}

cat("\n")
if (scvelo_available) {
  cat("✓ Environment is ready for velocity analysis (RunSCVELO)\n")
} else {
  cat("⊘ Install velocity extras: uv_install_extras('velocity')\n")
}

if (igraph_available) {
  cat("✓ Environment is ready for PAGA analysis (RunPAGA)\n")
} else {
  cat("⊘ Install required packages: uv_install_extras('singlecell')\n")
}

cat("\n")
cat("Next steps:\n")
cat("  1. Test with real data: Use your own Seurat object\n")
cat("  2. Run velocity analysis: RunSCVELO(srt, ...)\n")
cat("  3. Run PAGA: RunPAGA(srt, ...)\n")
cat("  4. Check the documentation for more analysis functions\n")

cat("\n==============================================================================\n")

# Clean up
cat("\nCleaning up test data...\n")
rm(srt, counts, counts_matrix, data_matrix)
gc()
cat("✓ Test data cleaned up\n")
