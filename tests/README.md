# SCP Package Tests

This directory contains comprehensive tests for the SCP (Single Cell Pipeline) package, including tests for all major workflow functions and compatibility with both Seurat V4 and V5.

## Test Structure

- `testthat/` - Contains all unit tests
  - `test-workflow-comprehensive.R` - Comprehensive tests for all major SCP functions
  - Other test files for specific functionality

- `testthat.R` - Standard R package test runner

## Running Tests

### Quick Tests
For a quick check of basic functionality:
```r
source("run_quick_tests.R")
```

### Comprehensive Tests
To run all comprehensive workflow tests:
```r
source("run_comprehensive_tests.R")
```

### Using testthat
To run all tests using the testthat framework:
```r
library(testthat)
library(SCP)
test_package("SCP")
```

### Running Specific Tests
To run only the comprehensive workflow tests:
```r
library(testthat)
test_file("tests/testthat/test-workflow-comprehensive.R")
```

## Test Coverage

The comprehensive test suite covers:

1. **Standard_SCP Pipeline**
   - Normalization
   - HVF finding
   - Scaling
   - Linear dimension reduction (PCA)
   - Non-linear dimension reduction (UMAP, t-SNE, etc.)
   - Clustering

2. **Integration_SCP**
   - Multiple integration methods (Uncorrected, Harmony, Seurat, etc.)
   - Batch effect correction
   - Both srtList and srtMerge inputs

3. **Seurat Compatibility**
   - Automatic detection of Seurat V4 vs V5
   - Proper data access for both versions
   - Version conversion utilities

4. **Analysis Functions**
   - RunDEtest (differential expression)
   - RunEnrichment (pathway enrichment)
   - RunCellQC (quality control)

5. **Visualization**
   - CellDimPlot
   - FeatureDimPlot
   - CellStatPlot
   - GroupHeatmap

6. **Cell Annotation**
   - RunSingleR
   - Cell type prediction

7. **Advanced Features** (if Python environment available)
   - RunPAGA
   - Imputation methods

## Requirements

- R >= 4.1.0
- Seurat (V4 or V5)
- SCP package and its dependencies
- Optional: Python environment for advanced features

## Notes

- Some tests are skipped on CRAN to reduce testing time
- Tests requiring Python packages are skipped if the environment is not available
- Integration method tests are skipped if required packages are not installed
- The test data is generated randomly for reproducibility