# Interactive Test Suite for SCP Python Integration

This directory contains interactive test scripts that can be run line-by-line in RStudio or Jupyter for debugging and verification of Python 3.10-3.13 compatibility and M-series Mac optimizations.

## Test Scripts

### 1. `test_python_environment.R`
**Purpose**: Test Python environment setup and package compatibility

**What it tests**:
- UV installation and configuration
- Python environment existence
- Python version (3.10-3.13 range)
- Core package imports (numpy, pandas, scipy, matplotlib, scanpy, anndata, igraph)
- Optional package groups (velocity, trajectory, deep learning, spatial)
- Reticulate configuration
- Python 3.10+ features (match statements, type hints)
- Environment management functions (ListEnv, VerifyEnv)

**Run this first** to ensure your environment is properly set up.

### 2. `test_m_series_mac.R`
**Purpose**: Test M-series MacBook detection and optimizations

**What it tests**:
- Platform detection (Darwin + arm64)
- M-series specific package versions
- NUMBA JIT configuration (should be disabled on M-series)
- PyTorch MPS backend availability
- Environment variables for thread control
- Threading configuration
- Memory handling

**Run this** if you're on an M-series Mac (M1, M2, M3, etc.) to verify optimizations.

### 3. `test_python_analysis.R`
**Purpose**: Test Python-based analysis functions

**What it tests**:
- Creating test Seurat objects
- Basic preprocessing pipeline
- Python module imports via reticulate
- SCVELO prerequisites
- Basis validation logic
- M-series configuration in Python code
- PAGA igraph compatibility
- Memory management
- Seurat to AnnData data format conversion

**Run this** to verify that analysis functions will work correctly.

## How to Use

### Option 1: RStudio (Recommended)

1. Open the test script in RStudio
2. Run line-by-line using `Ctrl+Enter` (Windows/Linux) or `Cmd+Enter` (Mac)
3. Read the output after each section
4. If a test fails, you can modify variables and re-run specific sections

Example workflow:
```r
# In RStudio, open: tests/interactive/test_python_environment.R
# Execute line-by-line or select a section and run

# Start from the beginning
source("tests/interactive/test_python_environment.R")

# Or run specific sections interactively
# ... (copy/paste sections as needed)
```

### Option 2: Command Line

```bash
# From the package root directory
Rscript tests/interactive/test_python_environment.R
Rscript tests/interactive/test_m_series_mac.R
Rscript tests/interactive/test_python_analysis.R
```

### Option 3: Jupyter Notebook

If you have an R kernel in Jupyter:

```r
# In a Jupyter notebook cell
source("tests/interactive/test_python_environment.R")
```

Or copy/paste sections into individual cells for step-by-step execution.

## Recommended Testing Order

1. **First time setup**:
   ```r
   # Step 1: Install UV (if not already installed)
   install_uv()

   # Step 2: Create Python environment with all extras
   PrepareEnv(python_version = "3.10", extras = "all")

   # Step 3: Configure reticulate
   use_uv_env()

   # Step 4: Run environment tests
   source("tests/interactive/test_python_environment.R")
   ```

2. **M-series Mac users** (after step 1):
   ```r
   source("tests/interactive/test_m_series_mac.R")
   ```

3. **Test analysis functions** (after steps 1-2):
   ```r
   source("tests/interactive/test_python_analysis.R")
   ```

4. **Regular verification** (after updates):
   ```r
   # Quick verification
   VerifyEnv(verbose = TRUE)

   # Full verification
   source("tests/interactive/test_python_environment.R")
   ```

## Interpreting Results

### Success Indicators
- `✓` = Test passed
- `⊘` = Test skipped (not applicable or prerequisite missing)

### Failure Indicators
- `✗` = Test failed
- `⚠` = Warning (may not be critical)

### Common Issues and Fixes

#### Issue: UV not installed
```
✗ UV is not installed
```
**Fix**:
```r
install_uv()
```

#### Issue: Environment not found
```
✗ UV environment not found
```
**Fix**:
```r
PrepareEnv(extras = "all")
```

#### Issue: Python version out of range
```
⚠ Python version outside 3.10-3.13 range
```
**Fix**:
```r
PrepareEnv(force = TRUE, python_version = "3.10", extras = "all")
```

#### Issue: Missing packages
```
✗ scvelo not installed
```
**Fix**:
```r
uv_install_extras("velocity")
# or for all packages
PrepareEnv(update = TRUE, extras = "all")
```

#### Issue: Reticulate using wrong Python
```
⚠ Using different Python (not SCP environment)
```
**Fix**:
```r
use_uv_env()
```

#### Issue: NUMBA JIT failing on M-series Mac
```
✗ NUMBA test failed
```
**Expected behavior**: This is expected on M-series Macs. The SCP Python functions automatically disable JIT compilation for stability. Verify this in the M-series test script.

## Debugging Tips

### 1. Run in sections
Instead of running the entire script, execute one test section at a time:
```r
# Just run TEST 1
# Copy lines from TEST 1 section and execute
```

### 2. Check intermediate variables
After a test fails, inspect the variables:
```r
# Check Python path
print(python_path)

# Check if file exists
file.exists(python_path)

# Try running Python directly
system2(python_path, "--version")
```

### 3. Verbose output
Enable verbose mode for functions:
```r
VerifyEnv(verbose = TRUE)
ListEnv(show_details = TRUE)
```

### 4. Manual Python testing
Test Python directly from R:
```r
# Import a module
reticulate::import("numpy")

# Run Python code
reticulate::py_run_string("
import numpy as np
print(np.__version__)
")
```

### 5. Reset if needed
If something is badly broken:
```r
# Remove environment and start fresh
RemoveEnv(prompt = FALSE)
PrepareEnv(force = TRUE, python_version = "3.10", extras = "all")
use_uv_env()
```

## Test Coverage

These tests cover:
- ✓ Python 3.10, 3.11, 3.12, 3.13 compatibility
- ✓ M-series Mac (M1, M2, M3) optimizations
- ✓ UV package manager integration
- ✓ Core scientific packages (numpy, pandas, scipy, etc.)
- ✓ Single-cell packages (scanpy, anndata, scvelo)
- ✓ igraph compatibility (including PAGA known issues)
- ✓ NUMBA JIT configuration
- ✓ Threading and memory management
- ✓ Reticulate integration
- ✓ Environment management functions
- ✓ Basis validation for velocity analysis
- ✓ Data format conversions (Seurat ↔ AnnData)

## Contributing

If you add new Python functionality to SCP:
1. Add corresponding tests to these scripts
2. Document any new dependencies
3. Test on both Intel and M-series Macs if possible
4. Test across Python 3.10-3.13

## Support

If tests fail and you can't resolve the issue:
1. Check the output carefully - it usually suggests fixes
2. Run `VerifyEnv()` for diagnostic information
3. Check the main SCP documentation
4. Open an issue on GitHub with the test output

---

**Note**: These are *interactive* tests, not automated testthat tests. They're designed for debugging and manual verification. For automated CI/CD testing, use the testthat suite in `tests/testthat/`.
