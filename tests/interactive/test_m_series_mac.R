#!/usr/bin/env Rscript
# Interactive Test Script: M-series Mac Compatibility
#
# This script tests M-series MacBook detection and optimization features.
# Run this line-by-line in RStudio for debugging.
#
# Author: SCP Development Team
# Date: 2025-01-14

# ==============================================================================
# SETUP
# ==============================================================================

cat("\n==============================================================================\n")
cat("SCP M-series Mac Compatibility Test Suite\n")
cat("==============================================================================\n\n")

# Load required packages
cat("Loading packages...\n")
library(SCP)

# ==============================================================================
# TEST 1: Platform Detection
# ==============================================================================

cat("\n--- TEST 1: Platform Detection ---\n")

# Get system information
os_name <- Sys.info()["sysname"]
machine <- Sys.info()["machine"]
release <- Sys.info()["release"]

cat("Operating System:", os_name, "\n")
cat("Machine:", machine, "\n")
cat("OS Release:", release, "\n")

# Determine if M-series Mac
is_m_series <- os_name == "Darwin" && machine == "arm64"
cat("\nM-series Mac detected:", is_m_series, "\n")

if (is_m_series) {
  cat("✓ This is an M-series MacBook\n")
  cat("  M-series specific optimizations will be applied\n")
} else {
  cat("⊘ Not an M-series Mac\n")
  if (os_name == "Darwin") {
    cat("  This is a Mac but not ARM architecture (Intel Mac?)\n")
  } else {
    cat("  This is not a Mac\n")
  }
}

# ==============================================================================
# TEST 2: Python Environment Detection
# ==============================================================================

cat("\n--- TEST 2: Python Environment ---\n")

env_exists <- uv_env_exists()
cat("UV environment exists:", env_exists, "\n")

if (env_exists) {
  # Get Python path
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

  if (file.exists(python_path)) {
    cat("Python executable:", python_path, "\n")

    # Configure reticulate
    tryCatch({
      use_uv_env()
      cat("✓ Reticulate configured\n")
    }, error = function(e) {
      cat("✗ Error configuring reticulate:", e$message, "\n")
    })
  }
} else {
  cat("⊘ Environment not found. Run PrepareEnv() first.\n")
}

# ==============================================================================
# TEST 3: M-series Specific Python Packages
# ==============================================================================

cat("\n--- TEST 3: M-series Optimized Packages ---\n")

if (env_exists && exists("python_path") && file.exists(python_path)) {

  m_packages <- list(
    "numba" = "JIT compiler (should be disabled on M-series)",
    "torch" = "Deep learning with MPS backend support",
    "numpy" = "Numerical computing (should be ARM-optimized)",
    "scipy" = "Scientific computing (should be ARM-optimized)"
  )

  cat("Checking M-series relevant packages...\n\n")

  for (pkg in names(m_packages)) {
    desc <- m_packages[[pkg]]
    cat(sprintf("Testing %s (%s):\n", pkg, desc))

    # Get package version
    version_cmd <- sprintf("import %s; print(%s.__version__)", pkg, pkg)
    version_result <- system2(python_path, c("-c", version_cmd),
                             stdout = TRUE, stderr = TRUE)

    if (is.null(attr(version_result, "status")) || attr(version_result, "status") == 0) {
      cat(sprintf("  ✓ %s version: %s\n", pkg, version_result[1]))

      # Package-specific checks
      if (pkg == "numba" && is_m_series) {
        cat("  Checking NUMBA configuration for M-series...\n")

        # Check if NUMBA_DISABLE_JIT is recommended
        numba_check <- system2(python_path,
          c("-c", "import os; print('NUMBA_DISABLE_JIT should be 1 for M-series stability')"),
          stdout = TRUE, stderr = TRUE)
        cat("  Note:", numba_check[1], "\n")
      }

      if (pkg == "torch" && is_m_series) {
        cat("  Checking PyTorch MPS backend...\n")

        # Check if MPS is available
        mps_check <- system2(python_path,
          c("-c", "import torch; print('MPS available:', torch.backends.mps.is_available())"),
          stdout = TRUE, stderr = TRUE)

        if (is.null(attr(mps_check, "status")) || attr(mps_check, "status") == 0) {
          cat(" ", mps_check[1], "\n")
        } else {
          cat("  ⊘ Could not check MPS availability\n")
        }
      }

    } else {
      cat(sprintf("  ✗ %s not installed\n", pkg))
    }
    cat("\n")
  }

} else {
  cat("Skipped (environment or Python not found)\n")
}

# ==============================================================================
# TEST 4: Environment Variables for M-series
# ==============================================================================

cat("\n--- TEST 4: M-series Environment Variables ---\n")

if (is_m_series && env_exists && exists("python_path") && file.exists(python_path)) {

  cat("Testing M-series specific environment variables...\n\n")

  # Environment variables that should be set for M-series
  m_series_env_vars <- c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "NUMBA_NUM_THREADS"
  )

  # Test setting environment variables
  test_script <- "
import os
import platform

# Detect M-series
is_m_series = platform.system() == 'Darwin' and platform.machine() == 'arm64'
print(f'M-series detected: {is_m_series}')

# Check important environment variables
env_vars = [
    'OMP_NUM_THREADS',
    'OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS',
    'VECLIB_MAXIMUM_THREADS',
    'NUMEXPR_NUM_THREADS',
    'NUMBA_NUM_THREADS',
    'NUMBA_DISABLE_JIT'
]

print('\\nEnvironment Variables:')
for var in env_vars:
    value = os.environ.get(var, 'NOT SET')
    print(f'  {var}: {value}')
"

  result <- system2(python_path, c("-c", test_script),
                   stdout = TRUE, stderr = TRUE)

  if (length(result) > 0) {
    cat(paste(result, collapse = "\n"), "\n")
  }

  cat("\nNote: These variables are typically set within SCP Python functions\n")
  cat("      (e.g., SCVELO, PAGA) to ensure M-series stability.\n")

} else {
  if (!is_m_series) {
    cat("Skipped (not an M-series Mac)\n")
  } else {
    cat("Skipped (environment not found)\n")
  }
}

# ==============================================================================
# TEST 5: NUMBA Configuration
# ==============================================================================

cat("\n--- TEST 5: NUMBA Configuration ---\n")

if (env_exists && exists("python_path") && file.exists(python_path)) {

  cat("Testing NUMBA JIT compilation...\n\n")

  # Test NUMBA with JIT enabled (default)
  numba_test_script <- "
import numba
import numpy as np

# Simple function to test NUMBA
@numba.jit(nopython=True)
def sum_array(arr):
    total = 0.0
    for i in range(arr.shape[0]):
        total += arr[i]
    return total

# Test it
test_arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
result = sum_array(test_arr)
print(f'NUMBA sum result: {result}')
print('✓ NUMBA JIT compilation works')
"

  cat("Testing NUMBA with default settings...\n")
  result <- system2(python_path, c("-c", numba_test_script),
                   stdout = TRUE, stderr = TRUE)

  if (is.null(attr(result, "status")) || attr(result, "status") == 0) {
    cat(paste(result, collapse = "\n"), "\n")
  } else {
    cat("✗ NUMBA test failed\n")
    if (is_m_series) {
      cat("  This is expected on M-series Macs. NUMBA JIT should be disabled.\n")
    }
  }

  # Test NUMBA with JIT disabled (M-series mode)
  if (is_m_series) {
    cat("\nTesting NUMBA with JIT disabled (M-series mode)...\n")

    numba_disabled_script <- "
import os
os.environ['NUMBA_DISABLE_JIT'] = '1'

import numba
import numpy as np

# Check NUMBA config
print(f'NUMBA version: {numba.__version__}')
print(f'NUMBA_DISABLE_JIT: {os.environ.get(\"NUMBA_DISABLE_JIT\", \"not set\")}')

# Simple function - will run without JIT
@numba.jit(nopython=True)
def sum_array(arr):
    total = 0.0
    for i in range(arr.shape[0]):
        total += arr[i]
    return total

# Test it (will use Python interpreter, not JIT)
test_arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
result = sum_array(test_arr)
print(f'Sum result: {result}')
print('✓ NUMBA works with JIT disabled (safer for M-series)')
"

    result <- system2(python_path, c("-c", numba_disabled_script),
                     stdout = TRUE, stderr = TRUE)

    if (length(result) > 0) {
      cat(paste(result, collapse = "\n"), "\n")
    }
  }

} else {
  cat("Skipped (environment or Python not found)\n")
}

# ==============================================================================
# TEST 6: Memory and Threading Test
# ==============================================================================

cat("\n--- TEST 6: Threading Configuration ---\n")

if (env_exists && exists("python_path") && file.exists(python_path)) {

  threading_test <- "
import os
import numpy as np
import threading

print(f'System threading info:')
print(f'  Active threads: {threading.active_count()}')

# Check environment variables for thread control
thread_vars = {
    'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS', 'not set'),
    'OPENBLAS_NUM_THREADS': os.environ.get('OPENBLAS_NUM_THREADS', 'not set'),
    'MKL_NUM_THREADS': os.environ.get('MKL_NUM_THREADS', 'not set'),
}

print('\\nThread control variables:')
for var, value in thread_vars.items():
    print(f'  {var}: {value}')

# Test numpy threading
print('\\nNumPy configuration:')
print(f'  NumPy version: {np.__version__}')

# Small computation test
arr = np.random.rand(1000, 1000)
result = np.dot(arr, arr.T)
print('  ✓ NumPy computation successful')
"

  cat("Testing threading configuration...\n")
  result <- system2(python_path, c("-c", threading_test),
                   stdout = TRUE, stderr = TRUE)

  if (length(result) > 0) {
    cat(paste(result, collapse = "\n"), "\n")
  }

  if (is_m_series) {
    cat("\nNote: For M-series Macs, single-threaded execution is recommended\n")
    cat("      for stability in scVelo and other compute-intensive tasks.\n")
  }

} else {
  cat("Skipped (environment or Python not found)\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("M-SERIES MAC TEST SUMMARY\n")
cat("==============================================================================\n\n")

if (is_m_series) {
  cat("✓ M-series Mac detected\n\n")

  cat("Recommendations for M-series:\n")
  cat("  1. NUMBA JIT should be disabled (NUMBA_DISABLE_JIT=1)\n")
  cat("  2. Single-threaded execution for stability (OMP_NUM_THREADS=1, etc.)\n")
  cat("  3. Use PyTorch with MPS backend for GPU acceleration\n")
  cat("  4. These settings are automatically applied in SCP Python functions\n\n")

  if (env_exists) {
    cat("✓ Python environment is set up\n")
    cat("  You can now use SCP's Python analysis functions\n")
    cat("  (RunSCVELO, RunPAGA, etc.) with M-series optimizations\n")
  } else {
    cat("⊘ Python environment not found\n")
    cat("  Run: PrepareEnv(extras = 'velocity')\n")
  }

} else {
  cat("⊘ Not an M-series Mac\n")
  cat("  M-series specific optimizations are not needed\n")
  cat("  Standard Python configuration will be used\n")
}

cat("\n==============================================================================\n")
