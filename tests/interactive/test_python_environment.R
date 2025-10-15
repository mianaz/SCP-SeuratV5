#!/usr/bin/env Rscript
# Interactive Test Script: Python Environment and Compatibility
#
# This script tests Python 3.10-3.13 compatibility and environment setup.
# Run this line-by-line in RStudio for debugging.
#
# Author: SCP Development Team
# Date: 2025-01-14

# ==============================================================================
# SETUP
# ==============================================================================

cat("\n==============================================================================\n")
cat("SCP Python Environment Compatibility Test Suite\n")
cat("==============================================================================\n\n")

# Load the SCP package
cat("Loading SCP package...\n")
library(SCP)

# ==============================================================================
# TEST 1: UV Installation Check
# ==============================================================================

cat("\n--- TEST 1: UV Installation ---\n")

# Check if UV is installed
uv_available <- check_uv()
cat("UV installed:", uv_available, "\n")

if (uv_available) {
  uv_version <- system2("uv", "--version", stdout = TRUE, stderr = TRUE)
  cat("UV version:", uv_version[1], "\n")
} else {
  cat("WARNING: UV not installed. Some tests will be skipped.\n")
  cat("Install with: install_uv()\n")
}

# ==============================================================================
# TEST 2: Environment Existence
# ==============================================================================

cat("\n--- TEST 2: Environment Existence ---\n")

# Check if UV environment exists
env_exists <- uv_env_exists()
cat("UV environment exists:", env_exists, "\n")

if (env_exists) {
  pkg_dir <- system.file("", package = "SCP")
  if (pkg_dir == "") {
    pkg_dir <- getwd()
  }
  venv_path <- file.path(pkg_dir, ".venv")
  cat("Environment path:", venv_path, "\n")

  # Get environment size
  env_size <- sum(file.info(list.files(venv_path, recursive = TRUE,
                                        full.names = TRUE, all.files = TRUE))$size,
                  na.rm = TRUE)
  cat("Environment size:", round(env_size / 1024^2, 2), "MB\n")
} else {
  cat("WARNING: Environment not found.\n")
  cat("Create with: PrepareEnv(extras = 'all')\n")
}

# ==============================================================================
# TEST 3: Python Executable and Version
# ==============================================================================

cat("\n--- TEST 3: Python Executable ---\n")

if (env_exists) {
  # Find Python executable
  if (Sys.info()["sysname"] == "Windows") {
    python_path <- file.path(venv_path, "Scripts", "python.exe")
  } else {
    python_path <- file.path(venv_path, "bin", "python")
  }

  python_exists <- file.exists(python_path)
  cat("Python executable exists:", python_exists, "\n")

  if (python_exists) {
    cat("Python path:", python_path, "\n")

    # Get Python version
    python_version <- system2(python_path, "--version", stdout = TRUE, stderr = TRUE)
    cat("Python version:", python_version[1], "\n")

    # Parse version to check if it's 3.10-3.13
    version_match <- regmatches(python_version[1],
                                regexec("Python (3\\.[0-9]+)", python_version[1]))
    if (length(version_match[[1]]) > 1) {
      py_version <- version_match[[1]][2]
      cat("Parsed version:", py_version, "\n")

      # Check if version is in supported range (3.10-3.13)
      version_num <- as.numeric(gsub("3\\.", "", py_version))
      is_supported <- version_num >= 10 && version_num < 14
      cat("Version in supported range (3.10-3.13):", is_supported, "\n")

      if (!is_supported) {
        cat("WARNING: Python version outside 3.10-3.13 range!\n")
      }
    }
  }
} else {
  cat("Skipped (environment not found)\n")
}

# ==============================================================================
# TEST 4: Core Python Packages Import
# ==============================================================================

cat("\n--- TEST 4: Core Python Packages ---\n")

if (env_exists && exists("python_path") && file.exists(python_path)) {
  core_packages <- c(
    "numpy", "pandas", "scipy", "matplotlib", "seaborn",
    "scikit-learn", "h5py", "numba", "anndata", "scanpy", "igraph"
  )

  cat("Testing core package imports...\n")

  package_results <- list()
  for (pkg in core_packages) {
    result <- system2(python_path,
                     c("-c", sprintf("import %s; print(%s.__version__)",
                                    pkg, pkg)),
                     stdout = TRUE,
                     stderr = TRUE)

    success <- !is.null(attr(result, "status")) && attr(result, "status") == 0 ||
               is.null(attr(result, "status"))

    if (success && length(result) > 0) {
      cat(sprintf("  ✓ %-20s %s\n", pkg, result[1]))
      package_results[[pkg]] <- list(installed = TRUE, version = result[1])
    } else {
      cat(sprintf("  ✗ %-20s NOT INSTALLED\n", pkg))
      package_results[[pkg]] <- list(installed = FALSE, version = NA)
    }
  }

  # Summary
  installed_count <- sum(sapply(package_results, function(x) x$installed))
  cat("\nPackages installed:", installed_count, "/", length(core_packages), "\n")

} else {
  cat("Skipped (environment or Python not found)\n")
}

# ==============================================================================
# TEST 5: Optional Packages (Velocity, Trajectory, etc.)
# ==============================================================================

cat("\n--- TEST 5: Optional Packages ---\n")

if (env_exists && exists("python_path") && file.exists(python_path)) {
  optional_packages <- list(
    "Velocity" = c("scvelo", "loompy"),
    "Trajectory" = c("phate", "palantir", "cellrank"),
    "Deep Learning" = c("torch"),
    "Spatial" = c("squidpy")
  )

  cat("Testing optional package groups...\n\n")

  for (group_name in names(optional_packages)) {
    cat(sprintf("%s:\n", group_name))
    pkgs <- optional_packages[[group_name]]

    for (pkg in pkgs) {
      result <- system2(python_path,
                       c("-c", sprintf("import %s; print(%s.__version__)",
                                      pkg, pkg)),
                       stdout = TRUE,
                       stderr = TRUE)

      success <- !is.null(attr(result, "status")) && attr(result, "status") == 0 ||
                 is.null(attr(result, "status"))

      if (success && length(result) > 0) {
        cat(sprintf("  ✓ %-20s %s\n", pkg, result[1]))
      } else {
        cat(sprintf("  ⊘ %-20s not installed\n", pkg))
      }
    }
    cat("\n")
  }

} else {
  cat("Skipped (environment or Python not found)\n")
}

# ==============================================================================
# TEST 6: Reticulate Configuration
# ==============================================================================

cat("\n--- TEST 6: Reticulate Configuration ---\n")

# Try to configure reticulate
if (env_exists && exists("python_path") && file.exists(python_path)) {
  cat("Configuring reticulate to use SCP environment...\n")

  tryCatch({
    use_uv_env()
    cat("✓ Reticulate configured successfully\n")

    # Get Python config
    py_config <- reticulate::py_config()
    cat("Python binary:", py_config$python, "\n")
    cat("Python version:", py_config$version, "\n")
    cat("numpy:", py_config$numpy, "\n")

    # Verify it's using the correct Python
    if (normalizePath(py_config$python, mustWork = FALSE) ==
        normalizePath(python_path, mustWork = FALSE)) {
      cat("✓ Using SCP UV environment\n")
    } else {
      cat("⚠ WARNING: Using different Python!\n")
      cat("  Expected:", python_path, "\n")
      cat("  Actual:", py_config$python, "\n")
    }

  }, error = function(e) {
    cat("✗ Error configuring reticulate:", e$message, "\n")
  })

} else {
  cat("Skipped (environment or Python not found)\n")
}

# ==============================================================================
# TEST 7: Environment Management Functions
# ==============================================================================

cat("\n--- TEST 7: Environment Management Functions ---\n")

# Test ListEnv()
cat("\nTesting ListEnv()...\n")
cat("----------------------------------------\n")
env_info <- ListEnv(show_details = TRUE)

# Test VerifyEnv()
cat("\nTesting VerifyEnv()...\n")
cat("----------------------------------------\n")
env_valid <- VerifyEnv(verbose = TRUE)

if (env_valid) {
  cat("\n✓ Environment verification PASSED\n")
} else {
  cat("\n✗ Environment verification FAILED\n")
}

# ==============================================================================
# TEST 8: Python 3.12 Specific Features
# ==============================================================================

cat("\n--- TEST 8: Python 3.12 Compatibility ---\n")

if (env_exists && exists("python_path") && file.exists(python_path)) {

  # Test match statement (Python 3.10+)
  cat("\nTesting Python 3.10+ match statement...\n")
  match_test <- system2(python_path,
    c("-c", "x = 1; match x:\n    case 1:\n        print('match works')"),
    stdout = TRUE,
    stderr = TRUE)

  if (!is.null(attr(match_test, "status")) && attr(match_test, "status") != 0) {
    cat("✗ Match statement failed (not Python 3.10+?)\n")
  } else {
    cat("✓ Match statement works (Python 3.10+)\n")
  }

  # Test type hints (should work in all versions)
  cat("\nTesting type hints...\n")
  type_test <- system2(python_path,
    c("-c", "from typing import Optional\ndef foo(x: int) -> Optional[str]:\n    return str(x)\nprint(foo(42))"),
    stdout = TRUE,
    stderr = TRUE)

  if (!is.null(attr(type_test, "status")) && attr(type_test, "status") != 0) {
    cat("✗ Type hints failed\n")
  } else {
    cat("✓ Type hints work\n")
  }

} else {
  cat("Skipped (environment or Python not found)\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("TEST SUMMARY\n")
cat("==============================================================================\n\n")

test_results <- list(
  "UV Installed" = uv_available,
  "Environment Exists" = env_exists,
  "Python Executable" = if (env_exists && exists("python_path")) file.exists(python_path) else FALSE,
  "Environment Valid" = if (exists("env_valid")) env_valid else FALSE
)

for (test_name in names(test_results)) {
  status <- if (test_results[[test_name]]) "✓ PASS" else "✗ FAIL"
  cat(sprintf("%-25s %s\n", test_name, status))
}

# Overall result
all_passed <- all(unlist(test_results))
cat("\n")
if (all_passed) {
  cat("✓✓✓ ALL TESTS PASSED ✓✓✓\n")
  cat("Your Python environment is ready for SCP!\n")
} else {
  cat("✗✗✗ SOME TESTS FAILED ✗✗✗\n")
  cat("Please review the output above for issues.\n")
  cat("\nQuick fixes:\n")
  cat("  1. Install UV: install_uv()\n")
  cat("  2. Create environment: PrepareEnv(extras = 'all')\n")
  cat("  3. Configure reticulate: use_uv_env()\n")
}

cat("\n==============================================================================\n")
