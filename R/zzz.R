.onAttach <- function(libname, pkgname) {
  options(future.globals.maxSize = Inf)

  # Check for Seurat version
  tryCatch({
    if (requireNamespace("Seurat", quietly = TRUE)) {
      # Verify the Seurat version
      seurat_version <- packageVersion("Seurat")
      is_v5 <- as.numeric(strsplit(as.character(seurat_version), "\\.")[[1]][1]) >= 5

      if (!is_v5) {
        packageStartupMessage(
          "\n[WARNING] You are running Seurat v", seurat_version,
          ". This fork is optimized for Seurat v5+",
          "\nPlease update your Seurat objects using UpdateSeuratObject()."
        )
      }
    }
  }, error = function(e) {
    # Silently fail if Seurat version check fails
  })

  # Apply macOS fixes if needed
  if (Sys.info()["sysname"] == "Darwin") {
    Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  }

  # Python environment initialization message
  packageStartupMessage(
    "SCP v", packageVersion("SCP"), " loaded.",
    "\nFor Python features, use PrepareEnv() to set up the environment."
  )
}

# UV Python environment functions

#' Get Python executable path from UV environment
#'
#' @return Path to Python executable in UV environment
#' @keywords internal
uv_python <- function() {
  pkg_dir <- system.file("", package = "SCP")
  if (pkg_dir == "") {
    pkg_dir <- getwd()
  }

  venv_path <- file.path(pkg_dir, ".venv")

  if (!dir.exists(venv_path)) {
    stop("UV environment not found. Please run PrepareEnv() first.")
  }

  # Find Python executable in venv
  if (Sys.info()["sysname"] == "Windows") {
    python_path <- file.path(venv_path, "Scripts", "python.exe")
  } else {
    python_path <- file.path(venv_path, "bin", "python")
  }

  if (!file.exists(python_path)) {
    stop("Python executable not found in UV environment")
  }

  invisible(python_path)
}

# macOS compatibility functions
is_macos <- function() {
  Sys.info()["sysname"] == "Darwin"
}

is_apple_silicon <- function() {
  if (!is_macos()) return(FALSE)

  # Check if the machine type contains ARM/M1/M2/M3 indicators
  grepl("arm64|aarch64", Sys.info()["machine"], ignore.case = TRUE)
}

apply_macos_fixes <- function() {
  # Set environment variable to prevent MKL library conflicts
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")

  # Set OpenMP threads to avoid conflicts on macOS
  Sys.setenv(OMP_NUM_THREADS = "1")

  invisible(NULL)
}