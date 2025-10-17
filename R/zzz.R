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
          "\n",
          "\033[33m", # Yellow color for warning
          "\u26A0 NOTICE: Running Seurat v", seurat_version, " (v4)",
          "\033[39m", # Reset color
          "\n  SCP supports both Seurat v4 and v5, but v5 is recommended for:",
          "\n  - Better performance with large datasets",
          "\n  - Improved memory efficiency",
          "\n  - Latest integration methods",
          "\n",
          "\n  To upgrade: install.packages('Seurat')",
          "\n  Then update objects: UpdateSeuratObject(your_object)",
          "\n"
        )
      } else {
        packageStartupMessage(
          "\n",
          "\033[32m", # Green color
          "\u2713 Running Seurat v5 - optimal performance enabled",
          "\033[39m", # Reset color
          "\n"
        )
      }
    }
  }, error = function(e) {
    # Silently fail if Seurat version check fails
  })

  # Set environment variables to suppress known warnings
  # These need to be set before Python is initialized
  Sys.setenv(PYTHONWARNINGS = "ignore::DeprecationWarning:pkg_resources")
  Sys.setenv(OMP_NUM_THREADS = "1")  # Avoid nested parallelism issues

  # Apply macOS fixes if needed
  if (Sys.info()["sysname"] == "Darwin") {
    Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  }

  # Show appropriate message based on environment status
  if (uv_env_exists()) {
    packageStartupMessage(
      "SCPNext v", packageVersion("SCPNext"), " loaded.",
      "\nPython environment found. Run use_uv_env() to activate."
    )
  } else {
    packageStartupMessage(
      "SCPNext v", packageVersion("SCPNext"), " loaded.",
      "\nFor Python features, use PrepareEnv() to set up the environment."
    )
  }
}

