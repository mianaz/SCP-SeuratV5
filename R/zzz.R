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

# Simplified Python environment functions
conda_python <- function(envname = "SCP_env", conda = "auto", ...) {
  envs <- reticulate::conda_list(conda = conda, ...)
  if (!envname %in% envs[["name"]]) {
    stop("conda environment ", envname, " is not available")
  }
  python_bin <- envs[envs[["name"]] == envname, ][["python"]][1]

  path_cut <- python_bin
  while (!file.exists(path_cut)) {
    path_cut_new <- dirname(path_cut)
    if (path_cut_new == path_cut) {
      break
    }
    path_cut <- path_cut_new
  }

  if (file.exists(path_cut)) {
    python_bin <- path_cut
  }

  invisible(python_bin)
}

get_envname <- function() {
  envname <- getOption("SCP_env_name", default = "SCP_env")
  return(envname)
}

find_conda <- function() {
  tryCatch(
    {
      conda <- reticulate::conda_binary()
      conda <- conda[which.min(nchar(conda))]
      return(conda)
    },
    error = function(error) {
      return(NULL)
    }
  )
}

env_exist <- function(conda = "auto", envname = "SCP_env", envs_dir = NULL) {
  tryCatch(
    {
      if (is.null(envs_dir)) {
        conda_envs <- reticulate::conda_list(conda = conda)
        env_exist <- envname %in% conda_envs[["name"]]
      } else {
        env_exist <- dir.exists(paste0(envs_dir, "/", envname))
      }
    },
    error = function(e) {
      env_exist <- FALSE
    }
  )
  return(env_exist)
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