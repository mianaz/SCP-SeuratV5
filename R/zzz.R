# Check for and potentially install critical dependencies
check_critical_dependencies <- function() {
  # Define critical R packages
  critical_deps <- list(
    cran_imports = c("Matrix", "Seurat", "SeuratObject", "reticulate", "rlang", "dplyr", "ggplot2"),
    cran_suggests = c("devtools", "withr"),
    bioc_imports = c("AnnotationDbi", "biomaRt", "BiocParallel", "ComplexHeatmap", "clusterProfiler"),
    bioc_suggests = c("GO.db", "GOSemSim", "slingshot")
  )
  
  # Check if each critical package is installed
  missing_deps <- list()
  for (type in names(critical_deps)) {
    missing_deps[[type]] <- character(0)
    for (pkg in critical_deps[[type]]) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        missing_deps[[type]] <- c(missing_deps[[type]], pkg)
      }
    }
  }
  
  # Check if BiocManager is installed
  has_biocmanager <- requireNamespace("BiocManager", quietly = TRUE)
  if (!has_biocmanager) {
    missing_deps$bioc_manager <- "BiocManager"
  } else {
    missing_deps$bioc_manager <- character(0)
  }
  
  # If there are missing dependencies, provide information about them
  if (length(unlist(missing_deps)) > 0) {
    cran_missing <- c(missing_deps$cran_imports, missing_deps$cran_suggests)
    bioc_missing <- c(missing_deps$bioc_imports, missing_deps$bioc_suggests)
    
    packageStartupMessage(
      "\n",
      "╔════════════════════════════════════════════════════════════╗\n",
      "║                  SCP Dependency Notice                     ║\n",
      "╚════════════════════════════════════════════════════════════╝\n"
    )
    
    if (length(missing_deps$cran_imports) > 0) {
      packageStartupMessage("Required CRAN packages missing: ", 
                           paste(missing_deps$cran_imports, collapse = ", "))
    }
    
    if (length(missing_deps$bioc_manager) > 0) {
      packageStartupMessage("BiocManager is required for Bioconductor packages")
    }
    
    if (length(missing_deps$bioc_imports) > 0) {
      packageStartupMessage("Required Bioconductor packages missing: ", 
                           paste(missing_deps$bioc_imports, collapse = ", "))
    }
    
    if (length(missing_deps$cran_suggests) > 0 || length(missing_deps$bioc_suggests) > 0) {
      packageStartupMessage("\nAdditional recommended packages:")
      if (length(missing_deps$cran_suggests) > 0) {
        packageStartupMessage("CRAN: ", paste(missing_deps$cran_suggests, collapse = ", "))
      }
      if (length(missing_deps$bioc_suggests) > 0) {
        packageStartupMessage("Bioconductor: ", paste(missing_deps$bioc_suggests, collapse = ", "))
      }
    }
    
    packageStartupMessage(
      "\nTo install all dependencies automatically, run:\n",
      "  SCP::install_all_dependencies()\n",
      "\nOr to install packages manually:\n"
    )
    
    if (length(cran_missing) > 0) {
      packageStartupMessage(
        sprintf('install.packages(c("%s"))', paste(cran_missing, collapse = '", "')),
        "\n"
      )
    }
    
    if (length(missing_deps$bioc_manager) > 0) {
      packageStartupMessage('install.packages("BiocManager")', "\n")
    }
    
    if (length(bioc_missing) > 0 && has_biocmanager) {
      packageStartupMessage(
        sprintf('BiocManager::install(c("%s"))', paste(bioc_missing, collapse = '", "')),
        "\n"
      )
    } else if (length(bioc_missing) > 0) {
      packageStartupMessage(
        "# After installing BiocManager, run:\n",
        sprintf('BiocManager::install(c("%s"))', paste(bioc_missing, collapse = '", "')),
        "\n"
      )
    }
    
    packageStartupMessage(
      "╔════════════════════════════════════════════════════════════╗\n",
      "║        Some SCP functionality may be limited                ║\n",
      "╚════════════════════════════════════════════════════════════╝\n"
    )
  }
}

.onAttach <- function(libname, pkgname) {
  options(future.globals.maxSize = Inf)
  
  # Check for critical dependencies first
  check_critical_dependencies()
  
  # Check for Seurat version
  tryCatch({
    if (requireNamespace("Seurat", quietly = TRUE)) {
      # Verify the Seurat version
      seurat_version <- packageVersion("Seurat")
      is_v5 <- as.numeric(seurat_version) >= 5
      
      if (!is_v5) {
        packageStartupMessage(
          "\n⚠️ IMPORTANT: You are running Seurat v", seurat_version, 
          ". SCP is optimized for Seurat v5+",
          "\nPlease update your Seurat objects using UpdateSeuratObject() before analysis.",
          "\nSee ?UpdateSeuratObject for details."
        )
      }
    }
  }, error = function(e) {
    # Silently fail if Seurat version check fails
  })
  
  # Check for macOS and apply necessary fixes
  tryCatch({
    if (requireNamespace("SCP", quietly = TRUE)) {
      # Check if our OS detection functions are available
      if (exists("is_macos", envir = asNamespace("SCP"), inherits = FALSE) &&
          exists("is_apple_silicon", envir = asNamespace("SCP"), inherits = FALSE) &&
          exists("apply_macos_fixes", envir = asNamespace("SCP"), inherits = FALSE)) {
        
        # Use the enhanced functions
        if (get("is_macos", envir = asNamespace("SCP"))()) {
          # Apply essential macOS fixes silently
          get("apply_macos_fixes", envir = asNamespace("SCP"))()
        }
      } else {
        # Fallback to basic detection if functions aren't available yet
        if (Sys.info()["sysname"] == "Darwin") {
          # Apply basic macOS-specific fixes silently at package load
          Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
        }
      }
    } else {
      # Fallback if namespace isn't loaded yet
      if (Sys.info()["sysname"] == "Darwin") {
        # Apply basic macOS-specific fixes silently at package load
        Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
      }
    }
  }, error = function(e) {
    # Fallback if anything goes wrong
    if (Sys.info()["sysname"] == "Darwin") {
      Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
    }
  })
  
  env <- FALSE
  # Only initialize Python environment if explicitly requested by setting SCP_env_init
  if (isTRUE(getOption("SCP_env_init", default = FALSE))) {
    tryCatch({
      conda <- find_conda()
      if (!is.null(conda)) {
        envs_dir <- reticulate:::conda_info(conda = conda)$envs_dirs[1]
        env <- env_exist(conda = conda, envname = get_envname(), envs_dir = envs_dir)
        if (isFALSE(env)) {
          packageStartupMessage("SCP python environment not found.")
        }
      } else {
        packageStartupMessage("Conda not found.")
      }
      if (isTRUE(env)) {
        Sys.unsetenv("RETICULATE_PYTHON")
        python_path <- conda_python(conda = conda)
        reticulate::use_python(python_path, required = TRUE)

        pyinfo <- utils::capture.output(reticulate::py_config())
        pyinfo_mesg <- c(
          "====================== SCP conda environment ======================",
          paste0("conda:          ", conda),
          paste0("environment:    ", paste0(envs_dir, "/", get_envname())),
          "======================== SCP python config ========================",
          pyinfo,
          "==================================================================="
        )
        invisible(lapply(pyinfo_mesg, packageStartupMessage))
        # Only load matplotlib and scanpy if needed - prevent potential crashes
        packageStartupMessage("Python environment initialized. Use SCP::run_Python() to execute Python commands.")
        packageStartupMessage("Conda path can be specified with the command `options(reticulate.conda_binary = \"/path/to/conda\")` before loading the package")
      } else {
        packageStartupMessage("If you have already created an SCP python environment using conda, you can specify the conda path by setting options(reticulate.conda_binary = \"/path/to/conda\", SCP_env_name = \"SCP_env\") before loading the package.")
      }
    }, error = function(e) {
      packageStartupMessage("Error initializing Python environment: ", e$message)
      packageStartupMessage("Python features will be disabled. Use PrepareEnv() to manually set up the environment.")
    })
  } else {
    packageStartupMessage("Python environment initialized.set options(SCP_env_init = TRUE) to disable this message.")
  }
}
