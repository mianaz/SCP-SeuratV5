.onAttach <- function(libname, pkgname) {
  options(future.globals.maxSize = Inf)
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
    packageStartupMessage("Python environment not initialized. Run PrepareEnv() to set up Python, or set options(SCP_env_init = TRUE) before loading SCP.")
  }
}
