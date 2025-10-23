# Python-Based Integration Methods
# Contains scVI and SCANVI integration functions
# These methods require Python environment with scvi-tools package

#' Setup Python Environment for scvi-tools (Internal)
#'
#' Common setup for Python-based integration methods (scVI, SCANVI).
#' Checks and configures Python environment, imports required modules.
#'
#' @param num_threads Integer. Number of threads for scvi-tools
#'
#' @return List with scvi and scipy Python modules
#' @keywords internal
.setup_scvi_env <- function(num_threads = 8) {
  # Windows-specific JAX installation for scvi-tools
  if (.Platform$OS.type == "windows" && !exist_Python_pkgs(packages = "scvi-tools")) {
    suppressWarnings(system2(
      command = get_uv_python_path(),
      args = "-m pip install jax[cpu]===0.3.20 -f https://whls.blob.core.windows.net/unstable/index.html --use-deprecated legacy-resolver",
      stdout = TRUE
    ))
  }

  # Activate Python environment
  use_uv_env()

  # Check for scvi-tools
  if (!reticulate::py_module_available("scvi")) {
    stop("Python module 'scvi-tools' is required. Install with: uv_install(packages = 'scvi-tools')")
  }

  # Import modules
  scvi <- import("scvi")
  scipy <- import("scipy")

  # Configure threads
  scvi$settings$num_threads <- as.integer(num_threads)

  return(list(scvi = scvi, scipy = scipy))
}


#' scVI Integration
#'
#' Integration using scVI (single-cell Variational Inference) or PeakVI for ATAC-seq data.
#' Requires Python environment with scvi-tools package.
#'
#' @inheritParams Integration_SCP
#' @param scVI_dims_use Integer vector. Dimensions from scVI latent space to use for downstream analysis.
#'   If NULL, uses all dimensions.
#' @param model Character. Model to use: "SCVI" for RNA-seq or "PEAKVI" for ATAC-seq. Default is "SCVI".
#' @param SCVI_params List. Parameters passed to scvi.model.SCVI(). Default is empty list.
#' @param PEAKVI_params List. Parameters passed to scvi.model.PEAKVI(). Default is empty list.
#' @param num_threads Integer. Number of threads for scvi-tools. Default is 8.
#'
#' @return A Seurat object with scVI integration results
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject CreateAssayObject CreateDimReducObject LayerData LayerData<-
#' @importFrom reticulate import
#' @export
scVI_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                           do_normalization = NULL, normalization_method = "LogNormalize",
                           do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst",
                           nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                           scVI_dims_use = NULL,
                           nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3),
                           nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                           neighbor_metric = "euclidean", neighbor_k = 20L,
                           cluster_algorithm = "louvain", cluster_resolution = 0.6,
                           model = "SCVI", SCVI_params = list(), PEAKVI_params = list(),
                           num_threads = 8, seed = 11) {

  # Validate parameters using framework helper (no linear_reduction for scVI)
  params <- .validate_reduction_params(
    linear_reduction = "pca",  # dummy value, not used by scVI
    nonlinear_reduction = nonlinear_reduction,
    cluster_algorithm = cluster_algorithm,
    linear_reduction_dims = 50,  # dummy value
    linear_reduction_dims_use = NULL,
    srtMerge = srtMerge
  )
  cluster_algorithm <- params$cluster_algorithm
  cluster_algorithm_index <- params$cluster_algorithm_index
  nonlinear_reduction <- params$nonlinear_reduction

  # Setup Python environment
  py_modules <- .setup_scvi_env(num_threads)
  scvi <- py_modules$scvi
  scipy <- py_modules$scipy

  set.seed(seed)

  # Input validation
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }

  # Save raw object for appending
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }

  # Prepare data
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  } else {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
  }

  # Convert to AnnData
  adata <- srt_to_adata(srtMerge, features = HVF, assay_X = DefaultAssay(srtMerge),
                        assay_layers = NULL, verbose = FALSE)
  adata[["X"]] <- scipy$sparse$csr_matrix(adata[["X"]])

  # Train model
  if (model == "SCVI") {
    scvi$model$SCVI$setup_anndata(adata, batch_key = batch)
    params_list <- list(adata = adata)
    for (nm in names(SCVI_params)) {
      params_list[[nm]] <- SCVI_params[[nm]]
    }
    trained_model <- invoke(.fn = scvi$model$SCVI, .args = params_list)
    trained_model$train()

    srtIntegrated <- srtMerge
    srtMerge <- NULL
    corrected <- t(as_matrix(trained_model$get_normalized_expression()))
    srtIntegrated[["scVIcorrected"]] <- CreateAssayObject(counts = corrected)
    DefaultAssay(srtIntegrated) <- "scVIcorrected"
    VariableFeatures(srtIntegrated[["scVIcorrected"]]) <- HVF

  } else if (model == "PEAKVI") {
    message("Assay is ChromatinAssay. Using PeakVI workflow.")
    scvi$model$PEAKVI$setup_anndata(adata, batch_key = batch)
    params_list <- list(adata = adata)
    for (nm in names(PEAKVI_params)) {
      params_list[[nm]] <- PEAKVI_params[[nm]]
    }
    trained_model <- invoke(.fn = scvi$model$PEAKVI, .args = params_list)
    trained_model$train()

    srtIntegrated <- srtMerge
    srtMerge <- NULL
  }

  # Extract latent representation
  latent <- as_matrix(trained_model$get_latent_representation())
  rownames(latent) <- colnames(srtIntegrated)
  colnames(latent) <- paste0("scVI_", seq_len(ncol(latent)))
  srtIntegrated[["scVI"]] <- CreateDimReducObject(embeddings = latent, key = "scVI_",
                                                  assay = DefaultAssay(srtIntegrated))

  if (is.null(scVI_dims_use)) {
    scVI_dims_use <- 1:ncol(srtIntegrated[["scVI"]]@cell.embeddings)
  }

  # Clustering using framework helper - note: scVI doesn't use linear_reduction name
  srtIntegrated <- .post_integration_clustering(
    srt = srtIntegrated,
    prefix = "scVI",
    reduction_name = "scVI",
    dims = scVI_dims_use,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    HVF = HVF,
    cluster_algorithm = cluster_algorithm,
    linear_reduction = ""  # scVI doesn't append linear_reduction to cluster name
  )

  # Nonlinear reductions using framework helper
  srtIntegrated <- .post_integration_reductions(
    srt = srtIntegrated,
    prefix = "scVI",
    reduction_use = "scVI",
    reduction_dims = scVI_dims_use,
    graph_use = "scVI_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed
  )

  # Finalize
  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["scVI_HVF"]] <- HVF

  # Append if requested
  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(
      srt_raw = srtMerge_raw, srt_append = srtIntegrated,
      pattern = paste0(assay, "|scVI|Default_reduction"),
      overwrite = TRUE, verbose = FALSE
    )
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}


#' SCANVI Integration
#'
#' Integration using SCANVI (single-cell ANnotation using Variational Inference).
#' Semi-supervised method that uses cell type labels for improved integration.
#' Requires Python environment with scvi-tools package.
#'
#' @inheritParams Integration_SCP
#' @param SCANVI_dims_use Integer vector. Dimensions from SCANVI latent space to use for downstream analysis.
#'   If NULL, uses all dimensions.
#' @param labels_key Character. Column name in metadata containing cell labels for semi-supervised training.
#'   Required parameter.
#' @param unlabeled_category Character. Category in labels_key representing unlabeled cells.
#'   Default is "Unknown".
#' @param SCVI_params List. Parameters passed to scvi.model.SCVI() for pretraining. Default is empty list.
#' @param SCANVI_params List. Parameters passed to SCANVI.from_scvi_model(). Default is empty list.
#' @param num_threads Integer. Number of threads for scvi-tools. Default is 8.
#'
#' @return A Seurat object with SCANVI integration results
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom SeuratObject CreateAssayObject CreateDimReducObject LayerData LayerData<-
#' @importFrom reticulate import
#' @export
SCANVI_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                             do_normalization = NULL, normalization_method = "LogNormalize",
                             do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst",
                             nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                             SCANVI_dims_use = NULL,
                             labels_key = NULL, unlabeled_category = "Unknown",
                             nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3),
                             nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                             neighbor_metric = "euclidean", neighbor_k = 20L,
                             cluster_algorithm = "louvain", cluster_resolution = 0.6,
                             SCVI_params = list(), SCANVI_params = list(),
                             num_threads = 8, seed = 11) {

  # Validate parameters using framework helper
  params <- .validate_reduction_params(
    linear_reduction = "pca",  # dummy value, not used by SCANVI
    nonlinear_reduction = nonlinear_reduction,
    cluster_algorithm = cluster_algorithm,
    linear_reduction_dims = 50,  # dummy value
    linear_reduction_dims_use = NULL,
    srtMerge = srtMerge
  )
  cluster_algorithm <- params$cluster_algorithm
  cluster_algorithm_index <- params$cluster_algorithm_index
  nonlinear_reduction <- params$nonlinear_reduction

  # Validate labels_key
  if (is.null(labels_key) || length(labels_key) != 1 || !is.character(labels_key)) {
    stop("'labels_key' must be a single character string specifying the metadata column with labels.")
  }

  # Setup Python environment
  py_modules <- .setup_scvi_env(num_threads)
  scvi <- py_modules$scvi
  scipy <- py_modules$scipy

  set.seed(seed)

  # Input validation
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }

  # Save raw object for appending
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }

  # Prepare data
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  } else {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
  }

  # Validate labels_key exists in metadata
  if (!labels_key %in% colnames(srtMerge@meta.data)) {
    stop(paste0("labels_key '", labels_key, "' was not found in the merged object's metadata."))
  }
  label_values <- srtMerge@meta.data[[labels_key]]
  if (!unlabeled_category %in% unique(label_values)) {
    stop(paste0("unlabeled_category '", unlabeled_category, "' was not found in column '", labels_key, "'."))
  }

  # Convert to AnnData
  adata <- srt_to_adata(srtMerge, features = HVF, assay_X = DefaultAssay(srtMerge),
                        assay_layers = NULL, verbose = FALSE)
  adata[["X"]] <- scipy$sparse$csr_matrix(adata[["X"]])

  # Train SCVI model first (pretraining)
  scvi$model$SCVI$setup_anndata(adata, batch_key = batch, labels_key = labels_key)
  scvi_params <- list(adata = adata)
  for (nm in names(SCVI_params)) {
    scvi_params[[nm]] <- SCVI_params[[nm]]
  }
  scvi_model <- invoke(.fn = scvi$model$SCVI, .args = scvi_params)
  scvi_model$train()

  # Train SCANVI model from SCVI
  scanvi_args <- list(
    scvi_model = scvi_model,
    unlabeled_category = unlabeled_category
  )
  for (nm in names(SCANVI_params)) {
    scanvi_args[[nm]] <- SCANVI_params[[nm]]
  }
  model <- invoke(.fn = scvi$model$SCANVI$from_scvi_model, .args = scanvi_args)
  model$train()

  srtIntegrated <- srtMerge
  srtMerge <- NULL

  # Get corrected expression
  corrected <- t(as_matrix(model$get_normalized_expression()))
  srtIntegrated[["SCANVIcorrected"]] <- CreateAssayObject(counts = corrected)
  DefaultAssay(srtIntegrated) <- "SCANVIcorrected"
  VariableFeatures(srtIntegrated[["SCANVIcorrected"]]) <- HVF

  # Extract latent representation
  latent <- as_matrix(model$get_latent_representation())
  rownames(latent) <- colnames(srtIntegrated)
  colnames(latent) <- paste0("SCANVI_", seq_len(ncol(latent)))
  srtIntegrated[["SCANVI"]] <- CreateDimReducObject(embeddings = latent, key = "SCANVI_",
                                                    assay = DefaultAssay(srtIntegrated))

  if (is.null(SCANVI_dims_use)) {
    SCANVI_dims_use <- 1:ncol(srtIntegrated[["SCANVI"]]@cell.embeddings)
  }

  # Clustering using framework helper
  srtIntegrated <- .post_integration_clustering(
    srt = srtIntegrated,
    prefix = "SCANVI",
    reduction_name = "SCANVI",
    dims = SCANVI_dims_use,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    HVF = HVF,
    cluster_algorithm = cluster_algorithm,
    linear_reduction = ""  # SCANVI doesn't append linear_reduction to cluster name
  )

  # Nonlinear reductions using framework helper
  srtIntegrated <- .post_integration_reductions(
    srt = srtIntegrated,
    prefix = "SCANVI",
    reduction_use = "SCANVI",
    reduction_dims = SCANVI_dims_use,
    graph_use = "SCANVI_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed
  )

  # Finalize
  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["SCANVI_HVF"]] <- HVF

  # Append if requested
  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(
      srt_raw = srtMerge_raw, srt_append = srtIntegrated,
      pattern = paste0(assay, "|SCANVI|Default_reduction"),
      overwrite = TRUE, verbose = FALSE
    )
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}
