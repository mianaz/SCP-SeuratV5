# Main Seurat Integration Function
# Supports both V4 (list-based) and V5 (layer-based) workflows with auto-detection

#' Seurat Integration with Auto-detection
#'
#' Main Seurat integration function that automatically detects Seurat version
#' and uses the appropriate workflow (V4 list-based or V5 layer-based).
#'
#' @inheritParams Integration_SCP
#' @param use_v5_workflow Logical or NULL. If NULL, auto-detect Seurat version.
#'   If TRUE, force V5 layer-based workflow. If FALSE, force V4 list-based workflow.
#' @param integration_method For V5 workflow, integration method to use. Options:
#'   "CCAIntegration", "RPCAIntegration", "HarmonyIntegration", "FastMNNIntegration", "JointPCAIntegration".
#'   For V4 workflow, uses traditional FindIntegrationAnchors + IntegrateData.
#' @param FindIntegrationAnchors_params List of parameters for FindIntegrationAnchors (V4 only)
#' @param IntegrateData_params List of parameters for IntegrateData (V4 only)
#' @param IntegrateEmbeddings_params List of parameters for IntegrateEmbeddings (V4 only)
#' @param IntegrateLayers_params List of parameters for IntegrateLayers (V5 only)
#'
#' @return A \code{Seurat} object with integrated data.
#'
#' @importFrom Seurat FindIntegrationAnchors IntegrateData IntegrateEmbeddings
#' @importFrom Seurat IntegrateLayers CCAIntegration RPCAIntegration HarmonyIntegration
#' @importFrom Seurat FastMNNIntegration JointPCAIntegration
#' @importFrom SeuratObject JoinLayers
#' @export
Seurat_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                                   do_normalization = NULL, normalization_method = "LogNormalize",
                                   do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst",
                                   nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                                   do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                                   scale_within_batch = FALSE,
                                   linear_reduction = "pca", linear_reduction_dims = 50,
                                   linear_reduction_dims_use = NULL, linear_reduction_params = list(),
                                   force_linear_reduction = FALSE,
                                   nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3),
                                   nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                                   neighbor_metric = "euclidean", neighbor_k = 20L,
                                   cluster_algorithm = "louvain", cluster_resolution = 0.6,
                                   use_v5_workflow = NULL, integration_method = "CCAIntegration",
                                   FindIntegrationAnchors_params = list(), IntegrateData_params = list(),
                                   IntegrateEmbeddings_params = list(), IntegrateLayers_params = list(),
                                   verbose = TRUE, seed = 11) {
  
  # Input validation
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  
  set.seed(seed)
  
  # Detect workflow
  workflow_info <- .detect_seurat_workflow(
    srt = if (!is.null(srtList)) srtList else srtMerge,
    use_v5_workflow = use_v5_workflow
  )
  
  cat(paste0("[", Sys.time(), "]", " Using Seurat ", workflow_info$version, " ", 
             workflow_info$workflow, "-based workflow", 
             if (workflow_info$auto_detected) " (auto-detected)" else " (forced)", "...\n"))
  
  # Prepare inputs based on workflow
  inputs <- .prepare_integration_inputs(srtList, srtMerge, batch, assay, workflow_info$workflow)
  srtList <- inputs$srtList
  srtMerge <- inputs$srtMerge
  
  # Get method-specific parameters
  method_params <- .get_integration_method_params("Seurat", workflow_info$workflow)
  
  # Check if method supports the detected workflow
  if (!is.null(method_params$workflow) && method_params$workflow != workflow_info$workflow) {
    cat(paste0("[", Sys.time(), "]", " Method doesn't support ", workflow_info$workflow, 
               " workflow, falling back to ", method_params$workflow, "...\n"))
    workflow_info$workflow <- method_params$workflow
    inputs <- .prepare_integration_inputs(srtList, srtMerge, batch, assay, workflow_info$workflow)
    srtList <- inputs$srtList
    srtMerge <- inputs$srtMerge
  }
  
  # Execute appropriate workflow
  if (workflow_info$workflow == "layers") {
    # V5 layer-based workflow
    srtIntegrated <- .seurat_integrate_v5_layers(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, normalization_method = normalization_method,
      do_HVF_finding = do_HVF_finding, HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      do_scaling = do_scaling, vars_to_regress = vars_to_regress,
      regression_model = regression_model, scale_within_batch = scale_within_batch,
      linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims,
      linear_reduction_dims_use = linear_reduction_dims_use,
      linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      nonlinear_reduction = nonlinear_reduction, nonlinear_reduction_dims = nonlinear_reduction_dims,
      nonlinear_reduction_params = nonlinear_reduction_params, force_nonlinear_reduction = force_nonlinear_reduction,
      neighbor_metric = neighbor_metric, neighbor_k = neighbor_k,
      cluster_algorithm = cluster_algorithm, cluster_resolution = cluster_resolution,
      integration_method = integration_method, IntegrateLayers_params = IntegrateLayers_params,
      verbose = verbose, seed = seed
    )
  } else {
    # V4 list-based workflow
    srtIntegrated <- .seurat_integrate_v4_list(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, normalization_method = normalization_method,
      do_HVF_finding = do_HVF_finding, HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      do_scaling = do_scaling, vars_to_regress = vars_to_regress, 
      regression_model = regression_model, scale_within_batch = scale_within_batch,
      linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims,
      linear_reduction_dims_use = linear_reduction_dims_use, 
      linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      nonlinear_reduction = nonlinear_reduction, nonlinear_reduction_dims = nonlinear_reduction_dims,
      nonlinear_reduction_params = nonlinear_reduction_params, force_nonlinear_reduction = force_nonlinear_reduction,
      neighbor_metric = neighbor_metric, neighbor_k = neighbor_k,
      cluster_algorithm = cluster_algorithm, cluster_resolution = cluster_resolution,
      FindIntegrationAnchors_params = FindIntegrationAnchors_params,
      IntegrateData_params = IntegrateData_params,
      IntegrateEmbeddings_params = IntegrateEmbeddings_params,
      verbose = verbose, seed = seed
    )
  }
  
  # Post-process results
  srtIntegrated <- .post_process_integration(srtIntegrated, workflow_info$workflow, assay)

  # TODO: Append functionality disabled - Integration_SCP_append function not implemented yet
  # if (isTRUE(append)) {
  #   srtIntegrated <- Integration_SCP_append(
  #     srtMerge = srtIntegrated,
  #     srtList = if (!is.null(srtList)) srtList else .convert_layers_to_list(srtIntegrated, batch, assay),
  #     batch = batch,
  #     integration_method = "Seurat",
  #     seed = seed
  #   )
  # }

  return(srtIntegrated)
}

#' V5 Layer-based Seurat Integration (Internal)
#'
#' @inheritParams Seurat_integrate
#' @return Seurat object with integrated layers
.seurat_integrate_v5_layers <- function(srtMerge, batch, assay, do_normalization, normalization_method,
                                      do_HVF_finding, HVF_source, HVF_method, nHVF, HVF_min_intersection, HVF,
                                      do_scaling, vars_to_regress, regression_model, scale_within_batch,
                                      linear_reduction, linear_reduction_dims, linear_reduction_dims_use,
                                      linear_reduction_params, force_linear_reduction,
                                      nonlinear_reduction, nonlinear_reduction_dims, nonlinear_reduction_params,
                                      force_nonlinear_reduction, neighbor_metric, neighbor_k,
                                      cluster_algorithm, cluster_resolution, integration_method,
                                      IntegrateLayers_params, verbose = TRUE, seed) {
  
  # Verify object is Seurat v5
  if (!IsSeurat5(srtMerge)) {
    stop("V5 layer workflow requires a Seurat v5 object. Use UpdateSeuratObject() to convert.")
  }
  
  # Check inputs (skip wasteful list processing for v5 layer workflow)
  checked <- check_srtMerge(
    srtMerge = srtMerge, batch = batch, assay = assay,
    do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = HVF_source, HVF_method = HVF_method,
    nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
    vars_to_regress = vars_to_regress, seed = seed,
    skip_list_processing = TRUE  # Skip redundant split->merge operations
  )
  srtMerge <- checked[["srtMerge"]]
  HVF <- checked[["HVF"]]
  assay <- checked[["assay"]]
  type <- checked[["type"]]
  
  # Split layers by batch if not already split
  layers <- SeuratObject::Layers(srtMerge[[assay]])
  data_layers <- layers[grepl("^data\\.", layers)]

  if (length(data_layers) <= 1) {
    # Not split yet, split now
    cat(paste0("[", Sys.time(), "]", " Splitting layers by batch: ", batch, "...\n"))
    srtMerge[[assay]] <- split(srtMerge[[assay]], f = srtMerge[[batch, drop = TRUE]])
  } else {
    cat(paste0("[", Sys.time(), "]", " Layers already split by batch (", length(data_layers), " batches)...\n"))
  }
  
  # Normalize per-layer if needed
  if (isTRUE(do_normalization) || is.null(do_normalization)) {
    if (normalization_method == "LogNormalize") {
      cat(paste0("[", Sys.time(), "]", " Normalizing data (LogNormalize) across layers...\n"))
      srtMerge <- NormalizeData(srtMerge, assay = assay, normalization.method = "LogNormalize", verbose = FALSE)
    } else if (normalization_method == "SCT") {
      cat(paste0("[", Sys.time(), "]", " Normalizing data (SCTransform) across layers...\n"))
      require_packages("glmGamPoi")
      srtMerge <- SCTransform(srtMerge, assay = assay, vst.flavor = "v2", verbose = FALSE)
    }
  }
  
  # Find variable features
  if (isTRUE(do_HVF_finding)) {
    cat(paste0("[", Sys.time(), "]", " Finding variable features across layers...\n"))
    srtMerge <- FindVariableFeatures(srtMerge, assay = assay, selection.method = HVF_method, nfeatures = nHVF, verbose = FALSE)
    HVF <- VariableFeatures(srtMerge)
  }
  
  # Scale data using optimized V5 function
  # Note: PCA requires scaled data, so we always scale before running PCA
  if (isTRUE(do_scaling) || is.null(do_scaling) || !.check_scaled_v5(srtMerge, HVF, assay)) {
    scale_start <- Sys.time()
    cat(paste0("[", scale_start, "]", " Scaling data across layers (optimized)...\n"))
    if (verbose) {
      layers <- SeuratObject::Layers(srtMerge[[assay]])
      data_layers <- layers[grepl("^data\\.", layers)]
      cat(paste0("[", scale_start, "]", "   Scaling ", length(HVF), " variable features across ", length(data_layers), " split layers...\n"))
      cat(paste0("[", scale_start, "]", "   Number of cells: ", ncol(srtMerge), "\n"))
    }
    srtMerge <- ScaleData(
      srtMerge,
      features = HVF,
      assay = assay,
      vars.to.regress = vars_to_regress,
      verbose = FALSE
    )
    scale_end <- Sys.time()
    scale_time <- difftime(scale_end, scale_start, units = "secs")
    cat(paste0("[", scale_end, "]", "   ScaleData completed in ", round(scale_time, 2), " seconds\n"))
  }

  # Run PCA
  cat(paste0("[", Sys.time(), "]", " Running PCA on the merged object...\n"))
  srtMerge <- RunPCA(srtMerge, assay = assay, features = HVF, npcs = linear_reduction_dims, verbose = FALSE)
  
  # Map integration method to function
  method_function <- switch(integration_method,
    "CCAIntegration" = CCAIntegration,
    "RPCAIntegration" = RPCAIntegration,
    "HarmonyIntegration" = HarmonyIntegration,
    "FastMNNIntegration" = FastMNNIntegration,
    "JointPCAIntegration" = JointPCAIntegration
  )
  
  # Integrate layers using IntegrateLayers
  cat(paste0("[", Sys.time(), "]", " Performing integration using ", integration_method, "...\n"))
  params <- list(
    object = srtMerge,
    method = method_function,
    orig.reduction = "pca",
    new.reduction = "integrated.dr",
    verbose = FALSE
  )
  
  # Merge user parameters
  for (nm in names(IntegrateLayers_params)) {
    params[[nm]] <- IntegrateLayers_params[[nm]]
  }
  
  srtIntegrated <- invoke(.fn = IntegrateLayers, .args = params)

  # Store HVF
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Seurat_HVF"]] <- HVF

  # Standard Seurat v5 workflow: Perform downstream analysis on integrated.dr
  # This follows official Seurat v5 best practices for layer-based integration

  # Determine dimensions to use
  dims_use <- if (!is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use
  } else {
    1:min(30, ncol(srtIntegrated@reductions[["integrated.dr"]]))
  }

  # Setup cluster algorithm
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  # Clustering on integrated reduction
  if (!is.null(cluster_resolution)) {
    if (verbose) cat(paste0("[", Sys.time(), "]", " Finding neighbors on integrated.dr...\n"))
    srtIntegrated <- FindNeighbors(
      object = srtIntegrated,
      reduction = "integrated.dr",
      dims = dims_use,
      annoy.metric = neighbor_metric,
      k.param = neighbor_k,
      graph.name = c("Seurat_integrated_nn", "Seurat_integrated_snn"),
      verbose = FALSE
    )

    if (verbose) cat(paste0("[", Sys.time(), "]", " Finding clusters (", cluster_algorithm, ")...\n"))
    srtIntegrated <- FindClusters(
      object = srtIntegrated,
      resolution = cluster_resolution,
      algorithm = cluster_algorithm_index,
      method = "igraph",
      graph.name = "Seurat_integrated_snn",
      verbose = FALSE
    )

    # Reorder clusters
    if (verbose) cat(paste0("[", Sys.time(), "]", " Reordering clusters...\n"))
    srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters")
    srtIntegrated[[paste0("Seurat", "clusters")]] <- Idents(srtIntegrated)
    srtIntegrated[["seurat_clusters"]] <- NULL
  }

  # UMAP on integrated reduction
  if (!is.null(nonlinear_reduction)) {
    for (nr in nonlinear_reduction) {
      if (!nr %in% c("umap", "tsne")) {
        warning(paste0("Nonlinear reduction '", nr, "' not supported in v5 workflow. Skipping."))
        next
      }

      for (n in nonlinear_reduction_dims) {
        if (verbose) cat(paste0("[", Sys.time(), "]", " Computing ", toupper(nr), " (", n, "D) on integrated.dr...\n"))

        reduction_name <- paste0("Seurat", toupper(nr), n, "D")

        if (nr == "umap") {
          srtIntegrated <- RunUMAP(
            object = srtIntegrated,
            reduction = "integrated.dr",
            dims = dims_use,
            reduction.name = reduction_name,
            reduction.key = paste0(reduction_name, "_"),
            n.components = n,
            verbose = FALSE,
            seed.use = seed
          )
        } else if (nr == "tsne") {
          srtIntegrated <- RunTSNE(
            object = srtIntegrated,
            reduction = "integrated.dr",
            dims = dims_use,
            reduction.name = reduction_name,
            reduction.key = paste0(reduction_name, "_"),
            dims.embed = n,
            verbose = FALSE,
            seed.use = seed
          )
        }
      }
    }
  }

  # Join layers at the end (following Seurat v5 best practice)
  # This makes the object ready for differential expression analysis
  if (verbose) cat(paste0("[", Sys.time(), "]", " Joining layers after downstream analysis...\n"))
  srtIntegrated[[assay]] <- JoinLayers(srtIntegrated[[assay]])

  return(srtIntegrated)
}

#' V4 List-based Seurat Integration (Internal)
#'
#' @inheritParams Seurat_integrate
#' @return Seurat object with integrated data
.seurat_integrate_v4_list <- function(srtList, batch, assay, do_normalization, normalization_method,
                                    do_HVF_finding, HVF_source, HVF_method, nHVF, HVF_min_intersection, HVF,
                                    do_scaling, vars_to_regress, regression_model, scale_within_batch,
                                    linear_reduction, linear_reduction_dims, linear_reduction_dims_use,
                                    linear_reduction_params, force_linear_reduction,
                                    nonlinear_reduction, nonlinear_reduction_dims, nonlinear_reduction_params,
                                    force_nonlinear_reduction, neighbor_metric, neighbor_k,
                                    cluster_algorithm, cluster_resolution,
                                    FindIntegrationAnchors_params, IntegrateData_params, IntegrateEmbeddings_params,
                                    verbose = TRUE, seed) {

  # Check inputs
  checked <- check_srtList(
    srtList = srtList, batch = batch, assay = assay,
    do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = HVF_source, HVF_method = HVF_method,
    nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
    vars_to_regress = vars_to_regress, seed = seed, check_v5 = FALSE
  )
  srtList <- checked[["srtList"]]
  HVF <- checked[["HVF"]]
  assay <- checked[["assay"]]
  type <- checked[["type"]]

  # Setup cluster algorithm
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  # Check minimum cell counts
  if (min(sapply(srtList, ncol)) < 50) {
    warning("The cell count in some batches is lower than 50, which may not be suitable for the current integration method.", immediate. = TRUE)
  }

  # Handle TFIDF normalization
  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'rlsi' integration workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
    FindIntegrationAnchors_params[["reduction"]] <- "rlsi"
    if (is.null(FindIntegrationAnchors_params[["dims"]])) {
      FindIntegrationAnchors_params[["dims"]] <- 2:min(linear_reduction_dims, 30)
    }
    if (!requireNamespace("Signac", quietly = TRUE)) {
      stop("Package 'Signac' is required for TFIDF normalization. Install it with: install.packages('Signac')")
    }
    srtMerge <- Reduce(merge, srtList)
    srtMerge <- Signac::RunTFIDF(object = srtMerge, assay = DefaultAssay(srtMerge), verbose = FALSE)
    srtMerge <- RunDimReduction(
      srtMerge,
      prefix = "", features = HVF, assay = DefaultAssay(srtMerge),
      linear_reduction = "svd", linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    srtMerge[["lsi"]] <- srtMerge[["svd"]]
    for (i in seq_along(srtList)) {
      srt <- srtList[[i]]
      cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (svd) on the data ", i, " ...\n"))
      srt <- RunDimReduction(
        srt,
        prefix = "", features = HVF, assay = DefaultAssay(srt),
        linear_reduction = "svd", linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
        verbose = FALSE, seed = seed
      )
      srt[["lsi"]] <- srt[["svd"]]
      srtList[[i]] <- srt
    }
  }

  # Handle RPCA reduction
  if (isTRUE(FindIntegrationAnchors_params[["reduction"]] == "rpca")) {
    cat(paste0("[", Sys.time(), "]", " Use 'rpca' integration workflow...\n"))
    for (i in seq_along(srtList)) {
      srt <- srtList[[i]]
      if (isTRUE(do_scaling) || (is.null(do_scaling) && !.check_scaled_v5(srt, HVF, DefaultAssay(srt)))) {
        cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data ", i, " ...\n"))
        srt <- ScaleData(
          srt,
          features = HVF,
          assay = DefaultAssay(srt),
          vars.to.regress = vars_to_regress,
          verbose = FALSE
        )
      }
      cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (pca) on the data ", i, " ...\n"))
      srt <- RunDimReduction(
        srt,
        prefix = "", features = HVF, assay = DefaultAssay(srt),
        linear_reduction = "pca", linear_reduction_dims = linear_reduction_dims,
        linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
        verbose = FALSE, seed = seed
      )
      srtList[[i]] <- srt
    }
  }

  if (is.null(names(srtList))) {
    names(srtList) <- paste0("srt_", seq_along(srtList))
  }

  # Perform integration based on normalization method
  if (normalization_method %in% c("LogNormalize", "SCT")) {
    cat(paste0("[", Sys.time(), "]", " Perform FindIntegrationAnchors on the data...\n"))
    params1 <- list(
      object.list = srtList,
      normalization.method = normalization_method,
      anchor.features = HVF,
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke(.fn = FindIntegrationAnchors, .args = params1)

    cat(paste0("[", Sys.time(), "]", " Perform integration(Seurat) on the data...\n"))
    params2 <- list(
      anchorset = srt_anchors,
      new.assay.name = "Seuratcorrected",
      normalization.method = normalization_method,
      features.to.integrate = HVF,
      verbose = FALSE
    )
    for (nm in names(IntegrateData_params)) {
      params2[[nm]] <- IntegrateData_params[[nm]]
    }
    srtIntegrated <- invoke(.fn = IntegrateData, .args = params2)

    DefaultAssay(srtIntegrated) <- "Seuratcorrected"
    VariableFeatures(srtIntegrated[["Seuratcorrected"]]) <- HVF

    if (isTRUE(do_scaling) || (is.null(do_scaling) && !.check_scaled_v5(srtIntegrated, HVF, DefaultAssay(srtIntegrated)))) {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
      srtIntegrated <- ScaleData(
        srtIntegrated,
        features = HVF,
        assay = DefaultAssay(srtIntegrated),
        vars.to.regress = vars_to_regress,
        verbose = FALSE
      )
    }

    cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
    srtIntegrated <- RunDimReduction(
      srtIntegrated,
      prefix = "Seurat", features = HVF, assay = DefaultAssay(srtIntegrated),
      linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims,
      linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- srtIntegrated@reductions[[paste0("Seurat", linear_reduction)]]@misc[["dims_estimate"]] %||% 1:linear_reduction_dims
    }
  } else if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " Perform FindIntegrationAnchors on the data...\n"))
    params1 <- list(
      object.list = srtList,
      normalization.method = "LogNormalize",
      anchor.features = HVF,
      reduction = "rlsi",
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    srt_anchors <- invoke(.fn = FindIntegrationAnchors, .args = params1)

    cat(paste0("[", Sys.time(), "]", " Perform integration(Seurat) on the data...\n"))
    params2 <- list(
      anchorset = srt_anchors,
      reductions = srtMerge[["lsi"]],
      new.reduction.name = "Seuratlsi",
      verbose = FALSE
    )
    for (nm in names(IntegrateEmbeddings_params)) {
      params2[[nm]] <- IntegrateEmbeddings_params[[nm]]
    }
    srtIntegrated <- invoke(.fn = IntegrateEmbeddings, .args = params2)

    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- 2:max(srtIntegrated@reductions[[paste0("Seurat", linear_reduction)]]@misc[["dims_estimate"]]) %||% 2:linear_reduction_dims
    }
    linear_reduction <- "lsi"
  }

  # Clustering and reductions
  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(
        object = srtIntegrated, reduction = paste0("Seurat", linear_reduction), dims = linear_reduction_dims_use,
        annoy.metric = neighbor_metric, k.param = neighbor_k,
        graph.name = paste0("Seurat_", c("KNN", "SNN")), verbose = FALSE
      )

      cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
      srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution,
                                   algorithm = cluster_algorithm_index, method = "igraph",
                                   graph.name = "Seurat_SNN", verbose = FALSE)
      cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", layer = "data")
      srtIntegrated[["seurat_clusters"]] <- NULL
      srtIntegrated[["Seuratclusters"]] <- Idents(srtIntegrated)
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srtIntegrated,
            prefix = "Seurat",
            reduction_use = paste0("Seurat", linear_reduction), reduction_dims = linear_reduction_dims_use,
            graph_use = "Seurat_SNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Seurat_HVF"]] <- HVF

  return(srtIntegrated)
}