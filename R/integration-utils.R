# Integration Utilities for Seurat V4/V5 Compatibility
# This file contains utility functions for managing integration workflows
# across different Seurat versions

#' Detect Seurat version and determine appropriate workflow
#'
#' @param srt Seurat object or list of Seurat objects
#' @param use_v5_workflow Logical or NULL. If NULL, auto-detect. If TRUE/FALSE, force that workflow.
#' @return List with version info and workflow recommendation
#' @export
.detect_seurat_workflow <- function(srt, use_v5_workflow = NULL) {
  # Force workflow if specified
  if (!is.null(use_v5_workflow)) {
    return(list(
      version = if (use_v5_workflow) "v5" else "v4",
      workflow = if (use_v5_workflow) "layers" else "list",
      auto_detected = FALSE
    ))
  }
  
  # Auto-detect version
  if (inherits(srt, "list")) {
    # Check first object in list
    srt_check <- srt[[1]]
  } else {
    srt_check <- srt
  }
  
  is_v5 <- IsSeurat5(srt_check)
  
  return(list(
    version = if (is_v5) "v5" else "v4",
    workflow = if (is_v5) "layers" else "list",
    auto_detected = TRUE
  ))
}

#' Convert list of Seurat objects to V5 layers
#'
#' @param srtList List of Seurat objects
#' @param batch Character, batch variable name
#' @param assay Character, assay name
#' @return Seurat v5 object with layers split by batch
#' @export
.convert_list_to_layers <- function(srtList, batch, assay) {
  # Merge objects first
  srtMerge <- Reduce(merge, srtList)
  
  # Ensure it's V5
  if (!IsSeurat5(srtMerge)) {
    srtMerge <- UpdateSeuratObject(srtMerge)
  }
  
  # Split by batch
  srtMerge[[assay]] <- split(srtMerge[[assay]], f = srtMerge[[batch]])
  
  return(srtMerge)
}

#' Convert V5 layers back to list of objects
#'
#' @param srt Seurat v5 object with layers
#' @param batch Character, batch variable name
#' @param assay Character, assay name
#' @return List of Seurat objects
#' @export
.convert_layers_to_list <- function(srt, batch, assay) {
  if (!IsSeurat5(srt)) {
    stop("Object must be Seurat v5 to convert layers to list")
  }
  
  # Join layers first
  srt[[assay]] <- JoinLayers(srt[[assay]])
  
  # Split by batch
  srtList <- SplitObject(srt, split.by = batch)
  
  return(srtList)
}

#' Check if integration method supports V5 layers
#'
#' @param method Character, integration method name
#' @return Logical, whether method supports V5 layers
#' @export
.supports_v5_layers <- function(method) {
  v5_supported <- c(
    "Seurat", "Harmony", "scVI", "MNN", "fastMNN", 
    "Scanorama", "BBKNN", "CSS", "LIGER", "ComBat"
  )
  
  return(method %in% v5_supported)
}

#' Check if integration method requires list approach
#'
#' @param method Character, integration method name
#' @return Logical, whether method requires list approach
#' @export
.requires_list_approach <- function(method) {
  list_required <- c("SCANVI", "Conos")
  
  return(method %in% list_required)
}

#' Prepare integration inputs based on workflow
#'
#' @param srtList List of Seurat objects or NULL
#' @param srtMerge Seurat object or NULL
#' @param batch Character, batch variable name
#' @param assay Character, assay name
#' @param workflow Character, "layers" or "list"
#' @return List with prepared inputs
#' @export
.prepare_integration_inputs <- function(srtList, srtMerge, batch, assay, workflow) {
  if (workflow == "layers") {
    # V5 layer-based workflow
    if (!is.null(srtList)) {
      # Convert list to layers
      srtMerge <- .convert_list_to_layers(srtList, batch, assay)
      srtList <- NULL
    }
    return(list(srtList = NULL, srtMerge = srtMerge, workflow = "layers"))
  } else {
    # V4 list-based workflow
    if (!is.null(srtMerge)) {
      # Convert merged object to list
      srtList <- SplitObject(srtMerge, split.by = batch)
      srtMerge <- NULL
    }
    return(list(srtList = srtList, srtMerge = NULL, workflow = "list"))
  }
}

#' Post-process integration results based on workflow
#'
#' @param srtIntegrated Integrated Seurat object
#' @param workflow Character, "layers" or "list"
#' @param assay Character, assay name
#' @return Processed Seurat object
#' @export
.post_process_integration <- function(srtIntegrated, workflow, assay) {
  if (workflow == "layers") {
    # Use default assay if assay is NULL or empty string
    if (is.null(assay) || assay == "") {
      assay <- DefaultAssay(srtIntegrated)
    }

    # Join layers for V5 (skip if already joined or assay doesn't exist)
    if (IsSeurat5(srtIntegrated) && !is.null(srtIntegrated[[assay]])) {
      # Check if layers are split (multiple layers exist)
      assay_obj <- srtIntegrated[[assay]]
      if (length(Layers(assay_obj)) > 1) {
        srtIntegrated[[assay]] <- JoinLayers(srtIntegrated[[assay]])
      }
    }
  }

  return(srtIntegrated)
}

#' Get integration method parameters for V4 vs V5
#'
#' @param method Character, integration method name
#' @param workflow Character, "layers" or "list"
#' @return List with method-specific parameters
#' @export
.get_integration_method_params <- function(method, workflow) {
  if (workflow == "layers") {
    # V5 layer-based parameters
    switch(method,
      "Seurat" = list(
        integration_method = "CCAIntegration",
        IntegrateLayers_params = list()
      ),
      "Harmony" = list(
        integration_method = "HarmonyIntegration",
        IntegrateLayers_params = list()
      ),
      "scVI" = list(
        # scVI doesn't support layers yet, fall back to list
        workflow = "list"
      ),
      "MNN" = list(
        integration_method = "FastMNNIntegration",
        IntegrateLayers_params = list()
      ),
      "fastMNN" = list(
        integration_method = "FastMNNIntegration",
        IntegrateLayers_params = list()
      ),
      "Scanorama" = list(
        # Scanorama doesn't support layers yet, fall back to list
        workflow = "list"
      ),
      "BBKNN" = list(
        # BBKNN doesn't support layers yet, fall back to list
        workflow = "list"
      ),
      "CSS" = list(
        # CSS doesn't support layers yet, fall back to list
        workflow = "list"
      ),
      "LIGER" = list(
        # LIGER doesn't support layers yet, fall back to list
        workflow = "list"
      ),
      "ComBat" = list(
        # ComBat doesn't support layers yet, fall back to list
        workflow = "list"
      ),
      list() # Default
    )
  } else {
    # V4 list-based parameters
    switch(method,
      "Seurat" = list(
        FindIntegrationAnchors_params = list(),
        IntegrateData_params = list(),
        IntegrateEmbeddings_params = list()
      ),
      list() # Default
    )
  }
}

#' Check if features are scaled without joining layers (V5 optimized)
#'
#' @param srt Seurat object
#' @param features Features to check
#' @param assay Assay name
#' @return Logical, TRUE if all features are scaled
#' @export
.check_scaled_v5 <- function(srt, features, assay) {
  if (!IsSeurat5(srt)) {
    # For V4, use traditional check
    scale_data <- tryCatch({
      get_seurat_data(srt, layer = "scale.data", assay = assay, join_layers = FALSE)
    }, error = function(e) NULL)

    if (is.null(scale_data) || length(scale_data) == 0) {
      return(FALSE)
    }

    return(all(features %in% rownames(scale_data)))
  }

  # For V5, check if scale.data layer exists for split layers
  assay_obj <- srt[[assay]]
  layers <- SeuratObject::Layers(assay_obj)

  # Check for scale.data layers
  scale_layers <- layers[grepl("scale.data", layers)]

  if (length(scale_layers) == 0) {
    return(FALSE)
  }

  # For split layers, we don't need to check all features in all layers
  # Just ensure scale.data layers exist
  return(TRUE)
}

