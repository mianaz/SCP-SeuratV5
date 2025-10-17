# SCPNext Refactoring Plan

> **Generated:** 2025-10-17
> **Target:** Reduce redundancy by ~6,500-7,000 lines (18-20%)
> **Goal:** Lean, efficient package leveraging R's built-in capabilities

## Executive Summary

This refactoring plan addresses significant code duplication in the SCPNext package, following the CLAUDE.md principle: "The goal should be a lean, efficient package that leverages R's built-in capabilities rather than reimplementing them."

### Key Metrics
- **Current Size:** 36,188 lines of code
- **Estimated Redundancy:** 6,500-7,000 lines (18-20%)
- **Target Size:** ~29,000-30,000 lines
- **Maintainability Improvement:** 40-50% (fewer locations to update)

---

## Phase 1: Critical Consolidation (Highest ROI)

### Priority 1A: Integration Methods Consolidation

**Current State:** 15 integration functions with 90%+ identical code
**Files:** SCP-workflow.R (lines 100-2900)
**Redundancy:** ~2,500 lines → Target: ~400 lines
**Impact:** 2,100 lines saved (33% of total redundancy)

#### Current Pattern (Repeated 15x)

```r
# Current: Each integration method is 150-200 lines of duplicated code
Harmony_integrate <- function(srtList, batch, append = FALSE, srtMerge = NULL,
                             assay = NULL, do_normalization = FALSE, ...) {
  # 1. Validate inputs (20-30 lines, duplicated)
  check <- check_srtList(
    srtList = srtList,
    batch = batch,
    assay = assay,
    do_normalization = do_normalization,
    ...
  )
  # Extract from check
  srtList <- check[["srtList"]]
  assay <- check[["assay"]]

  # 2. Merge objects (15-20 lines, duplicated)
  if (is.null(srtMerge)) {
    srtMerge <- merge(x = srtList[[1]], y = srtList[-1])
    srtMerge[[batch]] <- ...
  }

  # 3. Prepare features (30-40 lines, duplicated)
  VariableFeatures(srtMerge, assay = assay) <- ...
  DefaultAssay(srtMerge) <- assay

  # 4. Call integration method (5-10 lines, ONLY UNIQUE PART)
  srtMerge <- RunHarmony2(
    object = srtMerge,
    group.by.vars = batch,
    assay.use = assay,
    ...
  )

  # 5. Append results (30-40 lines, duplicated)
  if (isTRUE(append)) {
    srtMerge <- Integration_SCP_append(
      srtMerge = srtMerge,
      srtList = srtList,
      batch = batch,
      ...
    )
  }

  # 6. Return (duplicated)
  return(srtMerge)
}

# This pattern is repeated for:
# Seurat_integrate, scVI_integrate, MNN_integrate, fastMNN_integrate,
# Scanorama_integrate, BBKNN_integrate, CSS_integrate, LIGER_integrate,
# ComBat_integrate, Conos_integrate, SCANVI_integrate, Uncorrected_integrate
```

#### Proposed Solution: Factory Pattern with Method Dispatch

```r
#' Batch Correction and Integration
#'
#' Unified interface for all integration methods with method-specific dispatch
#'
#' @param srtList List of Seurat objects or merged Seurat object
#' @param batch Character, batch variable name
#' @param method Character, integration method. Options:
#'   - "Seurat": Seurat CCA integration
#'   - "Harmony": Harmony batch correction
#'   - "scVI": scVI deep learning integration
#'   - "MNN": Mutual nearest neighbors
#'   - "fastMNN": Fast MNN
#'   - "Scanorama": Scanorama integration
#'   - "BBKNN": BBKNN k-NN graph
#'   - "CSS": Cluster similarity spectrum
#'   - "LIGER": iNMF factorization
#'   - "ComBat": ComBat batch correction
#'   - "Conos": Conos integration
#'   - "SCANVI": Semi-supervised scVI
#'   - "Uncorrected": No correction (baseline)
#' @param ... Additional arguments passed to method-specific function
#'
#' @export
BatchCorrect <- function(srtList, batch, method = "Harmony",
                        append = FALSE, srtMerge = NULL,
                        assay = NULL, do_normalization = FALSE, ...) {

  # 1. Validate and prepare (single implementation)
  check <- check_srtList(
    srtList = srtList,
    batch = batch,
    assay = assay,
    do_normalization = do_normalization,
    ...
  )
  srtList <- check[["srtList"]]
  assay <- check[["assay"]]

  # 2. Merge if needed (single implementation)
  if (is.null(srtMerge)) {
    srtMerge <- .merge_srtList(srtList, batch, assay)
  }

  # 3. Prepare features (single implementation)
  srtMerge <- .prepare_integration_features(srtMerge, assay, check)

  # 4. Dispatch to method-specific handler (lookup table)
  integration_handler <- .get_integration_handler(method)
  srtMerge <- integration_handler(
    srtMerge = srtMerge,
    srtList = srtList,
    batch = batch,
    assay = assay,
    ...
  )

  # 5. Append if requested (single implementation)
  if (isTRUE(append)) {
    srtMerge <- .append_integration_results(
      srtMerge = srtMerge,
      srtList = srtList,
      batch = batch,
      method = method,
      ...
    )
  }

  return(srtMerge)
}

# Internal: Method dispatch lookup
.get_integration_handler <- function(method) {
  handlers <- list(
    "Harmony" = .integrate_harmony,
    "Seurat" = .integrate_seurat,
    "scVI" = .integrate_scvi,
    "MNN" = .integrate_mnn,
    "fastMNN" = .integrate_fastmnn,
    "Scanorama" = .integrate_scanorama,
    "BBKNN" = .integrate_bbknn,
    "CSS" = .integrate_css,
    "LIGER" = .integrate_liger,
    "ComBat" = .integrate_combat,
    "Conos" = .integrate_conos,
    "SCANVI" = .integrate_scanvi,
    "Uncorrected" = .integrate_uncorrected
  )

  method <- match.arg(method, names(handlers))
  return(handlers[[method]])
}

# Internal: Harmony-specific handler (only unique logic)
.integrate_harmony <- function(srtMerge, batch, assay, ...) {
  require_packages("harmony")
  srtMerge <- RunHarmony2(
    object = srtMerge,
    group.by.vars = batch,
    assay.use = assay,
    ...
  )
  return(srtMerge)
}

# Internal: scVI-specific handler (only unique logic)
.integrate_scvi <- function(srtMerge, srtList, batch, assay, ...) {
  use_uv_env()
  require_packages("reticulate")

  # Convert to AnnData
  adata <- srt_to_adata(srtMerge, assay = assay)

  # Run scVI
  scvi <- import("scvi")
  scvi$model$SCVI$setup_anndata(adata, batch_key = batch)
  model <- scvi$model$SCVI(adata, ...)
  model$train()

  # Get latent representation
  latent <- model$get_latent_representation()

  # Add to Seurat
  srtMerge[["scvi"]] <- CreateDimReducObject(
    embeddings = latent,
    key = "scvi_",
    assay = assay
  )

  return(srtMerge)
}

# ... (similar for other 11 methods, ~20-30 lines each)

# Backward compatibility: Keep old function names as thin wrappers
#' @rdname BatchCorrect
#' @export
Harmony_integrate <- function(...) {
  BatchCorrect(method = "Harmony", ...)
}

#' @rdname BatchCorrect
#' @export
scVI_integrate <- function(...) {
  BatchCorrect(method = "scVI", ...)
}

# ... (similar for other 13 methods)
```

**Benefits:**
- Single point of maintenance for common logic
- Consistent error handling and validation
- Easy to add new integration methods
- Backward compatible with existing code
- Reduces 2,500 lines to ~400 lines (84% reduction)

---

### Priority 1B: Dimension Reduction Consolidation

**Current State:** 11 Run* functions × 3-4 S3 methods = 33+ functions
**Files:** Seurat-function.R (1,871 lines)
**Redundancy:** ~800 lines → Target: ~200 lines
**Impact:** 600 lines saved (9% of total redundancy)

#### Current Pattern (Repeated 11x with 3-4 methods each)

```r
# Current: Generic + Seurat + Assay + default (repeated 11 times)

# Generic
RunNMF <- function(object, ...) {
  UseMethod("RunNMF")
}

# Seurat method
RunNMF.Seurat <- function(object, assay = NULL, slot = "data", ...) {
  assay <- assay %||% DefaultAssay(object)
  data <- GetAssayData(object, assay = assay, slot = slot)
  embeddings <- RunNMF.default(data, ...)

  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = embeddings,
    key = reduction.key,
    assay = assay
  )
  return(object)
}

# Default method
RunNMF.default <- function(object, nfeatures = 50, ...) {
  require_packages("RcppML")
  result <- RcppML::nmf(object, k = nfeatures, ...)
  return(result@w)
}

# This entire pattern repeats for:
# RunMDS, RunGLMPCA, RunDM, RunUMAP2, RunPaCMAP, RunPHATE,
# RunTriMap, RunLargeVis, RunFR, RunHarmony2
```

#### Proposed Solution: Unified Factory Function

```r
#' Dimensionality Reduction
#'
#' Unified interface for all dimensionality reduction methods
#'
#' @param object Seurat object or matrix
#' @param method Character, reduction method. Options:
#'   - "nmf": Non-negative matrix factorization
#'   - "mds": Multidimensional scaling
#'   - "glmpca": Generalized PCA for count data
#'   - "dm": Diffusion map
#'   - "umap": UMAP (enhanced)
#'   - "pacmap": PaCMAP
#'   - "phate": PHATE
#'   - "trimap": TriMap
#'   - "largevis": LargeVis
#'   - "fr": Force-directed graph
#' @param ... Method-specific parameters
#'
#' @return Seurat object with reduction or matrix of embeddings
#' @export
RunDimReduction <- function(object, method = "umap", ...) {
  UseMethod("RunDimReduction")
}

#' @export
RunDimReduction.Seurat <- function(object, method = "umap",
                                   assay = NULL, slot = "data",
                                   reduction.name = NULL,
                                   reduction.key = NULL,
                                   n_components = 50, ...) {

  # Validate method
  method <- tolower(method)
  valid_methods <- c("nmf", "mds", "glmpca", "dm", "umap",
                     "pacmap", "phate", "trimap", "largevis", "fr")
  method <- match.arg(method, valid_methods)

  # Get data
  assay <- assay %||% DefaultAssay(object)
  data <- GetAssayData(object, assay = assay, slot = slot)

  # Set reduction names
  reduction.name <- reduction.name %||% method
  reduction.key <- reduction.key %||% paste0(method, "_")

  # Dispatch to method
  embeddings <- .run_reduction_method(
    data = data,
    method = method,
    n_components = n_components,
    ...
  )

  # Add to object
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = embeddings,
    key = reduction.key,
    assay = assay
  )

  # Log command
  object <- LogSeuratCommand(object)

  return(object)
}

#' @export
RunDimReduction.default <- function(object, method = "umap",
                                    n_components = 50, ...) {
  method <- match.arg(tolower(method),
                     c("nmf", "mds", "glmpca", "dm", "umap",
                       "pacmap", "phate", "trimap", "largevis", "fr"))

  embeddings <- .run_reduction_method(
    data = object,
    method = method,
    n_components = n_components,
    ...
  )

  return(embeddings)
}

# Internal: Method dispatch
.run_reduction_method <- function(data, method, n_components, ...) {
  switch(method,
    "nmf" = .reduce_nmf(data, n_components, ...),
    "mds" = .reduce_mds(data, n_components, ...),
    "glmpca" = .reduce_glmpca(data, n_components, ...),
    "dm" = .reduce_dm(data, n_components, ...),
    "umap" = .reduce_umap(data, n_components, ...),
    "pacmap" = .reduce_pacmap(data, n_components, ...),
    "phate" = .reduce_phate(data, n_components, ...),
    "trimap" = .reduce_trimap(data, n_components, ...),
    "largevis" = .reduce_largevis(data, n_components, ...),
    "fr" = .reduce_fr(data, n_components, ...)
  )
}

# Internal: Method-specific implementations (only unique logic)
.reduce_nmf <- function(data, n_components, ...) {
  require_packages("RcppML")
  result <- RcppML::nmf(data, k = n_components, ...)
  return(result@w)
}

.reduce_umap <- function(data, n_components, ...) {
  require_packages("uwot")
  result <- uwot::umap(
    X = t(data),
    n_components = n_components,
    ...
  )
  rownames(result) <- colnames(data)
  return(result)
}

# ... (similar for other 8 methods, ~15-20 lines each)

# Backward compatibility: Keep old function names
#' @rdname RunDimReduction
#' @export
RunNMF <- function(object, ...) {
  RunDimReduction(object, method = "nmf", ...)
}

#' @rdname RunDimReduction
#' @export
RunUMAP2 <- function(object, ...) {
  RunDimReduction(object, method = "umap", ...)
}

# ... (similar for other 9 methods)
```

**Benefits:**
- Eliminates S3 method boilerplate duplication
- Consistent interface across all reduction methods
- Single point for parameter validation
- Easy to add new reduction methods
- Reduces 800 lines to ~200 lines (75% reduction)
- Backward compatible

---

### Priority 1C: Plotting Engine Consolidation

**Current State:** 60+ plotting functions with massive duplication
**Files:** SCP-plot.R (14,990 lines)
**Redundancy:** ~3,000 lines of duplicated logic
**Target:** Extract shared rendering engine (~500 lines), keep thin wrappers
**Impact:** 2,500 lines saved (38% of total redundancy)

#### Current Pattern Issues

```r
# Duplicated across CellDimPlot, FeatureDimPlot, CellDimPlot3D, FeatureDimPlot3D:

# 1. Color scaling logic (50-80 lines, repeated 4x)
if (is.numeric(data$color)) {
  p <- p + scale_color_gradientn(...)
} else {
  p <- p + scale_color_manual(...)
}

# 2. Faceting logic (40-60 lines, repeated 10x)
if (!is.null(split.by)) {
  if (is.numeric(split.by) || length(unique(...)) > facet_threshold) {
    p <- p + facet_wrap(...)
  } else {
    p <- p + facet_grid(...)
  }
}

# 3. Legend customization (30-40 lines, repeated 15x)
if (isFALSE(legend)) {
  p <- p + theme(legend.position = "none")
} else if (is.character(legend)) {
  p <- p + theme(legend.position = legend)
}

# 4. Label/mark additions (100+ lines, repeated 4x)
if (isTRUE(add_mark)) {
  # Complex mark addition logic
}

# 5. Density overlay (80+ lines, repeated 3x)
if (isTRUE(add_density)) {
  # Density calculation and overlay
}

# 6. Trajectory overlay (150+ lines, repeated 4x)
if (!is.null(lineages)) {
  # Trajectory line rendering
}

# 7. Velocity overlay (120+ lines, repeated 3x)
if (isTRUE(velocity)) {
  # Velocity arrow rendering
}
```

#### Proposed Solution: Centralized Plot Engine

```r
# File: R/plot-engine.R (new file)

#' Internal: Core plotting engine
#'
#' Handles all common plotting operations with consistent interface
.plot_engine <- function(data, plot_type, aesthetic_mapping, layers = list(), ...) {

  # 1. Create base plot
  p <- ggplot(data, aesthetic_mapping)

  # 2. Add base layer (scatter, hex, etc.)
  p <- .add_base_layer(p, plot_type, ...)

  # 3. Add optional layers in order
  for (layer in layers) {
    p <- .add_layer(p, layer, data, ...)
  }

  # 4. Apply scales
  p <- .apply_scales(p, data, ...)

  # 5. Apply faceting
  p <- .apply_facets(p, ...)

  # 6. Apply theme
  p <- .apply_theme(p, ...)

  # 7. Customize legend
  p <- .customize_legend(p, ...)

  return(p)
}

# Centralized color scaling (used by all plots)
.apply_scales <- function(p, data, color.by, palette = NULL,
                         colors = NULL, na_color = "grey80", ...) {

  if (is.numeric(data[[color.by]])) {
    # Continuous scale
    if (is.null(colors)) {
      colors <- palette_default("gradient")
    }
    p <- p + scale_color_gradientn(
      colors = colors,
      na.value = na_color,
      ...
    )
  } else {
    # Discrete scale
    if (is.null(colors)) {
      colors <- palette_dimplot(n = length(unique(data[[color.by]])))
    }
    p <- p + scale_color_manual(
      values = colors,
      na.value = na_color,
      ...
    )
  }

  return(p)
}

# Centralized faceting (used by all plots)
.apply_facets <- function(p, split.by = NULL, facet_threshold = 12,
                         nrow = NULL, ncol = NULL, ...) {

  if (is.null(split.by)) {
    return(p)
  }

  n_facets <- length(unique(p$data[[split.by]]))

  if (n_facets > facet_threshold) {
    p <- p + facet_wrap(
      as.formula(paste("~", split.by)),
      nrow = nrow,
      ncol = ncol,
      ...
    )
  } else {
    p <- p + facet_grid(
      as.formula(paste("~", split.by)),
      ...
    )
  }

  return(p)
}

# Centralized layer additions
.add_layer <- function(p, layer, data, ...) {
  switch(layer$type,
    "density" = .add_density_layer(p, layer, data, ...),
    "mark" = .add_mark_layer(p, layer, ...),
    "trajectory" = .add_trajectory_layer(p, layer, ...),
    "velocity" = .add_velocity_layer(p, layer, ...),
    "paga" = .add_paga_layer(p, layer, ...),
    "graph" = .add_graph_layer(p, layer, ...),
    p  # Return unchanged if unknown layer
  )
}

# ... (Other centralized functions)
```

#### Updated User-Facing Functions (Thin Wrappers)

```r
#' Cell Dimension Plot
#'
#' @export
CellDimPlot <- function(srt, group.by = NULL, reduction = NULL,
                        dims = c(1, 2), split.by = NULL,
                        add_density = FALSE, add_mark = FALSE,
                        add_trajectory = FALSE, lineages = NULL,
                        velocity = FALSE, paga = NULL, ...) {

  # 1. Prepare data
  plot_data <- .prepare_dimplot_data(
    srt = srt,
    group.by = group.by,
    reduction = reduction,
    dims = dims,
    split.by = split.by
  )

  # 2. Define aesthetic mapping
  aes_map <- aes(x = .data[[paste0("dim", dims[1])]],
                 y = .data[[paste0("dim", dims[2])]],
                 color = .data[[group.by]])

  # 3. Define layers to add
  layers <- list()
  if (add_density) layers <- c(layers, list(list(type = "density")))
  if (add_mark) layers <- c(layers, list(list(type = "mark")))
  if (add_trajectory) layers <- c(layers, list(list(type = "trajectory", lineages = lineages)))
  if (velocity) layers <- c(layers, list(list(type = "velocity")))
  if (!is.null(paga)) layers <- c(layers, list(list(type = "paga", paga = paga)))

  # 4. Call plot engine
  p <- .plot_engine(
    data = plot_data,
    plot_type = "scatter",
    aesthetic_mapping = aes_map,
    layers = layers,
    color.by = group.by,
    split.by = split.by,
    ...
  )

  return(p)
}

#' Feature Dimension Plot
#'
#' @export
FeatureDimPlot <- function(srt, features, reduction = NULL,
                          dims = c(1, 2), split.by = NULL,
                          add_density = FALSE, add_trajectory = FALSE,
                          lineages = NULL, velocity = FALSE, ...) {

  # Very similar to CellDimPlot, but different data preparation
  plot_data <- .prepare_feature_dimplot_data(
    srt = srt,
    features = features,
    reduction = reduction,
    dims = dims,
    split.by = split.by
  )

  aes_map <- aes(x = .data[[paste0("dim", dims[1])]],
                 y = .data[[paste0("dim", dims[2])]],
                 color = .data[["expression"]])

  layers <- list()
  if (add_density) layers <- c(layers, list(list(type = "density")))
  if (add_trajectory) layers <- c(layers, list(list(type = "trajectory", lineages = lineages)))
  if (velocity) layers <- c(layers, list(list(type = "velocity")))

  p <- .plot_engine(
    data = plot_data,
    plot_type = "scatter",
    aesthetic_mapping = aes_map,
    layers = layers,
    color.by = "expression",
    split.by = split.by,
    ...
  )

  return(p)
}

# CellDimPlot3D, FeatureDimPlot3D, StatPlot, etc. follow similar pattern
```

**Benefits:**
- Single point of maintenance for all plotting logic
- Consistent behavior across all plot types
- Easy to add new plot features globally
- Dramatically reduces SCP-plot.R size
- Users see no change in API
- Reduces duplication by ~2,500 lines

---

## Phase 2: Medium Priority Consolidation

### Priority 2A: Reference Mapping Functions

**Current State:** 5 mapping functions with shared structure
**Files:** SCP-projection.R (775 lines)
**Redundancy:** ~400 lines → Target: ~200 lines
**Impact:** 200 lines saved

#### Proposed Solution

```r
#' Reference Mapping
#'
#' Unified interface for reference-based cell projection
#'
#' @param srt Query Seurat object
#' @param srt_ref Reference Seurat object
#' @param method Mapping method: "knn", "pca", "seurat", "css", "symphony"
#' @param ... Method-specific parameters
#'
#' @export
ProjectCells <- function(srt, srt_ref, method = "knn", ...) {

  method <- match.arg(tolower(method),
                     c("knn", "pca", "seurat", "css", "symphony"))

  # Common validation
  .validate_projection_inputs(srt, srt_ref)

  # Dispatch to method
  result <- switch(method,
    "knn" = .project_knn(srt, srt_ref, ...),
    "pca" = .project_pca(srt, srt_ref, ...),
    "seurat" = .project_seurat(srt, srt_ref, ...),
    "css" = .project_css(srt, srt_ref, ...),
    "symphony" = .project_symphony(srt, srt_ref, ...)
  )

  return(result)
}

# Backward compatibility
#' @rdname ProjectCells
#' @export
RunKNNMap <- function(...) ProjectCells(method = "knn", ...)

#' @rdname ProjectCells
#' @export
RunPCAMap <- function(...) ProjectCells(method = "pca", ...)

# ... (similar for others)
```

---

### Priority 2B: Doublet Detection Consolidation

**Current State:** 4 doublet detection methods with duplicated validation
**Files:** SCP-cellqc.R (555 lines)
**Redundancy:** ~200 lines → Target: ~100 lines
**Impact:** 100 lines saved

#### Proposed Solution

```r
# Create shared validation wrapper
.run_doublet_method <- function(srt, method, ...) {
  # Common validation (extracted once)
  srt <- .validate_for_doublet_calling(srt)

  # Dispatch to method
  result <- switch(method,
    "scDblFinder" = .doublet_scdblfinder(srt, ...),
    "scds" = .doublet_scds(srt, ...),
    "Scrublet" = .doublet_scrublet(srt, ...),
    "DoubletDetection" = .doublet_detection(srt, ...)
  )

  # Common result formatting (extracted once)
  srt <- .format_doublet_results(srt, result, method)

  return(srt)
}

# Simplify db_* functions to method-specific logic only
db_scDblFinder <- function(srt, ...) {
  .run_doublet_method(srt, "scDblFinder", ...)
}
```

---

### Priority 2C: P-value Combination Methods

**Current State:** 6 p-value combination functions
**Files:** SCP-analysis.R (lines 1700-2000)
**Redundancy:** ~150 lines → Target: ~30 lines
**Impact:** 120 lines saved

#### Proposed Solution

```r
#' Combine P-values
#'
#' @param p Vector of p-values
#' @param method Combination method: "fisher", "stouffer", "wilkinson",
#'               "minimum", "maximum", "mean", "sum", "vote"
#' @param weights Optional weights for each p-value
#'
#' @export
CombinePvalues <- function(p, method = "fisher", weights = NULL) {

  method <- match.arg(method,
    c("fisher", "stouffer", "wilkinson", "minimum", "maximum",
      "mean", "sum", "vote"))

  # Remove NAs
  if (any(is.na(p))) {
    warning("Removing ", sum(is.na(p)), " NA p-values")
    p <- p[!is.na(p)]
  }

  # Dispatch
  switch(method,
    "fisher" = .combine_fisher(p, weights),
    "stouffer" = .combine_stouffer(p, weights),
    "wilkinson" = .combine_wilkinson(p, weights),
    "minimum" = min(p) * length(p),  # Bonferroni correction
    "maximum" = max(p),
    "mean" = mean(p),
    "sum" = sum(p) / length(p),
    "vote" = sum(p < 0.05) / length(p)
  )
}

# Keep old function names for compatibility
wilkinsonp <- function(p, ...) CombinePvalues(p, "wilkinson", ...)
minimump <- function(p, ...) CombinePvalues(p, "minimum", ...)
# ... etc
```

---

## Phase 3: File Organization & Modularization

### Reorganize utils.R (1,442 lines → 4 focused files)

**Current:** Mixed concerns in single file
**Target:** Logical separation

```
Current: R/utils.R (1,442 lines)
├── Python/UV management (500 lines)
├── Data access functions (400 lines)
├── Color/palette functions (200 lines)
├── Misc utilities (200 lines)
└── Plotting helpers (142 lines)

Proposed:
R/utils-python.R (500 lines) - Python/UV management
R/utils-data.R (400 lines) - Seurat V4/V5 data access
R/utils-colors.R (200 lines) - Color and palette functions
R/utils-misc.R (200 lines) - General utilities
R/plot-engine.R (500 lines) - NEW: Centralized plotting core
```

---

## Phase 4: Implementation Strategy

### Step 1: Create New Infrastructure (Non-Breaking)

1. Create `R/plot-engine.R` with centralized plotting functions
2. Create `R/integration-engine.R` with `BatchCorrect()` function
3. Create `R/reduction-engine.R` with unified `RunDimReduction()`
4. Add comprehensive tests for new functions

### Step 2: Migrate Existing Functions (Backward Compatible)

1. Update existing functions to call new centralized functions
2. Keep old function names as thin wrappers
3. Test thoroughly to ensure no behavior changes
4. Update documentation with deprecation warnings

### Step 3: Optimize and Clean

1. Remove redundant internal functions
2. Consolidate validation logic
3. Update examples and vignettes
4. Performance profiling and optimization

### Step 4: Documentation Update

1. Update all function documentation
2. Create migration guide for users
3. Add new examples showcasing unified interface
4. Update README with new patterns

---

## Migration Guide for Users

### For Integration Methods

```r
# Old way (still works)
srt <- Harmony_integrate(srtList, batch = "batch")

# New way (recommended)
srt <- BatchCorrect(srtList, batch = "batch", method = "Harmony")

# Or use the unified dispatcher
srt <- Integration_SCP(srtList, batch = "batch", method = "Harmony")
```

### For Dimension Reduction

```r
# Old way (still works)
srt <- RunUMAP2(srt, n_components = 30)

# New way (recommended)
srt <- RunDimReduction(srt, method = "umap", n_components = 30)
```

### For Reference Mapping

```r
# Old way (still works)
srt <- RunKNNMap(srt, srt_ref)

# New way (recommended)
srt <- ProjectCells(srt, srt_ref, method = "knn")
```

---

## Testing Strategy

### Unit Tests

```r
# Test each integration method produces same results
test_that("BatchCorrect produces same results as old functions", {
  # Load test data
  data("pancreas_sub")
  srtList <- SplitObject(pancreas_sub, split.by = "tech")

  # Run old way
  srt_old <- Harmony_integrate(srtList, batch = "tech")

  # Run new way
  srt_new <- BatchCorrect(srtList, batch = "tech", method = "Harmony")

  # Compare embeddings
  expect_equal(
    Embeddings(srt_old, "harmony"),
    Embeddings(srt_new, "harmony"),
    tolerance = 1e-6
  )
})
```

### Integration Tests

```r
# Test full pipeline still works
test_that("Standard workflow completes successfully", {
  srt <- Standard_SCP(pancreas_sub)
  srtList <- SplitObject(srt, split.by = "tech")
  srt <- BatchCorrect(srtList, batch = "tech", method = "Harmony")
  srt <- RunDimReduction(srt, method = "umap")
  p <- CellDimPlot(srt, group.by = "celltype")

  expect_s3_class(srt, "Seurat")
  expect_s3_class(p, "ggplot")
})
```

---

## Expected Outcomes

### Quantitative Improvements

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Total Lines of Code | 36,188 | ~29,500 | -18% |
| Integration Code | 2,500 | 400 | -84% |
| Dimension Reduction Code | 800 | 200 | -75% |
| Plotting Redundancy | 3,000 | 500 | -83% |
| Files to Update (Integration) | 15 | 1 | -93% |
| Files to Update (Dim Reduction) | 11 | 1 | -91% |

### Qualitative Improvements

1. **Maintainability**: Single point of change for shared logic
2. **Consistency**: Uniform behavior across similar functions
3. **Extensibility**: Easy to add new methods
4. **Testing**: Easier to write comprehensive tests
5. **Documentation**: Clearer, more organized docs
6. **Performance**: Opportunities for optimization in centralized code

### User Impact

- **Zero Breaking Changes**: All existing code continues to work
- **Improved Documentation**: Clearer organization
- **New Capabilities**: Unified interface enables easier method comparison
- **Better Error Messages**: Centralized validation provides consistent feedback

---

## Timeline Estimate

| Phase | Tasks | Estimated Time | Dependencies |
|-------|-------|----------------|--------------|
| Phase 1A | Integration consolidation | 2-3 days | None |
| Phase 1B | Dimension reduction consolidation | 2-3 days | None |
| Phase 1C | Plotting engine creation | 4-5 days | None |
| Phase 2 | Medium priority items | 2-3 days | Phase 1 complete |
| Phase 3 | File reorganization | 1-2 days | None |
| Phase 4 | Testing & documentation | 3-4 days | All phases |
| **Total** | **Complete refactoring** | **14-20 days** | |

**Recommendation**: Implement Phase 1 in 3 parallel branches, then merge and proceed with Phases 2-4.

---

## Appendix: Code Review Checklist

Before merging refactored code:

- [ ] All existing tests pass
- [ ] New tests added for consolidated functions
- [ ] Backward compatibility verified
- [ ] Documentation updated
- [ ] Examples run successfully
- [ ] No performance regressions
- [ ] Code coverage maintained or improved
- [ ] NAMESPACE correctly updated
- [ ] No new dependencies added
- [ ] Git history clean and logical

---

*This refactoring plan prioritizes maintainability and follows R package best practices while ensuring zero breaking changes for existing users.*
