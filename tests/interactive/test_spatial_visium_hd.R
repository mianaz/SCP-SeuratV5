#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(reticulate)
})

for (f in list.files('R', pattern='\\.R$', full.names = TRUE)) source(f, chdir = TRUE)

base <- Sys.getenv("SCP_TEST_VISIUM_BASE", unset = NA)
if (is.na(base) || !dir.exists(base)) {
  stop("Set SCP_TEST_VISIUM_BASE to Visium HD base (directory containing 'spatial' and matrix). Example: binned_outputs/square_008um")
}

EnsureEnv(required = c('scanpy','squidpy','anndata'))

set.seed(1)
srt <- LoadSpatialData(base, assay = 'Spatial')
if (ncol(srt) > 5000) srt <- subset(srt, cells = sample(colnames(srt), 5000))
srt <- PreprocessSpatial(srt, normalization_method = 'LogNormalize', npcs = 20, run_umap = TRUE, run_clustering = TRUE, verbose = TRUE)
srt <- RunSpatialFeatures(srt, method = 'markvariogram', verbose = TRUE)
srt <- RunSpatialGraph(srt, method = 'knn', n_neighbors = 6, verbose = TRUE)
srt <- RunSpatialEnrichment(srt, cluster_key = 'seurat_clusters', n_perms = 10, verbose = TRUE)
srt <- RunSpatialAutocorr(srt, mode = 'moran', n_perms = 10, n_jobs = 1, verbose = TRUE)
srt <- IdentifySpatialDomains(srt, method='leiden', resolution = 0.5, verbose = TRUE)

print("DONE: spatial pipeline on Visium HD subset")

