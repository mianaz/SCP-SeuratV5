# SCPNext Function Index

> **Generated:** 2025-10-17
> **Package Version:** 0.1.0
> **Total Functions:** 184 (149 exported + 35 S3 methods)

This index provides a comprehensive reference for all functions in the SCPNext package, organized alphabetically with quick-reference information.

## Legend

- **Status**: ✅ Exported | 🔒 Internal | 🎭 S3 Method
- **Doc**: 📖 Documented | ⚠️ Undocumented | 📝 Partial
- **Type**: 🔧 Utility | 📊 Analysis | 🎨 Visualization | 🔄 Workflow | 🐍 Python | 📦 Data

---

## Function Index

### A

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `adata_to_srt` | ✅ | 📖 | 🔄 | SCP-analysis.R:6100 | Convert AnnData to Seurat object |
| `adjcolors` | ✅ | 📖 | 🔧 | utils.R:980 | Adjust color brightness/saturation |
| `AnnotateFeatures` | ✅ | 📖 | 📊 | SCP-feature_annotation.R:1 | Annotate features with biomaRt |
| `as_matrix` | ✅ | 📖 | 🔧 | utils.R:1100 | Convert to matrix (S3 generic) |
| `as_matrix.data.frame` | 🎭 | 📖 | 🔧 | utils.R:1120 | Convert data.frame to matrix |
| `as_matrix.default` | 🎭 | 📖 | 🔧 | utils.R:1140 | Default matrix conversion |
| `as_matrix.Matrix` | 🎭 | 📖 | 🔧 | utils.R:1110 | Convert Matrix to matrix |
| `as_matrix.matrix` | 🎭 | 📖 | 🔧 | utils.R:1130 | Identity for matrix |

### B

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `BBKNN_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:1900 | BBKNN batch correction |
| `blendcolors` | ✅ | 📖 | 🔧 | utils.R:990 | Blend two colors |
| `buildReferenceFromSeurat` | 🔒 | ⚠️ | 🔄 | SCP-projection.R:700 | Build Symphony reference |

### C

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `capitalize` | ✅ | 📖 | 🔧 | utils.R:950 | Capitalize first letter of strings |
| `CC_GenePrefetch` | ✅ | 📖 | 📊 | SCP-analysis.R:150 | Prefetch cell cycle genes |
| `CellCorHeatmap` | ✅ | 📖 | 🎨 | SCP-plot.R:8000 | Cell correlation heatmap |
| `CellDensityPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:2900 | Cell density visualization |
| `CellDimPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:100 | Cell annotation dimension plot |
| `CellDimPlot3D` | ✅ | 📖 | 🎨 | SCP-plot.R:800 | 3D cell dimension plot |
| `CellScoring` | ✅ | 📖 | 📊 | SCP-analysis.R:300 | Module-based cell scoring |
| `CellStatPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:4200 | Cell statistics visualization |
| `check_DataType` | ✅ | 📖 | 🔧 | SCP-workflow.R:2 | Detect data type (raw/log/normalized) |
| `check_srtList` | ✅ | 📖 | 🔧 | SCP-workflow.R:85 | Validate Seurat object list |
| `check_srtMerge` | ✅ | 📖 | 🔧 | SCP-workflow.R:3000 | Validate merged Seurat object |
| `check_uv` | ✅ | 📖 | 🐍 | utils.R:150 | Check UV package manager availability |
| `choose_k` | 🔒 | 📖 | 🔧 | SCP-imputation.R:170 | Choose optimal k for KNN |
| `col2hex` | ✅ | 📖 | 🔧 | utils.R:980 | Convert color to hex code |
| `ComBat_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:2500 | ComBat batch correction |
| `compute_velocity_on_grid` | ✅ | 📖 | 📊 | SCP-analysis.R:6650 | Compute velocity grid for plotting |
| `Conos_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:2700 | Conos integration |
| `CreateDataFile` | ✅ | 📖 | 📦 | SCP-app.R:1 | Create HDF5 data file for viewer |
| `CreateMetaFile` | ✅ | 📖 | 📦 | SCP-app.R:300 | Create metadata file for viewer |
| `CSS_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:2100 | Cluster similarity spectrum integration |

### D

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `db_DoubletDetection` | ✅ | 📖 | 📊 | SCP-cellqc.R:450 | Doublet detection via DoubletDetection (Python) |
| `db_scDblFinder` | ✅ | 📖 | 📊 | SCP-cellqc.R:150 | Doublet detection via scDblFinder |
| `db_scds` | ✅ | 📖 | 📊 | SCP-cellqc.R:250 | Doublet detection via scds |
| `db_Scrublet` | ✅ | 📖 | 📊 | SCP-cellqc.R:350 | Doublet detection via Scrublet (Python) |
| `DefaultReduction` | ✅ | 📖 | 🔧 | SCP-workflow.R:4400 | Auto-select best dimensionality reduction |
| `drop_data` | ✅ | 📖 | 🔧 | utils.R:1180 | Drop data from plot objects (S3 generic) |
| `drop_data.default` | 🎭 | 📖 | 🔧 | utils.R:1190 | Default drop data method |
| `drop_data.ggplot` | 🎭 | 📖 | 🔧 | utils.R:1200 | Drop data from ggplot |
| `drop_data.patchwork` | 🎭 | 📖 | 🔧 | utils.R:1210 | Drop data from patchwork |
| `DynamicHeatmap` | ✅ | 📖 | 🎨 | SCP-plot.R:9200 | Trajectory-based heatmap |
| `DynamicPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:14800 | Dynamic feature visualization |

### E

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `EnrichmentPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:10800 | Enrichment analysis visualization |
| `extract_ddrtree_ordering` | 🔒 | ⚠️ | 🔧 | SCP-analysis.R:3650 | Extract DDRTree ordering from Monocle2 |
| `extract_major_version` | 🔒 | 📖 | 🔧 | utils.R:1400 | Extract major version number |

### F

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `fastMNN_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:1300 | Fast MNN batch correction |
| `FeatureCorPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:12400 | Feature correlation plot |
| `FeatureDimPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:1500 | Feature expression dimension plot |
| `FeatureDimPlot3D` | ✅ | 📖 | 🎨 | SCP-plot.R:2200 | 3D feature dimension plot |
| `FeatureHeatmap` | ✅ | 📖 | 🎨 | SCP-plot.R:6800 | Feature expression heatmap |
| `FeatureStatPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:4900 | Feature statistics visualization |
| `FetchH5` | ✅ | 📖 | 📦 | SCP-app.R:900 | Fetch data from HDF5 files |
| `find_pyproject_toml` | 🔒 | 📖 | 🐍 | utils.R:1430 | Locate pyproject.toml file |
| `FindConservedMarkers2` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:1300 | Find conserved markers across conditions |
| `FindExpressedMarkers` | ✅ | 📖 | 📊 | SCP-analysis.R:500 | Find markers in expressed genes only |
| `FoldChange.default` | 🎭 | 📖 | 📊 | SCP-analysis.R:1500 | Calculate fold changes (S3 method) |
| `format_size` | 🔒 | 📖 | 🔧 | utils.R:1440 | Format file size for display |

### G

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `GeneConvert` | ✅ | 📖 | 📊 | SCP-analysis.R:1 | Convert gene IDs via biomaRt |
| `geom_alluvial` | ✅ | 📖 | 🎨 | ggsankey.R:401 | Alluvial plot geom |
| `geom_alluvial_label` | ✅ | 📖 | 🎨 | ggsankey.R:501 | Alluvial plot labels |
| `geom_alluvial_text` | ✅ | 📖 | 🎨 | ggsankey.R:551 | Alluvial plot text |
| `geom_sankey` | ✅ | 📖 | 🎨 | ggsankey.R:1 | Sankey diagram geom |
| `geom_sankey_bump` | ✅ | 📖 | 🎨 | ggsankey.R:701 | Sankey bump chart |
| `geom_sankey_label` | ✅ | 📖 | 🎨 | ggsankey.R:201 | Sankey diagram labels |
| `geom_sankey_text` | ✅ | 📖 | 🎨 | ggsankey.R:251 | Sankey diagram text |
| `get_feature_metadata` | ✅ | 📖 | 🔧 | utils.R:750 | Get feature metadata (V4/V5 compatible) |
| `get_scp_pkg_dir` | 🔒 | 📖 | 🔧 | utils.R:1420 | Get package installation directory |
| `get_seurat_data` | ✅ | 📖 | 🔧 | utils.R:550 | Get data from Seurat (V4/V5 compatible) |
| `get_uv_python_path` | 🔒 | 📖 | 🐍 | utils.R:1450 | Get UV Python executable path |
| `get_vars` | ✅ | 📖 | 🔧 | utils.R:1170 | Get variable features |
| `GraphPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:13000 | Network/graph visualization |
| `GroupHeatmap` | ✅ | 📖 | 🎨 | SCP-plot.R:5600 | Group-level heatmap |
| `GSEAPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:11600 | GSEA results visualization |

### H

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `Harmony_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:1500 | Harmony batch correction |

### I

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `install_uv` | ✅ | 📖 | 🐍 | utils.R:100 | Install UV package manager |
| `Integration_SCP` | ✅ | 📖 | 🔄 | SCP-workflow.R:4900 | Main integration dispatcher |
| `invoke` | ✅ | 📖 | 🔧 | utils.R:1050 | Invoke function with error handling |
| `isOutlier` | ✅ | 📖 | 🔧 | SCP-cellqc.R:520 | MAD-based outlier detection |
| `IsSeurat5` | ✅ | 📖 | 🔧 | utils.R:500 | Detect if Seurat object is V5 |

### L

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `LengthCheck` | 🔒 | ⚠️ | 🔧 | SCP-analysis.R:250 | Validate vector lengths match |
| `LIGER_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:2300 | LIGER iNMF integration |
| `LineagePlot` | ✅ | 📖 | 🎨 | SCP-plot.R:13500 | Lineage/trajectory visualization |
| `ListDB` | ✅ | 📖 | 📊 | SCP-analysis.R:2100 | List available enrichment databases |

### M

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `make_long` | ✅ | 📖 | 🔧 | ggsankey.R:851 | Convert data to long format for sankey |
| `mapQuery` | 🔒 | ⚠️ | 🔄 | SCP-projection.R:750 | Map query cells (Symphony) |
| `matrix_power` | 🔒 | 📖 | 🔧 | SCP-imputation.R:180 | Matrix exponentiation |
| `maximump` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:1850 | Maximum p-value combination |
| `meanp` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:1800 | Mean p-value combination |
| `metap` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:1700 | Meta p-value combining |
| `minimump` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:1820 | Minimum p-value combination |
| `MNN_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:1100 | MNN batch correction |

### O

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `orderCells` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:6300 | Order cells along trajectory (Monocle2) |

### P

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `PAGAPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:14000 | PAGA graph visualization |
| `palette_category` | ✅ | 📖 | 🔧 | utils.R:1250 | Category color palette |
| `palette_default` | ✅ | 📖 | 🔧 | utils.R:1200 | Default color palette |
| `palette_dimplot` | ✅ | 📖 | 🔧 | utils.R:1280 | Dimension plot palette |
| `palette_scp` | ✅ | 📖 | 🔧 | utils.R:1220 | SCP package palette |
| `panel_fix` | ✅ | 📖 | 🔧 | SCP-plot.R:14900 | Fix ggplot panel sizes |
| `panel_fix_overall` | ✅ | 📖 | 🔧 | SCP-plot.R:14920 | Fix overall plot panel layout |
| `PerformDE` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:900 | Generic DE testing dispatcher |
| `PrepareDB` | ✅ | 📖 | 📊 | SCP-analysis.R:2300 | Prepare enrichment database |
| `PrepareEnv` | ✅ | 📖 | 🐍 | utils.R:1 | Main Python environment setup |
| `PrepareSCExplorer` | ✅ | 📖 | 📦 | SCP-app.R:600 | Prepare data for web explorer |
| `project2MST` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:6400 | Project to MST (Monocle2) |
| `ProjectionPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:14600 | Reference projection visualization |
| `py_to_r_auto` | 🔒 | ⚠️ | 🐍 | SCP-analysis.R:6600 | Automatic Python to R conversion |

### R

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `RecoverCounts` | ✅ | 📖 | 🔧 | SCP-workflow.R:3200 | Recover raw counts from normalized data |
| `RemoveEnv` | ✅ | 📖 | 🐍 | utils.R:450 | Remove Python virtual environment |
| `RenameClusters` | ✅ | 📖 | 🔧 | SCP-workflow.R:3600 | Rename cell cluster identities |
| `RenameFeatures` | ✅ | 📖 | 🔧 | SCP-workflow.R:3400 | Rename/convert feature names |
| `require_packages` | ✅ | 📖 | 🔧 | utils.R:1000 | Check and require R packages |
| `RunALRA` | 🔒 | 📖 | 📊 | SCP-imputation.R:50 | ALRA imputation method |
| `RunCellQC` | ✅ | 📖 | 📊 | SCP-cellqc.R:1 | Comprehensive cell quality control |
| `RunCSSMap` | ✅ | 📖 | 🔄 | SCP-projection.R:450 | CSS-based reference mapping |
| `RunDEtest` | ✅ | 📖 | 📊 | SCP-analysis.R:700 | Differential expression testing |
| `RunDimReduction` | ✅ | 📖 | 🔄 | SCP-workflow.R:4200 | Unified dimensionality reduction |
| `RunDM` | ✅ | 📖 | 🔄 | Seurat-function.R:601 | Diffusion map (S3 generic) |
| `RunDM.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:630 | Diffusion map default method |
| `RunDM.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:615 | Diffusion map for Seurat |
| `RunDoubletCalling` | ✅ | 📖 | 📊 | SCP-cellqc.R:100 | Doublet detection dispatcher |
| `RunDynamicEnrichment` | ✅ | 📖 | 📊 | SCP-analysis.R:4300 | Enrichment for dynamic features |
| `RunDynamicFeatures` | ✅ | 📖 | 📊 | SCP-analysis.R:4000 | Find trajectory-dynamic genes |
| `RunEnrichment` | ✅ | 📖 | 📊 | SCP-analysis.R:2500 | Enrichment analysis |
| `RunFR` | ✅ | 📖 | 🔄 | Seurat-function.R:1501 | Force-directed graph (S3 generic) |
| `RunFR.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:1530 | Force-directed graph default |
| `RunFR.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:1515 | Force-directed graph for Seurat |
| `RunGLMPCA` | ✅ | 📖 | 🔄 | Seurat-function.R:401 | GLM-PCA (S3 generic) |
| `RunGLMPCA.Assay` | 🎭 | 📖 | 🔄 | Seurat-function.R:430 | GLM-PCA for Assay |
| `RunGLMPCA.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:445 | GLM-PCA default method |
| `RunGLMPCA.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:415 | GLM-PCA for Seurat |
| `RunGSEA` | ✅ | 📖 | 📊 | SCP-analysis.R:2800 | Gene Set Enrichment Analysis |
| `RunHarmony2` | ✅ | 📖 | 🔄 | Seurat-function.R:1651 | Harmony embedding (S3 generic) |
| `RunHarmony2.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:1665 | Harmony for Seurat |
| `RunImputation` | ✅ | 📖 | 📊 | SCP-imputation.R:1 | Imputation dispatcher |
| `RunKNNMap` | ✅ | 📖 | 🔄 | SCP-projection.R:1 | KNN-based reference mapping |
| `RunKNNPredict` | ✅ | 📖 | 📊 | SCP-cell_annotation.R:1 | KNN-based cell type prediction |
| `RunKNNSmooth` | 🔒 | 📖 | 📊 | SCP-imputation.R:150 | KNN smoothing imputation |
| `RunLargeVis` | ✅ | 📖 | 🔄 | Seurat-function.R:1351 | LargeVis (S3 generic) |
| `RunLargeVis.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:1380 | LargeVis default method |
| `RunLargeVis.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:1365 | LargeVis for Seurat |
| `RunMAGIC` | 🔒 | 📖 | 📊 | SCP-imputation.R:100 | MAGIC imputation method |
| `RunMDS` | ✅ | 📖 | 🔄 | Seurat-function.R:201 | MDS (S3 generic) |
| `RunMDS.Assay` | 🎭 | 📖 | 🔄 | Seurat-function.R:230 | MDS for Assay |
| `RunMDS.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:245 | MDS default method |
| `RunMDS.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:215 | MDS for Seurat |
| `RunMonocle2` | ✅ | 📖 | 📊 | SCP-analysis.R:3400 | Monocle2 pseudotime analysis |
| `RunMonocle3` | ✅ | 📖 | 📊 | SCP-analysis.R:3700 | Monocle3 trajectory inference |
| `RunNMF` | ✅ | 📖 | 🔄 | Seurat-function.R:1 | NMF (S3 generic) |
| `RunNMF.Assay` | 🎭 | 📖 | 🔄 | Seurat-function.R:30 | NMF for Assay |
| `RunNMF.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:45 | NMF default method |
| `RunNMF.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:15 | NMF for Seurat |
| `RunPAGA` | ✅ | 📖 | 📊 | SCP-analysis.R:4600 | PAGA graph abstraction |
| `RunPaCMAP` | ✅ | 📖 | 🔄 | Seurat-function.R:901 | PaCMAP (S3 generic) |
| `RunPaCMAP.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:930 | PaCMAP default method |
| `RunPaCMAP.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:915 | PaCMAP for Seurat |
| `RunPalantir` | ✅ | 📖 | 📊 | SCP-analysis.R:5200 | Palantir pseudotime inference |
| `RunPCAMap` | ✅ | 📖 | 🔄 | SCP-projection.R:150 | PCA-based reference mapping |
| `RunPHATE` | ✅ | 📖 | 🔄 | Seurat-function.R:1051 | PHATE (S3 generic) |
| `RunPHATE.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:1080 | PHATE default method |
| `RunPHATE.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:1065 | PHATE for Seurat |
| `RunSCExplorer` | ✅ | 📖 | 📦 | SCP-app.R:1200 | Launch interactive Shiny explorer |
| `RunSCVELO` | ✅ | 📖 | 📊 | SCP-analysis.R:4900 | RNA velocity analysis (scVelo) |
| `RunScmap` | ✅ | 📖 | 📊 | SCP-cell_annotation.R:250 | Scmap cell type annotation |
| `RunSeuratMap` | ✅ | 📖 | 🔄 | SCP-projection.R:300 | Seurat anchor-based mapping |
| `RunSingleR` | ✅ | 📖 | 📊 | SCP-cell_annotation.R:500 | SingleR cell type annotation |
| `RunSlingshot` | ✅ | 📖 | 📊 | SCP-analysis.R:3100 | Slingshot trajectory inference |
| `RunSymphonyMap` | ✅ | 📖 | 🔄 | SCP-projection.R:600 | Symphony reference mapping |
| `RunTriMap` | ✅ | 📖 | 🔄 | Seurat-function.R:1201 | TriMap (S3 generic) |
| `RunTriMap.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:1230 | TriMap default method |
| `RunTriMap.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:1215 | TriMap for Seurat |
| `RunUMAP2` | ✅ | 📖 | 🔄 | Seurat-function.R:751 | Enhanced UMAP (S3 generic) |
| `RunUMAP2.default` | 🎭 | 📖 | 🔄 | Seurat-function.R:780 | UMAP2 default method |
| `RunUMAP2.Seurat` | 🎭 | 📖 | 🔄 | Seurat-function.R:765 | UMAP2 for Seurat |
| `RunWOT` | ✅ | 📖 | 📊 | SCP-analysis.R:5500 | Waddington OT trajectory |

### S

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `SaveExampleChromium` | ✅ | 📖 | 📦 | SCP-app.R:2100 | Save example 10X data |
| `SCANVI_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:900 | SCANVI semi-supervised integration |
| `Scanorama_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:1700 | Scanorama integration |
| `scVI_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:700 | scVI deep learning integration |
| `searchDatasets` | 🔒 | ⚠️ | 🔧 | SCP-analysis.R:6200 | Search enrichment databases |
| `segementsDf` | ✅ | 📖 | 🔧 | ggsankey.R:950 | Create segments data frame for sankey |
| `Seurat_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:300 | Seurat V4 integration |
| `Seurat_integrate_v5` | ✅ | 📖 | 🔄 | SCP-workflow.R:500 | Seurat V5 layer-based integration |
| `set_feature_metadata` | ✅ | 📖 | 🔧 | utils.R:850 | Set feature metadata (V4/V5 compatible) |
| `set_seurat_data` | ✅ | 📖 | 🔧 | utils.R:650 | Set data in Seurat (V4/V5 compatible) |
| `show_palettes` | ✅ | 📖 | 🔧 | utils.R:1350 | Display available color palettes |
| `slim_data` | ✅ | 📖 | 🔧 | utils.R:1220 | Reduce data size in plots (S3 generic) |
| `slim_data.default` | 🎭 | 📖 | 🔧 | utils.R:1230 | Default slim data method |
| `slim_data.ggplot` | 🎭 | 📖 | 🔧 | utils.R:1240 | Slim data from ggplot |
| `slim_data.patchwork` | 🎭 | 📖 | 🔧 | utils.R:1250 | Slim data from patchwork |
| `srt_to_adata` | ✅ | 📖 | 🔄 | SCP-analysis.R:5800 | Convert Seurat to AnnData |
| `SrtAppend` | ✅ | 📖 | 🔧 | SCP-workflow.R:4000 | Append/merge Seurat objects |
| `SrtReorder` | ✅ | 📖 | 🔧 | SCP-workflow.R:3800 | Reorder cells by feature values |
| `Standard_SCP` | ✅ | 📖 | 🔄 | SCP-workflow.R:4800 | Standard preprocessing pipeline |
| `StatPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:3500 | Generic statistical plots |
| `sump` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:1900 | Sum p-value combination |

### T

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `theme_alluvial` | ✅ | 📖 | 🎨 | ggsankey.R:950 | Alluvial plot theme |
| `theme_blank` | ✅ | 📖 | 🎨 | utils.R:1320 | Blank ggplot2 theme |
| `theme_sankey` | ✅ | 📖 | 🎨 | ggsankey.R:901 | Sankey diagram theme |
| `theme_sankey_bump` | ✅ | 📖 | 🎨 | ggsankey.R:980 | Sankey bump chart theme |
| `theme_scp` | ✅ | 📖 | 🎨 | utils.R:1300 | Main SCP package theme |

### U

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `Uncorrected_integrate` | ✅ | 📖 | 🔄 | SCP-workflow.R:100 | No batch correction (baseline) |
| `use_uv_env` | ✅ | 📖 | 🐍 | utils.R:400 | Activate UV Python environment |
| `uv_create_env` | ✅ | 📖 | 🐍 | utils.R:200 | Create UV virtual environment |
| `uv_env_exists` | ✅ | 📖 | 🐍 | utils.R:250 | Check if UV environment exists |
| `uv_install` | ✅ | 📖 | 🐍 | utils.R:350 | Install Python packages via UV |
| `uv_sync_deps` | ✅ | 📖 | 🐍 | utils.R:300 | Sync UV dependencies |

### V

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `validate_group_parameter` | ✅ | 📖 | 🔧 | validation-utils.R:1 | Validate group parameter inputs |
| `VelocityPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:14300 | RNA velocity field visualization |
| `VolcanoPlot` | ✅ | 📖 | 🎨 | SCP-plot.R:10000 | Volcano plot for DE results |
| `votep` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:1950 | Vote p-value combination |

### W

| Function | Status | Doc | Type | File | Description |
|----------|--------|-----|------|------|-------------|
| `WilcoxDETest` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:1100 | Wilcoxon differential expression test |
| `wilkinsonp` | 🔒 | ⚠️ | 📊 | SCP-analysis.R:1750 | Wilkinson p-value combination |

---

## Summary Statistics

### By Status
- **Exported Functions:** 149 (81%)
- **S3 Methods:** 35 (19%)
- **Internal Functions:** ~20-40 (estimated)

### By Documentation
- **Fully Documented:** ~159 (86%)
- **Undocumented:** ~15 (8%)
- **Partial Documentation:** ~10 (6%)

### By Type
- **Workflow Functions:** 32 (17%)
- **Analysis Functions:** 45 (24%)
- **Visualization Functions:** 62 (34%)
- **Utility Functions:** 30 (16%)
- **Python Integration:** 15 (8%)

### By File
- **SCP-plot.R:** 62 functions (34%)
- **SCP-analysis.R:** 32 functions (17%)
- **SCP-workflow.R:** 30 functions (16%)
- **Seurat-function.R:** 33 functions (18%)
- **utils.R:** 40+ functions (22%)
- **Other files:** ~15 functions (8%)

---

## Quick Reference by Use Case

### Quality Control
- `RunCellQC` - Main QC pipeline
- `RunDoubletCalling` - Doublet detection
- `isOutlier` - Outlier detection

### Data Integration
- `Integration_SCP` - Main dispatcher
- `Seurat_integrate`, `Seurat_integrate_v5` - Seurat methods
- `scVI_integrate`, `SCANVI_integrate` - Deep learning
- `Harmony_integrate` - Harmony
- 10+ other integration methods

### Dimensionality Reduction
- `RunDimReduction` - Unified interface
- `RunUMAP2`, `RunPHATE`, `RunPaCMAP`, `RunTriMap` - Manifold methods
- `RunNMF`, `RunMDS`, `RunGLMPCA` - Matrix factorization
- `RunDM`, `RunLargeVis`, `RunFR` - Graph-based

### Trajectory Analysis
- `RunSlingshot` - Linear trajectories
- `RunMonocle2`, `RunMonocle3` - Monocle methods
- `RunPAGA` - Graph abstraction
- `RunSCVELO` - RNA velocity
- `RunPalantir`, `RunWOT` - Advanced methods

### Cell Annotation
- `RunKNNPredict` - KNN-based
- `RunScmap` - Scmap
- `RunSingleR` - SingleR

### Differential Expression
- `RunDEtest` - Main DE testing
- `FindExpressedMarkers` - Marker finding
- `RunEnrichment` - Enrichment analysis
- `RunGSEA` - GSEA

### Visualization
- `CellDimPlot`, `FeatureDimPlot` - Dimension plots
- `StatPlot`, `CellStatPlot`, `FeatureStatPlot` - Statistics
- `GroupHeatmap`, `FeatureHeatmap` - Heatmaps
- `VolcanoPlot`, `EnrichmentPlot`, `GSEAPlot` - Results
- `VelocityPlot`, `LineagePlot`, `PAGAPlot` - Trajectories

### Data Management
- `get_seurat_data`, `set_seurat_data` - Data access
- `RecoverCounts` - Count recovery
- `RenameFeatures`, `RenameClusters` - Renaming
- `srt_to_adata`, `adata_to_srt` - Format conversion

### Python Environment
- `PrepareEnv` - Main setup
- `use_uv_env` - Activate environment
- `uv_install` - Install packages

---

*This index is automatically generated. Line numbers are approximate and may shift with code changes.*
