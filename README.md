# SCPNext: Next-Generation Single-Cell Pipeline

A next-generation evolution of the [SCP package](https://github.com/zhanghao-njmu/SCP) with full Seurat V5 support, Python 3.12+ compatibility, and modern dependency management.

[![R](https://img.shields.io/badge/R-%3E%3D4.1.0-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-%3E%3D4.2.0-green.svg)](https://satijalab.org/seurat/)
[![Python](https://img.shields.io/badge/Python-3.10--3.12-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Features

### 🎯 Seurat V5 Compatibility
- **Full V5 Support**: Native compatibility with Seurat V5's layer-based data structure
- **Version Agnostic**: Seamlessly handles both Seurat V4 and V5 objects
- **Automatic Detection**: Smart detection and conversion of legacy objects
- **Zero Migration Cost**: Existing pipelines work without modification

### ⚡ Fast Python Dependency Management
- **UV-Only**: Ultra-fast Python package manager (10-100x faster than pip/conda)
- **Modular Installation**: Install only the features you need
- **Cross-Platform**: Works on macOS, Linux, and Windows
- **Simple Setup**: One command to install all dependencies

### 🔬 Comprehensive Single-Cell Analysis
- RNA velocity analysis (scVelo, CellRank)
- Trajectory inference (Palantir, WOT)
- Advanced dimensionality reduction (UMAP, PHATE, PaCMAP, TriMAP)
- Batch correction (Scanorama, BBKNN, Harmony)
- Doublet detection (Scrublet, DoubletDetection)
- Deep learning methods (scvi-tools)
- And much more...

## Installation

### Install from GitHub

**Recommended: Fast installation with pak (10-100x faster)**

```r
# Install pak if needed
install.packages("pak")

# Install SCPNext (fast!)
pak::pak("mianaz/SCP-SeuratV5")
```

**Alternative: Traditional installation**

```r
# Install devtools if needed
install.packages("devtools")

# Install SCPNext
devtools::install_github("mianaz/SCP-SeuratV5")
```

### Python Environment Setup (One-Time Only)

SCP uses Python for advanced analysis features, managed exclusively via UV (ultra-fast package manager).

**First time setup:**
```r
library(SCPNext)

# Install core dependencies only (~5 seconds)
PrepareEnv()

# OR install all features (~60 seconds)
PrepareEnv(extras = "all")

# OR install specific features
PrepareEnv(extras = c("velocity", "trajectory", "spatial"))
```

**That's it!** In future R sessions, run `library(SCPNext)` and then `use_uv_env()` to activate Python features.

**Add features later:**
```r
library(SCPNext)
use_uv_env()
uv_install(extras = "deeplearning")
```

**Performance**: UV is 10-100x faster than pip/conda!

For detailed Python setup instructions, see [PYTHON_DEPENDENCIES.md](PYTHON_DEPENDENCIES.md).

## Quick Start

### Basic Workflow

```r
library(SCPNext)

# Load your Seurat object
srt <- readRDS("your_seurat_object.rds")

# Check if it's Seurat V5
IsSeurat5(srt)

# Convert to V5 if needed (using Seurat's built-in function)
if (!IsSeurat5(srt)) {
  srt <- UpdateSeuratObject(srt)
}

# Run standard analysis
srt <- RunStandardWorkflow(srt)

# Visualize
CellDimPlot(srt, group.by = "cell_type")
FeatureDimPlot(srt, features = c("CD3D", "CD8A", "CD4"))
```

### Working with Seurat V5

```r
# Version-agnostic data access
counts <- get_seurat_data(srt, layer = "counts")
data <- get_seurat_data(srt, layer = "data")

# Set data back
srt <- set_seurat_data(srt, data = normalized_data, layer = "data")

# Feature metadata
feature_meta <- get_feature_metadata(srt)
srt <- set_feature_metadata(srt, metadata = new_meta)

# Variable features (use Seurat's built-in function)
var_genes <- Seurat::VariableFeatures(srt)
```

### Python-Based Analysis

```r
# RNA Velocity with scVelo
srt <- RunScvelo(srt,
  group.by = "cell_type",
  linear_reduction = "pca",
  nonlinear_reduction = "umap"
)

# Trajectory with Palantir
srt <- RunPalantir(srt,
  start_cell = "CELL_001",
  num_waypoints = 500
)

# Cell fate with CellRank
srt <- RunCellRank(srt,
  group.by = "cell_type",
  linear_reduction = "pca"
)
```

## Requirements

### R Dependencies
- **R** >= 4.1.0
- **Seurat** >= 5.0.0
- **SeuratObject** >= 5.0.0
- Plus standard Bioconductor packages (automatically installed)

### Python Dependencies (Optional)
- **Python** 3.10, 3.11, or 3.12 (3.10 recommended)
- Managed via UV (ultra-fast package manager)
- Quick start:

```r
# First time only - install Python environment
library(SCPNext)
PrepareEnv(extras = c("spatial", "velocity"))

# Future sessions - load the package and activate Python
library(SCPNext)
use_uv_env()
```

See [PYTHON_DEPENDENCIES.md](PYTHON_DEPENDENCIES.md) for the complete list.

## Key Functions

### Version Management
- `IsSeurat5(srt)` - Check if object is Seurat V5
- Use Seurat's `UpdateSeuratObject(srt)` to convert to V5 if needed

### Data Access (Version-Agnostic)
- `get_seurat_data(srt, layer)` - Get data from any layer
- `set_seurat_data(srt, data, layer)` - Set data to a layer
- `get_feature_metadata(srt)` - Get feature metadata
- `set_feature_metadata(srt, metadata)` - Set feature metadata
- Use `Seurat::VariableFeatures(srt)` for variable features

### Python Environment
- `PrepareEnv(extras)` - Set up UV Python environment with optional feature groups (one-time)
- `uv_install(extras, packages)` - Install feature groups and/or individual Python packages to existing environment
- `RemoveEnv()` - Remove UV Python environment
- `check_uv()` - Check if UV is installed
- `uv_env_exists()` - Check if UV environment exists
- `use_uv_env()` - Activate Python environment (required before using Python-based functions)

### Analysis Functions
See the [original SCP documentation](https://github.com/zhanghao-njmu/SCP) for the full list of analysis functions.

### Examples/Tests
- Minimal scRNA smoke test: `tests/testthat/test_seurat_v5_basics.R`
- Optional Chromium test: set `SCP_TEST_10X_H5` to a 10x HDF5 path
- Optional Visium HD test: set `SCP_TEST_VISIUM_BASE` to a binned output dir (e.g., `.../square_008um`) and run `tests/interactive/test_spatial_visium_hd.R`

## Documentation

- **Python Dependencies**: [PYTHON_DEPENDENCIES.md](PYTHON_DEPENDENCIES.md)
- **UV Setup Guide**: [UV_SETUP.md](UV_SETUP.md)
- **Original SCP Docs**: [SCP Repository](https://github.com/zhanghao-njmu/SCP)

## Testing Your Installation

### Test Python Dependencies

```r
library(SCPNext)
use_uv_env()  # Activate Python environment

# Check if environment exists
uv_env_exists()

# Check Python configuration
py_config <- reticulate::py_config()
print(py_config)

# Test importing a module
reticulate::py_module_available("scanpy")
reticulate::py_module_available("anndata")
```

### Quick Test

```r
library(SCPNext)

# Test basic functionality
IsSeurat5(your_object)

# Test Python environment is active
reticulate::py_config()
```

## Troubleshooting

### Python Environment Issues

```r
# Remove and recreate environment
RemoveEnv()
PrepareEnv(force = TRUE, extras = "all")

# Verify installation
library(SCPNext)
use_uv_env()
reticulate::py_config()
reticulate::py_module_available("scanpy")
```

### Seurat V5 Conversion Issues

```r
# Check if object is V5
IsSeurat5(srt)

# Force update to V5 using Seurat's built-in function
srt <- UpdateSeuratObject(srt)

# Verify it's now V5
IsSeurat5(srt)
```

### Missing Optional Dependencies

```r
# Install specific feature groups to existing environment
uv_install(extras = c("velocity", "trajectory"))

# Or install individual packages
uv_install(packages = "scvelo")

# Or install both
uv_install(extras = "velocity", packages = "custom_package")

# Or reinstall with new extras
PrepareEnv(extras = c("velocity", "trajectory", "spatial"), update = TRUE)
```

## Performance Tips

1. **pak Installation**: Use `pak::pak()` instead of `devtools::install_github()` for 5-10x faster installation
2. **Modular Installation**: Install only specific `extras` you need instead of "all" for faster setup
3. **Seurat V5 Native**: Use V5 objects for better memory efficiency
4. **Apple Silicon**: Use `extras = "apple_silicon"` for Metal acceleration in deep learning
5. **UV Speed**: Python environment setup is 10-100x faster than traditional pip/conda

## Platform Support

| Platform | Python (UV) | Status |
|----------|-------------|--------|
| macOS (Intel) | ✅ | Full support |
| macOS (Apple Silicon) | ✅ | Full support + Metal acceleration |
| Linux | ✅ | Full support |
| Windows | ✅ | Full support |

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## Citation

If you use SCPNext in your research, please cite:

**SCPNext:**
```
Zhang, H. & Contributors (2025). SCPNext: Next-Generation Single Cell Pipeline.
GitHub: https://github.com/mianaz/SCP-SeuratV5
```

**Original SCP package:**
```
Zhang, H. (2024). SCP: Single Cell Pipeline.
GitHub: https://github.com/zhanghao-njmu/SCP
```

## License

GPL-3.0-or-later

## Acknowledgments

SCPNext is built upon the excellent [SCP package](https://github.com/zhanghao-njmu/SCP) by Hao Zhang. We are deeply grateful to the original author for creating such a comprehensive and well-designed toolkit for single-cell analysis.

**What SCPNext adds:**
- Full Seurat V5 compatibility and layer-based architecture support
- Python 3.12+ support with modern dependency management via UV
- Apple Silicon (M1/M2/M3/M4) optimization and testing
- Enhanced integration methods and updated API compatibility
- Modernized codebase for long-term maintainability

**Additional thanks to:**
- Seurat team for the Seurat V5 architecture
- Posit team for the pak R package manager
- Astral team for the UV Python package manager
- Single-cell community for feedback and contributions

## Related Projects

- [Seurat](https://satijalab.org/seurat/) - Single-cell RNA-seq analysis
- [scanpy](https://scanpy.readthedocs.io/) - Python-based single-cell analysis
- [scVelo](https://scvelo.readthedocs.io/) - RNA velocity analysis
- [pak](https://pak.r-lib.org/) - Fast R package manager (recommended for installation)
- [UV](https://docs.astral.sh/uv/) - Fast Python package manager

## Support

For issues and questions:
- **Bug Reports**: [GitHub Issues](https://github.com/mianaz/SCP-SeuratV5/issues)
- **Original SCP**: [SCP Issues](https://github.com/zhanghao-njmu/SCP/issues)
- **Python Dependencies**: See [PYTHON_DEPENDENCIES.md](PYTHON_DEPENDENCIES.md)

---

**Version**: 0.7.0
**Last Updated**: 2025-10-14
