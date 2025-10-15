# SCP - Single Cell Pipeline (Seurat V5 Compatible)

An adaptation of the [SCP package](https://github.com/zhanghao-njmu/SCP) with full Seurat V5 support and modern Python dependency management.

[![R](https://img.shields.io/badge/R-%3E%3D4.1.0-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-%3E%3D5.0.0-green.svg)](https://satijalab.org/seurat/)
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

```r
# Install devtools if needed
install.packages("devtools")

# Install SCP
devtools::install_github("mianaz/SCP-SeuratV5")
```

### Python Environment Setup

SCP uses Python for advanced analysis features, managed exclusively via UV (ultra-fast package manager).

```r
library(SCP)

# Install core dependencies only (~5 seconds)
PrepareEnv()

# Install all features (~60 seconds)
PrepareEnv(extras = "all")

# Install specific features
PrepareEnv(extras = c("velocity", "trajectory", "spatial"))

# Add more features to existing environment later
uv_install_extras("deeplearning")
```

**Performance**: UV is 10-100x faster than pip/conda!

For detailed Python setup instructions, see [PYTHON_DEPENDENCIES.md](PYTHON_DEPENDENCIES.md).

## Quick Start

### Basic Workflow

```r
library(SCP)

# Load your Seurat object (automatic V5 conversion if needed)
srt <- readRDS("your_seurat_object.rds")

# Check if it's Seurat V5
IsSeurat5(srt)

# Convert to V5 if needed
srt <- EnsureSeurat5(srt)

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

# Variable features
var_genes <- get_var_features(srt)
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
# One-time install
PrepareEnv(extras = c("spatial", "velocity"))

# In new sessions, activate environment
library(SCP)
use_uv_env()
```

See [PYTHON_DEPENDENCIES.md](PYTHON_DEPENDENCIES.md) for the complete list.

## Key Functions

### Version Management
- `IsSeurat5(srt)` - Check if object is Seurat V5
- `EnsureSeurat5(srt)` - Convert to V5 if needed
- `GetSeuratVersion(srt)` - Get detailed version info

### Data Access (Version-Agnostic)
- `get_seurat_data(srt, layer)` - Get data from any layer
- `set_seurat_data(srt, data, layer)` - Set data to a layer
- `get_feature_metadata(srt)` - Get feature metadata
- `set_feature_metadata(srt, metadata)` - Set feature metadata
- `get_var_features(srt)` - Get variable features

### Python Environment
- `PrepareEnv(extras)` - Set up UV Python environment with optional feature groups
- `use_uv_env()` - Configure R to use UV environment
- `check_uv()` - Check if UV is installed
- `uv_env_exists()` - Check if UV environment exists
- `uv_install_extras(extras)` - Add optional feature groups to existing environment
- `uv_install_packages(packages)` - Install individual Python packages
- `RemoveEnv()` - Remove UV Python environment
- `check_Python(packages)` - Verify Python packages are available

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
library(SCP)

# Activate UV environment
use_uv_env()

# Check specific packages
check_Python(c("numpy", "pandas", "scanpy", "anndata"))

# Check UV installation
check_uv()
uv_env_exists()
```

### Quick Test

```r
library(SCP)

# Test basic functionality
IsSeurat5(your_object)

# Test Python environment
use_uv_env()
py_config <- reticulate::py_config()
print(py_config)
```

## Troubleshooting

### Python Environment Issues

```r
# Remove and recreate environment
RemoveEnv()
PrepareEnv(force = TRUE)

# Or manually
uv_remove_env()
PrepareEnv(extras = "all")

# Verify installation
use_uv_env()  # Activate environment
check_Python(c("numpy", "scanpy", "anndata"))  # Check packages
```

### Seurat V5 Conversion Issues

```r
# Get detailed version info
GetSeuratVersion(srt)

# Force update to V5
srt <- UpdateSeuratObject(srt)
srt <- EnsureSeurat5(srt)
```

### Missing Optional Dependencies

```r
# Install specific feature groups to existing environment
uv_install_extras(c("velocity", "trajectory"))

# Or reinstall with new extras
PrepareEnv(extras = c("velocity", "trajectory", "spatial"), update = TRUE)
```

## Performance Tips

1. **Modular Installation**: Install only specific `extras` you need instead of "all" for faster setup
2. **Seurat V5 Native**: Use V5 objects for better memory efficiency
3. **Apple Silicon**: Use `extras = "apple_silicon"` for Metal acceleration in deep learning
4. **UV Speed**: Python environment setup is 10-100x faster than traditional pip/conda

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

If you use SCP in your research, please cite the original package:

```
Zhang, H. (2024). SCP: Single Cell Pipeline.
GitHub: https://github.com/zhanghao-njmu/SCP
```

## License

GPL-3.0-or-later

## Acknowledgments

- Original SCP package by [Hao Zhang](https://github.com/zhanghao-njmu)
- Seurat team for Seurat V5
- UV team at Astral for the amazing Python package manager

## Related Projects

- [Seurat](https://satijalab.org/seurat/) - Single-cell RNA-seq analysis
- [scanpy](https://scanpy.readthedocs.io/) - Python-based single-cell analysis
- [scVelo](https://scvelo.readthedocs.io/) - RNA velocity analysis
- [UV](https://docs.astral.sh/uv/) - Fast Python package manager

## Support

For issues and questions:
- **Bug Reports**: [GitHub Issues](https://github.com/mianaz/SCP-SeuratV5/issues)
- **Original SCP**: [SCP Issues](https://github.com/zhanghao-njmu/SCP/issues)
- **Python Dependencies**: See [PYTHON_DEPENDENCIES.md](PYTHON_DEPENDENCIES.md)

---

**Version**: 0.7.0
**Last Updated**: 2025-10-14
