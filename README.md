# SCPNext: Next-Generation Single-Cell Pipeline

A comprehensive R package for single-cell data analysis, built on the foundation of [SCP](https://github.com/zhanghao-njmu/SCP) with enhanced Seurat V5 support and streamlined Python integration.

[![R](https://img.shields.io/badge/R-%3E%3D4.1.0-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-%3E%3D4.2.0-green.svg)](https://satijalab.org/seurat/)
[![Python](https://img.shields.io/badge/Python-3.10--3.13-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## What's New in SCPNext

This release brings significant improvements to stability, usability, and performance:

### Major Enhancements

- **Full Seurat V5 Support**: Native compatibility with Seurat V5 objects and workflows
- **Streamlined Python Integration**: Faster, more reliable Python dependency management using [uv](https://docs.astral.sh/uv/)
- **Enhanced Stability**: Fixed critical bugs in parallel processing and enrichment analysis
- **Comprehensive Documentation**: New getting-started vignette with 94 working examples
- **Modular Installation**: Install only the Python features you need (velocity, trajectory, deep learning, etc.)

### Key Features

**Quality Control & Preprocessing**
- Automated cell QC with multiple metrics
- Flexible normalization and scaling
- Batch effect correction

**Dimensional Reduction & Clustering**
- Multiple methods: PCA, UMAP, t-SNE, PHATE
- Graph-based and hierarchical clustering
- Integration with Seurat V5 workflows

**Differential Expression & Enrichment**
- Multiple statistical tests
- GO/KEGG pathway enrichment
- GSEA support with parallel processing

**Advanced Analysis** (Python-powered)
- RNA velocity (scVelo)
- Trajectory inference (Slingshot, Palantir, CellRank)
- Cell annotation (automated and reference-based)
- Deep learning integration (scvi-tools)

**Visualization**
- Publication-ready plots
- Interactive exploration via SCExplorer
- Heatmaps, dot plots, violin plots, and more

## Installation

### Install from GitHub

**Fast installation with pak**

```r
# Install pak if needed
install.packages("pak")

# Install SCPNext (fast!)
pak::pak("mianaz/SCP-SeuratV5")
```

**Alternative**

```r
# Install devtools if needed
install.packages("devtools")

# Install SCPNext
devtools::install_github("mianaz/SCP-SeuratV5")
```

### Python Environment Setup

SCPNext uses Python for advanced analysis features. Python dependency management is handled via [uv](https://docs.astral.sh/uv/), an ultra-fast package manager.

**First-time setup:**

```r
library(SCPNext)

# Option 1: Basic installation (recommended) - ~30 seconds
# Includes: Essential single-cell tools for standard analysis
PrepareEnv(extras = "basic")

# Option 2: Full installation - ~90 seconds
# Includes everything: velocity, trajectory, deep learning, machine learning
PrepareEnv(extras = "all")
```

**What's included:**

- **`basic`** (recommended for most users):
  - Essential single-cell tools (UMAP, PAGA, Scanorama, BBKNN, Leiden, Louvain, Scrublet)
  - All features needed for standard QC, normalization, clustering, and visualization

- **`all`** (for advanced users):
  - Everything in `basic`
  - RNA velocity (scVelo, CellRank)
  - Trajectory inference (PHATE, Palantir)
  - Deep learning integration (scvi-tools, PyTorch)
  - Machine learning (XGBoost)

**Activating Python in future sessions:**

```r
library(SCPNext)
use_uv_env()  # Activate Python environment
```

**Adding packages later:**

```r
# Add specific Python packages
uv_install(packages = "scanpy")

# Upgrade from basic to all
uv_install(extras = "all")
```

## Quick Start

See our comprehensive [getting-started vignette](vignettes/getting-started.Rmd) for a complete walkthrough.

## Known limitations and To-Dos
- Most functions mentioned in [getting-started vignette](vignettes/getting-started.Rmd) should work ok now with a comprehensive debugging effort. 
- A refactoring for the integration module is underway, and have not been thoroughly tested.
- Some trajectory methods (Monocle, CellRank, Palatir, Monocle3, etc) have not been thoroughly tested.
- scANVI integration is new and has not been thoroughly tested.
- Added imputation methods (ALRA, MAGIC) but have not been thoroughly tested.
- Enrichment part has been rewritten to replace deprecated packages. RunGSEA, feature annotation functions should work fine. EnrichPlot will produce empty plots for unknown reasons. Also I was not able to fully reproduce EnrichPlots in original SCP. Needs investigation.

## Requirements

### R Dependencies
- **R** >= 4.1.0
- **Seurat** >= 4.0.0
- **SeuratObject** >= 4.0.0

## License

GPL-3.0-or-later

## Acknowledgments

SCPNext is built upon the excellent [SCP package](https://github.com/zhanghao-njmu/SCP) by Hao Zhang. We are deeply grateful to the original author for creating such a comprehensive and well-designed toolkit for single-cell analysis.


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

**Version**: 0.1.0
**Last Updated**: 2025-10-22
