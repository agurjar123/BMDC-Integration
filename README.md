# BMDC Integration Analysis

A comprehensive pipeline for integrating and analyzing single-cell RNA-seq data from bone marrow-derived cells (BMDCs) across different treatment conditions using scVI-tools.

## Overview

This repository contains tools for batch integration and analysis of single-cell RNA-seq data from three experimental conditions:
- Standard bone marrow (BM)
- TGF-β treated samples
- IL-10 treated samples

The pipeline leverages **scVI** (single-cell Variational Inference) for robust batch correction and integration, combined with **scanpy** for downstream analysis and visualization.

## Features

- **Automated data loading** from Seurat-exported count matrices
- **Batch integration** using scVI deep generative models
- **Quality control** and preprocessing workflows
- **Dimensionality reduction** with UMAP
- **Clustering analysis** using Leiden algorithm
- **Visualization** of integrated datasets

## Requirements

- Python 3.8+
- CUDA-compatible GPU (recommended for scVI training)

See `requirements.txt` for complete dependency list.

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/BMDC-Integration.git
cd BMDC-Integration
```

2. Create a conda environment (recommended):
```bash
conda create -n bmdc-integration python=3.9
conda activate bmdc-integration
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Data Structure

Your data should be organized as follows:

```
counts/
├── counts_1.mtx          # Count matrix for sample 1
├── features_1.tsv        # Gene features for sample 1
├── barcodes_1.tsv        # Cell barcodes for sample 1
├── metadata_1.csv        # Cell metadata for sample 1
├── counts_2.mtx          # Count matrix for sample 2
├── features_2.tsv        # Gene features for sample 2
├── barcodes_2.tsv        # Cell barcodes for sample 2
├── metadata_2.csv        # Cell metadata for sample 2
├── counts_3.mtx          # Count matrix for sample 3
├── features_3.tsv        # Gene features for sample 3
├── barcodes_3.tsv        # Cell barcodes for sample 3
└── metadata_3.csv        # Cell metadata for sample 3
```

## Usage

### Running the Python Script

```bash
python integrate_bmdc.py
```

This script will:
1. Load all three datasets from the `counts/` directory
2. Perform quality control and filtering
3. Train the scVI model for batch integration
4. Generate UMAP embeddings
5. Perform Leiden clustering
6. Save integrated data to `integrated_data.h5ad`

### Running the Jupyter Notebook

For interactive analysis:

```bash
jupyter notebook AG_082225.ipynb
```

The notebook provides step-by-step visualization and analysis of the integration process.

## Output Files

- `concatenated_datasets.h5ad` - Raw concatenated data
- `preprocessed.h5ad` - Preprocessed data with normalized counts
- `integrated_data.h5ad` - Final integrated dataset with scVI embeddings
- `figures/umap_integration.pdf` - UMAP visualization of batch integration

## Pipeline Steps

### 1. Data Loading
Load count matrices, features, barcodes, and metadata from Seurat exports.

### 2. Quality Control
- Filter cells with < 200 genes
- Filter genes expressed in < 3 cells
- Calculate mitochondrial gene percentages
- Remove low-quality cells

### 3. Normalization
- Total count normalization (10,000 counts per cell)
- Log transformation
- Highly variable gene selection

### 4. scVI Integration
- Train deep generative model with batch correction
- Extract batch-corrected latent representations
- Default: 30 latent dimensions, 2 hidden layers

### 5. Downstream Analysis
- Neighborhood graph construction
- UMAP embedding
- Leiden clustering

## Configuration

Key parameters in the pipeline:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_latent` | 30 | Dimensionality of latent space |
| `n_layers` | 2 | Number of hidden layers |
| `max_epochs` | 400 | Maximum training epochs |
| `patience` | 20 | Early stopping patience |
| `n_top_genes` | 2000 | Number of highly variable genes |
| `leiden_resolution` | 0.5 | Clustering resolution |

## Citation

If you use this pipeline, please cite:

**scVI-tools:**
```
Lopez, R., Regier, J., Cole, M. et al.
Deep generative modeling for single-cell transcriptomics.
Nat Methods 15, 1053–1058 (2018).
```

**Scanpy:**
```
Wolf, F., Angerer, P. & Theis, F.
SCANPY: large-scale single-cell gene expression data analysis.
Genome Biol 19, 15 (2018).
```

## License

MIT License - see LICENSE file for details

## Contact

For questions or issues, please open an issue on GitHub or contact [your email].

## Acknowledgments

This pipeline was developed for analyzing BMDC responses to immunomodulatory treatments.
