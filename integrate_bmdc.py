import scanpy as sc
import scvi
import pandas as pd
import numpy as np
from scipy.io import mmread
import anndata as ad
import os
import tempfile
from pathlib import Path

import anndata as ad
import mudata as md
import muon
import numpy as np
import pooch
import scanpy as sc
import scvi
import seaborn as sns
import torch

# Set up scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')
torch.set_float32_matmul_precision("high")

# Function to load each dataset
def load_seurat_export(i):
    # Load count matrix
    counts = mmread(f"counts/counts_{i}.mtx").T.tocsr()  # Transpose to cells x genes
    
    # Load features
    features = pd.read_csv(f"counts/features_{i}.tsv", sep="\t", header=None, names=["gene_ids", "gene_symbols"])
    
    # Load barcodes  
    barcodes = pd.read_csv(f"counts/barcodes_{i}.tsv", sep="\t", header=None, names=["barcode"])
    
    # Load metadata
    metadata = pd.read_csv(f"counts/metadata_{i}.csv", index_col=0)
    
    # Create AnnData object
    adata = ad.AnnData(X=counts, obs=metadata, var=features.set_index("gene_ids"))
    adata.var_names = features["gene_symbols"].values
    adata.obs_names = barcodes["barcode"].values
    
    return adata

# Load the three datasets
adata1 = load_seurat_export(1)
adata2 = load_seurat_export(2) 
adata3 = load_seurat_export(3)

# Add batch labels with actual condition names
adata1.obs['batch'] = 'Standard_BM'
adata2.obs['batch'] = 'TGF_beta_treated'
adata3.obs['batch'] = 'IL_10_treated'

# Concatenate datasets
adata = sc.concat([adata1, adata2, adata3], label='batch', keys=['Standard_BM', 'TGF_beta_treated', 'IL_10_treated'])
adata.write('concatenated_datasets.h5ad')


#preprocessing
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata  # keep full dimension safe
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",
    layer="counts",
    subset=True,
    batch_key="batch",  # Change depending on the batch key for your dataset
)

adata.write_h5ad("preprocessed.h5ad")


# Make a copy to avoid view issues with scvi
adata = adata.copy()

# Basic preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalize and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='batch')
adata.raw = adata
adata = adata[:, adata.var.highly_variable].copy()

# Setup scvi model
scvi.model.SCVI.setup_anndata(adata, layer=None, batch_key='batch')

# Train scvi model
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train(max_epochs=400, patience=20)

# Get integrated latent representation
adata.obsm["X_scvi"] = model.get_latent_representation()

# Compute neighborhood graph and UMAP
sc.pp.neighbors(adata, use_rep="X_scvi")
sc.tl.umap(adata, min_dist=0.3)

# Leiden clustering
sc.tl.leiden(adata, key_added="leiden_scvi", resolution=0.5)

# Plot results
sc.pl.umap(adata, color=['batch', 'leiden_scvi'], ncols=2, frameon=False, save='_integration.pdf')

# Save integrated data
adata.write('integrated_data.h5ad')