"""
Configuration file for BMDC integration analysis.

Modify these parameters according to your experimental setup and computational resources.
"""

# Data paths
DATA_DIR = "counts"
OUTPUT_DIR = "."
FIGURES_DIR = "figures"

# Sample information
SAMPLE_IDS = [1, 2, 3]
BATCH_LABELS = {
    1: 'Standard_BM',
    2: 'TGF_beta_treated',
    3: 'IL_10_treated'
}

# Quality control parameters
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 3
MAX_MT_PERCENT = 20  # Maximum percentage of mitochondrial genes

# Preprocessing parameters
TARGET_SUM = 1e4  # Normalization target sum
N_TOP_GENES = 2000  # Number of highly variable genes

# scVI model parameters
SCVI_CONFIG = {
    'n_layers': 2,          # Number of hidden layers
    'n_latent': 30,         # Dimensionality of latent space
    'gene_likelihood': 'nb' # Gene likelihood: 'nb' (negative binomial) or 'zinb'
}

# Training parameters
TRAINING_CONFIG = {
    'max_epochs': 400,      # Maximum number of training epochs
    'patience': 20,         # Early stopping patience
    'batch_size': 128,      # Training batch size
    'use_gpu': True         # Use GPU if available
}

# UMAP parameters
UMAP_CONFIG = {
    'min_dist': 0.3,        # Minimum distance in UMAP
    'n_neighbors': 15,      # Number of neighbors for UMAP
    'metric': 'euclidean'   # Distance metric
}

# Clustering parameters
LEIDEN_RESOLUTION = 0.5     # Resolution for Leiden clustering

# Visualization parameters
FIGURE_DPI = 80
FIGURE_FORMAT = 'pdf'       # Output format: 'pdf', 'png', 'svg'

# Random seed for reproducibility
RANDOM_SEED = 42
