import argparse
import scvi
import sys
import os
import scanpy as sc
import pandas as pd
import numpy as np
import torch

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--adata', help="Input adata.h5ad for integration",type=str)
parser.add_argument('--run_name', type=str)
parser.add_argument('--scvi_batch_key', type=str)
parser.add_argument("--scvi_categorical_covariate_keys", nargs="+", type=str)
parser.add_argument('--n_HVG', type=int)
parser.add_argument("--gene_groups_remove", nargs="+", type=str)
parser.add_argument('--n_latent', type=int)
parser.add_argument('--n_hidden', type=int)
parser.add_argument('--n_layers', type=int)
parser.add_argument('--dropout_rate', type=float)
parser.add_argument('--gene_likelihood', type=str)
args = parser.parse_args()

#if no categ covar, make true None
if args.scvi_categorical_covariate_keys[0] == 'None':
    args.scvi_categorical_covariate_keys = None

# Data
adata = sc.read(args.adata)

# HVGs
dataHVG=adata.copy()
sc.pp.normalize_total(dataHVG, target_sum=1e4)
sc.pp.log1p(dataHVG) 
# Remove select gene groups
for i in args.gene_groups_remove:
    dataHVG = dataHVG[:, (dataHVG.var[i] != True)]
# Scanpy default HVG
sc.pp.highly_variable_genes(dataHVG, n_top_genes=args.n_HVG, subset=True, batch_key=args.scvi_batch_key)
# Subset adata to HVG, use this for integration
adata_hvg = adata[:,dataHVG.var_names].copy()
adata.var['intHVG'] = adata.var['gene'].isin(dataHVG.var_names)
del dataHVG

# Settup AnnData object for scvi models.
scvi.model.SCVI.setup_anndata(
    adata_hvg, batch_key=args.scvi_batch_key, 
    categorical_covariate_keys=args.scvi_categorical_covariate_keys)

# Initialize scvi model 
scvi_model = scvi.model.SCVI(
    adata_hvg, n_latent=args.n_latent, n_hidden=args.n_hidden, n_layers=args.n_layers, 
    dropout_rate=args.dropout_rate, gene_likelihood=args.gene_likelihood)

# Train scvi model
scvi_model.to_device('cuda:0') # moves model to GPU 0, cuda:0 just moves to the first available cuda 
n_epoch_scvi=int(400*(20000/adata.shape[0]))
scvi_model.train(n_epoch_scvi)

# Obtain the latent space, build KNN graph/UMAP
adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=10) 
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15, key_added='n_neighbors15') 
sc.tl.umap(adata)
adata.obs['UMAP1'] = adata.obsm['X_umap'][:,0]
adata.obs['UMAP2'] = adata.obsm['X_umap'][:,1]

# Clustering
clustering_res=[0.6,0.8,1.0,1.2,1.5,2.0]
[ sc.tl.leiden(adata, resolution=i, key_added='leiden_'+str(i)) for i in clustering_res ]

# Create pymde 
from scvi.model.utils import mde
import pymde
adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])

# Save model
scvi_model.to_device('cpu') # back to GPU before saving
scvi_model.save(args.run_name, overwrite=True) 

# Prep AnnData for downstream
#normalise counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#other pl
adata.obs['log10_counts'] = np.log10(adata.obs['total_counts'])
adata.obs['log10_n_genes_by_counts'] = np.log10(adata.obs['n_genes_by_counts'])

# Save adata
adata.write(f"/{args.run_name}.h5ad",compression='gzip')
