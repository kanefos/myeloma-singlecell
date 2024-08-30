
import argparse
import sys
import os
import tempfile
import math
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

import numba
import scvi
import torch
import pynndescent

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--adata_ref_dir', type=str)
parser.add_argument('--model_ref_dir', type=str)
parser.add_argument('--adata_query_dir', type=str)
parser.add_argument('--run_name', type=str)
parser.add_argument("--label_keys", nargs="+", type=str)
parser.add_argument('--surgery_epochs', type=int)
parser.add_argument('--scanvi_label', type=str)
args = parser.parse_args()


# Weighted Prediction Function ################################################

def weighted_prediction(weights, ref_cats):
    """Get highest weight category."""
    N = len(weights)
    predictions = np.zeros((N,), dtype=ref_cats.dtype)
    uncertainty = np.zeros((N,))
    for i in range(N):
        obs_weights = weights[i]
        obs_cats = ref_cats[i]
        best_prob = 0
        for c in np.unique(obs_cats):
            cand_prob = np.sum(obs_weights[obs_cats == c])
            if cand_prob > best_prob:
                best_prob = cand_prob
                predictions[i] = c
                uncertainty[i] = max(1 - best_prob, 0)

    return predictions, uncertainty


# Import #######################################################################

adata_ref = sc.read(args.adata_ref_dir)
adata_ref = adata_ref[:,adata_ref.var['intHVG']]

model = scvi.model.SCVI.load(args.model_ref_dir, adata_ref.copy())
adata_query = sc.read(args.adata_query_dir)

# scANVI model #################################################################

#Try to obtain a better latent representation/predictions by using the labels to inform 
#the latent space. This is where scANVI comes in. scANVI uses semi-supervised learning 
#to improve the model learned with scVI, allowing us to transfer our cell type knowledge
#from the reference to the query data.

# scANVI
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model, unlabeled_category='Unknown', labels_key='pheno', adata=adata_ref.copy())

scanvi_model.to_device('cuda:0') # moves model to GPU 0, cuda:0 just moves to the first available cuda 
scanvi_model.train(args.surgery_epochs) # similar n

# Learn a neighbors index on reference latent space
X_train = adata_ref.obsm["X_scVI"]
ref_nn_index = pynndescent.NNDescent(X_train)
ref_nn_index.prepare()

# Build and train query model ##################################################

# Prep the reference for scvi
scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_model)
adata_query.obs["scanvi_label"] = "unlabeled"

# Build model for query dataset using reference model
adata_query.obs['pheno'] = 'Unknown'
adata_query.obs['scanvi_label'] = 'Unknown'
query_model = scvi.model.SCANVI.load_query_data(adata_query, scanvi_model)

# Move to cuda GPU
query_model.to_device('cuda:0') # moves model to GPU 0, cuda:0 just moves to the first available cuda 

# scArches/scANVI-specific query training arguments.
train_kwargs_surgery = {
    "early_stopping": True,
    "early_stopping_monitor": "elbo_train",
    "early_stopping_patience": 10,
    "early_stopping_min_delta": 0.001,
    "plan_kwargs": {"weight_decay": 0.0},
}

# Train query model
query_model.train(max_epochs=args.surgery_epochs, **train_kwargs_surgery)


# Query model output #############################################################

# Get latent representation of query model
query_emb = ad.AnnData(query_model.get_latent_representation())
query_emb.obs_names = adata_query.obs_names

# mde
from scvi.model.utils import mde
query_emb.obsm["X_mde"] = mde(query_emb.X)
query_emb.obs['X_mde1'] = query_emb.obsm['X_mde'][:,0]
query_emb.obs['X_mde2'] = query_emb.obsm['X_mde'][:,1]

# Nearest Neighbors Search of query within reference KNN
ref_neighbors, ref_distances = ref_nn_index.query(query_emb.X)
#find the nearest neighbors (ref_neighbors) and their corresponding distances (ref_distances) 
#for each point in the query set (query_emb.X) based on a reference dataset. The ref_nn_index 
#seems to be a precomputed nearest neighbors index.

# Convert distances to affinities
stds = np.std(ref_distances, axis=1)
stds = (2.0 / stds) ** 2
stds = stds.reshape(-1, 1)
ref_distances_tilda = np.exp(-np.true_divide(ref_distances, stds))
weights = ref_distances_tilda / np.sum(ref_distances_tilda, axis=1, keepdims=True)

# For each annotation level, get prediction and uncertainty
for l in args.label_keys:
    ref_cats = adata_ref.obs[l].cat.codes.to_numpy()[ref_neighbors]
    p, u = weighted_prediction(weights, ref_cats)
    p = np.asarray(adata_ref.obs[l].cat.categories)[p]
    query_emb.obs[l + "_pred"], query_emb.obs[l + "_uncertainty"] = p, u

# Filter our predictions on the uncertainty threshold
#Using in the HLCA manuscript:
#uncertainty_threshold = 0.2
#Set higher for now, can edit later manually
uncertainty_threshold = 0.5
for l in args.label_keys:
    mask = query_emb.obs[l + "_uncertainty"] > uncertainty_threshold
    print(f"{l}: {sum(mask)/len(mask)} unknown")
    query_emb.obs[l + "_pred"].loc[mask] = "Unknown"

# Save outputs ############################################################################
    
# Output
query_emb.obs.to_csv(f'~/Scratch/IntegratedBM/Surgery_output/{args.run_name}.csv')
adata_query.obs.to_csv(f'~/Scratch/IntegratedBM/Surgery_output/{args.run_name}-obs.csv')
