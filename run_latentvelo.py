import latentvelo as ltv
import numpy as np
import scanpy as sc
import scvelo as scv
import pandas as pd
from config import input_path,output_path,seed,gene_number
import os

np.random.seed(seed)

input_file = input_path
adata = sc.read_h5ad(input_path )

adata.var['velocity_genes'] = True
adata.var['velocity_genes'].value_counts()

adata = ltv.utils.standard_clean_recipe(adata)
adata.var['velocity_genes'] = True

spliced_key = 'spliced'
unspliced_key = 'unspliced'

spliced_library_sizes = adata.layers[spliced_key].sum(1)
unspliced_library_sizes = adata.layers[unspliced_key].sum(1)

if len(spliced_library_sizes.shape) == 1:
       spliced_library_sizes = spliced_library_sizes[:,None]
if len(unspliced_library_sizes.shape) == 1:
       unspliced_library_sizes = unspliced_library_sizes[:,None]

adata.obs['spliced_size_factor'] = spliced_library_sizes #spliced_all_size_factors
adata.obs['unspliced_size_factor'] = unspliced_library_sizes #unspliced_all_size_factors

model = ltv.models.VAE(observed = gene_number) # observed: number of genes
epochs, val_ae, val_traj = ltv.train(model,adata,name="latentvelo")

latent_adata, adata = ltv.output_results(model, adata, gene_velocity=True,embedding='umap')

adata.layers['velocity'] = adata.layers['velo_s']

adata.write_h5ad(output_path + "latentvelo.h5ad")
latent_adata.write_h5ad(output_path + "latentvelo_latent.h5ad")