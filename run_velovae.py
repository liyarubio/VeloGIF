import os
import scanpy as sc
import numpy as np
import sys
import torch
import scvelo as scv
from scipy.sparse import csr_matrix
import velovae as vv
from config import input_path, output_path, seed, device, embed, data_cluster, gene_number

np.random.seed(seed)
torch.manual_seed(seed)

input_file = input_path

adata = sc.read_h5ad(input_file)  

adata.layers['spliced'] = csr_matrix(adata.layers['spliced'])
adata.layers['unspliced'] = csr_matrix(adata.layers['unspliced'])

vv.preprocess(adata, n_gene=gene_number)

vae = vv.VAE(adata, 
             tmax=20, 
             dim_z=len(set(adata.obs[data_cluster])),
             device=device)

# default
config = {
    # You can change any hyperparameters here!
    # 'learning_rate': 1e-3,
    # 'learning_rate_ode': 2e-3,
    # 'learning_rate_post': 1e-3
}

vae.train(adata,
          config=config,
          plot=False,
          embed=embed)

vae.save_model(output_path + "/velovae", 'encoder_vae', 'decoder_vae')
vae.save_anndata(adata, 'vae', output_path, file_name="velovae.h5ad")
adata.layers['velocity'] = adata.layers['vae_velocity']
adata.write(output_path + "/velovae.h5ad")