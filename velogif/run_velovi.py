import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
import os
from velovi import preprocess_data, VELOVI
import matplotlib.pyplot as plt
from config import input_path,output_path,seed

torch.manual_seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(seed)

input_file = input_path
adata = sc.read_h5ad(input_file)

VELOVI.setup_anndata(adata, spliced_layer="Ms", unspliced_layer="Mu")
vae = VELOVI(adata)
vae.train()

latent_time = vae.get_latent_time()
velocities = vae.get_velocity(velo_mode='spliced')
velocity_u = vae.get_velocity(velo_mode='unspliced')

t = latent_time
scaling = 20 / t.max(0)

adata.layers["velocity"] = velocities / scaling
adata.layers['velocity_u'] = velocity_u / scaling
adata.layers["latent_time_velovi"] = latent_time

adata.write_h5ad(output_path + "velovi.h5ad")