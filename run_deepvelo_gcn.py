import numpy as np
import scvelo as scv
import torch
import scanpy as sc
import os
from deepvelo.utils.scatter import scatter
from deepvelo.utils.preprocess import autoset_coeff_s
from deepvelo.utils.plot import statplot, compare_plot
from deepvelo import train, Constants
from deepvelo.utils import (
    velocity,
    velocity_confidence,
    continuity_confidence,
    update_dict,
    cross_boundary_correctness,)
from config import input_path,output_path,seed
os.environ['CUDA_LAUNCH_BLOCKING'] = '1'

torch.manual_seed(seed)
np.random.seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

input_file = input_path
print(f"Reading file: {input_file}")

adata = sc.read_h5ad(input_file)
print(adata)

configs = {
    "name": "DeepVelo_GB", # name of the experiment
    'n_gpu': 0, # whether use gpu
    "loss": {"args": {"coeff_s": autoset_coeff_s(adata)}},
    "arch":{'args': {'pred_unspliced': True}},
    "trainer": {"verbosity": 0}, # increase verbosity to show training progress
}
configs = update_dict(Constants.default_configs, configs)

# initial velocity
velocity(adata, mask_zero=False) # use scvelo stochastic
trainer = train(adata, configs)

print(adata)
adata.write_h5ad(output_path + "deepvelo_gcn.h5ad")