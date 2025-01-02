import scanpy as sc
import scvelo as scv
import numpy as np
import os
from config import input_path , output_path ,seed ,n_job

np.random.seed(seed)

# stochastic
input_file = input_path 
adata = sc.read_h5ad(input_file)
scv.tl.velocity(adata,n_jobs=n_job) # default: stochastic
adata.write_h5ad(output_path + "scvelo.sto.h5ad")

# dynamic
input_file = input_path 
adata = sc.read_h5ad(input_file)
scv.tl.recover_dynamics(adata, n_jobs=n_job,var_names= "all")
scv.tl.velocity(adata, mode='dynamical',n_jobs=n_job)
adata.write_h5ad(output_path + "scvelo.dyn.h5ad")