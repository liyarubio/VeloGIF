import pandas as pd
import scanpy as sc
import numpy  as np
import os
import sys
import dynamo as dyn
from config import input_path , output_path ,seed ,n_job,data_cluster

np.random.seed(seed)

input_file = input_path if input_path.endswith('.h5ad') else os.path.join(input_path, 'redeem_young.h5ad')

adata = sc.read_h5ad(input_file)

dyn.pp.recipe_monocle(adata)

try:
    dyn.tl.dynamics(adata, cores=n_job) # model auto:stochastic 
except Exception as e:
    print("stochastic mode did not converge, change to deterministic mode")
    dyn.tl.dynamics(adata, cores=n_job,model="deterministic")

print(adata)
adata.var = adata.var.iloc[:,0:10]
adata.uns['cell_phase_genes'] = None
adata.obs = adata.obs[[data_cluster]] # delete other columns to avoid save error
adata.layers['velocity'] = np.nan_to_num(adata.layers['velocity_S'].toarray())

adata.write_h5ad(output_path + "dynamo.h5ad")