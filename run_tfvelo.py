import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use('AGG')
import os, sys
import random
import shutil
from config import input_path, output_path, seed, data_cluster, embed, n_job,tfvelo_path
np.random.seed(seed)
random.seed(seed)

# os.chdir(tfvelo_path)
sys.path.insert(0,tfvelo_path)
import TFvelo as TFv

def check_data_type(adata):
    for key in list(adata.var):
        if adata.var[key][0] in ['True', 'False']:
            adata.var[key] = adata.var[key].map({'True': True, 'False': False})
    return

def data_type_tostr(adata, key):
    if key in adata.var.keys():
        if adata.var[key][0] in [True, False]:
            adata.var[key] = adata.var[key].map({True: 'True', False:'False'})
    return

def run_tfvelo(adata, key):
    # os.chdir(tfvelo_path) # need same database path
    adata.obs['clusters'] = adata.obs[key].copy()
    adata.layers["total"] = adata.layers["spliced"].todense() + adata.layers["unspliced"].todense()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    n_cells, n_genes = adata.X.shape

    TFv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    adata.X = adata.layers["total"].copy()
    gene_names = []
    for tmp in adata.var_names:
        gene_names.append(tmp.upper())
    TFv.pp.moments(adata, n_pcs=30, n_neighbors=15)
    TFv.pp.get_TFs(adata, databases='ENCODE ChEA')
    print(adata)
    adata.uns['genes_pp'] = np.array(adata.var_names)
    flag = TFv.tl.recover_dynamics(adata, n_jobs=n_job, max_iter=20, var_names='all',
                                   WX_method="lsq_linear", WX_thres=20, max_n_TF=99, n_top_genes=2000,fit_scaling=True, use_raw=0, init_weight_method="correlation",
                                   n_time_points=1000)
    return

# script_path = os.path.abspath(__file__)
# script_dir = os.path.dirname(script_path)
# os.chdir(script_dir)

adata = sc.read_h5ad(input_path)

# copy database
try:
    shutil.copytree(tfvelo_path + 'ENCODE/', os.getcwd()+ '/ENCODE/')
except Exception as e:
    print(e)

try:
    shutil.copytree(tfvelo_path +'ChEA/', os.getcwd()+ '/ChEA/')
except Exception as e:
    print(e)

run_tfvelo(adata,data_cluster)

n_cells = adata.shape[0]
expanded_scaling_y = np.expand_dims(np.array(adata.var['fit_scaling_y']), 0).repeat(n_cells, axis=0)
adata.layers['velocity'] = np.nan_to_num(adata.layers['velo_hat'] / expanded_scaling_y)

adata.write(output_path +'tfvelo.h5ad')
