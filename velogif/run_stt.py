from tools import sctt as st
from tools import networks
import scanpy as sc
import numpy as np
import os
from config import input_path , output_path  ,seed , data_cluster,embed

np.random.seed(seed)

input_file = input_path

adata = sc.read_h5ad(input_file)
adata_orig = adata.copy()

adata.obs['attractor'] = adata.obs[data_cluster].values
n_states = len(adata.obs['attractor'].values.unique())
adata_aggr = st.dynamical_iteration(adata,return_aggr_obj=True,n_states = n_states)

adata_aggr.obs_names = adata.obs_names

gene1 = set(adata_aggr.var_names)
gene2 = set(adata_orig.var_names)
gene = list(gene1 & gene2)

adata_aggr= adata_aggr[:,gene ]
adata_orig = adata_orig[:,gene]

adata_aggr.obs_names = adata_orig.obs_names
adata_aggr.obs[data_cluster] = adata_orig.obs[data_cluster]
emb = 'X_'+embed
adata_aggr.obsm[emb] = adata_orig.obsm[emb]

adata_aggr.write(output_path +'stt.h5ad')
