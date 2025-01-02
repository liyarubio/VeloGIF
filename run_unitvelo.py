import pandas as pd
import numpy as np
import scanpy as sc
import scvelo as scv
import unitvelo as utv
import tensorflow as tf
import os
os.environ['TF_USE_LEGACY_KERAS'] = 'True'
from config import input_path , output_path  ,seed ,device , data_cluster

np.random.seed(seed)
tf.random.set_seed(seed)

adata = sc.read_h5ad(input_path)
adata.var['highly_variable'] = True # use all genes

velo_config = utv.config.Configuration()
adata = utv.run_model(adata,label=data_cluster,config_file=velo_config)

adata.write_h5ad(output_path +"unitvelo.h5ad")