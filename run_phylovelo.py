import phylovelo as pv
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import re
import numpy as np
import pickle
import scvelo as scv
from config import output_path,pv_adata_path,pv_pkl_path

with open(pv_pkl_path, 'rb') as f:
    sd = pickle.load(f)

exp = sd.x_normed
depth = sd.depths
sel_xdr = sd.Xdr

pv.velocity_inference(sd,sd.depths, cutoff=0.95, target='x_normed')
pv.velocity_embedding(sd, target='x_normed')

adata = sc.read_h5ad(pv_adata_path)

v = pd.DataFrame(sd.velocity,columns=['i'])
vv = pd.concat([v['i']] * exp.shape[0], axis=1)
vv = vv.T

scv.pp.moments(adata) # for visualization
adata.layers['velocity'] = vv.values
adata.write_h5ad(output_path+'phylovelo.h5ad')