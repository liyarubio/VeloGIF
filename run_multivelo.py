import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt
import re

from config import  mv_rna_path,mv_atac_path,output_path,n_job,seed

np.random.seed(seed)

adata_rna = sc.read_h5ad(mv_rna_path)
adata_atac = sc.read_h5ad(mv_atac_path)

adata_result = mv.recover_dynamics_chrom(adata_rna, 
                                         adata_atac,
                                         n_jobs=n_job)

adata_result.layers['velocity'] = adata_result.layers['velo_s']
adata_result.write(output_path  + "multivelo.h5ad")