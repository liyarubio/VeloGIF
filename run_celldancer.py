import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import celldancer as cd
import celldancer.utilities as cdutil
import os
from celldancer.utilities import export_velocity_to_dynamo
from config import input_path , output_path  ,seed  , data_cluster ,embed , n_job

np.random.seed(seed)

input_file = input_path

print(f"Reading file: {input_file}")

adata = sc.read_h5ad(input_file)
print(adata)
emb = 'X_' + embed

cdutil.adata_to_df_with_embed(adata,
                              us_para=['Mu','Ms'],
                              # cell_type_para='cell_type',
                              cell_type_para=data_cluster,
                              embed_para=emb,
                              save_path='res/celldancer_input.csv'
                             )


df = pd.read_csv('res/celldancer_input.csv')
loss_df, cellDancer_df=cd.velocity(df,n_jobs=n_job,
                                   speed_up = False)


adata_cd = export_velocity_to_dynamo(cellDancer_df,adata)
print(adata_cd)

adata_cd.layers["velocity"] = adata_cd.layers["velocity_S"].toarray()

adata_cd.write_h5ad(output_path + "/celldancer.h5ad")