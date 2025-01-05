import velocyto as vcy
import scanpy as sc
import scvelo as scv
import numpy as np
import os
from config import input_path,output_path,seed

np.random.seed(seed)

input_file = input_path if input_path.endswith('.h5ad') else os.path.join(input_path, 'adata.h5ad')
print(f"Reading file: {input_file}")

adata = sc.read_h5ad(input_file)


vlm = scv.utils.convert_to_loom(adata)

vlm.fit_gammas(fit_offset=False)
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift()
vlm.extrapolate_cell_at_t()
velocity_s = np.transpose(vlm.velocity)

adata.layers['velocity'] = velocity_s
adata.write_h5ad(output_path+'velocyto.h5ad')