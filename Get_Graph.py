import os
import scanpy as sc
import numpy as np
import scvelo as scv
import pandas as pd
import math
import argparse
import matplotlib.pyplot as plt
from tools.caculate_eval import evaluate ,get_global
from config import output_path,Methods_name,edges,data_cluster,embed,figures_path,evals_path,input_path,n_job,velocity_layer,celltype_palette

eval_all = pd.DataFrame(columns=['Method','GDC','CBDir','ICCoh'])
method_to_draw = []
for method in Methods_name.keys():
   if os.path.exists(output_path + method + '.h5ad'):
        method_to_draw.append(method)

print(f'A total of {str(len(method_to_draw))} .h5ad files of RNA velocity results were detected')
        
nrow = 4
ncol = math.ceil(len(method_to_draw)/nrow)
fig, axs = plt.subplots(nrow, ncol, figsize=(5 * ncol, 5 * nrow))
f = 0

for method in method_to_draw:
    try :
        r = int(f/ncol)
        c = f % ncol
        ax = axs[r, c]

        adata = sc.read_h5ad(output_path + method + '.h5ad')

        scv.pp.neighbors(adata)

        all_nan_columns = np.where(np.isnan(adata.layers[velocity_layer]).all(axis=0))[0].tolist()
        if all_nan_columns:
            adata._inplace_subset_var(np.delete(adata.var_names, all_nan_columns))

        scv.tl.velocity_graph(adata, basis=embed,vkey=velocity_layer,n_jobs=n_job)

        # Visualization 
        # merge
        scv.pl.velocity_embedding_stream(adata,basis=embed,vkey=velocity_layer,color=data_cluster,
                                        title=Methods_name[method],palette=celltype_palette,
                                        ax = ax,frameon=False,show=False
                                        )

        # save each figure
        scv.pl.velocity_embedding_stream(adata,basis=embed,vkey=velocity_layer,color=data_cluster,
                                        title=Methods_name[method],palette=celltype_palette,
                                        save=figures_path + Methods_name[method] + '.svg'
                                        )
                                        
        ax.set_title(ax.get_title(), fontsize=16) 
        f = f+1

        # Evaluation
        cell_number = adata.obs.shape[0]
        cluster_edges, k_cluster, k_velocity, x_emb=edges,data_cluster,velocity_layer,'X_'+embed
        list, eval_list  = evaluate(adata, cluster_edges, k_cluster, k_velocity, x_emb=x_emb, verbose=True)
        in_cluster = eval_list["ic_coh"]
        cross_boundary =eval_list["crs_bdr_crc"]
        eval_global = get_global(in_cluster,cross_boundary,cell_number)

        new_row = {'Method':Methods_name[method],'GDC':eval_global,'CBDir':list[0],'ICCoh':list[1]}
        eval_all = eval_all.append(new_row, ignore_index=True)
    
    except Exception as e:
        print("run {}  get an error:{}".format(method,e))

blank_fig = nrow * ncol - len(method_to_draw)
if blank_fig > 0:
    for i in range(0,blank_fig):
        r = int(f/ncol)
        c = f % ncol
        ax = axs[r, c]
        ax.set_visible(False)
        ax.set_facecolor('none')
        f = f+1

plt.savefig(figures_path + 'Merge.svg',format='svg')
eval_all.to_csv(evals_path + 'Eval.csv')