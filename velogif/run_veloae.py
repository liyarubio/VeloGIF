import pickle, os
import numpy as np
import scvelo as scv
import pandas as pd
import anndata as ad
import torch
from veloproj import *
import scanpy as sc
from config import input_path , output_path,seed ,device, embed,n_job

np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
np.random.seed(seed)

parser = get_parser()
args = parser.parse_args(args=['--device', device])

torch.backends.cudnn.deterministic = True
device = torch.device(args.device if args.device.startswith('cuda') and torch.cuda.is_available() else "cpu")

input_file = input_path
adata = sc.read_h5ad(input_file)

# run scvelo first, get initial velocity
scv.tl.velocity(adata,n_jobs=n_job)

def main_AE(args, adata):
    spliced = adata.layers['Ms']
    unspliced = adata.layers['Mu']
    tensor_s = torch.FloatTensor(spliced).to(device)
    tensor_u = torch.FloatTensor(unspliced).to(device)
    tensor_x = torch.FloatTensor(adata.X.toarray()).to(device)
    tensor_v = torch.FloatTensor(adata.layers['velocity']).to(device)

    model = init_model(adata, args, device)

    inputs = [tensor_s, tensor_u]
    xyids = [0, 1]
    if args.use_x:
        inputs.append(tensor_x)

    model = fit_model(args, adata, model, inputs, tensor_v, xyids, device)
    return tensor_s, tensor_u, tensor_x  

tensor_s, tensor_u, tensor_x = main_AE(args, adata)

model = init_model(adata, args, device)
model.load_state_dict(torch.load(args.model_name))
model = model.to(device)
model.eval()
with torch.no_grad():
    x = model.encoder(tensor_x)
    s = model.encoder(tensor_s)
    u = model.encoder(tensor_u)
    
    v = estimate_ld_velocity(s, u, device=device, perc=[5, 95], 
                                norm=args.use_norm, fit_offset=args.fit_offset_pred, 
                                use_offset=args.use_offset_pred).cpu().numpy()
    x = x.cpu().numpy()
    s = s.cpu().numpy()
    u = u.cpu().numpy()


emb = 'X_' +embed
adata_new = new_adata(adata, 
                  x, s, u, v, 
                  n_nb_newadata=args.n_nb_newadata, # default 30 
                  X_emb_key=emb)
adata_new.layers['velocity'] = adata_new.layers['new_velocity']
adata_new.write_h5ad(output_path + "veloae.h5ad")