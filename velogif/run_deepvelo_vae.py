import os
import scvelo as scv
import scanpy as sc
import numpy as np
import sklearn
import pandas as pd
import seaborn as sns
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import seaborn as sns
from config import input_path , output_path  ,seed ,device,n_job
from tools import vae_from_deepvelo_SA as dv

np.random.seed(seed)
tf.random.set_seed(seed)

input_file = input_path
adata = sc.read_h5ad(input_file)

# run scvelo first, get initial velocity
scv.tl.velocity(adata,n_jobs=n_job)

adata.var["velocity_genes"] = True
adata_raw = adata.copy()

X = np.tile(adata_raw.X.A[:, adata.var["velocity_genes"]], (5, 1))
Y = np.tile(adata.layers["velocity"][:, adata.var["velocity_genes"]], (5, 1))
noise_sigma = (adata_raw.X.A.std()/70)**2
X[adata_raw.shape[0]:, :] += \
    np.random.normal(0, noise_sigma, X[adata_raw.shape[0]:, :].shape)

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, 
                                                    Y, 
                                                    test_size=0.1, 
                                                    random_state=seed # set 2024
                                                    )

encoder = dv.create_encoder(X.shape[1])
decoder = dv.create_decoder(X.shape[1])

autoencoder = dv.VAE(encoder, decoder)
opt = keras.optimizers.Adam(learning_rate = 0.00005) # default: learning_rate 0.001; in deepvelo tutorial 0.00005
autoencoder.compile(optimizer=opt)

es = keras.callbacks.EarlyStopping(monitor='val_loss', patience=3)  # as tutorial figure2 set
autoencoder.fit(X_train, y_train,
        epochs=100, # as tutorial figure2 set
        batch_size=2, # as tutorial figure2 set
        shuffle=True, # as tutorial figure2 set
        validation_data=(X_test, y_test),
        callbacks=[es])

X = adata_raw.X.A[:, adata.var["velocity_genes"]]
velocity_deepvelo = autoencoder.predict(X)
print(velocity_deepvelo.shape)
adata.layers['velocity'] = velocity_deepvelo

adata.write_h5ad(output_path + "deepvelo_vae.h5ad")