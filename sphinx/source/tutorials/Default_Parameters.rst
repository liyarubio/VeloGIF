.. role:: raw-html(raw)
   :format: html

Default Parameters
==================

velocyto
--------
The pipeline was employed as per the guideline provided on the velocyto website via
:raw-html:`<a href="http://velocyto.org/velocyto.py/tutorial/analysis.html" target="_blank">http://velocyto.org/velocyto.py/tutorial/analysis.html</a>`.

scVelo
------
We followed the guideline provided on the scVelo website via
:raw-html:`<a href="https://scvelo.readthedocs.io/en/stable/VelocityBasics.html" target="_blank">https://scvelo.readthedocs.io/en/stable/VelocityBasics.html</a>`
when running "stochastic" mode.

Additionally, we followed the guideline via
:raw-html:`<a href="https://scvelo.readthedocs.io/en/stable/DynamicalModeling.html" target="_blank">https://scvelo.readthedocs.io/en/stable/DynamicalModeling.html</a>`
when running "dynamical" mode. In this mode, ``scv.tl.recover_dynamics()`` was required to run first,
and the parameter ``mode='dynamical'`` was set when running ``scv.tl.velocity()``.

VeloAE
------
We followed the pipeline provided on the VeloAE GitHub repository:
:raw-html:`<a href="https://github.com/qiaochen/VeloAE/blob/main/notebooks/dentategyrus/model-dentategyrus_gat.ipynb" target="_blank">https://github.com/qiaochen/VeloAE/blob/main/notebooks/dentategyrus/model-dentategyrus_gat.ipynb</a>`.

The ``--g-rep-dim`` and ``--k-dim`` parameters were set to 10 in simulation datasets from dyngen
due to the fewer number of genes in these datasets. In other datasets, these two parameters were
set to 100, which is the default value.

Dynamo
------
We followed the guideline on the dynamo website via
:raw-html:`<a href="https://dynamo-release.readthedocs.io/en/latest/notebooks/zebrafish.html" target="_blank">https://dynamo-release.readthedocs.io/en/latest/notebooks/zebrafish.html</a>`
for conventional scRNA-seq datasets. The default "stochastic model" was used in most datasets,
while the alternative model, "deterministic model," was used when the "stochastic model" did not converge.

Additionally, we followed the guideline via
:raw-html:`<a href="https://dynamo-release.readthedocs.io/en/latest/notebooks/tutorial_hsc_velocity.html" target="_blank">https://dynamo-release.readthedocs.io/en/latest/notebooks/tutorial_hsc_velocity.html</a>`
for metabolic labeling datasets.

VeloVAE
-------
We followed the tutorial on the VeloVAE GitHub repository:
:raw-html:`<a href="https://github.com/welch-lab/VeloVAE/blob/master/notebooks/velovae_example.ipynb" target="_blank">https://github.com/welch-lab/VeloVAE/blob/master/notebooks/velovae_example.ipynb</a>`.

The default model, "continuous veloVAE" model, was employed, and ``tmax`` was set to 20. As suggested
in the tutorial, the ``dim_z`` was set to the number of cell types in datasets with known cell types,
while it was set to 5 in datasets with unknown cell types.

MultiVelo
---------
We followed the tutorial on the MultiVelo GitHub repository:
:raw-html:`<a href="https://github.com/welch-lab/MultiVelo/blob/main/Examples/MultiVelo_Demo.ipynb" target="_blank">https://github.com/welch-lab/MultiVelo/blob/main/Examples/MultiVelo_Demo.ipynb</a>`.

UniTVelo
--------
We followed the tutorial on the UniTVelo GitHub repository:
:raw-html:`<a href="https://github.com/StatBiomed/UniTVelo/blob/main/notebooks/Figure2_ErythroidMouse.ipynb" target="_blank">https://github.com/StatBiomed/UniTVelo/blob/main/notebooks/Figure2_ErythroidMouse.ipynb</a>`.

The "Unified-time mode" was employed by keeping the ``FIT_OPTION`` to 1, which is the default parameter.

DeepVelo (VAE-based)
--------------------
We followed the pipeline on the DeepVelo GitHub repository:
:raw-html:`<a href="https://github.com/aicb-ZhangLabs/deepvelo/blob/main/code/Figure2.ipynb" target="_blank">https://github.com/aicb-ZhangLabs/deepvelo/blob/main/code/Figure2.ipynb</a>`.

We used the VAE framework rather than the normal AE framework, and both frameworks were available
in the tutorial. The parameters ``epochs=100``, ``batch_size=2``, and ``shuffle=True`` were used,
consistent with the tutorial. The RNA velocity matrix calculated by scVelo's stochastic mode was
used to train the model.

PhyloVelo
---------
We followed the tutorial on the PhyloVelo website:
:raw-html:`<a href="https://phylovelo.readthedocs.io/en/latest/notebook/KPTracer-3726-t1.html" target="_blank">https://phylovelo.readthedocs.io/en/latest/notebook/KPTracer-3726-t1.html</a>`.

cellDancer
----------
We followed the guideline on the cellDancer website:
:raw-html:`<a href="https://guangyuwanglab2021.github.io/cellDancer_website/notebooks/case_study_gastrulation.html" target="_blank">https://guangyuwanglab2021.github.io/cellDancer_website/notebooks/case_study_gastrulation.html</a>`.

The parameter ``speed_up`` was set to ``False`` to use all cells to train the model, allowing for
subsequent quantitative comparisons.

veloVI
------
We followed the tutorial on the veloVI website via
:raw-html:`<a href="https://velovi.readthedocs.io/en/latest/tutorial.html" target="_blank">https://velovi.readthedocs.io/en/latest/tutorial.html</a>`.

LatentVelo
----------
We followed the tutorial on the LatentVelo GitHub repository:
:raw-html:`<a href="https://github.com/Spencerfar/LatentVelo/blob/main/paper_notebooks/Dentate%20gyrus.ipynb" target="_blank">https://github.com/Spencerfar/LatentVelo/blob/main/paper_notebooks/Dentate%20gyrus.ipynb</a>`.

The "LatentVelo VAE model" was used, and the parameters were kept at their default values.

DeepVelo (GCN-based)
--------------------
We followed the tutorial on the DeepVelo GitHub repository:
:raw-html:`<a href="https://github.com/bowang-lab/DeepVelo/blob/main/examples/figure2.ipynb" target="_blank">https://github.com/bowang-lab/DeepVelo/blob/main/examples/figure2.ipynb</a>`.

TFvelo
------
We followed the tutorial on the TFvelo GitHub repository via
:raw-html:`<a href="https://github.com/xiaoyeye/TFvelo/blob/main/TFvelo_demo.ipynb" target="_blank">https://github.com/xiaoyeye/TFvelo/blob/main/TFvelo_demo.ipynb</a>`.

STT
---
We followed the tutorial on the STT GitHub repository:
:raw-html:`<a href="https://github.com/cliffzhou92/STT/blob/release/example_notebooks/example-mouse_eryth.ipynb" target="_blank">https://github.com/cliffzhou92/STT/blob/release/example_notebooks/example-mouse_eryth.ipynb</a>`.

The parameter ``n_states`` was set to the number of cell types in datasets with known cell types,
while it was set to 2 in datasets with unknown cell types.