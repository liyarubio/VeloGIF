.. raw:: html

    <style>
    /* 表格容器样式 */
    .wy-table-responsive {
        overflow-x: auto;
        position: relative;
    }

    /* 表格基础样式 */
    .wy-table-responsive table {
        margin-bottom: 24px;
        width: 100%;
        border-collapse: separate;
        border-spacing: 0;
    }

    /* 首列固定样式 */
    .wy-table-responsive table th:first-child,
    .wy-table-responsive table td:first-child {
        position: sticky;
        left: 0;
        background-color: white;
        z-index: 2;
        border-right: 2px solid #e1e4e5;
    }

    /* 表头首列样式 */
    .wy-table-responsive table th:first-child {
        background-color: #f3f6f6;
        z-index: 3;
    }

    /* 首列阴影效果 */
    .wy-table-responsive table th:first-child::after,
    .wy-table-responsive table td:first-child::after {
        content: "";
        position: absolute;
        top: 0;
        right: -2px;
        bottom: 0;
        width: 4px;
        box-shadow: inset -4px 0 4px -2px rgba(0, 0, 0, 0.1);
    }
    </style>

Methods Introduction
==============================



.. list-table::
   :header-rows: 1

   * - Method
     - Year
     - Publication
     - Based splicing
     - Key Words / Major innovation
     - Core Algorithm
     - Website
     - GPU acceleration
     - Deep Learning Framework
     - Need origin RNA velocity
   * - `velocyto <http://velocyto.org/velocyto.py/index.html>`_
     - 2018
     - Nature
     - √
     - steady state
     - extreme expression quantiles regression
     - √
     - 
     - 
     - 
   * - `scVelo <https://scvelo.readthedocs.io/en/stable/>`_
     - 2020
     - Nature Biotechnology
     - √
     - dynamic, EM, latent time
     - expectation-maximization (EM)
     - √
     - 
     - 
     - 
   * - `VeloAE <https://github.com/qiaochen/VeloAE>`_
     - 2021
     - Proceedings of the National Academy of Sciences
     - √
     - autoencoder, GCN, attention
     - autoencoder (AE)
     - 
     - √
     - Torch
     - optional: for velocity constraint
   * - `Dynamo <https://dynamo-release.readthedocs.io/en/latest/index.html>`_
     - 2022
     - Cell
     - splicing / new + exp / new + u + s
     - GMM, transcriptomic vector fields
     - generalized method of moments (GMM), negative binomial (NB) method or curve fitting.
     - √
     - 
     - 
     - 
   * - `VeloVAE <https://github.com/welch-lab/VeloVAE>`_
     - 2022
     - International Conference on Machine Learning
     - √
     - variational mixture of ODEs, latent state
     - variational autoencoder (VAE)
     - 
     - √
     - Torch
     - 
   * - `MultiVelo <https://github.com/welch-lab/MultiVelo/>`_
     - 2022
     - Nature Biotechnology
     - c + u + s
     - chromatin accessibility
     - EM
     - √
     - √
     - Torch
     - 
   * - `UniTVelo <https://github.com/StatBiomed/UniTVelo>`_
     - 2022
     - Nature Communications
     - √
     - temporally unified, top-down
     - radial basis function (RBF)
     - 
     - √
     - TensorFlow
     - 
   * - `DeepVelo (VAE-based) <https://github.com/aicb-ZhangLabs/deepvelo>`_
     - 2022
     - Science Advances
     - √
     - neural differential equation
     - VAE
     - 
     - √
     - TensorFlow
     - for training
   * - `cellDancer <https://github.com/GuangyuWangLab2021/cellDancer>`_
     - 2023
     - Nature Biotechnology
     - √
     - single-cell resolution inference of velocity kinetics
     - deep neural networks (DNN)
     - √
     - 
     - Torch
     - 
   * - `PhyloVelo <https://github.com/kunwang34/PhyloVelo>`_
     - 2023
     - Nature Biotechnology
     - lineage + exp
     - phylogenetic time, monotonically expressed genes
     - diffusion process of gene expressions
     - √
     - 
     - 
     - 
   * - `veloVI <https://velovi.readthedocs.io/en/latest/>`_
     - 2023
     - Nature Methods
     - √
     - variational inference, velocity uncertainty
     - VAE
     - √
     - √
     - Torch
     - 
   * - `LatentVelo <https://github.com/Spencerfar/LatentVelo>`_
     - 2023
     - Cell Reports Methods
     - √
     - neural ODE, batch, variational autoencoder
     - VAE
     - 
     - √
     - Torch
     - 
   * - `DeepVelo (GCN-based) <https://github.com/bowang-lab/DeepVelo>`_
     - 2024
     - Genome Biology
     - √
     - GCN, time-dependent kinetics and multiple lineages
     - convolutional networks (GCN)
     - 
     - √
     - Torch
     - 
   * - `TFvelo <https://github.com/xiaoyeye/TFvelo>`_
     - 2024
     - Nature Communications
     - regulatory + exp
     - gene regulatory
     - Least Absolute Shrinkage and Selection Operator (LASSO) regression, generalized EM
     - 
     - 
     - 
     - 
   * - `STT <https://github.com/cliffzhou92/STT>`_
     - 2024
     - Nat Methods
     - √
     - multiscale dynamical model, spatial transcriptomics
     - stochastic differential equation (SDE), generalized Perron cluster cluster analysis (GPCCA)
     - 
     - 
     - 
     - 

.. raw:: html

    <script>
        document.addEventListener('DOMContentLoaded', function() {
            // 获取所有链接
            var links = document.getElementsByTagName('a');
            // 为每个链接添加target="_blank"属性
            for(var i = 0; i < links.length; i++) {
                links[i].setAttribute('target', '_blank');
                links[i].setAttribute('rel', 'noopener noreferrer');
            }
        });
    </script>

    .. raw:: html

    <script>
        document.addEventListener('DOMContentLoaded', function() {
            var links = document.getElementsByTagName('a');
            for(var i = 0; i < links.length; i++) {
                links[i].setAttribute('target', '_blank');
                links[i].setAttribute('rel', 'noopener noreferrer');
            }
        });
    </script>
