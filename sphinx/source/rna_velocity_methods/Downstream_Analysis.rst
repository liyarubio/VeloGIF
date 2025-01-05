Downstream Analysis
============================================

.. note:: All "√" symbols in the table are clickable links that will redirect you to the corresponding documentation, tutorials, or example notebooks for each tool.

.. raw:: html

    <style>
        /* 表格容器样式 */
    .wy-table-responsive {
        margin: 20px 0;
        overflow: visible;  /* 改为visible以避免出现双滚动条 */
        position: relative;
    }

    /* 表格主体样式 */
    .custom-table-container {
        max-height: 80vh;
        overflow-y: auto;
        overflow-x: auto;
        position: relative;
        margin-bottom: 24px;
        border: 1px solid #ddd;
        border-radius: 4px;
    }

    /* 美化滚动条 */
    .custom-table-container::-webkit-scrollbar {
        width: 10px;
        height: 10px;
    }

    .custom-table-container::-webkit-scrollbar-track {
        background: #f1f1f1;
        border-radius: 5px;
    }

    .custom-table-container::-webkit-scrollbar-thumb {
        background: #888;
        border-radius: 5px;
    }

    .custom-table-container::-webkit-scrollbar-thumb:hover {
        background: #555;
    }

    /* 保持其他原有样式 */
    .custom-table {
        width: 100%;
        border-collapse: collapse;
        position: relative;
    }

    .custom-table th, .custom-table td {
        text-align: center;
        vertical-align: middle;
        padding: 8px;
        border: 1px solid #ddd;
        min-width: 100px;
    }

    /* 固定表头和第一列的样式保持不变 */
    .custom-table thead tr th {
        position: sticky;
        top: 0;
        background: white;
        z-index: 2;
    }

    .custom-table td:first-of-type {
        position: sticky;
        left: 0;
        background: white;
        z-index: 1;
        border: 1px solid #ddd;
        min-width: 150px;
        max-width: 200px;
        white-space: normal;
        box-shadow: 2px 0 5px rgba(0,0,0,0.1);
    }

    </style>
    <div class="custom-table-container">
    <table class="custom-table">
        <thead>
            <tr>
                <th></th>
                <th></th>
                <th>velocyto</th>
                <th>scVelo</th>
                <th>VeloAE</th>
                <th>Dynamo</th>
                <th>VeloVAE</th>
                <th>MultiVelo</th>
                <th>UniTVelo</th>
                <th>DeepVelo (VAE)</th>
                <th>cellDancer</th>
                <th>PhyloVelo</th>
                <th>veloVI</th>
                <th>LatentVelo</th>
                <th>DeepVelo (GCN)</th>
                <th>TFvelo</th>
                <th>STT</th>
            </tr>
        </thead>
        <tbody>
            <!-- Basic Function 部分 (第2至11行) -->
            <tr>
                <th rowspan="10">Basic Function</th>
                <td>compute u & s exp matrix</td>
                <td><a href="http://velocyto.org/velocyto.py/tutorial/cli.html#running-velocyto">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
            </tr>
            <tr>
                <td>filter and normalize</td>
                <td><a href="http://velocyto.org/velocyto.py/tutorial/analysis.html#start-a-new-analysis-preliminary-filtering">√</a></td>
                <td><a href="https://scvelo.readthedocs.io/en/stable/VelocityBasics.html#Preprocess-the-Data">√</a></td>
                <td></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html#preprocess-data">√</a></td>
                <td><a href="https://github.com/welch-lab/VeloVAE/blob/master/notebooks/velovae_example.ipynb">√</a></td>
                <td><a href="https://github.com/welch-lab/MultiVelo/blob/main/Examples/MultiVelo_Demo.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://phylovelo.readthedocs.io/en/latest/notebook/C.elegans_demo.html">√</a></td>
                <td></td>
                <td><a href="https://github.com/Spencerfar/LatentVelo/blob/3cdf02946e53411aad68700ea0bba1254c48ec54/paper_notebooks/Dentate%20gyrus.ipynb">√</a></td>
                <td></td>
                <td><a href="https://github.com/xiaoyeye/TFvelo/blob/main/TFvelo_run_demo.py">√</a></td>
                <td></td>
            </tr>
            <tr>
                <td>compute cell x gene velocity matrix</td>
                <td><a href="http://velocyto.org/velocyto.py/tutorial/analysis.html#gamma-fit-and-extrapolation">√</a></td>
                <td><a href="https://scvelo.readthedocs.io/en/stable/VelocityBasics.html#Estimate-RNA-velocity">√</a></td>
                <td></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html#velocity-vectors">√</a></td>
                <td><a href="https://github.com/welch-lab/VeloVAE/blob/master/notebooks/velovae_example.ipynb">√</a></td>
                <td><a href="https://github.com/welch-lab/MultiVelo/blob/main/Examples/MultiVelo_Demo.ipynb">√</a></td>
                <td><a href="https://github.com/StatBiomed/UniTVelo/blob/main/notebooks/Figure2_ErythroidMouse.ipynb">√</a></td>
                <td><a href="https://github.com/aicb-ZhangLabs/deepvelo/blob/main/code/Figure2.ipynb">√</a></td>
                <td><a href="https://guangyuwanglab2021.github.io/cellDancer_website/notebooks/case_study_gastrulation.html">√</a></td>
                <td><a href="https://phylovelo.readthedocs.io/en/latest/notebook/C.elegans_demo.html">√</a></td>
                <td><a href="https://velovi.readthedocs.io/en/latest/tutorial.html">√</a></td>
                <td><a href="https://github.com/Spencerfar/LatentVelo/blob/3cdf02946e53411aad68700ea0bba1254c48ec54/paper_notebooks/Dentate%20gyrus.ipynb">√</a></td>
                <td><a href="https://github.com/bowang-lab/DeepVelo/blob/main/examples/figure2.ipynb">√</a></td>
                <td><a href="https://github.com/xiaoyeye/TFvelo/blob/main/TFvelo_run_demo.py">√</a></td>
                <td><a href="https://github.com/cliffzhou92/STT/blob/release/example_notebooks/example-bone-marrow.ipynb">√</a></td>
            </tr>
            <tr>
                <td>compute unsplice velocity</td>
                <td></td>
                <td><a href="https://scvelo.readthedocs.io/en/stable/DynamicalModeling.html">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://guangyuwanglab2021.github.io/cellDancer_website/notebooks/case_study_gastrulation.html">√</a></td>
                <td></td>
                <td><a href="https://velovi.readthedocs.io/en/latest/tutorial.html">√</a></td>
                <td><a href="https://github.com/Spencerfar/LatentVelo/blob/3cdf02946e53411aad68700ea0bba1254c48ec54/paper_notebooks/Dentate%20gyrus.ipynb">√</a></td>
                <td><a href="https://github.com/bowang-lab/DeepVelo/blob/main/examples/figure2.ipynb">√</a></td>
                <td></td>
                <td><a href="https://github.com/cliffzhou92/STT/blob/release/example_notebooks/example-bone-marrow.ipynb">√</a></td>
            </tr>
            <tr>
                <td>simulate data</td>
                <td></td>
                <td><a href="https://github.com/theislab/scvelo_notebooks/blob/master/perspectives/Perspectives_parameters.ipynb">√</a></td>
                <td></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/API.html#simulation-sim">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://guangyuwanglab2021.github.io/cellDancer_website/celldancer.simulation.simulate.html">√</a></td>
                <td><a href="https://phylovelo.readthedocs.io/en/latest/autoapi/phylovelo/sim_utils/index.html">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/xiaoyeye/TFvelo/blob/main/simulation/TFvelo_synthetic_demo.py">√</a></td>
                <td><a href="https://github.com/cliffzhou92/STT/blob/1a1e61e4848002bf127670a198844cb11fd067a5/data/emt_sim/functions_emt.py#L78">√</a></td>
            </tr>
            <tr>
                <td>compute transition probability</td>
                <td><a href="http://velocyto.org/velocyto.py/tutorial/analysis.html#projection-of-velocity-onto-embeddings">√</a></td>
                <td><a href="https://scvelo.readthedocs.io/en/stable/scvelo.utils.get_transition_matrix.html#scvelo.utils.get_transition_matrix">√</a></td>
                <td></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html#velocity-vectors">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://guangyuwanglab2021.github.io/cellDancer_website/notebooks/case_study_gastrulation.html">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
            </tr>
            <tr>
                <td>compute low-dim velocity</td>
                <td><a href="http://velocyto.org/velocyto.py/tutorial/analysis.html#projection-of-velocity-onto-embeddings">√</a></td>
                <td><a href="https://scvelo.readthedocs.io/en/stable/VelocityBasics.html#Estimate-RNA-velocity">√</a></td>
                <td><a href="https://github.com/qiaochen/VeloAE/tree/main/notebooks">√</a></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html#velocity-vectors">√</a></td>
                <td><a href="https://github.com/welch-lab/VeloVAE/blob/master/notebooks/velovae_example.ipynb">√</a></td>
                <td><a href="https://github.com/welch-lab/MultiVelo/blob/main/Examples/MultiVelo_Demo.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td><a href="https://guangyuwanglab2021.github.io/cellDancer_website/notebooks/case_study_gastrulation.html">√</a></td>
                <td><a href="https://phylovelo.readthedocs.io/en/latest/notebook/C.elegans_demo.html">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/xiaoyeye/TFvelo/blob/main/TFvelo_demo.ipynb">√</a></td>
                <td></td>
            </tr>
            <tr>
                <td>compute pseudotime / latent time</td>
                <td><a href="https://github.com/velocyto-team/velocyto-notebooks/blob/master/python/hgForebrainGlutamatergic.ipynb">√</a></td>
                <td><a href="https://scvelo.readthedocs.io/en/stable/DynamicalModeling.html#Latent-time">√</a></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/welch-lab/VeloVAE/blob/master/notebooks/velovae_example.ipynb">√</a></td>
                <td><a href="https://github.com/welch-lab/MultiVelo/blob/main/Examples/MultiVelo_Demo.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td><a href="https://guangyuwanglab2021.github.io/cellDancer_website/notebooks/case_study_gastrulation.html">√</a></td>
                <td><a href="https://phylovelo.readthedocs.io/en/latest/notebook/getting_start.html#5.-PhyloVelo-pseudotime-analysis">√</a></td>
                <td><a href="https://velovi.readthedocs.io/en/latest/tutorial.html">√</a></td>
                <td></td>
                <td><a href="https://github.com/bowang-lab/DeepVelo/blob/main/examples/figure3.ipynb">√</a></td>
                <td><a href="https://github.com/xiaoyeye/TFvelo/blob/main/TFvelo_demo.ipynb">√</a></td>
                <td></td>
            </tr>
            <tr>
                <td>identify key genes/ pathway</td>
                <td></td>
                <td><a href="https://scvelo.readthedocs.io/en/stable/scvelo.tl.rank_velocity_genes.html#scvelo.tl.rank_velocity_genes">√</a></td>
                <td></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/_autosummary/dynamo.vf.rank_genes.html#dynamo.vf.rank_genes">√</a></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/StatBiomed/UniTVelo/blob/main/notebooks/Figure2_ErythroidMouse.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/YosefLab/velovi_reproducibility/blob/main/notebooks/case_studies/dentategyrus.ipynb">√</a></td>
                <td><a href="https://github.com/Spencerfar/LatentVelo/blob/main/paper_notebooks/Gastrulation.ipynb">√</a></td>
                <td><a href="https://github.com/bowang-lab/DeepVelo/blob/main/examples/incorporate_cellrank.ipynb">√</a></td>
                <td></td>
                <td><a href="https://github.com/cliffzhou92/STT/blob/release/example_notebooks/example-emt.ipynb">√</a></td>
            </tr>
            <tr>
                <td>evaluate velocity</td>
                <td></td>
                <td><a href="https://scvelo.readthedocs.io/en/stable/VelocityBasics.html#Speed-and-coherence">√</a></td>
                <td><a href="https://github.com/qiaochen/VeloAE/tree/main/notebooks">√</a></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/_autosummary/dynamo.tl.cell_wise_confidence.html#dynamo.tl.cell_wise_confidence">√</a></td>
                <td><a href="https://github.com/welch-lab/VeloVAE/blob/master/notebooks/velovae_example.ipynb">√</a></td>
                <td></td>
                <td><a href="https://github.com/StatBiomed/UniTVelo/blob/main/notebooks/Figure2_ErythroidMouse.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://velovi.readthedocs.io/en/latest/tutorial.html">√</a></td>
                <td><a href="https://github.com/Spencerfar/LatentVelo/blob/main/paper_notebooks/Hindbrain.ipynb">√</a></td>
                <td><a href="https://github.com/bowang-lab/DeepVelo/blob/main/examples/figure2.ipynb">√</a></td>
                <td></td>
                <td></td>
            </tr>

            <!-- Extended Function 部分 (第12至17行) -->
            <tr>
                <th rowspan="6">Extended Function</th>
                <td>find root / end points</td>
                <td><a href="https://github.com/velocyto-team/velocyto-notebooks/blob/master/python/DentateGyrus.ipynb">√</a></td>
                <td><a href="https://scvelo.readthedocs.io/en/stable/scvelo.tl.terminal_states.html#scvelo.tl.terminal_states">√</a></td>
                <td></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/_autosummary/dynamo.pd.andecestor.html#dynamo.pd.andecestor">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
            </tr>
            <tr>
                <td>predict in silico perturbation </td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/notebooks/perturbation_tutorial/perturbation_tutorial.html">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/aicb-ZhangLabs/deepvelo/blob/main/code/Figure2.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                
            </tr>
            <tr>
                <td>predict transition path / trajectory</td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/notebooks/lap_tutorial/lap_tutorial.html">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/aicb-ZhangLabs/deepvelo/blob/main/code/Figure2.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/Spencerfar/LatentVelo/blob/main/paper_notebooks/Gastrulation.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/cliffzhou92/STT/blob/1a1e61e4848002bf127670a198844cb11fd067a5/README.md?plain=1#L34">√</a></td>
            </tr>
            <tr>
                <td>reconstruct vector field</td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html#vector-field-reconstruction">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
            </tr>
            <tr>
                <td>correct batch effect</td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/Spencerfar/LatentVelo/blob/main/paper_notebooks/Bifurcation%20batch%20correction.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td></td>
            </tr>
            <tr>
                <td>epigenome–transcriptome interactions</td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td><a href="https://github.com/welch-lab/MultiVelo/blob/main/Examples/MultiVelo_Demo.ipynb">√</a></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
                <td></td>
            </tr>

            <!-- Mode 行 (第18行) -->
            

            <tr>
                <th rowspan="1">Mode</th>
                <td></td>
                <td></td>
                <td>
                    <a href="https://scvelo.readthedocs.io/en/stable/getting_started.html#velocity-tools">deterministic</a>
                    <a href="https://scvelo.readthedocs.io/en/stable/getting_started.html#velocity-tools">stochastic</a>
                    <a href="https://scvelo.readthedocs.io/en/stable/DynamicalModeling.html">dynamic</a>
                </td>
                <td></td>
                <td>
                    <a href="https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html#learn-dynamics">deterministic</a>
                    <a href="https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html#learn-dynamics">stochastic</a>
                </td>
                <td>
                    <a href="https://github.com/welch-lab/VeloVAE/tree/master/notebooks/paper_notebooks/continuous">continuous veloVAE</a>
                    <a href="https://github.com/welch-lab/VeloVAE/tree/master/notebooks/paper_notebooks/continuous">continuous Full VB</a>
                    <a href="https://github.com/welch-lab/VeloVAE/tree/master/notebooks/paper_notebooks/discrete">discrete VAE</a>
                    <a href="https://github.com/welch-lab/VeloVAE/tree/master/notebooks/paper_notebooks/discrete">discrete FullVB</a>
                    <a href="https://github.com/welch-lab/VeloVAE/tree/master/notebooks/paper_notebooks/branching_ode">branching ODE</a>
                </td>
                <td></td>
                <td>
                    <a href="https://github.com/StatBiomed/UniTVelo/blob/main/notebooks/Figure2_ErythroidMouse.ipynb">Unified-time</a>
                    <a href="https://github.com/StatBiomed/UniTVelo/blob/main/notebooks/SuppFig6_Pancreas.ipynb">Independent mode</a>
                </td>
                <td>
                    <a href="https://github.com/aicb-ZhangLabs/deepvelo/blob/main/code/Figure2.ipynb">normal AE</a>
                    <a href="https://github.com/aicb-ZhangLabs/deepvelo/blob/main/code/Figure2.ipynb">VAE</a>
                </td>
                <td></td>
                <td></td>
                <td></td>
                <td>
                    <a href="https://github.com/Spencerfar/LatentVelo/blob/main/paper_notebooks/Mouse%20erythroid.ipynb">VAE</a>
                    <a href="https://github.com/Spencerfar/LatentVelo/blob/main/paper_notebooks/Hindbrain.ipynb">annotated VAE</a>
                    <a href="https://github.com/Spencerfar/LatentVelo/blob/3cdf02946e53411aad68700ea0bba1254c48ec54/latentvelo/models/refine_model.py#L101">RefineODE</a>
                    <a href="https://github.com/Spencerfar/LatentVelo/blob/main/paper_notebooks/Mouse%20brain%20ATAC.ipynb">ATACRegMode</a>
                </td>
                <td></td>
                <td></td>
                <td></td>
            </tr>
        </tbody>
    </table>
    <script>
        document.addEventListener('DOMContentLoaded', function() {
            document.querySelectorAll('.custom-table a').forEach(function(link) {
                link.setAttribute('target', '_blank');
                link.setAttribute('rel', 'noopener noreferrer');
            });
        });
    </script>
    </div>