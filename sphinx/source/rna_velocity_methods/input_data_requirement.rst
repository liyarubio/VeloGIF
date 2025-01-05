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

    /* 首列固定样式 - 调整宽度 */
    .wy-table-responsive table th:first-child,
    .wy-table-responsive table td:first-child {
        position: sticky;
        left: 0;
        background-color: white;
        z-index: 2;
        border-right: 2px solid #e1e4e5;
        min-width: 200px;  /* 为首列设置最小宽度 */
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

    /* 其他列的宽度 */
    .wy-table-responsive table th,
    .wy-table-responsive table td {
        min-width: 80px;  /* 为其他列设置合适的宽度 */
        text-align: center;  /* 居中对齐 */
    }
    </style>

Input Data Requirement
========================

.. list-table::
   :header-rows: 1

   * - Data Type
     - velocyto
     - scVelo
     - VeloAE
     - Dynamo
     - VeloVAE
     - MultiVelo
     - UniTVelo
     - DeepVelo(VAE)
     - cellDancer
     - PhyloVelo
     - veloVI
     - LatentVelo
     - DeepVelo(GCN)
     - TFvelo
     - STT
   * - conventional scRNA-seq
     - √
     - √
     - √
     - √
     - √
     - 
     - √
     - √
     - √
     - 
     - √
     - √
     - √
     - √
     - √
   * - Without unsplice information
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - √
     - 
     - 
     - 
     - √
     - 
   * - scRNA-seq with metabolic label
     - 
     - 
     - 
     - √
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
   * - spatial transcriptome
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - √
   * - scRNA-seq with scATAC-seq
     - 
     - 
     - 
     - 
     - 
     - √
     - 
     - 
     - 
     - 
     - 
     - √
     - 
     - 
     - 
   * - scRNA-seq with lineage data
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - √
     - 
     - 
     - 
     - 
     - 


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