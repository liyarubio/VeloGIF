.. _getting_started:

Getting Started
==========================

1. Start the Container and Mount Local Directories
----------------------------------------------------------------------------------------------------------------------

- Start the container using the following command, mounting the GPU and specifying the directory:

  .. code-block:: bash

      docker run -d --name <container name> -v <your host path>:<container path> <image>

- eg. With GPU:

  .. code-block:: bash

      docker run -d --gpus all --name velogif -v /Your_Path/VeloGIF-main:/velogif velogif-registry.ap-southeast-1.cr.aliyuncs.com/velogif/velogif:latest

- eg. Without GPU:

  .. code-block:: bash

      docker run -d --name velogif -v /Your_Path/VeloGIF-main:/velogif velogif-registry.ap-southeast-1.cr.aliyuncs.com/velogif/velogif:latest


2. Access the Container Console
---------------------------------------------------------------------------------------------

- Enter the interactive terminal of the container with the following command:

  .. code-block:: bash

      docker exec -it velogif /bin/bash

- Navigate to the folder we just mounted:

  .. code-block:: bash

      cd /velogif

3. Custom Parameter
--------------------------------------------------------------------------------

- By modifying the ``config.py`` file, define the file locations for input data and output results, select the algorithm to run, as well as other parameters.

Input Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Input data should be provided in the form of |h5ad|, and the layers contain splice and unsplice matrices. Data preprocessing can refer to |velocyto| and |scvelo|.

.. |h5ad| raw:: html

   <a href="https://anndata.readthedocs.io/en/stable/" target="_blank">.h5ad</a>

.. |velocyto| raw:: html

   <a href="http://velocyto.org/velocyto.py/tutorial/cli.html#running-velocyto" target="_blank">velocyto</a>

.. |scvelo| raw:: html

   <a href="https://scvelo.readthedocs.io/en/stable/VelocityBasics.html#Preprocess-the-Data" target="_blank">scVelo</a>

- Note that |multivelo| need chromatin accessibility information, |phylovelo| need lineage information, and |tfvelo| need gene regulatory database.

.. |multivelo| raw:: html

   <a href="https://github.com/welch-lab/MultiVelo/" target="_blank">Mutivelo</a>

.. |phylovelo| raw:: html

   <a href="https://phylovelo.readthedocs.io/en/latest" target="_blank">PhyloVelo</a>

.. |tfvelo| raw:: html

   <a href="https://github.com/xiaoyeye/TFvelo" target="_blank">TFvelo</a>

- We provide |demo_data|, which is sampled 1000 cells from |redeem|. ReDeeM dataset include splice, unsplice, lineage, and chromatin accessibility information. Please put the downloaded ``Demo_data_for_VeloGIF.zip`` into the currently mounted ``VeloGIF-main`` folder, and unzip it to form the ``Input_data`` folder.

.. |demo_data| raw:: html

   <a href="https://figshare.com/ndownloader/files/51461450" target="_blank">demo data</a>

.. |redeem| raw:: html

   <a href="https://doi.org/10.1038/s41586-024-07066-z" target="_blank">ReDeeM dataset</a>

  .. code-block:: bash

      unzip Demo_data_for_VeloGIF.zip

Select Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Users can choose algorithms based on our comprehensive |benchmark|, the |characteristics|, and the |input_data| required for the different algorithms. Note that VeloGIF provides the running environment for all 15 algorithms.

.. |benchmark| raw:: html

   <a href="https://sysomics.com/velogif/benchmark/Overall_Performance.html" target="_blank">benchmark</a>

.. |characteristics| raw:: html

   <a href="https://sysomics.com/velogif/rna_velocity_methods/Methods_Introduction.html" target="_blank">characteristics</a>

.. |input_data| raw:: html

   <a href="https://sysomics.com/velogif/rna_velocity_methods/input_data_requirement.html" target="_blank">input data</a>

  .. code-block:: python

      Methods =['velocyto',
                'scvelo',
                'veloae',
                'dynamo',
                'velovae',
                'unitvelo',
                'deepvelo_vae',
                'celldancer',
                'velovi',
                'latentvelo',
                'deepVelo_gcn',
                'stt',
                'tfvelo', 
                'multivelo',
                'phylovelo']

Detailed Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Default common parameters:

  .. code-block:: python

      n_job = 10 # Number of parallel jobs.
      device = 'cuda:0' # GPU
      seed = 2024 # random seed
      embed = 'umap' # Key for embedding
      data_cluster = 'CellType' # Key for annotations of observations/cells, a column included in adata.obs
      gene_number = 2000 # Gene number
      velocity_layer = 'velocity'# Key for velocity to be visualized

- Each algorithm requires different parameters, and we keep the default parameters of the algorithm, as detailed on |default_params|. Users also can customize the parameters of each algorithm by modifying ``run_X.py``.

.. |default_params| raw:: html

   <a href="https://sysomics.com/velogif/tutorials/Default_Parameters.html" target="_blank">Default Parameters</a>

Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- VeloGIF visualizes all results by default. Users can also select the result to draw by modifying ``Methods_name`` dictionary.

Evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Users can quantitatively evaluate RNA velocity results by customizing ``edges`` list and defining cell transfer directions based on prior knowledge. eg. The transformation from Hematopoietic stem cells (HSC) to multipotent progenitor (MPP).

4. Execute the Scripts
------------------------------------------------------------------------------------

- After entering the container, navigate to the application directory and run the calculation script:

  .. code-block:: bash

      python Run_all_methods.py

5. Exploring the Output
------------------------------------------------------------------------------------

.. code-block:: text

    result
    ├── evals
    │   └── Eval.csv                    # GDC, CBDir, and ICCoh value of each method
    ├── figures
    │   ├── Merge.svg                   # Velocity stream of all methods
    │   ├── cellDancer.svg              # Velocity stream of each method
    │   ├── DeepVelo (GCN-based).svg
    │   ├── DeepVelo (VAE-based).svg
    │   ├── Dynamo.svg
    │   ├── LatentVelo.svg
    │   ├── MultiVelo.svg
    │   ├── scVelo (dynamic).svg
    │   ├── scVelo (stochastic).svg
    │   ├── STT.svg
    │   ├── TFvelo.svg
    │   ├── UniTVelo.svg
    │   ├── veloAE.svg
    │   ├── velocyto.svg
    │   ├── veloVAE.svg
    │   └── veloVI.svg
    ├── execution_log.txt               # Log file for running all methods
    ├── celldancer.h5ad                 # .h5ad files contain results of each methods
    ├── deepvelo_gcn.h5ad
    ├── deepvelo_vae.h5ad
    ├── dynamo.h5ad
    ├── latentvelo.h5ad
    ├── multivelo.h5ad
    ├── phylovelo.h5ad
    ├── scvelo.dyn.h5ad
    ├── scvelo.sto.h5ad
    ├── stt.h5ad
    ├── tfvelo.svg
    ├── unitvelo.h5ad
    ├── veloae.h5ad
    ├── velocyto.h5ad
    ├── velovae.h5ad
    └── velovi.h5ad


