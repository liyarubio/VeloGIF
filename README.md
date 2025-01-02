# VeloGIF
<img src="https://github.com/liyarubio/VeloGIF/blob/main/Figure/VeloGIF.png" width="50%" height="50%">
More details in our website : https://sysomics.com/velogif

## Project Overview
Our [benchmark of RNA velocity methods](https://sysomics.com/velogif/benchmark/Framework_and_Datasets.html) indicated a significant complementarity among the various methods, leading us to propose an integration tool, ```VeloGIF```. The ```Graphics Interchange Format (GIF)``` is a prevalent image file format widely utilized to create live photos due to its excellent compatibility. Our integrated tool employed multiple RNA velocity methods to expand static transcriptional snapshots to dynamic photos, and it is also highly compatible, hence we named it VeloGIF. ```VeloGIF``` allows users to concurrently obtain the computation results of 15 various RNA velocity methods, as well as visualization and evaluation results, and could select the outcomes that align most closely with the expectations based on prior knowledge. VeloGIF employed [Docker container](https://www.docker.com) to wrap each method, which avoids dependency issues.

## Installation

### 1. Install Docker (Required)

- Visit the [Docker official website](https://www.docker.com) and download Docker.
  
- Add docker permissions to the current user.
   ```
   sudo groupadd docker 
   sudo gpasswd -a ${USER} docker
   sudo service docker restart
   ```
   
- Verify that Docker is running correctly:
   ```
   docker --version
   ```
   
### 2. Install NVIDIA Container Toolkit (Optional, For CUDA Acceleration)

If you need to use GPU and CUDA, install the NVIDIA Container Toolkit.

- Ensure that your system has NVIDIA drivers installed and that CUDA is properly configured.

   - You can check if the NVIDIA drivers are available using the following command:

   ```
   nvidia-smi
   ```
   
- Install the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html):

### 3. Set the Image Repository and Pull the Image

- configure the provided image repository address and account information, then pull the required Docker image:
```
docker pull your_repository/velogif:latest
```

## Getting started
### 1. Start the Container (Mount GPU and Specify Directory)

-Start the container using the following command, mounting the GPU and specifying the directory:

```
docker run -d --name <container name> -v <your host path>:<container path> <image>
```
- eg. With GPU:
```
docker run -d --gpus all --name velogif -v /Your_Path:/velogif your_repository/velogif:latest
```

- eg. Without GPU:
```
docker run -d --name velogif -v /Your_Path:/velogif your_repository/velogif:latest
```

### 2. Access the Container Console

- Enter the interactive terminal of the container with the following command:

```
docker exec -it velogif /bin/bash
```

- Navigate to the folder we just mounted
```
cd /velogif
```

### 3. Custom parameter
- By modifying the ```config.py``` file, define the file locations for input data and output results, select the algorithm to run, as well as other parameters.

#### Input data
- Input data should be provided in the form of [.h5ad](https://anndata.readthedocs.io/en/stable/), and the layers contain splice and unsplice matrices. Data preprocessing can refer to [velocyto](http://velocyto.org/velocyto.py/tutorial/cli.html#running-velocyto) and [scVelo](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html#Preprocess-the-Data).
  
- Note that [Mutivelo](https://github.com/welch-lab/MultiVelo/) need chromatin accessibility information, [PhyloVelo](https://phylovelo.readthedocs.io/en/latest) need lineage information, and [TFvelo](https://github.com/xiaoyeye/TFvelo) need gene regulatory database. 

- We provide demo data, which is sampled 1000 cells from [ReDeeM dataset](https://doi.org/10.1038/s41586-024-07066-z). ReDeeM dataset include splice, unsplice, lineage, and chromatin accessibility information.

#### Select algorithms
- Users can choose algorithms based on our comprehensive [benchmark](https://sysomics.com/velogif/benchmark/Overall_Performance.html), the [characteristics](https://sysomics.com/velogif/rna_velocity_methods/Methods_Introduction.html), and the [input data](https://sysomics.com/velogif/rna_velocity_methods/input_data_requirement.html) required for the different algorithms. Note that VeloGIF provides the running environment for all 15 algorithms.
  
```
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
```

#### Detailed parameters
- Default common parameters
```
n_job = 10 # Number of parallel jobs.
device = 'cuda:0' # GPU
seed = 2024 # random seed
embed = 'umap' # Key for embedding
data_cluster = 'CellType' # Key for annotations of observations/cells, a column included in adata.obs
gene_number = 2000 # Gene number 
```

- Each algorithm requires different parameters, and we keep the default parameters of the algorithm, as detailed on [Default Parameters](https://sysomics.com/velogif/tutorials/Default_Parameters.html). Users also can customize the parameters of each algorithm by modifying ```run_X.py```.

#### Visulization
- VeloGIF visualizes all results by default. Users can also select the result to draw by modifying ```Methods_name```dictionary.

#### Evaluation
- Users can quantitatively evaluate RNA velocity results by customizing ```edges``` list and defining cell transfer directions based on prior knowledge, for example, The transformation from Hematopoietic stem cells (HSC) to multipotent progenitor (MPP).

### 4. Execute the Scripts
- After entering the container, navigate to the application directory and run the calculation script:
```
python Run_all_methods.py
```
### 5. Exploring the Output 
```
result
├── evals
│   └── Eval.csv                    # GDC, CBDir, and ICCoh value of each method
├── figures
│   ├── Merge.svg                   # Velocity stream of all methods
│   ├── cellDancer.svg              # Velocity stream of each method
│   ├── DeepVelo (GCN-based).svg
│   ├── DeepVelo (VAE-based).svg
│   ├── Dynamo.svg
│   ├── LatentVelo.svg
│   ├── MultiVelo.svg
│   ├── scVelo (dynamic).svg
│   ├── scVelo (stochastic).svg
│   ├── STT.svg
│   ├── TFvelo.svg
│   ├── UniTVelo.svg
│   ├── veloAE.svg
│   ├── velocyto.svg
│   ├── veloVAE.svg
│   └── veloVI.svg
├── execution_log.txt               # Log file for running all methods
├── celldancer.h5ad                 # .h5ad files contain results of each methods, the RNA velocity result in adata.layers['velocity']
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
```

## Contact Us
If you have any questions or suggestions, please contact the project maintainers:
- Email: bio_liyaru@163.com
- GitHub Repository: https://github.com/liyarubio/VeloGIF

### Thank you for your support and usage!
