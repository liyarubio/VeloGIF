# VeloGIF
![image](https://github.com/liyarubio/VeloGIF/blob/main/Figure/Fig1A.png)

## Project Overview
Our benchmark of RNA velocity indicated a significant complementarity among the various methods, leading us to propose an integration tool, VeloGIF. VeloGIF allows users to concurrently obtain the computation results of 15 various RNA velocity methods, as well as visualization and evaluation results, and could select the outcomes that align most closely with the expectations based on prior knowledge. VeloGIF employed Docker container to wrap each method, which avoids dependency issues.

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

-Enter the interactive terminal of the container with the following command:

```
docker exec -it velogif /bin/bash
```
### 3. Custom parameter
- By modifying the ```config.py``` file, define the file locations for input data and output results, select the algorithm to run, as well as other parameters.

#### Input data
- Input data should be provided in the form of [.h5ad](https://anndata.readthedocs.io/en/stable/), and the layers contain splice and unsplice matrices. Data preprocessing can refer to [velocyto](http://velocyto.org/velocyto.py/tutorial/cli.html#running-velocyto) and [scVelo](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html#Preprocess-the-Data).
  
- Note that [Mutivelo](https://github.com/welch-lab/MultiVelo/) need chromatin accessibility information, [PhyloVelo](https://phylovelo.readthedocs.io/en/latest) need lineage information, and [TFvelo](https://github.com/xiaoyeye/TFvelo) need gene regulatory database. 

- We provide demo data, which is sampled 1000 cells from [ReDeeM dataset](https://doi.org/10.1038/s41586-024-07066-z). ReDeeM dataset include splice, unsplice, lineage, and chromatin accessibility information.

#### Select algorithms
- Users can choose algorithms based on our comprehensive [benchmark](https://sysomics.com/velogif/benchmark/Overall_Performance.html), the [characteristics], and the [input data](https://sysomics.com/velogif/rna_velocity_methods/Methods_Introduction.html) required for the different algorithms. Note that VeloGIF provides the running environment for all 15 algorithms.

```
Methods =[
    # conventional scRNA-seq
    'velocyto',
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
    'tfvelo', # with GRN
    'multivelo', # with scATAC
    'phylovelo'  # with lineage info.
]
```

### 4. Execute the Scripts for Demo Data in the Container

- After entering the container, navigate to the application directory and run the calculation script:

```
cd /velogif
python Run_all_methods.py
```




## Notes

1. Ensure your host environment has sufficient GPU resources and that CUDA dependencies are properly installed.
2. If you encounter any issues during usage, check whether Docker and the NVIDIA Container Toolkit are installed correctly, or review the relevant logs to troubleshoot.
3. If you need to customize scripts or extend functionality, modify or add files in the `/app` directory and rebuild the Docker image as needed.

## Contact Us

If you have any questions or suggestions, please contact the project maintainers:

- Email: [support@velogif.com](mailto:support@velogif.com)
- GitHub Repository: https://github.com/your_repository/velogif

Thank you for your support and usage!
