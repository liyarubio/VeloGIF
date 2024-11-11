# VeloGIF

## Project Overview

This project provides a pre-configured Docker image containing a fully functional Conda environment, which supports 12 different execution sub-environments. Additionally, it includes convenient scripts for one-click computation, automated figure generation, and packaging, allowing users to quickly set up and use the environment.

![image](https://github.com/liyarubio/VeloGIF/blob/main/Figure/Fig1A.png)


## Prerequisites

Please ensure that you have installed the latest version of Docker and have configured it properly.

## Steps to Use

### 1. Configure the Image Repository and Pull the Image

First, set up the provided image repository address and account information, then pull the required Docker image.

### 2. Start the Container (with GPU and Directory Mounting)

Use the following command to start the container with GPU support and mount the specified directory:

```
docker run -d --gpus all --name rnagifv2 -v /opt/pysc:/app rnagifv2:latest
```

### 3. Access the Container Console

To access the interactive terminal within the container, use the following command:

```
docker exec -it rnagifv2 /bin/bash
```

### 4. Execute the Run Script within the Container

Once inside the container, navigate to the application directory and execute the computation script:

```
cd /app
./run_all.sh
```

After the computation completes, you may proceed with the figure generation and packaging scripts. Ensure that the environment and input/output directories for figures are set correctly:

```
./get_all.sh
./package.sh
```

Upon completion, a ZIP file containing the packaged output will be available in the specified directory.
