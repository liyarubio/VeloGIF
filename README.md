# VeloGIF
![image](https://github.com/liyarubio/VeloGIF/blob/main/Figure/Fig1A.png)

## Project Overview
Our benchmark of RNA velocity indicated a significant complementarity among the various methods, leading us to propose an integration tool, VeloGIF. VeloGIF allows users to concurrently obtain the computation results of 15 various RNA velocity methods, as well as visualization and evaluation results, and could select the outcomes that align most closely with the expectations based on prior knowledge. VeloGIF employed Docker container to wrap each method, which avoids dependency issues.

## Installation

### 1. Install Docker (Required)

1. Visit the [Docker official website](https://www.docker.com) and download Docker.
2. Add docker permissions to the current user.
   ```
   sudo groupadd docker 
   sudo gpasswd -a ${USER} docker
   sudo service docker restart
   ```
2. Verify that Docker is running correctly:
   ```
   docker --version
   ```
   
### 2. Install NVIDIA Container Toolkit (Optional, For CUDA Acceleration)

If you need to use GPU and CUDA, install the NVIDIA Container Toolkit.

1. Ensure that your system has NVIDIA drivers installed and that CUDA is properly configured.

   - You can check if the NVIDIA drivers are available using the following command:

   ```
   nvidia-smi
   ```
   
2. Install the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html):


## Getting started

### 1. Set the Image Repository and Pull the Image

First, configure the provided image repository address and account information, then pull the required Docker image:

```
docker pull your_repository/velogif:latest
```

### 2. Start the Container (Mount GPU and Specify Directory)

Start the container using the following command, mounting the GPU and specifying the directory:

- Mount the host directory to the container using the -v option.

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

### 3. Access the Container Console

Enter the interactive terminal of the container with the following command:

```
docker exec -it velogif /bin/bash
```

### 4. Execute the Scripts for Demo Data in the Container

After entering the container, navigate to the application directory and run the calculation script:
- The demo data is sampled 1000 cells from [Redeem dataset](https://doi.org/10.1038/s41586-024-07066-z)

```
cd /velogif
python Run_all_methods.py
```

## Input data



## Notes

1. Ensure your host environment has sufficient GPU resources and that CUDA dependencies are properly installed.
2. If you encounter any issues during usage, check whether Docker and the NVIDIA Container Toolkit are installed correctly, or review the relevant logs to troubleshoot.
3. If you need to customize scripts or extend functionality, modify or add files in the `/app` directory and rebuild the Docker image as needed.

## Contact Us

If you have any questions or suggestions, please contact the project maintainers:

- Email: [support@velogif.com](mailto:support@velogif.com)
- GitHub Repository: https://github.com/your_repository/velogif

Thank you for your support and usage!
