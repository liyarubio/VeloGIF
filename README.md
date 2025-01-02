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

- With GPU:
- Mount the host directory to the container using the -v option.
```
docker run -d --name <container name> -v <your host path>:<container path> <image>
```
- eg.
```
docker run -d --gpus all --name velogif -v /Your_Path:/velogif your_repository/velogif:latest
```

- Without GPU:
```
docker run -d --name velogif -v /Your_Path:/velogif your_repository/velogif:latest
```

### 3. Access the Container Console

Enter the interactive terminal of the container with the following command:

```
docker exec -it velogif /bin/bash
```

### 4. Execute the Scripts in the Container

After entering the container, navigate to the application directory and run the calculation script:

```
cd /app
python run_all.py
```

Once the computation is complete, you can generate charts by running the following script:

```
python Get_graph.py
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









---

# velogif 项目

## 项目概述

本项目提供了一个预配置的 Docker 镜像，其中包含了一个完整的 Conda 环境，支持 15 种不同方法的执行子环境。此外，项目还提供了一键运行的计算脚本和一键生成图表的脚本，方便用户快速启动和使用。

## 准备工作

### 1. 安装 Docker（必须）

1. 访问 [Docker 官方网站](https://www.docker.com/products/docker-desktop/) 并下载最新版 Docker。
2. 根据操作系统（Windows、macOS 或 Linux）的具体指南完成安装：
   - **Windows**: 下载 Docker Desktop 安装程序，运行后根据提示完成安装，安装完成后需要启用 WSL 2 功能。
   - **macOS**: 下载对应 Apple Silicon 或 Intel 芯片的 Docker Desktop 版本，并运行安装程序完成安装。
   - **Linux**: 执行以下命令安装 Docker（以 Ubuntu 为例）：
     ```
     sudo apt-get update
     sudo apt-get install -y ca-certificates curl gnupg
     sudo install -m 0755 -d /etc/apt/keyrings
     curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
     echo \
       "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
       $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
     sudo apt-get update
     sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
     ```
3. 启动 Docker 服务（仅限 Linux）：
   ```
   sudo systemctl start docker
   sudo systemctl enable docker
   ```
4. 验证 Docker 是否正常运行：
   ```
   docker --version
   ```

### 2. 安装 NVIDIA Container Toolkit（以支持 CUDA，加速计算过程）

如果您需要使用 GPU 和 CUDA，请安装 NVIDIA Container Toolkit。

1. 确保系统已安装 NVIDIA 驱动程序，且 CUDA 已正确配置。
   - 您可以通过以下命令检查 NVIDIA 驱动程序是否可用：

   ```
   nvidia-smi
   ```

2. 安装 NVIDIA Container Toolkit：
   - 添加 NVIDIA 的软件包仓库：

     ```
     distribution=$(. /etc/os-release;echo $ID$VERSION_ID) \
     && curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
     && curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | \
        sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
     ```

   - 更新包列表并安装 NVIDIA Container Toolkit：

     ```
   sudo apt-get update
   sudo apt-get install -y nvidia-container-toolkit
     ```

   - 配置容器运行时并重启 Docker：

     ```
     sudo nvidia-ctk runtime configure --runtime=docker
     sudo systemctl restart docker
     ```

3. 测试 GPU 是否在容器中正常运行：

   ```
   docker run --rm --gpus all nvidia/cuda:12.2.0-base nvidia-smi
   ```

如果上述命令输出了 GPU 信息，则说明安装成功。

## 使用步骤

### 1. 设置镜像仓库并拉取镜像

首先，配置提供的镜像仓库地址和账户信息，然后拉取所需的 Docker 镜像：

```
docker pull your_repository/velogif:latest
```

### 2. 启动容器（挂载 GPU 与指定目录）

使用以下命令启动容器，并挂载 GPU 与指定目录：

```
docker run -d --gpus all --name velogif -v /opt/pysc:/app velogif:latest
```

### 3. 进入容器控制台

使用以下命令进入容器的交互式终端：

```
docker exec -it velogif /bin/bash
```

### 4. 在容器中执行运行脚本

进入容器后，切换到应用目录并运行计算脚本：

```
cd /app
python run_all.py
```

计算完成后，您可以继续执行图表生成脚本：

```
python Get_graph.py
```

## 注意事项

1. 请确保您的主机环境具备足够的 GPU 资源，以及 CUDA 的相关依赖已正确安装。
2. 如果在使用过程中遇到任何问题，请检查是否正确安装了 Docker 和 NVIDIA Container Toolkit，或查阅相关日志以定位问题。
3. 如果需要自定义脚本或扩展功能，请在 `/app` 目录下修改或新增文件，并根据需要重新构建 Docker 镜像。

## 联系我们

如有任何疑问或建议，请联系项目维护者：

- 邮箱：support@velogif.com
- GitHub 仓库：[[Issues · liyarubio/VeloGIF](https://github.com/liyarubio/VeloGIF/issues)](https://github.com/your_repository/velogif)

感谢您的使用与支持！
