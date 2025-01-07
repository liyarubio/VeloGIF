.. _main_installation:



Installation
===============================


1. Install Docker (Required)
-------------------------------------------------------------------------------------------

• Visit the |docker_link| and download Docker.

.. |docker_link| raw:: html

   <a href="https://www.docker.com/" target="_blank">Docker official website</a>

• Add Docker permissions to the current user:

   .. code-block:: bash

       sudo groupadd docker
       sudo gpasswd -a ${USER} docker
       sudo service docker restart

• Verify that Docker is running correctly:

   .. code-block:: bash

       docker --version

2. Install NVIDIA Container Toolkit (Optional, For CUDA Acceleration)
-------------------------------------------------------------------------------------------

If you need to use GPU and CUDA, install the |nvidia_link|. Ensure that your system has NVIDIA drivers installed and that CUDA is properly configured.

.. |nvidia_link| raw:: html

   <a href="https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html" target="_blank">NVIDIA Container Toolkit</a>

• You can check if the NVIDIA drivers are available using the following command:

   .. code-block:: bash

       nvidia-smi

Install the |nvidia_toolkit_link|.

Note that CUDA is required for LatentVelo.

.. |nvidia_toolkit_link| raw:: html

   <a href="https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html" target="_blank">NVIDIA Container Toolkit</a>

3. Set the Image Repository and Pull the Image
-------------------------------------------------------------------------------------------

• Configure the provided image repository address and account information, then pull the required Docker image:

   .. code-block:: bash

       docker pull velogif-registry.ap-southeast-1.cr.aliyuncs.com/velogif/velogif:latest

• Check local images and ensure that the VeloGIF image is successfully pulled:

   .. code-block:: bash

       docker images

4. Download the Scripts and Extract Files
-------------------------------------------------------------------------------------------

• Download the scripts from the |velogif_link|.

.. |velogif_link| raw:: html

   <a href="https://github.com/liyarubio/VeloGIF/archive/refs/heads/main.zip" target="_blank">VeloGIF GitHub repository</a>

• Extract the downloaded files:

   .. code-block:: bash

       unzip VeloGIF-main.zip
       cd VeloGIF-main/velogif/
       unzip tools.zip