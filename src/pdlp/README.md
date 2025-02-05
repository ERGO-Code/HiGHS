## Building HiGHS to run PDLP on a GPU

### Installation

#### CUDA-Toolkit 

https://developer.nvidia.com/cuda-toolkit

#### Check and further installation

* Before installing anything more, run the following in a terminal

``
nvidia-smi
``

This should print a table of GPU-related stats. If it doesn't, you should install the “nvidia-utils” package which will also contain the “nvidia-smi” tool inside it. To install this package, run the command in the terminal:

``
sudo apt install nvidia-utils-515
``

* Run the following in a terminal

``
nvcc --version
``

This should give you details about the cuda compiler driver. If it doesn't then run

``
apt install nvidia-cuda-toolkit
``

When both work, you should check that `/usr/local/cuda` exists. 

#### Build

* Checkout the HiGHS codebase from Github, and use the branch `cuda-updates`
* From your HiGHS build directory, run

``
cmake -DCUPDLP_GPU=ON -DALL_TESTS=ON ..; make
``






