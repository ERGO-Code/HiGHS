## Building HiGHS to run PDLP on a GPU

### Installation

#### CUDA-Toolkit 

https://developer.nvidia.com/cuda-toolkit

#### nvidia-smi check drivers

To utilize the maximum potential of any hardware, relevant drivers for that hardware are installed onto the system. There exists a tool named “nvidia-smi” which helps to manage NVidia GPU hardware all through the terminal on the Linux operating system. 

* Before installing anything, run the following in a terminal to check if nvidia-smi is present and can see the GPUs

```
nvidia-smi
```

This should print a table of GPU-related stats. If it doesn't don't worry: installing the Cuda Toolkit from https://developer.nvidia.com/cuda-toolkit should help!

### Install CUDA and Driver from the Nvidia website

Download Cuda Tooklit and the Nvidia driver from https://developer.nvidia.com/cuda-toolkit for the appropriate OS and architecture.

Note, that installing the `nvidia-cuda-toolkit` package from apt will work on many Linux installations, *but* if the drivers are missing and you attempt to install them with
```
apt install nvidia-cuda-toolkit
apt install nvidia-utils-515
```
or `nvidia-utils-some-other-number`, this may lead to errors with versions of cuda libraries on versions of Lunux that are resolved by cleaning up the cuda drivers and toolkit from the machine and installing again from the website. 

Check again that the drivers are present by running 
```
nvidia-smi
```

* Run the following in a terminal

```
nvcc --version
```

If you see a version, this means that the cuda compiler nvcc is up and running.

When both work, on Linux, you should check that `/usr/local/cuda` exists. 

#### Build

```
cmake -S . -B build -DCUPDLP_GPU=ON -DALL_TESTS=ON ..
cmake --build build --parallel
```

