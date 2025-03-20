# [HiGHS on a GPU](@id highs-gpu)

## Using HiGHS on a GPU

HiGHS can run `cuPDLP-C` on Nvidia GPUs.

### Requiremets

CUDA Toolkit and CMake. 

A [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit) installation is required, along with the matching NVidia driver. Please install both following the instructions on Nvidia's website.

HiGHS must be build locally with CMake. 

Make sure the cuda compiler `nvcc` is installed by running 

```
nvcc --version
```

### Build HiGHS with GPU support

HiGHS must be built, from the root directory, with 

```
cmake -S. -Bbuild -DCUPDLP_GPU=ON
cmake --build build --parallel
```

This uses [FindCUDAToolkit](https://cmake.org/cmake/help/latest/module/FindCUDAToolkit.html) to find a cuda installation locally.

#### Find CUDA

If cuda is not found automatically, there is an extra option `-DCUPDLP_FIND_CUDA=ON`, to be used with `-DCUPDLP_GPU=ON`, which instead uses `cuPDLP-C`'s `FindCUDAConf.cmake`. 

This requires the environment variable `CUDA_HOME` to be set to the directory with the cuda installation. Having set this, run 

```
cmake -S. -Bbuild -DCUPDLP_GPU=ON -DCUPDLP_FIND_CUDA=ON
cmake --build build --parallel
```

to build HiGHS. 