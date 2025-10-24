# [GPU acceleration](@id gpu)

From HiGHS v1.10.0, its first order primal-dual LP (PDLP) solver
[cuPDLP-C](https://github.com/COPT-Public/cuPDLP-C) can be run on an
NVIDIA GPU under Linux and Windows. However, to achieve this, CUDA
utilities must be installed and HiGHS must be built locally using
CMake, as described below.

### PDLP: A health warning

First order solvers for LP are still very much "work in
progress". Although impressive results have been reported, these are
often to lower accuracy than is achieved by simplex and interior point
solvers, have been obtained using top-of-the-range GPUs, and not
achieved for all problem classes. Note that, due to PDLP using
relative termination conditions, a solution deemed optimal by PDLP may
not be accepted as optimal by HiGHS. The user should consider the
infeasibility data returned by [HighsInfo](@ref HighsInfo) to decide
whether the solution is acceptable to them.

#### Termination criteria

Although the PDLP solver may report that it has terminated with an
optimal solution, HiGHS may identify that the solution returned by
PDLP is not optimal. As discussed in [HiGHS feasibility and optimality
tolerances](@ref kkt), this is due to PDLP using relative termination
criteria and (unlike interior point solvers) not satisfying
feasibility to high accuracy.

If you use the HiGHS PDLP solver, in the first instance it is
recommended that you increase the feasibility and optimality
tolerances to `1e-4`, since this will result in the algorithm
terminating much sooner. There are multiple feasibility and optimality
tolerances, but all will be set to the value of the
[`kkt_tolerance`](@ref option-kkt-tolerance) option (if it differs
from its default value of `1e-4`) so this is recommended in the first
instance.

### Requirements

CUDA Toolkit and CMake. 

A [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit)
installation is required, along with the matching NVIDIA
driver. Please install both following the instructions on NVIDIA's
website.

HiGHS must be build locally with CMake. 

Make sure the CUDA compiler `nvcc` is installed by running 

```
nvcc --version
```

### Build HiGHS with GPU support

See [Building HiGHS with NVidia GPU support](@ref gpu-build).