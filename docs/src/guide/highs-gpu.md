# [HiGHS on a GPU](@id highs-gpu)

From HiGHS v1.10.0, its first order primal-dual LP (PDLP) solver [cuPDLP-C](https://github.com/COPT-Public/cuPDLP-C) can be run on an NVIDIA GPU under Linux and Windows. However, to achieve this, CUDA utilities must be installed and HiGHS must be built locally using `cmake`, as described below.

### PDLP: A health warning

First order solvers for LP are still very much ``work in progress''. Although impressive results have been reported, these are often to lower accuracy than is achieved by simplex and interior point solvers, have been obtained using top-of-the-range GPUs, and not achieved for all problem classes.

#### Termination criteria

Although the PDLP solver may report that it has terminated with an optimal solution, HiGHS may identify that the solution returned by PDLP is not optimal. To explain how this can occur, and allow users to decide whether to accept such a non-optimal solution, it is necessary to discuss the termination criteria used by PDLP.

For the LP
```math
\begin{aligned}
\min                \quad & c^T\! x        \\
\textrm{subject to} \quad & Ax = b  \\
                          & x \ge 0,
\end{aligned}
```
the optimality conditions are that, at a point ``x``, there exist ``y`` and ``s`` such that
```math
\begin{aligned}
Ax=b&\textrm{Primal~equations}\\
A^Ty+s=c&\textrm{Dual~equations}\\
x\ge0&\textrm{Primal~feasibility}\\
s\ge0&\textrm{Dual~feasibility}\\
c^Tx-b^Ty=0&\textrm{Primal-dual~gap}
\end{aligned}
```
The PDLP algorithm determines values of ``x\ge0`` and ``y``, and chooses ``s`` to be the non-negative values of ``c-A^Ty``. Hence it guarantees primal and dual feasibility by construction. PDLP terminates when
```math
\begin{aligned}
\|Ax-b\|_2&\le \epsilon_P(1+\|b\|_2)\\
\|c-A^Ty-s\|_2&\le \epsilon_D(1+\|c\|_2)\\
|c^Tx-b^Ty|&\le \epsilon_{PD}(|c^Tx|+|b^Ty|)
\end{aligned}
```
where the value of ``\epsilon_P`` used is the HiGHS option [`primal_feasibility_tolerance`](primal_feasibility_tolerance), the value of ``\epsilon_D`` used is the HiGHS option [`dual_feasibility_tolerance`](dual_feasibility_tolerance), and the value of ``\epsilon_{PD}`` used is the HiGHS option [`pdlp_d_gap_tol`](p_gap_tol). 

If you use the HiGHS PDLP solver, in the first instance it is recommended that you increase the primal and dual feasibility tolerances to `1e-4`, since this will result in the algorithm terminating much sooner.

### Requirements

CUDA Toolkit and CMake. 

A [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit) installation is required, along with the matching NVIDIA driver. Please install both following the instructions on NVIDIA's website.

HiGHS must be build locally with CMake. 

Make sure the CUDA compiler `nvcc` is installed by running 

```
nvcc --version
```

### Build HiGHS with GPU support

HiGHS must be built, from the root directory, with 

```
cmake -S. -Bbuild -DCUPDLP_GPU=ON
cmake --build build --parallel
```

This uses [FindCUDAToolkit](https://cmake.org/cmake/help/latest/module/FindCUDAToolkit.html) to find a CUDA installation locally.

#### Find CUDA

If CUDA is not found automatically, there is an extra option `-DCUPDLP_FIND_CUDA=ON`, to be used with `-DCUPDLP_GPU=ON`, which instead uses `cuPDLP-C`'s `FindCUDAConf.cmake`. 

This requires the environment variable `CUDA_HOME` to be set to the directory with the CUDA installation. Having set this, run 

```
cmake -S. -Bbuild -DCUPDLP_GPU=ON -DCUPDLP_FIND_CUDA=ON
cmake --build build --parallel
```

to build HiGHS. 
