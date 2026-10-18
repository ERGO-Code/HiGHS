# Install HiGHS

## Compile from source

HiGHS uses CMake as build system, and requires at least version
3.15. Details about building from source using CMake can be found in `HiGHS/cmake/README.md`.

### Bazel build

Alternatively, building with Bazel is supported for Bazel-based projects. To build HiGHS, from the root directory, run

```
bazel build //...
```

## Install via a package manager

HiGHS can be installed using a package manager in the cases of
[`Julia`](@ref HiGHS.jl), [`Python`](@ref python-getting-started), [`CSharp`](@ref nuget) and [`Rust`](@ref Rust).

Note, that HiGHS is available via apt on Linux. For simplex, ipx and the MIP and QP solvers, the execution should be as expected. We advise users not to use HiPO from the apt installation, the Metis version linked there is not thread safe. If you consider using HiPO, please use the binaries linked below, compilation from source or the python wrapper.

## Precompiled Binaries

Precompiled static binaries are available at https://github.com/ERGO-Code/HiGHS/releases.

Additionally, there is one package containing shared libraries for Windows x64.

The `*-mit` binary packages contain HiGHS and are MIT-licenced.
The `*-apache` binary packages contain HiGHS with HiPO and are Apache-licenced, due to the licensing of the dependencies of HiPO. For more information, see [THIRD_PARTY_NOTICES.md](https://github.com/ERGO-Code/HiGHS/blob/master/THIRD_PARTY_NOTICES.md).

If you have any questions or requests for more platforms and binaries, please get in touch with us at hello@highs.dev.

To install a precompiled binary, download and extract the archive corresponding to your Operating System and architecture, the executable is located at `/bin/highs`.

## HiGHS with HiPO

HiGHS does not have any external dependencies, however, the new interior point solver HiPO uses BLAS. At the moment HiPO is optional and can be enabled via CMake.

### External ordering heuristics

HiPO also relies on a fill-reducing ordering heuristic. HiGHS includes the source code of Metis, AMD and RCM, three open-source ordering heuristics. Their source code is already part of the HiGHS library, so there is no need to link them. In particular, there is no need to have Metis installed separately, as in previous versions of HiPO. These source codes can be found in extern/metis, extern/amd, extern/rcm, together with the respective license files. Notice that the HiGHS source code is MIT licensed. However, if you build HiGHS with HiPO support, then libhighs and the HiGHS executables are licensed Apache 2.0, due to the presence of Metis and AMD.

### BLAS

On MacOS no BLAS installation is required because HiPO uses [Apple Accelerate](https://developer.apple.com/accelerate/) by default.

On Windows and Linux, you can either compile OpenBLAS at configure time using the option `-DBUILD_OPENBLAS=ON` (`OFF` by default) or compile BLAS using the instructions below.

#### MacOS

To build HiPO on MacOS, run
```
cmake -S. -B build -DHIPO=ON
```

#### Linux and Windows: Compile OpenBLAS at configure time

```
cmake -S. -B build -DHIPO=ON -DBUILD_OPENBLAS=ON
```

#### Linux and Windows: Link with BLAS installatied on your machine

On Linux, libblas and libopenblas are supported. We recommend libopenblas for its better performance, and it is found by default if available on the system. Install with

```
sudo apt update
sudo apt install libopenblas-dev
```

To build HiPO, run
```
cmake -S. -B build -DHIPO=ON
```

On Windows, OpenBLAS is required. It could be installed via [vcpkg](https://learn.microsoft.com/en-us/vcpkg/get_started/overview) with

```
vcpkg install openblas[threads]:x64-windows-static
```

Note, that `[threads]` is required for HiPO.

On Windows, you also need to specify the path to OpenBLAS. If it was installed with vcpkg as suggested above, add the path to `vcpkg.cmake` to the CMake flags, e.g.
```
cmake -S. -B build -DHIPO=ON -DCMAKE_TOOLCHAIN_FILE="C:/vcpkg/scripts/buildsystems/vcpkg.cmake"
```

##### Path to BLAS

To specify explicitly which BLAS vendor to look for, `BLA_VENDOR` could be set in CMake, e.g. `-DBLA_VENDOR=Apple` or `-DBLA_VENDOR=OpenBLAS`. Alternatively, to specify which BLAS library to use, set `BLAS_LIBRARIES` to the full path of the library e.g. `-DBLAS_LIBRARIES=/path_to/libopenblas.so`.


## [Building HiGHS with NVidia GPU support](@id gpu-build)

HiGHS must be built, from the root directory, with

```
cmake -S. -Bbuild -DCUPDLP_GPU=ON
cmake --build build --parallel
```

This uses [FindCUDAToolkit](https://cmake.org/cmake/help/latest/module/FindCUDAToolkit.html) to find a CUDA installation locally. For more details on HiGHS with CMake, see `HiGHS/cmake/README.md`.


#### Find CUDA

If CUDA is not found automatically, there is an extra option `-DCUPDLP_FIND_CUDA=ON`, to be used with `-DCUPDLP_GPU=ON`, which instead uses `cuPDLP-C`'s `FindCUDAConf.cmake`.

This requires the environment variable `CUDA_HOME` to be set to the directory with the CUDA installation. Having set this, run

```
cmake -S. -Bbuild -DCUPDLP_GPU=ON -DCUPDLP_FIND_CUDA=ON
cmake --build build --parallel
```

to build HiGHS.

### Bazel build with CUDA

Alternatively, for Bazel run

```
bazel build //... --//:cupdlp_gpu
```

It may be necessary to also specify the architecture, e.g.

```
bazel build //... --//:cupdlp_gpu --@rules_cuda//cuda:archs=sm_89
```

## [Building HiGHS with AMD GPU support](@id gpu-build-amd)

The native HiPDLP solver can also run on an AMD GPU using
[ROCm](https://rocm.docs.amd.com/) / HIP. This requires a ROCm
installation providing the HIP compiler and the hipBLAS and hipSPARSE
libraries. Make sure the HIP compiler is available by running

```
hipcc --version
```

ROCm 10.0 or newer is recommended (this is the version the HIP backend
is tested against). On such a ROCm, the supported GPU architectures are
`gfx908` and newer (e.g. `gfx908`/MI100, `gfx90a`/MI200,
`gfx942`/MI300, and recent RDNA cards).

Then build HiGHS, from the root directory, with

```
cmake -S. -Bbuild -DHIPDLP_HIP=ON
cmake --build build --parallel
```

CMake must be able to find ROCm. If it is not installed in the default
location, point it there, for example

```
export PATH=/opt/rocm/bin:$PATH
export CMAKE_PREFIX_PATH=/opt/rocm
```

By default the HIP device code is compiled for a generic set of GPU
architectures. To target the specific GPU on the build machine (which
also speeds up compilation and linking), set `CMAKE_HIP_ARCHITECTURES`
to its `gfx` target, for example

```
cmake -S. -Bbuild -DHIPDLP_HIP=ON -DCMAKE_HIP_ARCHITECTURES=gfx90a
```

You can find the `gfx` identifier of the installed GPU with `rocminfo`
(look for the `gfx` name, e.g. `gfx90a` for MI200-class cards or
`gfx942` for MI300). Multiple architectures may be given as a
semicolon-separated list, e.g. `-DCMAKE_HIP_ARCHITECTURES="gfx90a;gfx942"`.

By default the host C/C++ sources are compiled with the system compiler
(e.g. GCC) and only the HIP device code with ROCm's compiler. To build
the whole of HiGHS with the ROCm toolchain instead, for a uniform
Clang-based build, point CMake at `amdclang` / `amdclang++`

```
cmake -S. -Bbuild -DHIPDLP_HIP=ON \
  -DCMAKE_C_COMPILER=amdclang -DCMAKE_CXX_COMPILER=amdclang++
```

The HIP backend compiles the same HiPDLP source as the CUDA backend,
selecting the AMD implementation at build time. Once built, the solver
is selected at run time by setting the [__solver__](@ref
option-solver) option to "hipdlp".

To check the ROCm / HIP backend on the local machine, run the example
`call_highs_hipdlp` (also registered as the ctest
`cxx_examples_call_highs_hipdlp`), which solves a small LP with
`solver = "hipdlp"` and verifies the result. A successful run is a
quick end-to-end sanity check of the GPU backend.

To confirm the work is actually running on the GPU, watch `rocm-smi`
(for example `watch -n 0.1 rocm-smi`) while the solve runs and check
that GPU utilisation and memory usage rise.
