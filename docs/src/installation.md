# Install HiGHS

## Compile from source

HiGHS uses CMake as build system, and requires at least version
3.15. Details about building from source using CMake can be found in `HiGHS/cmake/README.md`.

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
vcpkg install openblas[threads]
```

Note, that `[threads]` is required for HiPO.

On Windows, you also need to specify the path to OpenBLAS. If it was installed with vcpkg as suggested above, add the path to `vcpkg.cmake` to the CMake flags, e.g.
```
cmake -S. -B build -DHIPO=ON -DCMAKE_TOOLCHAIN_FILE="C:/vcpkg/scripts/buildsystems/vcpkg.cmake"
```

##### Path to BLAS

To specify explicitly which BLAS vendor to look for, `BLA_VENDOR` coud be set in CMake, e.g. `-DBLA_VENDOR=Apple` or `-DBLA_VENDOR=OpenBLAS`. Alternatively, to specify which BLAS library to use, set `BLAS_LIBRARIES` to the full path of the library e.g. `-DBLAS_LIBRARIES=/path_to/libopenblas.so`.


## Bazel build

Alternatively, building with Bazel is supported for Bazel-based projects. To build HiGHS, from the root directory, run

```
bazel build //...
```

## Install via a package manager

HiGHS can be installed using a package manager in the cases of
[`Julia`](@ref HiGHS.jl), [`Python`](@ref python-getting-started), [`CSharp`](@ref nuget) and [`Rust`](@ref Rust).

## Precompiled Binaries

From v1.13.0 onwards, precompiled static binaries are available at https://github.com/ERGO-Code/HiGHS/releases.

Additionally, there is one package containing shared libraries for Windows x64.

The `*-mit` binary packages contain HiGHS and are MIT-licenced.
The `*-apache` binary packages contain HiGHS with HiPO and are Apache-licenced, due to the licensing of the dependencies of HiPO. For more information, see [THIRD_PARTY_NOTICES.md](https://github.com/ERGO-Code/HiGHS/blob/master/THIRD_PARTY_NOTICES.md).

If you have any questions or requests for more platforms and binaries, please get in touch with us at hello@highs.dev.

To install a precompiled binary, download and extract the archive corresponding to your Operating System and architecture, the executable is located at `/bin/highs`.

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

### Bazel build with Cuda

Alternatively, for Bazel run

```
bazel build //... --//:cupdlp_gpu
```