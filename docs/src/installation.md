# Install HiGHS

## Compile from source

HiGHS uses CMake as build system, and requires at least version
3.15. Details about building from source using CMake can be found in `HiGHS/cmake/README.md`.

## HiGHS with HiPO

HiGHS does not have any external dependencies, however, the new interior point solver HiPO uses BLAS and Metis. At the moment HiPO is optional and can be enabled via CMake. To build HiPO, you need to have Metis and BLAS installed on your machine. Please follow the instructions below.

#### BLAS

On Linux, libblas and libopenblas are supported. We recomment libopenblas for its better performance, and it is found by default if available on the system. Install with

```
sudo apt update
sudo apt install libopenblas-dev
```

On MacOS no BLAS installation is required because HiPO uses [Apple Accelerate](https://developer.apple.com/accelerate/) by default.

On Windows, OpenBLAS is required. It could be installed via [vcpkg](https://learn.microsoft.com/en-us/vcpkg/get_started/overview) with

```
vcpkg install openblas[threads]
```
Note, that `[threads]` is required for HiPO.

To specify explicitly which BLAS vendor to look for, `BLA_VENDOR` coud be set in CMake, e.g. `-DBLA_VENDOR=Apple` or `-DBLA_VENDOR=OpenBLAS`. Alternatively, to specify which BLAS library to use, set `BLAS_LIBRARIES` to the full path of the library e.g. `-DBLAS_LIBRARIES=/path_to/libopenblas.so`.

#### Metis
There are some known issues with Metis so the recommented version is in [this fork](https://github.com/galabovaa/METIS/tree/510-ts), branch 510-ts. This is version 5.10 with several patches for more reliable build and execution. Clone the repository with
```
git clone https://github.com/galabovaa/METIS.git
cd METIS
git checkout 510-ts
```

Then build with
```
cmake -S. -B build
-DGKLIB_PATH=/path_to_METIS_repo/GKlib
-DCMAKE_INSTALL_PREFIX=path_to_installs_dir
cmake --build build
cmake --install build
```

On Windows, do not forget to specify configuration type
```
cmake --build build --config Release
```

### HiPO

To install HiPO, on Linux and MacOS, run
```
cmake -S. -B build -DHIPO=ON -DMETIS_ROOT=path_to_installs_dir
```
On Windows, you also need to specify the path to OpenBLAS. If it was installed with vcpkg as suggested above, add the path to `vcpkg.cmake` to the CMake flags, e.g.
```
-DCMAKE_TOOLCHAIN_FILE="C:/vcpkg/scripts/buildsystems/vcpkg.cmake"
```

## Bazel build

Alternatively, building with Bazel is supported for Bazel-based projects. To build HiGHS, from the root directory, run

```
bazel build //...
```

## Install via a package manager

HiGHS can be installed using a package manager in the cases of
[`Julia`](@ref HiGHS.jl), [`Python`](@ref python-getting-started), [`CSharp`](@ref nuget) and [`Rust`](@ref Rust).

## Precompiled Binaries

_These binaries are provided by the Julia community and are not officially
supported by the HiGHS development team. If you have trouble using these
libraries, please open a GitHub issue and tag `@odow` in your question._

Precompiled static executables are available for a variety of platforms at

 * [https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases](https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases)

Multiple versions are available. Each version has the form `vX.Y.Z`. In
general, you should choose the most recent version.

To install a precompiled binary, download the appropriate `HiGHSstatic.vX.Y.Z.[platform-string].tar.gz`
file and extract the executable located at `/bin/highs`.

Do not download the file starting with `HiGHSstatic-logs`. These files contain
information from the automated compilation system. Click "Show all N assets"
to see more files.

### Platform strings

The GitHub releases contain precompiled binaries for a number of different
platforms. These are indicated by the platform-specific string in each
filename.

 * For Windows users: choose the file ending in `x86_64-w64-mingw32-cxx11.tar.gz`
 * For M1 macOS users: choose the file ending in `aarch64-apple-darwin.tar.gz`
 * For Intel macOS users: choose the file ending in `x86_64-apple-darwin.tar.gz`

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