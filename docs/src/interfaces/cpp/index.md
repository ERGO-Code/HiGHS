# [Getting started](@id cpp-getting-started)

HiGHS can be cloned from [GitHub](https://github.com/ERGO-Code/HiGHS) with the command

``` bash
git clone https://github.com/ERGO-Code/HiGHS.git
```

### Building HiGHS from source code

HiGHS uses CMake (minimum version 3.15) as a build system, and can use the following compilers

- Clang ` clang `
- GNU ` g++ `
- Intel ` icc `
- Microsoft ` MSVC `

The simplest setup is to create a build folder (within the folder into
which HiGHS has been downloaded) and then build HiGHS within it. The
name of the build folder is arbitrary but, assuming it is HiGHS/build,
the full sequence of commands required is as follows

``` bash
cd HiGHS
mkdir build
cd build
cmake -DFAST_BUILD=ON ..
cmake --build .
```

This creates the [executable](@ref Executable) `build/bin/highs`.

### Test build

To perform a quick test to see whether the compilation was successful, run `ctest` from within the build folder.

``` bash
ctest
```

### Install

The default installation location may need administrative
permissions. To install, after building and testing, run

``` bash
cmake --install .
```

To install in a specified installation directory run CMake with the
`CMAKE_INSTALL_PREFIX` flag set:

``` bash
cmake -DFAST_BUILD=ON -DCMAKE_INSTALL_PREFIX=/path/to/highs_install ..
cmake --build .
cmake --install .
```
