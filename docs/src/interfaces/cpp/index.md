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

Instructions for building HiGHS from source code are in `HiGHS/cmake/README.md`.

The simplest setup is to build HiGHS in a build directory within the root directory. The
name of the build folder is arbitrary but, assuming it is `build`,
the sequence of commands is as follows

``` bash
cd HiGHS
cmake -S. -B build 
cmake --build build --parallel
```