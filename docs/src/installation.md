# Install HiGHS

## Install via a package manager

In most programming languages supported by HiGHS, such as `Julia`, `Python`, and
`Rust`, you can install HiGHS using the languages package manager. Conssult the
corresponding interface documentation for details.

## Precompiled Binaries

_These binaries are provided by the Julia community and are not officially
supported by the HiGHS development team. If you have trouble using these
libraries, please open a GitHub issue and tag `@odow` in your question._

Precompiled static executables are available for a variety of platforms at

 * [https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases](https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases)

Each download includes library files for linking to external projects and a
stand-alone executable.

To install a precompiled binary, download the appropriate `.tar.gz` file and
extract the executable located at `/bin/highs`.

 * For Windows users: if in doubt, choose the file ending in `x86_64-w64-mingw32-cxx11.tar.gz`
 * For M1 macOS users: choose the file ending in `aarch64-apple-darwin.tar.gz`
 * For Intel macOS users: choose the file ending in `x86_64-apple-darwin.tar.gz`

## Compile from source

HiGHS uses CMake as build system, and requires at least version 3.15. First
setup a build folder and call CMake as follows:

```bash
$ mkdir build
$ cd build
$ cmake -DFAST_BUILD=ON ..
```

Then compile the code using:

```bashs
$ cmake --build .
```

To test whether the compilation was successful, run

```bash
$ ctest
```

HiGHS is installed using the command

```bash
$ cmake --install .
```
