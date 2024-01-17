# Install HiGHS

## Install via a package manager

HiGHS can be installed using a package manager in the cases of
[`Julia`](@ref HiGHS.jl), [`Python`](@ref python-getting-started), and
[`Rust`](@ref Rust).

## Precompiled Binaries

_These binaries are provided by the Julia community and are not officially
supported by the HiGHS development team. If you have trouble using these
libraries, please open a GitHub issue and tag `@odow` in your question._

Precompiled static executables are available for a variety of platforms at

 * [https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases](https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases)

Multiple versions are available. Each version has the form `vX.Y.Z`. In
general, you should choose the most recent versinon.

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

## Compile from source

HiGHS uses CMake as build system, and requires at least version
3.15. After extracting HiGHS from
[GitHub](https://github.com/ERGO-Code/HiGHS), setup a build folder and
call CMake as follows:

```bash
$ mkdir build
$ cd build
$ cmake ..
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

This installs the library in `lib/`, as well as all header files in `include/highs/`. For a custom
installation in `install_folder` run

```bash
$ cmake -DCMAKE_INSTALL_PREFIX=install_folder .
```

and then

```bash
$ cmake --install .
```

To use the library from a CMake project use

`find_package(HiGHS)`

and add the correct path to HIGHS_DIR.

## Windows 

By default, CMake builds the debug version of the binaries. These are generated in a directory `Debug`. To build a release version, add the option `--config Release`

```bash
    cmake -S . -B build
    cmake --build build --config Release
```

It is also possible to specify a specific Visual studio version to build with:
```bash
    cmake -G "Visual Studio 17 2022" -S . -B build
    cmake --build build
```

When building under Windows, some extra options are available.  One is building a 32 bit version or a 64 bit version. The default build is 64 bit. To build 32 bit, the following commands can be used from the `HiGHS/` directory:

```bash
    cmake -A Win32 -S . -DFAST_BUILD=OFF -B buildWin32
    cmake --build buildWin32
```

Another thing specific for windows is the calling convention, particularly important for the HiGHS dynamic library (dll). The default calling convention in windows is cdecl calling convention, however, dlls are most often compiled with stdcall. Most applications which expect stdcall, can't access dlls with cdecl and vice versa. To change the default calling convention from cdecl to stdcall the following option can be added
```bash
    cmake -DSTDCALL=ON -S . -DFAST_BUILD=OFF -B build
    cmake --build build
```
An extra note. With the legacy `-DFAST_BUILD=OFF`, under windows the build dll is called `highs.dll` however the exe expects `libhighs.dll` so a manual copy of `highs.dll` to `libhighs.dll` is needed. Of course all above options can be combined with each other.

