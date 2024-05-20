# HiGHS CMake Build Instructions 

| OS       | C++   | Python   | .NET   |
|:-------- | :---: | :------: | :----: |
| Linux    | [![Status][linux_cpp_svg]][linux_cpp_link] | [![Status][linux_python_svg]][linux_python_link] | [![Status][linux_dotnet_svg]][linux_dotnet_link] |
| MacOS    | [![Status][macos_cpp_svg]][macos_cpp_link] | [![Status][macos_python_svg]][macos_python_link] | [![Status][macos_dotnet_svg]][macos_dotnet_link] |
| Windows  | [![Status][windows_cpp_svg]][windows_cpp_link] | [![Status][windows_python_svg]][windows_python_link] | [![Status][windows_dotnet_svg]][windows_dotnet_link] |

[linux_cpp_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/cmake-linux-cpp.yml/badge.svg
[linux_cpp_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/cmake-linux-cpp.yml
[macos_cpp_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/cmake-macos-cpp.yml/badge.svg
[macos_cpp_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/cmake-macos-cpp.yml
[windows_cpp_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/cmake-windows-cpp.yml/badge.svg
[windows_cpp_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/cmake-windows-cpp.yml

[linux_python_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-python-ubuntu.yml/badge.svg
[linux_python_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-python-ubuntu.yml
[macos_python_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-python-macos.yml/badge.svg
[macos_python_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-python-macos.yml
[windows_python_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-python-win.yml/badge.svg
[windows_python_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-python-win.yml

[linux_dotnet_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-csharp-ubuntu.yml/badge.svg
[linux_dotnet_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-csharp-ubuntu.yml
[macos_dotnet_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-csharp-macos.yml/badge.svg
[macos_dotnet_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-csharp-macos.yml
[windows_dotnet_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-csharp-win.yml/badge.svg
[windows_dotnet_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/test-csharp-win.yml

<!--# ?branch=main -->

## Introduction
<nav for="cmake"> |
<a href="#requirement">Requirement</a> |
<a href="#build">C++</a> |
<a href="#cmake-options">Options</a> |
<a href="#integrating-highs-in-your-cmake-project">Integration</a> |
</nav>

HiGHS can be built from source using CMake: <http://www.cmake.org/>. CMake works by generating native Makefiles or build projects that can be used in the compiler environment of your choice.

HiGHS can be built as a standalone project or it could be incorporated into an existing CMake project.

## Requirement
You'll need:

* `CMake >= 3.15`.
* A C++11 compiler

#### Supported compilers 

Here is a list of the supported compilers:

* Clang `clang`
* GNU `g++`
* Intel `icc`
* Microsoft `MSVC`

# Build

To build the C++ library and executable run

``` bash
cd HiGHS
cmake -S. -B build 
cmake --build build --parallel
```

This creates the [executable](@ref Executable) `build/bin/highs`. To perform a quick test to see whether the compilation was successful, run `ctest` from within the build folder.

``` bash
ctest 
```

On Windows, the configuration type must be specified:
``` bash
ctest -C Release
```

## Install

The default installation location may need administrative
permissions. To install, after building and testing, run

``` bash
cmake --install build 
```

form the root directory. 

To install in a specified installation directory run CMake with the
`CMAKE_INSTALL_PREFIX` flag set:

``` bash
cmake -S. -B build -DCMAKE_INSTALL_PREFIX=/path/to/highs_install 
cmake --build build --parallel
cmake --install build
```

# CMake Options

There are several options that can be passed to CMake to modify how the code
is built.<br>
To set these options and parameters, use `-D<Parameter_name>=<value>`.

All CMake options are passed at configure time, i.e., by running <br>
`cmake -S. -B<your_chosen_build_directory> -DOPTION_ONE=ON -DOPTION_TWO=OFF ...` <br>
before running `cmake --build <your_chosen_build_directory>`<br>

For example, to generate build files in a new
subdirectory called 'build', run:

```sh
cmake -S. -Bbuild 
```
and then build with:

```sh
cmake --build build
```

Following is a list of available options:
| CMake Option | Default Value | Note |
|:-------------|:--------------|:-----|
| `CMAKE_BUILD_TYPE` | Release | see CMake documentation [here](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html) |
| `BUILD_SHARED_LIBS` | ON(*) |  Build shared libraries (.so or .dyld). * OFF by default on Windows |
| `BUILD_CXX` | ON | Build C++ |
| `FORTRAN` | OFF | Build Fortran interface |
| `CSHARP` | OFF | Build CSharp wrapper |
| `BUILD_DOTNET` | OFF | Build .Net package |
| `PYTHON_BUILD_SETUP` | OFF | Build Python bindings. Called at `pip install` from pyproject.toml |
| `ZLIB` | ON | Use ZLIB if available |
| `ALL_TESTS` | OFF | Run unit tests and extended instance test set |

<!-- Following is a list of available options, for the full list run:

```sh
cmake -S. -Bbuild -LH
``` -->

HiGHS can be integrated into other CMake-based projects. 

# Integrating HiGHS in your CMake Project

If you already have HiGHS installed on your system, you can use `find_package()` to include HiGHS in your C++ CMake project. 

```
project(LOAD_HIGHS LANGUAGES CXX)

set(HIGHS_DIR path_to_highs_install/lib/cmake/highs)

find_package(HIGHS REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main highs::highs)
```

The line 
```
set(HIGHS_DIR path_to_highs_install/lib/cmake/highs)
```
adds the HiGHS installation path to `HIGHS_DIR`. This is equivalent to building this project with
```
cmake -DHIGHS_DIR=path_to_highs_install/lib/cmake/highs ..
```

Alternatively, if you wish to include the code of HiGHS within your project, FetchContent is also available as follows: 

```
project(LOAD_HIGHS LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(
    highs
    GIT_REPOSITORY "https://github.com/ERGO-Code/HiGHS.git"
    GIT_TAG        "latest"
)

FetchContent_MakeAvailable(highs)

add_executable(main call_from_cpp.cc)
target_link_libraries(main highs::highs)
```
