## Download 

HiGHS can be cloned from the Edinburgh Group in Research and Optimization ([ERGO](https://www.maths.ed.ac.uk/ERGO/)) [GitHub repo](https://www.github.com/ERGO-COde/HiGHS).

``` bash
git clone https://github.com/ERGO-Code/HiGHS.git
```

Precompiled executables are available for a variety of platforms at https://github.com/JuliaBinaryWrappers/HiGHS_jll.jl/releases

Note that HiGHS is still pre-1.0, so the version numbers in the releases do not match versions of HiGHS in this repository.

For Windows users: if in doubt, choose the `x86_64-w64-mingw32-cxx11.tar.gz` file

For Mac users: choose the `x86_64-apple-darwin.tar.gz` file.

### Build HiGHS from source code

HiGHS uses CMake as a build system. The simplest setup is to create a build folder (within the folder into which HiGHS has been downloaded) and then build HiGHS within it. The name of the build folder is arbitrary but, assuming it is HiGHS/build, the full sequence of commands required is as follows

``` bash
cd HiGHS
mkdir build
cd build
cmake ..
cmake --build . 
```

This creates the executable `build/bin/highs`.

### Test Build

To perform a quick test to see whether the compilation was successful, run ctest from within the build folder.

``` bash
ctest 
```

## Install 

The default installation location may need administrative permissions. To install, after building and testing, run 

``` bash
cmake --install . 
```

To install in a specified installation directory run CMake with the CMAKE_INSTALL_PREFIX flag set: 

``` bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/highs_install ..
cmake --build .
cmake --install . 
```
