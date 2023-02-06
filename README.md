# HiGHS - Linear optimization software

[![Build Status](https://github.com/ERGO-Code/HiGHS/workflows/build/badge.svg)](https://github.com/ERGO-Code/HiGHS/actions?query=workflow%3Abuild+branch%3Amaster)
[![PyPi](https://img.shields.io/pypi/v/highspy.svg)](https://pypi.python.org/pypi/highspy)
[![PyPi](https://img.shields.io/pypi/dm/highspy.svg)](https://pypi.python.org/pypi/highspy)

HiGHS is a high performance serial and parallel solver for large scale sparse
linear optimization problems of the form

    Minimize (1/2) x^TQx + c^Tx subject to L <= Ax <= U; l <= x <= u

where Q must be positive semi-definite and, if Q is zero, there may be a requirement that some of the variables take integer values. Thus HiGHS can solve linear programming (LP) problems, convex quadratic programming (QP) problems, and mixed integer programming (MIP) problems. It is mainly written in C++, but also has some C. It has been developed and tested on various Linux, MacOS and Windows installations using both the GNU (g++) and Intel (icc) C++ compilers. Note that HiGHS requires (at least) version 4.9 of the GNU compiler. It has no third-party dependencies.

HiGHS has primal and dual revised simplex solvers, originally written by Qi Huangfu and further developed by Julian Hall. It also has an interior point solver for LP written by Lukas Schork, an active set solver for QP written by Michael Feldmeier, and a MIP solver written by Leona Gottwald. Other features have been added by Julian Hall and Ivet Galabova, who manages the software engineering of HiGHS and interfaces to C, C#, FORTRAN, Julia and Python.

Although HiGHS is freely available under the MIT license, we would be pleased to learn about users' experience and give advice via email sent to highsopt@gmail.com.

Reference
---------
If you use HiGHS in an academic context, please acknowledge this and cite the following article.
Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

http://www.maths.ed.ac.uk/hall/HuHa13/

Wikipedia
---------

The project has an entry on Wikipedia: https://en.wikipedia.org/wiki/HiGHS_optimization_solver

Documentation
-------------

The rest of this file gives brief documentation for HiGHS. Comprehensive documentation is available via https://www.highs.dev.

Download
--------

Precompiled static executables are available for a variety of platforms at:
https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases

_These binaries are provided by the Julia community and are not officially supported by the HiGHS development team. If you have trouble using these libraries, please open a GitHub issue and tag `@odow` in your question._

**Installation instructions**

To install, download the appropriate file and extract the executable located at `/bin/highs`.

 * For Windows users: if in doubt, choose the file ending in `x86_64-w64-mingw32.tar.gz`
 * For M1 macOS users: choose the file ending in `aarch64-apple-darwin.tar.gz`
 * For Intel macOS users: choose the file ending in `x86_64-apple-darwin.tar.gz`

**Shared libaries**

For advanced users, precompiled executables using shared libraries are available for a variety of platforms at:
https://github.com/JuliaBinaryWrappers/HiGHS_jll.jl/releases.

Similar download instructions apply.

 * These files link against `libstdc++`. If you do not have one installed, download the platform-specific libraries from: 
   https://github.com/JuliaBinaryWrappers/CompilerSupportLibraries_jll.jl/releases/tag/CompilerSupportLibraries-v0.5.1%2B0
   and copy all the libraries into the same folder as the `highs` executable.
 * Unless using the FORTRAN interface, any of versions libgfortran3-libgfortran5 should work.
   If in doubt, Windows users should choose the `x86_64-w64-mingw32-libgfortran5.tar.gz`.

Compilation
-----------

HiGHS uses CMake as build system. First setup
a build folder and call CMake as follows

    mkdir build
    cd build
    cmake ..

Then compile the code using

    make

This installs the executable `bin/highs`.
The minimum CMake version required is 3.15.

Testing
-------

To perform a quick test whether the compilation was successful, run

    ctest

Run-time options
----------------

In the following discussion, the name of the executable file generated is
assumed to be `highs`.

HiGHS can read plain text MPS files and LP files and the following command
solves the model in `ml.mps`

    highs ml.mps

HiGHS options
-------------
Usage:
    highs [OPTION...] [file]
    
      --model_file arg        File of model to solve.
      --presolve arg          Presolve: "choose" by default - "on"/"off" are alternatives.
      --solver arg            Solver: "choose" by default - "simplex"/"ipm" are alternatives.
      --parallel arg          Parallel solve: "choose" by default - "on"/"off" are alternatives.
      --run_crossover arg     Run crossover after IPM: "on" by default - "choose"/"off" are alternatives.
      --time_limit arg        Run time limit (seconds - double).
      --options_file arg      File containing HiGHS options.
      --solution_file arg     File for writing out model solution.
      --write_model_file arg  File for writing out model.
      --random_seed arg       Seed to initialize random number generation.
      --ranging arg           Report cost, bound, RHS and basic solution ranging in any solution file: "off" by default - "on" is alternatives.
      --read_solution_file    File of solution to be read
      
      --version               Print version number      
  -h, --help                  Print help.
  
  Note:
  
  * If the file defines a quadratic term in the objective (so the problem is a QP or MIQP) and "simplex" or "ipm" is selected for the solver option, then the quadratic term will be ignored.
  * If the file constrains some variables to take integer values and defines a quadratic term in the objective, then the problem is MIQP and cannot be solved by HiGHS

Language interfaces and further documentation
---------------------------------------------

There are HiGHS interfaces for C, C#, FORTRAN, and Python in HiGHS/src/interfaces, with example driver files in HiGHS/examples. 
Documentation is availble via https://www.highs.dev/, and we are happy to give a reasonable level of support via
email sent to highsopt@gmail.com.

Parallel code
-------------

Parallel computation within HiGHS is limited to the dual simplex solver.
However, performance gain is unlikely to be significant at present. 
For the simplex solver, at best, speed-up is limited to the number of memory channels, rather than the number of cores. 

HiGHS will identify the number of available threads at run time, and restrict their use to the value of the HiGHS option `threads`.

If run with `threads=1`, HiGHS is serial. The `--parallel` run-time
option will cause the HiGHS parallel dual simplex solver to run in serial. Although this
could lead to better performance on some problems, performance will typically be
diminished.

If multiple threads are available, and run with `threads>1`, HiGHS will use multiple threads. 
Although the best value will be problem and architecture dependent, for the simplex solver `threads=8` is typically a
good choice. 
Although HiGHS is slower when run in parallel than in serial for some problems, it can be faster with multiple threads.

HiGHS Library
-------------

HiGHS is compiled in a shared library. Running

`make install`

from the build folder installs the library in `lib/`, as well as all header files in directories rooted at `include/highs/`. For a custom installation in `install_folder` run

`cmake -DCMAKE_INSTALL_PREFIX=install_folder ..`

and then

`make install`

To use the library from a CMake project use

`find_package(HiGHS)`

and add the correct path to HIGHS_DIR.

Compiling and linking without CMake
-----------------------------------

An executable defined in the file `use_highs.cpp` (for example) is linked with the HiGHS library as follows. After running the code above, compile and run with

`g++ -o use_highs use_highs.cpp -I install_folder/include/highs -L install_folder/lib/ -lhighs`

`LD_LIBRARY_PATH=install_folder/lib/ ./use_highs`

Interfaces
==========

Julia
-----

A Julia interface is available at https://github.com/jump-dev/HiGHS.jl.

Rust
----

HiGHS can be used from rust through the [`highs` crate](https://crates.io/crates/highs). The rust linear programming modeler [**good_lp**](https://crates.io/crates/good_lp) supports HiGHS. 

R
------

An R interface is available through the [`highs` R package](https://cran.r-project.org/package=highs).

Javascript
----------

HiGHS can be used from javascript directly inside a web browser thanks to [highs-js](https://github.com/lovasoa/highs-js). See the [demo](https://lovasoa.github.io/highs-js/) and the [npm package](https://www.npmjs.com/package/highs).

Node.js
-------

HiGHS has a [native Node.js](https://www.npmjs.com/package/highs-solver) interface.

C#
--

Here are observations on calling HiGHS from C#:

- [highs_csharp_api.cs](https://github.com/ERGO-Code/HiGHS/blob/master/src/interfaces/highs_csharp_api.cs) contains all the PInvoke you need. Copy it into your C# project.
- Make sure that the native HiGHS library (highs.dll, libhighs.dll, libhighs.so, ... depending on your platform) can be found at runtime. How to do this is platform dependent, copying it next to your C# executable should work in most cases. You can use msbuild for that. At least on linux installing HiGHS system wide should work, too.
- Make sure that all dependencies of the HiGHS library can be found, too. E.g. if HiGHS was build using Visual C++ make sure that the MSVCRuntime is installed on the machine you want to run your application on.
- Depending on the name of your HiGHS library it might be necessary to change the constant "highslibname", see [document](https://learn.microsoft.com/en-us/dotnet/standard/native-interop/cross-platform) on writing cross platform P/Invoke code if necessary.
- Call the Methods in highs_csharp_api.cs and have fun with HiGHS.

This is the normal way to call plain old C from C# with the great simplification that you don't have to write the PInvoke declarations yourself.

Python
------

There are two ways to build the Python interface to HiGHS. 

__From PyPi__

HiGHS has been added to PyPi, so should be installable using the command 

pip install highspy

The installation can be tested using the example [minimal.py](https://github.com/ERGO-Code/HiGHS/blob/master/examples/minimal.py), yielding the output

    Running HiGHS 1.2.2 [date: 2022-09-04, git hash: 8701dbf19]
    Copyright (c) 2022 ERGO-Code under MIT licence terms
    Presolving model
    2 rows, 2 cols, 4 nonzeros
    0 rows, 0 cols, 0 nonzeros
    0 rows, 0 cols, 0 nonzeros
    Presolve : Reductions: rows 0(-2); columns 0(-2); elements 0(-4) - Reduced to empty
    Solving the original LP from the solution after postsolve
    Model   status      : Optimal
    Objective value     :  1.0000000000e+00
    HiGHS run time      :          0.00

or the more didactic [call_highs_from_python.py](https://github.com/ERGO-Code/HiGHS/blob/master/examples/call_highs_from_python.py). 

__Directly__

In order to build the Python interface, build and install the HiGHS
library as described above, ensure the shared library is in the
`LD_LIBRARY_PATH` environment variable, and then run

`pip install ./`

from `src/interfaces/highspy` (there should be a `setup.py` file there).

You may also require

* `pip install pybind11`
* `pip install pyomo`

The Python interface can then be tested as above
