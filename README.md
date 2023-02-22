# HiGHS - Linear optimization software

[![Build Status](https://github.com/ERGO-Code/HiGHS/workflows/build/badge.svg)](https://github.com/ERGO-Code/HiGHS/actions?query=workflow%3Abuild+branch%3Amaster)
[![PyPi](https://img.shields.io/pypi/v/highspy.svg)](https://pypi.python.org/pypi/highspy)
[![PyPi](https://img.shields.io/pypi/dm/highspy.svg)](https://pypi.python.org/pypi/highspy)

## Table of Contents

*   [About HiGHS](#about)
<!-- *   [Codemap](#codemap) -->
*   [Documentation](#documentation)
*   [Precompiled Binaries](#precompiled-binaries)
*   [Compilation](#compilation)
*   [Licence](#licence)
*   [Reference](#reference)

## About HiGHS 

HiGHS is a high performance serial and parallel solver for large scale sparse
linear optimization problems of the form

    Minimize (1/2) x^TQx + c^Tx subject to L <= Ax <= U; l <= x <= u

where Q must be positive semi-definite and, if Q is zero, there may be a requirement that some of the variables take integer values. Thus HiGHS can solve linear programming (LP) problems, convex quadratic programming (QP) problems, and mixed integer programming (MIP) problems. It is mainly written in C++, but also has some C. It has been developed and tested on various Linux, MacOS and Windows installations using both the GNU (g++) and Intel (icc) C++ compilers. Note that HiGHS requires (at least) version 4.9 of the GNU compiler. It has no third-party dependencies.

HiGHS has primal and dual revised simplex solvers, originally written by Qi Huangfu and further developed by Julian Hall. It also has an interior point solver for LP written by Lukas Schork, an active set solver for QP written by Michael Feldmeier, and a MIP solver written by Leona Gottwald. Other features have been added by Julian Hall and Ivet Galabova, who manages the software engineering of HiGHS and interfaces to C, C#, FORTRAN, Julia and Python.

Find out more about HiGHS at https://www.highs.dev.

Although HiGHS is freely available under the MIT license, we would be pleased to learn about users' experience and give advice via email sent to highsopt@gmail.com.

## Documentation

Documentation is available at https://ergo-code.github.io/HiGHS/. Executables

## Precompiled binaries

Precompiled static executables are available for a variety of platforms at:
https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases

_These binaries are provided by the Julia community and are not officially supported by the HiGHS development team. If you have trouble using these libraries, please open a GitHub issue and tag `@odow` in your question._

See https://ergo-code.github.io/HiGHS/cpp/get-started.html#Precompiled-executables.

## Compilation

HiGHS uses CMake as build system. First setup a build folder and call CMake as follows

    mkdir build
    cd build
    cmake -DFAST_BUILD=ON ..

Then compile the code using

    cmake --build . 

This installs the executable `bin/highs`.
The minimum CMake version required is 3.15.

To perform a quick test whether the compilation was successful, run

    ctest

HiGHS can read plain text MPS files and LP files and the following command
solves the model in `ml.mps`

    highs ml.mps

Language interfaces and further documentation
---------------------------------------------

There are HiGHS interfaces for C, C#, FORTRAN, and Python in HiGHS/src/interfaces, with example driver files in HiGHS/examples. 
We are happy to give a reasonable level of support via email sent to highsopt@gmail.com.

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

Webdemo
-------

Alternatively, HiGHS can directly be compiled into a single HTML file and used
in a browser. This requires `emscripten` to be installed from their website
(unfortunately, e.g. `sudo apt install emscripten` in Ubuntu Linux is broken):

    https://emscripten.org/docs/getting_started/downloads.html

Then, run

    sh build_webdemo.sh

This will create the file `build_webdemo/bin/highs.html`. For fast edit
iterations run

    find src app | entr -rs 'make -C build_webdemo highs; echo'

This will rebuild `highs.html` every time a source file is modified (e.g.
from Visual Studio Code).

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

## Licence

MIT License

Copyright (c) 2022 HiGHS


## Reference

If you use HiGHS in an academic context, please acknowledge this and cite the following article.
Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

http://www.maths.ed.ac.uk/hall/HuHa13/
