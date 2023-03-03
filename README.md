# HiGHS - Linear optimization software

[![Build Status](https://github.com/ERGO-Code/HiGHS/workflows/build/badge.svg)](https://github.com/ERGO-Code/HiGHS/actions?query=workflow%3Abuild+branch%3Amaster)
[![PyPi](https://img.shields.io/pypi/v/highspy.svg)](https://pypi.python.org/pypi/highspy)
[![PyPi](https://img.shields.io/pypi/dm/highspy.svg)](https://pypi.python.org/pypi/highspy)

## Table of Contents

*   [About HiGHS](#about-highs)
*   [Documentation](#documentation)
*   [Precompiled binaries](#precompiled-binaries)
*   [Compilation](#compilation)
*   [Interfaces](#interfaces)
*   [Python](#python)
*   [Example](#google-colab-example)
*   [Reference](#reference)

About HiGHS 
-----------

HiGHS is a high performance serial and parallel solver for large scale sparse
linear optimization problems of the form

$$ \min \quad \dfrac{1}{2}x^TQx + c^Tx \qquad \textrm{s.t.}~ \quad L \leq Ax \leq U; \quad l \leq x \leq u $$

where Q must be positive semi-definite and, if Q is zero, there may be a requirement that some of the variables take integer values. Thus HiGHS can solve linear programming (LP) problems, convex quadratic programming (QP) problems, and mixed integer programming (MIP) problems. It is mainly written in C++, but also has some C. It has been developed and tested on various Linux, MacOS and Windows installations. No third-party dependencies are required.

HiGHS has primal and dual revised simplex solvers, originally written by Qi Huangfu and further developed by Julian Hall. It also has an interior point solver for LP written by Lukas Schork, an active set solver for QP written by Michael Feldmeier, and a MIP solver written by Leona Gottwald. Other features have been added by Julian Hall and Ivet Galabova, who manages the software engineering of HiGHS and interfaces to C, C#, FORTRAN, Julia and Python.

Find out more about HiGHS at https://www.highs.dev.

Although HiGHS is freely available under the MIT license, we would be pleased to learn about users' experience and give advice via email sent to highsopt@gmail.com.

Documentation
-------------

Documentation is available at https://ergo-code.github.io/HiGHS/.

Precompiled binaries
--------------------

Precompiled static executables are available for a variety of platforms at
https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases

_These binaries are provided by the Julia community and are not officially supported by the HiGHS development team. If you have trouble using these libraries, please open a GitHub issue and tag `@odow` in your question._

See https://ergo-code.github.io/HiGHS/binaries.html.

Compilation
-----------

HiGHS uses CMake as build system, and requires at least version 3.15. First setup a build folder and call CMake as follows

    mkdir build
    cd build
    cmake -DFAST_BUILD=ON ..

Then compile the code using

    cmake --build . 

This installs the executable `bin/highs`.

To test whether the compilation was successful, run

    ctest

HiGHS can read MPS files and (CPLEX) LP files, and the following command
solves the model in `ml.mps`

    highs ml.mps
    
HiGHS is installed using the command

    cmake --install .

with the optional setting of `--prefix <prefix>  = The installation prefix CMAKE_INSTALL_PREFIX` if it is to be installed anywhere other than the default location.

Interfaces 
----------

There are HiGHS interfaces for C, C#, FORTRAN, and Python in [HiGHS/src/interfaces](https://github.com/ERGO-Code/HiGHS/blob/master/src/interfaces), with example driver files in [HiGHS/examples](https://github.com/ERGO-Code/HiGHS/blob/master/examples). More on language and modelling interfaces can be found at https://ergo-code.github.io/HiGHS/interfaces.html.

We are happy to give a reasonable level of support via email sent to highsopt@gmail.com.

Python
------

There are two ways to build the Python interface to HiGHS. 

__From PyPi__

HiGHS has been added to PyPi, so should be installable using the command 

    pip install highspy

The installation can be tested using the example [minimal.py](https://github.com/ERGO-Code/HiGHS/blob/master/examples/minimal.py), yielding the output

    Running HiGHS 1.5.0 [date: 2023-02-22, git hash: d041b3da0]
    Copyright (c) 2023 HiGHS under MIT licence terms
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

    pip install ./

from `src/interfaces/highspy` (there should be a `setup.py` file there).

You may also require

* `pip install pybind11`
* `pip install pyomo`

The Python interface can then be tested as above.

Google Colab Example
-----------------------------
The [Google Colab Example Notebook](https://colab.research.google.com/drive/1JmHF53OYfU-0Sp9bzLw-D2TQyRABSjHb?usp=sharing) demonstrates how to call HiGHS via the Python interface `highspy`.

Reference
---------

If you use HiGHS in an academic context, please acknowledge this and cite the following article.

Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: [10.1007/s12532-017-0130-5](https://link.springer.com/article/10.1007/s12532-017-0130-5)
