# HiGHS - Linear optimization software

[![Build Status](https://github.com/ERGO-Code/HiGHS/workflows/build/badge.svg)](https://github.com/ERGO-Code/HiGHS/actions?query=workflow%3Abuild+branch%3Amaster)
[![Conan Center](https://img.shields.io/conan/v/highs)](https://conan.io/center/recipes/highs)
[![PyPi](https://img.shields.io/pypi/v/highspy.svg)](https://pypi.python.org/pypi/highspy)
[![PyPi](https://img.shields.io/pypi/dm/highspy.svg)](https://pypi.python.org/pypi/highspy)

## Table of Contents

- [HiGHS - Linear optimization software](#highs---linear-optimization-software)
  - [Table of Contents](#table-of-contents)
  - [About HiGHS](#about-highs)
  - [Documentation](#documentation)
  - [Installation](#installation)
    - [Precompiled binaries](#precompiled-binaries)
    - [Compilation](#compilation)
    - [Meson](#meson)
    - [Python](#python)
  - [Interfaces](#interfaces)
    - [Python](#python-1)
      - [From PyPi](#from-pypi)
      - [Build directly from Git](#build-directly-from-git)
      - [Testing](#testing)
      - [Google Colab Example](#google-colab-example)
  - [Reference](#reference)

## About HiGHS
-----------

HiGHS is a high performance serial and parallel solver for large scale sparse
linear optimization problems of the form

$$ \min \quad \dfrac{1}{2}x^TQx + c^Tx \qquad \textrm{s.t.}~ \quad L \leq Ax \leq U; \quad l \leq x \leq u $$

where Q must be positive semi-definite and, if Q is zero, there may be a requirement that some of the variables take integer values. Thus HiGHS can solve linear programming (LP) problems, convex quadratic programming (QP) problems, and mixed integer programming (MIP) problems. It is mainly written in C++, but also has some C. It has been developed and tested on various Linux, MacOS and Windows installations. No third-party dependencies are required.

HiGHS has primal and dual revised simplex solvers, originally written by Qi Huangfu and further developed by Julian Hall. It also has an interior point solver for LP written by Lukas Schork, an active set solver for QP written by Michael Feldmeier, and a MIP solver written by Leona Gottwald. Other features have been added by Julian Hall and Ivet Galabova, who manages the software engineering of HiGHS and interfaces to C, C#, FORTRAN, Julia and Python.

Find out more about HiGHS at https://www.highs.dev.

Although HiGHS is freely available under the MIT license, we would be pleased to learn about users' experience and give advice via email sent to highsopt@gmail.com.

## Documentation

Documentation is available at https://ergo-code.github.io/HiGHS/.

## Installation

There are various ways to install the HiGHS library. These are detailed below.

### Precompiled binaries
--------------------

Precompiled static executables are available for a variety of platforms at
https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases

_These binaries are provided by the Julia community and are not officially supported by the HiGHS development team. If you have trouble using these libraries, please open a GitHub issue and tag `@odow` in your question._

See https://ergo-code.github.io/HiGHS/stable/installation/#Precompiled-Binaries.

### Compilation
---------------

HiGHS uses CMake as build system, and requires at least version 3.15. First setup a build folder and call CMake as follows

    mkdir build
    cd build
    cmake ..

Then compile the code using

    cmake --build .

This installs the executable `bin/highs`.

As an alternative it is also possible to let `cmake` create the build folder and thus build everything from the HiGHS directory, as follows

    cmake -S . -B build
    cmake --build build


To test whether the compilation was successful, run

    ctest

HiGHS can read MPS files and (CPLEX) LP files, and the following command
solves the model in `ml.mps`

    highs ml.mps

HiGHS is installed using the command

    cmake --install .

with the optional setting of `--prefix <prefix>  = The installation prefix CMAKE_INSTALL_PREFIX` if it is to be installed anywhere other than the default location.

### Meson
-----

HiGHs can also use the `meson` build interface:

``` sh
meson setup bbdir -Dwith_tests=True
meson test -C bbdir
```


### Python
-----

Installing from PyPI through your Python package manager of choice (e.g., `pip`) will also 
install the HiGHS library if not already present. HiGHS is available as `highspy` on [PyPi](https://pypi.org/project/highspy/).

If `highspy` is not already installed, run:

```bash
$ pip install highspy
```

## Interfaces
There are HiGHS interfaces for C, C#, FORTRAN, and Python in [HiGHS/src/interfaces](https://github.com/ERGO-Code/HiGHS/blob/master/src/interfaces), with example driver files in [HiGHS/examples](https://github.com/ERGO-Code/HiGHS/blob/master/examples). More on language and modelling interfaces can be found at https://ergo-code.github.io/HiGHS/stable/interfaces/other/.

We are happy to give a reasonable level of support via email sent to highsopt@gmail.com.

### Python

There are two ways to install the Python interface. Building directly 
from Git assumes that you have already installed the HiGHS library. 
Installing from PyPI through your Python package manager of choice (e.g., `pip`) will also install the HiGHS library if not already present. 

#### From PyPi

HiGHS is available as `highspy` on [PyPi](https://pypi.org/project/highspy/).
This will not only install the Python interface, but also the HiGHS library 
itself.

If `highspy` is not already installed, run:

```bash
$ pip install highspy
```

#### Build directly from Git

In order to build the Python interface, build and install the HiGHS
library as described above, ensure the shared library is in the
`LD_LIBRARY_PATH` environment variable, and then run

    meson setup bbdir -Dwith_pybind11=True
    meson compile -C bbdir
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/bbdir
    pip install ./

from the HiGHS directory.

The Python interface can then be tested as below.

#### Testing

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

#### Google Colab Example

The [Google Colab Example Notebook](https://colab.research.google.com/drive/1JmHF53OYfU-0Sp9bzLw-D2TQyRABSjHb?usp=sharing) demonstrates how to call HiGHS via the Python interface `highspy`.


## Reference


If you use HiGHS in an academic context, please acknowledge this and cite the following article.

Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: [10.1007/s12532-017-0130-5](https://link.springer.com/article/10.1007/s12532-017-0130-5)
