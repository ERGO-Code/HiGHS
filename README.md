# HiGHS - Linear optimization software

[![Build Status](https://github.com/ERGO-Code/HiGHS/workflows/build/badge.svg)](https://github.com/ERGO-Code/HiGHS/actions?query=workflow%3Abuild+branch%3Amaster)
[![Conan Center](https://img.shields.io/conan/v/highs)](https://conan.io/center/recipes/highs)
[![PyPi](https://img.shields.io/pypi/v/highspy.svg)](https://pypi.python.org/pypi/highspy)
[![PyPi](https://img.shields.io/pypi/dm/highspy.svg)](https://pypi.python.org/pypi/highspy)

- [HiGHS - Linear optimization software](#highs---linear-optimization-software)
  - [About HiGHS](#about-highs)
  - [Documentation](#documentation)
  - [Installation](#installation)
    - [Build from source using CMake](#build-from-source-using-cmake)
    - [Precompiled binaries](#precompiled-binaries)
  - [Interfaces](#interfaces)
      - [Python](#python)
  - [Reference](#reference)

## About HiGHS

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

### Build from source using CMake

HiGHS uses CMake as build system, and requires at least version 3.15. To generate build files in a new subdirectory called 'build', run:

```sh
    cmake -S . -B build
    cmake --build build
```
This installs the executable `bin/highs` and the library `lib/highs`.

To test whether the compilation was successful, change into the build directory and run

```sh
    ctest
```

HiGHS can read MPS files and (CPLEX) LP files, and the following command
solves the model in `ml.mps`

```sh
    highs ml.mps
```
HiGHS is installed using the command

```sh
    cmake --install build
```

with the optional setting of `--prefix <prefix>`, or the cmake option `CMAKE_INSTALL_PREFIX` if it is to be installed anywhere other than the default location.


As an alternative, HiGHS can be installed using the `meson` build interface:
``` sh
meson setup bbdir -Dwith_tests=True
meson test -C bbdir
```
_The meson build files are provided by the community and are not officially supported by the HiGHS development team._

### Precompiled binaries

Precompiled static executables are available for a variety of platforms at
https://github.com/JuliaBinaryWrappers/HiGHSstatic_jll.jl/releases

_These binaries are provided by the Julia community and are not officially supported by the HiGHS development team. If you have trouble using these libraries, please open a GitHub issue and tag `@odow` in your question._

See https://ergo-code.github.io/HiGHS/stable/installation/#Precompiled-Binaries.


## Interfaces

There are HiGHS interfaces for C, C#, FORTRAN, and Python in [HiGHS/src/interfaces](https://github.com/ERGO-Code/HiGHS/blob/master/src/interfaces), with example driver files in [HiGHS/examples](https://github.com/ERGO-Code/HiGHS/blob/master/examples). More on language and modelling interfaces can be found at https://ergo-code.github.io/HiGHS/stable/interfaces/other/.

We are happy to give a reasonable level of support via email sent to highsopt@gmail.com.

#### Python

The python package `highspy` is a thin wrapper around HiGHS and is available on [PyPi](https://pypi.org/project/highspy/). It can be easily installed via `pip` by running

```sh
$ pip install highspy
```

Alternatively, `highspy` can be built from source.  Download the HiGHS source code and run 

```sh
pip install . 
```
from the root directory. 

The HiGHS C++ library no longer needs to be separately installed. The python package `highspy` depends on the `numpy` package and `numpy` will be installed as well, if it is not already present. 

The installation can be tested using the small example [call_highs_from_python_highspy.py](https://github.com/ERGO-Code/HiGHS/blob/master/examples/call_highs_from_python_highspy.py).

The [Google Colab Example Notebook](https://colab.research.google.com/drive/1JmHF53OYfU-0Sp9bzLw-D2TQyRABSjHb?usp=sharing) also demonstrates how to call `highspy`.

## Reference

If you use HiGHS in an academic context, please acknowledge this and cite the following article.

Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: [10.1007/s12532-017-0130-5](https://link.springer.com/article/10.1007/s12532-017-0130-5)
