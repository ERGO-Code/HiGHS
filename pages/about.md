---
title: About HiGHS
permalink: /about/
---

## HiGHS

[![Build Status](https://travis-ci.org/ERGO-Code/HiGHS.svg?branch=master)](https://travis-ci.org/ERGO-Code/HiGHS)

HiGHS is a high performance serial and parallel solver for large-scale sparse linear programming (LP) models developed in C++11, with interfaces to C, C#, FORTRAN, Julia and Python.

HiGHS is freely available under the MIT licence.
We would be pleased to learn about users' experience and give advice, to get in touch with the team email [highsopt@gmail.com](mailto:highsopt@gmail.com).

More information about HiGHS can be found at the [website](http://www.highs.dev) and in the README.md in the [GitHub repo](https://www.github.com/ERGO-COde/HiGHS).

#### Release

Latest Release Version 1.0.

#### Requirements

##### OS
HiGHS can be used on Windows, Linux and MacOS.

##### Compilers

HiGHS has been tested and used with the following compilers: 

* GNU ` g++ ` version 6.0 or later, 
* Intel ` icc ` version
* Clang ` clang ` version

##### Dependencies

No third party sortware is required by HiGHS.

In order to build HiGHS from source CMake 3.15 is required. For precompiled executables and libraries please contact us at [highsopt@gmail.com](mailto:highsopt@gmail.com).

#### Specification

HiGHS solves large scale sparse linear programming (LP) problems of the form

``` 

    Minimize c^Tx subject to L <= Ax <= U; l <= x <= u

```

The HiGHS core solver implementes the dual revised simplex method in parallel C++ using OpenMP directives. An interior point solver is available as an optional feature.

##### Reference

Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

http://www.maths.ed.ac.uk/hall/HuHa13/

### Performance

The performance of HiGHS relative to some commercial and open-source simplex solvers may be assessed via the Mittelmann benchmarks on http://plato.asu.edu/ftp/lpsimp.html

### Documentation

Documentation can be found at todo: add url.
More information about HiGHS can be found at the [website](http://www.highs.dev).

#### Quick Guide

This is an example of building HiGHS on Linux using CMake. For more information refer to our [documentation](todo: add url).

If you have not already downloaded HiGHS, run

``` bash
git clone https://github.com/ERGO-Code/HiGHS.git
```

Enter the HiGHS directory and create a build subdirectory. Enter `build/` and generate Makefiles using CMake:

``` bash
    mkdir build
    cd build
    cmake ..
```

Then compile and test the HiGHS code using

``` bash
    make
    ctest
```

If the compilation and linking were successful, the executable `highs` is now located in the `bin/` subdirectory. To solve an LP model specified in `problem.mps` run 

``` bash
./bin/highs /path_to_file_dir/problem.mps
```

For usage and runtime option information see 

``` bash
./bin/highs --help
```

More details on HiGHS options and how to use the library see todo: docsy link

##### Parallel code

At the moment the parallel option is temporarily unavailable due to a large
refactoring in progress. Once the refactoring is complete the updated documentation will be published.

