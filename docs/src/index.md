# HiGHS - High Performance Optimization Software

<!-- 
[![Build Status](https://github.com/ERGO-Code/HiGHS/workflows/build/badge.svg)](https://github.com/ERGO-Code/HiGHS/actions?query=workflow%3Abuild+branch%3Amaster) -->

![Build Status][fast_build_svg]][fast_build_link] 
[![Build Status][linux_build_svg]][linux_build_link] 
[![Build Status][macos_build_svg]][macos_build_link] 
[![Build Status][windows_build_svg]][windows_build_link] 
\
[![Conan Center](https://img.shields.io/conan/v/highs)](https://conan.io/center/recipes/highs)
\
[![PyPi](https://img.shields.io/pypi/v/highspy.svg)](https://pypi.python.org/pypi/highspy)
[![PyPi](https://img.shields.io/pypi/dm/highspy.svg)](https://pypi.python.org/pypi/highspy)
\
[![NuGet version](https://img.shields.io/nuget/v/Highs.Native.svg)](https://www.nuget.org/packages/Highs.Native)
[![NuGet download](https://img.shields.io/nuget/dt/Highs.Native.svg)](https://www.nuget.org/packages/Highs.Native)

[fast_build_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-fast.yml/badge.svg
[fast_build_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-fast.yml
[linux_build_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-linux.yml/badge.svg
[linux_build_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-linux.yml
[macos_build_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-macos.yml/badge.svg
[macos_build_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-macos.yml
[windows_build_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-windows.yml/badge.svg
[windows_build_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-windows.yml

!!! warning
    This HiGHS documentation is a work in progress.

HiGHS is software for the definition, modification and solution of large scale
sparse linear optimization models.

HiGHS is freely available from [GitHub](https://github.com/ERGO-Code/HiGHS)
under the MIT licence and has no third-party dependencies.

## Specification

HiGHS can solve linear programming (LP) models of the form:
```math
\begin{aligned}
\min                \quad & c^T\! x        \\
\textrm{subject to} \quad & L \le Ax \le U  \\
                          & l \le x \le u,
\end{aligned}
```
as well as mixed integer linear programming (MILP) models of the same form, for
which some of the variables must take integer values.

HiGHS also solves quadratic programming (QP) models, which contain an additional
objective term ``\frac{1}{2}x^T\! Q x``, where the Hessian matrix ``Q`` is
positive semi-definite. HiGHS cannot solve QP models where some of the variables
must take integer values.

Read the [Terminology](@ref) section for more details.

## Using HiGHS

HiGHS can be used as a standalone executable on Windows, Linux and MacOS. There
is also a C++11 library that can be used within a C++ project or, via its C, C#,
FORTRAN, Julia, and Python interfaces.

Get started by following [Install HiGHS](@ref).

## Overview

The standalone [Executable](@ref) allows models to be solved from
[MPS](https://en.wikipedia.org/wiki/MPS_(format)) or (CPLEX)
[LP](https://web.mit.edu/lpsolve/doc/CPLEX-format.htm) files, with full control
of the HiGHS run-time options, and the solution can be written to files in human
and computer-readable formats.

The HiGHS shared library allows models to be loaded, built and modified. It can
also be used to extract solution data and perform other operations relating to
the incumbent model. The basic functionality is introduced via a [`Guide`](@ref guide-basic),
with links to examples of its use in the `Python` interface `highspy`. This makes use of the C++
structures and enums, and is as close as possible to the native C++ library
calls. These can be studied via the [C++ header file](https://github.com/ERGO-Code/HiGHS/blob/master/src/Highs.h).

The C interface cannot make use of the C++ structures and enums, and its methods are documented [explicitly](@ref c-api).

## Solution algorithms

For LPs, HiGHS has implementations of both the revised simplex and interior
point methods. MIPs are solved by branch-and-cut, and QPs by active set.

## Citing HiGHS

If you use HiGHS in an academic context, please cite the following article:

Parallelizing the dual revised simplex method,
Q. Huangfu and J. A. J. Hall,
_Mathematical Programming Computation_, 10 (1), 119-142, 2018.
DOI: [10.1007/s12532-017-0130-5](https://link.springer.com/article/10.1007/s12532-017-0130-5)

## Performance benchmarks

The performance of HiGHS relative to some commercial and open-source simplex
solvers may be assessed via the Mittlemann benchmarks:

 * [LP (find primal-dual feasible point)](http://plato.asu.edu/ftp/lpfeas.html)
 * [LP (find optimal basic solution)](http://plato.asu.edu/ftp/lpopt.html)
 * [MILP benchmarks](http://plato.asu.edu/ftp/milp.html).

## Feedback

Your comments or specific questions on HiGHS would be greatly appreciated, so
please send an email to `highsopt@gmail.com` to get in touch with the
development team.
