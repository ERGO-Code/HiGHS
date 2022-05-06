# HiGHS - High Performance Optimization Software
[![Build Status](https://github.com/ERGO-Code/HiGHS/workflows/build/badge.svg)](https://github.com/ERGO-Code/HiGHS/actions?query=workflow%3Abuild+branch%3Amaster)

HiGHS is a high performance serial and parallel solver for
large-scale sparse linear programming (LP) and mixed-integer
programming (MIP) models developed in C++11, with interfaces to C, C#,
FORTRAN, Julia and Python.

HiGHS is freely available under the MIT licence. 

HiGHS can be used as a standalone executable on Windows, Linux and MacOS. There is
a C++11 library which can be used within a C++ project or, via one of
the interfaces, to a project written in other languages.

Installing HiGHS from source code requires CMake minimum version 3.15, but no other third-party utilities. 

Your comments or specific questions on HiGHS would be greatly appreciated, so please send an email to  [highsopt@gmail.com](mailto:highsopt@gmail.com) to get in touch with the team.

### LP and MIP 
Linear programming problems (LP) in general bounded form are defined as
```math
\textrm{min} \qquad c^Tx \qquad \textrm{subject to} \qquad L <= Ax <= U; \qquad l <= x <= u.
```
and mixed integer programming (MIP) problems of the same form, for which some of the variables must take integer values. 

More on the theory background can be found at [www.highs.dev](http://www.highs.info). 