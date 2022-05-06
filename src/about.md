## Specification

HiGHS solves large scale sparse linear programming (LP) problems of the form
```math
\textrm{min} \qquad c^Tx \qquad \textrm{subject to} \qquad L <= Ax <= U; \qquad l <= x <= u.
```

The HiGHS core solver implementes the dual revised simplex method in parallel C++. An interior point solver is available as an optional feature.

## OS
HiGHS can be used on Windows, Linux and MacOS.

## Compilers

HiGHS can be used with the following compilers:

- Clang ` clang `
- GNU ` g++ ` 
- Intel ` icc `

## Dependencies

No third party sortware is required by HiGHS, except for the Threads library.

In order to build HiGHS from source CMake 3.15 is required. For precompiled executables and libraries please contact us at [highsopt@gmail.com](mailto:highsopt@gmail.com).

## Reference

Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

http://www.maths.ed.ac.uk/hall/HuHa13/

### Performance

The performance of HiGHS relative to some commercial and open-source simplex solvers may be assessed via the Mittelmann benchmarks on http://plato.asu.edu/ftp/lpsimp.html
