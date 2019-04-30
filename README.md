# HiGHS - Linear optimization software

[![Build Status](https://travis-ci.org/ERGO-Code/HiGHS.svg?branch=master)](https://travis-ci.org/ERGO-Code/HiGHS)

HiGHS is a high performance serial and parallel solver for large scale sparse
linear programming (LP) problems of the form

    Maximize c^Tx subject to L <= Ax <= U; l <= x <= u

It is written in C++ with OpenMP directives. It is based on the dual revised
simplex method implemented in HSOL.

Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

http://www.maths.ed.ac.uk/hall/HuHa13/

HSOL was originally written by Qi Huangfu, with features such as presolve,
crash and advanced basis start added by Julian Hall and Ivet Galabova and
further work by Michael Feldmeier.

HSOL has been developed and tested on various linux installations
using both the GNU (g++) and Intel (icc) C++ compilers.

Compilation
-----------

HiGHS uses CMake as build system. To compile the run you need to setup
a build directory and call

    mkdir build
    cd build
    cmake ..

Then compile the code using

    make

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

      --file arg             Filename of LP to solve.
      --presolve arg         Use presolve: off by default.
      --crash arg            Use crash to start simplex: off by default.
      --simplex arg          Use simplex solver: on by default.
      --ipm arg              Use interior point method solver: off by
                             default.
      --parallel arg         Use parallel solve: off by default.
      --time_limit arg       Use time limit.
      --iteration_limit arg  Use iteration limit (integer).
      --options_file arg     File containing HiGHS options.
      --parser arg           Mps parser type: swap back to fixed format
                             parser.
  -h, --help                 Print help.


Parallel code
-------------

At the moment the parallel option is temporarily unavailable due to a large
refactoring in progress. This document will be updated once we have completed
the interface currently being developed.

In order to use OpenMP if available, set`-DOPENMP=ON` during the configuration
step (`cmake ..`).

When compiled with the parallel option on, the number of threads used at run
time is the value of the environment variable `OMP_NUM_THREADS`. For example,
to use HiGHS with eight threads to solve `ml.mps` execute

    export OMP_NUM_THREADS=8
    highs --parallel ml.mps

If `OMP_NUM_THREADS` is not set, either because it has not been set or due to
executing the command

    unset OMP_NUM_THREADS

then all available threads will be used.

If run with `OMP_NUM_THREADS=1`, HiGHS is serial. The `--parallel` run-time
option will cause HiGHS to use serial minor iterations and, although this
could lead to better performance on some problems, performance will typically be
diminished.

When compiled with the parallel option and `OMP_NUM_THREADS>1` or unset, HiGHS
will use multiple threads. If `OMP_NUM_THREADS` is unset, HiGHS will try to use
all available threads so performance may be very slow. Although the best value
will be problem and architecture dependent, `OMP_NUM_THREADS=8` is typically a
good choice. Although HiGHS is slower when run in parallel than in serial for
some problems, it is typically faster in parallel.

HiGHS Library
-------------

Highs is compiled in a shared library. Running

`make install`

installs the highs executable in the bin/ and the highs library in the
lib/ folder, as well as all header files in include/. For a custom
installation in `install_folder` run

`cmake -DCMAKE_INSTALL_PREFIX=install_folder ..`

and then

`make install`

To use the library from a cmake project use

`find_package(HiGHS)`

and add the correct path to HIGHS_DIR.

Compiling and linking without cmake
Suppose we want to link an executable defined in file `use_highs.cpp` with the
highs library. After running the code above compile and run with

`g++ -o use_highs use_highs.cpp -I install_folder/include/ -L install_folder/lib/ -lhighs`

`LD_LIBRARY_PATH=install_folder/lib/ ./use_highs`

Interfaces
----------

GAMS
----

Set custom options with `-D<option>=<value>` during the configuration step (`cmake ..`):

- `GAMS_ROOT`:
    path to GAMS system: enables building of GAMS interface

If build with GAMS interface, then HiGHS can be made available as solver
in GAMS by adding an entry for HiGHS to the file gmscmpun.txt in the GAMS
system directory (gmscmpnt.txt on Windows):
```
HIGHS 11 5 0001020304 1 0 2 LP RMIP
gmsgenus.run
gmsgenux.out
/path/to/libhighs.so his 1 1
```
OSI
---
- `OSI_ROOT`:
    path to COIN-OR/Osi build/install directory (OSI_ROOT/lib/pkg-config/osi.pc should exist)
