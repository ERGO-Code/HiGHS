# HiGHS - Linear optimization software

[![Build Status](https://travis-ci.org/ERGO-Code/HiGHS.svg?branch=master)](https://travis-ci.org/ERGO-Code/HiGHS)

HiGHS is a high performance serial and parallel solver for large scale
sparse linear programming (LP) problems of the form

    Maximize c^Tx subject to L <= Ax <= U; l <= x <= u

It is written in C++ with OpenMP directives. It is based on the dual
revised simplex method implemented in HSOL, exploiting parallelism using either "parallel
minor iterations" (PAMI) or "single iteration parallelism" (SIP). A
full technical reference to PAMI and SIP is

Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: 10.1007/s12532-017-0130-5

http://www.maths.ed.ac.uk/hall/HuHa13/

HSOL was originally written by Qi Huangfu, with features such as
presolve, crash and advanced basis start added by Julian Hall and Ivet
Galabova.

HSOL has been developed and tested on various linux installations
using both the GNU (g++) and Intel (icc) C++ compilers.

Compilation
-----------

HiGHS uses CMake as build system. To compile the run you need to setup
a build directory and define your build configuration:

    mkdir build
    cd build
    cmake .. [add otional parameters here]

Afterwards you may compile the code wrt the set configuration using your
defined build generator, e.g. `make`:

    make

Useful options
--------------

Set custom options with `-D<option>=<value>` during the configuration step (`cmake ..`):

- `OLD_PARSER`: 
    on: Uses the original fixed-format MPS file reader
    off: Uses the new free-format MPS file reader (which requires Boost)
- `OPENMP`:
    on: Causes OpenMP to be used if available - cmake checks for this
- `SCIP_DEV`:
    on: Suppresses all output so that HiGHS runs clean within SCIP
- `HiGHSDEV`:
    on: Includes a lot of testing and development code which should not be used in "production" or when running optimized code
- `HiGHSRELEASE`:
    on: Defined when CMAKE_BUILD_TYPE=Release

Testing
-------

To perform a quick test whether the compilation was successful, run

    ctest

Run-time options
----------------

In the following discussion, the name of the executable file generated
is assumed to be `HiGHS`.

HiGHS can only read plain text MPS files, and the following command
solves the model in `ml.mps`

    highs -f ml.mps

Usage
-----

```
usage: highs [options] -f fName.mps 

Options:
    -p mode  : use presolve mode. Values:
             : Off On
    -c mode  : use crash mode to mode. Values:
             : Off LTSSF LTSSF1 LTSSF2 LTSSF3 LTSSF4 LTSSF5 LTSSF6 LTSSF7
    -e edWt  : set edge weight to edWt. Values:
             : Dan Dvx DSE DSE0 DSE1
    -s       : use option sip
    -m [cut] : use pami. Cutoff optional double value.
    -t fName : use pami with partition file fName
    -d       : debug mode on
```

Run-time options `-p` and `-s` direct HiGHS to use PAMI or SIP.

When compiled with the OpenMP directives invoked, the number of
threads used at run time is the value of the environment variable
`OMP_NUM_THREADS`. For example, to use HiGHS with PAMI and eight
threads to solve `ml.mps` execute

    export OMP_NUM_THREADS=8
    highs -m -f ml.mps

If `OMP_NUM_THREADS` is not set, either because it has not been set or
due to executing the command

    unset OMP_NUM_THREADS

then all available threads will be used.

Observations
------------

When compiled without the OpenMP directives, or if run with
`OMP_NUM_THREADS=1`, HiGHS is serial. The `-sip` run-time option will not
affect performance. The `-pami` run-time option will cause HiGHS to use
serial minor iterations and, although this could lead to better
performance on some problems, performance will typically be
diminished.

When compiled with the OpenMP directives and `OMP_NUM_THREADS>1` or
unset, HiGHS will use multiple threads if the `-pami` or `-sip` run-time
option is specified. If `OMP_NUM_THREADS` is unset, HiGHS will try to use
all available threads so performance may be very slow. Although the
best value will be problem and architecture dependent,
`OMP_NUM_THREADS=8` is typically a good choice. Although HiGHS is slower
when run in parallel than in serial for some problems, it is typically
faster, with the `-pami` option usually faster than the `-sip` option.
