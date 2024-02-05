# cuPDLP-C observations

This directory contains files from [cuPDLP-C v0.3.0](https://github.com/COPT-Public/cuPDLP-C/tree/v0.3.0). Below are some issues experienced when integrating them into HiGHS.

## Preprocessing issue

The following line is not recognised by g++, 

> #if !(CUPDLP_CPU)

so I've had to replace all ocurrences by

> #ifndef CUPDLP_CPU

This yields a compiler warning about "extra tokens at end of #ifndef
directive" in the case of the following, but it's not a problem for
now, as CUPDLP_CPU is set

> #ifndef CUPDLP_CPU & USE_KERNELS

## cmake issues

CUPDLP_CPU and CUPDLP_DEBUG should both set when building. However, they are not recognised so are forced by the following lines in cupdlp_defs.h

#define CUPDLP_CPU
#define CUPDLP_DEBUG (1)

## Macro definitions

When definitions in [glbopts.h](https://github.com/ERGO-Code/HiGHS/blob/add-pdlp/src/pdlp/cupdlp/glbopts.h) such as the following are used in [CupdlpWrapper.cpp](https://github.com/ERGO-Code/HiGHS/blob/add-pdlp/src/pdlp/CupdlpWrapper.cpp) there is a g++ compiler error, because `typeof` isn't recognised

> #define CUPDLP_INIT(var, size)                                  \
  {                                                             \
    (var) = (typeof(var))malloc((size) * sizeof(typeof(*var))); \
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }

Hence there is a set of type-specific definitions in `CupdlpWrapper.h`, such as 

>#define cupdlp_init_double(var, size)\
   {\
     (var) = (double*)malloc((size) * sizeof(double));\
   }

## C methods not picked up by g++

Three methods
* `double infNorm(double *x, cupdlp_int n);`
* `void cupdlp_haslb(cupdlp_float *haslb, const cupdlp_float *lb, const cupdlp_float bound, const cupdlp_int len);`
* `void cupdlp_hasub(cupdlp_float *hasub, const cupdlp_float *ub, const cupdlp_float bound, const cupdlp_int len);`

are declared in [cupdlp_linalg.h](https://github.com/ERGO-Code/HiGHS/blob/add-pdlp/src/pdlp/cupdlp/cupdlp_linalg.h) and defined in [cupdlp_linalg.c](https://github.com/ERGO-Code/HiGHS/blob/add-pdlp/src/pdlp/cupdlp/cupdlp_linalg.c) but not picked up by g++. Hence duplicate methods are declared and defined in [CupdlpWrapper.h](https://github.com/ERGO-Code/HiGHS/blob/add-pdlp/src/pdlp/CupdlpWrapper.h) and [CupdlpWrapper.cpp](https://github.com/ERGO-Code/HiGHS/blob/add-pdlp/src/pdlp/CupdlpWrapper.cpp).


## Overflow when setting nIterLim

The line

settings->nIterLim = INFINITY;

in `cupdlp_utils.c` yields a compiler warning

overflow in conversion from ‘float’ to ‘cupdlp_int’ {aka ‘int’} changes value from ‘+Inff’ to ‘2147483647’

and results in non-deterministic behaviour. If nothing else, `nIterLim` is sometimes negative!

Fixed by introducing the following to glbopts.h, and using it to set nIterLim

#define I_INFINITY 2147483647

## Values of row duals

Dual values returned from cuPDLP-c seem always to be non-negative, even if they correspond to a pure upper-bounded constraint that has been negated. Since `PDHG_PostSolve` converts the solution to the problem solved by cuPDLP-c into a solution for the original problem, "un-permuting" `y` according to the reording of the constraints, it should negate the duals for pure upper-bounded constraints.

## Problem with sys/time.h

The HiGHS branch add-pdlp compiles and runs fine on @jajhall's Linux machine, but CI tests on GitHub fail utterly due to `sys/time.h` not being found. Since HiGHS won't be using the cuPDLP-c timing, this can be commented out using a compiler directive.

## Handling infeasible or unbounded problems

cuPDLP-c fails to terminate with the infeasible and unbounded LPs in unit tests `pdlp-infeasible-lp` and `pdlp-unbounded-lp` in `highs/check/TestPdlp.cpp`. In both cases the primal and dual step sizes grow steadily - eventually heading to infinity. Presumably, once they have reached a tolerance, cuPDLP-c should terminate so that infeasibility and unboundedness can be deduced according to whether the current iterate is primal feasible (as it is for `pdlp-unbounded-lp`).

## To be done

- Remove cuPDLP-c timing using a compiler directive
- Make cuPDLP-c less chatty
- Create HiGHS options to feed cuPDLP-c
- Return iteration count etc from cuPDLP-c





