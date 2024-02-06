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


## Problem with sys/time.h

The HiGHS branch add-pdlp compiles and runs fine on @jajhall's Linux machine, but CI tests on GitHub fail utterly due to `sys/time.h` not being found. Since HiGHS won't be using the cuPDLP-c timing, this can be commented out using a compiler directive.

## Making cuPDLP-c less chatty

As a research code, `cuPDLP-c` naturally produces a lot of output. Within HiGHS it should produce less output, and it should be possible to make it run silently. The simplest way to do this is introduce a logging level parameter into `cuPDLP-c` that, when zero, yields no output, when 1 yields just summary logging at the end, and when 2 or more produces the logging that you would wish to see. I guess that this would be added to `CUPDLP_INT_USER_PARAM_INDEX`.

A related issue is the use of `fp` and `fp_sol`. HiGHS won't be using these, so I set them to null pointers. `cuPDLP-c` already doesnt print the solution if `fp_sol` is a null pointer, so (like me) I suggest that the call to `writeJson(fp, pdhg);` is conditional on `if (fp)`. 

## Handling infeasible or unbounded problems

cuPDLP-c now terminates with status `INFEASIBLE_OR_UNBOUNDED` for the infeasible and unbounded LPs in unit tests `pdlp-infeasible-lp` and `pdlp-unbounded-lp` in `highs/check/TestPdlp.cpp`. In the case of the unbounded LP, PDLP identifies a primal feasible point, so unboundedness can be deduced. This is done in `HighsSolve.cpp:131.

## Returning the iteration count

The cuPDLP-c iteration count is held in `pdhg->timers->nIter`, but `pdhg` is destroyed in `LP_SolvePDHG`, so add `cupdlp_int* num_iter` to the parameter list of this method.

## To be done

- Make cuPDLP-c less chatty





