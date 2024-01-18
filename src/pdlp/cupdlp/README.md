# cuPDLP-C observations

This directory contains files from [cuPDLP-C v0.3.0](https://github.com/COPT-Public/cuPDLP-C/tree/v0.3.0). Below are some issues expereinced when integrating them into HiGHS.

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






