
## [Fortran](@id fortran-api)

The interface is in
[`HiGHS/src/interfaces/highs_fortran_api.f90`]. Its
methods are simply bindings to the [`C` API](@ref c-api)

To include in the build, switch the Fortran CMake parameter on:
``` bash
cmake -DFORTRAN=ON ..
```