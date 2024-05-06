# Other Interfaces

!!! note
    Some of the interfaces listed on this page are not officially supported by
    the HiGHS development team and are contributed by the community.

## AMPL

HiGHS can be used via AMPL, see the [AMPL Documentation](https://dev.ampl.com/solvers/highs/index.html).

## C#

Here are observations on calling HiGHS from C#:

 * The file [`highs_csharp_api.cs`](https://github.com/ERGO-Code/HiGHS/blob/master/src/interfaces/highs_csharp_api.cs)
   contains all the PInvoke you need. Copy it into your C# project.
 * Make sure, that the native HiGHS library (`highs.dll`, `libhighs.dll`,
   `libhighs.so`, ... depending on your platform) can be found at runtime. How
   to do this is platform dependent, copying it next to your C# executable
   should work in most cases. You can use msbuild for that. On linux, installing
   HiGHS system wide should work.
 * Make sure that all dependencies of the HiGHS library can be found, too. For
   example, if HiGHS was build using `Visual C++` make sure that the
   `MSVCRuntime` is installed on the machine you want to run your application
   on.
 * Depending on the name of your HiGHS library, it might be necessary to change
   the constant "highslibname". See [document](https://learn.microsoft.com/en-us/dotnet/standard/native-interop/cross-platform)
   on writing cross platform P/Invoke code if necessary.
 * Call the Methods in `highs_csharp_api.cs` and have fun with HiGHS.

This is the normal way to call plain old C from C# with the great simplification
that you don't have to write the PInvoke declarations yourself.

## [Fortran](@id fortran-api)

The interface is in
[`highs_fortran_api.f90`](https://github.com/ERGO-Code/HiGHS/blob/master/src/interfaces/highs_fortran_api.f90). Its
methods are simply bindings to the [`C` API](@ref c-api)

To include in the build, switch the Fortran CMake parameter on:
```
cmake -DFORTRAN=ON ..
```

## GAMS

The interface is available at [GAMSlinks](https://github.com/coin-or/GAMSlinks/),
including [pre-build libraries](https://github.com/coin-or/GAMSlinks/releases).

## Javascript

 * HiGHS can be used from javascript directly inside a web browser thanks to
   [highs-js](https://github.com/lovasoa/highs-js). See the [demo](https://lovasoa.github.io/highs-js/)
   and the [npm package](https://www.npmjs.com/package/highs).
 * Alternatively, HiGHS also has a [native Node.js](https://www.npmjs.com/package/highs-solver)
   interface.

## R

 * An R interface is available through the [`highs` R package](https://cran.r-project.org/web/packages/highs/index.html).

## Rust

 * HiGHS can be used from rust through the [`highs` crate](https://crates.io/crates/highs).
   The rust linear programming modeler [**good_lp**](https://crates.io/crates/good_lp)
   supports HiGHS.
