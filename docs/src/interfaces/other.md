# Other Interfaces

!!! note
    Some of the interfaces listed on this page are not officially supported by
    the HiGHS development team and are contributed by the community.

## AMPL

HiGHS can be used via AMPL, see the [AMPL Documentation](https://dev.ampl.com/solvers/highs/index.html).


## GAMS

HiGHS is included in GAMS distributions from GAMS 40 onward, see the
[GAMS/HiGHS documentation](https://www.gams.com/latest/docs/S_HIGHS.html).
GAMS provides the HiGHS LP and MIP solvers, including the CPU and NVIDIA GPU
versions of PDLP.

## Javascript

 * HiGHS can be used from Javascript directly inside a web browser thanks to
   [highs-js](https://github.com/lovasoa/highs-js). See the [demo](https://lovasoa.github.io/highs-js/)
   and the [npm package](https://www.npmjs.com/package/highs).
 * Alternatively, HiGHS also has a [native Node.js](https://www.npmjs.com/package/highs-solver)
   interface. The last release was in August 2025, so it may not include recent
   versions of HiGHS.

## MATLAB

* [HiGHSMEX](https://github.com/savyasachi/HiGHSMEX) is a MATLAB interface that provides all the functionality of HiGHS, except the following: Reading problem data from a model file; Setting names for the rows and columns of the model, or setting name for the objective; Advanced features such as solution of systems using the current basis matrix.

  The interface is avalailable for Windows, MacOS and Linux, and has been tested on Windows and MacOS. Pre-built binaries (mex files) for Windows and MacOS are available in the repository, which also includes instructions for building from source in [README.md](https://github.com/savyasachi/HiGHSMEX/blob/main/README.md).

* The HiGHS MIP and dual simplex LP solvers have been used _within_ MATLAB (so for all architectures) by default since release 2024a.

## R

 * An R interface is available through the [`highs` R package](https://cran.r-project.org/web/packages/highs/index.html).

## Rust

 * HiGHS can be used from rust through the [`highs` crate](https://crates.io/crates/highs).
   The rust linear programming modeler [**good_lp**](https://crates.io/crates/good_lp)
   supports HiGHS.
