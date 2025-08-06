# Other Interfaces

!!! note
    Some of the interfaces listed on this page are not officially supported by
    the HiGHS development team and are contributed by the community.

## AMPL

HiGHS can be used via AMPL, see the [AMPL Documentation](https://dev.ampl.com/solvers/highs/index.html).


## GAMS

The interface is available at [GAMSlinks](https://github.com/coin-or/GAMSlinks/),
including [pre-build libraries](https://github.com/coin-or/GAMSlinks/releases).

## Javascript

 * HiGHS can be used from Javascript directly inside a web browser thanks to
   [highs-js](https://github.com/lovasoa/highs-js). See the [demo](https://lovasoa.github.io/highs-js/)
   and the [npm package](https://www.npmjs.com/package/highs).
 * Alternatively, HiGHS also has a [native Node.js](https://www.npmjs.com/package/highs-solver)
   interface.

## MATLAB

* [HiGHSMEX](https://github.com/savyasachi/HiGHSMEX) is a MATLAB interface for Windows to provide all the functionality of HiGHS, except the following: Reading problem data from a model file; Setting names for the rows and columns of the model, or setting name for the objective; Advanced features such as solution of systems using the current basis matrix.

* The HiGHS MIP and dual simplex LP solvers have been used _within_ MATLAB (so for all architectures) by default since release 2024a.

## R

 * An R interface is available through the [`highs` R package](https://cran.r-project.org/web/packages/highs/index.html).

## Rust

 * HiGHS can be used from rust through the [`highs` crate](https://crates.io/crates/highs).
   The rust linear programming modeler [**good_lp**](https://crates.io/crates/good_lp)
   supports HiGHS.
