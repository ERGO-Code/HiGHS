## Build changes

### HiGHS on nixpkgs

HiGHS now has a `flake.nix` to build the binary, allowing `nix` users to try it out

### Python build update

Highspy with setuptools from v1.7.0 only worked on Python 3.12
For v1.7.0 we have dropped setuptools and switched to scikit-build-core

### Windows versions

Fixed version info of shared library
Added version info to executable

## Code changes

Inserting `pdlp_iteration_count` into various structs (for v1.7.0) breaks the C API, so it has been moved to the end of those structs

`setBasis` has been added to `highspy`

`writePresolvedModel` has been added

Saved MIP solution pool is populated when presolve reduces MIP to empty

Compilation date has been removed improve build reproducibility. Methods to print compilation dates are deprecated

Logging and error return when user-supplied solution or basis is rejected on vector dimensions

Memory allocation errors in presolve are caught and `Highs::run()` returns `HighsStatus::kError` with `model_status_ = HighsModelStatus::kMemoryLimit`

QP solver logging is now neater and quieter

Any Hessian for the incumbent model is modified with zero entries when adding columns to the model, and rows/columns are removed when columns are deleted from the model.

Minor bug fix in MIP presolve

QP solver will now hot start given a basis and solution


