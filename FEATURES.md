## Build changes

## Code changes

Added `int64_t mip_total_lp_iterations` to `HighsCallbackDataOut` and modified accessor function

`Highs::writeSolution` and `Highs::writeBasis` now being done via `HighsIO` logging, so can be redirected to logging callback.

Introduced `const double kHighsUndefined` as value of undefined values in a user solution. It's equal to `kHighsInf`

Added `Highs::setSolution(const HighsInt num_entries, const HighsInt* index, const double* value);` to allow a sparse primal solution to be defined. When a MIP is solved to do this, the value of (new) option `mip_max_start_nodes` is used for `mip_max_nodes` to avoid excessive cost

Added options `write_presolved_model_to_file` and `write_presolved_model_file` so that presolved model can be written via a command line option

Added `Highs::feasibilityRelaxation` to solve the problem of minimizing a (possibly weighted) sum of (allowable) infeasibilities in an LP/MIP.



