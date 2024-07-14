## Build changes

## Code changes

Added `int64_t mip_total_lp_iterations` to `HighsCallbackDataOut` and modified accessor function

`Highs::writeSolution` and `Highs::writeBasis` now being done via `HighsIO` logging, so can be redirected to logging callback.

Introduced `const double kHighsUndefined` as value of undefined values in a user solution. It's equal to `kHighsInf`

Added `Highs::setSolution(const HighsInt num_entries, const HighsInt* index, const double* value);` to allow a sparse primal solution to be defined.


