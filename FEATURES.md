## Build changes

## Code changes

Forcing column reduction now checks the bound on the column dual rather than whether the dual row activity is zero fixing [#2409](https://github.com/ERGO-Code/HiGHS/issues/2409)

Now handling correctly the case where an infeasible MIP has a feasible relaxation, so no ray is computed fixing [#2415](https://github.com/ERGO-Code/HiGHS/issues/2415)

Fixed minor bug exposed by [#2441](https://github.com/ERGO-Code/HiGHS/issues/2441) in `Highs::setSolution()` for a sparse user solution when the moidel is empty, and only clearing the dual data before solving with modified objective in `Highs::multiobjectiveSolve()` so that user-supplied solution is not cleared.

The irreducible infeasibility system (IIS) facility now detects infeasibility due to bounds on constraint activity values (implied by variable bounds) being incompatible with constraint bounds. A `kIisStrategyLight mode` for the `iis_strategy` option has been introduced so that only infeasibility due to incompatible variable/constraint bounds and constraint activity values is checked for. The LP corresponding to any known IIS is now formed and held as a data member of the `HighsIis` class. It can be obtained as a const reference using `Highs::getIisLp()`, and written to a file using `Highs::writeIisModel(const std::string& filename = "")`

Prompted by [#2463](https://github.com/ERGO-Code/HiGHS/issues/2463), the HiGHS solution and basis files now match data to any column and row names in the model, only assuming that the data are aligned with column and row indices if there are no names in the model. This requires a new version (v2) of the HiGHS basis file. Basis files from v1 are still read, but deprecated. Now, when writing out a model, basis or solution, column and row names are added to the model - previously they were created temporarily and inconsistently on the fly. If the model has existing names, then distinctive names are created to replace any blank names, but names with spaces or duplicate names yield an error status return.

