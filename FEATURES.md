## Build changes

## Code changes

Forcing column reduction now checks the bound on the column dual rather than whether the dual row activity is zero fixing [#2409](https://github.com/ERGO-Code/HiGHS/issues/2409)

Now handling correctly the case where an infeasible MIP has a feasible relaxation, so no ray is computed fixing [#2415](https://github.com/ERGO-Code/HiGHS/issues/2415)

Fixed minor bug exposed by [#2441](https://github.com/ERGO-Code/HiGHS/issues/2441) in `Highs::setSolution()` for a sparse user solution when the moidel is empty, and only clearing the dual data before solving with modified objective in `Highs::multiobjectiveSolve()` so that user-supplied solution is not cleared.

The irreducible infeasibility system (IIS) facility now detects infeasibility due to bounds on constraint activity values (implied by variable bounds) being incompatible with constraint bounds. A `kIisStrategyLight mode` for the `iis_strategy` option has been introduced so that only infeasibility due to incompatible variable/constraint bounds and constraint activity values is checked for. The LP corresponding to any known IIS is now formed and held as a data member of the `HighsIis` class. It can be obtained as a const reference using `Highs::getIisLp()`, and written to a file using `Highs::writeIisModel(const std::string& filename = "")`

Prompted by [#2463](https://github.com/ERGO-Code/HiGHS/issues/2463), the HiGHS solution and basis files now match data to any column and row names in the model, only assuming that the data are aligned with column and row indices if there are no names in the model. This requires a new version (v2) of the HiGHS basis file. Basis files from v1 are still read, but deprecated. Now, when writing out a model, basis or solution, column and row names are added to the model - previously they were created temporarily and inconsistentyly on the fly. If the model has existing names, then distinctive names are created to replace any blank names, but names with spaces or duplicate names yield an error status return.

Refactored strong branching to minimize duplicated code

Only for LPs is there a choice of solver. Previously, when setting the `solver` option to anything other than "choose", any incumbent model was solved as an LP, using that LP solver. This has caused confusiuon for users, and is unnecessary now that there is the `solve_relaxation` option. Now, if the incumbent model is a QP or MIP, it is solved as such (unless `solve_relaxation` is true for a MIP), and the value of the `solver` option only determines what solver is used to solve an LP. If the value of `solver` is "choose", then HiGHS will use what it expects to be the best solver for the problem; if value of `solver` is "ipm", then HiGHS will use what it expects to be the better IPM solver (of HiPO and IPX) for the problem; if value of `solver` is "hipo", then HiGHS will use the HiPO IPM solver (if available in the build); if value of `solver` is "ipx", then HiGHS will use the IPX IPM solver; if value of `solver` is "pdlp", then HiGHS will use the PDLP first-order solver. The option `mip_lp_solver` has been introduced to define which LP solver is used when solving LPs in the MIP solver for which an advanced basis is not known - typically the "root node" LP. Note that The PDLP solver cannot be used to solve such LPs, since it does not yield a basic solution. If an interior point solver fails to obtain a basic solution, the simplex solver will then be used. The option `mip_ipm_solver` has been introduced to define which IPM solver is used when solving LPs in the MIP solver for which IPM is mandatory - typically the analytic centre calculation. When LPs are to be solved by an IPM solver, the HiPO solver is used (if available in the build) unless IPX has been specified explicitly. 

As per [#2487](https://github.com/ERGO-Code/HiGHS/issues/2487), trivial heuristics now run before feasibility jump (FJ), and FJ will use any existing incumbent. FJ will clip any finite variable values in the incumbent to lower and upper bounds, and falls back to the existing logic (lower bound if finite, else upper bound if finite, else 0) for any infinite values in the incumbent.

Prompted by [#2460](https://github.com/ERGO-Code/HiGHS/issues/2460), the options `user_cost_scale` and `user_bound_scale` apply uniform (power-of-two) scaling to the objective and bounds of a model, and now respect the following restrictions
- For a MIP, column bounds cannot be scaled, so the scaling is achieved by scaling the cost and constraint matrix column
- For a QP, Hessian entries must be scaled down (up) when bounds are scaled up (down) so that all terms in the objective are scaled by a constant.

After HiGHS determines and logs the coefficient ranges and warns about extreme values, it recommends values of `user_cost_scale` and `user_bound_scale` if
- All the objective coefficients (bound values) are smaller than the "excessively small constant" `kExcessivelySmallObjectiveCoefficient` (`kExcessivelySmallBoundValue`) - both of which are 1e-4 - suggesting that they are scaled up so that the largest value becomes (just over) the excessively small constant. Since the smallest value can be arbitrarily small, scaling so that this becomes (just over) the small constant is inadvisable, as the largest value could then be scaled up to an extremely large value.
- All the objective coefficients (bound values) are larger than the "excessively large constant" `kExcessivelyLargeObjectiveCoefficient` (`kExcessivelyLargeBoundValue`) - both of which are 1e6 - suggesting that they are scaled down so that the largest value becomes (just over) the excessively large constant.

The recommended objective scaling is determined with respect to the recommended bound scaling
Users can obtain the recommended values of user_cost_scale and user_bound_scale by calling `Highs::getObjectiveBoundScaling`

The irreducible infeasibility system (IIS) facility now detects infeasibility due to bounds on constraint activity values being incompatible with constraint bounds. A `kIisStrategyLight` mode for the `iis_strategy` option has been introduced so that only infeasibility due to incompatible variable/constraint bounds and constraint activity values is checked for. The model corresponding to any known IIS is now formed and held as a data member of the `HighsIis` class. The `HighsIis` class is available via `highspy`, and its data members are available via the C API.

Prompted by [#2463](https://github.com/ERGO-Code/HiGHS/issues/2463), when HiGHS writes out a solution or basis for a model without column or row names, it creates names. This avoids a mis-match between the ordering of variables when such a model is written out as a .lp file, and then this and a solution or a basis is read in.

Prompted by [#2487](https://github.com/ERGO-Code/HiGHS/issues/2487), the feasibility jump MIP heuristic starts from any incumbent solution, allowing it to improve on an existing feasible solution.

Prompted by [#2528](https://github.com/ERGO-Code/HiGHS/issues/2528), the logging for the HiGHS interior point method (IPM) solvers and PDLP solver has been standardised, and logging has been added to the IPM solver during time-consuming computational phases.

Prompted by [#2557](https://github.com/ERGO-Code/HiGHS/issues/2557), `Highs::getFixedLp` has been added so that, after solving a MIP, the LP with discrete variables fixed at their optimal values can be formed, allowing it to be passed to HiGHS and solved as an LP. 

Prompted by [#2581](https://github.com/ERGO-Code/HiGHS/issues/2581), the QP example in `call_highs_from_csharp.cs` has been corrected

Prompted by [#2582](https://github.com/ERGO-Code/HiGHS/issues/2582), the C API constants are declared `static const` to prevent multiple definition linker errors

