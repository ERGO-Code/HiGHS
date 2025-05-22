## Build changes

## Code changes

Fixed incorrect assertion in `HighsMipSolver::solutionFeasible()` (fixing [#2204](https://github.com/ERGO-Code/HiGHS/issues/2204))

As part of [#2251](https://github.com/ERGO-Code/HiGHS/issues/2251) cuPDLP-C will start from the incumbent solution in HiGHS. For a model that has been changed, the user must supply a starting solution via a call to `Highs::setSolution`

getColIntegrality now returns `HighsVarType::kContinuous` when `model_.lp_.integrality_` is empty (fixing [#2261](https://github.com/ERGO-Code/HiGHS/issues/2261))

Now ensuring that when solving a scaled LP with useful but unvalidated basis, it does not lose its scaling after validation, since the scaling factors will be applied to the solution (fixing [#2267](https://github.com/ERGO-Code/HiGHS/issues/2267))

By setting non-empty values of options `read_solution_file`, `read_basis_file`, `write_model_file` (with extension `.lp` or `.mps`), `write_solution_file`, `solution_file`, `write_basis_file`, these files will be read or written when calling `Highs::run()`. Hence options previously only available via the command line interface can be use (for example) by modelling languages that only call `Highs::run()` (fixing [#2269](https://github.com/ERGO-Code/HiGHS/issues/2269)).

Bug [#2273](https://github.com/ERGO-Code/HiGHS/issues/2273) fixed

As part of [#2286](https://github.com/ERGO-Code/HiGHS/pull/2286), the root of the HiGHS source files is now `highs/`, rather than `src/`

ZI rounding and shifting MIP primal heuristics have been added (see [#2287](https://github.com/ERGO-Code/HiGHS/pull/2287)). They are off by default, but can be activated by setting the options `mip_heuristic_run_zi_round` and `mip_heuristic_run_shifting` to be true. Options `mip_heuristic_run_rins`, `mip_heuristic_run_rens` and `mip_heuristic_run_root_reduced_cost` to run the RINS, RENS and rootReducedCost heuristics have been added. These are true by default, but setting them to be false can accelerate the MIP solver on easy problems.

Added `Highs_changeRowsBoundsByRange` to C API, fixing [#2296](https://github.com/ERGO-Code/HiGHS/issues/2296)

Corrected docstrings for `Highs_getReducedRow`, motivated by [#2312](https://github.com/ERGO-Code/HiGHS/issues/2312)

LP file reader no longer fails when there is no objective section. Fix is [#2316](https://github.com/ERGO-Code/HiGHS/pull/2316), but this exposes code quality issue [#2318](https://github.com/ERGO-Code/HiGHS/issues/2318)

Added a max scale factor (+1024) when scaling up coefficients in `preprocessBaseInequality` and `postprocessCut`. Fix is [#2337](https://github.com/ERGO-Code/HiGHS/pull/2337)

Corrected the bounds used in when strengthening coefficients in `HPresolve::rowPresolve`, fixing [#1517](https://github.com/ERGO-Code/HiGHS/issues/1517)

Fixed numerical error in `highs/mip/HighsCliqueTable.cpp`, closing [#2320](https://github.com/ERGO-Code/HiGHS/issues/2320)

Fixed bug in `highs/mip/HighsFeasibilityJump.cpp`, closing [#2331](https://github.com/ERGO-Code/HiGHS/issues/2331)

Tightened CMIR cuts, leading to small performance gain,  closing [#2333](https://github.com/ERGO-Code/HiGHS/issues/2333)

Scaling the tolerance in forcing row reduction to avoid use of rows with small coefficients and bounds,  closing [#2290](https://github.com/ERGO-Code/HiGHS/issues/2290)

Fixed bug when calculating a coefficient in one of the cuts in `separateImpliedBounds` in `highs/mip/HighsImplications.cpp`

Added `CSECTION` to the exceptions for keywords that are followed by text, and thus cannot be used as names of columns, RHS, ranges, bounds etc.

Introduced the following KKT error measures to `HighsInfo`: `num_relative_primal_infeasibilities`; `max_relative_primal_infeasibility`; `num_relative_dual_infeasibilities`; `max_relative_dual_infeasibility`; `num_primal_residual_errors`; `max_primal_residual_error`; `num_dual_residual_errors`; `max_dual_residual_error`; `num_relative_primal_residual_errors`; `max_relative_primal_residual_error`; `num_relative_dual_residual_errors`; `max_relative_dual_residual_error`; `num_complementarity_violations`; `max_complementarity_violation`; `primal_dual_objective_error.` The relative values are used to assess whether a solution deemed to be optimal by the first order LP solver `cuPDLP-C` or interior point solver `IPX` (without crossover) is truly acceptable. They also enable users to determine whether a solution corresponding to `HighsModelStatus::kUnknown` is acceptable to them as optimal. Also introduced options `complementarity_tolerance` used to assess whether the (relative) primal-dual objective error is acceptable, and `kkt_tolerance` which, if set to a value other than `kDefaultKktTolerance = 1e-7` is used as the tolerance for all the KKT error measures. The HiGHS documentation has been updated to reflect the new options and `HighsInfo` data, and logging messages indicate when KKT error measures are not satisfied, despite the solver considering the LP solution to be optimal.





Added a max scale factor (+1024) when scaling up coefficients in `preprocessBaseInequality` and `postprocessCut`. Fix is [#2337](https://github.com/ERGO-Code/HiGHS/pull/2337).

