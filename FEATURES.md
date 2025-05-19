## Build changes

## Code changes

Fixed incorrect assertion in `HighsMipSolver::solutionFeasible()` (fixing [#2204](https://github.com/ERGO-Code/HiGHS/issues/2204)))

getColIntegrality now returns `HighsVarType::kContinuous` when `model_.lp_.integrality_` is empty (fixing [#2261](https://github.com/ERGO-Code/HiGHS/issues/2261))

Now ensuring that when solving a scaled LP with useful but unvalidated basis, it does not lose its scaling after validation, since the scaling factors will be applied to the solution (fixing [#2267](https://github.com/ERGO-Code/HiGHS/issues/2267))

By setting non-empty values of options `read_solution_file`, `read_basis_file`, `write_model_file` (with extension `.lp` or `.mps`), `write_solution_file`, `solution_file`, `write_basis_file`, these files will be read or written when calling `Highs::run()`. Hence options previously only available via the command line interface can be use (for example) by modelling languages that only call `Highs::run()` (fixing [#2269](https://github.com/ERGO-Code/HiGHS/issues/2269)).

Bug [#2273](https://github.com/ERGO-Code/HiGHS/issues/2273)) fixed

As part of [#2286](https://github.com/ERGO-Code/HiGHS/pull/2286)), the root of the HiGHS source files is now `highs/`, rather than `src/`

ZI rounding and shifting MIP primal heuristics have been added (see [#2287](https://github.com/ERGO-Code/HiGHS/pull/2287)). They are off by default, but can be activated by setting the options `mip_heuristic_run_zi_round` and `mip_heuristic_run_shifting` to be true. Options `mip_heuristic_run_rins`, `mip_heuristic_run_rens` and `mip_heuristic_run_root_reduced_cost` to run the RINS, RENS and rootReducedCost heuristics have been added. These are true by default, but setting them to be false can accelerate the MIP solver on easy problems.

Added `Highs_changeRowsBoundsByRange` to C API, fixing [#2296](https://github.com/ERGO-Code/HiGHS/issues/2296))

Corrected docstrings for `Highs_getReducedRow`, motivated by [#2312](https://github.com/ERGO-Code/HiGHS/issues/2312))

LP file reader no longer fails when there is no objective section. Fix is [#2316](https://github.com/ERGO-Code/HiGHS/pull/2316)), but this exposes code quality issue [#2318](https://github.com/ERGO-Code/HiGHS/issues/2318))

Added a max scale factor (+1024) when scaling up coefficients in `preprocessBaseInequality` and `postprocessCut`. Fix is [#2337](https://github.com/ERGO-Code/HiGHS/pull/2337).

