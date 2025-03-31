## Build changes

## Code changes

Fixed incorrect assertion in `HighsMipSolver::solutionFeasible()` (fixing #2204)

getColIntegrality now returns `HighsVarType::kContinuous` when `model_.lp_.integrality_` is empty (fixing #2261)

