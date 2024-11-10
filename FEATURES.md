## Build changes

## Code changes

When primal infeasiblity is detected in presolve, no dual ray is available so, previously, the `has_dual_ray` parameter of `Highs::getDualRay` returned false and that was it. Now, if a null pointer is not passed for `dual_ray_value`, `Highs::getDualRay` will compute a dual ray - at the cost of solving the feasiblility LP without presolve. The same is now true for `Highs::getPrimalRay`. `Highs::getDualUnboundednessDirection` has been introduced to determine the product between the constraint matrix and the dual ray, forcing the calculation of the latter if necessary. Once a dual ray is known for the incumbent model in HiGHS, subsequent calls to `Highs::getDualRay` and `Highs::getDualUnboundednessDirection` will be vastly cheaper

The method `Highs::getDualObjectiveValue` now exitsts to compute the dual objective value, returning `HighsStatus::kError` if it is not possible.

The method `Highs::getStandardFormLp` now exists to return the incumbent LP in standard form - overlooking any integrality or Hessian. To determine the sizes of the vectors of data, the method is called without specifying pointers to the data arrays.

Added documentation on the use of presolve when solving an incumbent model, and clarifying the use of the method `Highs::presolve`.

HiGHS will now read a `MIPLIB` solution file

Added time limit check to `HPresolve::strengthenInequalities`

Added `getColIntegrality` to `highspy`

Now computing the primal-dual integral, reporting it, and making it available as `HighsInfo::primal_dual_integral`

Trivial primal heuristics "all zero", "all lower bound", "all upper bound", and "lock point" added to the MIP solver














