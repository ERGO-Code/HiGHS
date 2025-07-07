## Build changes

## Code changes

Forcing column reduction now checks the bound on the column dual rather than whether the dual row activity is zero fixing [#2409](https://github.com/ERGO-Code/HiGHS/issues/2409)

Now handling correctly the case where an infeasible MIP has a feasible relaxation, so no ray is computed fixing [#2415](https://github.com/ERGO-Code/HiGHS/issues/2415)

Fixed minor bug exposed by [#2441](https://github.com/ERGO-Code/HiGHS/issues/2441) in Highs::setSolution() for a sparse user solution when the moidel is empty, and only clearing the dual data before solving with modified objective in Highs::multiobjectiveSolve() so that user-supplied solution is not cleared.


