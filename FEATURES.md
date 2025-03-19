## Build changes

Added code coverage report.

Replaced command line parsing library with CLI11. Removed C++17 reference with cxxopts, which is no longer in HiGHS. 

## Code changes

Any LP offset is communicated to the IPM solver, and used in logging and primal/dual objective calculations. 

If there is a valid basis when Highs::run() is called, presolve isn't skipped unless the solver option is "simplex" or "choose" (when simplex will always be chosen if there is an advanced basis).

Added basis solve methods to highspy

Added methods to get primal/dual ray and dual unboundedness direction to highspy

When a presolved LP has model status kUnknown, rather than returning this to the user, it performs postsolve and then uses the basis to solve the original LP

Fixed bug in presolve when pointers stored in HighsMatrixSlice get invalidated when the coefficient matrix is reallocated (e.g. when non-zeros are added in HPresolve::addToMatrix)

Primal and dual residual tolerances - applied following IPM or PDLP solution - now documented as options

Highs::getCols (Highs::getRows) now runs in linear time if the internal constraint matrix is stored column-wise (row-wise). Added ensureColwise/Rowwise to the Highs class, the C API and highspy so that users can set the internal constraint matrix storage orientation

When columns and rows are deleted from the incumbent LP after a basic solution has been found, HiGHS no longer invalidates the basis. Now it maintains the basic and nonbasic status of the remaining variables and constraints. When the model is re-solved, this information is used to construct a starting basis.

Fixed bugs in presolve

When running from the command line, changes to default option values are reported

Added callback to allow users to supply integer feasible solutions to the MIP solver during execution

Bug fix for primal heuristics in the MIP solver

Model status is set appropriately when a solver's claimed optimality doesn't satify the general HiGHS primal/dual feasibilily tolerances. Affects IPM without crossover and PDLP

Command line parsing now done with pure C++11 code

Added command line flags to read basis from and write basis to a file 

Bug fixes in records of primal/dual rays

MPS read utility improved. Error logging is now less verbose; inability to handle USERCUTS section is properly logged

Implemented lifting for probing as described by Achterberg et al in _Presolve Reductions in Mixed Integer Programming._ INFORMS Journal on Computing 32(2):473-506 (2019). Not used by default, but option mip_lifting_for_probing allows it to be used with two levels of modification

Propagated updates from cuPDLP-C

Added GPU support for cuPDLP-C