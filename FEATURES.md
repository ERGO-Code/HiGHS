The highlights of v1.16 are that the parallel MIP solver is now used by
default, a rewritten LP file reader, a Hessian oracle for QP, support
for AMD GPUs in HiPDLP, Windows ARM64 Python wheels and NuGet packages,
and a large number of correctness fixes in the MIP solver and presolve.

## Code changes

### MIP solver

Following [PR #3406](https://github.com/ERGO-Code/HiGHS/pull/3406), the
parallel MIP solver introduced in v1.15 is now enabled by default: the
branch-and-bound tree is explored using multiple threads unless the
`parallel` option is set to "off". The worker/node-search design was
reworked in [PR #3192](https://github.com/ERGO-Code/HiGHS/pull/3192),
deterministic early termination when a worker proves that the search
can stop was added in
[PR #3124](https://github.com/ERGO-Code/HiGHS/pull/3124), the timeout
is now enforced while creating sub-MIPs
([PR #3237](https://github.com/ERGO-Code/HiGHS/pull/3237)), and the
parallel status line in the MIP log was updated in
[PR #3416](https://github.com/ERGO-Code/HiGHS/pull/3416).

Following [PR #3137](https://github.com/ERGO-Code/HiGHS/pull/3137),
HiGHS now uses one thread per physical core that the program is allowed
to run on, respecting the OS affinity mask, rather than halving the
logical core count. This fixes the use of only half of the cores on
Apple silicon, reported in
[#3221](https://github.com/ERGO-Code/HiGHS/issues/3221), where the
hard-coded cache-line size was also replaced by a constant
([PR #3281](https://github.com/ERGO-Code/HiGHS/pull/3281)). The
parallel scheduler documentation, particularly for the C API, was
updated in [PR #3148](https://github.com/ERGO-Code/HiGHS/pull/3148) and
[PR #3449](https://github.com/ERGO-Code/HiGHS/pull/3449), and the
`resetGlobalScheduler` calls were removed from the C API in
[PR #3141](https://github.com/ERGO-Code/HiGHS/pull/3141).

New MIP techniques include a single machine scheduling separator
([PR #2990](https://github.com/ERGO-Code/HiGHS/pull/2990)), dual fixing
embedded into probing
([PR #3244](https://github.com/ERGO-Code/HiGHS/pull/3244)), aggregated
implications ([PR #3264](https://github.com/ERGO-Code/HiGHS/pull/3264)),
implications kept between restarts
([PR #3209](https://github.com/ERGO-Code/HiGHS/pull/3209)), enhanced
clique extraction
([PR #3286](https://github.com/ERGO-Code/HiGHS/pull/3286)), and the
clique table maintained through presolve
([PR #3312](https://github.com/ERGO-Code/HiGHS/pull/3312)). The strong
branching budget is reduced before two restarts
([PR #3213](https://github.com/ERGO-Code/HiGHS/pull/3213)), dive length
is limited ([PR #3336](https://github.com/ERGO-Code/HiGHS/pull/3336)),
and the number of nodes in the RENS/RINS depth-first search is limited
([PR #3445](https://github.com/ERGO-Code/HiGHS/pull/3445)).

Several bugs that led to a suboptimal solution being reported as
optimal have been fixed. Prompted by
[#3170](https://github.com/ERGO-Code/HiGHS/issues/3170) and
[#3171](https://github.com/ERGO-Code/HiGHS/issues/3171), invalid cuts
from redundant variable bounds and from inconsistent substitutions in
the path-mixing transform are no longer generated
([PR #3177](https://github.com/ERGO-Code/HiGHS/pull/3177) and
[PR #3179](https://github.com/ERGO-Code/HiGHS/pull/3179)). Conflict
analysis no longer derives empty conflicts from objective-based bound
changes ([#3273](https://github.com/ERGO-Code/HiGHS/issues/3273), fixed
in [PR #3275](https://github.com/ERGO-Code/HiGHS/pull/3275)), and no
longer stores a weaker bound as the strongest
([#3433](https://github.com/ERGO-Code/HiGHS/issues/3433), fixed in
[PR #3443](https://github.com/ERGO-Code/HiGHS/pull/3443)). Prompted by
[#3451](https://github.com/ERGO-Code/HiGHS/issues/3451) and
[#3452](https://github.com/ERGO-Code/HiGHS/issues/3452), variables are
now fixed using the analytic centre only when it is well centred
([PR #3456](https://github.com/ERGO-Code/HiGHS/pull/3456) and
[PR #3408](https://github.com/ERGO-Code/HiGHS/pull/3408)).

Prompted by [#3190](https://github.com/ERGO-Code/HiGHS/issues/3190),
crashes and hangs in conflict analysis under extreme numerical ranges
were fixed ([PR #3196](https://github.com/ERGO-Code/HiGHS/pull/3196)),
and a crash in `HighsConflictPool::syncConflictPool` was fixed in
[PR #3439](https://github.com/ERGO-Code/HiGHS/pull/3439). The
feasibility jump heuristic was made more robust, fixing a segfault
reported in [#3262](https://github.com/ERGO-Code/HiGHS/issues/3262)
([PR #3269](https://github.com/ERGO-Code/HiGHS/pull/3269)), and the
jump value is now computed correctly for rows with negative
coefficients ([#3302](https://github.com/ERGO-Code/HiGHS/issues/3302)
and [#3303](https://github.com/ERGO-Code/HiGHS/issues/3303), fixed in
[PR #3310](https://github.com/ERGO-Code/HiGHS/pull/3310)). An
out-of-bounds access in symmetry detection
([#3166](https://github.com/ERGO-Code/HiGHS/issues/3166)) was fixed in
[PR #3182](https://github.com/ERGO-Code/HiGHS/pull/3182), and symmetry
handling now always retains the branch as a reason
([PR #3357](https://github.com/ERGO-Code/HiGHS/pull/3357)).

Prompted by [#2866](https://github.com/ERGO-Code/HiGHS/issues/2866),
MIPs with very large column bounds are now handled more robustly
([PR #3139](https://github.com/ERGO-Code/HiGHS/pull/3139),
[PR #3184](https://github.com/ERGO-Code/HiGHS/pull/3184) and
[PR #3267](https://github.com/ERGO-Code/HiGHS/pull/3267)), and a warning
is given for badly scaled models with large bounds on integer variables
([PR #3243](https://github.com/ERGO-Code/HiGHS/pull/3243)). A MIP
solution is no longer discarded when its repair LP fails to solve
([PR #3181](https://github.com/ERGO-Code/HiGHS/pull/3181), fixing
[#3180](https://github.com/ERGO-Code/HiGHS/issues/3180)), and stalling
in `HighsRedcostFixing::addRootRedcost` was avoided in
[PR #3241](https://github.com/ERGO-Code/HiGHS/pull/3241).

Prompted by [#3360](https://github.com/ERGO-Code/HiGHS/issues/3360), an
unnecessary copy of the LP relaxation when setting up the MIP solver was
removed ([PR #3363](https://github.com/ERGO-Code/HiGHS/pull/3363)),
recovering part of a performance regression on small MIPs. Prompted by
[#3316](https://github.com/ERGO-Code/HiGHS/issues/3316), MIP solves are
no longer slower when row and column names are supplied through the
API ([PR #3335](https://github.com/ERGO-Code/HiGHS/pull/3335)).

Prompted by [#3414](https://github.com/ERGO-Code/HiGHS/issues/3414),
MIP feasibility is now checked with the MIP tolerances after removing
user scaling, via the new methods
`HighsOptions::getMipPrimalFeasibilityTolerance()` and
`HighsOptions::getMipDualFeasibilityTolerance()`
([PR #3432](https://github.com/ERGO-Code/HiGHS/pull/3432)). Prompted by
[#3178](https://github.com/ERGO-Code/HiGHS/issues/3178), the LP offset
and any dual bound in `HighsInfo` are now scaled consistently with user
cost and bound scaling
([PR #3186](https://github.com/ERGO-Code/HiGHS/pull/3186)).

### Presolve

Following [PR #3167](https://github.com/ERGO-Code/HiGHS/pull/3167), a
"light" presolve is available via the new `presolve_light` option,
performing an initial sweep of low-cost presolve rules without forming
the full presolve data structures
([#3040](https://github.com/ERGO-Code/HiGHS/issues/3040)). The initial
sweep was refactored in
[PR #3210](https://github.com/ERGO-Code/HiGHS/pull/3210) and
[PR #3300](https://github.com/ERGO-Code/HiGHS/pull/3300), and a crash
in it was fixed in
[PR #3442](https://github.com/ERGO-Code/HiGHS/pull/3442).

Following [PR #3255](https://github.com/ERGO-Code/HiGHS/pull/3255),
Fourier-Motzkin elimination is available in presolve, building on
infrastructure for adding rows and columns in presolve
([PR #3218](https://github.com/ERGO-Code/HiGHS/pull/3218) and
[PR #3250](https://github.com/ERGO-Code/HiGHS/pull/3250)). Continuous
singletons are now removed from double-sided rows
([PR #2962](https://github.com/ERGO-Code/HiGHS/pull/2962)), zero-cost
singletons are handled for implied integers
([PR #3455](https://github.com/ERGO-Code/HiGHS/pull/3455)),
substitution candidates are detected when a ranged row becomes dual
implied free ([PR #3437](https://github.com/ERGO-Code/HiGHS/pull/3437)),
and dual substitution has an early stop and a row size limit
([PR #3404](https://github.com/ERGO-Code/HiGHS/pull/3404)). The
`removeSlacks` reduction, a duplicate of the zero-cost singleton
reduction, was removed in
[PR #3309](https://github.com/ERGO-Code/HiGHS/pull/3309).

Following [PR #3334](https://github.com/ERGO-Code/HiGHS/pull/3334),
presolve tracks whether basis or primal-dual postsolve may be required,
so that reductions that cannot be postsolved accordingly are not used
([PR #3329](https://github.com/ERGO-Code/HiGHS/pull/3329),
[PR #3342](https://github.com/ERGO-Code/HiGHS/pull/3342) and
[PR #3348](https://github.com/ERGO-Code/HiGHS/pull/3348)). As a result,
the options `lp_presolve_requires_basis_postsolve` and
`use_implied_bounds_from_presolve` are now redundant
([PR #3311](https://github.com/ERGO-Code/HiGHS/pull/3311)).

Several presolve bugs leading to incorrect results were fixed.
Fractional bounds on integer variables are now rounded
([#3271](https://github.com/ERGO-Code/HiGHS/issues/3271) and
[#3274](https://github.com/ERGO-Code/HiGHS/issues/3274), fixed in
[PR #3272](https://github.com/ERGO-Code/HiGHS/pull/3272) and
[PR #3276](https://github.com/ERGO-Code/HiGHS/pull/3276)), and
integrality is no longer lost when upgrading bounds from variable
bounds ([PR #3268](https://github.com/ERGO-Code/HiGHS/pull/3268) and
[PR #3320](https://github.com/ERGO-Code/HiGHS/pull/3320)). Small
contributions are now handled correctly in
`HPresolve::strengthenInequalities`
([#2900](https://github.com/ERGO-Code/HiGHS/issues/2900) and
[#3461](https://github.com/ERGO-Code/HiGHS/issues/3461), fixed in
[PR #3234](https://github.com/ERGO-Code/HiGHS/pull/3234) and
[PR #3468](https://github.com/ERGO-Code/HiGHS/pull/3468)), and
worst-case bounds are rounded correctly in dual fixing
([#3389](https://github.com/ERGO-Code/HiGHS/issues/3389), fixed in
[PR #3411](https://github.com/ERGO-Code/HiGHS/pull/3411)). Prompted by
[#3216](https://github.com/ERGO-Code/HiGHS/issues/3216), clique table
fixes in [PR #3220](https://github.com/ERGO-Code/HiGHS/pull/3220) allow
feasible incumbents to be reconstructed after presolve. Prompted by
[#3391](https://github.com/ERGO-Code/HiGHS/issues/3391), dual
infeasibility is now checked in `HPresolve::rowPresolve`
([PR #3392](https://github.com/ERGO-Code/HiGHS/pull/3392)), and
[PR #3367](https://github.com/ERGO-Code/HiGHS/pull/3367) fixes stale row
positions in parallel-row detection.

Prompted by [#3140](https://github.com/ERGO-Code/HiGHS/issues/3140), a
regression in dual fixing for single equations was fixed
([PR #3232](https://github.com/ERGO-Code/HiGHS/pull/3232)), and a
further presolve regression was fixed in
[PR #3318](https://github.com/ERGO-Code/HiGHS/pull/3318). Dual-implied
bounds are now reset when column bounds change or equalities are
upgraded ([PR #3270](https://github.com/ERGO-Code/HiGHS/pull/3270) and
[PR #3332](https://github.com/ERGO-Code/HiGHS/pull/3332)), and redundant
row flags are no longer reset
([PR #3341](https://github.com/ERGO-Code/HiGHS/pull/3341)). Prompted by
[#3118](https://github.com/ERGO-Code/HiGHS/issues/3118), the solver no
longer reports an infeasible problem as optimal when a starting
solution is supplied
([PR #3135](https://github.com/ERGO-Code/HiGHS/pull/3135)), and edge-case
assertion failures in presolve and domain propagation were fixed in
[PR #3149](https://github.com/ERGO-Code/HiGHS/pull/3149).

Postsolve fixes include primal-dual postsolve for sparsify
([PR #3324](https://github.com/ERGO-Code/HiGHS/pull/3324)), parallel
row postsolve in the basis case
([PR #3313](https://github.com/ERGO-Code/HiGHS/pull/3313)), dominated
and duplicate columns
([PR #3285](https://github.com/ERGO-Code/HiGHS/pull/3285) and
[PR #3290](https://github.com/ERGO-Code/HiGHS/pull/3290)), parallel
column merging ([PR #3261](https://github.com/ERGO-Code/HiGHS/pull/3261)),
and the dual flip after postsolve
([PR #3248](https://github.com/ERGO-Code/HiGHS/pull/3248)). Prompted by
[#3323](https://github.com/ERGO-Code/HiGHS/issues/3323), the basis is
now checked for consistency after postsolve
([PR #3350](https://github.com/ERGO-Code/HiGHS/pull/3350)). Following
[PR #3299](https://github.com/ERGO-Code/HiGHS/pull/3299), the active
costs for an LP solution are identified using effective costs after
free column singletons have been substituted. Presolve/postsolve debug
utilities were added in
[PR #3219](https://github.com/ERGO-Code/HiGHS/pull/3219), and
`HPresolve.cpp` was split up in
[PR #3398](https://github.com/ERGO-Code/HiGHS/pull/3398).

### LP and QP solvers

Rather than provide an explicit Hessian matrix, users can provide an
"oracle" callback that provides the necessary data associated with the
Hessian ([PR #3153](https://github.com/ERGO-Code/HiGHS/pull/3153)).
Callbacks for the first feasible point and user interrupts were added
to the QP active set method in
[PR #3214](https://github.com/ERGO-Code/HiGHS/pull/3214).

Prompted by [#3322](https://github.com/ERGO-Code/HiGHS/issues/3322),
the QP active set solver no longer reports convex QPs as non-convex
depending on column order
([PR #3325](https://github.com/ERGO-Code/HiGHS/pull/3325)), and NaNs and
a false "Unbounded" status caused by zero-sized Cholesky storage
([#3453](https://github.com/ERGO-Code/HiGHS/issues/3453)) were fixed in
[PR #3458](https://github.com/ERGO-Code/HiGHS/pull/3458).

Following [PR #3466](https://github.com/ERGO-Code/HiGHS/pull/3466),
`Highs::crossover` has been generalised so that it can be run from an
approximate primal-dual solution, such as one from PDLP, IPX or HiPO,
including for maximization and LPs with ranged rows
([#3459](https://github.com/ERGO-Code/HiGHS/issues/3459)).

Prompted by [#3397](https://github.com/ERGO-Code/HiGHS/issues/3397),
IPX now terminates when its iterates stall
([PR #3401](https://github.com/ERGO-Code/HiGHS/pull/3401)), a timeout
was added to forming the initial basis and factorization in IPX
([#3120](https://github.com/ERGO-Code/HiGHS/issues/3120), fixed in
[PR #3185](https://github.com/ERGO-Code/HiGHS/pull/3185)), empty rows
are accounted for in the IPX solution
([PR #3319](https://github.com/ERGO-Code/HiGHS/pull/3319)), and the
primal-dual gap computed by IPX and HiPO now accounts for the objective
offset ([PR #3427](https://github.com/ERGO-Code/HiGHS/pull/3427)).
Prompted by [#3413](https://github.com/ERGO-Code/HiGHS/issues/3413),
IPM is no longer used to try to get a basic solution when computing rays
([PR #3435](https://github.com/ERGO-Code/HiGHS/pull/3435)).

Prompted by [#2300](https://github.com/ERGO-Code/HiGHS/issues/2300), a
false unboundedness error in dual simplex phase 1 was fixed
([PR #3353](https://github.com/ERGO-Code/HiGHS/pull/3353)), excessive
primal value rebuilds in simplex are now handled
([PR #3351](https://github.com/ERGO-Code/HiGHS/pull/3351)), and an error
when the objective is zero and `output_flag` is false was fixed in
[PR #3226](https://github.com/ERGO-Code/HiGHS/pull/3226). Prompted by
[#3383](https://github.com/ERGO-Code/HiGHS/issues/3383),
`kUnboundedOrInfeasible` is returned when this status cannot be
resolved, regardless of the `allow_unbounded_or_infeasible` option
setting ([PR #3421](https://github.com/ERGO-Code/HiGHS/pull/3421)).

Prompted by [#2924](https://github.com/ERGO-Code/HiGHS/issues/2924) and
[#3150](https://github.com/ERGO-Code/HiGHS/issues/3150), the calls to
simplex, the QP solver and the MIP solver are now wrapped in try-catch
blocks ([PR #3168](https://github.com/ERGO-Code/HiGHS/pull/3168)).

### HiPO and HiPDLP

HiPO now performs parallel solves
([PR #3200](https://github.com/ERGO-Code/HiGHS/pull/3200)), constructs
the normal equations in parallel
([PR #3187](https://github.com/ERGO-Code/HiGHS/pull/3187)) and uses a
simplified factorisation
([PR #3159](https://github.com/ERGO-Code/HiGHS/pull/3159)). The new
options `hipo_parallel_force` and `hipo_parallel_forbid` control the
parallel techniques used, and `hipo_serial_blas` forces OpenBLAS to run
in serial when used by HiPO. A C interface to the HiPO linear solver
was added in [PR #3066](https://github.com/ERGO-Code/HiGHS/pull/3066).

Prompted by [#3236](https://github.com/ERGO-Code/HiGHS/issues/3236),
the HiPO status is no longer discarded if IPX refinement fails
([PR #3284](https://github.com/ERGO-Code/HiGHS/pull/3284)), and prompted
by [#3369](https://github.com/ERGO-Code/HiGHS/issues/3369), HiPO no
longer cycles because of a zero slack for removed empty rows
([PR #3376](https://github.com/ERGO-Code/HiGHS/pull/3376)). Further
small HiPO bugs were fixed in
[PR #3448](https://github.com/ERGO-Code/HiGHS/pull/3448).

Prompted by
[#3467](https://github.com/ERGO-Code/HiGHS/issues/3467), HiPDLP now
handles maximization and removes free rows rather than adding free
slacks ([PR #3469](https://github.com/ERGO-Code/HiGHS/pull/3469)). The
time limit for cuPDLP-C was fixed in
[PR #3450](https://github.com/ERGO-Code/HiGHS/pull/3450).

### File I/O

The LP file reader has been rewritten as a parser for a formal grammar
of the LP file format
([PR #3355](https://github.com/ERGO-Code/HiGHS/pull/3355)). Invalid
files are now rejected, rather than being read silently as a different
model ([PR #3343](https://github.com/ERGO-Code/HiGHS/pull/3343)), and
syntax errors are reported with the line and column of the error. The
reader is also faster and uses less memory
([PR #3446](https://github.com/ERGO-Code/HiGHS/pull/3446)). A warning
is now given for constant LHS entries in LP file constraints
([PR #3283](https://github.com/ERGO-Code/HiGHS/pull/3283)), and
semi-integer variables are now written correctly to LP files
([PR #3380](https://github.com/ERGO-Code/HiGHS/pull/3380)).

Floating point precision when writing to file was increased in
[PR #3340](https://github.com/ERGO-Code/HiGHS/pull/3340), and semi
fixed bounds in MPS files are now read correctly
([PR #3390](https://github.com/ERGO-Code/HiGHS/pull/3390)). Prompted by
[#3317](https://github.com/ERGO-Code/HiGHS/issues/3317), the HiGHS
solution file definition has been extended to accommodate a partial
solution ([PR #3321](https://github.com/ERGO-Code/HiGHS/pull/3321)),
and prompted by [#3417](https://github.com/ERGO-Code/HiGHS/issues/3417),
the objective bound and MIP relative gap are now written to the pretty
and raw solution files
([PR #3431](https://github.com/ERGO-Code/HiGHS/pull/3431)).

### Interfaces and API

Following [PR #3257](https://github.com/ERGO-Code/HiGHS/pull/3257),
`highspy` has new `HighsError` and `HighsStatusError` exceptions,
supports minimize/maximize with a constant objective, and no longer
makes repeated calls to the same callback. Prompted by
[#3144](https://github.com/ERGO-Code/HiGHS/issues/3144), the definition
of `PYBIND11_DETAILED_ERROR_MESSAGES` was removed from the `highspy`
bindings, as it broke the ODR for other code using pybind
([PR #3346](https://github.com/ERGO-Code/HiGHS/pull/3346)).

Prompted by [#3288](https://github.com/ERGO-Code/HiGHS/issues/3288), the
C API methods `Highs_getBasis` and `Highs_getSolution` now return
`kHighsStatusError` if there is no valid basis or primal solution
([PR #3291](https://github.com/ERGO-Code/HiGHS/pull/3291)), and the enum
constants in `highs_c_api.h` were updated in
[PR #3403](https://github.com/ERGO-Code/HiGHS/pull/3403).

Prompted by [#3476](https://github.com/ERGO-Code/HiGHS/issues/3476),
`n` calls to `changeColsIntegrality` on an LP now have linear rather
than quadratic cost ([PR #3477](https://github.com/ERGO-Code/HiGHS/pull/3477)).
Prompted by [#3129](https://github.com/ERGO-Code/HiGHS/issues/3129),
`addRows` is fast again for matrices with many explicit zeros
([PR #3169](https://github.com/ERGO-Code/HiGHS/pull/3169)), and
extra copies of `HighsSparseMatrix` are avoided when building the
opposite format ([PR #3428](https://github.com/ERGO-Code/HiGHS/pull/3428)).
Following [PR #3386](https://github.com/ERGO-Code/HiGHS/pull/3386), the
incumbent model is not modified when the data supplied yields an error.

Prompted by [#3368](https://github.com/ERGO-Code/HiGHS/issues/3368),
`Highs::feasibilityRelaxation` now minimizes the violation when the
model is a maximization
([PR #3370](https://github.com/ERGO-Code/HiGHS/pull/3370)), and the
`iis_strategy` option description was corrected in
[PR #3436](https://github.com/ERGO-Code/HiGHS/pull/3436).

Prompted by [#3164](https://github.com/ERGO-Code/HiGHS/issues/3164), row
and column names are now given when logging illegal data in
`assessBounds` and `assessMatrix`
([PR #3165](https://github.com/ERGO-Code/HiGHS/pull/3165)). User logging
for the dependent equations search and parallel simplex was improved
in [PR #3238](https://github.com/ERGO-Code/HiGHS/pull/3238), and the
message given when the value of the `solver` option is not valid was
fixed in [PR #3162](https://github.com/ERGO-Code/HiGHS/pull/3162).

Memory and thread safety issues detected by sanitizers were fixed in
[PR #3306](https://github.com/ERGO-Code/HiGHS/pull/3306),
[PR #3430](https://github.com/ERGO-Code/HiGHS/pull/3430) (prompted by
[#3419](https://github.com/ERGO-Code/HiGHS/issues/3419)),
[PR #3472](https://github.com/ERGO-Code/HiGHS/pull/3472),
[PR #3475](https://github.com/ERGO-Code/HiGHS/pull/3475) and
[PR #3481](https://github.com/ERGO-Code/HiGHS/pull/3481), and
`std::vector<bool>` was removed in
[PR #3197](https://github.com/ERGO-Code/HiGHS/pull/3197).

## Build changes

Following [PR #3423](https://github.com/ERGO-Code/HiGHS/pull/3423),
Windows ARM64 Python wheels and NuGet packages are now available, with
OpenBLAS built with optimised kernels. Free-threaded Python (cp314t)
wheels are now built
([PR #3158](https://github.com/ERGO-Code/HiGHS/pull/3158)), and
`highspy` can be built from the sdist with HiPO off
([PR #3480](https://github.com/ERGO-Code/HiGHS/pull/3480)).

OpenBLAS is built with `-mfma` for its AVX2 `DYNAMIC_ARCH` kernels on
Linux ([PR #3399](https://github.com/ERGO-Code/HiGHS/pull/3399)), and
with locking on and threading off on Windows ARM64
([PR #3482](https://github.com/ERGO-Code/HiGHS/pull/3482)). There is an
option to not set the number of OpenBLAS threads to 1
([PR #3402](https://github.com/ERGO-Code/HiGHS/pull/3402)), and on Apple
users can choose the BLAS library used for the HiPO build
([PR #3242](https://github.com/ERGO-Code/HiGHS/pull/3242)).

A CUDA build fix was merged in
[PR #3378](https://github.com/ERGO-Code/HiGHS/pull/3378), the build was
fixed for Xcode 26.6 in
[PR #3228](https://github.com/ERGO-Code/HiGHS/pull/3228), the precise
floating point model is used with the Intel compiler
([PR #3447](https://github.com/ERGO-Code/HiGHS/pull/3447)),
`isDuplicate` was fixed for 32-bit
([PR #3349](https://github.com/ERGO-Code/HiGHS/pull/3349)), and Bazel
fixes were merged in
[PR #3147](https://github.com/ERGO-Code/HiGHS/pull/3147).

Following [PR #3082](https://github.com/ERGO-Code/HiGHS/pull/3083),
HiPDLP can now run on AMD GPUs using ROCm/HIP. This requires a ROCm
installation providing the HIP compiler and the hipBLAS and hipSPARSE
libraries, and is enabled with the CMake option `-DHIPDLP_HIP=ON`;
ROCm 9.0 or newer is recommended. A missing GPU is now handled
gracefully ([PR #3387](https://github.com/ERGO-Code/HiGHS/pull/3388)),
and the ROCm headers of the HIP compiler are now preferred over a
system HIP install
([PR #3485](https://github.com/ERGO-Code/HiGHS/pull/3486)). See the
[installation documentation](https://ergo-code.github.io/HiGHS/dev/installation/#gpu-build-amd)
for details.


CI now runs more sanitizer tests on PRs
([PR #3175](https://github.com/ERGO-Code/HiGHS/pull/3175) and
[PR #3474](https://github.com/ERGO-Code/HiGHS/pull/3474)), uses Visual
Studio 18 ([PR #3337](https://github.com/ERGO-Code/HiGHS/pull/3337)),
and the debug-mode tests were sped up in
[PR #3330](https://github.com/ERGO-Code/HiGHS/pull/3330). The Julia
`build_tarballs.jl` `BLA_VENDOR` setting was updated in
[PR #3161](https://github.com/ERGO-Code/HiGHS/pull/3161), and the
documentation was updated to Documenter@1 and Clang@0.19 in
[PR #3229](https://github.com/ERGO-Code/HiGHS/pull/3229).
