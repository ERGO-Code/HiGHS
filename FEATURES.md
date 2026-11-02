The highlights of v1.16 are the use of parallelism by default in the
MIP solver, a rewritten LP file reader, a Hessian oracle for QP,
support for AMD GPUs in HiPDLP, Windows ARM64 Python wheels and NuGet
packages, and a large number of correctness fixes in the MIP solver
and presolve.

## Code changes

### Model file I/O

The LP file reader has been rewritten as a parser for a formal grammar
of the LP file format ([PR
#3355](https://github.com/ERGO-Code/HiGHS/pull/3355)). Invalid files
are now rejected, rather than being read silently as a different
model, and syntax errors are reported with the line and column of the
error. The reader is also faster and uses less memory ([PR
#3446](https://github.com/ERGO-Code/HiGHS/pull/3446)).

The floating point precision when writing to an MPS file was increased
([PR #3340](https://github.com/ERGO-Code/HiGHS/pull/3340)).

### Presolve

Following [PR #3167](https://github.com/ERGO-Code/HiGHS/pull/3167), a
"light" presolve is available via the new `presolve_light` option,
performing an initial sweep of low-cost presolve rules without forming
the full presolve data structures. Fourier-Motzkin elimination was added
to presolve in [PR #3255](https://github.com/ERGO-Code/HiGHS/pull/3255),
and new MIP techniques include a single machine scheduling separator
([PR #2990](https://github.com/ERGO-Code/HiGHS/pull/2990)), dual fixing
within probing ([PR #3244](https://github.com/ERGO-Code/HiGHS/pull/3244))
and aggregated implications
([PR #3264](https://github.com/ERGO-Code/HiGHS/pull/3264)). The options
`lp_presolve_requires_basis_postsolve` and
`use_implied_bounds_from_presolve` are now redundant
([PR #3334](https://github.com/ERGO-Code/HiGHS/pull/3334) and
[PR #3311](https://github.com/ERGO-Code/HiGHS/pull/3311)).

### Parallelism

Following [PR #3406](https://github.com/ERGO-Code/HiGHS/pull/3406), the
parallel MIP solver introduced in v1.15 is now enabled by default: the
branch-and-bound tree is explored using multiple threads unless the
`parallel` option is set to "off". The parallel design was reworked in
[PR #3192](https://github.com/ERGO-Code/HiGHS/pull/3192), with
deterministic early termination added in
[PR #3124](https://github.com/ERGO-Code/HiGHS/pull/3124).

Following [PR #3137](https://github.com/ERGO-Code/HiGHS/pull/3137),
HiGHS now uses one thread per physical core that the program is
allowed to run on, respecting the OS affinity mask, which fixes the
use of only half of the cores on Apple silicon. The parallel scheduler
documentation was updated in [PR
#3449](https://github.com/ERGO-Code/HiGHS/pull/3449).

### QP solvers

Rather than provide an explicit Hessian matrix, users can provide an
"oracle" callback that provides the necessary data associated with the
Hessian ([PR #3153](https://github.com/ERGO-Code/HiGHS/pull/3153)).
Callbacks for the first feasible point and user interrupts were added
to the QP active set method in [PR
#3214](https://github.com/ERGO-Code/HiGHS/pull/3214).  The QP active
set solver no longer reports some convex QPs as non-convex ([PR
#3325](https://github.com/ERGO-Code/HiGHS/pull/3325)).

### Simplex

A false unboundedness error in dual simplex was fixed ([PR
#3353](https://github.com/ERGO-Code/HiGHS/pull/3353)). Following [PR
#3421](https://github.com/ERGO-Code/HiGHS/pull/3421),
`kUnboundedOrInfeasible` is returned when this status cannot be
resolved, regardless of the `allow_unbounded_or_infeasible` option
setting.

### HiPDLP

HiPDLP now also handles maximization ([PR
#3469](https://github.com/ERGO-Code/HiGHS/pull/3469)).

### HiPO

HiPO now performs parallel solves
([PR #3200](https://github.com/ERGO-Code/HiGHS/pull/3200)), constructs
the normal equations in parallel
([PR #3187](https://github.com/ERGO-Code/HiGHS/pull/3187)) and uses a
simplified factorisation
([PR #3159](https://github.com/ERGO-Code/HiGHS/pull/3159)). The new
options `hipo_parallel_force`, `hipo_parallel_forbid` and
`hipo_serial_blas` control its use of parallelism, and a C interface to
the HiPO linear solver was added in
[PR #3066](https://github.com/ERGO-Code/HiGHS/pull/3066).

### IPX

IPX now
terminates when its iterates stall
([PR #3401](https://github.com/ERGO-Code/HiGHS/pull/3401)) and checks
the time limit when forming its initial basis
([PR #3185](https://github.com/ERGO-Code/HiGHS/pull/3185)).

### Crossover

Following [PR #3466](https://github.com/ERGO-Code/HiGHS/pull/3466),
`Highs::crossover` can be run from an approximate primal-dual solution,
such as one from PDLP, IPX or HiPO, including for maximization and LPs
with ranged rows.

### Interfaces

Following [PR #3257](https://github.com/ERGO-Code/HiGHS/pull/3257),
`highspy` has new `HighsError` and `HighsStatusError` exceptions and
supports minimize/maximize with a constant objective, and
`PYBIND11_DETAILED_ERROR_MESSAGES` is no longer defined in the bindings
([PR #3346](https://github.com/ERGO-Code/HiGHS/pull/3346)).

The C API methods `Highs_getBasis` and `Highs_getSolution` now return
`kHighsStatusError` if there is no valid basis or primal solution ([PR
#3291](https://github.com/ERGO-Code/HiGHS/pull/3291)).

### Other

Performance issues in `changeColsIntegrality`, `addRows` and when
names are supplied through the API were fixed in [PR
#3477](https://github.com/ERGO-Code/HiGHS/pull/3477), [PR
#3169](https://github.com/ERGO-Code/HiGHS/pull/3169) and [PR
#3335](https://github.com/ERGO-Code/HiGHS/pull/3335), and
`Highs::feasibilityRelaxation` now minimizes the violation for
maximization models ([PR
#3370](https://github.com/ERGO-Code/HiGHS/pull/3370)).

### HiGHS solution file

The HiGHS solution file can now hold a partial solution
([PR #3321](https://github.com/ERGO-Code/HiGHS/pull/3321)), the
objective bound and MIP relative gap are written to the pretty and raw
solution files ([PR #3431](https://github.com/ERGO-Code/HiGHS/pull/3431)),
and 

### Bug fixes

Many bugs that could lead to a suboptimal solution being reported as
optimal, or a feasible model being reported as infeasible, were fixed.
In the MIP solver these include invalid cuts
([PR #3177](https://github.com/ERGO-Code/HiGHS/pull/3177),
[PR #3179](https://github.com/ERGO-Code/HiGHS/pull/3179)), conflict
analysis ([PR #3275](https://github.com/ERGO-Code/HiGHS/pull/3275),
[PR #3443](https://github.com/ERGO-Code/HiGHS/pull/3443)) and fixings
from the analytic centre
([PR #3456](https://github.com/ERGO-Code/HiGHS/pull/3456)). In presolve
they include fractional bounds on integer variables
([PR #3272](https://github.com/ERGO-Code/HiGHS/pull/3272),
[PR #3276](https://github.com/ERGO-Code/HiGHS/pull/3276)),
`HPresolve::strengthenInequalities`
([PR #3234](https://github.com/ERGO-Code/HiGHS/pull/3234),
[PR #3468](https://github.com/ERGO-Code/HiGHS/pull/3468)), dual fixing
([PR #3411](https://github.com/ERGO-Code/HiGHS/pull/3411)), the clique
table ([PR #3220](https://github.com/ERGO-Code/HiGHS/pull/3220)) and a
number of postsolve fixes. Crashes in conflict analysis, feasibility
jump and symmetry detection were fixed in
[PR #3196](https://github.com/ERGO-Code/HiGHS/pull/3196),
[PR #3269](https://github.com/ERGO-Code/HiGHS/pull/3269) and
[PR #3182](https://github.com/ERGO-Code/HiGHS/pull/3182), and MIPs with
very large column bounds are handled more robustly
([PR #3267](https://github.com/ERGO-Code/HiGHS/pull/3267)). The
performance of feasibility jump was improved by fixing two bugs: jump
candidates with negative constraint coefficients are no longer skipped,
and the objective contribution now has the correct sign
([PR #3310](https://github.com/ERGO-Code/HiGHS/pull/3310)).

## Build changes

Following [PR #3423](https://github.com/ERGO-Code/HiGHS/pull/3423),
Windows ARM64 Python wheels and NuGet packages are now available, with
OpenBLAS built with optimised kernels. Free-threaded Python (cp314t)
wheels are now built
([PR #3158](https://github.com/ERGO-Code/HiGHS/pull/3158)), and
`highspy` can be built from the sdist with HiPO off
([PR #3480](https://github.com/ERGO-Code/HiGHS/pull/3480)).

Following [PR #3083](https://github.com/ERGO-Code/HiGHS/pull/3083),
HiPDLP can now run on AMD GPUs using ROCm/HIP. This requires a ROCm
installation providing the HIP compiler and the hipBLAS and hipSPARSE
libraries, and is enabled with the CMake option `-DHIPDLP_HIP=ON`;
ROCm 10.0 or newer is recommended. A missing GPU is now handled
gracefully ([PR #3388](https://github.com/ERGO-Code/HiGHS/pull/3388)).
See the
[installation documentation](https://ergo-code.github.io/HiGHS/dev/installation/#gpu-build-amd)
for details.

OpenBLAS is built with `-mfma` for its AVX2 `DYNAMIC_ARCH` kernels on
Linux ([PR #3399](https://github.com/ERGO-Code/HiGHS/pull/3399)), and on
Apple users can choose the BLAS library used for the HiPO build
([PR #3242](https://github.com/ERGO-Code/HiGHS/pull/3242)). CI now runs
more sanitizer tests on PRs
([PR #3474](https://github.com/ERGO-Code/HiGHS/pull/3474)).