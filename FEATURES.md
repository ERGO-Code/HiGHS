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

Following [PR #3167](https://github.com/ERGO-Code/HiGHS/pull/3167),
presolve starts with an initial sweep of low-cost presolve rules
without forming the full presolve data structures. The new
`presolve_mode` option ([PR
#3536](https://github.com/ERGO-Code/HiGHS/pull/3536)) selects the
presolve that is performed: 0 (default) uses all rules, unless PDLP is
used, in which case the "light" rules are used; 1 uses only the
initial sweep; 2 uses the "light" rules; and 3 uses all rules. The
rules can be modified further using `presolve_rule_off`. Fourier-Motzkin elimination was added
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
[PR #3066](https://github.com/ERGO-Code/HiGHS/pull/3066). HiPO log
output is now shown when a log callback is used with `log_to_console`
set to false ([PR #3592](https://github.com/ERGO-Code/HiGHS/pull/3592)),
and IPX is always used to compute the analytic centre for the MIP
solver, since HiPO does not yet perform centring
([PR #3522](https://github.com/ERGO-Code/HiGHS/pull/3522)).

### IPX

IPX now
terminates when its iterates stall
([PR #3401](https://github.com/ERGO-Code/HiGHS/pull/3401)) and checks
the time limit when forming its initial basis
([PR #3185](https://github.com/ERGO-Code/HiGHS/pull/3185)). False
infeasibility and solve errors on small, badly scaled LPs were fixed in
[PR #3513](https://github.com/ERGO-Code/HiGHS/pull/3513).

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

Following [PR #3372](https://github.com/ERGO-Code/HiGHS/pull/3372),
`Highs::passModel` and `Highs::passHessian` trim any vectors that are
larger than the model dimensions require, and return
`HighsStatus::kWarning`. User-supplied solutions are now assessed in
the same way in all contexts ([PR
#3556](https://github.com/ERGO-Code/HiGHS/pull/3556)), and the HiGHS
log header is written by more methods, including `writeModel` and
`writeSolution` ([PR #3496](https://github.com/ERGO-Code/HiGHS/pull/3496)
and [PR #3517](https://github.com/ERGO-Code/HiGHS/pull/3517)).

**Breaking change:** The C API methods `Highs_getBasis` and
`Highs_getSolution` now return `kHighsStatusError`, rather than
`kHighsStatusOk`, if there is no valid basis or primal solution ([PR
#3291](https://github.com/ERGO-Code/HiGHS/pull/3291)). This includes
an LP found to be infeasible, for which `Highs_run` still returns
`kHighsStatusOk`. The output arrays are not written when an error is
returned, so check the return value before reading them. The soname of
the shared library remains `libhighs.so.1`, so existing binaries get
this behaviour without being rebuilt.

### Loading `highs_extras`

When HiPO is built into the separate shared library `highs_extras`,
HiGHS loads it at run time. Following [PR
#3547](https://github.com/ERGO-Code/HiGHS/pull/3547), the locations
from which it is loaded are restricted, so that HiGHS cannot be made to
load a different library placed in a directory that it searches.

- The new methods `Highs::initializeExtras(directory, search_if_empty)`
  and `Highs_initializeExtras` load the library from a given
  directory. If `directory` is empty, the default locations are
  searched only if `search_if_empty` is true. Call one of these
  methods once, before HiPO is first used. Later calls return the
  result of the first call. Neither method is needed when the extras
  are linked statically.
- If neither method is called, HiGHS searches the default locations
  when HiPO is first used. On Windows these are now only the folder of
  the application's `.exe` and the system folder: the `PATH` and the
  current folder are no longer searched. Any DLL that `highs_extras.dll`
  depends on, such as `openblas.dll`, must be in one of these folders.
- `highspy` loads the library only from the installed
  `highspy-extras` package.

If the library cannot be loaded, the reason is now reported
([PR #3555](https://github.com/ERGO-Code/HiGHS/pull/3555)).

### Other

Performance issues in `changeColsIntegrality`, `addRows` and when
names are supplied through the API were fixed in [PR
#3477](https://github.com/ERGO-Code/HiGHS/pull/3477), [PR
#3169](https://github.com/ERGO-Code/HiGHS/pull/3169) and [PR
#3335](https://github.com/ERGO-Code/HiGHS/pull/3335), and
`Highs::feasibilityRelaxation` now minimizes the violation for
maximization models ([PR
#3370](https://github.com/ERGO-Code/HiGHS/pull/3370)). Two checks
that ran on every simplex solve but whose results were not used, taking
about 2% of MIP solution time, are now only run in debug mode ([PR
#3515](https://github.com/ERGO-Code/HiGHS/pull/3515)), the enumeration
presolve makes far fewer row comparisons ([PR
#3529](https://github.com/ERGO-Code/HiGHS/pull/3529)), and changing
the orientation of a sparse matrix no longer copies it ([PR
#3485](https://github.com/ERGO-Code/HiGHS/pull/3485)).

The primal-dual integral is now bounded above by the solution time
([PR #3521](https://github.com/ERGO-Code/HiGHS/pull/3521)).

### HiGHS solution file

The HiGHS solution file can now hold a partial solution
([PR #3321](https://github.com/ERGO-Code/HiGHS/pull/3321)), the
objective bound and MIP relative gap are written to the pretty and raw
solution files ([PR #3431](https://github.com/ERGO-Code/HiGHS/pull/3431)),
and an infinite loop when reading an invalid solution file was fixed
([PR #3504](https://github.com/ERGO-Code/HiGHS/pull/3504)).

### Bug fixes

Many bugs that could lead to a suboptimal solution being reported as
optimal, or a feasible model being reported as infeasible, were fixed.
In the MIP solver these include invalid cuts
([PR #3177](https://github.com/ERGO-Code/HiGHS/pull/3177),
[PR #3179](https://github.com/ERGO-Code/HiGHS/pull/3179),
[PR #3507](https://github.com/ERGO-Code/HiGHS/pull/3507),
[PR #3591](https://github.com/ERGO-Code/HiGHS/pull/3591)), conflict
analysis ([PR #3275](https://github.com/ERGO-Code/HiGHS/pull/3275),
[PR #3443](https://github.com/ERGO-Code/HiGHS/pull/3443)), conflict
and clique propagation
([PR #3571](https://github.com/ERGO-Code/HiGHS/pull/3571),
[PR #3540](https://github.com/ERGO-Code/HiGHS/pull/3540)), domain
propagation with infinite column ranges
([PR #3525](https://github.com/ERGO-Code/HiGHS/pull/3525)) and fixings
from the analytic centre
([PR #3456](https://github.com/ERGO-Code/HiGHS/pull/3456)). In presolve
they include fractional bounds on integer variables
([PR #3272](https://github.com/ERGO-Code/HiGHS/pull/3272),
[PR #3276](https://github.com/ERGO-Code/HiGHS/pull/3276)),
`HPresolve::strengthenInequalities`
([PR #3234](https://github.com/ERGO-Code/HiGHS/pull/3234),
[PR #3468](https://github.com/ERGO-Code/HiGHS/pull/3468)), dual fixing
([PR #3411](https://github.com/ERGO-Code/HiGHS/pull/3411)), the clique
table ([PR #3220](https://github.com/ERGO-Code/HiGHS/pull/3220)),
variable bound aggregation
([PR #3562](https://github.com/ERGO-Code/HiGHS/pull/3562)), parallel
row merging ([PR #3520](https://github.com/ERGO-Code/HiGHS/pull/3520)),
loss of precision in comparisons with the feasibility tolerance
([PR #3554](https://github.com/ERGO-Code/HiGHS/pull/3554),
[PR #3575](https://github.com/ERGO-Code/HiGHS/pull/3575)), the
objective offset after the initial sweep
([PR #3548](https://github.com/ERGO-Code/HiGHS/pull/3548)) and a
number of postsolve fixes. The zero cost singleton column rule is now
switched off when dual postsolve may be required
([PR #3498](https://github.com/ERGO-Code/HiGHS/pull/3498)). Crashes in conflict analysis, feasibility
jump, symmetry detection, `removeDoubletonEquations` and the cut pool
when `mip_pool_age_limit` is zero were fixed in
[PR #3196](https://github.com/ERGO-Code/HiGHS/pull/3196),
[PR #3269](https://github.com/ERGO-Code/HiGHS/pull/3269),
[PR #3182](https://github.com/ERGO-Code/HiGHS/pull/3182),
[PR #3566](https://github.com/ERGO-Code/HiGHS/pull/3566) and
[PR #3574](https://github.com/ERGO-Code/HiGHS/pull/3574). Symmetry
detection now respects the time limit and no longer stalls on hash
collisions ([PR #3535](https://github.com/ERGO-Code/HiGHS/pull/3535)),
and MIPs with
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
([PR #3158](https://github.com/ERGO-Code/HiGHS/pull/3158),
[PR #3502](https://github.com/ERGO-Code/HiGHS/pull/3502)), as are
wheels for Python 3.15
([PR #3542](https://github.com/ERGO-Code/HiGHS/pull/3542)), and
`highspy` can be built from the sdist with HiPO off
([PR #3480](https://github.com/ERGO-Code/HiGHS/pull/3480)).

Following [PR #3083](https://github.com/ERGO-Code/HiGHS/pull/3083),
HiPDLP can now run on AMD GPUs using ROCm/HIP. This requires a ROCm
installation providing the HIP compiler and the hipBLAS and hipSPARSE
libraries, and is enabled with the CMake option `-DHIPDLP_HIP=ON`;
ROCm 10.0 or newer is recommended. A missing GPU is now handled
gracefully ([PR #3388](https://github.com/ERGO-Code/HiGHS/pull/3388)),
and the ROCm headers of the HIP compiler are preferred over a system
HIP installation
([PR #3486](https://github.com/ERGO-Code/HiGHS/pull/3486)).
See the
[installation documentation](https://ergo-code.github.io/HiGHS/dev/installation/#gpu-build-amd)
for details.

Following [PR #3409](https://github.com/ERGO-Code/HiGHS/pull/3409),
HiPO is built with clang-cl on Windows x64, as it already was on ARM64,
for the release binaries and Python wheels. This makes HiPO
substantially faster on Windows, matching its performance on Linux and
MacOS. vcpkg triplets for clang-cl are provided in
`cmake/vcpkg-triplets`
([PR #3530](https://github.com/ERGO-Code/HiGHS/pull/3530),
[PR #3541](https://github.com/ERGO-Code/HiGHS/pull/3541)).

`highs.pc` now lists the private dependencies of a static build, so
`pkg-config --static` gives C++ applications the zlib and `highs_extras`
libraries needed to link it, and OpenBLAS or Apple Accelerate if
`highs_extras` uses one of them
([PR #3572](https://github.com/ERGO-Code/HiGHS/pull/3572)).
`find_package(highs)` now works for a HiPO build with every BLAS route,
including bundled and static vcpkg OpenBLAS
([PR #3567](https://github.com/ERGO-Code/HiGHS/pull/3567)). An invalid
OpenBLAS CMake file, such as the one in Ubuntu 26.04, is now ignored in
favour of the standard BLAS search
([PR #3489](https://github.com/ERGO-Code/HiGHS/pull/3489)). The
LGPL-licensed sorting routines in the bundled GKlib were replaced by
`std::sort` ([PR #3483](https://github.com/ERGO-Code/HiGHS/pull/3483)).

OpenBLAS is built with `-mfma` for its AVX2 `DYNAMIC_ARCH` kernels on
Linux ([PR #3399](https://github.com/ERGO-Code/HiGHS/pull/3399)), and on
Apple users can choose the BLAS library used for the HiPO build
([PR #3242](https://github.com/ERGO-Code/HiGHS/pull/3242)). CI now runs
more sanitizer tests on PRs
([PR #3474](https://github.com/ERGO-Code/HiGHS/pull/3474)).