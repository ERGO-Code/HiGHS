Patch release v1.15.1 fixes:

Apache release binaries not including HiPO fixed in [PR #3112](https://github.com/ERGO-Code/HiGHS/pull/3112).

Build for users with C++20 fixed in [PR #3114](https://github.com/ERGO-Code/HiGHS/pull/3114).

OpenBLAS setup issue affecting performance fixed in [PR #3122](https://github.com/ERGO-Code/HiGHS/pull/3122).

---


The highlights of v1.15 are the first variant of the parallel MIP solver
and the addition of HiPO in Python.

## Code changes

Following [PR #2886](https://github.com/ERGO-Code/HiGHS/pull/2886), our
prototype multithreaded MIP solver is available, refactoring
the worker/node-search logic of the branch-and-cut solver to run
across multiple threads. This is the first release of the parallel MIP
solver.

A fix to the parallel scheduler was merged in
[PR #3087](https://github.com/ERGO-Code/HiGHS/pull/3087).

Following [PR #2994](https://github.com/ERGO-Code/HiGHS/pull/2994),
HiPO can now be used to solve LPs and QPs from `highspy`, with
documentation added in
[PR #3024](https://github.com/ERGO-Code/HiGHS/pull/3024) and
[PR #3060](https://github.com/ERGO-Code/HiGHS/pull/3060). HiPO
performance was further improved by faster handling of free variables
([PR #3013](https://github.com/ERGO-Code/HiGHS/pull/3013)) and faster
triangular solves
([PR #3014](https://github.com/ERGO-Code/HiGHS/pull/3014)).

Prompted by [#2957](https://github.com/ERGO-Code/HiGHS/issues/2957),
square Hessians that are only slightly asymmetric are now accepted by
`Highs::qFormatOk`, rather than being rejected outright
([PR #2965](https://github.com/ERGO-Code/HiGHS/pull/2965),
[PR #2984](https://github.com/ERGO-Code/HiGHS/pull/2984) and
[PR #3003](https://github.com/ERGO-Code/HiGHS/pull/3003)). Identification
of non-convexity in the QP solvers was also improved, with extended
documentation added in
[PR #3056](https://github.com/ERGO-Code/HiGHS/pull/3056), and a further
bug affecting some QPs was fixed in
[PR #3069](https://github.com/ERGO-Code/HiGHS/pull/3069). QP hot start
was added in [PR #3089](https://github.com/ERGO-Code/HiGHS/pull/3089).

Following [PR #2982](https://github.com/ERGO-Code/HiGHS/pull/2982),
several bugs in the computation of an Irreducible Infeasible Subsystem
(IIS) via `Highs::getIis` were fixed.

Following [PR #3046](https://github.com/ERGO-Code/HiGHS/pull/3046),
zero-cost singleton columns are no longer fixed to an infinite bound
during `HPresolve::dualFixing`, and the size of shifts applied during
reduced-cost fixing was reduced in
[PR #2986](https://github.com/ERGO-Code/HiGHS/pull/2986). Strengthened
variable bounds found during presolve are now retained
([PR #3010](https://github.com/ERGO-Code/HiGHS/pull/3010)), a bug in
the computation of the fractional value of a variable with finite
lower and upper bounds was fixed
([PR #3005](https://github.com/ERGO-Code/HiGHS/pull/3005)), and the
handling of infeasibilities for semi-continuous and semi-integer
variables was corrected
([PR #3018](https://github.com/ERGO-Code/HiGHS/pull/3018)).

Following [PR #3042](https://github.com/ERGO-Code/HiGHS/pull/3042),
parallel simplex is no longer used when solving the LP relaxations
that arise within the MIP solver, since it isn't safe to do so when
the MIP solver is itself running in parallel.

Following [PR #2971](https://github.com/ERGO-Code/HiGHS/pull/2971),
`highspy` now supports scalar division, and significant improvements
were made to `highspy`'s static typing in
[PR #2983](https://github.com/ERGO-Code/HiGHS/pull/2983). Context
manager support was added to `highspy` in
[PR #3093](https://github.com/ERGO-Code/HiGHS/pull/3093).

Following [PR #3039](https://github.com/ERGO-Code/HiGHS/pull/3039),
`Highs::setBasis` and `Highs::setLogicalBasis` have been implemented,
and additional methods were added to the C# wrapper in
[PR #3041](https://github.com/ERGO-Code/HiGHS/pull/3041).

Prompted by [#3044](https://github.com/ERGO-Code/HiGHS/issues/3044),
a further bug was fixed in
[PR #3092](https://github.com/ERGO-Code/HiGHS/pull/3092), and
[PR #3091](https://github.com/ERGO-Code/HiGHS/pull/3091) corrects the
re-checking of an implied bound, fixing
[#3090](https://github.com/ERGO-Code/HiGHS/issues/3090).

Following [PR #2981](https://github.com/ERGO-Code/HiGHS/pull/2981),
C-heap memory leaks in solver state cleanup were fixed, and
[PR #3081](https://github.com/ERGO-Code/HiGHS/pull/3081) fixes a bug
in column stuffing.

Following [PR #3016](https://github.com/ERGO-Code/HiGHS/pull/3016),
`Highs::presolve()` now logs and returns failure appropriately when
called inconsistently with the state of the incumbent model.

Prompted by [#3007](https://github.com/ERGO-Code/HiGHS/issues/3007), a
potential typo in `FactorHiGHSSolver::chooseNla()` was corrected, and
several other minor fixes were collected in
[PR #3017](https://github.com/ERGO-Code/HiGHS/pull/3017) and
[PR #3004](https://github.com/ERGO-Code/HiGHS/pull/3004) (the latter
for `HighsCliqueTable`).

Following [PR #3057](https://github.com/ERGO-Code/HiGHS/pull/3057),
lines in HiGHS options files that contain only spaces are now ignored,
and the text of the first line containing an error is printed in the
resulting message.

## Build changes

The Python build has been updated and an additional python package is
available. The HiPO dependencies are linked via the optional
`highspy-extras`, e.g. Metis, for all platforms, and OpenBLAS, for
Windows and Linux. The `highspy-extras` package is automatically
consumed by `highspy` and does not need to be imported manually. Note,
that `highspy-extras` is distributed under the Apache 2.0 license, due
to the dependencies' licensing.

On Linux, `libblas` is no longer supported, in favour of OpenBLAS. MKL
support would be considered at a later stage.

Prompted by [#3000](https://github.com/ERGO-Code/HiGHS/issues/3000),
[PR #3027](https://github.com/ERGO-Code/HiGHS/pull/3027) fixes a
recent `alpine:edge` compilation issue.

Following [PR #3025](https://github.com/ERGO-Code/HiGHS/pull/3025),
all remaining thread sanitizer data race warnings are cleared.

GPU support was added for local Python pip installs in
[PR #3071](https://github.com/ERGO-Code/HiGHS/pull/3071) and the BLAS
library used by HiPO is now checked at runtime
([PR #3080](https://github.com/ERGO-Code/HiGHS/pull/3080)).

The `cibuildwheel` workflow and `pyproject` configuration were updated
to `cibuildwheel` 4.0.0
([PR #3061](https://github.com/ERGO-Code/HiGHS/pull/3061)), and the
NuGet publishing workflow now uses trusted publishing rather than an
API key
([PR #3062](https://github.com/ERGO-Code/HiGHS/pull/3062),
following [PR #2959](https://github.com/ERGO-Code/HiGHS/pull/2959)).

The macOS GitHub Actions build times were significantly improved in
[PR #3033](https://github.com/ERGO-Code/HiGHS/pull/3033), and the C#
wrapper build was fixed to statically link the MSVC STL and correctly
install on win32
([PR #3070](https://github.com/ERGO-Code/HiGHS/pull/3070)).

Following [PR #3001](https://github.com/ERGO-Code/HiGHS/pull/3001),
the install location of the readme and license files was corrected,
and explicit `-fexceptions` copts were added to the Bazel config
([PR #3072](https://github.com/ERGO-Code/HiGHS/pull/3072) and
[PR #3075](https://github.com/ERGO-Code/HiGHS/pull/3075)).
