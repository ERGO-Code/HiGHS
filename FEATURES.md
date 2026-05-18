The highlights of v1.15 are the first variant of the parallel MIP solver and the addition of HiPO in Python.

## Code changes


## Build changes

The Python build has been updated and an additional python package is available. The HiPO dependencies are linked via the optional `highspy-extras`, e.g. Metis, for all platforms, and OpenBLAS, for Windows and Linux. The `highspy-extras` package is automatically consumed by `highspy` and does not need to be imported manually. Note, that `highspy-extras` is distributed under the Apache 2.0 license, due to the dependencies\` licensing.

On Linux, `libblas` is no longer supported, in favour of OpenBLAS. MKL support would be considered at a later stage.

Prompted by [#2849](https://github.com/ERGO-Code/HiGHS/issues/3000), [PR #3027](https://github.com/ERGO-Code/HiGHS/pull/3027) fixes a recent `alpine:edge` compilation issue.

Following [PR #3025](https://github.com/ERGO-Code/HiGHS/pull/3025), all remaining thread sanitizer data race warnings are cleared.