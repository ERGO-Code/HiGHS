# [Parallelism](@id parallelism)

## Generally

HiGHS has increasing opportunities for exploiting parallel
computing. When using a CPU, these are currently restricted to the
dual simplex solver for LP, the factorisation-based interior point solver (HiPO),
and the MIP solver. Details of these and future plans are set out below.
HiGHS has an implementation of a first order method (PDLP) for solving LPs
that can exploit the availability of a [GPU](@ref gpu).

The value of the [threads](@ref option-threads) option determines the 
number of threads used by HiGHS. By default it is zero, in which case 
HiGHS will use half the available threads on a machine. 
If it is set to one, then HiGHS will never use more than one thread. 
The maximum value that is advantageous is machine-dependent, 
but it is unlikely to be more than eight due to most computation in 
HiGHS being memory-bound.

Note that since parallelism is controlled by a static scheduler,
concurrent Highs instances must use the same value of `threads`.

## Dual simplex

By default, the HiGHS dual simplex solver runs in serial. However, it
has a variant allowing concurrent processing. This variant is used
when the [parallel](@ref option-parallel) option is set to "on" and
[simplex\_strategy](@ref option-simplex-strategy) is set to
`kSimplexStrategyDualTasks` (2) or `kSimplexStrategyDualMulti` (3).

The concurrency used will be the value of
[simplex\_max\_concurrency](@ref option-simplex-max-concurrency). If
this is fewer than the number of threads available, parallel
performance may be less than anticipated.

The speed-up achieved using the dual simplex solver is normally
bounded by the number of memory channels in the architecture, and
typically less than the values achieved by [Huangfu and
Hall](https://link.springer.com/article/10.1007/s12532-017-0130-5). This
is because enhancements to the serial dual simplex solver in recent
years have not been propagated to the parallel solver.

Unless an LP has significantly more variables than constraints, the
parallel dual simplex solver is unlikely to be worth using.

## IPM

The interior point solver HiPO uses multi-threading in various phases:

* To run multiple orderings heuristics on different Newton system approaches, 
  in order to choose the best one.
* To construct the normal equations matrix. 
* To process the elimination tree during the multifrontal factorisation 
  (_tree level_ parallelism).
* To perform the dense factorisation of the frontal matrices (_node level_ parallelism).
* To perform the triangular solves.

Running multiple orderings and building the normal equations in parallel is always 
advantageous, so is performed regardless of the value of the 
[parallel](@ref option-parallel) option.

If the [parallel](@ref option-parallel) option is set "on" or "choose", the parallelism 
during the factorisation and the triangular solves is controlled by a heuristic. 
Otherwise, these phases do not exploit parallelism.

If the [parallel](@ref option-parallel) option is set "on", the level of parallelism 
in the factorisation can be refined using the 
[hipo\_parallel\_type](@ref option-hipo-parallel-type) option, which can be "tree" for 
tree level only, "node" for node level only, "both" for both levels, or "choose" to 
leave the selection to the solver.

The options [hipo\_parallel\_force](@ref option-hipo-parallel-force) and 
[hipo\_parallel\_forbid](@ref option-hipo-parallel-forbid) can be used to override the 
default parallelism settings. They are bit maps in which the various parallel phases are 
controlled by the following values:
-   1: analyse phase
-   2: reordering of normal equations
-   4: reordering of augmented system
-   8: building of normal equations structure
-  16: building of normal equations values
-  32: processing of the elimination tree
-  64: dense factorisation of frontal matrices
- 128: forward solve with factorisation
- 256: diagonal solve with factorisation
- 512: backward solve with factorisation

Setting `hipo_parallel_force` (resp. `hipo_parallel_forbid`) to one of these values, 
or a sum of values, forces (resp. forbids) the use of parallelism in the corresponding 
phases. These options override any other behaviour enforced by other options. If a given 
phase is both forced and forbidden, the default behaviour is used instead. 
For instance, setting `hipo_parallel_force` to 81 = 1+16+64 and `hipo_parallel_forbid` 
to 68 = 4+64  forces the use of parallelism in the analyse phase and for building the 
normal equations values, and forbids it for the reordering of augmented system. 

By default, `hipo_parallel_type` is set to `choose`, `hipo_parallel_force` and 
`hipo_parallel_forbid` are set to 0, meaning that the choice of parallelism is left 
completely to the internal heuristics. These heuristics are very much likely to be 
correct, but the options provide ways of overriding them.

The extent to which parallelism is used in HiPO depends on the value of the
[threads](@ref option-threads) option (see above).

### [BLAS threads](@id parallel-blas-threads)

HiPO calls BLAS routines for the dense linear algebra in the factorisation
and solves. Many BLAS libraries are multi-threaded, and their threads are
independent of the HiGHS [threads](@ref option-threads) option. Since HiPO
already exploits parallelism at a higher level, also running OpenBLAS with
multiple threads oversubscribes the cores and can slow HiPO down
considerably, as well as making its behaviour non-deterministic. HiPO
performs best when OpenBLAS runs in serial.

When HiGHS is built with OpenBLAS, HiPO therefore sets the number of OpenBLAS
threads to one before solving, and stops with an error if OpenBLAS does not
then report using a single thread. This also applies when OpenBLAS is given
directly via `BLAS_LIBRARIES` (see [Installation](@ref install-blas)), as long as
the library name contains "openblas".

The number of OpenBLAS threads is a __process-wide__ setting, so it affects
any other component of the same process that uses the same OpenBLAS library.
If HiGHS is one part of a larger application that relies on multi-threaded
OpenBLAS, set the
[hipo\_serial\_blas](@ref option-hipo-serial-blas) option to false. HiPO
then leaves the number of OpenBLAS threads unchanged, and it is up to the
user to choose it, for instance with the `OPENBLAS_NUM_THREADS` environment
variable or by calling `openblas_set_num_threads`.

HiPO only changes the number of OpenBLAS threads if it is not already one, and
HiGHS instances in the same process never change it concurrently. However,
`openblas_set_num_threads` is not thread safe, so other code in the process
should not change the number of OpenBLAS threads while HiPO is running.

For other BLAS libraries, HiPO does not change the number of threads:

- Apple Accelerate: HiPO adjusts its own parallelism to suit Accelerate, so no
  action is needed.
- libblastrampoline (used by HiGHS.jl, see [HiGHS.jl](@ref HiGHS.jl)): if using OpenBLAS, call
  `BLAS.set_num_threads(1)` in Julia before solving with HiPO.

## MIP

If the [parallel](@ref option-parallel) option is not set to "off", the MIP solver
will explore the branch-and-bound tree using multiple threads.
This exploration includes cuts and heuristics that are not run from the root node.

In addition, the MIP solver utilises parallelism when performing 
symmetry detection on the model, when querying
clique tables, and when the interior point solver is used to compute
the analytic centre. This parallelism is always advantageous, so is
performed regardless of the value of the [parallel](@ref option-parallel) option.

The extent to which parallelism is used in the MIP solver depends on the value of the 
[threads](@ref option-threads) option (see above).

## Future plans

First-order solvers for LP are still very much in their infancy, and
are not robust. Hence the availability of a PDLP solver for LP is
unlikely to be used to enhance other solvers in HiGHS in the short or
medium term.

## Alternative

Given the limited scope for parallelisation in HiGHS, if you have
multiple independent instances to solve, then assign one instance per
worker (cpu core, thread or machine) so that multiple instances are
solved concurrently.
