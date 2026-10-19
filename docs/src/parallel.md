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
HiGHS uses one thread per physical core available to your program (see [here](@ref scheduler-details)).
If it is set to one, then HiGHS will never use more than one thread. 
The maximum value that is advantageous is machine-dependent, 
but it is unlikely to be more than eight due to most computation in 
HiGHS being memory-bound.

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

## [MIP](@id parallel-mip)

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



## Parallel scheduler: initialization, shutdown and thread safety

### [Terminology](@id parallel-terminology)

- **Application thread**: a thread you create yourself (with `std::thread`, OpenMP, pthreads, Julia, Python, and so on), or the main thread of your program.
- **Scheduler**: the pool of threads HiGHS uses to run parallel work.
- **Worker thread**: a thread created by HiGHS when it creates a scheduler. Worker threads only run HiGHS tasks.
- **Owning thread**: the application thread that caused a scheduler to be created. It takes part in the parallel work together with the worker threads.

### [Each application thread has its own scheduler](@id scheduler-details)

A scheduler is tied to the application thread that created it, not to a `Highs` instance.

- The first time an application thread calls `run()` or `presolve()`, HiGHS creates a scheduler for that thread. The thread becomes its owning thread.
- All `Highs` instances used on the same application thread share that thread's scheduler.
- Application threads never share a scheduler. If several application threads each use HiGHS, each gets its own scheduler and its own worker threads. They do not interfere with each other.
- By default, HiGHS uses one thread per physical core that your program is allowed to run on (its CPU affinity). Hyperthreads (extra "virtual" cores on the same physical core) are not counted. The thread that calls HiGHS is included in this number, so a scheduler using `N` threads creates `N - 1` worker threads.
- The number of worker threads is fixed when the scheduler is created. Changing the `threads` option later does not resize an existing scheduler. To change it, reset the scheduler first (see [Shutdown](#shutdown)).

Because every scheduler adds worker threads, running HiGHS from many application threads with the default thread count can create far more threads than you have cores, which slows everything down. In that case set the `threads` option explicitly.

### Shutdown

A scheduler is shut down when its owning thread exits, or when you call `Highs::resetGlobalScheduler(bool blocking)`. Without an explicit call, a scheduler stays alive between solves, so that it can be reused.

The `blocking` argument controls what happens to the worker threads:

- `resetGlobalScheduler(true)` waits until all worker threads have finished before returning. Use it when you need the threads to be gone, for example before unloading a library or exiting a program that checks for leftover threads.
- `resetGlobalScheduler(false)` (the default) returns immediately. The worker threads are told to stop and finish on their own shortly afterwards, so a few of them may still be running when the call returns.

When a scheduler is shut down because its owning thread exits, the workers are stopped without waiting for them, as with `false`.

Resetting only affects the scheduler of the calling thread. If the calling thread has no scheduler, the call does nothing. Other threads' schedulers are not affected.

### Thread safety

Call `resetGlobalScheduler()` only:

- from the thread that owns the scheduler (the one that ran the solve), and
- when no `run()` or `presolve()` is in progress on that thread, including from inside a callback.

Do not call it from inside a HiGHS parallel task or from a worker thread. That would stop the whole scheduler while other workers may still be running tasks.

After a reset, the next `run()` or `presolve()` on that thread creates a new scheduler.

### Using the C API

The C API follows the same rules. Each application thread that calls `Highs_run` (or any other function that solves a model) gets its own scheduler, and `Highs_destroy` does not shut it down. Instances created and used on the same thread share that thread's scheduler. To control the number of threads, set the `threads` option with `Highs_setIntOptionValue(highs, "threads", n)` before the first solve on a thread. To free the scheduler explicitly, call `Highs_resetGlobalScheduler(1)` to wait for the worker threads to finish, or `Highs_resetGlobalScheduler(0)` to return immediately. The same restrictions as in C++ apply: call it from the thread that ran the solve, when no solve is in progress on that thread (including from inside a callback), and never from a worker thread or a HiGHS task.

