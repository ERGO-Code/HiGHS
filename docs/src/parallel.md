# Parallelism

## Generally

HiGHS currently has limited opportunities for exploiting parallel
computing. These are currently restricted to the dual simplex solver
for LP, and the MIP solver. Details of these and future plans are set
out below.

By default, when running in parallel, HiGHS will use half the
available threads on a machine. This number can be modified by setting
the value of the
[threads](@ref)
option.

## Dual simplex

By default, the HiGHS dual simplex solver runs in serial. However, it
has a variant allowing concurrent processing. This variant is used
when the
[parallel](@ref)
option is set "on", by specifying `--parallel` when running the
[executable](@ref Executable) via
the command line, or by setting it via a library call in an
application.

The concurrency used will be the value of
[simplex\_max\_concurrency](@ref). If
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

## MIP

The only parallel computation currently implemented in the MIP solver
occurs when performing symmetry detection on the model, and when
querying clique tables. This parallelism is always advantageous, so is
performed regardless of the value of the
[parallel](@ref) option.

## Future plans

The MIP solver has been written with parallel tree seach in mind, and
it is hoped that this will be implemented before the end of 2024. The
parallel LP solver will also enhance the MIP solver performance by
speeding up the solution of the root node.

Development of a parallel interior point solver will start in 2023,
and is expected to be completed by the end of 2024.



