
In order to use OpenMP if available, set `-DOPENMP=ON` during the configuration
step ( `cmake ..` ).

When compiled with the parallel option on, the number of threads used at run
time is the value of the environment variable `OMP_NUM_THREADS` . For example, 
to use HiGHS with eight threads to solve `ml.mps` execute

    export OMP_NUM_THREADS=8
    highs --parallel ml.mps

If `OMP_NUM_THREADS` is not set, either because it has not been set or due to
executing the command

    unset OMP_NUM_THREADS

then all available threads will be used.

If run with `OMP_NUM_THREADS=1` , HiGHS is serial. The `--parallel` run-time
option will cause HiGHS to use serial minor iterations and, although this
could lead to better performance on some problems, performance will typically be
diminished.

When compiled with the parallel option and `OMP_NUM_THREADS>1` or unset, HiGHS
will use multiple threads. If `OMP_NUM_THREADS` is unset, HiGHS will try to use
all available threads so performance may be very slow. Although the best value
will be problem and architecture dependent, `OMP_NUM_THREADS=8` is typically a
good choice. Although HiGHS is slower when run in parallel than in serial for
some problems, it is typically faster in parallel.