# [Solvers](@id solvers)

## Introduction

HiGHS has implementations of the three main solution techniques for LP
(simplex, interior point and primal-dual hybrid gradient), two
solution techniques for QP (active set and interior point), and one
MIP solver. By default HiGHS will choose the most appropriate
technique for a given problem, but this can be over-ridden by setting
the option [__solver__](@ref option-solver).

What follows below is an overview of all the solvers and key
[options](@ref option-definitions) that define their algorithmic
features. The relevant [__solver__](@ref option-solver) settings for
each problem type, and their interpretation when incompatible with the
problem type, are then [summarised](@ref solver-option).

## LP

#### Simplex

HiGHS has efficient implementations of both the primal and dual
simplex methods, although the dual simplex solver is likely to be
faster and is more robust, so is used by default. Similarly, for
reasons discussed in the [parallelism](@ref) section, the parallel
dual simplex solver is unlikely to be worth using. The novel features
of the dual simplex solver are described in

_Parallelizing the dual revised simplex method_, Q. Huangfu and
J. A. J. Hall, Mathematical Programming Computation, 10 (1), 119-142,
2018 [DOI:
10.1007/s12532-017-0130-5](https://link.springer.com/article/10.1007/s12532-017-0130-5).

* Setting the option [__solver__](@ref option-solver) to "simplex" forces the simplex solver to be used
* The option [__simplex\_strategy__](@ref option-simplex-strategy)
  determines whether the primal solver or one of the parallel solvers is
  to be used.

#### [Interior point](@id solvers-lp-ipm)

HiGHS has two interior point (IPM) solvers:

* IPX is based on the preconditioned conjugate gradient method, as discussed in

  _Implementation of an interior point method with basis
  preconditioning_, Mathematical Programming Computation, 12, 603-635,
  2020. [DOI:
  10.1007/s12532-020-00181-8](https://link.springer.com/article/10.1007/s12532-020-00181-8).

  This solver is serial.

* HiPO is based on a direct factorisation, as discussed in 

  _A factorisation-based regularised interior point method using the
  augmented system_, F. Zanetti and J. Gondzio, 2025, [available on
  arxiv](https://arxiv.org/abs/2508.04370)

  This solver is parallel.

  The [__hipo\_system__](@ref option-hipo-system) option can be used to
  select the approach to use when solving the Newton systems within
  the interior point solver: select "augmented" to force the solver to
  use the augmented system, "normaleq" for normal equations, or
  "choose" to leave the choice to the solver.

  The option [__hipo\_ordering__](@ref option-hipo-ordering) can be used
  to select the fill-reducing heuristic to use during the
  factorisation:
  
  * Nested dissection, obtained setting the option
    [__hipo\_ordering__](@ref option-hipo-ordering) to "metis".
  
  * Approximate minimum degree, obtained setting the option
    [__hipo\_ordering__](@ref option-hipo-ordering) to "amd".
  
  * Reverse Cuthill-McKee, obtained setting the option
    [__hipo\_ordering__](@ref option-hipo-ordering) to "rcm".

For small LPs, IPX is often faster than HiPO. However, as problem size
grows, HiPO becomes more efficient, and its advantage can be more
than an order of magnitude.

#### Primal-dual hybrid gradient method

HiGHS has a primal-dual hybrid gradient implementation for LP (PDLP)
that is run on an NVIDIA [GPU](@ref gpu) if CUDA is installed. If this
is not possible, the PDLP solver is run on a CPU, but it is unlikely
to be competitive with the HiGHS interior point or simplex solvers.

## QP

HiGHS has two solvers for convex QP: a primal active set method, and
an interior point method. The active set implementation uses a dense
Cholesky factorization of the reduced Hessian, and the the limit on
its dimension is determined by the option
[__qp\_nullspace\_limit__](@ref option-qp-nullspace-limit). The
interior point solver is HiPO, so see [above](@ref solvers-lp-ipm) for
the key algorithmic options.

## MIP

MIPs are solved using a sophisticated branch-and-cut solver that is
still largely single-threaded. Users can choose the solver for LP
sub-problems as follows

* LPs where an advanced basis is not known. This is generally the case
  at the root node of the branch-and-bound tree, but can occur at
  other times in the MIP solution process. By default the simplex
  solver is used, but the [__mip\_lp\_solver__](@ref
  option-mip-lp-solver) option allows the user to determine whether
  one of the IPM solvers is used.

* LPs where an interior point solver must be used. This is generally
  when computing the analytic centre of the constraints for the
  centralised rounding heuristic. By default IPX is used, but the
  [__mip\_ipm\_solver__](@ref option-mip-ipm-solver) option allows the
  user to determine whether HiPO is used.

Otherwise, for LPs where an advanced basis is known, only the simplex
solver can be used.

## [Summary](@id solver-option)

The option [__solver__](@ref option-solver) can be set to:
* "simplex", which selects the simplex solver.
* "ipm", which selects the HiPO solver (or IPX if HiPO is not available in the build).
* "ipx", which selects the IPX solver.
* "hipo", which selects the HiPO solver, for both LP and QP.
* "pdlp", which selects the PDLP solver.
* "qpasm", which selects the QP active-set method.
* "choose", which selects the default solver for the given problem ("simplex" for LP, "qpasm" for QP).

The option [__solver__](@ref option-solver) is ignored and the default solver is used if:
* The problem is an LP and solver is set to "qpasm".
* The problem is a QP and solver is set to "simplex", "ipx" or "pdlp".
* The problem is a MIP and solver is not set to "choose".

