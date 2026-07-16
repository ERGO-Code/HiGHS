# HiPO on PyPSA energy models: investigation notes and changes on this branch

This branch documents an investigation into why HiGHS's HiPO interior-point solver
underperforms Gurobi on one specific family of problems — full-year, hourly-resolution
PyPSA energy-system LPs — together with two verified, behaviour-preserving changes that
came out of it. Measurements were taken on a 4-core VM with the instances from the
public benchmark bucket of
[open-energy-transition/solver-benchmark](https://github.com/open-energy-transition/solver-benchmark)
(openenergybenchmark.org).

## Motivation

From the published benchmark results (best run per instance/solver, 16-vCPU VMs):

- HiPO is within a **geometric-mean factor 2.08 of Gurobi** across 40 co-solved PyPSA
  instances, and **beats Gurobi by 3–8x on sector-coupled models**.
- On **hourly-resolution electricity models** it collapses: 12–150x behind, with 9 of 12
  transmission-expansion (`trex`) variants timing out on instances Gurobi solves in
  minutes. This family is the gap.

## Diagnosis (measured on `pypsa-de-elec-10-1h`, 4.5M rows raw / 1.17M presolved)

- The normal-equations form fails outright on this structure (forcing it reports an
  integer overflow in the factor); the augmented-system fallback works as designed.
- `metis` ordering (the default choice) is right: `amd` gives 8x the flops, `rcm`
  exhausts memory.
- The elimination tree is chain-shaped (year-long time coupling), so supernodes are
  small: the factorisation runs at only ~1.3 GFLOP/s, and tree-level parallelism has
  little to work with (2 -> 4 threads: factorisation 165s -> 126s, solves unchanged).
- **The triangular solve phase is 55–65% of the numerical work at every problem size**
  (solve/factorise time ratio ~1.5–1.8 from 95k to 1.17M rows) and is entirely serial —
  an Amdahl ceiling of ~2x on multicore hardware until it is parallelised. The HiPO
  paper (Zanetti & Gondzio, arXiv:2508.04370) lists exactly this as future work.
- There are ~8 triangular solves per IPM iteration. With `HIPO_TIMING_LEVEL 2`, the
  solve phase splits into: dense BLAS ~44%, **pivot-swap permutations ~22%**, sparse
  scatter ~4%, and ~30% overhead (dominated by per-block heap allocations).
- Iteration counts are inflated by data scaling: these instances carry 1e10 bounds, and
  `user_bound_scale=-14` reached a 2.5x smaller duality gap in identical wall time with
  refinement residuals improving from ~1.0 to ~1e-13.

Suggested priority for closing the hourly-family gap:
1. scaling defaults for energy-model data (an existing option; also fixable upstream in
   PyPSA/linopy by emitting true infinities instead of 1e10 big-M bounds);
2. parallelising the triangular solves (the acknowledged engineering gap);
3. supernode-amalgamation and regularisation tuning (compounds with both).

## Changes on this branch

### 1. `193d2bcf` — expose FactorHiGHS tuning constants as advanced runtime options

The tuning parameters in `highs/ipm/hipo/factorhighs/FactorHiGHSSettings.h` (supernode
amalgamation thresholds, Bunch-Kaufman pivot threshold, parallel-gemm block sizes,
dynamic regularisation coefficient, Metis seed) were compile-time constants, so
experimenting required a rebuild per candidate value. They now live in a
`hipo::HipoTuning` struct populated once during single-threaded setup from twelve new
advanced options (`hipo_sn_*`, `hipo_spops_weight`, `hipo_alpha_bk`, `hipo_block_*`,
`hipo_min_consecutive_sums`, `hipo_dynamic_reg_coeff`, `hipo_metis_seed`).

Verified: defaults reproduce the previous behaviour exactly (same iterations, objective
and flop count); changed values demonstrably reach the factorisation (amalgamation flops
move); out-of-range values are rejected by option bounds.

### 2. `8107f6e4` — serial hygiene in the triangular solve phase

Three behaviour-preserving fixes in `highs/ipm/hipo/factorhighs/`:

- **Workspace hoist**: the solve handlers allocated a `std::vector<double>` per column
  block per supernode per solve (millions of allocations per solve at 1M+ rows). A
  single workspace, sized once per factorisation in `Numeric::finaliseFactor()`, is now
  reused.
- **Per-block swap-skip flags**: both `permuteWithSwaps` passes ran on every block of
  every supernode on every solve, even when pivoting never swapped anything. A per-block
  flag computed once per factorisation skips them.
- **`solveAS` buffer reuse**: the augmented-system rhs concatenation reused no memory;
  it now uses a persistent buffer.

Verified **bit-identical**: iterate-by-iterate logs match the baseline exactly on a full
solve of `pypsa-eur-elec-op-2-3h` (35 iterations) and `pypsa-eur-elec-op-10-3h`
(98 iterations); on `pypsa-de-elec-10-1h` the 41-iteration baseline prefix is identical.
Effect: solve time -7% (small), -14% (medium); on the hourly instance, per-solve cost
-12% and **48 vs 41 iterations within the same 600s budget**. Per-factorisation time
unchanged (3.71s -> 3.70s), confirming the change is solve-local.

## References

- F. Zanetti, J. Gondzio, *A factorisation-based regularised interior point method
  using the augmented system*, arXiv:2508.04370 (the HiPO paper; benchmarks include 45
  PyPSA-Eur problems and name solve-phase parallelisation as future work).
- J. Hogg, J. Reid, J. Scott, *A DAG-based sparse Cholesky solver for multicore
  architectures*, SIAM J. Sci. Comput., 2010.
- J. Park, M. Smelyanskiy, P. Dubey, *Sparsifying Synchronization for High-Performance
  Shared-Memory Sparse Triangular Solver*, ISC 2014.
