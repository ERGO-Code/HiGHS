# PyPSA SciGRID-DE end-to-end baseline

End-to-end LP solve times for the parametric PyPSA SciGRID-DE workload
(`check/perf/scripts/build_scigrid_mps.py`). The 24h variant comes
verbatim from `pypsa.examples.scigrid_de()`; longer horizons replicate the
24h snapshot index to scale up the time-coupled LP.

Build: `cmake --build build --target highs-bin` (FAST_BUILD, Release, no LTO).
Host: Linux 6.18, x86_64.
HiGHS: 1.14.0, git hash `7df0786` (branch `claude/highs-cpu-optimization-review-XyO7z` at the perf-harness commit, no source-level optimisations applied yet).

| horizon | nvars     | ncons     | mps size | iters    | HiGHS time | wall time | objective       |
| ------- | --------- | --------- | -------- | -------- | ---------- | --------- | --------------- |
| 24h     | 59,640    | 142,968   |  14 MB   |  16,121  |   ~4.6 s   |   5.2 s   | 6.6848173236e+06 |
| 168h    | 417,480   | 1,000,776 |  99 MB   | 102,493  |  45.58 s   |  51.05 s  | 6.6848173236e+06 |
| 504h    | 1,252,440 | 3,002,328 | 298 MB   | 302,974  | 343.79 s   | 362.26 s  | 6.6848173236e+06 |

Reproduce:

```sh
# One-shot (fast, useful for iteration):
check/perf/scripts/run_endtoend.sh 168

# Full real-workload measurement:
check/perf/scripts/run_endtoend.sh 504
```

Use 504h as the "before / after" target when evaluating CPU optimisations
whose benefit only shows up over a long simplex run (~300k iterations of
FTRAN/BTRAN/PRICE/saxpy). 168h is fine for shake-down.

The objective value is identical across horizons because we replicate the
same 24h pattern — that's intentional (it makes correctness regressions
trivial to spot) and does *not* invalidate the timing signal: the simplex
work scales with `nvars × iters`, both of which grow with the horizon.

## Where the time is spent (168h SciGRID-DE, full solve)

Source: HiGHS internal timers via `highs_analysis_level=60`, full log at
`scigrid_168h_analysis.log`. Wall time of the analysis run was 48.8 s
(45.6 s in HiGHS proper).

### Simplex inner-loop breakdown (% of full solve)

| Bucket                    | Wall % | Calls   | µs/call | What it is                                  |
| ------------------------- | -----: | ------: | ------: | ------------------------------------------- |
| **INVERT (refactorise)**  | 17.4%  |      74 |  94 200 | Basis matrix LU factorisation               |
| **FTRAN_DSE**             | 14.2%  | 102 493 |    55.5 | Solve `B·x = a` for the entering column     |
| **FTRAN**                 | 13.8%  | 102 493 |    53.9 | Solve `B·x = a` (used in pivot/edge calc)   |
| **FTRAN_BFRT**            |  9.4%  |  60 808 |    62.2 | FTRAN for the bound-flipping ratio test     |
| **BTRAN**                 |  7.0%  | 102 502 |    27.3 | Solve `Bᵀ·y = c` (basis update)             |
| **PRICE**                 |  3.2%  | 102 502 |    12.7 | Compute `aᵀ·N` (reduced costs)              |
| **UPDATE_PRIMAL**         |  2.2%  | 409 951 |     2.2 | Apply the pivot to the primal solution      |
| **CHUZR_DUAL**            |  1.8%  | 102 496 |     6.9 | Choose pivot row (ratio test)               |
| **COMPUTE_DUAL**          |  1.3%  |      91 |  5 793  | Recompute dual values (rebuild)             |
| Everything else accounted |  6.3%  | -       | -       | CHUZC, UPDATE_*, COMPUTE_*, ratio cleanup   |
| **Total accounted**       | 76.6%  | -       | -       | Of full simplex `40.16 s`                   |

**Roll-ups:**

- **HFactor solves dominate.** FTRAN + FTRAN_DSE + FTRAN_BFRT + BTRAN = **44.4 %** of total wall time. Add INVERT at 17.4% and HFactor accounts for **>61 %** of the simplex run.
- **PRICE is small (3.2 %).** This is the dual simplex; PRICE is far cheaper than I had guessed pre-data. Optimisations on `priceByRow` will not move the needle much on this workload.
- **CHUZC family combined ≈ 1.9 %.** Tiny. Don't bother.

### HFactor breakdown (all factor calls, 168h)

| Sub-clock                    | Wall % | Calls   | µs/call |
| ---------------------------- | -----: | ------: | ------: |
| **FTRAN Upper**              | 32.4%  | 265 881 |    48.9 |
| **INVERT Kernel** (Markowitz)| 13.3%  |      74 |  72 400 |
| **FTRAN Upper FT** (FT update)|  8.3%  | 265 881 |    12.5 |
| **FTRAN Upper Hyper2**       | 10.9%  | 127 026 |    34.6 |
| **FTRAN Upper Hyper1**       |  8.0%  |  40 772 |    78.4 |
| **FTRAN Lower**              |  4.1%  | 265 881 |     6.3 |
| **FTRAN Lower Hyper**        |  3.4%  | 264 888 |     5.2 |
| **BTRAN Lower Hyper**        |  5.4%  | 102 502 |    21.3 |
| **INVERT Simple**            |  3.2%  |      74 |  17 200 |
| **INVERT Finish**            |  0.9%  |      74 |   4 700 |

**Surprises worth flagging:**

- **`FTRAN Upper` is 8× the cost of `FTRAN Lower`** (32.4 % vs 4.1 %). The U-factor of SciGRID-DE bases is much denser than L. Optimising `HFactor::ftranU` and the `solveHyper` it calls is **the single highest-impact target** on this workload.
- **`FTRAN Upper FT`** (Forrest-Tomlin update, 8.3 %) is its own substantial bucket. This is the patch-up step that lives in `HFactor::ftranFT` (`HFactor.cpp:1848+`), called between INVERT cycles. Worth its own attention.
- **The hyper-sparse paths dominate** (`FTRAN Upper Hyper1/2/3` = 21.3 % combined; `BTRAN Lower Hyper` = 5.4 %). The dense-fallback paths we benchmarked in `highs_perf` rarely fire on real workloads. **SIMD on the dense fallback wouldn't help; prefetch/`__restrict__` on the hyper-sparse loops would.**

### Function-level (callgrind, 24h full solve, ~24 G instructions)

Source: `valgrind --tool=callgrind`, full annotation at
`scigrid_24h_callgrind.txt`. Top self-time (instructions retired):

| % Ir   | Function                                                  |
| -----: | --------------------------------------------------------- |
| 14.4%  | `HFactor::ftranU` (combined call-stack contexts)          |
|  6.21% | `HFactor::ftranFT`                                        |
|  9.5%  | `solveHyper` (combined contexts)                          |
|  5.7%  | `HFactor::ftranL`                                         |
|  4.9%  | `HFactor::btranL` + `HFactor::btranU`                     |
|  3.76% | `HFactor::buildKernel` (INVERT)                           |
|  3.34% | `HEkkDualRHS::updatePrimal`                               |
|  3.13% | `HighsSparseMatrix::priceByRowWithSwitch`                 |
|  3.08% | `presolve::HPresolve::dualFixing` (one-shot startup cost) |
|  2.03% | `HFactor::buildSimple` (INVERT)                           |
|  1.59% | `HighsSparseMatrix::priceByColumn`                        |

(File-parse of MPS adds another ~3% via `__memchr_avx2` and string ops; not
relevant for repeated solves of an in-memory model.)

## Implication for the original optimisation review

The data **promotes** these items in the priority list:

1. **`HFactor::ftranU` and the hyper-sparse path inside `solveHyper`.** ~14% of all instructions; ~33% of wall time. This is now the #1 target.
2. **`HFactor::ftranFT`.** 6% of instructions; 8% of wall time. Different code path, same kind of inner loop.
3. **`HFactor::buildKernel`** (Markowitz pivoting). 13% of factor wall time. Different optimisation profile (dense work, scoring).

It **demotes** these:

- `priceByRow` and `priceByColumn` — dual simplex barely uses them (~3% each).
- The dense SpMV fallbacks — almost never hit.
- `HVectorBase::saxpy` standalone — embedded inside `solveHyper`, gets covered by the same fix.

LTO + `-march=native` (the cheap items) still help every function uniformly — keep them.

