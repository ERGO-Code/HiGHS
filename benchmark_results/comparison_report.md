# HiGHS Original vs Cache-Optimized Solver Benchmark Comparison

This report presents a head-to-head comparison between the original HiGHS MIP solver and the L1/L2 cache-optimized version (prefetch + loop unrolling) on lightweight MIPLIB instances (solved with a 60-second time limit).

## Comparison Table

| Instance | Solver Version | Status | Primal Bound (Obj) | Dual Bound | Nodes | LP Iters | Time (s) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **markshare_4_0** | Original | Time | 1 | 0 | 1041394 | 2751681 | 60.00 |
| | Optimized | Time | 1 | 0 | 1039030 | 2744680 | 60.00 |
| | **Speedup** | | | | | | **1.00x** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **pk1** | Original | Time | 11 | 3.73843291027 | 233638 | 2345601 | 60.00 |
| | Optimized | Time | 16 | 3.36725306104 | 242758 | 2218315 | 60.00 |
| | **Speedup** | | | | | | **1.00x** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **neos5** | Original | Optimal | 14.9999999998 | 14.9987951807 | 216997 | 2663178 | 56.26 |
| | Optimized | Unknown | N/A | N/A | 0 | 0 | 0.00 |
| | **Speedup** | | | | | | **0.00x** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **mas76** | Original | Unknown | N/A | N/A | 0 | 0 | 0.00 |
| | Optimized | Unknown | N/A | N/A | 0 | 0 | 0.00 |
| | **Speedup** | | | | | | **0.00x** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **assign1-5-8** | Original | Unknown | N/A | N/A | 0 | 0 | 0.00 |
| | Optimized | Unknown | N/A | N/A | 0 | 0 | 0.00 |
| | **Speedup** | | | | | | **0.00x** |
| --- | --- | --- | --- | --- | --- | --- | --- |

## Saved Artifacts
The `.sol` (solution) and `.log` (run log) files have been saved to `/home/yunus/Masaüstü/Projects/HiGHS/benchmark_results/` for all runs.
