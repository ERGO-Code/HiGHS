# HiGHS Original vs Lagromory Solver Benchmark Comparison

This report presents a head-to-head comparison between the original HiGHS MIP solver and the new native Lagromory (Relax-and-Cut) enhanced solver on lightweight MIPLIB instances (solved with a 60-second time limit).

## Comparison Table

| Instance | Solver Version | Status | Primal Bound (Obj) | Dual Bound | Nodes | LP Iters | Time (s) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **markshare_4_0** | Original | Time | 1 | 0 | 1028230 | 2716222 | 60.00 |
| | Lagromory | Time | 1 | 0 | 988186 | 2504532 | 60.00 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **pk1** | Original | Time | 11 | 3.73624327435 | 233305 | 2341806 | 60.00 |
| | Lagromory | Time | 13.999999 | 3.60580235899 | 239564 | 2291166 | 60.00 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **neos5** | Original | Optimal | 14.9999999998 | 14.9987951807 | 216997 | 2663178 | 56.25 |
| | Lagromory | Time | 14.9999999993 | 14.1436035258 | 184909 | 2944017 | 60.00 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **mas76** | Original | Time | 40005.054141 | 39691.7546079 | 275261 | 1942578 | 60.00 |
| | Lagromory | Time | 40005.054142 | 39670.1085081 | 282289 | 1933943 | 60.00 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **assign1-5-8** | Original | Time | 212 | 198 | 32494 | 1152774 | 60.00 |
| | Lagromory | Time | 212 | 198 | 54164 | 1380274 | 60.00 |
| --- | --- | --- | --- | --- | --- | --- | --- |

## Saved Artifacts
The `.sol` (solution) and `.log` (run log) files have been saved to `/home/yunus/Masaüstü/Projects/HiGHS/benchmark_results/` for all runs.
