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
