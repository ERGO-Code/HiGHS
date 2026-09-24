# HiPDLP Phase 7 Cross-Solver Validation Report

Date: 2026-09-04

Status: **COMPLETE. Phase 7 correctness validation is closed.**

## Findings

1. The canonical validation matrix is complete: 39 smoke rows, 824 Netlib
   regression rows, 912 full-Netlib rows, 32 energy-regression rows, and 72
   full-energy rows. `summary.csv` contains 1,879 unique
   scope/solver/instance rows.
2. No new rewrite correctness defect was found while filling the missing
   matrix. Every executed rewrite `OPTIMAL` result has complete finite
   postsolve/KKT data and independently recomputed original-LP primal, dual,
   and gap metrics at or below `1e-4`.
3. Full Netlib accepts rewrite CPU on 111/114 and rewrite CUDA on 107/114.
   Every rejection is a 300-second `TIME_LIMIT`. The four CPU/CUDA terminal
   differences (`d2q06c`, `fit1d`, `greenbeb`, and `sierra`) are CUDA time
   limits against CPU optimal results and are classified
   `TIME_OR_ITERATION_LIMIT`.
4. Energy regression accepts rewrite CPU and CUDA on 4/4. Their statuses and
   iteration counts match on every instance, and their largest paired
   objective relative difference is `5.84e-14`.
5. Full energy was run for all nine instances and all eight solver
   configurations. Rewrite CPU and CUDA each accept 6/9. Both time out on
   `zen-garden-eur-PI-28-2ts` and `pypsa-eur-elec_16_1h_2050`; both hit the
   harness process timeout on `pypsa-de-elec-trex_vopt_50_1h`. Missing timeout
   metrics are rejected, never counted as success.
6. No run was refused for memory or disk. The preflight found 3.696 TB disk
   available, 120,252,540 kB host memory available, 16,777,212 kB free swap,
   and 83.138 GiB of 121.689 GiB CUDA-visible unified memory free. The largest
   observed rewrite process used about 8.1 GiB RSS on CPU and 5.4 GiB RSS on
   CUDA. `RESOURCE_LIMIT` and `NOT_RUN_RESOURCE_LIMIT` counts are zero.
7. The previously fixed false-`INFEASIBLE` claims on `fffff800`, `gfrd-pnc`,
   `greenbea`, and `vtp.base` remain resolved. The root cause was acceptance of
   an approximate Farkas direction without charging stationarity error at
   bound scale. The focused
   `hi-pdlp-farkas-accounts-for-stationarity-error-at-bound-scale` regression
   and the post-fix CPU/CUDA artifacts remain the evidence; no solver change
   was needed in this completion campaign.
8. Baseline first-order solvers still have many solver-reported `OPTIMAL`
   results rejected by the common KKT/reference gates. cuOpt's reported
   reduced costs remain preserved and inconsistent; all 227 parseable cuOpt
   Netlib/energy canonical results use `derived_from_row_dual` for KKT.

## Repository Identity

The campaign began from branch `feature/hipdlp-rewrite` at
`db17f8cfbea6e7ce79b0c76239570c3e32cca635` with the recorded dirty tree: 21
modified and 8 untracked files. During the campaign the checkout advanced
outside this agent to `5a396aaadbcd5a76b7cfa688cbf2588ddd666f1b`
(`Add HiPDLP Phase 7 validation`) and became clean. This agent made no commit.
The final reporting pass leaves only this report and the two Phase 7 plan files
modified in the working tree; they remain uncommitted.

All benchmark results use the fixed CPU and CUDA builds named in the plan.
Run metadata records executable, manifest, and solver-configuration identities,
so the evidence does not rely on the final repository `HEAD` alone.

No Phase 8 implementation work was performed. Legacy remains present, and no
legacy CUDA comparison is claimed because no distinct `CUPDLP_GPU=ON` legacy
build exists.

## Reused CPU Netlib Evidence

No rewrite CPU Netlib suite, shard, or resume run was executed during this
completion campaign, and no new `netlib-*hipdlp-rewrite-cpu*` directory was
created. The canonical report reuses these immutable artifacts:

| Scope | Artifact | Rows | Accepted / statuses |
| --- | --- | ---: | --- |
| smoke | `smoke-recheck-20260903-02/rewrite-cpu` | 13 | 13; 13 `OPTIMAL` |
| Netlib regression | `netlib-regression-recheck-20260903-01/rewrite-cpu` | 103 | 97; 97 `OPTIMAL`, 6 `TIME_LIMIT` |
| full Netlib | `netlib-full-hipdlp-rewrite-cpu-20260903-01` | 114 | 111; 111 `OPTIMAL`, 3 `TIME_LIMIT` |

All 221 reused rewrite CPU `OPTIMAL` rows pass the independent original-LP KKT
gate in their respective scopes. The full-Netlib CPU time limits are `pilot`,
`fit2d`, and `greenbea`.

## Result Counts

### Smoke And Netlib Regression

| Solver | Smoke accepted | Netlib regression accepted | Regression statuses |
| --- | ---: | ---: | --- |
| rewrite CPU | 13/13 | 97/103 | 97 `OPTIMAL`, 6 `TIME_LIMIT` |
| rewrite CUDA | 13/13 | 96/103 | 96 `OPTIMAL`, 7 `TIME_LIMIT` |
| legacy | 5/13 | 28/103 | 81 `OPTIMAL`, 22 process timeouts |
| in-tree cuPDLP | n/a | 40/103 | 95 `OPTIMAL`, 4 `TIME_LIMIT`, 4 `UNKNOWN` |
| cuPDLPx | n/a | 49/103 | 99 `OPTIMAL`, 3 `TIME_LIMIT`, 1 `SOLVE_ERROR` |
| cuOpt 26.08 | n/a | 45/103 | 99 `OPTIMAL`, 3 `TIME_LIMIT`, 1 parser failure |
| HiGHS simplex | n/a | 102/103 | 102 `OPTIMAL`, 1 `TIME_LIMIT` |
| HiGHS IPM | n/a | 103/103 | 103 `OPTIMAL` |

The regression `pilot4` difference remains fully classified. At 300 seconds in
full Netlib both rewrite backends return independently verified `OPTIMAL`; the
CUDA result takes 2,905,720 iterations and 229.69 seconds.

### Full Netlib

| Solver | Accepted | Statuses |
| --- | ---: | --- |
| rewrite CPU | 111/114 | 111 `OPTIMAL`, 3 `TIME_LIMIT` |
| rewrite CUDA | 107/114 | 107 `OPTIMAL`, 7 `TIME_LIMIT` |
| legacy | 30/114 | 92 `OPTIMAL`, 22 process timeouts |
| in-tree cuPDLP | 43/114 | 110 `OPTIMAL`, 4 `UNKNOWN` |
| cuPDLPx | 52/114 | 112 `OPTIMAL`, 1 `TIME_LIMIT`, 1 `SOLVE_ERROR` |
| cuOpt 26.08 | 47/114 | 112 `OPTIMAL`, 1 `TIME_LIMIT`, 1 parser failure |
| HiGHS simplex | 113/114 | 113 `OPTIMAL`, 1 `TIME_LIMIT` |
| HiGHS IPM | 114/114 | 114 `OPTIMAL` |

Rewrite CUDA time limits are `d2q06c`, `fit1d`, `fit2d`, `greenbea`,
`greenbeb`, `pilot`, and `sierra`. Baseline acceptance is lower than reported
`OPTIMAL` counts because every solver is held to the same independent KKT and
trusted-objective gates.

### Energy Regression

| Solver | Accepted | Statuses |
| --- | ---: | --- |
| rewrite CPU | 4/4 | 4 `OPTIMAL` |
| rewrite CUDA | 4/4 | 4 `OPTIMAL` |
| legacy | 0/4 | 3 `OPTIMAL`, 1 process timeout |
| in-tree cuPDLP | 1/4 | 4 `OPTIMAL` |
| cuPDLPx | 2/4 | 4 `OPTIMAL` |
| cuOpt 26.08 | 1/4 | 4 `OPTIMAL` |
| HiGHS simplex | 4/4 | 4 `OPTIMAL` |
| HiGHS IPM | 4/4 | 4 `OPTIMAL` |

### Full Energy

| Solver | Accepted | Statuses |
| --- | ---: | --- |
| rewrite CPU | 6/9 | 6 `OPTIMAL`, 2 `TIME_LIMIT`, 1 process timeout |
| rewrite CUDA | 6/9 | 6 `OPTIMAL`, 2 `TIME_LIMIT`, 1 process timeout |
| legacy | 0/9 | 3 `OPTIMAL`, 6 process timeouts |
| in-tree cuPDLP | 1/9 | 6 `OPTIMAL`, 2 `TIME_LIMIT`, 1 process timeout |
| cuPDLPx | 2/9 | 7 `OPTIMAL`, 1 `TIME_LIMIT`, 1 process timeout |
| cuOpt 26.08 | 2/9 | 7 `OPTIMAL`, 1 `TIME_LIMIT`, 1 process timeout |
| HiGHS simplex | 7/9 | 7 `OPTIMAL`, 1 `TIME_LIMIT`, 1 process timeout |
| HiGHS IPM | 7/9 | 7 `OPTIMAL`, 1 `TIME_LIMIT`, 1 process timeout |

All 72 planned full-energy rows were executed. There are no hidden resource
skips. Process timeouts have blank reported solver status and unavailable
metrics in raw data, as required; canonical rows retain those missing fields
and `timed_out=true`.

## Disagreement Classes

`disagreements.csv` contains 626 rows:

| Classification | Rows |
| --- | ---: |
| `BASELINE_DEFECT` | 501 |
| `TIME_OR_ITERATION_LIMIT` | 108 |
| `REFERENCE_OR_MODEL_ISSUE` | 5 |
| `UNSUPPORTED_STATUS_OR_CERTIFICATE` | 8 |
| resolved historical `REWRITE_DEFECT` | 4 |
| `RESOURCE_LIMIT` | 0 |
| `UNRESOLVED` | 0 |

The four `REWRITE_DEFECT` rows are the resolved historical Farkas cases, not
new failures. The five reference/model rows are the corrected historical
`forplan` reference plus the cuPDLPx/cuOpt `forplan` parser results in regression
and full Netlib. The eight unsupported-status rows are in-tree cuPDLP's four
unverified `UNKNOWN` results represented in both Netlib scopes.

No canonical run reports `INFEASIBLE`, `UNBOUNDED`, or
`UNBOUNDED_OR_INFEASIBLE`, so there is no unverified certificate claim in the
completion matrix.

## Machine-Readable Evidence

The result root is
`/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase7`:

- `summary.csv`: 1,879 unique canonical rows;
- `disagreements.csv`: 626 classified rows, zero `UNRESOLVED`;
- `evidence-audit.csv`: 83 run-level or special evidence directories;
- `generate_canonical_summary.py`: cardinality, unique-key, rewrite-KKT, and
  zero-`UNRESOLVED` assertions;
- `full-energy-resource-check-20260904-01/resource-check.json`: resource
  decision and observed huge-instance pressure.

Raw commands, stdout/stderr, adapter `result.json`, metadata, and parsed rows
remain in their immutable run directories. `forplan` uses the trusted objective
`-664.2189612722`. cuOpt raw reduced costs remain unchanged.

All timings are descriptive single runs. They are not the formal
five-repetition performance gate.

## Remaining Non-Phase-7 Gates

The two historical pre-implementation target-absence RED invocations remain
unverifiable because no artifact was retained. A distinct legacy CUDA baseline
was not configured and is not claimed. Neither item is a Phase 7 correctness
completion blocker.

Formal paired performance gates and the legacy stabilization/removal decision
remain outside Phase 7. Phase 7 now satisfies the correctness prerequisite to
begin Phase 8 work, but legacy must not be deleted until those separate gates
authorize removal.