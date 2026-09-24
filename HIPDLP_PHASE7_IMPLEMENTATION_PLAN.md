# HiPDLP Phase 7 Cross-Solver Validation Implementation Plan

> **For agentic workers:** Execute this plan task by task. Keep solver fixes
> separate from harness work, and stop a comparison as soon as a correctness
> failure has a reproducible case.

**Goal:** Establish whether the rewritten HiPDLP CPU and CUDA solvers return
correct, finite, independently verified results across the existing smoke,
Netlib, and energy suites when compared with legacy HiPDLP, in-tree cuPDLP,
external cuPDLPx, NVIDIA cuOpt, and high-accuracy HiGHS references.

**Architecture:** Every solver runs in a separate process and emits one stable
JSON summary per instance. HiPDLP and HiGHS-family adapters call the in-tree KKT
and certificate verifiers directly. The cuOpt adapter writes a complete
HiGHS-compatible raw solution, then invokes one small in-tree verifier so cuOpt
is judged by the same original-LP calculations as HiPDLP.

**Tech Stack:** CMake/Ninja HiGHS builds, Python 3.12, the existing
`lp-benchmarks` Python harness, CUDA 13, cuOpt 26.08, cuPDLPx, CSV/JSON.

**Spec:** `HIPDLP_REWRITE_PLAN.md`, especially Sections 6 and 8 Phase 7.

## Global Constraints

- Correctness comes before timing. Do not tune the solver during a validation
  run.
- Use `presolve=off`, FP64 where selectable, a 60-second smoke limit, and
  primal, dual, gap, and independent KKT tolerances of `1e-4`.
- A solver's own `OPTIMAL` status is not sufficient. Accept it only when the
  original-LP verifier reports all three relative metrics at or below `1e-4`.
- Compare objectives to trusted references with mixed denominator
  `1 + abs(reference)` and a limit of `1e-3`.
- Treat a legacy disagreement as data, not as proof that the rewrite is wrong.
  High-accuracy HiGHS plus the independent verifier arbitrate disagreements.
- Keep the internal legacy selector until Phase 7 is complete. Do not delete or
  rewrite legacy code in this phase.
- Do not implement adaptive step size, Malitsky-Pock, CUDA Graph solving,
  OpenMP, ROCm, mixed precision, or multistart here.
- Do not run standalone hash, dependency, symbol, binary-string, or source-tree
  audits. Existing manifest-integrity and metadata fields in the benchmark
  harness may remain, but they are not Phase 7 work products.
- Do not run five-repetition performance gates in Phase 7. This phase records
  comparable single-run timings only; formal performance gates remain separate.
- Do not commit unless the user explicitly requests it.

## Fixed Inputs

| Component | Path or version |
| --- | --- |
| HiPDLP source | `/home/yanyu/github/hipdlp2` |
| CUDA build | `/home/yanyu/builds/hipdlp2-m7b-cuda-20260826` |
| CPU build | `/home/yanyu/builds/hipdlp2-m7b-cpu-20260826` |
| Benchmark harness | `/home/yanyu/github/lp-benchmarks` |
| Benchmark data | `/home/yanyu/benchmark-data/lp` |
| Results root | `/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase7` |
| cuPDLPx | `/home/yanyu/github/cuPDLPx/build-baseline/cupdlpx` |
| cuOpt source | `/home/yanyu/github/cuopt`, tag `v26.08.00` |
| cuOpt commit | `400863c1a92914bd47f5d68245f9a00e2d0d2df4` |
| cuOpt environment | `/home/yanyu/github/cuopt-26.08-env` |
| cuOpt executable | `/home/yanyu/github/cuopt-26.08-env/bin/cuopt_cli` |

The cuOpt environment stores `PYTHONNOUSERSITE=1`. Always enter it with
`micromamba run -p /home/yanyu/github/cuopt-26.08-env`; do not call the system
Python for cuOpt.

## Common Result Schema

Each adapter must write `result.json` with these top-level fields:

```json
{
  "schema": "lp-cross-solver-result/1",
  "solver": "hipdlp-rewrite-cuda",
  "instance": "/absolute/path/model.mps",
  "status": "OPTIMAL",
  "return_status": "ok",
  "solver_time_sec": 0.0,
  "iterations": 0,
  "primal_objective": 0.0,
  "relative_primal_residual": 0.0,
  "relative_dual_residual": 0.0,
  "relative_objective_gap": 0.0,
  "kkt_input_valid": true,
  "kkt_pass": true,
  "certificate_kind": "none",
  "certificate_pass": false
}
```

Use JSON `null` for an unavailable scalar. Never write zero as a substitute for
an unavailable metric. The benchmark runner must reject `OPTIMAL` when any KKT
field is missing, non-finite, or above `1e-4`.

## Task 1: Freeze the Phase 7 Execution Point

**Files:**

- Modify only if needed to fix a failed readiness check:
  `highs/pdlp/hipdlp/api/options.cc`
- Test: `check/TestPdlpHi.cpp`
- Test: `check/TestHiPdlpSolver.cpp`

- [x] Rebuild only the existing CPU and CUDA validation targets:

```bash
cmake --build /home/yanyu/builds/hipdlp2-m7b-cpu-20260826 \
  --target unit_tests hipdlp_migration_adapter -j 2
cmake --build /home/yanyu/builds/hipdlp2-m7b-cuda-20260826 \
  --target unit_tests hipdlp_migration_adapter -j 2
```

- [x] Run the narrow routing and solver checks:

```bash
/home/yanyu/builds/hipdlp2-m7b-cpu-20260826/bin/unit_tests '[migration]'
/home/yanyu/builds/hipdlp2-m7b-cuda-20260826/bin/unit_tests \
  '[phase7][solver][backend][cuda]'
/home/yanyu/builds/hipdlp2-m7b-cuda-20260826/bin/unit_tests \
  '[phase7][solver][scaling][sctest]'
```

- [x] Record the current commit and `git status --short` in the human report.
  The working tree is currently dirty, so do not describe results using `HEAD`
  alone and do not calculate an extra diff hash.

**Gate:** Default migration selects rewrite, explicit CUDA runs CUDA, explicit
PID reaches the solve loop, and the focused postsolve/scaling test passes.

**Evidence (2026-09-03):** On branch `feature/hipdlp-rewrite` at
`db17f8cfbea6e7ce79b0c76239570c3e32cca635` with the dirty working tree
recorded in the Phase 7 audit, both prescribed builds completed successfully.
The CPU `[migration]` filter passed 193 assertions in 2 cases; CUDA
`[phase7][solver][backend][cuda]` passed 61 assertions in 1 case; and CUDA
`[phase7][solver][scaling][sctest]` passed 10 assertions in 1 case. The rebuilt
migration-adapter SHA-256 values are
`eb3f260371809c95c1b53382a56384202e32ad3575e6620bc21d605e4ddd8300`
for `/home/yanyu/builds/hipdlp2-m7b-cpu-20260826` and
`5f538652e15746c1f64bc43b2b105de8139196013f92b497717125560e62e571`
for `/home/yanyu/builds/hipdlp2-m7b-cuda-20260826`.

## Task 2: Add One Adapter for HiGHS-Family Solvers

**Files:**

- Create: `check/LpSolverValidationAdapter.cpp`
- Modify: `check/CMakeLists.txt`
- Modify: `BUILD.bazel`
- Modify: `check/meson.build`

**Interface:**

```text
lp_solver_validation_adapter
  --instance FILE
  --solver simplex|ipm|pdlp
  --presolve on|off
  --time-limit SECONDS
  --kkt-tolerance TOLERANCE
  --output-dir DIRECTORY
```

- [ ] Write a focused test invocation around `check/instances/afiro.mps` that
  initially fails because the target does not exist.
- [x] Implement the adapter with the public `Highs` API. Read the MPS before
  starting the solver wall timer, set the requested solver and options, call
  `run()`, and pass the returned `HighsSolution` to
  `hipdlp::verification::verifyKkt()`.
- [x] For infeasible or unbounded statuses, expose certificate fields only when
  a ray is available and independently verified. Do not turn a missing ray into
  a successful certificate.
- [x] Write `result.json` using the common schema. Map HiGHS model statuses to
  the same uppercase strings used by `HiPdlpMigrationAdapter.cpp`.
- [x] Register one executable target, not one target per solver.
- [x] Build it in the CPU build and run `afiro` through `simplex`, `ipm`, and
  `pdlp`. Each run must produce parseable JSON; simplex and IPM must pass KKT.

```bash
cmake --build /home/yanyu/builds/hipdlp2-m7b-cpu-20260826 \
  --target lp_solver_validation_adapter -j 2
```

**Gate:** The same executable produces complete, independently verified JSON
for simplex, IPM, and in-tree cuPDLP without changing solver implementation.

**Evidence (2026-09-03):** The CPU target build succeeded and one
`lp_solver_validation_adapter` target is registered in each supported build
file. Fresh one-instance harness runs are preserved in
`task2-afiro-simplex-20260903-01`, `task2-afiro-ipm-20260903-01`, and
`task2-afiro-pdlp-20260903-01` below the Phase 7 results root. All three are
accepted `OPTIMAL` results with independent KKT pass: simplex used 22
iterations with relative metrics `2.27e-16 / 1.04e-16 / 6.11e-17`, IPM used
7 with `2.84e-17 / 3.47e-17 / 0`, and in-tree cuPDLP used 200 with
`3.00e-5 / 2.39e-5 / 8.65e-5`. The common executable SHA-256 is
`8b036dcfd03663e9cdae29eb7adb6e440897894434f3ee20d7edcb405b71f9af`.
The focused `[phase7][certificate]` filter passes 7 assertions in 1 case, and
source inspection confirms that absent rays remain `none`/`false` while only
independently verified rays set `certificate_pass=true`. No retained command
or session artifact proves the historical pre-implementation missing-target
RED invocation, so that checkbox remains open.

## Task 3: Add the External-Solution Verifier

**Files:**

- Create: `check/ExternalLpSolutionVerifier.cpp`
- Modify: `check/CMakeLists.txt`
- Modify: `BUILD.bazel`
- Modify: `check/meson.build`

**Interface:**

```text
external_lp_solution_verifier
  --instance FILE
  --solution HIGHS_RAW_SOLUTION
  --solver NAME
  --status STATUS
  --solver-time SECONDS
  --iterations COUNT
  --tolerance TOLERANCE
  --output FILE
```

- [x] Generate one valid HiGHS raw solution with the existing HiGHS CLI and
  one deliberately perturbed solution for `afiro`.
- [ ] Confirm the new target is absent before implementation.
- [x] Implement the verifier with `Highs::readModel()` and
  `Highs::readSolution()`. Recompute row activities from the original LP before
  calling `verifyKkt()`; do not trust row activities stored by an external
  solver.
- [x] Preserve the supplied solver status as reported data, but set `kkt_pass`
  only from the independent verifier. Emit the common `result.json` schema.
- [x] Verify that the valid solution passes and the perturbed solution fails.

**Gate:** External primal, row-dual, and reduced-cost vectors are checked by the
same C++ verifier used for HiPDLP, with no copied Python KKT implementation.

**Evidence (2026-09-03):** The target builds in the CPU tree. The preserved
artifacts and both verifier logs are under
`task3-external-verifier-20260903-01` in the Phase 7 results root. HiGHS CLI
generated a complete raw AFIRO solution; the first unsupported CLI-option
attempt is also retained as a failed log. The valid solution passes with
relative primal, dual, and gap metrics `2.00e-16 / 5.17e-14 / 2.06e-14`.
Changing only primal `X01` from `80` to `81` raises the relative primal
residual to `2.12e-3`; the supplied `OPTIMAL` status remains reported while
`kkt_pass` becomes false. Both results use `reduced_cost_source=reported`, and
the valid result passes the reported reduced-cost consistency check. No
retained pre-implementation artifact proves the historical absent-target
check, so that checkbox remains open.

## Task 4: Add the Pinned cuOpt Adapter

**Files:**

- Create: `/home/yanyu/github/lp-benchmarks/adapters/cuopt_adapter.py`
- Create: `/home/yanyu/github/lp-benchmarks/tests/test_cuopt_adapter.py`

**Interface:**

```text
cuopt_adapter.py
  --instance FILE
  --time-limit SECONDS
  --tolerance TOLERANCE
  --verifier FILE
  --output-dir DIRECTORY
```

- [x] Write a test using cuOpt's pinned `afiro_original.mps`. Assert that the
  adapter emits a HiGHS raw solution containing 32 columns, 27 rows, all column
  reduced costs, and all row duals.
- [x] Implement the adapter with `Problem.read()` and `SolverSettings` from
  cuOpt 26.08. Set these exact parameters:

```python
settings.set_parameter("method", 1)
settings.set_parameter("pdlp_precision", 1)
settings.set_parameter("presolve", 0)
settings.set_parameter("time_limit", args.time_limit)
settings.set_parameter("absolute_primal_tolerance", args.tolerance)
settings.set_parameter("relative_primal_tolerance", args.tolerance)
settings.set_parameter("absolute_dual_tolerance", args.tolerance)
settings.set_parameter("relative_dual_tolerance", args.tolerance)
settings.set_parameter("absolute_gap_tolerance", args.tolerance)
settings.set_parameter("relative_gap_tolerance", args.tolerance)
settings.set_parameter("log_to_console", False)
```

- [x] Read primal values from `Variable.Value`, reduced costs from
  `Variable.ReducedCost`, and row duals from `Constraint.DualValue`. Compute row
  activities from each constraint's variable/coefficient data rather than from
  its reported slack.
- [x] Use `problem.SolveTime`, `problem.SolutionStats.nb_iterations`, and
  `problem.Status.name`. Write a complete HiGHS raw solution and invoke
  `external_lp_solution_verifier` to produce `result.json`.
- [x] Run the test through the pinned environment:

```bash
/home/yanyu/.local/bin/micromamba run \
  -p /home/yanyu/github/cuopt-26.08-env \
  python -m unittest \
  /home/yanyu/github/lp-benchmarks/tests/test_cuopt_adapter.py
```

**Gate:** cuOpt is forced to FP64 GPU PDLP, its full solution is independently
verified, and AFIRO is accepted with finite metrics at or below `1e-4`.

**Evidence (2026-09-03):** cuOpt is pinned at tag `v26.08.00`, commit
`400863c1a92914bd47f5d68245f9a00e2d0d2df4`. The command printed above fails
before test discovery because `unittest` treats an absolute file path as a
module name; that failure is retained. The equivalent pinned invocation
`micromamba run -p /home/yanyu/github/cuopt-26.08-env python
/home/yanyu/github/lp-benchmarks/tests/test_cuopt_adapter.py` passes 1 test in
1.322 seconds. A persistent rerun is under
`task4-cuopt-afiro-20260903-01` in the Phase 7 results root: AFIRO is
`OPTIMAL` in 460 iterations, with relative metrics
`5.74e-6 / 1.20e-5 / 8.26e-5` and `kkt_pass=true`. Its raw solution contains
both `# Columns 32` and `# Rows 27` sections for primal and dual data.

cuOpt 26.08 reports 32 finite but incorrect all-zero reduced costs on this
fixture; the independent consistency check fails with maximum absolute error
`10`. The raw values are preserved and classified as a baseline defect. KKT
uses `reduced_cost_source=derived_from_row_dual` in the C++ verifier, so this
evidence validates cuOpt's primal and row-dual solution but does not claim its
reported reduced costs are correct.

## Task 5: Add Solver Configurations to the Harness

**Files:**

- Create: `/home/yanyu/github/lp-benchmarks/solvers/hipdlp-rewrite-cpu.json`
- Create: `/home/yanyu/github/lp-benchmarks/solvers/hipdlp-rewrite-cuda.json`
- Create: `/home/yanyu/github/lp-benchmarks/solvers/hipdlp-legacy.json`
- Create: `/home/yanyu/github/lp-benchmarks/solvers/highs-simplex.json`
- Create: `/home/yanyu/github/lp-benchmarks/solvers/highs-ipm.json`
- Create: `/home/yanyu/github/lp-benchmarks/solvers/cupdlp-in-tree.json`
- Create: `/home/yanyu/github/lp-benchmarks/solvers/cuopt-26.08.json`
- Retain: `/home/yanyu/github/lp-benchmarks/solvers/cupdlpx.json`

- [x] Give every configuration the same acceptance limits: uppercase terminal
  status, all three relative KKT limits `1e-4`, and reference objective error
  `1e-3`.
- [x] HiPDLP rewrite configurations must pass `--step-size pid`,
  `--presolve off`, `--scaling auto`, and `--output-dir {output_dir}`.
- [x] The legacy configuration must pass `--implementation legacy` and
  `--backend auto`. Do not label this CPU or CUDA; the build determines the
  legacy backend.
- [x] cuOpt configuration must execute `adapters/cuopt_adapter.py` with an
  environment whose `PATH` begins with
  `/home/yanyu/github/cuopt-26.08-env/bin` and whose
  `PYTHONNOUSERSITE` is `1`.
- [x] Add only the smallest generic runner extension needed to supply the
  external verifier path to the cuOpt adapter. Prefer one repeatable
  `--command-variable NAME=VALUE` option over hard-coding a HiPDLP build path in
  the checked-in solver configuration. Cover placeholder substitution with a
  standard-library `unittest`.
- [x] Run `--dry-run` for `afiro` with all seven new configurations and inspect
  the exact commands. A dry run must not create a result directory.

**Gate:** All solver commands point to real executables, use comparable options,
and produce `result.json` in the run directory.

**Evidence (2026-09-03):** `python3 tests/test_run_benchmark.py` passes 6
tests. Fresh AFIRO dry runs for rewrite CPU, rewrite CUDA, legacy, HiGHS
simplex, HiGHS IPM, in-tree cuPDLP, and cuOpt all exit zero. A structured
configuration check confirms identical `OPTIMAL`, `1e-4` KKT, and `1e-3`
reference-objective gates, and confirms that zero dry-run output directories
were created. The exact expanded commands and unittest output are preserved in
`task5-config-validation-20260903-01` below the Phase 7 results root.

## Task 6: Run the 13-Instance HiPDLP Smoke Gate

Use these exact instances from `check/instances`:

```text
25fv47 adlittle afiro avgas blending chip e226 scrs8 sctest shell stair
standata standgub
```

- [x] Run rewrite CPU once per instance with PID, auto scaling, presolve off,
  and a 60-second limit.
- [x] Run rewrite CUDA once per instance with the same options.
- [x] Run legacy once per instance through the internal legacy selector for
  characterization only.
- [x] Require rewrite CPU and CUDA to report `OPTIMAL`, finite values, complete
  postsolved solutions, and independent KKT pass on all 13 instances.
- [x] Compare CPU and CUDA status and objective. Tolerate floating-point metric
  differences, but investigate any different terminal decision.
- [x] Keep known legacy false-optimal or timeout outcomes as legacy defects;
  do not weaken the rewrite gate to match them.

**Gate:** Rewrite CPU and CUDA are both 13/13 independently accepted. If either
is not, stop the larger Phase 7 matrix, preserve the failing command and JSON,
and fix only the demonstrated solver defect before rerunning this task.

**Evidence (2026-09-03):** The provenance-stable rerun is under
`smoke-recheck-20260903-02` in the Phase 7 results root. Its external manifest
contains the exact 13 repository fixtures and verifies every input SHA-256;
the manifest SHA-256 is
`bfb0f6251d2ef4883d8a16f4e5bb0ebc706a0877fdea4c39adc5498f58096e47`.
Rewrite CPU and CUDA are both accepted `OPTIMAL` on 13/13 with finite complete
KKT metrics. Their status and iteration count agree on every instance; the
maximum paired objective relative difference is `1.97e-12` on `scrs8`.
CPU and CUDA executable SHA-256 values are `eb3f260371809c95c1b53382a56384202e32ad3575e6620bc21d605e4ddd8300`
and `5f538652e15746c1f64bc43b2b105de8139196013f92b497717125560e62e571`.
The before/after provenance snapshots are byte-identical.

Legacy is accepted on 5/13. It reports `OPTIMAL` but fails the independent KKT
gate on `25fv47`, `adlittle`, `avgas`, `standata`, and `standgub`; `scrs8`,
`sctest`, and `stair` hit the process timeout. These are retained as legacy
characterization failures. The earlier `smoke-recheck-20260903-01` attempt is
also retained; the harness correctly rejected the mismatched Netlib manifest
before solving, and it is not counted as evidence.

## Task 7: Run Cross-Solver Correctness Suites

Run one measured execution and no warm-up per solver. Start with
`role=regression`; only proceed to full suites when all unexplained failures are
resolved.

```bash
cd /home/yanyu/github/lp-benchmarks

python3 run_benchmark.py \
  --manifest manifests/netlib-feasible.csv \
  --data-root /home/yanyu/benchmark-data/lp/netlib/feasible \
  --solver-config solvers/SOLVER.json \
  --executable EXECUTABLE \
  --output-dir /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase7/RUN \
  --select role=regression \
  --warmups 0 --repetitions 1 --time-limit 60
```

- [x] Execute the regression selection for rewrite CPU, rewrite CUDA, legacy,
  in-tree cuPDLP, cuPDLPx, cuOpt 26.08, HiGHS simplex, and HiGHS IPM.
- [x] Triage every rejected result before expanding the suite.
- [x] Run the full Netlib feasible manifest with a 300-second per-instance
  limit for the same available solvers.
- [x] Run the energy regression selection with a 300-second limit.
- [x] Run the full energy manifest only after checking available disk, RAM, and
  GPU memory against the huge instances. A resource refusal is not a solver
  correctness failure; record it as `NOT_RUN_RESOURCE_LIMIT`.
- [x] Do not claim a legacy CUDA comparison unless a distinct legacy CUDA build
  was actually configured with `CUPDLP_GPU=ON`. The current rewrite CUDA build
  has cuPDLP CUDA disabled and is not a legacy CUDA baseline.

**Gate:** No supported, executed instance has an unexplained wrong status,
non-finite result, failed certificate, failed KKT gate, or reference-objective
error above `1e-3`.

**Completion evidence (2026-09-04):** The canonical summary contains 824
Netlib regression, 912 full-Netlib, 32 energy-regression, and 72 full-energy
rows, with the requested eight solvers represented in every matrix. Full
Netlib accepts rewrite CPU on 111/114 and rewrite CUDA on 107/114; all rewrite
rejections are 300-second time limits. Energy regression accepts both rewrite
backends on 4/4. Full energy accepts each on 6/9, with the same two solver time
limits and one harness process timeout. Every rewrite `OPTIMAL` has complete
finite postsolve data and passes all three original-LP `1e-4` KKT gates.

The full-energy resource artifact records 3.696 TB available disk,
120,252,540 kB available host memory, 16,777,212 kB free swap, and 83.138 GiB
of 121.689 GiB CUDA-visible unified memory free before execution. All 72 rows
were run; no row is `NOT_RUN_RESOURCE_LIMIT` or `RESOURCE_LIMIT`. Exact commands
and logs remain in the immutable run directories. The legacy evidence uses the
CPU build; no legacy CUDA result is claimed.

## Task 8: Triage Disagreements

Classify each disagreement as exactly one of:

```text
REWRITE_DEFECT
BASELINE_DEFECT
REFERENCE_OR_MODEL_ISSUE
TIME_OR_ITERATION_LIMIT
UNSUPPORTED_STATUS_OR_CERTIFICATE
RESOURCE_LIMIT
UNRESOLVED
```

- [x] Reproduce a disagreement with one direct command outside the harness.
- [x] Compare the original-LP KKT metrics and high-accuracy simplex/IPM status
  before reading iteration traces.
- [x] For rewrite CPU/CUDA disagreement, compare the first failing convergence
  boundary, restart decision, and PID outcome. Do not patch the CUDA backend
  when CPU and CUDA exhibit the same trajectory.
- [x] Add a focused regression test before fixing a rewrite defect, make the
  smallest solver change, rerun that test, then rerun only the affected
  instance across CPU and CUDA.
- [x] Leave baseline defects unchanged and document them.

**Gate:** `UNRESOLVED` is empty for every instance included in the Phase 7
completion claim.

**Regression-scope evidence (2026-09-03):** `pilot4` is the only terminal
CPU/CUDA disagreement in the fresh 103-instance rewrite runs. High-accuracy
simplex and IPM both return independently verified `OPTIMAL` near
`-2581.139258883`. At the common 788,840-iteration budget, CPU and CUDA both
return `ITERATION_LIMIT` with independently KKT-passing iterates; CPU takes
10.1 seconds and CUDA 60.7 seconds. The first restart-event split is at
iteration 295,840: both candidates have nearly identical original-LP KKT
metrics, but `current_error - previous_trial_error` is `-5.30e-14` on CPU and
`+1.32e-13` on CUDA, so only CUDA takes the strict necessary-decay branch.
CPU takes the same controlled PID restart one 40-iteration cadence later.
This is classified `TIME_OR_ITERATION_LIMIT`, not a CUDA backend defect, and no
solver change was made. Raw references, equal-budget results, all 493 sampled
boundaries, and debugger captures are under `pilot4-direct-triage-20260903-01`.

Historical false infeasible claims on `fffff800`, `gfrd-pnc`, `greenbea`, and
`vtp.base` are classified resolved `REWRITE_DEFECT`; the focused certificate
regression was observed red before the bound-scale stationarity fix, passes,
and current runs no longer make those claims. The corrected
`forplan` reference is `-664.2189612722`. cuOpt's all-zero reported reduced
costs and external first-order KKT failures remain documented baseline defects;
no upstream solver source was modified. The completed full Netlib and energy
scopes add no new rewrite defect, and every rejected row has a final
classification.

## Task 9: Produce the Phase 7 Report

**Files:**

- Create: `docs/hipdlp-phase7-cross-solver-report.md`
- Create:
  `/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase7/summary.csv`
- Create:
  `/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase7/disagreements.csv`

- [x] Report one row per solver/instance with status, objective, iterations,
  solver time, the three independent metrics, KKT/certificate verdict, and
  acceptance reason.
- [x] Put findings first: rewrite defects, unresolved disagreements, and tests
  not run. Keep raw command logs in result directories rather than pasting them
  into the report.
- [x] Include separate counts for smoke, Netlib regression/full, and energy
  regression/full. Do not hide resource-skipped instances in a pass rate.
- [x] State plainly that timings are descriptive single runs, not the formal
  five-repetition performance gate.
- [x] Update `HIPDLP_REWRITE_PLAN.md` with a short Phase 7 outcome and remaining
  blockers. Do not append assertion totals, dependency scans, or hash audits.

**Evidence (2026-09-04):** `summary.csv` contains 1,879 unique canonical
scope/solver/instance rows: 39 smoke, 824 Netlib regression, 912 full Netlib,
32 energy regression, and 72 full energy. `disagreements.csv` contains 626
classified rows with zero `UNRESOLVED`. `evidence-audit.csv` covers 83
run-level or special evidence directories. The findings-first human report is
`docs/hipdlp-phase7-cross-solver-report.md`; it records complete scope counts,
resource evidence, descriptive single-run timing boundaries, and the remaining
non-Phase-7 performance/stabilization gates. The rewrite plan carries the same
completion outcome.

## Phase 7 Completion Criteria

Phase 7 is complete only when:

- rewrite CPU and CUDA pass the 13-instance smoke gate 13/13;
- every rewrite `OPTIMAL` result in the executed Netlib and energy scope passes
  the original-LP `1e-4` primal, dual, and objective-gap gates;
- every rewrite infeasible or unbounded claim has a verified certificate;
- no rewrite result is non-finite or success-shaped after an error;
- all status/objective disagreements have a non-`UNRESOLVED` classification;
- cuOpt 26.08, cuPDLPx, in-tree cuPDLP, and HiGHS references are represented in
  the regression comparison, or a concrete unavailable reason is reported;
- `summary.csv`, `disagreements.csv`, raw results, and the human report agree.

Formal performance gates, legacy deletion, and the end of the stabilization
window are not Phase 7 completion criteria.

**Completion record (2026-09-04):** All Phase 7 completion criteria above are
satisfied. Phase 7 clears the correctness prerequisite to begin Phase 8 work;
the separate performance and stabilization gates still control legacy removal.

## Phase 8 Handoff Update (2026-09-08)

Phase 7 remains the historical correctness campaign above. Phase 8 implementation
and new acceptance evidence are recorded separately in
`docs/hipdlp-phase8-gpu-performance-report.md`; they do not overwrite the Phase 7
matrix. True legacy CUDA has now been frozen in an independent pre-rewrite
`CUPDLP_GPU=ON` build outside the source tree.

Actual-option verification found that the Phase 7 migration adapter set final
`kkt_tolerance=1e-4` but retained the independent
`pdlp_optimality_tolerance=1e-7` default. Phase 8 protocol r2 reproduces those
effective rewrite settings. The first new driver's explicit `1e-4` iteration
option failed the larger energy reference-objective preflight and is retained
as rejected configuration evidence, not as a speedup. Final KKT and reference
acceptance thresholds remain `1e-4` and `1e-3`, respectively. Formal paired
performance execution for that historical four-case denominator retained all
24 jobs and 144 records, with 126 accepted and 18 rejected legacy records. CPU
and CUDA legacy gates were NOT_PASS; accepted-subset diagnostic ratios must not
replace the full declared denominator. The current focused CPU build matrix and CUDA regressions pass,
but the performance, full-release and stabilization/removal requirements remain
unsatisfied. See the Phase 8 report for complete before/after timing, measured
communication/allocation costs, and the unresolved >20% legacy regressions.

The later Phase 8 legacy-qualified rerun corrected the four-case selection and
accuracy mismatch without changing this Phase 7 record. Phase 7 supplied only
30 full-Netlib legacy outcomes that were accepted within 60 seconds and no
eligible energy outcomes. Fresh true legacy CPU/CUDA qualification at uniform
`1e-4` froze a 23-instance strict intersection before rewrite execution. The
formal 552-record run leaves CPU at 22/23 complete pairs and CUDA at 20/23 due
to rewrite time limits, so both performance gates remain NOT_PASS. The complete
selection and result chain is documented in the Phase 8 report and
`20260908T205450Z-legacy-qualified-1e4` artifacts.