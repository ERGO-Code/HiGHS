# HiPDLP Two-SpMV Step And CPU Projection Fusion

Date: 2026-09-14

Follow-up: the [2026-09-22 execution optimizations](#cpu-clearing-cuda-fusion-and-production-replay-2026-09-22)
below add CPU clearing removal, dedicated CUDA fusion, and production replay.
Earlier measurements and the historical scope statements remain separate.

Status: **Implemented and verified for the scoped CPU/CUDA regression tests.**
The nine-case CPU production comparison completed all 18 jobs and accepted all
108 records. The geometric mean of optimized/before solver-time ratios is
**0.8420**, with a paired-instance 95% bootstrap interval of **[0.7794, 0.9597]**.
The primary campaign retains a `blend` regression; separate same-core
diagnostics are reported below, not substituted into that result.

The primary campaign is an engineering comparison between two rewrite
versions, not a Phase 8 release acceptance pass. A subsequent fixed-core
rerun against true legacy is reported separately below. Existing baseline
artifacts and the original before/after results are unchanged.

## Implementation

The production Halpern-PDHG step now computes only the products required for
the update: `A * reflected_primal` and `A^T * next_dual`. It no longer computes
`A * next_primal`, and the restart transition likewise no longer refreshes an
unused current-primal product. The `current_ax` pointer was removed from
`HalpernPdhgState`, so callers cannot mistake an outdated value for a live cache.
The production solver keeps its separate iteration-zero product; candidate
products required for convergence and final verification remain unchanged.

The CPU backend implements `fusedProjectedPrimalUpdateAndReflect` and
`fusedProjectedDualUpdateAndReflect`. Each generates the projected candidate
and its reflection in one loop, replacing separate projection, copy, scale,
and axpy calls. Both outputs remain available to convergence, FPE, and restart.
The implementation checks projected, doubled, and reflected values explicitly;
intermediate doubling overflow is rejected even if a hypothetical later
subtraction could cancel it. Invalid handles, aliases, dimensions, and step
sizes are rejected before the loop. Arithmetic errors return failure and do
not promise transactional preservation of all scratch output elements.

The backend-neutral default composes the original operations. In the
2026-09-14 implementation, CUDA used this compatibility implementation;
no dedicated fused CUDA kernel was added in that change. Its CUDA benefit
was limited to the shared removal of the unused SpMV, without a CUDA timing
claim. The 2026-09-22 follow-up below supersedes this implementation limitation.

Explicit trace snapshots and the standalone kernel benchmark compute `Ax`
when they request a snapshot. Thus existing diagnostic outputs remain valid
without making every production iteration compute an unused quantity.

Code references:

- [Shared step](../highs/pdlp/hipdlp/algorithm/halpern_pdhg_step.cc#L29)
- [Backend contract and fallback](../highs/pdlp/hipdlp/backend/backend.hpp#L153)
- [CPU fused operations](../highs/pdlp/hipdlp/backend/cpu/cpu_backend.cc#L476)
- [Trace snapshot preparation](../highs/pdlp/hipdlp/algorithm/fixed_trace.cc#L61)
- [Regression tests](../check/TestHiPdlpBackend.cpp#L33)

### Numerical Boundary

Removing the unused `Ax` also removes error detection that occurred solely
while calculating that unused quantity. The optimized step still validates
the actual update products and projected/reflected outputs; requested snapshots
and convergence candidates still calculate and validate their required
products. This is not a claim that every pathological overflow of the discarded
diagnostic quantity is reported at the old iteration. No stopping tolerance,
PID/restart rule, certificate acceptance, or final original-model KKT check was
relaxed.

The internal backend virtual interface and step-state layout changed. Internal
consumers must be rebuilt with matching headers and libraries; archived drivers
must keep using their archived library bundle.

## Verification

| Validation | Cases | Assertions | Result |
| --- | ---: | ---: | --- |
| Clang CPU, narrow step-fusion tests | 3 | 291 | Pass |
| Clang CPU, affected trace/FPE/restart/PID/solver union | 79 | 4564 | Pass |
| GCC/NVCC CUDA, same affected union | 105 | 19243 | Pass |
| Production CPU benchmark smoke, two instances and two libraries | 4 solves | Not Catch-based | Pass |
| Standalone kernel benchmark, old-compatible and new branches | 2 runs | Not Catch-based | Pass |

The direct step test asserts exactly one forward product, one transpose
product, and zero standalone reflection copy/scale/axpy calls on CPU. Numerical
tests cover bounded and free primal variables, equality and inequality dual
projection, alias/argument rejection, intermediate doubling overflow, and
exact equality with the composed operations for lengths 0, 1, and 257.

The affected regression command was:

```sh
unit_tests '[step-fusion],[phase3c],[phase5g],[phase7h],[phase7m],[phase7n],[phase7o],[phase7p],[solver]'
```

CPU-only builds correctly report the CUDA-only certificate comparison as
unavailable. The CUDA-configured run executes that comparison and the CUDA
fusion-contract fallback. These are scoped regressions, not a claim that the
entire repository test suite or all release build matrices were run.

The standalone CPU kernel benchmark's 400-step `afiro` state file matches
byte-for-byte between the frozen-header compatibility branch and the new
snapshot path, including `Ax`. Production smoke checks for `afiro` and `degen3`
also preserve the frozen result's iteration count and raw solution bytes.

## Primary CPU Measurements

Both sides are native rewrite with `auto` scaling, PID, Halpern restart, FP64,
presolve off, one thread, and all four solver tolerances explicitly `1e-4`.
No initialization override, preload hook, scalar trace, or phase profiling is
loaded. The before library is the previously frozen production rewrite; the
after library is a separate production build. Both use Clang 18, Release
`-O3 -DNDEBUG -g`, no LTO, no testing hooks, and the same frozen runtime
dependencies. The new build has `BUILD_TESTING=OFF`.

There are nine predeclared cases, two roles, one warmup plus five measured
fresh solves per process: 18 serial jobs and 108 records. Both roles are
remeasured in this campaign, with deterministic shuffled job order. Every
repeat has a 60-second solve budget and a separate 90-second watchdog. Failed
or numerically changed runs cannot be dropped from the denominator.

Each returned solution must independently pass original-model KKT at `1e-4`
and reference-objective relative error at most `1e-3`. Additionally, its raw
solution bytes and iteration count must equal the pre-optimization reference.
All 108 records satisfy these conditions: 54/54 before and 54/54 after.

Timing uses the identical `steady_clock` boundary around `solveLpHiPdlp`.
Model input and external verification are outside this timer. Setup, norm
estimation, iteration work, candidate checks, postsolve, and local teardown
remain inside it. Warmups do not enter the medians but must pass verification.

Values below are five-measurement medians in **seconds**. A ratio below one
means faster after optimization; the same iteration count applies to both.

| Case | Before seconds | After seconds | After / before | Iterations |
| --- | ---: | ---: | ---: | ---: |
| afiro | 0.000341 | 0.000290 | 0.8501 | 400 |
| agg2 | 0.038690 | 0.030697 | 0.7934 | 3760 |
| beaconfd | 0.025386 | 0.019675 | 0.7750 | 3640 |
| blend | 0.003812 | 0.005308 | 1.3922 | 2240 |
| boeing1 | 0.079001 | 0.060461 | 0.7653 | 9320 |
| brandy | 0.119471 | 0.092315 | 0.7727 | 24320 |
| degen2 | 0.027723 | 0.022391 | 0.8077 | 2640 |
| degen3 | 0.165411 | 0.128608 | 0.7775 | 3320 |
| greenbeb | 57.139066 | 44.989109 | 0.7874 | 773440 |

All nine pairs are valid. The paired geometric mean is 0.8420299569, or
15.8% less time on that basis. The 10,000-resample paired-instance bootstrap
interval is [0.7793827846, 0.9596624472]. Eight cases improve in the primary
campaign, typically by approximately 19-23% outside `afiro`; `greenbeb` saves
12.150 seconds. Both optimization changes are included together, so these
measurements do not separately attribute time savings to fusion versus removal
of the unused SpMV.

Full-precision numbers: [result CSV](hipdlp-phase8-step-fusion.csv).

## Blend Follow-Up

The primary `blend` result is a 39.2% regression and remains in the table and
geometric mean. The optimized measured times were 6.272, 6.307, 5.308, 3.170,
and 3.107 ms, while before was approximately 3.8 ms. These nonstationary samples
prompted a separate, explicitly conditional diagnostic rather than selective
reruns replacing the original measurements.

Five unpinned interleaved before/after blocks still exhibited two time bands;
even before ranged from approximately 3.8 to 7.2 ms across process medians.
This did not establish a stable performance improvement or regression by itself.

The machine exposes CPUs with different maximum-frequency ranges (CPU 0:
2808 MHz; CPU 5: 3900 MHz). A further diagnostic fixed both roles to the same
allowed CPU using `taskset` on the owned solver subprocess only. No clock,
power, privilege, or other-user process settings were changed. Each CPU used
three interleaved paired blocks, each role again with one warmup and five
measurements.

| Fixed CPU | Before process-median range, ms | After process-median range, ms | Paired after/before range |
| --- | --- | --- | --- |
| CPU 0 | 7.214-7.237 | 6.246-6.258 | 0.8630-0.8664 |
| CPU 5 | 3.810-3.840 | 3.095-3.119 | 0.8105-0.8187 |

All six same-core pairs improve. These results support cross-core scheduling
as a confound in the unpinned comparison; the original samples did not record
CPU identity, so their exact scheduling history is not proven. CPU affinity
also does not control frequency or competing load. The follow-ups preserve
all solution bytes and iteration counts, with independent KKT rechecks, and
do not alter the primary nine-case estimate or its denominator. No universal
`blend` speedup is claimed from its selectively investigated samples.

## Controlled Rerun With Legacy

At the user's request, a new experiment repeated `blend` with three roles:
before rewrite, optimized rewrite, and the true frozen pre-rewrite CPU solver.
Each of CPU 0 and CPU 5 ran three balanced blocks, rotating the three-role
order. Every role/block used one warmup and five measured fresh solves. All
108 records passed independent KKT at `1e-4` and reference-objective error at
most `1e-3`; every raw solution and iteration count matched its role's frozen
reference. All solver subprocesses were pinned with `taskset`, without
frequency, power, or competing-load changes.

The following ranges are the three process medians on each core, in ms:

| Fixed CPU | Before rewrite | Optimized rewrite | Legacy | Optimized / before |
| --- | --- | --- | --- | --- |
| CPU 0 | 7.221-7.234 | 6.235-6.239 | 1.308-1.321 | 0.8619-0.8636 |
| CPU 5 | 3.808-3.834 | 3.070-3.087 | 0.582-0.588 | 0.8040-0.8060 |

The optimization improves all six same-core three-way blocks: approximately
13.6-13.8% less time on CPU 0 and 19.4-19.6% on CPU 5. The apparent primary
39.2% regression does not reproduce with this control. The two unpinned time
bands are consistent with the distinct per-core results, supporting scheduling
as the explanation. Because original samples did not record CPU identity,
that historical assignment remains unproven; the primary regression is not
erased or replaced by this rerun.

### Nine-Case Optimized Versus Legacy

The same nine cases were also rerun as an optimized/legacy comparison,
entirely on CPU 5. Case order was frozen before execution, the two roles were
adjacent for each case with alternating role order, and every job again used
one warmup and five measurements. All 18 jobs and 108 records completed and
passed the independent acceptance checks. Both denominators were measured
in this new campaign, not taken from historical timing data.

The optimized role uses **native rewrite auto scaling and initialization**;
legacy uses its native scaling and initialization. Neither role uses an
initialization override, preload hook, or trajectory logger. This differs from
the earlier supplied-legacy-initial-condition experiment. Both roles keep PID,
FP64, one thread, presolve off, a 60-second budget, and the four solver and
independent KKT tolerances at `1e-4`.

Times are five-measurement medians in **seconds**:

| Case | Legacy seconds | Optimized seconds | Optimized / legacy | Legacy iterations | Optimized iterations |
| --- | ---: | ---: | ---: | ---: | ---: |
| afiro | 0.000080 | 0.000289 | 3.5952 | 280 | 400 |
| agg2 | 0.013045 | 0.030712 | 2.3543 | 3080 | 3760 |
| beaconfd | 0.010615 | 0.019603 | 1.8468 | 3400 | 3640 |
| blend | 0.000585 | 0.003086 | 5.2777 | 840 | 2240 |
| boeing1 | 0.058280 | 0.060388 | 1.0362 | 15520 | 9320 |
| brandy | 0.018434 | 0.092363 | 5.0106 | 8480 | 24320 |
| degen2 | 0.011304 | 0.022346 | 1.9769 | 2520 | 2640 |
| degen3 | 0.074296 | 0.128216 | 1.7257 | 3240 | 3320 |
| greenbeb | 23.869795 | 44.981208 | 1.8844 | 656680 | 773440 |

With 9/9 valid pairs, the optimized/legacy geometric mean is **2.4110780085**,
and the 10,000-resample paired-instance bootstrap 95% interval is
**[1.7448003023, 3.3735514650]**. This is a descriptive fixed-nine CPU result,
not a new release-gate verdict. The interval does not quantify machine-state
uncertainty or generalize to the full 23-case CPU/CUDA acceptance set.

For `blend`, the remaining ratio of 5.2777 decomposes descriptively into
2240/840 = 2.6667 times as many iterations and 1.9792 times the amortized
complete-solve time per iteration. The latter includes setup, checks,
postsolve, and teardown; it is not a pure kernel timer. Thus the optimization
did not make the mathematical trajectory worse, but it did not remove the
native scaling/initialization convergence difference or all execution costs.
`boeing1` is now close to legacy in this campaign (1.0362x), while `degen3`
and `greenbeb` remain approximately 1.73x and 1.88x slower.

The newly pinned legacy `blend` median on CPU 5 is 0.584706 ms, compared with
approximately 1.3 ms in earlier unpinned campaigns. Ratios using those different
denominators must not be mixed. In particular, comparing this 2.4111 GM with
the earlier matched-initial GM would conflate changed initialization, CPU
placement, and measurement campaigns.

Full precision: [fixed-core optimized/legacy CSV](hipdlp-phase8-optimized-vs-legacy-cpu.csv).

## Remaining Scope

- CPU certificate screening and repeated original-model certificate work were
  not optimized in this change.
- A dedicated CUDA project/reflection kernel remains future work; this change
  only exercises its compatible composed implementation and shared two-SpMV step.
- This campaign does not combine the new optimization with captured legacy
  initialization. It compares unchanged native rewrite initial conditions
  before and after the execution-path changes.
- The primary CPU campaign did not fix affinity or frequency. Its bootstrap
  samples instances, not unobserved machine-state changes. The same-core
  follow-up explains a risk without retroactively repairing the primary design.
- There is no new legacy performance gate or permission to remove legacy code.

## Artifacts

The independent study root contains the fixed plan, source snapshots, matching
build and loader records, copied optimized library, raw solutions, per-repeat
independent verifier outputs, process logs, and all timing samples:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260914-step-fusion
```

The primary recheck verifies 314 input hashes and 613 sealed output hashes,
recomputes acceptance rather than trusting stored flags, compares stdout with
the individual records, and confirms solution/iteration parity. The separate
unpinned follow-up contains 60 identity-matching solves and 150 output hashes;
the same-core follow-up contains 72 identity-matching solves and 180 output
hashes. Each follow-up independently re-verifies one identical solution per
process, rather than counting six copies as six independent verifier executions.

```sh
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260914-step-fusion/compare.mjs check
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260914-step-fusion/blend-followup.mjs check
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260914-step-fusion/blend-affinity.mjs check
```

The new three-way `blend` rerun and fixed-CPU-5 nine-case comparison have a
separate frozen plan with 352 input hashes. Each scope has 613 sealed output
hashes, and their combined 216/216 records are individually reverified, with
driver/stdout/solution/verifier bindings checked from raw artifacts:

```sh
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260914-step-fusion-legacy-rerun/rerun.mjs check blend
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260914-step-fusion-legacy-rerun/rerun.mjs check batch
```

## CPU Clearing, CUDA Fusion, And Production Replay (2026-09-22)

The three user-requested execution changes are implemented and independently
validated. Scaling, the certified norm guard, PID, restart, FP64 arithmetic,
certificate acceptance, and original-model KKT tolerances are unchanged.
No legacy replacement or full Phase 8 performance-gate claim is made.

### Implementation

- CPU SpMV no longer clears its output before writing every row. Empty rows
  still explicitly receive positive zero, and accumulation order is unchanged.
- CUDA overrides primal and dual projection/reflection with dedicated fused
  kernels. Each operation uses one primitive enqueue and no D2D copy, replacing
  four enqueues and one copy. The original projection expressions and separate
  rounded doubling/reflection are retained, including intermediate overflow
  detection. Standalone projection operations remain available as an oracle.
- Production solving records steps 2 through 40 once and replays the same
  39-step buffer. Step 1 remains direct for post-restart FPE; short iteration
  tails remain direct. The algorithm, not the backend, generates each step
  through the existing shared Halpern-PDHG function.
- A preallocated 117-double vector carries the current primal/dual steps and
  per-step blend weights. Every replay uploads 936 bytes on the backend stream;
  there are no new iteration-time device/pinned allocations, D2H copies, or
  explicit synchronization calls in the replay step. This is not zero-upload
  execution or a process-wide vendor-allocation guarantee.
- The internal `SolveContext::use_command_buffers` defaults to true. CPU
  retains direct execution; setting this internal control false preserves a
  direct CUDA comparison path. It is not a new public option.
- Iteration limits remain exact. Time is checked before submission and after
  completed blocks; an already submitted 39-step block cannot be interrupted.
  Positive-budget tests inject expiry at completed steps rather than relying
  on a nondeterministic wall-clock race.
- Exceptions during recording cancel capture before resource teardown.
  A real injected host allocation failure reproduced a stuck-capture error
  before this fix and successfully recovered on the next solve afterward.

The backend vtable, shared-step signature, `SolveContext`, and
`SolveDiagnostics` changed. Rebuild internal consumers with matching headers
and libraries; archived drivers must keep their archived bundles.

### CPU Microbenchmark

The same benchmark links the archived before/after SpMV objects, using Clang
18 Release, CPU5, three alternating paired blocks, one warmup plus five
measurements per shape. All 216 samples, including warmups, have identical
output checksums. Six synthetic CSR shapes use 257, 8192, or 65536 rows and
one or eight entries per nonempty row; every seventh row is empty.

| Rows | Entries per nonempty row | After / before |
| ---: | ---: | ---: |
| 257 | 1 | 0.956674 |
| 257 | 8 | 0.984237 |
| 8192 | 1 | 0.952980 |
| 8192 | 8 | 0.977344 |
| 65536 | 1 | 0.953689 |
| 65536 | 8 | 0.980149 |

This is approximately 1.6-4.7% less SpMV time on these synthetic shapes,
not a measured whole-solver CPU improvement. Some individual blocks vary;
no frequency/power control or general confidence interval is claimed.

### GPU Whole-Solve Comparison

The predeclared cohort is `afiro`, `e226`, and `stair`. Each stage has a
separately frozen production library and matching driver, built with GCC 13,
NVCC 13, architecture 121, Release, no LTO, and no test hooks. All use native
auto scaling, PID, Halpern restart, FP64, one thread, presolve off, and all
four solver tolerances plus independent KKT tolerance at 1e-4.

Three balanced blocks rotate the role order. Every job is pinned to CPU5,
with one warmup and five measurements, a 10-second per-solve budget, and a
five-minute campaign budget. All 27 jobs and 162 planned records completed.
Every returned solution independently passes original-model KKT and reduced
cost consistency; raw solution hashes and iteration counts are identical
across all three stages. No failed case was removed.

Times below are milliseconds: the median of three process medians, each
process median using its five measured solves. Timing encloses the same
solver wrapper, including setup, graph construction, uploads, checks,
postsolve, in-solver KKT, and teardown. External verification and model input
are outside the timer. This comparison does not add a trusted-reference
objective gate or compare against legacy.

| Case | Iterations | Before | Fusion only | Fusion + replay | Fusion / before | Replay / fusion | Final / before |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| afiro | 400 | 27.925 | 23.542 | 16.331 | 0.8431 | 0.6937 | 0.5848 |
| e226 | 5720 | 398.500 | 315.878 | 189.819 | 0.7927 | 0.6009 | 0.4763 |
| stair | 5560 | 436.389 | 356.694 | 231.537 | 0.8174 | 0.6491 | 0.5306 |

Fusion alone reduces time by 15.7-20.7%; replay then reduces fusion-stage
time by another 30.6-39.9%. Combined time reduction is 41.5-52.4% on these
three cases. These results do not establish performance on large models,
other GPUs, or the full 23-case release cohort. Small/coarse sample coverage
and uncontrolled clocks/load remain limitations.

### Final Validation And Evidence

| Validation | Cases | Assertions | Result |
| --- | ---: | ---: | --- |
| CPU affected regression union | 148 | 21587 | Pass |
| Actual CUDA affected regression union | 220 | 56407 | Pass |
| CUDA fusion/replay Compute Sanitizer memcheck | 10 | 1255 | Zero errors; zero leaked bytes |
| GPU production paired records | 162 records | Independent KKT and exact solution identity | Pass |

The broader initial union exposed four preexisting failures in
`hi-pdlp-phase5h-shared-step-ownership-boundary`: token bans conflict with the
resident-check implementation already in HEAD. The final union explicitly
excludes this unchanged structure test; it is not reported as passing.
CPU-only runs also report the existing unavailable CUDA-only comparison.

Read-only review found one important capture-exception issue, reproduced
and fixed with RED/GREEN evidence. Follow-up review found no remaining
Critical/Important issue in that fix. All touched source files have no
editor diagnostics. Full repository/release matrices were not run.

Artifacts are retained under:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260922-execution
```

`cpu-measurements` and `gpu-measurements-final` contain frozen plans, hashes,
raw per-process/per-repeat records, and full-precision summaries. Final GPU
libraries are `production-before`, `production-fused`, and
`production-recorded-final`. Matching source archives and build logs are
retained. The earlier `gpu-measurements` run was interrupted and is incomplete;
it is not pooled with the final campaign. Stale-timestamp build attempts are
also preserved separately and were not used as measurement inputs.