# HiPDLP Phase 8 GPU Performance Acceptance

Date: 2026-09-08

Later experiment: [2026-09-22 large-instance CPU/CUDA screening](hipdlp-phase8-large20-screen.md)
compares the latest rewrite with true legacy on 20 size-selected models.
It is a separate single-run cohort, not a replacement for this acceptance gate.

Status: **Legacy-qualified acceptance finished; Phase 8 BLOCKED.** The current
denominator is a 23-instance strict intersection selected only after fresh
frozen legacy CPU/CUDA qualification at uniform `1e-4` accuracy. All 92 formal
jobs and 552 planned records ran; 528 records pass and 24 candidate records
reach the 60-second limit. CPU has 22/23 valid pairs and CUDA 20/23, so the
predeclared no-denominator-shrink rule reports formal GM/CI as unavailable and
both backend gates as NOT_PASS. Release coverage, strict process-wide
allocation freedom and legacy removal remain unestablished. This report does
not replace the Phase 7 correctness record.

## Current Legacy-Qualified Rerun

The earlier four-instance performance conclusion documented under
"Historical Four-Case Formal Gate" below is retained as historical evidence,
but it is no longer the current comparison denominator. A new legacy-only
qualification and formal campaign completed under
`20260908T205450Z-legacy-qualified-1e4`.

Selection did not inspect rewrite outcomes. Phase 7 full Netlib contained 30
unique legacy results that were independently accepted and completed within 60
seconds; legacy energy contained none. Therefore the available data could not
honestly provide approximately 50 unique common cases. Fresh qualification ran
all 30 through both independent frozen legacy builds with PID, presolve off,
one thread, FP64, and all solver/KKT tolerances exactly `1e-4`. CPU accepted
29/30 and true legacy CUDA accepted 23/30. Their strict intersection froze 23
instances before any candidate execution:

```text
afiro, agg2, beaconfd, blend, boeing1, bore3d, brandy, degen2,
degen3, e226, fit1d, fit2d, greenbeb, grow7, kb2, pds-10, qap8,
recipe, scsd1, seba, vtp.base, wood1p, woodw
```

The excluded legacy qualification cases are `cycle`, `forplan`, `ganges`,
`scagr7`, `shell`, `stocfor2`, and `truss`. They were not replaced after
candidate results became available. Qualification materialized all 60 planned
records, with 52 accepted and eight rejected. `selected.json` is read-only and
pins both legacy records, processes, verifier outputs, raw solutions and models.

The formal campaign then ran the exact 23-instance set for candidate/legacy on
CPU/CUDA: 92 serialized solver-instance jobs, each with one warmup and five
measurements, giving 552/552 materialized records. Solver stopping, primal,
dual, final KKT and independent KKT tolerances are all `1e-4`; reference
objective relative error remains `1e-3`. This corrects the prior mismatch where
rewrite iteration tolerance was `1e-7`.

Formal acceptance is 528/552. Frozen legacy CPU and CUDA each accept 138/138
records. Candidate CPU accepts 132/138 and candidate CUDA accepts 120/138.
Every rejection is a fixed-set candidate time limit:

| Backend | Valid pairs / planned | Rejected candidate instances | Formal GM / CI | Verdict |
| --- | ---: | --- | --- | --- |
| CPU | 22 / 23 | `fit2d` | N/A | NOT_PASS |
| CUDA | 20 / 23 | `fit1d`, `fit2d`, `greenbeb` | N/A | NOT_PASS |

The predeclared rule reports GM and bootstrap CI only with all 23 pairs. It does
not shrink the denominator after a candidate failure. Consequently the formal
GM is not the geometric mean of the remaining successes; it is unavailable and
the backend verdict is NOT_PASS. For diagnosis only, recomputing a geometric
mean from the sealed accepted-pair ratios gives 2.3907 on CPU (22 pairs) and
2.1477 on CUDA (20 pairs). These values are deliberately absent from the formal
GM field and are not gate substitutes. There are 34 accepted-pair regressions
above 20%.

Notable accepted-pair medians include CPU `fit1d` 41.9455 s versus legacy
0.27165 s (154.41x), CPU `recipe` 0.07927 s versus 0.000615 s (128.82x), CUDA
`recipe` 2.9062 s versus 0.02156 s (134.82x), and CUDA `wood1p` 3.4920 s versus
0.30098 s (11.60x). Rewrite is faster on several cases, including CPU `grow7`
(0.0229x), CPU `bore3d` (0.0561x), CUDA `grow7` (0.0285x), and CUDA `bore3d`
(0.1503x), but those gains do not offset fixed-set failures or widespread
regressions.

The campaign output is fail-closed: records are re-assessed from numeric fields,
JSONL must equal individual record artifacts, raw solutions are hashed before
and after verification, campaign outputs are sealed, and completion is
independently re-derived. `formal-finalization.json` pins the protocol, harness,
selection, 552 records, completion and summary. Full details are in
`qualification-report.md`, `formal/summary.json`, `phase8-result-snippets.md`,
and the per-job logs under that run root.

This rerun confirms the user's concern: the earlier four cases were too small
and included unsuitable legacy denominators. The corrected common-qualified
test is broader and uses equal `1e-4` accuracy, but Phase 8 still fails because
rewrite has incomplete 60-second coverage and is much slower on many accepted
pairs. No threshold was relaxed and no slow case was replaced.

## Evidence Roots

All paths below are under
`/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/`:

- `20260908T000926Z-gpu-resident`: pre-edit rewrite CPU/CUDA binaries, complete
  dependency bundles, source snapshot, development RED/GREEN logs, resource
  audits, and independent pre-rewrite legacy CPU/CUDA builds.
- `20260908T163737Z-acceptance`: current production source and dependency
  freeze, matched compiler builds, fixed-workload measurements, latest Compute
  Sanitizer log, and the rejected first full-solve driver configuration.
- `20260908T163737Z-acceptance-r2`: corrected driver/protocol, twelve accepted
  production correctness runs, and the formal paired campaign. No old result
  directory or executable was overwritten.
- `20260908T205450Z-legacy-qualified-1e4`: legacy-only 30-instance candidate
  pool, 60 fresh qualification jobs, read-only 23-instance common selection,
  and the current 92-job/552-record all-`1e-4` formal campaign.

Within r2, `formal-paired` contains all commands, original records, verifier
outputs and provenance; `strict-results` contains independently revalidated
records without overwriting them. `final-tables/solver-medians.csv`,
`paired-ratios.csv` and `traffic-400-steps.csv` retain full numeric precision.
The campaign ran from 17:25:33 to 19:39:46 UTC on 2026-09-08. All 24 driver
processes exited normally; a solver time-limit status is still rejected.

Source HEAD is `5a396aaadbcd5a76b7cfa688cbf2588ddd666f1b`, plus the preserved
uncommitted changes. The candidate source archive SHA256 is
`3bd3f587820f5a8d591c3937f8e0b003e3ab8cbebc852aafc9fa2c660268d4b4`.
Archive contents and library/driver manifests, not HEAD alone, identify this
experiment. The later changes to acceptance tools are recorded separately;
they do not modify the frozen solver libraries.

## Implementation

- The shared production solve loop obtains backend-neutral raw convergence
  metrics and optional FPE/anchor-distance primitives from `CheckResult`.
  Stopping, restart and PID decisions remain in the algorithm layer.
- CUDA uses deterministic FP64 reductions, preallocated scratch and a pinned
  result packet. A normal packet is 112 bytes including error/protocol fields,
  with one D2H copy and one stream wait. Numeric scalar counts are separate
  from copies and bytes.
- Original bounds, objective offsets, scale recovery, robust norms and
  nonfinite handling follow the serial raw oracle. Arithmetic uncertainty and
  intervals crossing a stopping threshold request the original host oracle.
- Certificate screening uses necessary conditions with outward bounds.
  Suspected certificates still download and run the unchanged independent
  original-model verifier. The tight arithmetic enclosure is currently
  restricted to the checked AArch64/glibc 2.39 environment; other environments
  conservatively retain candidates and/or request the host oracle.
- Optional FPE failures cannot preempt mandatory optimality/certificate/limit
  handling. PID distances are computed before anchors move. Reference FPE
  checks after restart remain special, separately measured work.
- `FixedPointErrorResult::computed` distinguishes calculated numbers from
  uncomputed/failed results; the old three-readback count no longer determines
  numerical success.
- Candidate downloads are batched, pre-reserved and transactional. Terminal
  postsolve and original-model verification remain mandatory. Special paths
  may allocate; this is not a zero-communication or whole-solve-zero-allocation
  claim. The shared Halpern step and precision are unchanged; no solve graph,
  second PDHG algorithm, mixed precision or deferred checks were introduced.

## Configuration Correction

The first acceptance driver incorrectly set `pdlp_optimality_tolerance=1e-4`.
The actual Phase 7 migration path left that independent iteration option at
its default `1e-7`, although its final KKT option was `1e-4`. The rewrite option
translator does not replace the iteration option with `kkt_tolerance`.

This mistake was detected, not used as a speedup: both the frozen rewrite CPU
and candidate CPU stopped the larger energy model at 9,000 iterations and
objective `31756489973.752686`, passing KKT but missing the trusted reference
by `0.0036357736091352116`. Candidate CUDA behaved the same. The Phase 7
simplex and IPM references agree near `31872370698.00065`.

The rejected first protocol and all its raw outputs remain in the first
acceptance root. A real-driver test failed before the fix and passed afterward.
Protocol r2 restores nominal primal/dual and PDLP iteration options to `1e-7`,
retains independent KKT at `1e-4` and reference relative error at `1e-3`, and
reports values from the actual options rather than hard-coded JSON constants.
Frozen legacy still applies its historical KKT override to the stopping
tolerance; equal nominal options do not mean identical legacy/rewrite stopping
logic. This difference is preserved and is not hidden by relabeling timing.

## Builds And Timing

The six role/backend drivers call the same `solveLpHiPdlp` wrapper boundary.
Each owns a fresh solve state and started `HighsTimer`. Solver time includes
wrapper setup, transformation/scaling, backend initialization/upload,
operator-norm work, iterations/checks, final downloads, verification where the
implementation performs it, and local teardown. Model reading and the external
verifier subprocess are outside that timer.

Per-repeat driver wall includes setup and raw-solution writing. Process wall
includes loading, one model read, all six solves, output and exit; it is not
divided into a fictitious per-solve time. First-call and warmed solver times
are reported separately, without claiming flushed caches.

Production builds use `BUILD_TESTING=OFF`, Release `-O3 -DNDEBUG -g`, no test
hooks and no LTO. CPU uses Clang 18.1.3. All compared CUDA builds use GCC 13.3,
NVCC 13.0.88, `sm_121`, FP64 and optimized line information. CUDA true legacy
uses `CUPDLP_GPU=ON` at pre-rewrite revision `e29a6c04`; it is not the rewrite
build with `CUPDLP_GPU=OFF`. Its external source has only a recorded CMake
compatibility patch, not numerical changes. Each driver uses its matching
source headers, library, loader and frozen dependency directory.

The older architecture-75/audit builds and an earlier LTO-enabled rewrite
production build are not formal timing denominators. Full commands and hashes
are in `20260908T163737Z-acceptance/production-build-report.md` and the six
driver `build-provenance.json` files.

## Fixed Workload

The constant-space benchmark reuses `enqueueHalpernPdhgStep`, with sixteen
iteration vectors and no per-step snapshots. It restores the same start for
each 40/400-step workload; restart/PID updates and stopping are absent here,
so these numbers are not full-solve times. Check timing includes candidate
products, raw metrics, FPE and PID distances. Candidate CUDA also includes its
certificate screen, but separate host certificate validation is outside the
check timer. Any numerical host fallback is inside the timer.

Each completed measurement has one warmup and five repetitions. The table
uses the 400-step workload and medians, in milliseconds. Before is the
optimization-free rewrite, not legacy HiPDLP.

| CUDA instance | 400 steps before / after | Check before / after | Four-vector download before / after |
| --- | ---: | ---: | ---: |
| afiro | 21.791 / 22.185 | 0.196 / 0.234 | 0.0397 / 0.0206 |
| e226 | 23.137 / 23.110 | 0.472 / 0.369 | 0.0403 / 0.0212 |
| pypsa_medium | 28.043 / 28.449 | 16.163 / 0.407 | 0.0713 / 0.0383 |
| pypsa-eur-elec_50_168h | Not accepted | Not accepted | Not accepted |

Communication reduction does not guarantee time reduction: afiro's check is
19.4% slower in this microbenchmark; e226's is 21.8% faster and
pypsa_medium's is 97.5% faster. Shared step and SpMV times are largely unchanged.
These process jobs used a fixed reference-first order, not the randomized
formal solve order, and are not substitutes for the legacy gate.

Before CUDA has no Event interface, so the comparison above is synchronized
wall/wall only. Candidate Events are recorded on the backend stream and the
stop Event is completed before reading elapsed time. Candidate check Event
medians are 0.225, 0.361 and 0.398 ms respectively. Event intervals can include
host-induced submission gaps, not only kernel execution. The SpMV field
explicitly includes the existing finite-check kernel.

Additional disaggregated CUDA measurements, milliseconds, before / after:

| Instance | Initialization/upload, single call | Operator norm/steps, single call | Host postsolve/KKT, median | Host certificates, median |
| --- | ---: | ---: | ---: | ---: |
| afiro | 311.726 / 268.617 | 1.539 / 1.577 | 0.0102 / 0.0100 | 0.00474 / 0.00459 |
| e226 | 247.742 / 254.800 | 3.804 / 3.184 | 0.2057 / 0.2047 | 0.04344 / 0.04309 |
| pypsa_medium | 278.315 / 254.578 | 170.788 / 168.677 | 5.392 / 5.374 | 1.824 / 1.816 |
| pypsa-eur-elec_50_168h | 296.279 / 307.803 | 1046.783 / 1035.976 | Unavailable | Unavailable |

Setup observations are single-call data with context/cache effects, not five-
measurement estimates of startup improvement. The incomplete large-energy
microbenchmark still retains its setup observations, but no steady-state
acceptance claim is made from them.

The larger energy CPU/CUDA state check at `1e-9` failed in BOTH frozen rewrite
and candidate CUDA. Maximum scaled error is `7.774563701157255e-9`; all four
before/after CUDA state vectors are exactly equal and finite, with identical
step sizes. The failure threshold was not changed. Those two runs stopped
after warmup, so they do not supply five-measurement medians. All other fourteen
jobs completed; repeated states were identical, and smaller CUDA reference
errors were at most `9.12e-12`. Detailed samples and failed outputs are in
`20260908T163737Z-acceptance/kernel-runs`.

## Communication And Allocation

The production audit classifies a boundary as ordinary only when it continues
without download, host verification or FPE/PID fallback. Cadence packet counters
alone are not this classification. Earlier current-code 400-step audits observed:

| CUDA instance | Ordinary boundaries | D2H bytes / copies / waits | Explicit device / pinned allocations |
| --- | ---: | ---: | ---: |
| afiro | 8 | 896 / 8 / 8 | 0 / 0 |
| e226 | 5 | 560 / 5 / 5 | 0 / 0 |
| pypsa_medium | 5 | 560 / 5 / 5 | 0 / 0 |

This is measured constant-size communication on the real solve path. It is
not a claim that every cadence boundary takes the ordinary path; numerical
oracle, certificate, terminal and restart-reference work remains separate.

The new paired software audit executed all 16 production-wrapper calls after
formal timing ended. Each before/after pair used exactly 400 iterations, PID,
auto scaling, the r2 tolerances and the same instance, and returned
`ITERATION_LIMIT`. Counts cover setup through wrapper-local teardown, not just
the loop. No intercepted failures or in-flight calls remained at scope end.
Before means the frozen preoptimization rewrite, not the legacy denominator.

| CUDA instance | D2H bytes before / after | D2H API copies before / after | Explicit waits before / after | H2D bytes before / after |
| --- | ---: | ---: | ---: | ---: |
| afiro | 11,848 / 4,800 | 303 / 165 | 153 / 85 | 3,479 / 7,594 |
| e226 | 90,104 / 50,448 | 263 / 178 | 133 / 87 | 74,475 / 132,134 |
| pypsa_medium | 5,329,508 / 2,908,308 | 273 / 184 | 138 / 90 | 1,711,337 / 3,761,418 |
| pypsa-eur-elec_50_168h | 35,138,452 / 19,167,756 | 273 / 188 | 138 / 92 | 10,843,682 / 24,178,840 |

D2H bytes fall by 59.5%, 44.0%, 45.4% and 45.5%, respectively. H2D API copies
increase from 14 to 31 in each case because the prepared check/screen owns
additional data. Explicit device allocations increase from 42 to 68 and
pinned allocations from 0 to 10, all including setup. This is a communication
trade-off, not an overall memory-allocation reduction.

| CUDA instance | Host allocation API calls before / after | Host requested bytes before / after | Device requested bytes before / after | Added pinned bytes |
| --- | ---: | ---: | ---: | ---: |
| afiro | 29,269 / 29,632 | 50,129,201 / 51,144,183 | 9,382 / 15,453 | 2,004 |
| e226 | 29,270 / 29,642 | 51,157,709 / 52,113,631 | 124,770 / 195,833 | 16,276 |
| pypsa_medium | 29,297 / 29,699 | 85,047,623 / 80,710,757 | 4,794,028 / 7,487,809 | 968,884 |
| pypsa-eur-elec_50_168h | 29,557 / 30,084 | 276,411,895 / 241,712,621 | 31,088,569 / 48,820,659 | 6,388,692 |

CPU has zero CUDA API traffic; its wrapper host-allocation call counts are
262 / 263 for each smaller instance and 263 / 264 for the larger energy model.
Requested sizes and all free/failure counters are retained in the CSV.

These are diagnostic-only, process-wide interposable glibc/CUDA API counts.
Async success means accepted submission; completion is not inferred from an
enqueue time. Hidden vendor allocators, direct driver calls and implicit waits
are outside the measured surface. Nested/aliased API totals can overlap and
must not be summed into physical-memory usage. The shim perturbs execution and
its times are excluded from every formal timing ratio. See
`traffic-audit/runs/main-serial-v1/summary.json` and `traffic-audit/handoff.md`.

Explicit kernel preload removes first-use runtime allocation bundles from
ordinary operations. The process-wide glibc audit still observes recurring
72-byte allocation/free pairs on the vendor `cuda-EvtHandlr` thread, roughly
every 100 ms. Depending on timing, these appear in step, ordinary or reference
phases. They are not filtered out or waived. The solver thread/backend-owned
resource contract and strict process-wide zero allocation are different claims;
the latter remains blocked. See the original run's `lazy-allocation-handoff.md`
and machine-readable audit logs. Detailed diagnostics are not performance runs.

## Correctness And Gates

Protocol r2 production preflight accepted 12/12 original-LP verified results:
four PID instances on CPU and CUDA, plus afiro/e226 on the public
adaptive-to-fixed mapping for each backend. The larger energy result again
uses 115,280 iterations and passes the trusted-objective gate. These are
preflight calls, not the five-repetition performance aggregate.

The latest CUDA memcheck passed 36 cases / 18,752 assertions with zero errors
and zero leaked bytes, covering packets, certificate screening, downloads,
rounding edge cases and Event resources. Final current-source checks passed:

| Configuration/check | Cases | Assertions | Result |
| --- | ---: | ---: | --- |
| CMake Clang CPU, seven-tag union | 175 | 23,722 | Pass |
| CMake GCC CPU, seven-tag union | 175 | 23,722 | Pass |
| Bazel CPU, seven corresponding file targets | 175 | 23,722 | Pass |
| Meson CPU, seven corresponding file targets | 175 | 23,722 | Pass |
| CMake CUDA, seven-tag union | 236 | 60,159 | Pass |
| Hidden CUDA factory/fallback only | 1 | 28 | Pass |
| Actual GCC production-library oracle probe | 6 | Not Catch-based | Pass |

The production probe verifies exact finite norm/FMA fixtures, inclusive
stopping and unpoisoned mandatory downloads, using the actual frozen GCC
library. Its first attempt had a malformed test CSC start array and failed
before any backend call; the corrected separate binary passes all six cases.
The first CPU refresh stopped at a regeneration guard, not a compilation/test
failure. Routine regeneration of the new development GCC directory preserved
all options; the unchanged matrix script then passed all sixteen invocations.
Both original failure logs remain beside the successful follow-up evidence.

An unfiltered Bazel Backend target still has two recorded source-runfile
failures outside the selected seven-tag scope; this is not a fully green
unfiltered suite. No hidden-device run attempted the GPU-dependent suite.

### Historical Four-Case Formal Gate

The historical r2 experiment retained four instances per backend, seed `20260908`,
one warmup and five measurements, 300-second solve budgets, a separate
330-second per-repetition watchdog, per-instance medians and 10,000 paired
bootstrap resamples. CPU and CUDA are compared separately with their frozen
legacy counterparts. Geometric mean <= 1.05 is the stated gate; the 95% interval
is reported rather than added as a new pass condition. Rejected/missing cases
and >20% regressions cannot be hidden by reducing the denominator.

Historical result: **NOT_PASS on CPU and CUDA**. All 144 records are present;
`incomplete=true` in the summary means accepted-pair coverage is incomplete,
not that jobs were omitted. Both rewrite roles each accept 48/48 records.
Legacy CPU accepts 12/24, legacy CUDA 18/24. All solver libraries, source
snapshots, protocol inputs and dependencies pass the post-run identity check.

| Backend | Accepted legacy pairs / planned | Full four-pair gate | Diagnostic GM on accepted pairs only | Diagnostic 95% CI |
| --- | ---: | --- | ---: | --- |
| CPU | 2 / 4 | NOT_PASS | 3.9903 | [1.1672, 13.6416] |
| CUDA | 3 / 4 | NOT_PASS | 5.2256 | [2.3185, 21.0873] |

The diagnostic geometric means are **not** substitutes for the full four-pair
gate. Failed legacy cases remain in the declared denominator and block a pass.
The accepted subsets already exceed 1.05 substantially.

Rejected records are exactly six repeats each for legacy CPU pypsa_medium
(relative primal residual 0.0110072), legacy CUDA pypsa_medium (relative primal
0.00688722 and dual 0.000139624), and legacy CPU larger energy (`TIME_LIMIT`).
There were no missing results or abnormal driver process exits. The first two
sets were solver-reported `OPTIMAL`, but the independent verifier rejected
them. They are not accepted timing denominators.

Complete warmed solver medians, seconds; ratios use only accepted legacy data:

| Backend / instance | Preoptimization rewrite | Candidate | Frozen legacy | Candidate / before | Candidate / accepted legacy |
| --- | ---: | ---: | ---: | ---: | ---: |
| CPU afiro | 0.001090 | 0.001094 | 0.0000802 | 1.0031 | 13.6416 |
| CPU e226 | 0.047320 | 0.047964 | 0.041092 | 1.0136 | 1.1672 |
| CPU pypsa_medium | 11.5245 | 11.8123 | 81.8352, KKT rejected | 1.0250 | Unavailable |
| CPU larger energy | 168.3157 | 179.5850 | 300.0133, time limit | 1.0670 | Unavailable |
| CUDA afiro | 0.059564 | 0.062225 | 0.021320 | 1.0447 | 2.9186 |
| CUDA e226 | 0.677750 | 0.675715 | 0.291449 | 0.9970 | 2.3185 |
| CUDA pypsa_medium | 20.1207 | 20.5409 | 31.4432, KKT rejected | 1.0209 | Unavailable |
| CUDA larger energy | 246.3246 | 251.9912 | 11.9499 | 1.0230 | 21.0873 |

The four-instance candidate/preoptimization solver-time geometric means are
**1.026876** on CPU (95% CI [1.008340, 1.053360]) and **1.021250** on CUDA
(95% CI [1.003437, 1.038673]). This experiment does not establish a full-solve
speedup. The intervals bootstrap instances, not unobserved clock/load changes.

Candidate first-call and end-to-end boundaries, seconds:

| Backend / instance | First solver call | Warm solver median | Warm repeat wall median | Entire six-call process wall |
| --- | ---: | ---: | ---: | ---: |
| CPU afiro | 0.001254 | 0.001094 | 0.001194 | 0.011102 |
| CPU e226 | 0.063222 | 0.047964 | 0.048195 | 0.309818 |
| CPU pypsa_medium | 11.8340 | 11.8123 | 11.8277 | 71.0837 |
| CPU larger energy | 179.2613 | 179.5850 | 179.7017 | 1078.6766 |
| CUDA afiro | 0.353751 | 0.062225 | 0.062329 | 0.883762 |
| CUDA e226 | 0.895964 | 0.675715 | 0.676052 | 4.496203 |
| CUDA pypsa_medium | 20.7402 | 20.5409 | 20.5586 | 123.9275 |
| CUDA larger energy | 251.5021 | 251.9912 | 252.1207 | 1513.1743 |

Equivalent before/legacy cold and wall fields are in `solver-medians.csv`.
The process wall includes all repetitions and is not a per-solve denominator.

All >20% legacy regressions remain explicitly unresolved: CPU afiro (13.64x),
CUDA afiro (2.92x), CUDA e226 (2.32x), and CUDA larger energy (21.09x). Afiro
rewrite performs 760 iterations versus legacy's 280, partly reflecting the
preserved stopping/scaling differences; that does not explain or waive the
whole ratio. E226 rewrite performs 8,440 iterations versus legacy CPU 15,880
and CUDA 9,320, yet is slower. Larger energy rewrite performs 115,280 versus
legacy CUDA 129,640, yet takes 21.09x as long. Thus iteration count alone is
not a defensible explanation. The unchanged step cost, serial reference/FPE
special paths, extra setup storage and conservative fallback costs need a
separate full-solve optimization task, not relaxed acceptance settings.

This representative experiment cannot close full-release coverage, the
stabilization window, or legacy removal. No failures were waived or slow
instances deleted to manufacture a passing gate.

The numerical re-review found no additional concrete raw-metric/error-path
defect after the two exact overflow/threshold regressions were fixed. It does
not certify every FPE/restart/PID threshold: the existing production fixtures
compare final status/iteration and aggregate restart/PID counts, while the
parallel primitive-to-policy behavior at every strict window edge is not yet
exhaustively covered. The separate production GCC oracle probe now passes,
without conflating it with the Clang audit build. Remaining policy-window
coverage is not hidden by a large total assertion count.

## Environment And Limitations

Host: Linux AArch64, NVIDIA GB10, driver 580.159.03, CUDA 13.0.88. The single
accessible GPU query reported SM clock 214 MHz at that instant, memory clock
and power limit `N/A`. It does not establish sustained clocks or exclusive
load. No sudo, driver/permission/capability change, lock clock, power change,
GPU reset or other user's process termination was performed.

Nsight Systems with CUDA/NVTX only, CPU sampling/context switches disabled and
GPU metrics disabled, did not start its target and produced no timeline. The
agent-owned collection was stopped; the timeline is `unavailable` and was not
retried. NCU counters and CPU sampling were not used. Nonprofiler wall time,
completed CUDA Events and software audits remain the evidence. No privileged
information is required to continue acceptance, and its absence is not evidence
of a performance pass.

GPU jobs are serialized. Host-side editing, small harness tests and read-only
review occurred during the longer campaign; no machine-exclusive or controlled
background-load claim is made. Frozen solver sources, libraries and protocol
remain unchanged during formal execution. Review-driven harness hardening for
driver exit/timeout status and prerequisite-content validation is tested in the
worktree passed its 27 lightweight tests and independently rechecked all 144
saved records and the twelve-case prerequisite without overwriting them.
Strict revalidation preserves every original rejection and additionally checks
driver exit/timeout status. Both original and strict verdicts are NOT_PASS;
the frozen run's missing process-status guard did not change this campaign's
outcome because all driver processes exited normally. The strict checker hash
and zero post-run drift list are saved in `strict-results/summary.json`.

Reproduction entry points (use a new protocol/run directory for new solver
measurements; the paths here are immutable execution records):

```sh
R=/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260908T163737Z-acceptance-r2
python3 -B "$R/frozen/driver/hipdlp_phase8_benchmark.py" run \
  --protocol "$R/protocol.json" --output-dir "$R/formal-paired" \
  --execute --correctness-evidence "$R/production-correctness/summary.json"
python3 -B "$R/strict_revalidate.py"
python3 -B "$R/export_acceptance_tables.py"
node "$R/traffic-audit/run_audits.mjs" main-after-formal main-serial-v1
bash "$R/final-current-build-matrix.sh" --formal-finished
```

The first command was invoked with an outer progress-only wrapper whose exact
code and output are retained in the terminal/command record. No measurement
parameters were changed by that wrapper. Each script refuses output collisions;
do not reuse a completed output directory or regenerate a frozen protocol.

Legacy code and its migration selector remain in-tree. Removal needs satisfied
correctness/performance/build gates and an explicitly defined, evidenced
stabilization window; none is inferred from elapsed calendar time or from the
word deprecated. Independent cuPDLP and `solver="pdlp"` are unchanged.