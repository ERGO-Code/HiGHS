# HiPDLP check scheduling and fused FPE reduction

## Behavior

- Restart trials remain every 40 iterations.
- Certificate screening runs every 200 iterations on both backends. Initialization
  and iteration/time limits still check certificates. A time limit that expires
  during a convergence-only check also forces host certificate verification.
- Termination checks start at iteration 40. Let `r` be the largest of the primal,
  dual and gap residuals divided by their respective convergence tolerances. The
  next interval is 40 for `r <= 10`, 80 for `r <= 100`, and 200 otherwise.
  Certificate boundaries also check convergence. The intervals do not affect
  the restart trial cadence.
- Between termination checks, evaluate only the fixed-point error (FPE). If PID
  requests a restart, obtain fresh convergence metrics and anchor distances
  before updating its step sizes. Pilot observations retain their 40-step cadence.
- The CUDA FPE reduction uses two kernels: one traversal for the two stable norms
  and the interaction dot product, followed by a block-partial reduction. Its
  workspace is allocated before iteration. It retains the existing scalar
  readback and error checking.

This adds one check-input flag and a bounded interval rule. There is no new
public option or alternative solver loop.

## Measurements

Local RTX 4070 Ti SUPER; existing Release CUDA build. Baseline is the shared
library saved before the edits. Each configuration runs four times in one
process; numbers below are medians of the final three runs. Baseline and modified
runs execute serially without tests or compilation in the background.

`25fv47` and `80bau3b` use 2,000 iterations, auto scaling and `1e-12` optimality/KKT
tolerances. They hit the iteration limit, with eight restarts in both builds.
`adlittle` uses `1e-6` tolerances and a 20,000-iteration budget. Both builds reach
verified optimality at iteration 2,120 (fixed) and 2,840 (PID), with respectively
10 and 15 restarts.

| Model | Backend / steps | Iteration phase, ms | Reduction | Whole solve, ms | Full checks |
|---|---|---:|---:|---:|---:|
| 25fv47 | cpu / fixed | 137.57 → 121.36 | 11.8% | 152.46 → 136.39 | 50 → 11 |
| 25fv47 | cpu / pid | 153.57 → 110.09 | 28.3% | 168.66 → 124.78 | 50 → 18 |
| 25fv47 | cuda / fixed | 53.69 → 44.79 | 16.6% | 86.74 → 77.42 | 50 → 11 |
| 25fv47 | cuda / pid | 54.48 → 45.73 | 16.1% | 86.58 → 77.51 | 50 → 18 |
| 80bau3b | cuda / fixed | 57.77 → 48.35 | 16.3% | 111.73 → 103.17 | 50 → 11 |
| 80bau3b | cuda / pid | 58.64 → 49.08 | 16.3% | 112.78 → 103.44 | 50 → 18 |
| adlittle | cuda / fixed | 39.48 → 32.98 | 16.5% | 53.93 → 47.46 | 53 → 17 |
| adlittle | cuda / pid | 53.38 → 45.62 | 14.5% | 68.30 → 60.36 | 71 → 41 |

For the CPU 2,000-step runs, host certificate verifications fall from 51 to 11
(including initialization). GPU host certificate verifications remain at one;
the periodic device screen now runs ten times instead of fifty. FPE evaluation
counts and restart counts agree in these fixed-budget comparisons.

These are small local measurements, not a benchmark-suite speedup claim. Startup
includes preprocessing and backend setup; the first process warmup is excluded.
Phase timers measure host wall time, including waits. Moving FPE work outside the
boundary phase changes that phase's accounting, so the table compares the entire
iteration phase and whole solve instead.

## Validation

Build: `cmake --build build-scaling-cuda --target unit_tests -j 6`.
Tests: `build-scaling-cuda/bin/unit_tests '[hipdlp]'`.

The suite has 344 passing cases and nine baseline failures. All nine failures
were reproduced with the saved pre-change library:

- CPU/CUDA oracle fused-product and dual-product overflow fixtures (four cases).
- The phase5h shared-step ownership source-text test (one case).
- Fixed-trace restart cadence resource assertions (two cases).
- Legacy optimal-KKT and scrs8 characterization (two cases).

New coverage checks CPU/GPU certificate gating, convergence metrics with gating
disabled, the 200-step certificate cadence, a terminal check between certificate
boundaries, unchanged iterates and restart counts with extra observations, empty
and uneven GPU vector sizes, and norms at magnitudes `1e-200` and `1e200`.
CUDA Compute Sanitizer memcheck also passes the new empty/uneven/extreme-vector
FPE test with zero errors.

Existing tests cover numerical error/terminal precedence, CPU/GPU agreement,
scaling, graph replay, and allocation behavior.

## Remaining GPU synchronization

FPE-only trials still synchronize every 40 iterations. Removing that wait needs
device ownership of the Halpern iteration count and restart history, conditional
anchor/current-vector updates, PID state and step sizes, and graph scalar
parameters. Moving only the restart comparison does not remove the host
roundtrip. This change leaves that broader state migration for a separate
implementation and measurement.
