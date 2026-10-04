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

## Overlapped original-KKT verification (2026-10-03)

On the CUDA backend with recorded blocks (one block is the 39 graph steps
between two check boundaries, about 40 iteration times), a candidate due for
original-model KKT verification is still downloaded at its boundary, but its
postsolve and verification run on the host only after the next block has been
enqueued, so the device computes while the host verifies.

- An accepted candidate is returned with its own iteration count, solution and
  work counters (restarts, FPE evaluations, replays and so on are restored to
  their values at that boundary; timings keep counting). The discarded block is
  allowed to finish before teardown.
- A rejected candidate updates the backoff before the next boundary decides.
- Iteration and time limits, certificate boundaries and the scaling pilot keep
  the synchronous verification, so their precedence is unchanged. A failure
  inside the overlapped region first resolves the pending candidate, which keeps
  an acceptable candidate ahead of a later numerical error.
- Solves with a phase observer keep the synchronous verification, so phase
  audits still see the accepted boundary as the terminal one. The CPU backend is
  unchanged.

Measured against `35b83fd9c` (same machine and protocol as below), the solution
checksum, iteration count and verification count are identical on all nine
models:

| Model (CUDA) | Iterations | Loop time, s | Ratio |
|---|---:|---:|---:|
| 25fv47 | 9800 | 0.311 to 0.304 | 0.978 |
| scorpion | 16560 | 0.478 to 0.486 | 1.018 |
| pds-20 | 6640 | 0.568 to 0.561 | 0.988 |
| pypsa-eur-240t5n | 2440 | 0.241 to 0.234 | 0.967 |
| nug20 | 4240 | 0.374 to 0.360 | 0.962 |
| mod2 | 140880 | 7.197 to 7.038 | 0.978 |
| a2864 | 8480 | 11.53 to 11.37 | 0.986 |
| neso-2245 | 6880 | 6.349 to 6.166 | 0.971 |
| neso-2005 | 30000 (budget) | 63.92 to 63.13 | 0.988 |
| Geometric mean | | | 0.982 |

(25fv47 and scorpion verify for under 2 ms in total, so their rows are noise.)
Each verification can hide at most one block, while a host verification takes
two to five blocks on the larger models (pypsa 9.3 ms against a 4.0 ms block,
neso-2245 121 ms against 37 ms, neso-2005 418 ms against 85 ms). About 5-9% of
loop time on those models is still the device waiting for host verification.

## Original-model KKT screen on CUDA (2026-10-03)

The host verifier (`verifyKkt`) differs from the internal stopping test: it
uses largest entries rather than norms, normalizes by the bounds that are within
`sqrt(tolerance)` of being active and by the costs whose reduced costs are that
small, uses a midpoint rule for boxed variables, and adds a complementarity
test. The internal metrics therefore often pass long before it does (mod2: 613
rejected verifications).

The CUDA backend now bounds the verifier's four relative metrics from below on
the device, from the scaled candidate: row activities with the original
row-wise matrix (in 64-entry row segments, so long rows are summed by several
threads) and reduced costs from the scaled transpose as postsolve forms them.
Every quantity the host computes differently (summation order, long-double
unscaling, the scaling pipeline's roundings) carries an error bound, and each
`sqrt(tolerance)` threshold and midpoint test that the bound straddles takes
the side that makes the metric smaller. When a lower bound already exceeds the
tolerance the candidate is rejected without the host verifier and without
advancing the backoff; otherwise the host verifier decides as before, so every
`OPTIMAL` is still certified by the same independent verifier. The screen runs
only where internal stopping passes and the backoff has elapsed. The CPU
backend has no screen.

`hi-pdlp-cuda-kkt-screen-bounds-the-host-verifier` checks ten Netlib models
with equality, inequality, ranged and free rows in both senses, from the
optimum to heavily perturbed candidates, at tolerances `1e-4` and `1e-7`. Of 960
cases, 480 pass and 480 fail the verifier; the screen rejected all 480 failing
ones and none of the passing ones, and its bounds never exceeded the verifier's
metrics (smallest ratio 0.99998).

A first version with one warp per row cost 826 us per screen on mod2, because
the scalar FP64 tail ran on one lane while the warp still paid for 32 on this
FP64-weak GPU. One thread per row or column, plus short row segments for long
rows, brought it to 119 us.

Against the overlapped verification above (same protocol):

| Model (CUDA) | Iterations | Host verifications | Screens (rejected) | Loop time, s | Ratio |
|---|---:|---:|---:|---:|---:|
| 25fv47 | 9800 to 9560 | 6 to 1 | 26 (25) | 0.307 to 0.298 | 0.973 |
| scorpion | 16560 to 16520 | 5 to 1 | 15 (14) | 0.482 to 0.492 | 1.022 |
| pds-20 | 6640 to 6360 | 6 to 1 | 21 (20) | 0.563 to 0.507 | 0.901 |
| pypsa-eur-240t5n | 2440 | 4 to 1 | 7 (6) | 0.235 to 0.220 | 0.936 |
| nug20 | 4240 to 3840 | 7 to 1 | 38 (37) | 0.362 to 0.307 | 0.849 |
| mod2 | 140880 to 140520 | 43 to 1 | 614 (613) | 7.166 to 6.969 | 0.973 |
| a2864 | 8480 to 8440 | 5 to 1 | 13 (12) | 11.42 to 10.89 | 0.953 |
| neso-2245 | 6880 to 6680 | 7 to 1 | 21 (20) | 6.186 to 5.492 | 0.888 |
| neso-2005 | 30000 (budget) | 10 to 1 | 70 (70) | 62.74 to 61.35 | 0.978 |
| Geometric mean | | | | | 0.940 |

Acceptance no longer waits for the backoff, so the CUDA iteration counts now
equal the CPU ones, which verify every boundary. Together with the overlap, loop
time is about 0.92 of `35b83fd9c`.

## Follow-up experiments not adopted (2026-10-03)

Recorded so the same ideas are not re-measured. Machine: RTX A4500 (sm_86) and
i9-14900K, one pinned P-core per stream. Release builds from `git archive` of
this commit plus each variant. All solver and KKT tolerances `1e-4`, presolve
off, PID, auto scaling. Variants ran back to back with alternating order; small
CUDA models are best of two runs, large ones single runs. Single-run noise on
small models is about ±4%.

### Independent check of this change

Against the parent `5225aa077` on 20 models (CUDA) and 18 (CPU), iteration and
FPE counts are identical on all 38 runs; loop-time geometric mean is 0.852 on
CUDA and 0.896 on CPU. The fused FPE is bit-identical to the packet FPE (same
thread mapping, accumulation and merge order), so switching between the two does
not perturb restart decisions.

### Rate-based termination interval: no gain

Replacing the `r <= 10 / 100` rule by an extrapolated horizon (decay rate of
`log10(r)` over the second half of the run, next check at a quarter of the
remaining horizon, still capped by the 200-iteration certificate cadence) cut
full checks two to three times further (bnl2 4878 to 1645, 25fv47 223 to 90),
but loop time moved only to 0.987 (CUDA) and 1.003 (CPU) and some runs stopped
later (blendB +14% iterations). With certificates every 200 iterations a full
check costs only three to four iterations more than an FPE-only trial.

### Where CUDA loop time goes after this change (nsys)

Between the first and last graph launch:

| Model | PDHG steps in graphs | Other GPU work at boundaries | GPU idle |
|---|---:|---:|---:|
| 25fv47 | 78% | 17% | 5% |
| scorpion | 79% | 16% | 5% |
| mod2 | 83% | 15% | 2% |
| pds-20 | 72% | 19% | 9% |
| pypsa-eur-240t5n | 68% | 19% | 14% |
| nug20 | 70% | 17% | 14% |
| a2864 | 94% | 6% | 0.1% |

Boundary GPU work is dominated by `checkPartialsKernel`/`checkPacketKernel` and
the certificate kernels. Idle time is 45-70 us per boundary of synchronization
on small models; the 9-14% cases are host original-KKT verifications (about
5 ms each) during which the device waits.

### CPU original-KKT verification backoff: about 1%, not adopted

The CPU repeated rejected verifications every boundary (mod2 614, bnl2 707,
world 1001), 6-14% of CPU loop time. Two schedules after a rejection:

| Model (CPU) | Base: iterations / verifications | Doubling 1,2,4,8,16: iterations / verifications / time | `g = floor(sqrt(2 L c) / 40)`, `c = 15`: iterations / verifications / time |
|---|---|---|---|
| 25fv47 | 9560 / 26 | 9800 / 6 / 0.99 | 9680 / 15 / 0.95 |
| scorpion | 16520 / 15 | 16560 / 5 / 0.92 | 16560 / 11 / 1.00 |
| pds-20 | 6360 / 21 | 6640 / 6 / 1.09 | 6520 / 13 / 1.03 |
| pypsa-eur-240t5n | 2440 / 7 | 2440 / 4 / 0.99 | 2480 / 5 / 1.00 |
| nug20 | 3840 / 38 | 4240 / 7 / 1.09 | 4000 / 18 / 1.00 |
| mod2 | 140520 / 614 | 140880 / 43 / 0.90 | 140760 / 53 / 0.97 |
| a2864 | 8440 / 13 | 8480 / 5 / 1.01 | 8520 / 11 / 0.99 |
| Geometric mean | | 0.99 | 0.99 |

`L` is iterations since the first rejection and `c` the verification cost in
iterations (measured 10-25 on CPU). The time spent verifying cannot be removed
without delaying acceptance, because the passing candidate arrives inside the
rejection phase. The root cause is that the internal (scaled) stopping metrics
pass long before the original-model KKT check; aligning the two would remove
both the waste and the delay.

### CUDA fast-path check: about 2%, not adopted

Prototype: a plain-FP64 kernel for the stopping metrics (residual squares,
objectives and their absolute sums with warp-shuffle reductions); when they show
the candidate is clearly not optimal (residuals above tolerance by a 1e-6
relative margin, or the gap above it after a `64 eps n |sum|` rounding bound)
the certified kernels are skipped, otherwise they run as before. Certificate
boundaries always used the certified path.

Keeping the FPE and anchor distances bit-identical requires the robust scaled
norms (FP64 divisions per element and per merge); that version gained nothing.
With plain FP64 squares for those norms as well ("lean"):

| Model (CUDA) | Iterations | Loop time, s (off to lean) | Ratio |
|---|---:|---:|---:|
| 25fv47 | 9800 | 0.309 to 0.302 | 0.978 |
| scorpion | 16560 | 0.475 to 0.458 | 0.964 |
| pds-20 | 6640 | 0.574 to 0.551 | 0.961 |
| pypsa-eur-240t5n | 2440 | 0.243 to 0.230 | 0.947 |
| nug20 | 4240 | 0.376 to 0.370 | 0.985 |
| mod2 | 140880 | 7.193 to 6.939 | 0.965 |
| a2864 | 8480 | 11.53 to 11.57 | 1.003 |
| neso-2245 | 6880 | 6.367 to 6.343 | 0.996 |
| neso-2005 | 30000 (budget) | 63.21 to 63.69 | 1.008 |
| Geometric mean | | | 0.978 |

Iterations, restarts and PID updates were unchanged on all nine. On 25fv47 the
check cost per boundary fell from 164 us to 100 us; what remains is mostly the
certificate boundaries (certified check plus screen, about 380 us each), and the
idle synchronization above.
