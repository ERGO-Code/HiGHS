# HiPDLP Correctness-First Original-Model Stopping

Date: 2026-09-22

Same-day follow-up: [the energy50 GPU reduction repair](#energy50-gpu-reduction-repair)
below explains and fixes the CPU/GPU timeout discrepancy under the v2 contract.
Earlier diagnostic timings remain historical observations, not current results.

## Scope And Outcome

The rewrite's mandatory final verifier and the external solution verifier now
share `hipdlp-original-kkt/2`. The change fixes objective-offset dilution and
cancelling complementarity defects; it does not imitate legacy's stopping
iterations or relax the requested tolerance. The frozen legacy solver and
all previous benchmark artifacts remain unchanged.

This aligns final acceptance, not every internal metric. Internal transformed
L2 checks still trigger candidate validation and feed PID/restart. The optional
400-step early original-KKT trigger stays default-off. Thus the change does
not promise the earliest original-acceptable exit or identical trajectories
under arbitrary equivalent formulations.

## Original-Model Contract

The existing original-model primal feasibility, dual sign, row-activity, and
stationarity checks remain required. `Ax` and `A^T y` are recomputed from the
original LP and supplied solution. The column multiplier used for optimality
is the recomputed reduced cost; its difference from the supplied reduced cost
remains part of stationarity.

Let `p0` and `d0` be primal and dual objective contributions excluding the
additive objective constant. The objective scale is

$$
D = 1 + |p_0| + |d_0|.
$$

The gap is `abs(p0-d0)/D`. For each original column and row, choose the finite
bound indicated by the objective-sense-adjusted multiplier. Accumulate the
absolute multiplier times bound-slack product; a nonzero multiplier toward a
missing bound contributes its absolute magnitude, with the existing sign
violation also checked. The additional acceptance metric is

$$
r_{\mathrm{comp}} =
\frac{\sum_j |z_j(x_j-b_j)| + \sum_i |y_i((Ax)_i-b_i)|}{D}.
$$

The implementation handles free/fixed/one-sided/ranged bounds and min/max
sense explicitly. `OPTIMAL` requires finite metrics and all four normalized
metrics (primal, dual/stationarity, gap, complementarity) within the requested
KKT tolerance. Raw maximum complementarity remains available for diagnosis.

Gap and complementarity normalization divide through a common finite scale
to avoid denominator overflow. Objective sums retain ordinary addition order
for ordinary magnitudes and rescale by powers of two if an intermediate sum
overflows. Finite cancellation can therefore recover a finite objective;
nonfinite products, nonfinite aggregate complementarity, or final values
outside FP64 still fail closed. This is not arbitrary-precision arithmetic
or a guarantee against all summation-rounding effects.

Reported objectives still include the objective constant. External reports
now use the same verified primal objective rather than a second, offset-first
summation. If complete KKT metrics are nonfinite, the external objective is
reported as unavailable rather than as a separately trusted value.

### Important Accuracy Limit

Offset-free normalization is invariant to additive objective constants, but
is **not monotonically stricter than the old normalization**. For example,
minimize `x - 1e8` with `x >= 1e8`: the candidate `x = 1e8 + 1` has old relative
gap `1/2`, but v2 relative gap and complementarity approximately `5e-9`.
This is the declared tolerance relative to the original linear objective,
not an absolute error promise for a nearly cancelled reported objective.

An independently known reference objective, where available, remains a
separate benchmark quality check. No reference objective is passed to the
solver. The external absolute reduced-cost consistency field also remains
a separate diagnostic, not a hidden second definition of `kkt_pass`.
Historical benchmark gates that explicitly required it are not rewritten.

## Reproduced Defects

1. For minimize `x`, `0 <= x <= 1`, candidate `x=0.25` and reduced cost 1,
   adding `1e12` made the old normalized gap pass; adding `1e300` could round
   the computed absolute gap to zero. V2 retains absolute gap 0.25 and relative
   gap 0.2 across positive/negative offsets and both objective senses.
2. With costs `(1, 131072)`, lower bounds zero, and candidate
   `(1, -1/131072)`, primal violation is below 1e-4 and total gap is zero.
   Two opposite complementary defects nevertheless have absolute sum 2.
   V2 rejects this; column/row and upper/ranged variants are covered.
3. An internal-optimal iteration-zero candidate rejected by original KKT
   could previously execute one step even with `iteration_limit=0`. Explicit
   original rejection now preserves iteration/time limits. Verified optimal
   and certificate precedence remain, with zero-budget regressions.
4. Review identified finite objective cancellation after intermediate overflow,
   e.g. fixed variables at 1 with costs `(1e308,1e308,-1e308)` and offset
   `-1e308`. The rescaled accumulators accept the exact finite optimum while
   still rejecting true overflow. A second review identified the external
   report's divergent objective sum; real CLI tests reproduced `0` versus
   the correct `+2`/`-2` and now pass.

## Retained-Solution Recheck

The 34 records previously reported OPTIMAL in the fixed large20 screen were
reverified without rerunning or modifying their solutions. All 34 verifier
jobs completed. Counts below are KKT-only, not the historical full gate:

| Frozen role | OPTIMAL records rechecked | Old KKT passes | V2 KKT passes |
| --- | ---: | ---: | ---: |
| Legacy CPU | 8 | 1 | 1 |
| Rewrite CPU | 6 | 6 | 5 |
| Legacy CUDA | 14 | 2 | 1 |
| Rewrite CUDA | 6 | 6 | 5 |

All three removed KKT passes are the 50-node, 168-hour energy model. Rewrite's
9000-iteration CPU/CUDA solutions have relative complementarity approximately
0.00582621; legacy CUDA's has approximately 0.00994162. Their old primal/dual/
gap checks were individually below 1e-4. Rewrite's five historical full-gate
successes continue to pass v2.

## Bounded Correctness Solves

One diagnostic solve per backend reran only
`pypsa-eur-elec_50_168h`, with the existing PID/auto/FP64 settings, all solver
tolerances 1e-4, 60 seconds, and at most 100000 iterations. These use test
configurations and are **not performance samples**. Each returned solution
was reverified with the final independent CPU verifier after review fixes.

| Backend | Status | Iterations | Relative complementarity | Reference-objective error |
| --- | --- | ---: | ---: | ---: |
| CPU | OPTIMAL | 31200 | 1.21554e-6 | 3.15925e-6 |
| CUDA | TIME_LIMIT | 26240 | 1.58323e-4 | 3.38300e-4 |

The CPU solution passes all v2 metrics and the separate 1e-3 reference gate,
instead of accepting at 9000 iterations with reference error 0.00363577.
The CUDA solution does not pass v2 and is correctly returned as limited;
its reference error alone does not justify OPTIMAL. The solve runs preceded
the objective-overflow/reporting review fixes; final-verifier rechecks
confirm these exact stored solutions still have the shown outcomes. No
post-review timing or new complete-solve campaign is claimed.

## Huge GPU Initialization Diagnosis

A separate phase observer links matching headers and the frozen large20 CUDA
library, not the modified solver. On `pypsa-eur-elec_16_1h_2050`, with 60 seconds
and a one-iteration cap, it reproduces TIME_LIMIT with zero iterations:

| Phase | Diagnostic seconds |
| --- | ---: |
| Model read, outside solve | 18.027 |
| Preprocessing | 16.837 |
| Backend setup, total | 0.354 |
| Operator-norm estimate | 61.240 |
| Wrapper solve, total | 80.248 |

These are one instrumented run, not a phase speed comparison. Norm estimation
is the dominant measured phase. At that frozen checkpoint, CUDA `normKernel` and `dotKernel`
each use one thread to scan their vector; power iteration repeatedly invokes
them, matrix products, and scalar readbacks. This identifies a concrete next
profiling/optimization target, but does not assign the entire 61.240 seconds
to an individual kernel. The larger trex watchdog was not rerun and its phase
remains unknown.

No operator-norm bound, initialization step size, or reduction ordering was
changed here. Next work should preserve the certified safety bound while
investigating parallel deterministic reductions and bounded initialization;
do not recover speed by weakening norm safety or optimality checks.

## Verification And Compatibility

| Final verification | Cases | Assertions | Result |
| --- | ---: | ---: | --- |
| CPU affected union | 153 | 22148 | Pass |
| Actual CUDA affected union | 224 | 56818 | Pass |
| CUDA stopping/replay memcheck | 7 | 1078 | Zero errors, zero leaked bytes |
| External-verifier CLI objective regression | 1 per build, two signed fixtures | Python unittest | CPU/CUDA builds pass |
| Final independent saved-solution verification | 36 solutions | 34 historical plus 2 diagnostic | All verifier processes completed |

New tests exercised RED before the fixes. The old early-KKT test's offset-
dependent acceptance was replaced with assertions for rejection, continuation,
and offset-invariant original verification. No test tolerance was loosened.
The preexisting source-token structure failure remains outside the affected
union; a complete repository/release matrix was not run.

`KktMetrics` and `KktVerificationResult` layouts grew. Rebuild internal binary
consumers with matching headers and libraries. External JSON preserves
`lp-cross-solver-result/1` and existing fields, adding the explicit
`kkt_contract=hipdlp-original-kkt/2`, tolerance, finite/pass component flags,
and absolute/relative complementarity. Do not silently merge old and new
`kkt_pass` values as if their contracts were identical.

All evidence is under:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260922-original-stopping
```

`rechecked-solutions-final` seals the final verifier, source, 36 input solution
hashes, outputs, and summary. Earlier rechecks and diagnostic solves remain
separate; their pre-review CPU verifier/source artifacts are retained under
`pre-review-fix`. Initialization evidence is
`frozen-energy16-cuda-phases.jsonl`. No commit, branch, legacy solver edit,
new scaling default, or rewritten historical performance verdict is part of
this change.

## Energy50 GPU Reduction Repair

The question was why CPU reached an accepted solution at 31200 iterations
while CUDA did not finish within 60 seconds on `pypsa-eur-elec_50_168h`.
Matching the iteration budget showed that CUDA also passes KKT v2 at 31200
iterations; the discrepancy was predominantly execution cost, not evidence
that CUDA needed additional iterations.

### Matched Diagnostic

The pre-change CPU and CUDA phase probes used the same model, auto scaling,
PID, FP64, 31200-step cap, and a 120-second diagnostic time limit. The
post-change CUDA probe completed the same work within a 60-second limit.
All three have 24 restarts and 24 PID updates. The two CUDA runs have the
same final metrics and byte-identical raw solutions, not merely equal
iteration counts. This is instance-specific evidence, not bitwise equality
of every intermediate policy scalar or a universal CPU/CUDA parity claim.

| Diagnostic phase | CPU seconds | CUDA before seconds | CUDA after seconds |
| --- | ---: | ---: | ---: |
| Complete solve | 39.868 | 72.789 | 22.117 |
| Operator norm | 0.013609 | 1.036057 | 0.006662 |
| Fixed-point fallback | 0.000 | 49.948 | 1.807 |
| Restart-reference fixed-point work | 0.011103 | 1.497324 | 0.008160 |
| Host verification | 16.162 | 16.327 | 16.344 |

CUDA records 780 periodic checks and 775 host-oracle checks. The latter
includes iteration zero, so 774 periodic boundaries requested the host
oracle. Both CUDA versions retain these counts. A packet requiring the host
oracle is not used for fixed-point values; the fallback computes vector
differences, two norms, and a dot product. Its prior single-thread CUDA
norm/dot loops repeatedly traversed the full vectors, accounting for the
measured dominant fallback phase. The precise reason for each packet's
conservative fallback was not separately instrumented.

These phase observations are instrumented single runs, not production
medians. Host verification and fallback are nested within boundary time;
do not add boundary totals to these components. Asynchronous step submission
timing is not a GPU kernel execution timer.

### Implementation And Correctness

Only the CUDA backend's `dot`/`norm2` execution and its cuBLAS result pointer
mode change. Nonempty vectors up to `INT_MAX` use FP64 `cublasDdot` and
`cublasDnrm2` on the existing backend stream. Device pointer mode directs
results into the existing preallocated device scalar. A following finite
check preserves the sticky numerical error flag, and the host output is
published only after the unchanged readback/error synchronization succeeds.
Empty and larger-index cases retain the old guarded kernel paths.

There is no new library dependency. CPU kernels, `differenceNorm2`, the
mandatory host oracle, the certified operator-norm guard, step-size formula,
PID/restart rules, graph sequence, and all KKT/reference tolerances remain
unchanged. The snapshot comparison confirms that the only production solver
file differing between before and after is the CUDA backend implementation.

Parallel reduction changes floating-point association in general. Tests
cover small contributions that serial accumulation lost, an analytic
mixed-sign dot product, empty vectors, extreme finite norms, true overflow,
unchanged output on failure, and inherited sticky errors. Repeat equality is
tested only in the current environment; cross-device/version bitwise
determinism is not promised. Explicit cuBLAS status-failure injection was not
added, and internal audit counters do not cover vendor-internal allocations.

### Production Comparison

Separate immutable production builds contain the same final KKT v2 code and
driver source, with `BUILD_TESTING=OFF`, Release, no LTO, and no test hooks.
CPU uses Clang 18; both CUDA roles use GCC 13 / NVCC 13, sm_121. Solver
processes run serially on CPU5 with no clock/power changes. One warmup and
five measured fresh solves are performed per role, in the predeclared order
CPU, CUDA-before, CUDA-after. No trace observer or phase diagnostics run here.

Every role uses a 120-second budget and 31200-step cap so the original CUDA
implementation can complete. The wrapper timer includes setup, norm
initialization, graph construction, iterations/checks, in-solver KKT,
postsolve, and teardown; model reading, solution output, and external
verification are excluded. These are rewrite-before/after comparisons, not
new legacy benchmark denominators.

| Production role | Cold warmup seconds | Five-run median seconds | Measured range seconds | Iterations |
| --- | ---: | ---: | --- | ---: |
| CPU control | 39.619 | 39.552 | 39.489-39.812 | 31200 |
| CUDA before | 72.775 | 73.030 | 72.995-73.104 | 31200 |
| CUDA after | 22.204 | 21.662 | 21.551-21.721 | 31200 |

All 18/18 records, including warmups, pass the separate CPU original-model
KKT v2 verifier, reported reduced-cost consistency, and reference-objective
relative error <=1e-3. The reference error is approximately 3.16e-6. All six
GPU before/after solutions have the same raw SHA256; the CPU solution differs
in floating-point bits but independently meets the same acceptance criteria.

The CUDA after/before median ratio is 0.296620: 70.3% less time, or 3.37x
speedup for this instance. After/CPU is 0.547680, versus 1.846405 before/CPU.
The protocol has only one instance and a fixed process order; it supplies no
across-instance confidence interval or control for unobserved machine drift.
It does not establish the full 20-instance or release performance gate.

A separate production check restores the original 60-second time budget
with a nonbinding 100000-step limit. It returns OPTIMAL at 31200 iterations
in 22.140 seconds and produces the same independently accepted GPU solution.
This check is not pooled into the five-run medians.

### Final Checks And Remaining Work

- Actual CUDA affected union: 242 cases, 57641 assertions, all passed.
- Focused reduction/policy/replay memcheck: 19 cases, 1340 assertions,
   zero errors and zero leaked bytes.
- Production comparison: 18/18 accepted, input identities rechecked against
   1424 frozen hashes; raw records, solutions, and verifier outputs agree.
- Separate 60-second production check: accepted without changing accuracy.
- All touched editor diagnostics and diff whitespace checks are clean.

Host verification still takes approximately 16 seconds in the phase probe.
Its fallback frequency, certificate checks, and `differenceNorm2` remain
separate follow-up work. The new reductions also serve operator-norm
estimation, but the huge energy16/trex models were not rerun here; no measured
speed claim for those initialization paths is made. CPU code is unchanged;
the completed CPU production control is not a new full CPU regression suite.

Evidence is retained under:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260922-energy-gpu
```

The diagnostic JSON/solution files are separate from `production/`, which
contains both source snapshots, three build/loader logs and drivers, the
predeclared plan and digest, per-role sealed outputs, full-precision summary,
and the separate `budget60` result. Revalidate without running any solves:

```sh
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260922-energy-gpu/production/measure.mjs check
```

No commit, branch, legacy modification, tolerance relaxation, or replacement
of historical benchmark artifacts is included in this repair.

## Follow-Up CUDA Check Optimizations

This follow-up completes the reasonable reduction and execution candidates
identified after the energy50 repair. It supersedes the preceding remaining-work
list, not the archived measurements or their interpretation. The original KKT v2
contract, mandatory host oracle, scaling, PID/restart formulas, certified norm
bound, legacy solver, and all acceptance tolerances are unchanged.

### Retained Changes

- Nonempty norm/dot reductions use FP64 cuBLAS on the backend stream. CUDA 12+
   uses `cublasDnrm2_64`/`cublasDdot_64`; older toolkits retain ordinary cuBLAS
   for int-sized vectors and guarded serial fallbacks above that range. The
   installed CUDA 13 configuration is the one built and tested here.
- `differenceNorm2` computes checked parallel differences and reduces them with
   cuBLAS, never using the cancellation-prone squared-norm identity. Its scratch
   is allocated while vectors are created, not during checks. Persistent extra
   device memory is eight bytes times the largest vector created by the backend;
   growing the buffer temporarily retains the old allocation as well.
- Fixed-point fallback queues two norms and one dot, then performs one numeric
   readback and one error-flag readback: two copies, 28 bytes, one wait instead
   of six copies and three waits. The compatibility convergence API similarly
   batches its two norms into two copies, 20 bytes, one wait. It is not the main
   resident-check hot path.
- Raw stopping metrics, FPE metrics, anchor distances, prepared FPE products,
   and certificate screening now have independent validity. The resident packet
   remains 112 bytes with one readback/wait. A raw-metric fallback no longer
   discards independently safe FPE or certificate-screen results. Uncertain
   certificate screens still leave candidates possible for host verification;
   they cannot accept a certificate themselves.
- An FPE fallback reuses existing deltas and `A^T delta_y` only when their
   preparation is valid. Invalid optional products are not reused. Optional
   errors remain deferred until after optimal/certificate/limit precedence.
   Both normal and optional FPE preparation replace two copy-plus-axpy pairs
   with two checked difference kernels, removing two D2D copies and enqueues.
- Primal/dual blend uses one kernel instead of two. The forward SpMV output
   is checked inside its immediate dual projection consumer, before projection
   can hide Inf. Other SpMVs, including the end-of-step transpose product, keep
   their standalone checks. This removes two primitive enqueues per nonempty
   production step, from eight to six. Direct and Graph execution still use
   the same shared algorithm, including dynamic step/blend parameters.

Already fused checked scale/axpy kernels and the parallel resident/certificate
reductions are retained; blind cuBLAS replacement would add finite-check scans
or undermine their interval-safety contract. Empty reductions remain defined
and checked. New difference/blend grid loops use a wide unsigned index with a
compile-time overshoot bound. No near-INT_MAX allocation test was attempted;
wide cuBLAS entry points are not a claim that every backend operation has been
validated on such vectors. Explicit vendor-status fault injection remains a
test gap. Audit counters do not include allocations internal to cuBLAS/cuSPARSE.

The internal backend vtable, FPE evaluator signature, and check-result structure
changed; internal consumers must be rebuilt against matching headers.

### Production A/B

Two separate frozen production builds use GCC 13 / NVCC 13, sm_121, Release
`-O3`, no LTO, `BUILD_TESTING=OFF`, and no test hooks or phase instrumentation.
Both contain identical original-model verification and benchmark-driver code.
Runs are serial on CPU5, with one warmup and five measurements per case/role,
in fixed order: afiro before/after, then energy50 before/after. The timer includes
setup, checks, in-solver verification, and teardown, but excludes model reading,
solution output, and the separate CPU verification process.

The predeclared budget is 120 seconds and 100000 iterations, with FP64, auto
scaling, PID, and all four configured tolerances at 1e-4. Neither final solution
bits nor stopping iterations were required to equal the baseline. Every run
must independently pass original-model KKT v2, reported reduced-cost consistency,
and reference-objective relative error <=1e-3.

| Instance | Before median seconds | After median seconds | After/before | Iterations before/after |
| --- | ---: | ---: | ---: | ---: |
| afiro | 0.017405984 | 0.016233597 | 0.932645 | 400 / 400 |
| pypsa-eur-elec_50_168h | 21.771453159 | 9.372434466 | 0.430492 | 31200 / 31200 |

The warmed time reduction is 6.7% for afiro and 57.0% for energy50 (2.32x
speedup). Energy50 measured ranges are 21.622-22.270 seconds before and
9.287-9.401 seconds after. Cold warmups are 22.216 and 9.532 seconds respectively;
afiro cold warmups are 0.318050 and 0.317093 seconds. Do not mix these cold
observations into the five-run medians.

All 24/24 records, including warmups, pass independent acceptance. Afiro raw
solution bytes are unchanged. Energy50 raw bytes differ, as expected to be
possible when parallel reduction association changes; each role is individually
repeatable in this environment. Energy50 relative reference error is about
3.15925e-6 for both roles. The after solution has relative primal residual
2.29654e-11, dual residual 9.12955e-5, objective gap 8.02852e-5, and
complementarity 1.21554e-6. No tolerance was relaxed to obtain the speedup.

This is a two-case rewrite-before/after result, not a new CPU performance
comparison, legacy denominator, full 20-instance rerun, or release gate.
Fixed role order and unobserved machine drift limit statistical interpretation.
The combined gain is measured; individual contributions were not isolated with
separate production ablations. The much larger energy16/trex cases were not rerun.

### Verification And Evidence

- CPU-only Clang build: 170 affected cases, 22893 assertions passed.
- CUDA affected union: 247 cases, 57722 assertions passed; the subsequent
   reviewed grid-index repair passes 18 focused cases, 1584 assertions.
- Compute Sanitizer after that repair: 18 cases, 1855 assertions, zero errors
   and zero leaked bytes, including new reductions, optional errors, dynamic
   parameters, and production Graph paths.
- Independent review found and prompted the new grid-index overflow bound.
   Numerical input/output, optional error ordering, and raw fallback guards
   retain their focused regressions.
- Production evidence is sealed against 2118 source, driver, build, input,
   verifier, and resolved-library hashes. All 584 delivered `highs` source files
   match the measured after snapshot. Failed RED/build attempts are retained
   separately and are not counted as passing checks.

Evidence root:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260922-cuda-check-optimizations
```

`production/` contains both source snapshots, build logs, matching drivers,
loader identities, the frozen plan/digest, sealed runs, and full-precision
summary. The independent verifier uses its separately frozen CPU production
library. Recheck all input/output identities and acceptance without solving:

```sh
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260922-cuda-check-optimizations/production/measure.mjs check
```

Existing dirty files and historical artifacts are preserved. No commit or
branch was created.