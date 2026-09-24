# CPU Performance With Matched Initial Conditions

Date: 2026-09-14

Status: **Complete: 36/36 jobs, 216/216 records independently accepted.**
Matching the starting point and initial steps largely removes the iteration
count difference, but does not remove the CPU runtime gap. The matched
rewrite/legacy time geometric mean is **2.2664**, with nine of nine valid
pairs. This diagnostic result does not change the Phase 8 release verdict.

## Question And Scope

This experiment asks whether giving rewrite the old legacy solver's starting
point and initial step sizes removes the slowdown on the previous nine-case
CPU batch: `afiro`, `agg2`, `beaconfd`, `blend`, `boeing1`, `brandy`, `degen2`,
`degen3`, and `greenbeb`.

It is an independent diagnostic comparison, not a new Phase 8 release gate.
Production code, default options, prior benchmarks, and frozen libraries are
unchanged. All four variants use the same original-model acceptance criteria.

## Results

Every variant accepts 54/54 records. There are no solver timeouts, process
failures, missing records, or rejected solutions in this campaign. The run
started at 18:58:43 UTC and finished at 19:17:31 UTC on 2026-09-14. Times below
are the medians of five measurements after one warmup, in **seconds**. All
legacy denominators were remeasured in this campaign; no historical timing
sample is substituted.

| Case | Legacy | Rewrite auto | Rewrite scaling only | Rewrite matched | Matched / legacy |
| --- | ---: | ---: | ---: | ---: | ---: |
| afiro | 0.000188 | 0.000674 | 0.000261 | 0.000508 | 2.704 |
| agg2 | 0.013150 | 0.038686 | 0.038465 | 0.029291 | 2.228 |
| beaconfd | 0.010606 | 0.025453 | 0.017180 | 0.023463 | 2.212 |
| blend | 0.001302 | 0.003821 | 0.001757 | 0.003201 | 2.457 |
| boeing1 | 0.058414 | 0.078759 | 0.123644 | 0.123647 | 2.117 |
| brandy | 0.018556 | 0.119549 | 0.043418 | 0.041800 | 2.253 |
| degen2 | 0.011314 | 0.027795 | 0.027412 | 0.026330 | 2.327 |
| degen3 | 0.075082 | 0.166078 | 0.163857 | 0.160844 | 2.142 |
| greenbeb | 23.832823 | 57.078577 | 55.341168 | 48.310953 | 2.027 |

The iteration counts were identical across all six repeats of each variant:

| Case | Legacy | Rewrite auto | Rewrite scaling only | Rewrite matched |
| --- | ---: | ---: | ---: | ---: |
| afiro | 280 | 400 | 280 | 280 |
| agg2 | 3080 | 3760 | 4200 | 3080 |
| beaconfd | 3400 | 3640 | 2360 | 3400 |
| blend | 840 | 2240 | 880 | 840 |
| boeing1 | 15520 | 9320 | 14880 | 14880 |
| brandy | 8480 | 24320 | 8720 | 8360 |
| degen2 | 2520 | 2640 | 2640 | 2520 |
| degen3 | 3240 | 3320 | 3320 | 3240 |
| greenbeb | 656680 | 773440 | 755120 | 656680 |

All paired comparisons retain all nine cases. Ratios are numerator time
divided by denominator time, so a value below one means the numerator is
faster. The confidence interval is a 10,000-resample paired-instance bootstrap,
not a repeat-level or machine-load confidence interval.

| Numerator / denominator | Valid pairs | Time geometric mean | Bootstrap 95% interval |
| --- | ---: | ---: | --- |
| Rewrite auto / legacy | 9 / 9 | 2.7296 | [2.1091, 3.5769] |
| Rewrite scaling only / legacy | 9 / 9 | 2.0116 | [1.6959, 2.3573] |
| Rewrite matched / legacy | 9 / 9 | 2.2664 | [2.1589, 2.4015] |
| Rewrite matched / rewrite auto | 9 / 9 | 0.8303 | [0.6393, 1.0481] |
| Rewrite matched / rewrite scaling only | 9 / 9 | 1.1267 | [0.9313, 1.3976] |

The point estimate against rewrite-auto corresponds to 17.0% less time on the
geometric-mean basis. Its interval includes one, so this batch does not support
a robust across-instance improvement claim. Matching initialization is also
not uniformly better than using legacy scaling with rewrite's own initial
steps: the matched/scaling-only point estimate is 1.1267, and its interval also
includes one.

### Interpretation

- **Matching works as an experimental control.** All actual starting vectors
  and initial step scalars match. Seven cases now terminate at exactly legacy's
  iteration count; `boeing1` and `brandy` terminate at slightly fewer iterations.
  Equal terminal counts are not proof of equal intermediate trajectories.
- **Initialization is not the whole performance problem.** Despite those
  matched counts, every matched solve takes 2.03-2.70 times legacy's time.
  `degen3` is 3240 iterations on both sides but 0.160844 versus 0.075082 seconds;
  `greenbeb` is 656680 on both sides but 48.310953 versus 23.832823 seconds.
  These comparisons leave substantial execution/check/setup work to explain.
- **Improvements are case-dependent.** `brandy` drops from 0.119549 to 0.041800
  seconds relative to rewrite-auto, largely alongside its reduced iteration
  count. `boeing1` rises from 0.078759 to 0.123647 seconds. Neither case is
  omitted from the geometric mean.
- **Legacy's initial step is not automatically a better step.** On `beaconfd`,
  the scaling-only rewrite takes 2360 iterations and 0.017180 seconds; matching
  legacy's initial steps restores 3400 iterations and takes 0.023463 seconds.
  This control is intended to compare implementations, not select an optimum.

The previous [CPU diagnosis](hipdlp-phase8-cpu-diagnosis.md) identified three
SpMVs per rewrite step versus two in an ordinary legacy step, additional
vector/finite-value passes, and repeated certificate work. The present matched
comparison is consistent with a remaining implementation-cost gap, but does
not measure each primitive's causal share. No CPU-kernel optimization was
performed here.

### Timing Limitations

The small cases show a remaining reversal that must not be explained away:
`afiro` takes 0.000508 seconds matched versus 0.000261 with scaling only, despite
280 iterations in both; `blend` takes 0.003201 versus 0.001757 seconds despite
slightly fewer matched iterations. Their within-job ranges do not overlap:
`afiro` is 0.501-0.524 ms matched versus 0.258-0.278 ms scaling-only, and `blend`
is 3.166-3.211 ms versus 1.740-1.835 ms. This is not one discarded outlier.
The matched arm has additional initialization injection, and the experiment
does not separately time its contribution or distinguish it from differences
in check work and between-job machine state. No exact attribution is claimed.

Only one randomized job campaign was performed, without CPU-frequency control
or repeated interleaved process-level blocks. For example, this campaign's
legacy `afiro` median is 0.187728 ms, compared with 0.080432 ms in the earlier
campaign; the same iteration count does not make those timing denominators
interchangeable. The paired-instance bootstrap does not resolve this source
of uncertainty. The table reports the measured experimental wrappers, not an
estimate with their costs subtracted.

The prior [scaling study](hipdlp-phase8-scaling-mechanism.md) also showed that
completion under the current internal stopping rule is not necessarily the
first original-model-acceptable candidate. This experiment preserves both
solvers' stopping semantics and checks the returned solutions independently;
it does not establish first-hit times under a new common stopping policy.

Full-precision medians, iteration counts, initial steps, and ratios are in
[the result CSV](hipdlp-phase8-matched-initial-cpu.csv).

## Four Variants

| Variant | Solver library | Scaling | Initial point and steps |
| --- | --- | --- | --- |
| legacy | Frozen pre-rewrite CPU | Native Ruiz + Pock-Chambolle | Native |
| rewrite-auto | Frozen rewrite CPU | Curtis-Reid + Ruiz + Pock-Chambolle | Native rewrite |
| rewrite-scaling | Same rewrite CPU | Ruiz + Pock-Chambolle | Native rewrite |
| rewrite-matched | Same rewrite CPU | Ruiz + Pock-Chambolle | Captured legacy primal/dual point and exact tau, sigma, omega, eta |

The scaling-only arm separates the effect of disabling Curtis-Reid from the
additional effect of matching initialization. Scalar step sizes alone would
not represent the same effective updates in different scaling coordinates.
Even the corresponding legacy-scaling matrices are not bitwise identical:
the measured maximum relative coefficient difference ranges from `3.51e-16`
to `1.23e-15`. The sparsity coordinates and dimensions agree. No claim is made
that all later trajectories or solver policies are identical.

## Initial-Condition Verification

Both solvers initialize zero dual variables and project zero primal variables
onto the transformed bounds. Passive probes recorded their actual state after
this projection and before the first convergence check. All nine matched
preflights compare the actual primal and dual vectors, and each of the four
step-size scalars, exactly against legacy.

`greenbeb` has nonzero projected starting values. Native rewrite and legacy
differed in the last floating-point bits of some entries. The matched arm
therefore loads the captured legacy primal vector immediately after the
initial projection; original rewrite code subsequently copies the anchors and
computes the matrix-product caches. The zero dual starting point already
matches. The same one-time point-loading mechanism is used for every matched
case, not only `greenbeb`.

The initialization overlay calls the original rewrite scalar initializer and
then installs the captured `tau`, `sigma`, `omega`, and `eta`. The original
operator-norm estimation and initialization work remain inside solver time.
Subsequent PID weight updates, Halpern restart decisions, convergence checks,
and original-model KKT verification are unchanged.

Each timed matched solve asserts exactly one starting-point load and one
step-size override. The three rewrite arms all load the same primal-projection
forwarding overlay, so its extra per-step dispatch is common to them; the two
unmatched arms perform no point replacement. This experimental dispatch means
these runs are not claimed to measure completely uninstrumented rewrite.

All nine matched preflight solutions passed independent original-model KKT at
`1e-4` and reference-objective relative error at most `1e-3`.

### Safety Limitation

Legacy's larger initial steps bypass rewrite's certified norm upper bound in
this experiment only. Capturing a value from legacy is not a proof that the
step satisfies rewrite's general convergence guarantee. No production safety
guard is removed. Numerical failures, timeouts, or independent verification
failures must remain in the planned denominator.

## Measurement Protocol

- Nine fixed cases and four variants: 36 serialized solver-instance jobs.
- One warmup and five measured fresh solves in each process: 216 records.
- PID, Halpern restart, FP64, presolve off, one thread, solver seed `20260908`.
- Stopping, primal feasibility, dual feasibility, and final KKT tolerances all
  explicitly `1e-4`; independent original-model KKT also `1e-4`.
- A 60-second solver budget for every repeat and a separate 90-second watchdog
  reset at each ordered repeat event. No retries or replacement cases.
- Deterministic shuffled job order frozen before execution, using order seed
  `20260914`. Both legacy and rewrite are remeasured in this campaign.
- Primary time is `steady_clock` around the same `solveLpHiPdlp` wrapper.
  Model reading, raw solution writing, and external verification are outside
  this interval. Separate repeat-wall and whole-process-wall values are kept.
- Fixture reading happens once before timed solves. The prior legacy run used
  to acquire fixtures is not charged to rewrite. The one-time starting-point
  allocation/upload is inside solver time. This is a supplied-initial-condition
  experiment, not a deployable automatic initialization procedure.
- Per-instance comparisons use the five measured solver-time median. The
  geometric mean and 10,000-resample paired bootstrap interval are reported
  only when every one of the nine planned pairs is accepted. Warmup is excluded
  from timing medians but included in correctness and completeness checks.

## Evidence

The independent post-run recheck verifies 162/162 pinned input hashes and
1262/1262 sealed output hashes, re-derives 216/216 accepted records, and checks
all solution identities against their preflight results. The frozen input
plan, probe source/build records, exact initial states and fixtures, raw
solutions, per-repeat verifier results, process records, completion seal, and
full-precision summary are retained under:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260914-matched-initial-cpu
```

The original acceptance libraries and dependencies are those pinned by the
2026-09-08 legacy-qualified all-`1e-4` campaign. Diagnostic preflight runs used
a separate 180-second budget and are not performance samples.

The comparison driver has no scalar-trajectory or phase-profiler output.
Before measurement, all four variants passed a two-solve `afiro` smoke check
with the expected raw solution bytes and per-solve initialization counts.
Acceptance self-tests cover missing records, invalid configuration and
initialization counts, process failures, timeouts, nonfinite metrics, KKT
failures, and reference-objective failures. Summary self-tests cover warmup
exclusion, medians, duplicate/missing/rejected coverage, ratios, and bootstrap.

The independent summary command checks the frozen input
and output hashes, recomputes acceptance from numeric values, compares JSONL
with individual driver/verifier records, and verifies raw solution identities:

```sh
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260914-matched-initial-cpu/summarize.mjs --check
```