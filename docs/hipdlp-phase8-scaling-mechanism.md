# Scaling Mechanism And An Auto-Selection Strategy

Date: 2026-09-09

## Scope And Main Result

This CPU diagnostic follow-up studies `blend`, `brandy`, and `boeing1`.
Every comparison below uses the **same frozen rewrite library**.
`legacy` denotes rewrite's scaling mode, not the old HiPDLP solver.

The extra Curtis-Reid stage changes a diagonal preconditioner, not just a
single step size. Controlled runs show that the unfavorable direction for
`blend/brandy` persists when primal weight and restart times are held fixed.
However, the earlier conclusion that it helps `boeing1` needs qualification:
it reduces completion iterations under the current internal stopping rule,
but the less-scaled run already has an independently acceptable original-model
solution at iteration 4800, well before its reported finish at 14880.

No production code, default scaling policy, or frozen acceptance result was
changed. The auto-selection recommendations below are proposals, not a
validated implementation.

## What The Extra Stage Changes

The compared pipelines are `Ruiz -> Pock-Chambolle` and
`Curtis-Reid -> Ruiz -> Pock-Chambolle`. All other solver settings remain
fixed. Both initialize `eta=0.998` and the same primal weight for each model.

The implementation's Curtis-Reid stage alternates row and column averages of
coefficient logarithms. Its local target is balancing matrix coefficients,
not minimizing original-model KKT error or predicting the final active set.
It sees `A`, not the objective, right-hand sides, bounds, or optimum.
See [computeCurtisReidScaling](../highs/pdlp/hipdlp/model/scaling.cc#L232).

With positive diagonal scales and the convention

$$
\widetilde A=D_r A D_c,\qquad x=D_c\widetilde x,\qquad
y=D_r\widetilde y,
$$

the preconditioner entries in unscaled, transformed coordinates are

$$
T_j=\tau (D_c)_{jj}^{2},\qquad
\Sigma_i=\sigma (D_r)_{ii}^{2},\qquad
\tau=\eta/\omega,\quad \sigma=\eta\omega.
$$

These formulas describe the directional update weights; actual movement also
depends on residuals, projections, coupling, and Halpern averaging. A large
`T_j` or `Sigma_i` is not by itself a guarantee of faster convergence.

The ratios below compare extra scaling with no extra scaling at the identical
initial `eta/omega`. Only original variable columns are included; `boeing1`
has 384 original columns and 473 transformed columns including slacks.

| Case | Original-column T ratio: min / median / max | Row Sigma ratio: min / median / max |
| --- | --- | --- |
| blend | 0.0145 / 0.986 / 44.9 | 0.00723 / 0.772 / 34.1 |
| brandy | 0.00737 / 0.543 / 100 | 0.0100 / 1.00 / 75.8 |
| boeing1 | 0.199 / 53.2 / 471 | 0.00137 / 0.186 / 2.51 |

Thus this is not uniformly stronger preconditioning. Some directions become
much faster while others become much slower. Global PID weight can rebalance
primal versus dual motion, but cannot undo an arbitrary diagonal distortion.

## Where The Primal Tail Lives

The following snapshots locate the largest contributions to the unscaled,
transformed primal residual's squared norm. They are not original-model
relative KKT residuals. Row numbers in the accompanying CSV are zero-based.

| Case / iteration | Selected slow direction | Sigma without extra scaling | Sigma with extra scaling | Absolute row violation: without / with |
| --- | --- | ---: | ---: | --- |
| blend / 800 | Row 13, 46.7% of the extra-scaled primal squared residual | 6.499e-3 | 9.006e-4 | 2.904e-4 / 3.012e-1 |
| brandy / 8480 | Original row 144, 82.9% of the extra-scaled primal squared residual | 1.678e-4 | 8.779e-9 | 5.429e-3 / 8.057 |
| boeing1 / 9320 | Original row 92, 99.9986% of the less-scaled primal squared residual | 8.143e-6 | 1.349e-3 | 1.800e-1 / 4.223e-5 |

In `blend`, rows 13 and 11 together account for 84.8% of the extra-scaled
primal squared residual at iteration 800. Both receive smaller effective
dual updates despite the extra-scaled run's larger global `sigma`.

In `brandy`, row scaling and later weight evolution combine especially badly
for the dominant direction. At iteration 8480, the extra-scaled row scale is
0.0067133 instead of 0.0393499, and `omega` is 0.00019519 instead of 0.108563.
The resulting `Sigma` is about 19,109 times smaller. The first PID update was
already different: weight changes from the common 0.0101613 to 0.0076833 with
extra scaling, versus 0.0192528 without it.

In `boeing1`, the later global balance goes the other way. At iteration 9320,
extra scaling has `omega=1.00056`, versus 0.00332046 without it. For the
dominant slow row, the effective dual update is about 166 times larger and
the transformed primal violation is much smaller.

These snapshots explain the update geometry and identify bottlenecks. They
do not prove that changing only the selected row's scale would remove the
delay. The effects of the coupled system are tested more directly below.

## Controlled Experiments

Each cell is **without extra scaling / with extra scaling**. The same
rewrite, FP64, initialization bound, and four solver tolerances of `1e-4`
are retained. Completed solves also pass the frozen independent KKT verifier
and reference-objective relative error at most `1e-3`.

Fixed weight means the common initial weight is held constant. Fixed restart
means restart at iteration 40 and positive multiples of 160, unless the solver
terminates at that boundary first. This is an arbitrary common diagnostic schedule, not a recommended
production policy. The limit is 100000 iterations; capped cases remain failed
to complete and are not used as successful timing denominators.

| Case | Native PID + FPE restart | Fixed weight, native restart | Native PID, fixed restart | Fixed weight + fixed restart |
| --- | --- | --- | --- | --- |
| blend | 880 / 2240 | 1760 / 7720 | 1280 / 3760 | 1680 / 8400 |
| brandy | 8720 / 24320 | 11960 / 27840 | cap / cap | 22080 / 63320 |
| boeing1 | 14880 / 9320 | cap / 67840 | cap / cap | cap / 67560 |

All 24 runs materialized: 18 completed and passed independent verification,
and six reached the iteration cap with nonpassing KKT. The six baseline
control runs reproduce the earlier scalar traces and final solutions exactly.

The results disconfirm a simple "PID is the culprit" explanation. In `blend`,
PID actually reduces the relative disadvantage of extra scaling. The
`blend/brandy` direction persists even with both control mechanisms fixed.
The `boeing1` direction under the existing stopping rule also persists with
both fixed. There is no single multiplicative allocation of cause among
scaling, PID, and restart: they interact, and the fixed-weight/schedule choices
are not claimed to be optimal for either pipeline.

## Important Boeing1 Correction

The additional passive original-model KKT samples revealed an earlier valid
candidate in the run **without** the extra stage. Its raw solution was then
exported without changing the ongoing solve and passed the frozen external
verifier independently:

| Iteration | Original primal | Original dual | Original gap | Reference-objective error | Internal primal |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4800 | 6.1009e-5 | 2.9462e-6 | 4.0533e-5 | 2.5084e-5 | 1.0355e-3 |
| 8000 | 6.0945e-5 | 7.4744e-8 | 3.5547e-6 | 1.8608e-5 | 8.5458e-4 |

Both candidates satisfy the original acceptance criteria, but the internal
primal metric remains above `1e-4`. The current code only attempts terminal
original-model verification after internal stopping reports optimal, so this
run continues to 14880. The additional sampler observes valid candidates but
does not terminate or change the original trace.

The metrics use different norms, normalization, and transformed versus
original constraints. Equal numeric tolerances therefore do not imply equal
tests. The internal measure is not uniformly stricter: earlier work found
the opposite situation in `agg2`, `beaconfd`, and `brandy`, where an internal
pass preceded an original-model rejection.

Consequently, **9320 versus 14880 does not establish that extra scaling reaches
original-model accuracy sooner**. We have established only that the
less-scaled trajectory has an acceptable candidate by 4800. Sampling was not
exhaustive, so the first acceptable iteration of either trajectory has not
been determined. This is not a newly measured end-to-end speedup or a change
to the formal Phase 8 result.

Source: [internal residual metric](../highs/pdlp/hipdlp/algorithm/convergence.cc#L79),
[terminal verification trigger](../highs/pdlp/hipdlp/algorithm/solver.cc#L1201),
and [original-model relative primal metric](../highs/pdlp/hipdlp/verification/kkt_verifier.cc#L255).

## What Effective Auto Scaling Should Mean

The objective should be total time to an accepted original-model solution,
including scaling and any selection work. It should not be coefficient
balance, a scaled residual, one FPE value, or iteration count under mismatched
stopping surrogates.

Three approaches have different tradeoffs:

1. A conservative fixed pipeline is cheap and predictable, but is not truly
   instance-adaptive. `Ruiz -> Pock-Chambolle` is the natural baseline to test
   here; this three-case study does not establish it as globally best.
2. A feature-based selector is cheap at solve time, but coefficient range,
   sparsity, and scale-factor statistics alone do not reveal the active set.
   It needs held-out evaluation and must not encode instance names.
3. A small, budgeted portfolio can observe the actual solver dynamics, at the
   cost of setup, memory, and extra iterations. It requires retaining candidate
   solver states and a common evaluation metric. It is more than adding a flag
   to the current fixed pipeline.

The recommended design to evaluate is a conservative baseline plus an
optional **two-candidate** portfolio, not an ever-longer chain of scalings:

- First establish a consistent original-model stopping contract. In
  particular, evaluate whether near-convergence periodic original-model checks
  can detect cases like the 4800-step candidate without relying solely on the
  internal test. This changes stopping semantics and requires separate
  correctness review; do not silently remove existing requirements.
- Retain both current scaling candidates initially. Use scale safety,
  overflow/underflow risk, and memory/setup estimates to reject unsafe or
  uneconomical trials, not to declare a numerical winner.
- For sufficiently large solves, allocate a small, explicit, capped pilot
  budget. Evaluate several checkpoints and windows using original-model
  metrics, balancing primal, dual, and gap progress against elapsed cost.
  A possible merit is

  $$
  E_k=\max\left\{
  \frac{r_p^{\mathrm{orig}}(k)}{\epsilon_p},
  \frac{r_d^{\mathrm{orig}}(k)}{\epsilon_d},
  \frac{g^{\mathrm{orig}}(k)}{\epsilon_g}\right\}.
  $$

- Avoid eliminating a candidate from a single observation. Compare stable
  progress over multiple windows, account for oscillation and stalled
  components, and use a conservative fallback when evidence is inconclusive.
  If a candidate already meets the actual acceptance contract, prefer that
  fact to an extrapolated score.
- Continue from the selected state; charge all trial work to solver time.
  Tiny cases may not justify any portfolio overhead. If switching scaling
  mid-solve, primal/dual vectors, norm bounds, anchors, and PID/FPE history need
  consistent conversion or reset; scalar `omega` adjustment alone is not enough.
- Freeze selection parameters before testing on held-out cases. Evaluate
  CPU/CUDA separately, count failures, and include setup, worst regressions,
  memory use, and decision overhead. Three chosen examples are not enough to
  validate an automatic policy.

Even original-model early ranking is not infallible. At iteration 40 in
`brandy`, the maximum original KKT component is about 308 with extra scaling
versus 507 without: a one-checkpoint chooser selects the ultimately slower
pipeline. In `boeing1`, rankings reverse multiple times. Multiple-window
selection reduces a known risk but has not been demonstrated to solve it
reliably here. There is no validated universal selection rule in this study.

## Evidence And Recheck

The study used six passive geometry runs, 24 control runs, and one additional
unchanged solve that exported the two earlier `boeing1` candidates. All work
was CPU-only with independent outputs and 180-second diagnostic budgets.
No diagnostic time is presented as a formal performance measurement.

- [Control results](hipdlp-phase8-scaling-controls.csv): all 24 outcomes,
  including the six capped cases.
- [Residual directions](hipdlp-phase8-scaling-directions.csv): scale, weight,
  effective dual step, and row violations for nine bottleneck directions.
- Full logs, scale vectors, sampled iterates, original KKT samples, raw
  solutions, verifier output, build commands, and the 247-artifact hash index:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260909-scaling-mechanism
```

Read-only recheck, including the six formal finalization seals:

```sh
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260909-scaling-mechanism/analyze-mechanism.mjs --check
```