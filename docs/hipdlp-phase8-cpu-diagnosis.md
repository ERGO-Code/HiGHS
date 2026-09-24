# Phase 8 CPU Slowdown Diagnosis: Nine Cases

Date: 2026-09-08

This is a diagnostic study, not a solver change or a new acceptance campaign.
The frozen Phase 8 verdict remains unchanged. The requested cases are `afiro`,
`agg2`, `beaconfd`, `blend`, `boeing1`, `brandy`, `degen2`, `degen3`, and
`greenbeb`.

Follow-up: [scaling mechanism and controlled experiments](hipdlp-phase8-scaling-mechanism.md)
qualify the `boeing1` interpretation below. The rewrite legacy-scaling run
contains an externally verified acceptable original-model candidate at
iteration 4800, although its internal stopping test makes it run to 14880.
Reported completion counts are therefore not necessarily the first time each
trajectory reaches the common original-model accuracy contract.

## Main Findings

1. There are two distinct contributors: different convergence trajectories
   and a more expensive CPU implementation. Neither iteration counts nor one
   step-size scalar explains all nine regressions.
2. Initial primal weight is the same to floating-point precision in all nine
   cases. Weight trajectories diverge after the first restart because the
   scaled iterates and anchor distances differ.
3. Changing only rewrite's scaling option from `auto` to `legacy` removes most
   of the extra iterations for `afiro`, `blend`, and `brandy`, and improves
   `beaconfd`. It worsens `agg2` and `boeing1`, and changes nothing in the full
   scalar trajectories of `degen2` and `degen3`. A global scaling switch is not
   justified by these results.
4. `degen2`, `degen3`, and `boeing1` expose implementation cost particularly
   clearly. Rewrite has similar or fewer iterations, but is still slower.
5. Original-model KKT verification rejects intermediate rewrite candidates in
   `agg2`, `beaconfd`, and `brandy`. Removing those checks would violate the
   requested accuracy, not provide a valid speedup.

## Evidence And Method

- Formal CPU times below are the existing five-measurement medians after one
  warmup, using the exact frozen libraries from the legacy-qualified campaign.
- All stopping, primal, dual, KKT, and independent KKT tolerances are `1e-4`.
  PID, Halpern restart, FP64, presolve off, one thread, and seed `20260908` are
  retained. Reference-objective relative error must remain at most `1e-3`.
- The diagnostic programs link the same frozen Clang CPU libraries and loader
  bundles. Function wrappers call the original exported functions and record
  their existing scalar inputs and results; no production algorithm is patched.
- Twenty-seven trajectory runs cover legacy, rewrite-auto, and
  rewrite-legacy-scaling for every case. All 27 pass independent original-model
  KKT and reference-objective verification. All 18 baseline trajectory runs
  reproduce the formal iteration count and raw solution bytes exactly.
- Nine additional rewrite runs enable the existing `SolveDiagnostics` and
  record original-model KKT attempts. All nine preserve the formal solution
  bytes and iteration counts. These phase timings are single instrumented cold
  samples, not matched performance medians.
- Diagnostic solves have a separate 180-second budget to avoid censoring a
  trajectory due to observation overhead. Their times do not replace the
  formal 60-second acceptance measurements. No GPU runs were performed.

The formal source root is:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260908T205450Z-legacy-qualified-1e4
```

Diagnostic sources, compile commands, complete logs, solutions, verifier
outputs, every-40-step JSON traces, KKT rejections, and the hash index are under:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260908-cpu-nine-diagnostics
```

Repository-sized derived data:

- [Summary CSV](hipdlp-phase8-cpu-diagnostic-summary.csv): full-precision timing,
  step-size, weight, restart, residual, and phase-summary fields for nine cases.
- [Checkpoint CSV](hipdlp-phase8-cpu-diagnostic-checkpoints.csv): 333 exact
  iteration checkpoints across the three variants, without interpolation.
  Complete trajectories remain in the diagnostic root.

## Time Versus Iterations

Let `R` be rewrite-auto and `L` be frozen legacy. The exact descriptive
decomposition is

$$
\frac{T_R}{T_L} = \frac{N_R}{N_L}
  \frac{T_R/N_R}{T_L/N_L}.
$$

`T/N` is amortized complete-solve time per iteration, including initialization,
checks, postsolve, and teardown. It is not a pure step or SpMV timer, especially
for `afiro`.

| Case | Legacy ms | Rewrite ms | Legacy iterations | Rewrite iterations | Time ratio | Iteration ratio | Amortized T/N ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| afiro | 0.08043 | 0.67186 | 280 | 400 | 8.353 | 1.429 | 5.847 |
| agg2 | 12.99172 | 38.62794 | 3080 | 3760 | 2.973 | 1.221 | 2.436 |
| beaconfd | 10.59527 | 25.37936 | 3400 | 3640 | 2.395 | 1.071 | 2.237 |
| blend | 1.30562 | 3.82119 | 840 | 2240 | 2.927 | 2.667 | 1.098 |
| boeing1 | 58.40362 | 78.79652 | 15520 | 9320 | 1.349 | 0.601 | 2.247 |
| brandy | 18.36744 | 119.48834 | 8480 | 24320 | 6.505 | 2.868 | 2.268 |
| degen2 | 11.37177 | 27.67880 | 2520 | 2640 | 2.434 | 1.048 | 2.323 |
| degen3 | 74.29333 | 165.83441 | 3240 | 3320 | 2.232 | 1.025 | 2.178 |
| greenbeb | 23729.52026 | 56687.18673 | 656680 | 773440 | 2.389 | 1.178 | 2.028 |

`blend` is primarily a convergence regression. `brandy` combines roughly 2.87x
as many iterations with 2.27x amortized cost. `boeing1` is the opposite:
rewrite's convergence improvement is outweighed by implementation cost.

## Step Size And Primal Weight

The quantities recorded before a restart are the values used to generate that
boundary's candidate:

$$
\tau=\eta/\omega,\qquad \sigma=\eta\omega,\qquad
\eta=\sqrt{\tau\sigma},\qquad \omega=\sqrt{\sigma/\tau}.
$$

Weight-update events separately record the new values installed for later
iterations. The collector verifies both identities and checks that an updated
weight appears at the next boundary.

Both implementations use safety factor `0.998`. Legacy uses its power-method
estimate for the squared norm. Rewrite uses the maximum of its estimate and a
certified upper bound; the final Pock-Chambolle stage gives the bound `1` in
these runs. Consequently rewrite starts with `eta=0.998` in every case.

| Case | Initial omega, both implementations | Legacy eta | Rewrite eta |
| --- | ---: | ---: | ---: |
| afiro | 0.01317475934 | 1.021225674 | 0.998 |
| agg2 | 0.0002694409631 | 1.009983925 | 0.998 |
| beaconfd | 0.04749354479 | 0.998247444 | 0.998 |
| blend | 0.2788982491 | 1.036464818 | 0.998 |
| boeing1 | 0.6833586562 | 1.004797837 | 0.998 |
| brandy | 0.01016130404 | 1.004461855 | 0.998 |
| degen2 | 4.732044461 | 1.073452351 | 0.998 |
| degen3 | 1.420162030 | 1.022663600 | 0.998 |
| greenbeb | 1343.317562 | 1.012909126 | 0.998 |

For example, `afiro` starts with legacy `(tau, sigma) =
(77.5138010412, 0.01345440249)` versus rewrite
`(75.7509093151, 0.01314840982)`. Initial weight is identical. The scalar bound
difference is small relative to the observed 8.35x full-solve ratio.

The first PID update already shows trajectory-dependent differences:

| Case | Legacy omega after iteration 40 | Rewrite omega after iteration 40 |
| --- | ---: | ---: |
| afiro | 1.70191e-3 | 2.94211e-3 |
| agg2 | 1.80923e-4 | 5.93555e-4 |
| beaconfd | 8.29453e-3 | 3.73544e-2 |
| blend | 4.33358e-1 | 5.10960e-1 |
| boeing1 | 6.01413e-2 | 7.06588e-2 |
| brandy | 1.93454e-2 | 7.68330e-3 |
| degen2 | 6.89858 | 6.90855 |
| degen3 | 2.02945 | 2.02017 |
| greenbeb | 6.87541e2 | 1.64256e3 |

Different diagonal scaling means that equal scalar `tau` or `omega` would not
imply equal effective updates in original variable coordinates. These
differences establish divergence, not a defect in the PID controller. The
norm-bound contribution has not been isolated by overriding the safety guard.

Source: [step initialization](../highs/pdlp/hipdlp/algorithm/solver.cc#L957),
[scalar step policy](../highs/pdlp/hipdlp/algorithm/step_size_policy.cc#L71), and
[legacy initialization/PID](../highs/pdlp/hipdlp/pdhg.cc#L1987).

## Residuals At Matched Iterations

The following are the solver-internal relative `(primal, dual, gap)` values.
The legacy table's raw `pinf/dinf` have been normalized using the actual
unscaled problem norms, not compared directly with rewrite's relative values.
These are not the separate original-model KKT metrics.

Each row uses the first iteration at which either baseline implementation
finishes. Thus `boeing1` is sampled at rewrite's finish, while the others are
sampled at legacy's finish.

| Case | Iteration | Legacy: primal / dual / gap | Rewrite: primal / dual / gap |
| --- | ---: | --- | --- |
| afiro | 280 | 9.132e-6 / 4.122e-6 / 2.433e-5 | 2.480e-3 / 1.159e-5 / 2.774e-3 |
| agg2 | 3080 | 4.352e-6 / 4.486e-7 / 3.824e-5 | 9.312e-4 / 1.035e-5 / 4.308e-4 |
| beaconfd | 3400 | 4.906e-8 / 5.826e-6 / 2.015e-5 | 1.224e-3 / 5.878e-10 / 3.022e-4 |
| blend | 840 | 5.773e-5 / 1.465e-5 / 9.929e-5 | 7.038e-3 / 1.997e-3 / 3.342e-2 |
| boeing1 | 9320 | 8.542e-4 / 3.033e-7 / 5.640e-6 | 3.003e-5 / 1.073e-7 / 9.945e-5 |
| brandy | 8480 | 5.675e-5 / 9.165e-5 / 1.130e-5 | 3.724e-2 / 2.407e-5 / 9.492e-4 |
| degen2 | 2520 | 1.277e-5 / 2.564e-5 / 1.240e-6 | 4.747e-4 / 2.584e-5 / 5.905e-5 |
| degen3 | 3240 | 9.748e-5 / 5.562e-5 / 1.357e-5 | 1.211e-4 / 4.166e-5 / 1.508e-5 |
| greenbeb | 656680 | 9.868e-5 / 1.061e-7 / 3.018e-8 | 5.399e-2 / 1.054e-6 / 4.922e-6 |

`brandy` and `greenbeb` are dominated by primal feasibility at this point, not
the duality gap. `degen3` is much closer: primal residual is only slightly over
the threshold, accounting for two additional 40-step blocks rather than its
2.23x runtime. At `greenbeb` iteration 656680, weight is approximately `0.04351`
in legacy and `0.007446` in rewrite; this is an observed trajectory difference,
not proof that a particular weight update caused the delay.

## Scaling-Only Control

Only rewrite's supported `hipdlp_scaling` option is changed. Tolerances,
library, initialization safety rule, and PID policy remain fixed.

| Case | Legacy iterations | Rewrite auto | Rewrite legacy scaling | Restarts: legacy / auto / scaling control |
| --- | ---: | ---: | ---: | --- |
| afiro | 280 | 400 | 280 | 6 / 8 / 6 |
| agg2 | 3080 | 3760 | 4200 | 15 / 15 / 16 |
| beaconfd | 3400 | 3640 | 2360 | 15 / 15 / 15 |
| blend | 840 | 2240 | 880 | 7 / 11 / 7 |
| boeing1 | 15520 | 9320 | 14880 | 18 / 15 / 17 |
| brandy | 8480 | 24320 | 8720 | 15 / 22 / 15 |
| degen2 | 2520 | 2640 | 2640 | 10 / 10 / 10 |
| degen3 | 3240 | 3320 | 3320 | 9 / 9 / 9 |
| greenbeb | 656680 | 773440 | 755120 | 58 / 81 / 78 |

This control supports a scaling-driven convergence regression for `afiro`,
`blend`, and `brandy`. It also shows why neither a universal scaling rollback
nor restart count alone is an adequate explanation. `beaconfd` changes from
3640 to 2360 iterations with the same 15 restarts. Both `degen2` and `degen3`
have byte-identical complete scalar traces and final solutions under the two
rewrite scaling modes.

## Original-Model KKT Rejections

Rewrite deliberately verifies the recovered solution after the internal
stopping test passes. A failed original-model check continues the solve.

| Case | First internal pass | Final accepted iteration | Rejected KKT attempts | First rejection on original model |
| --- | ---: | ---: | ---: | --- |
| agg2 | 2800 | 3760 | 9 | dual residual 3.23721e-4 |
| beaconfd | 3600 | 3640 | 1 | primal residual 1.47035e-4 |
| brandy | 23680 | 24320 | 13 | primal residual 2.10611e-4 |

No intermediate original-model KKT rejection occurred in the other six
rewrite-auto runs. In `agg2`, for example, stopping at 2800 would appear faster
than legacy but return a solution outside the independent `1e-4` dual tolerance.
In `brandy`, these checks account for the final 640-iteration tail, not the
much larger preceding convergence difference.

Source: [optimal-candidate verification](../highs/pdlp/hipdlp/algorithm/solver.cc#L1201).

## CPU Execution Cost

The source and the separate phase samples point to concrete costs:

- [The shared rewrite step](../highs/pdlp/hipdlp/algorithm/halpern_pdhg_step.cc#L31)
  computes `A * reflected_primal`, `A * next_primal`, and `A^T * next_dual`:
  three SpMVs per step.
- [Legacy's ordinary CPU step](../highs/pdlp/hipdlp/pdhg.cc#L961) computes only
  `A * reflected_primal` and `A^T * current_dual`. It combines projection and
  reflection in the same loops and only stores the candidate at major steps.
- Rewrite separately performs projection, copy, scale, axpy, and blend. The
  [CPU backend](../highs/pdlp/hipdlp/backend/cpu/cpu_backend.cc#L346) additionally
  validates handles/dimensions, checks scale/axpy results, and scans outputs for
  nonfinite values. These are real extra passes, not just function-call syntax.
- Repeated certificate verification is substantial CPU work. For `greenbeb`,
  the diagnostics count 19,336 certificate probes, and host verification takes
  25.6% of the instrumented solve. That phase includes
  [ray and Farkas verification](../highs/pdlp/hipdlp/algorithm/solver.cc#L403),
  not just the final optimality KKT check.

| Case | Setup % | Step % | Boundary/check % | Certificate probes |
| --- | ---: | ---: | ---: | ---: |
| afiro | 28.0 | 42.4 | 22.1 | 10 |
| agg2 | 12.8 | 63.5 | 22.8 | 94 |
| beaconfd | 12.3 | 67.2 | 18.9 | 91 |
| blend | 13.3 | 56.4 | 24.6 | 56 |
| boeing1 | 4.5 | 72.0 | 22.3 | 233 |
| brandy | 1.8 | 70.3 | 25.9 | 608 |
| degen2 | 13.9 | 62.7 | 22.4 | 66 |
| degen3 | 15.3 | 65.2 | 19.2 | 83 |
| greenbeb | 0.1 | 71.1 | 28.8 | 19336 |

Setup combines preprocessing, backend setup, and norm estimation. The columns
need not total 100% because loop/control and other work remain. Host
verification is nested inside boundary time and must not be added to it again.
Instrumentation, cold caches, and clock reads affect these single-run shares;
they localize work but do not measure a speedup or apportion the entire legacy
regression exactly. In particular, `afiro` has only a 0.591 ms absolute formal
runtime difference despite its large ratio.

## Next Focused Experiments

1. Use `degen3` as the CPU execution-cost case: scaling does not change its
   trace, and convergence counts are already close. Measure the extra Ax and
   vector passes separately before changing fusion or freshness rules.
2. Use `blend` and `brandy` for scaling/PID/restart trajectory analysis. The
   scaling-only control is already informative. Compare the first divergent
   restart, pre-update anchor distances, FPE/reference ratios, and installed
   weights while holding scaling fixed.
3. Use `greenbeb` for the long primal-feasibility tail and repeated certificate
   costs. Its small scaling-only improvement and zero KKT rejections rule out
   two simple explanations.
4. Isolate the norm-initialization difference only in a separately verified
   experiment. Do not remove the certified bound, final original-model KKT,
   or finite-value checks merely to match legacy time.

The data does not yet establish a PID formula bug, a specific restart-threshold
bug, or the exact percentage attributable to each CPU primitive. The complete
traces retain actual restart iterations, FPE values, and pre/post weights for
that next comparison. Legacy's optimized restart-condition function was not
dynamically observable, so its actual restart events were captured at the PID
restart entry; unobserved legacy condition reasons are not invented.

## Recheck

The analysis script revalidates six formal seals, all baseline solution
identities, the 27 independent trajectory verifications, the nine profile
solutions, and 198 diagnostic artifact hashes before recomputing the CSVs:

```sh
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260908-cpu-nine-diagnostics/analyze-cpu.mjs --check
```

The study leaves all frozen acceptance outputs and production solver code
unchanged. Only this report and the two derived CSVs are added to the repository.