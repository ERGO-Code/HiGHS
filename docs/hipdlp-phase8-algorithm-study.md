# HiPDLP Bounded Algorithm Study

Date: 2026-09-15

## Outcome

**Implemented and verified, but the default acceleration hypothesis failed
held-out evaluation. Default scaling, PID, eta, restart and stopping remain
unchanged.** The new, default-off `hipdlp_early_kkt_check` option retains a
bounded original-model verification trigger for explicit experiments. It is
not a generally recommended performance flag.

The candidate evaluated was Ruiz -> Pock-Chambolle scaling plus an additional
near-convergence original-KKT check every 400 iterations. It is now requested
explicitly with `hipdlp_scaling=legacy` and `hipdlp_early_kkt_check=true`.
The proposed change to `auto` was withdrawn after the frozen validation group
regressed. No further scaling, PID or trigger parameters were tuned using that
group. The six validation cases are not fresh holdouts for future tuning.

Both formal campaigns accepted 168/168 records, including warmups. The first
campaign is the held-out evaluation; the second only confirms the final public
interface reproduces the already-frozen strategy. Their combined execution
time was 34.830 seconds, well within a single 600-second final-run budget.
There were no solver timeouts, numerical failures, rejected solutions or
unexecuted formal records. The improvement target, GM <= 0.90 on the fixed
validation group, was **not met**.

| Campaign | Group | Pairs | Complete-time GM | Paired bootstrap 95% interval | Solver-time GM |
| --- | --- | ---: | ---: | --- | ---: |
| First frozen evaluation | Debug | 8/8 | 0.794608 | [0.643627, 0.943097] | 0.693975 |
| First frozen evaluation | Validation | 6/6 | 2.051281 | [0.759470, 6.587841] | 2.174114 |
| Retained-interface replication | Debug | 8/8 | 0.803722 | [0.648132, 0.958485] | 0.690917 |
| Retained-interface replication | Validation | 6/6 | 1.964229 | [0.750645, 6.203898] | 2.168836 |

The intervals resample paired instances, not timing repetitions or machine
states. Six validation instances give a wide interval and limited generality;
neither a favorable debug GM nor an interval including one rescues the failed
predeclared validation point-estimate target.

## Protocol And Execution Ledger

This study uses the current rewrite at `085a959da38c35f52e531f08253af6c9be8a67dc`
as its primary before baseline. Existing benchmarks and user changes are not
modified. The experiment root is
`/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260915-algorithm-study`.

- [x] Inspect the actual worktree and the five requested prior evidence files.
- [x] Confirm all 14 model files and that CPU5 belongs to the allowed CPU set.
- [x] Freeze input identities, references, build identities and baseline timings.
- [x] Test the stopping-tail hypothesis with passive original-model KKT samples.
- [x] Run bounded single-factor scaling, PID and stopping experiments on debug cases.
- [x] Freeze the candidate and parameters before any held-out candidate run.
- [x] Run CPU/CUDA compatibility and equivalent-LP correctness regressions.
- [x] Interleave one warmup plus five measurements per role on CPU5.
- [x] Revalidate raw evidence and report medians, fixed-denominator GM, intervals,
      failures and every regression greater than 20 percent.

The fixed debug group is afiro, agg2, beaconfd, blend, boeing1, brandy, degen2,
degen3. The held-out group is e226, bore3d, grow7, kb2, qap8, vtp.base. Derived
LPs stay with their parent group and are correctness tests, not independent
held-out performance observations.

Fast experiments use at most 2 seconds or 100000 iterations per solve and a
60-second round budget. Final measurements use at most 10 seconds per solve
and a 600-second campaign budget. Missing, failed and limited runs stay in
the denominator; incomplete paired coverage makes the group GM unavailable.
Both sides use FP64, one thread, Release, presolve off, all solver and original
KKT tolerances 1e-4, and reference objective relative error at most 1e-3.
Reference objectives are used only for external acceptance, never by the
solver or candidate selection policy. No unguarded increase of eta is allowed.

## First Falsifiable Hypothesis

Some baseline time is spent after the first original-model KKT-acceptable
candidate because `algorithm/solver.cc` only attempts original KKT validation
after the internal stopping test passes. The cheap discriminator is passive
sampling at every existing 40-step convergence boundary on a few debug cases.
It records internal and original metrics separately and exports the first
original-acceptable and first internal-acceptable raw candidates for independent
verification. If the tail is small or additional checking costs outweigh it,
the stopping-trigger candidate will not be retained on that evidence.

The first probe changes no trajectory, tolerance, restart, PID, scaling,
operator-norm guard or verifier. Its diagnostic clock readings are not
performance measurements. Algorithm trials will change one factor at a time.

## Development Evidence

The baseline preflight accepted all 14 cases in less than the 60-second round
budget. Its slowest solve was degen3 at 0.129577 seconds. The baseline records
and 1110 input/source/runtime hashes are frozen before any candidate run.

The passive four-case probe exactly reproduced baseline final solution bytes
and iterations. Under auto scaling, the first original KKT pass equalled final
exit on afiro, blend, boeing1 and brandy. The first internal pass for brandy
was earlier, at 23680 versus final 24320; removing verification would be wrong.
Thus a stopping-only explanation of these current auto trajectories is rejected.

The first 8-case single-factor round accepted 64/64 records. These are one
measurement after one warmup, diagnostic comparisons rather than final medians:

| Factor | Solver-time GM | Complete-time GM | Decision |
| --- | ---: | ---: | --- |
| Ruiz + PC, without extra Curtis-Reid | 0.781342 | 0.885370 | Continue to interaction test; boeing1 regresses |
| PID log gains multiplied by 0.5 | 1.035231 | 1.052917 | Reject as default candidate |
| Near-convergence KKT probe every 40 steps | 1.162601 | 1.141127 | Reject due to verification overhead |

A second stopping-only candidate keeps the same thresholds but allows an
extra probe only every 400 steps, limiting full recovery/verification work.
It leaves all 8 exit iteration counts unchanged and is considered only as
a possible interaction with scaling. No eta or restart rule is changed.
All sampled initial eta values are 0.998 with certified squared norm bound 1.
PID never falls back in the four initial diagnostic trajectories, so the
uninitialized best-gap observation does not explain their performance.

## Retained Code

- The public boolean `hipdlp_early_kkt_check` defaults to false. Option
    translation passes it to the backend-independent solver; it selects no
    scaling pipeline and changes no PID or restart parameters.
- `shouldVerifyOriginalKkt` requests an extra check only at a positive multiple
    of 400 iterations, while internal stopping says continue, with the internal
    gap at its requested tolerance and both internal feasibility metrics at most
    the square root of their requested tolerances. These are **trigger**
    thresholds, not acceptance thresholds.
- The extra path is disabled if the original KKT tolerance is looser than any
    requested internal tolerance. Internal-optimal verification, numerical-error
    handling, iteration/time-limit precedence and certificate verification remain
    intact. The unchanged original-model verifier must accept before `OPTIMAL`
    can be returned; a rejected periodic candidate continues normally.
- `auto` still means Curtis-Reid -> enabled Ruiz/Pock-Chambolle stages. The
    scaling implementation is byte-identical to the before source. The explicit
    `legacy` scaling mode was already available; no instance-name, dimensions,
    reference objective or early winner selection enters the algorithm.
- The certified operator-norm bound, safety factor 0.998, finite-value checks,
    PID update formula and restart policy are unchanged. There is no unguarded
    eta increase or online portfolio cost to exclude from timing.

The added option changes `HighsOptions` and internal `SolverOptions` layouts.
Internal binary consumers must be rebuilt with matching headers and libraries.
Archived drivers continue using their own archived bundles.

Relevant code: [verification trigger](../highs/pdlp/hipdlp/algorithm/convergence.cc#L77),
[solver integration](../highs/pdlp/hipdlp/algorithm/solver.cc#L1199),
[public option](../highs/lp_data/HighsOptions.h#L1470),
[solver regressions](../check/TestHiPdlpSolver.cpp#L377).

## Timing Contract

Before is the current two-SpMV/fused-CPU rewrite, not legacy. Both sides were
remeasured in each campaign on allowed CPU5, with adjacent before/after jobs
and alternating role order by case. Each job performs one warmup and five
fresh measured solves in one process. Warmups must pass correctness but do not
enter medians. FP64, one thread, presolve off, PID, Halpern restart, Clang 18,
Release `-O3 -DNDEBUG -g`, no LTO and no testing hooks are matched.

The predeclared primary **complete time** is solver-wrapper wall time plus
independent verifier process wall time. It charges scaling, algorithm setup,
norm estimation, iteration/check work, solution recovery, in-solver final
verification and local teardown, then the external verifier's startup, model
read, checks and output. Initial model reading by the solver and raw solution
serialization are outside that primary interval and retained separately as
read/repeat/process wall fields. External verification costs are deliberately
not subtracted; they dominate sub-millisecond cases. Solver time is reported
separately so this overhead cannot conceal algorithmic regressions.

No diagnostics or preload hooks run in formal measurements. There is no
frequency/power control or guarantee about competing machine load. The second
campaign is not pooled with the first, and no old benchmark times serve as
denominators. Final per-solve limits are 10 seconds with a 600-second campaign
watchdog. Fast rounds use 2 seconds or 100000 iterations and 60 seconds per
round; these limits are not expanded after observing a slow case.

## Per-Instance Results

The following are five-measurement medians from the **retained-interface
replication**, in milliseconds. `Opt-in` is the explicit scaling-plus-trigger
combination, not the retained default. Iterations are before/opt-in and were
constant across all six repeats. Every row has 12/12 independently accepted
records and zero failures. The first campaign's full-precision values remain
separate in [the paired CSV](hipdlp-phase8-algorithm-study.csv).

| Group | Case | Before complete | Opt-in complete | Ratio | Before solver | Opt-in solver | Iterations |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Debug | afiro | 7.511 | 7.628 | 1.016 | 0.287 | 0.221 | 400/280 |
| Debug | agg2 | 40.183 | 38.789 | 0.965 | 30.724 | 29.451 | 3760/4200 |
| Debug | beaconfd | 28.164 | 22.030 | 0.782 | 19.566 | 13.322 | 3640/2360 |
| Debug | blend | 10.875 | 9.262 | 0.852 | 3.083 | 1.461 | 2240/880 |
| Debug | boeing1 | 69.346 | 42.600 | 0.614 | 60.617 | 33.344 | 9320/4800 |
| Debug | brandy | 100.887 | 42.910 | 0.425 | 92.487 | 34.555 | 24320/8720 |
| Debug | degen2 | 31.671 | 31.546 | 0.996 | 22.375 | 22.336 | 2640/2640 |
| Debug | degen3 | 149.402 | 153.058 | 1.024 | 128.461 | 136.437 | 3320/3320 |
| Validation | e226 | 41.957 | 65.910 | 1.571 | 30.046 | 57.288 | 5720/13400 |
| Validation | bore3d | 36.800 | 1041.189 | 28.293 | 28.459 | 1033.091 | 7640/297200 |
| Validation | grow7 | 67.955 | 24.938 | 0.367 | 59.236 | 16.156 | 14720/3600 |
| Validation | kb2 | 20.564 | 20.719 | 1.008 | 13.002 | 13.193 | 18600/18400 |
| Validation | qap8 | 28.355 | 27.790 | 0.980 | 16.308 | 16.409 | 600/600 |
| Validation | vtp.base | 19.553 | 69.721 | 3.566 | 11.442 | 61.783 | 4280/25160 |

All formal regressions above 20 percent, with no omissions:

| Campaign | Case | Complete ratio | Solver ratio |
| --- | --- | ---: | ---: |
| First evaluation | e226 | 1.752474 | 1.920703 |
| First evaluation | bore3d | 29.006374 | 36.508501 |
| First evaluation | vtp.base | 3.813058 | 5.432113 |
| Interface replication | e226 | 1.570891 | 1.906658 |
| Interface replication | bore3d | 28.293328 | 36.300876 |
| Interface replication | vtp.base | 3.565807 | 5.399897 |

The correctness table gives maxima over all six opt-in repetitions. Before
metrics, both campaigns, reference errors and failure counts are also in the
paired CSV. No verifier or acceptance tolerance was changed.

| Case | Original primal | Original dual | Original gap | Reference error | Failed records |
| --- | ---: | ---: | ---: | ---: | ---: |
| afiro | 3.741e-6 | 8.494e-6 | 1.014e-5 | 8.355e-6 | 0 |
| agg2 | 1.625e-6 | 4.714e-7 | 1.565e-5 | 3.143e-5 | 0 |
| beaconfd | 8.989e-7 | 1.040e-5 | 1.281e-5 | 7.813e-8 | 0 |
| blend | 5.548e-5 | 9.782e-6 | 4.289e-5 | 6.391e-5 | 0 |
| boeing1 | 6.101e-5 | 2.946e-6 | 4.053e-5 | 2.508e-5 | 0 |
| brandy | 5.097e-5 | 9.022e-5 | 1.047e-5 | 1.807e-5 | 0 |
| degen2 | 4.989e-5 | 4.286e-5 | 1.248e-5 | 3.130e-5 | 0 |
| degen3 | 5.741e-5 | 6.872e-5 | 1.044e-5 | 1.428e-5 | 0 |
| e226 | 3.491e-5 | 9.052e-6 | 7.240e-5 | 1.795e-4 | 0 |
| bore3d | 7.619e-5 | 3.812e-8 | 8.362e-5 | 4.428e-6 | 0 |
| grow7 | 3.726e-10 | 1.748e-9 | 1.617e-9 | 4.021e-9 | 0 |
| kb2 | 8.489e-5 | 2.396e-6 | 2.390e-7 | 2.660e-7 | 0 |
| qap8 | 1.096e-5 | 1.196e-5 | 6.667e-7 | 2.432e-8 | 0 |
| vtp.base | 6.753e-7 | 9.964e-5 | 1.562e-5 | 6.904e-5 | 0 |

## Single-Factor Results

Each factor was first checked on afiro, blend, boeing1 and brandy, then on
the fixed eight-case debug group. These quick results use one measured solve
after one warmup, not the final five-repeat protocol. All 144 records across
the three full debug rounds were accepted; all 72 measured rows are retained
in [the ablation CSV](hipdlp-phase8-algorithm-ablations.csv).

| Factor / round | Solver GM | Complete GM | Interpretation |
| --- | ---: | ---: | --- |
| Scaling only | 0.781342 | 0.885370 | Helps several debug cases, hurts boeing1 |
| PID log gains x0.5 only | 1.035231 | 1.052917 | Rejected; damping is not a general improvement |
| Near-convergence probe every 40 steps only | 1.162601 | 1.141127 | Rejected; recovery/verification overhead dominates |
| Same probe every 400 steps only | 1.019161 | 0.988362 | No reduced exit counts; no default speedup claim |
| Scaling only, combination-control round | 0.780104 | 0.874821 | Reproduces boeing1 regression |
| Scaling + 400-step probe | 0.695998 | 0.798877 | Chosen on debug evidence, then rejected as default on validation |

Quick-round regressions above 20 percent in either metric are preserved:
scaling-only boeing1 (solver +55.7%, complete +50.6%; control round +55.1%,
+48.3%); PID damping beaconfd (+44.9%, +35.4%) and brandy (+20.8% solver);
40-step probing brandy (+40.8%, +38.5%), degen2 (+21.5% solver) and degen3
(+65.3%, +58.5%). Neither the 400-step-only nor combined debug round has a
greater-than-20-percent regression. Cold exploratory samples remain in their
original directories and are not substituted for these warmed comparisons.

## Mechanism And First Hits

The 32 passive trajectories record eta, tau, sigma, omega, pre/post PID state,
anchor distances, FPE, restart condition/reason and interval at the existing
40-step boundaries. They also record internal and original-model metrics
separately, export the first original/internal candidates, and independently
verify the original-KKT first-hit solutions. Diagnostic elapsed times are
retained in `analysis.json`, not used as performance times.

The entries below are first original-KKT passing boundary / first internal
passing boundary / final exit. A dash means the run terminated before any
internal pass. These are first **sampled 40-step boundaries**, not claims about
unobserved intervening iterates.

| Debug case | Before auto | Scaling only | Combined candidate | Restarts before / combined |
| --- | --- | --- | --- | ---: |
| afiro | 400/400/400 | 280/280/280 | 280/280/280 | 8/6 |
| agg2 | 3760/2800/3760 | 4200/4200/4200 | 4200/4200/4200 | 15/16 |
| beaconfd | 3640/3600/3640 | 2360/2360/2360 | 2360/2360/2360 | 15/15 |
| blend | 2240/2240/2240 | 880/880/880 | 880/880/880 | 11/7 |
| boeing1 | 9320/9320/9320 | 4640/14880/14880 | 4640/-/4800 | 15/12 |
| brandy | 24320/23680/24320 | 8720/8720/8720 | 8720/8720/8720 | 22/15 |
| degen2 | 2640/2640/2640 | 2640/2640/2640 | 2640/2640/2640 | 10/10 |
| degen3 | 3240/3320/3320 | 3240/3320/3320 | 3240/3320/3320 | 9/9 |

This distinguishes two effects. Scaling changes the path and restart/PID
observations; periodic verification only discovers an already-acceptable
candidate. In boeing1, the first combined candidate is observed at 4640 and
the bounded trigger finds a passing candidate at 4800. In degen3, the passing
3240 candidate is between 400-step checks, so the retained trigger does not
save iterations and still adds verification work.

All sampled initial eta values are 0.998 with certified squared norm bound 1.
The trace checker verifies tau = eta/omega and sigma = eta*omega for every
recorded initialization/update. Large omega jumps are real but not by
themselves evidence of a defective controller: halving the log gains worsened
the debug GM. No PID fallback occurred in the native before/scaling/combined
traces. Restart counts are observations, not an additive attribution of cause.

## Correctness And Compatibility

| Final retained-code check | Cases | Assertions | Result |
| --- | ---: | ---: | --- |
| CPU relevant algorithm/solver regression union | 141 | 6957 | Pass |
| CUDA same union, including actual device execution | 156 | 15772 | Pass |
| Native default preservation against frozen before | 14 LPs | Iterations and all raw solution bytes | Pass |
| Two formal campaigns | 336 records | Independent KKT and reference objective | Pass |

Tests include scalar trigger boundaries, finite inputs, limits, stricter
internal tolerance preservation and public option translation. A deterministic
two-row LP proves that an internal-continue candidate can pass periodic KKT
at iteration 400; stricter KKT rejects the same probe and continues to the
401-step limit. Disabling the feature preserves the old behavior. Before
failed four integration assertions; the enabled implementation passes on
both CPU and CUDA.

Six small equivalent-LP variants cover unchanged data, positive row scaling,
positive column scaling, row/column permutations, objective scaling and a
combined transformation. Their recovered solutions and status pass original
KKT on both backends. These are correctness fixtures, not extra independent
validation instances. CUDA performance on large problems, the full HiGHS
test suite and the full release build matrix were not evaluated.

## Evidence And Recheck

Old benchmark outputs, plans, reports and user modifications were preserved;
no commit, branch change or destructive git operation was performed. New
source/runtime snapshots, compilation commands, references, model hashes,
every raw solution, driver/verifier record, diagnostic event and rejected
experimental candidate remain under the study root above.

The formal harness requires exact 8+6 case coverage, native roles, frozen
parameters and source/driver/library hash coverage. Reported iteration counts
must obey the configured cap; a failed, missing or rejected record makes its
group GM unavailable. Acceptance self-tests reject denominator shrinkage,
changed frozen parameters, missing candidate pins, diagnostics in formal runs,
over-budget iterations and numeric/protocol failures.

The independent recheck verifies 336/336 formal records, 5240 hash comparisons,
14 unchanged-default cases and 32 passive trajectories. Historical candidate
source hashes are checked against the frozen source snapshots, because the
working tree deliberately withdrew that default change. Retained-interface
hashes are additionally checked against the current production sources.

```sh
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260915-algorithm-study/analyze.mjs
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260915-algorithm-study/retained/study.mjs self-test
```

## Unresolved

There is no supported general scaling selector, damped-PID default or safe
larger-eta strategy from this study. The validation regressions are not
attributed solely to PID or to a verifier mismatch, and their root mechanisms
were not used to retune on the holdout. The retained trigger can add recovery
and verification cost, especially on CUDA; enabling it globally is not
recommended by these data. It remains opt-in because its stopping semantics
are correct and its benefit is reproducible on specific trajectories, not
because cross-instance performance has been established.

Future algorithm selection work needs a new held-out set and a cheap,
numerically motivated rule that survives cases such as bore3d. No claim is
made of universal LP acceleration, an improved default, or a Phase 8 release
performance pass.