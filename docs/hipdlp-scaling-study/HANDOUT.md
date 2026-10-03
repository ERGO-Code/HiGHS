# HiPDLP scaling study — handout (2026-10-01/02)

Branch `feature/hipdlp-rewrite`. Benchmarks: `/srv/mps_da` (read-only, shared; copy/gunzip to a scratch dir).
Hardware used: RTX A4500 + i9-14900K. Tolerance 1e-4 everywhere, fixed step size, presolve off.

## 1. What is committed

- **Curtis-Reid guard** (`highs/pdlp/hipdlp/model/scaling.{hpp,cc}`, `problem.{hpp,cc}`,
  `algorithm/solver.cc`, test in `check/TestHiPdlpModel.cpp`, docs in `docs/src/solvers.md`,
  `HIPDLP_REWRITE_PLAN.md`). In `applyScalingPipeline`, `auto` drops CR when
  max over {costs, finite row lower bounds, finite column bounds} of
  `(1+||x'_CR+Ruiz+PC||2)/(1+||x'_Ruiz+PC||2) > kCurtisReidInflationLimit (1e3)`;
  sets `ScalingMetadata::curtis_reid_rejected` and logs it.
  - Validated: neso-2245-filtered auto -> 106,080 its / 108 s optimal (was: no convergence in 714k its).
  - 232-instance CPU sweep: 229/232 identical; 3 differ only at the 30 s time-limit boundary.
    Guard trips only on neso-* and ns1688926 (unsolved by every config).
  - `[hipdlp]` unit tests: new test passes; 7 pre-existing failures (TestHiPdlpBackend 1718/1883/3714,
    TestHiPdlpTrace 862) fail identically on unmodified HEAD.
  - Threshold calibrated on one family (next-largest inflation 10^2.75, neso 10^3.9-10^4.8).
- `tools/` — scratch drivers/analysis (not built by CMake). `adapter_x.cpp` = migration adapter +
  `--scaling-mask --ruiz-iters --cr-iters --drop-tol --dump-solution`. Build against a build dir, e.g.
  `g++ -O2 -std=c++17 -I highs -I build -I highs/lp_data -I highs/util -I highs/io -I highs/pdlp -I extern
   -I highs/interfaces tools/X.cpp -L build/lib -lhighs -Wl,-rpath,build/lib`. Python scripts have a
  hard-coded `S=` scratch path — edit it.
- `data/` — sweep results (`cpu_all.csv`: 232 inst x 8 configs, CPU, 30 s, 14 parallel), features,
  guard decisions, complexity-factor proxies (true solution and pilots).
  Config codes: A=CR+Ruiz+PC (auto), L=Ruiz+PC (legacy), P=PC, R=Ruiz, C=CR, O=off, AP=CR+PC, AR=CR+Ruiz.

## 2. Findings

1. **neso failure is CR, not tiny entries.** neso-2245 and -filtered are the same LP after
   `readModel` (HiGHS drops |a|<1e-9). Median |a| 1.3e-4; rows with ~2000 nnz and one unit entry.
   CR's geometric mean gives rows R~1.4e5, their columns C~1e-8; stopping test is in unscaled space,
   so 4 columns carry 97% of unscaled dual violation (amplified 3.6e4x). Dropping |a|<1e-5 doesn't help;
   robust CR (rho=1e-2) helps but not enough. Ruiz+PC 106k its, PC 90k, auto never.
2. **Legacy vs rewrite on neso is a stopping-standard difference.** Both internal criteria are
   cuPDLP-style (2-norm/(1+||b||,||c||), boxed columns get a free dual slack; gap binds): legacy stops
   at 21k, rewrite internal passes at 31,440. The rewrite then requires `verifyKkt` (max-norm,
   boxed sign by midpoint): 35 boxed zero-cost interior columns have rc ~1e-3, so it is rejected until
   a restart cascade at 105,600/105,840; accepted at 106,080. Legacy at tol 3e-6 stops at 104,720.
3. **Why scaling matters.** Complexity factor `Kw = sigma(A')*(w*||(x*-x0)/C||^2+||y*/R||^2/w)/2`,
   w=(1+||c||)/(1+||b||) unscaled (fixed-step weight). On neso it orders all 8 configs over 14 decades.
   Over 122 instances, with true x*,y*, it picks the faster of L/A 63/76. Pipelines ending in PC have
   sigma~1, so the difference is the distance term. More stages is not monotone (blend: L 1800 vs A 7920;
   bandm: L 82880 vs A 5080).
4. **Auto policy evaluation (122 solvable instances, deterministic CPU iterations):**
   A 21,945 SGM / 111 solved; L 23,880 / 115; P 38,422 / 107; oracle(L,A) 13,179 / 122.
   PC-only default rejected. Static features do not predict L vs A.
   **Best candidate (not implemented):** two 200-iteration pilots (A and L), estimate Kw of each
   candidate from both pilot iterates, take the geometric mean, pick min -> GM 0.80 vs always-A
   (95% CI 0.68-0.94), 118/122 solved, picks L on both neso. Single pilot is self-biased.
   Worst regression: stair 18x (uninvestigated).

## 3. Open next steps

- Prototype the two-pilot selector in `algorithm/solver.cc` (needs two solver states or a restart),
  validate on held-out instances + neso on CUDA; consider warm-starting the chosen scaling.
- Decide whether `verifyKkt`'s boxed-column sign rule is the intended acceptance bar (3.4x iterations on neso).
- Investigate stair; PID step size changes the w term of Kw.

## 4. Continuation (2026-10-02)

An opt-in two-pilot prototype now exists in the working tree as
`hipdlp_scaling_selector=true`. See [the continuation report](SELECTOR_STUDY.md)
and `data/2026-10-02-selector/` for guard reproduction, HEAD/HEAD~1 test parity,
new held-out CPU validation and stair's slow dual-mode diagnosis. The new
held-out batch does not confirm the estimated 20% gain; the default stays off.

The subsequent [twenty-large-LP timing study](LARGE20_STUDY.md) compares auto,
Ruiz+PC, the selector and PC-only on CUDA at equal 1e-4 tolerances. PC-only has
the lowest observed aggregate time within a 180-second per-instance budget,
with 9/20 verified solves (Ruiz+PC and selector also 9/20; auto 8/20). A repeat
on the nine PC/Ruiz common solves confirms the order, although almost all of
the difference comes from neso-2245. Eleven all-policy timeouts prevent a
conclusion about total time to converge on every instance. Full copied-input
payload checks and source hashes are retained in `data/2026-10-02-large20/`.

## 5. Next direction: one 200-iteration pilot + a small model (2026-10-03)

User objective: choose a suitable scaling with the simplest effective method
to reduce **aggregate elapsed time**, especially on large LPs. This is the
next experiment, not an implemented policy; keep changes small and preserve
the existing opt-in prototype and results.

- Start with a fixed default (PC-only or Ruiz+PC). PC-only is the current
  large-CUDA timing baseline; Ruiz+PC solved more in the earlier 122-instance
  sample. Neither is established as the universal default.
- Run that default for 200 iterations. Collect cheap static features
  (coefficient/row/column norm distributions, bound types, cost/RHS scales)
  and progress at several pilot checkpoints (unscaled residuals and gap,
  primal/dual growth, fixed-point error, restarts, time per iteration).
  Account for any extra original-model KKT feature evaluation in elapsed time.
- Try a shallow decision tree or another small model with candidates
  P=PC, L=Ruiz+PC, A=CR+Ruiz+PC. Predict time saved relative to the default,
  allowing an uncertain decision to continue the default. Switch only when
  predicted savings justify extra setup/switch costs and regression risk.
  Keeping the default continues its existing state: the first 200 iterations
  are useful work. Cross-scaling warm starts are a separate validation item;
  do not blindly reuse anchors/restart history in different coordinates.
- Keep the CR guard before any CR candidate is run. Kw may be a feature,
  not the sole decision rule: stair's true-solution Kw also prefers the
  slower scaling. One-pilot Kw bias does not rule out learned progress features.
- Optional later fallback: a second pilot only for difficult, uncertain
  instances, with a time budget; avoid paying two-pilot costs on every LP.

Evaluate total time including scaling, features, inference, pilots and
switching, plus solved counts and worst regressions. Train with asymmetric
time-loss penalties and split by instance family; retain held-out names
outside `data/instances.txt`. The nine common large-CUDA solves leave only
about 1.7% headroom over fixed PC even for a free perfect P/L/A selector;
the eleven all-policy timeouts need longer, equal-budget measurements before
being assigned a winning label. Do not treat timeout times as convergence.

Use CPU iteration counts for policy comparisons and serial CUDA timings for
elapsed-time validation. Keep both tolerances equal (currently 1e-4), copy
read-only benchmark LPs to scratch, and **do not change `verifyKkt` without
asking the user**. This handoff update makes no solver/default changes.


## 6. One-pilot small-model experiment completed (2026-10-03)

See [ONE_PILOT_STUDY.md](ONE_PILOT_STUDY.md) and
`data/2026-10-03-onepilot/`. A separate study driver starts P, observes
40/80/120/160/200 iterations and evaluates a depth-two tree. Retaining P
preserves state; L/A switches are cold and charge the discarded 200 iterations
and all setup to the global budgets. Guarded A and final `verifyKkt` acceptance
are preserved. No production default or public learned-selector option is added;
the existing opt-in two-pilot prototype and previous results remain intact.

- Training: 52 instances / 36 name-derived family groups; 13 family-validation
  names; external held-out and large20 families excluded. Only verified P/L/A
  CPU iteration triples are labels; 3x regression penalty and conservative
  support/savings gates. The tree selects only three static features, and a
  static-only ablation produces identical choices.
- CPU: family-validation M reduces total iterations 74.87% against P, but
  fixed A is better and the worst M/P regression is lseu 6.109x. The 26-name
  external batch solves 23/26 for every policy; M reduces common-solve total
  iterations only 3.96% against P. Time-limited qap18/small_pre are censored.
  stair is a separate training-sample diagnostic and still selects slow L.
- Serial CUDA, nine large known solves: P 205.184 s, M 210.887 s (+2.78%),
  both 9/9. Repeat: P 203.148 s, M 210.035 s (+3.39%), both 9/9. The tree
  itself proposes P on every one; no runtime-gate rejection is responsible.
  Feature extraction costs 4.375 s in the primary run, exceeding the free
  measured P/L/A oracle headroom of 3.707 s (1.81%).
- Four independent large heldouts: observing the 200-iteration P pilot costs
  31.331 s versus 23.860 s plain (+31.31%); 0/4 converge at that pilot budget.
  Longer 240 s P/M runs newly solve cont11 at about 233 s; neso-2005 remains
  a timeout for both. L/A are not measured at 240 s, so no winning P/L/A label
  is assigned to either. The remaining previous all-policy timeouts are still
  unresolved.
- CPU/CUDA continuation checks and cold-switch/global-budget driver checks
  pass. Targeted unit tests: 80 assertions; full `[hipdlp]`: 341/350 cases pass,
  with the same inherited 9 failed cases / 16 assertions.

Conclusion: this one-pilot tree does not reduce aggregate large-CUDA time.
Keep it experimental. Next reduce observation cost using existing scaling/
resident-check data, and collect longer equal-budget serial CUDA convergence
labels on hard LPs before learning an elapsed-time policy. Do not infer a
universal default or a timeout winner from this sample. `verifyKkt` remains
unchanged; changing it still requires asking the user.
