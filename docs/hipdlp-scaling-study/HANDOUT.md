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
