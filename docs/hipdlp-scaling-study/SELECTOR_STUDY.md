# Two-pilot scaling selector continuation — 2026-10-02

The opt-in prototype is implemented on `feature/hipdlp-rewrite` as
`hipdlp_scaling_selector=true` (study adapter: `--scaling-selector`). The new
held-out results do **not** validate the previously estimated 20% improvement.
Keep the selector disabled by default. `verifyKkt` is unchanged.

## Implementation and comparison method

A = guarded CR+Ruiz+PC (`auto`), L = Ruiz+PC (`legacy`), S = selector.
The Curtis-Reid inflation guard runs first. If it rejects CR, S uses L directly
and skips both pilots. Otherwise S runs A for 200 iterations and then L for
200 iterations, estimates both candidates' Kw from both pilot solutions,
averages the log estimates, and continues the lower-scoring candidate.
Converged pilots can terminate early. Ties and invalid scores prefer A.

The two solver states live on nested call stacks: A remains suspended while L
runs. Continuing either winner preserves its iterates, anchors, restart state,
PID history, and CUDA command buffers. The discarded pilot counts against the
global iteration limit; all setup and pilots share the time budget. Budgets of
400 iterations or less skip selection. The prototype can retain two complete
backend states and does extra setup and host spectral estimation; iteration
metrics are not wall-time speedup claims.

Scoring matches `tools/proxy2.cpp`: a deterministic 60-step spectral estimate,
initial unscaled weight `(1+||c||2)/(1+||b||2)`, and synthetic slack columns
reconstructed from row activities. The geometric mean is evaluated in log
space. The solver's conservative step-size bounds and final KKT acceptance
are preserved. PID still scores with the *initial* weight, not the subsequent
controller weights.

All performance comparisons below use the deterministic serial CPU backend,
fixed steps, presolve off, and **both** `kkt_tolerance=1e-4` and
`pdlp_optimality_tolerance=1e-4`. CUDA is used separately for neso validation.
All input files were copied from `/srv/mps_da` into `/tmp` before solving.
Manifests record original source hashes; all 26 source hashes in the main
held-out batch were checked again after the runs and were unchanged.

## Build and committed baseline

Hardware: Intel Core i7-14700, NVIDIA GeForce RTX 4070 Ti SUPER. Release,
`ALL_TESTS=ON`, rewrite CUDA enabled. The system-default CUDA 12.0 compiler
cannot build the existing unconditional `cusparseSpMV_preprocess` call; using
the already installed CUDA 12.8 toolkit builds the unmodified commits.

```sh
cmake -S . -B build-scaling-cuda -DCMAKE_BUILD_TYPE=Release \
  -DALL_TESTS=ON -DHIGHS_HIPDLP_CUDA=ON \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.8/bin/nvcc \
  -DCUDAToolkit_ROOT=/usr/local/cuda-12.8
cmake --build build-scaling-cuda -j 12
build-scaling-cuda/bin/unit_tests '[hipdlp]'
```

HEAD `043320530` and its parent `b6dcd91f2` were built separately, with the
parent in a detached `/tmp` worktree. Both unmodified suites fail the same
**9 test cases / 16 assertions**. HEAD has 348 cases; the parent has 347.
Failures are in `TestHiPdlpBackend` (1718, 1884, 3714, 3743, 3744),
`TestHiPdlpTrace` (863), and two legacy CPU characterization cases in `TestPdlpHi`
(970/976 and 1042/1044). Thus this machine reproduces two additional legacy
failures beyond the handout's recorded seven, equally on both commits.

The final prototype suite has 349 cases, 340 passing and the same 9 failing
cases / 16 assertions. Its added selector test passes 32 assertions, including
bitwise equality of winning solution vectors against standalone A and L runs,
the 200-iteration overhead, a non-cadence global limit of 441, and short-budget
fallback. The existing observer-controlled time-limit test also passes.

The committed CUDA baseline logs **“Curtis-Reid rejected”** for the scratch
`neso-2245-filtered` and converges at **106,080 iterations**, with KKT passing,
in 123.1 solver seconds. S skips its pilots and reproduces the same iteration
count, objective and KKT metrics (124.2 solver seconds). This is guard
validation, not a selector speedup. The filtered scratch MPS removes 7,871,092
sub-1e-9 matrix entries; objective coefficients and bounds are retained.

`neso-2005` also rejects CR. Its extended A run reaches the 400,000-iteration
limit without passing KKT (680.2 solver seconds). Equal 2,000-iteration CUDA
A/S runs validate that S skips its pilots and follows the guarded pipeline;
no convergence or speedup claim is made for this instance.

## Held-out CPU results

The 26 main names in `data/heldout-instances.txt` are absent from the original
232-name `data/instances.txt`, including known filename aliases. The batch is
an exploratory sample, with several small LPs and related QAP formulations;
it is not a representative or independent-family performance guarantee.
Budget: 500,000 iterations, with a 600-second wall-time safety limit.

| Policy | KKT-verified solves | Iteration GM, 24 common solves | PAR2 SGM |
|---|---:|---:|---:|
| A | 24/26 | 627.39 | 2,155.67 |
| L | 24/26 | 589.69 | 2,096.43 |
| S | 24/26 | 658.21 | 2,207.65 |

S/A geometric mean iteration ratio is **1.0491**, i.e. **4.91% worse**, on
common solves. The handout-style shifted PAR2 ratio is **1.0165**, also worse.
`forest6` exhausts the same iteration budget under all three policies.
`small_pre` hits the wall-time safety limit under all three; it is included in
solve counts but **censored from all iteration comparisons**, so differing
iteration counts reached at 600 seconds cannot masquerade as policy gains.
PAR2 assigns a fixed-budget failure twice 500,000 iterations; SGM uses +1000.
Ordinary GM floors a zero-iteration solve at one iteration.

Worst common-solve regressions against A:

| Instance | A iterations | S iterations | S/A |
|---|---:|---:|---:|
| afiro | 240 | 400 | 1.667 |
| qap06lp | 400 | 600 | 1.500 |
| sc50a | 440 | 640 | 1.455 |
| sc50b | 480 | 600 | 1.250 |
| qap08lp | 920 | 1,120 | 1.217 |

These are predominantly pilot overhead on cheap or similarly scaled problems.
The larger held-out screen uses `L1_sixm250`, `istefanov_med`, `support19` and
`square41`, each with the same **2,000-iteration** budget for A/L/S. All three
policies solve **0/4**; each run reaches exactly 2,000 iterations. There are no
common solves and thus no solved-iteration GM; PAR2 SGM is 4,000 for each
policy. These results are saved separately. The original optional
20,000-iteration screen was stopped after observing expensive runs; its partial
observations are not mixed into either fixed-budget comparison.

## Why stair regresses

Fresh fixed-step CPU counts reproduce the saved result exactly:

| Policy | Fixed iterations | PID iterations |
|---|---:|---:|
| A | 15,600 | 5,560 |
| L | 287,400 | 29,080 |
| S | 287,600 | 29,280 |

The fixed selector regression is **18.4359x**. Its geometric-mean pilot Kw is
67,626 for A versus 25,054 for L, so it chooses L. Preserving the selected
state is working: S is exactly L plus the discarded 200 iterations.

Two effects explain why the score is misleading:

1. **The short L pilot has not developed its dual solution.** Its own pilot
   dual norm is 2.09 at iteration 200 versus about 37.53 at the simplex
   reference. The dual contribution to L's distance term is about 166 versus
   53,715 at the reference. It remains underestimated at 4,000 iterations
   (dual norm 9.88) and 16,000 iterations (15.84). A's dual norm is already
   close to its reference 4.69 by iteration 4,000.
2. **Kw omits the active problem's directional conditioning.** Even with a
   simplex solution instead of pilots, Kw still incorrectly prefers L:
   approximately 120,120 versus A's 152,812. Full operator condition numbers
   are similar, but on the reference active rows and strictly interior
   columns, the smallest singular value is 0.000454 for L versus 0.003061
   for A; condition numbers are 2,042 versus 310. Including weakly active
   zero-reduced-cost columns gives the same conclusion (2,067 versus 310).
   At iteration 200, **99.53% of L's active dual error** lies in the smallest
   singular direction of that latter operator. This supports a slow dual-mode
   explanation; these local spectral diagnostics do not derive the exact
   nonlinear restart convergence rate.

Longer pilots alone cannot fix the fact that the true-solution Kw ordering is
wrong. PID improves both pipelines and reduces the regression to **5.2662x**,
but the selector still chooses L. Its fixed initial-weight score is especially
limited when PID changes weights substantially; the L run's omega rises to
about 0.233 by iteration 28,800.

The evidence favors leaving this prototype opt-in. A follow-up should test an
independent safeguard sensitive to late dual progress or active conditioning;
calibrating a one-off stair threshold would not establish generalization.

## Reproduction and saved evidence

Data and baseline/unit logs are in `data/2026-10-02-selector/`. Raw full logs,
scratch LPs, commands, reference dumps and pilot solutions remain under
`/tmp/hipdlp-scaling-study`, `/tmp/hipdlp-scaling-heldout-baseline`,
`/tmp/hipdlp-scaling-heldout-selector`, and
`/tmp/hipdlp-scaling-medium-heldout-2000` on this machine. No LPs are added to
the repository.

Build `tools/adapter_x.cpp` against the Release library using the include/link
command in the handout. The following drivers accept paths rather than using
the prior scratch location:

```sh
python3 docs/hipdlp-scaling-study/tools/selector_sweep.py \
  --adapter /tmp/hipdlp-scaling-study/selector-adapter \
  --scratch /tmp/hipdlp-scaling-heldout-repro \
  --names docs/hipdlp-scaling-study/data/heldout-instances.txt \
  --held-out --configs A,L,S --iterations 500000 --time-limit 600 --jobs 6
python3 docs/hipdlp-scaling-study/tools/selector_report.py \
  /tmp/hipdlp-scaling-heldout-repro/results.json
```

Build `tools/stair_diag.cpp` with the same include/link flags, then:

```sh
/tmp/hipdlp-scaling-study/stair-diag \
  /tmp/hipdlp-scaling-study/lp/stair.mps \
  /tmp/hipdlp-scaling-study/quick/stair-reference
OPENBLAS_NUM_THREADS=1 python3 docs/hipdlp-scaling-study/tools/stair_analyze.py \
  --reference-prefix /tmp/hipdlp-scaling-study/quick/stair-reference \
  --trajectories /tmp/hipdlp-scaling-study/stair \
  --output /tmp/hipdlp-scaling-study/stair/analysis.json
```
