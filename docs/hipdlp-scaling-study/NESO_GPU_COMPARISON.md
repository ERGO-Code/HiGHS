# NESO GPU comparison: HPR-LP, cuPDLPx, and rewrite HiPDLP

Experiment dates: 2026-10-03 to 2026-10-04 (Europe/London).
GPU: NVIDIA GeForce RTX 4070 Ti SUPER, 16 GB.

This study runs both NESO instances and separates actual CR use, ingestion differences,
solver stopping criteria, and implementation overhead. The numerical solver sources
and production defaults in hipdlp2 are unchanged. All new files are study tools and
evidence. Official C solver copies under `/tmp/neso-gpu-comparison` have passive
solution-export patches only.

## Inputs and acceptance

The source instances are `/srv/mps_da/neso-{2245,2005}.mps.gz`.
Streaming decompression and SHA-256 checks confirm that both are byte-for-byte
identical to the existing `/tmp/hipdlp-scaling-study/lp/neso-{2245,2005}.mps`
copies used in this study. The provenance checks are retained in
[input-provenance.json](data/2026-10-03-neso-gpu/input-provenance.json).
Canonical scratch inputs explicitly remove matrix entries with `abs(a) <= 1e-9`
and declare variables appearing only in BOUNDS using a zero objective entry.
Objective coefficients, row declarations, and variable/row bounds are preserved.
This coefficient deletion is a model change; results refer to these ingested LPs,
not a proof of equivalence to the unfiltered mathematical models.

`neso_input_audit.cpp` checks that canonical and original inputs, when ingested
by HiGHS with its default coefficient threshold, have exactly identical matrices,
costs, bounds, names, dimensions, objective sense, and offset. Both audits pass.

| Instance | Rows | Columns | Retained matrix entries | Removed entries | BOUNDS-only columns explicitly declared |
|---|---:|---:|---:|---:|---:|
| neso-2245 | 45,715 | 93,151 | 12,980,419 | 7,871,092 | 29,762 |
| neso-2005 | 60,083 | 108,047 | 36,876,017 | 16,761,762 | 29,460 |

cuPDLPx's MPS reader ignores a BOUNDS entry for an unknown column. Its initial
exploratory neso-2245 run therefore loaded only 63,389 columns. That run is excluded
from the common-input comparison. An original-input neso-2005 exploration was
stopped before completion after discovering this discrepancy.
These omitted columns have no matrix or objective entries. Eliminating such
independent columns can preserve the optimal value when their bound intervals
are feasible, but it changes vector dimensions and can alter initialization and
weight/restart trajectories. The comparison therefore keeps explicit declarations.

All compared GPU solves run serially with presolve disabled. The C baselines use
600-second solver budgets; Julia's diagnostic CR screen uses 180 seconds per run.
Rewrite uses PID and `kkt_tolerance = pdlp_optimality_tolerance = 1e-4`.
External native tolerances are recorded per run and are not assumed equivalent.

Acceptance is `hipdlp-original-kkt/2` at `1e-4`, using the same `verifyKkt` routine
as rewrite. Both reconstructed and reported reduced-cost checks are retained when
available. Reported bound duals with full FP64 exports are the primary comparison;
stationarity is checked by that contract. cuPDLPx's default 10-digit text export
can introduce amplified errors at large bounds, so a passive CLI patch exports
17 digits and the reported reduced costs. Older rounded exports are identified
separately and cannot establish that the full-precision native solution fails.

Native solver seconds have different scopes. C HPR's solve time includes analysis,
power estimation and iterations but excludes preceding scaling/setup; rewrite's
adapter solve time includes preprocessing, backend setup and verification. Process
wall time includes input reading and export, and is saved separately. These single
runs are diagnostic evidence, not statistically controlled end-to-end speed ratios.
cuPDLPx measures its native runtime and budget with process CPU `clock()`, whereas
the process wrapper and rewrite use wall time. Its preprocessing field also uses
CPU `clock()` around asynchronous GPU calls and is not an end-to-end setup timer.

## Version-specific scaling findings

GitHub revisions used:

- Julia HPR-LP: `0f8f1501bcc7013b53fec6822e8da91929a39d2e`.
- HPR-LP-C: `79dbaa3fb082a1b826785d13dcd30d9ba350e056`.
- cuPDLPx: `083b90e482f9f0e8229eac965d189a1ba386ca20` (v0.3.0).
- hipdlp2 source: `35b83fd9c0f77dce61f32891be4d7e3d81f1ac9b`.

Julia HPR-LP and rewrite use the same basic CR iteration: alternate row/column
means of negative coefficient logarithms for 20 sweeps, floor coefficients at
`1e-300`, and clamp final multipliers to `[1e-30, 1e30]`. This is coefficient
balancing; it does not optimize active-set conditioning or original KKT progress.
Both then apply 10 Ruiz sweeps and L1 Pock–Chambolle scaling. Rewrite's L2 stage
is a separate optional stage; PC itself is L1 in all three implementations.

Current HPR-LP-C additionally rejects CR if either candidate scale spread exceeds
`1e8`. Its CR-enabled runs reject CR on both canonical NESO models; disabling CR
reproduces the same iterates/iteration counts. Thus those C results do not show
that CR helped NESO. The current Julia source lacks this guard, so its CR switch
is tested separately.

Rewrite auto computes both pipelines and rejects CR if it excessively inflates
scaled costs, right-hand sides or bounds relative to Ruiz+PC. Diagnostic maximum
norm inflation is about `10^3.902` (7,980x) for 2245 and `10^4.832` (67,920x) for
2005. Auto rejects CR on both and reproduces legacy's iteration counts and solution.

Current cuPDLPx defaults to a different 12-sweep geometric scaling, based on
`sqrt(min(abs(a)) * max(abs(a)))` with alternating rows and columns. This is the
Tomlin min/max geometric scaling cited in its source, not the CR log-average
algorithm. Its local older checkout (v0.2.2) did not have this default stage;
the experiments use the pinned current source instead.

Pinned sources:
[Julia CR](https://github.com/PolyU-IOR/HPR-LP/blob/0f8f1501bcc7013b53fec6822e8da91929a39d2e/src/utils.jl#L159),
[C HPR guard](https://github.com/PolyU-IOR/HPR-LP-C/blob/79dbaa3fb082a1b826785d13dcd30d9ba350e056/src/solver/scaling.cu#L237),
[cuPDLPx geometric scaling](https://github.com/MIT-Lu-Lab/cuPDLPx/blob/083b90e482f9f0e8229eac965d189a1ba386ca20/src/preconditioner.cu#L212).

<!-- GPU_RESULTS_BEGIN -->
## GPU results

All KKT entries below use full-precision exports and reported bound duals.
Native seconds retain each implementation's timing scope described above;
Julia process seconds include a fresh process, GPU initialization, JIT and warmup.
External native tolerances were tightened as needed for common acceptance.

| Instance | Solver/configuration | Native tolerance | Iterations | Native seconds | Process wall seconds | Common KKT |
|---|---|---:|---:|---:|---:|---|
| neso-2005 | cuPDLPx, default geometric scaling | 1e-4 | 414,800 | 600.074 | 625.51 | FAIL (time limit) |
| neso-2005 | cuPDLPx, geometric scaling off | 1e-8 | 111,200 | 161.122 | 179.61 | PASS |
| neso-2005 | HPR-C, CR requested/rejected | 1e-8 | 39,750 | 43.679 | 91.58 | PASS |
| neso-2005 | Julia HPR, CR on | 1e-8 | 118,816 | 180.114 | 381.22 | FAIL (time limit) |
| neso-2005 | Julia HPR, CR off | 1e-8 | 86,700 | 132.241 | 329.84 | PASS |
| neso-2005 | Rewrite, auto (CR rejected) | 1e-4 | 34,520 | 124.180 | 183.13 | PASS |
| neso-2005 | Rewrite, Ruiz+PC | 1e-4 | 34,520 | 111.736 | 160.53 | PASS |
| neso-2245 | cuPDLPx, default geometric scaling | 1e-4 | 402,800 | 216.078 | 223.08 | FAIL |
| neso-2245 | cuPDLPx, geometric scaling off | 1e-8 | 10,000 | 5.401 | 12.21 | PASS |
| neso-2245 | HPR-C, CR requested/rejected | 1e-8 | 11,000 | 1.994 | 15.77 | PASS |
| neso-2245 | Julia HPR, CR on | 1e-8 | 324,768 | 180.078 | 363.68 | FAIL (time limit) |
| neso-2245 | Julia HPR, CR off | 1e-8 | 8,200 | 6.042 | 202.86 | PASS |
| neso-2245 | Rewrite, auto (CR rejected) | 1e-4 | 6,880 | 22.817 | 43.31 | PASS |
| neso-2245 | Rewrite, Ruiz+PC | 1e-4 | 6,880 | 18.581 | 35.64 | PASS |

The initial C HPR 2245 runs at native `1e-4` stopped at 1,900 iterations
and failed common KKT (dual residual about 0.063). Their CR-on/off iterates
match because CR is rejected. Tightening to `1e-8` gives the accepted
11,000-iteration result above. C HPR CR-off 2005 at `1e-8` also matches
the accepted CR-requested result (39,750 iterations, 43.705 seconds).

At native `1e-4`, cuPDLPx without geometric scaling stopped at 2,200/5,400
iterations. The complete 2005 export fails common KKT (primal residual
0.000521, dual residual 0.214). The initial 2245 cuPDLPx exports used ten
digits; their derived-cost failures cannot diagnose full-precision solver
accuracy. The default 2245 rerun and tight ablations above use seventeen
digits and reported reduced costs. The full-precision default 2245 rerun still
fails common KKT: dual residual `5.81e-4`, gap `0.02185`, complementarity
`0.02185`, and stationarity residual `0.02326`. Thus text precision alone does
not explain that default run's rejection. Its native `OPTIMAL` label is not
sufficient for the original-coordinate accuracy required here.

Rewrite auto and legacy produce the same iteration counts and solutions
on both instances. Auto spends extra preprocessing work evaluating the CR
candidate it ultimately rejects. This is a performance opportunity in
the selection procedure, not evidence that the rejected CR pipeline
would converge better.

Complete numerical metrics and original argv are in
[verified-results.json](data/2026-10-03-neso-gpu/verified-results.json).
<!-- GPU_RESULTS_END -->

## Where rewrite differs from cuPDLPx

1. **Bounds/cost scaling and weight initialization.** HPR and cuPDLPx additionally
   normalize scaled bounds and costs. cuPDLPx initializes its scaled primal weight
   to one with this option enabled. Rewrite has matrix diagonal scaling and
   initializes its weight from the unscaled cost/RHS norm ratio. These lead to
   different effective primal/dual weights in original coordinates. A scalar
   weight adjustment does not undo a directionally different CR preconditioner.
2. **Conservative step bound.** Rewrite caps the PC-ended matrix's norm bound at
   at least one, yielding base step about `0.998`. cuPDLPx uses its power estimate:
   observed baseline/no-geometric steps are approximately 1.0003/1.0073. This is
   a modest scalar difference; it alone does not establish a cause for large
   iteration regressions.
3. **Restart parameters and cadence.** Rewrite checks restart every 40 iterations
   and uses necessary decay 0.8. Current cuPDLPx uses evaluation frequency 200 and
   necessary decay 0.5. Both default PID gains are P=0.99, I=0.01, D=0 and smoothing
   0.3, but different coordinates and restart times produce different histories.
4. **Model formulation.** Rewrite converts ranged rows to boxed slack variables;
   current HPR/cuPDLPx project directly against row lower/upper bounds. CR can
   therefore see a different operator when ranged constraints are present.
5. **GPU work and host verification.** C HPR has adaptive matrix reduction,
   specialized update backends, and deferred empty-row updates. Its tight 2245
   run continues on a much smaller set of free columns. Rewrite's actual phase
   timings are measured separately below; raw iteration count alone cannot
   identify a kernel performance difference.
6. **Stopping criteria.** HPR/cuPDLPx native residuals normalize L2 norms by global
   bound/cost norms. The common verifier measures individual row/variable
   violations, dual signs, gap, stationarity and complementarity. HPR's 2245
   `1e-4` stop at 1,900 iterations fails the common contract even with reported
   bound duals; its `1e-8` stop at 11,000 passes.

Useful rewrite sources: `model/scaling.cc`, `model/problem.cc`,
`algorithm/step_size_policy.cc`, `algorithm/restart_policy.cc`,
`algorithm/solver.cc`, and `verification/kkt_verifier.cc`, under
`highs/pdlp/hipdlp/`.

## Rewrite phase profile (2245, Ruiz+PC/PID)

An additional diagnostic run reproduces 6,880 iterations and passes common KKT.
It reports 17.915 solver wall seconds, including the following nested phases:

| Phase | Seconds |
|---|---:|
| Model preprocessing, including scaling | 11.027 |
| Backend setup | 0.632 |
| Operator norm estimation | 0.023 |
| Iteration loop | 6.141 |
| Check boundaries, inside the loop | 5.309 |
| Original-model host KKT verification, inside boundaries | 2.179 |
| Fixed-point fallback, nested | 0.851 |
| Candidate vector downloads, nested | 0.010 |

There are 141 check boundaries, 7 original KKT verifications (6 rejected),
7 candidate downloads, 21 restarts/weight updates, and 140 ordinary check
synchronizations. Ordinary check D2H traffic is only 15,680 bytes in total;
the record does not show a full-vector CPU download on every iteration.
The phase counters overlap and must not be added as independent costs.
Asynchronous GPU work can be charged to a synchronization boundary, so subtracting
boundary time from iteration time does not measure pure GPU kernel cost.

This identifies model preprocessing and the check/verification path as concrete
rewrite optimization targets. The native timing-scope difference also explains
part of the apparent C HPR/rewrite solver-seconds gap. C HPR's reduced matrix
contains about 18,800 free columns late in its accepted 2245 run, versus 93,151
columns in the complete model. A remaining kernel-speed difference cannot be
quantified from these wall phases alone. No analogous phase breakdown was run
for 2005.

The raw profile is retained as
[profile.json](data/2026-10-03-neso-gpu/runs/neso-2245-canonical/rewrite-profile/profile.json).

## Interpretation and next implementation work

The evidence supports keeping the CR guard for these NESO models. Julia HPR's
unguarded CR runs hit their 180-second main-solve limits on both instances
(2245: native KKT error 0.336; 2005: approximately 1). C HPR's enabled CR is
rejected before scaling. Thus there is no contradiction between CR being
present in HPR's source and CR being harmful here. This is evidence for these
models and parameter settings, not a claim that CR always hurts LPs.

The likely mechanism is that coefficient balancing does not account for costs,
bounds, or the submatrix controlling the active constraints. On these instances
CR substantially worsens the transformed data magnitudes despite subsequent
Ruiz/PC equilibration. This mechanism is an inference supported by the scaling
diagnostics and ablations; it has not been isolated into a proof of which
coordinate dominates convergence.

For rewrite/cuPDLPx parity, first benchmark Ruiz+PC without the extra geometric
stage and require the same original-model KKT contract. Next isolate bounds/cost
normalization together with its induced primal weight, then restart cadence and
decay. Measure changes separately on both instances. Implementing HPR-C-style
adaptive matrix reduction is a separate performance project. A GPU-resident
original-model verification path is worthwhile only to the extent shown by the
phase profile. These changes are proposed follow-up experiments; none have been
applied to the production solver in this study.

## Reproduction commands

Scratch root: `/tmp/neso-gpu-comparison`. Canonicalization and audit commands:

```sh
python3 docs/hipdlp-scaling-study/tools/neso_canonicalize.py \
  /tmp/hipdlp-scaling-study/lp /tmp/neso-gpu-comparison/canonical-checked
/tmp/neso-gpu-comparison/input-audit SOURCE.mps CANONICAL.mps
```

Build all C programs with existing CUDA 12.8, Release, without changing solver
sources. HPR-C uses `CMAKE_CUDA_ARCHITECTURES=89`, `BUILD_GPU_PRESOLVER=OFF`,
`BUILD_SHARED_LIB=OFF`; cuPDLPx uses its required PSLP v0.0.11 and
`CUPDLPX_BUILD_SHARED_LIB=OFF`. The Julia project is instantiated in its scratch
checkout; warmup is enabled for actual solves.

```sh
python3 docs/hipdlp-scaling-study/tools/neso_gpu_compare.py \
  --scratch /tmp/neso-gpu-comparison \
  --models /tmp/neso-gpu-comparison/canonical-checked \
  --names neso-2245-canonical,neso-2005-canonical \
  --solvers hpr-cr,hpr-no-cr,rewrite-auto,rewrite-legacy,cupdlpx
python3 docs/hipdlp-scaling-study/tools/neso_verify_runs.py \
  --scratch /tmp/neso-gpu-comparison \
  --verifier /home/yzhou/Documents/Github/hipdlp2/build-scaling-cuda/bin/external_lp_solution_verifier \
  --reported-verifier /tmp/neso-gpu-comparison/neso-kkt
```

Use `--hpr-tolerance 1e-8` for tight HPR checks;
`--solvers julia-cr,julia-no-cr --hpr-tolerance 1e-8 --time-limit 180` for Julia;
`--solvers cupdlpx-no-gm-tight --cup-tolerance 1e-8` for the tight cuPDLPx ablation.
Each run retains the exact argv, model hash, process result, native summary,
solution vectors, logs and independent verifier output. LPs and full vectors
remain in scratch, not in the repository.
