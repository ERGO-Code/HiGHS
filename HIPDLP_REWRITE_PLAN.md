# HiPDLP Rewrite Implementation Plan

## 1. Objective

Rewrite HiPDLP as a modular, high-performance solver that remains integrated
with HiGHS while separating:

- HiGHS model and option handling;
- problem transformation and scaling;
- the PDHG/Halpern algorithm;
- restart, step-size, and convergence policies;
- linear algebra and device-memory management;
- deterministic serial CPU and NVIDIA CUDA implementations;
- explicit unavailable placeholders for OpenMP and AMD ROCm;
- correctness verification, profiling, and benchmarking.

The first complete rewrite will preserve the current Halpern-PDHG behavior and
the public `solver = "hipdlp"` contract. It will provide deterministic reference
CPU and CUDA backends. OpenMP and ROCm remain explicit unavailable placeholders
in the runtime factory and build configuration during this effort. Low-precision
GPU SpMV and multistart warm-up are follow-on research projects, but the initial
interfaces must accommodate them without another solver rewrite.

The rewritten path uses deterministic Curtis-Reid scaling by default before
any enabled legacy Ruiz, Pock-Chambolle, or L2 stages. The implementation is
based on the algorithm and behavior in the MIT-licensed HPR-LP repositories,
with pinned revisions and independent HiGHS-native tests rather than a runtime
dependency.

## 2. Confirmed Decisions

| Area | Decision |
| --- | --- |
| Core ownership | Keep the solver coupled to HiGHS and its model/types |
| Backend selection | Conditionally compile any subset of backends and select one through a runtime factory |
| Initial backends | Deterministic serial CPU and CUDA; OpenMP and ROCm are unavailable placeholders |
| Initial algorithm | Reproduce current Halpern-PDHG behavior, then expose research policies |
| Default scaling | Curtis-Reid first, then enabled legacy scaling stages; deterministic host computation for every backend |
| Compatibility | Preserve solver name, relevant options, statuses, solution fields, logging, and iteration reporting |
| Migration | Run legacy and rewritten implementations side by side until cutover |
| Correctness | Independently recomputed KKT metrics plus consensus with in-tree cuPDLP, external cuPDLPx, cuOpt, and HiGHS reference solvers |
| Accuracy target | Relative primal residual, dual residual, and objective gap at or below `1e-4` |
| Performance | Rewritten CPU and CUDA geomean solver time must be within 5% of their frozen legacy HiPDLP baselines |
| CPU optimization | Defer OpenMP and architecture tuning; establish deterministic owned serial kernels first |
| Build systems | Full backend support in CMake; keep Bazel and Meson CPU builds working |
| AMD scope | Factory/build placeholder only; implementation and compile validation are deferred |
| Benchmark cadence | Smoke tests per change, representative Netlib/energy per PR, full Netlib/energy nightly, curated Mittelmann/custom release runs |

## 3. Current-State Problems to Remove

The current implementation has useful algorithmic work, but its boundaries
prevent safe extension:

- `highs/pdlp/hipdlp/pdhg.cc` contains preprocessing, postprocessing, algorithm
  control, convergence, CPU operations, CUDA allocation, CUDA execution, and
  reporting in one file.
- `highs/pdlp/hipdlp/pdhg.hpp` exposes CUDA library types and owns all host and
  device buffers directly in `PDLPSolver`.
- `CUPDLP_GPU` controls both in-tree cuPDLP and HiPDLP GPU behavior, preventing
  independent backend configuration.
- GPU setup and teardown use raw handles and pointers rather than RAII.
- CPU and GPU paths duplicate algorithm decisions, creating parity risk.
- `linalg.*`, `scaling.*`, and `restart.*` operate directly on mutable global
  solver state or `HighsLp`, limiting isolated tests.
- The current `check/TestPdlpHi.cpp` covers one feasible model and does not
  independently verify KKT conditions.
- Several error paths assert or terminate the process instead of returning a
  structured solver error.
- Current timing does not cleanly separate parsing, transformation, transfer,
  solve, verification, and end-to-end costs.
- The global default `pdlp_step_size_strategy = adaptive` is unsupported by
  legacy HiPDLP and is silently converted to PID. On `scrs8` with presolve off
  and `kkt_tolerance = 1e-4`, this path reached 10,381,800 iterations and the
  60-second limit with a relative gap above `1e10`, while CPU cuPDLP-C solved
  the same model in 0.31 seconds and 39,840 iterations.
- Forcing legacy HiPDLP to fixed step size makes `scrs8` terminate in 0.15
  seconds and 21,040 iterations, but it reports `OPTIMAL` with relative primal
  infeasibility `1.76e-4` and objective `907.714` versus the trusted objective
  near `904.297`. This is a stopping/validation defect, not an acceptable
  resolution of the PID failure.

The rewrite must remove these structural issues rather than moving the same
coupling into smaller files.

## 4. Target Architecture

### 4.1 Proposed source layout

```text
highs/pdlp/
  HiPdlpWrapper.{h,cpp}             # Stable HiGHS entry point
  hipdlp/
    api/
      options.{hpp,cc}              # Validated HiGHS option translation
      status.{hpp,cc}               # Internal-to-HiGHS status mapping
    model/
      problem.{hpp,cc}              # Transformed/scaled HiGHS LP data
      transform.{hpp,cc}            # Preprocess/postprocess mapping
      scaling.{hpp,cc}              # Scaling algorithms and scale metadata
    algorithm/
      solver.{hpp,cc}               # Backend-independent solve loop
      state.hpp                     # Iterate and scalar algorithm state
      halpern_pdhg.{hpp,cc}         # Initial algorithm implementation
      convergence.{hpp,cc}          # Residuals, gap, stopping decisions
      restart_policy.{hpp,cc}       # Restart interface and Halpern policy
      step_size_policy.{hpp,cc}     # Step-size interface and initial policy
    backend/
      backend.hpp                   # Coarse-grained backend contract
      backend_factory.{hpp,cc}      # Availability and runtime selection
      types.hpp                     # Opaque buffers, matrix views, capabilities
      cpu/
        cpu_backend.{hpp,cc}        # Deterministic serial FP64 backend
        sparse_kernels.{hpp,cc}
        vector_kernels.{hpp,cc}
      cuda/
        cuda_backend.{hpp,cu}
        cuda_kernels.{cuh,cu}
      # OpenMP and ROCm are factory capabilities only in this effort.
    verification/
      kkt_verifier.{hpp,cc}         # Independent verification on original LP
    diagnostics/
      iteration_log.{hpp,cc}
      timings.{hpp,cc}
      profile_events.hpp
    legacy/                         # Temporary migration-only adapter
```

Do not move legacy files into `legacy/` until the new path builds. Initially,
leave them unchanged and add an adapter that gives the wrapper an explicit
legacy/rewrite choice. Delete the legacy implementation only after all cutover
gates pass.

### 4.2 Dependency direction

Dependencies must flow in one direction:

```text
HiPdlpWrapper
  -> options + transform
  -> algorithm
  -> backend contract
  -> selected backend implementation
```

The algorithm may depend on backend-neutral handles and result structures. It
must not include CUDA, HIP, cuBLAS, cuSPARSE, rocBLAS, or rocSPARSE headers.
Backend implementations may depend on algorithm data contracts, but may not
make restart or termination decisions.

The independent KKT verifier consumes the original `HighsLp` and final
unscaled solution. It must not call the backend's convergence implementation.

### 4.3 Problem representation

Create a solver-owned `Problem` with:

- dimensions, objective sense, offset, costs, and variable bounds;
- normalized row bounds and equality/inequality metadata;
- CSR data for `A`;
- transpose data built once in the format preferred by the backend;
- preprocessing maps needed to recover the original solution;
- row and column scale vectors;
- unscaled norms used by stopping criteria;
- explicit index-width metadata.

The representation remains based on HiGHS scalar and index types. Transformation
must not mutate the caller's `HighsLp`. Construct the transformed problem once,
then upload or bind it to the selected backend once.

Empty rows/columns, fixed variables, free variables, boxed rows, infinite
bounds, objective offsets, maximization, and 32/64-bit vendor sparse-index
constraints must have explicit behavior and tests.

### 4.3.1 Default Curtis-Reid scaling

The rewrite adds a HiPDLP-specific scaling selector:

```text
hipdlp_scaling = auto | off | curtis_reid | legacy
hipdlp_curtis_reid_iterations = <positive integer, default 20>
```

Selection rules:

- `auto` is the rewrite default. It applies Curtis-Reid first, followed by the
  Ruiz, Pock-Chambolle, and L2 stages selected by the existing
  `pdlp_scaling_mode` bits;
- `curtis_reid` applies only Curtis-Reid scaling;
- `legacy` omits Curtis-Reid and retains the existing HiPDLP scaling sequence;
- `off`, or the existing `pdlp_features_off` scaling bit, disables all scaling;
- the legacy implementation keeps its current defaults during migration so
  differential tests can distinguish algorithm changes from scaling changes.

Curtis-Reid scaling is computed once on the solver-owned transformed problem,
after row normalization/slack construction and before backend binding or
upload. The same deterministic host implementation supplies scale vectors to
the serial CPU and CUDA backends; CUDA must not independently recompute
different scales.

For each stored nonzero $a_{ij}$, let $r_i$ and $c_j$ be row and column log
scales. Starting from zero, one sweep computes

$$
r_i = \frac{1}{|N(i)|}\sum_{j \in N(i)}\left(-\log|a_{ij}| - c_j\right),
\qquad
c_j = \frac{1}{|N(j)|}\sum_{i \in N(j)}\left(-\log|a_{ij}| - r_i\right).
$$

After the configured number of sweeps, define
$D_r = \operatorname{diag}(\exp(r))$ and
$D_c = \operatorname{diag}(\exp(c))$. Apply

$$
A_s = D_r A D_c,\quad
\ell^r_s = D_r\ell^r,\quad u^r_s = D_r u^r,\quad
c_s = D_c c,\quad
\ell^x_s = D_c^{-1}\ell^x,\quad u^x_s = D_c^{-1}u^x.
$$

Postsolve composes all enabled scale stages and recovers
$x = D_c x_s$, $y = D_r y_s$, and bound reduced costs
$z = D_c^{-1}z_s$ before applying row permutations and objective-sense
mapping.

Use the HPR-LP safeguards as minimum requirements: treat empty rows and columns
as scale one, floor $|a_{ij}|$ at `1e-300` before taking a logarithm, clamp
final multiplicative scales to `[1e-30, 1e30]`, reject non-finite input or
intermediate values with a structured error, and never change the sparse
pattern. Accumulate in a fixed row/column order for reproducibility.

Behavioral references are pinned to:

- `PolyU-IOR/HPR-LP` commit
  `0f8f1501bcc7013b53fec6822e8da91929a39d2e`, especially
  `src/utils.jl`;
- `PolyU-IOR/HPR-LP-C` commit
  `358295ca9af3a9f1413174f2f63e5bdf3032c548`, especially
  `src/scaling.cu`.

Both repositories are MIT licensed. Implement the algorithm in HiGHS-native
C++ and record source provenance in `THIRD_PARTY_NOTICES.md` if implementation
code is derived rather than independently written from the algorithm.

### 4.4 Backend contract

Use a small runtime-selected interface with coarse and fused operations. Do not
place virtual dispatch inside element loops. The contract should cover:

1. capability and device queries;
2. problem upload/binding and workspace creation;
3. vector initialization and host/device transfer;
4. `A*x` and `A^T*y`;
5. fused projected primal and dual updates;
6. fused Halpern/average updates;
7. dot products, norms, differences, and scalar reductions;
8. convergence-metric reduction into a small host result;
9. synchronization and timing boundaries;
10. solution download.

Backend-owned `ProblemHandle`, `VectorHandle`, and `Workspace` objects must be
move-only RAII types. Device APIs and raw pointers stay private to the backend.
All allocation is completed before entering the steady-state iteration loop.
Per-iteration heap allocation and implicit host/device transfer are forbidden.

The interface must carry a precision descriptor with separate storage,
compute, and reduction precision. The initial implementation supports only
FP64 for all three, but this boundary enables later low-precision SpMV with
FP64 accumulation or correction.

Reserve a capability for batched vectors and SpMM, but do not implement a
multistart algorithm in the rewrite. A backend must reject unsupported
capabilities explicitly rather than silently falling back.

### 4.5 Backend selection

Add a HiPDLP-specific option equivalent to:

```text
hipdlp_backend = auto | cpu | cuda | openmp | rocm
```

Rules:

- `auto` selects CUDA only when that backend was built and a compatible device
  is available; otherwise it selects deterministic serial CPU.
- `cpu` is deterministic, single-threaded, and intended for correctness,
  diagnosis, and the initial CPU implementation.
- `openmp` and `rocm` are reserved names that return a clear unavailable error
  in this effort; `auto` never selects them;
- requesting an unavailable backend returns a clear `HighsStatus::kError`;
  there is no silent fallback for an explicit request;
- backend availability is reported once in the startup log;
- the runtime factory contains selection only, not algorithm logic.

Use an independent CMake option instead of `CUPDLP_GPU`:

```text
HIGHS_HIPDLP_CUDA
```

Do not add empty OpenMP or ROCm build options that could imply support. Exact
future names should follow the prevailing HiGHS option convention when those
backends are implemented. In-tree cuPDLP configuration must remain independent.

### 4.6 Algorithm policies

Keep one backend-independent solve loop. Separate these replaceable policies:

- `RestartPolicy`;
- `StepSizePolicy`;
- `ConvergencePolicy`;
- `AveragingPolicy`.

The first policy set reproduces the current Halpern restart, fixed/PID behavior
that is reachable through HiGHS options, scaling sequence, iterate selection,
and convergence cadence. Unsupported legacy option combinations must continue
to produce a clear warning or error consistent with existing behavior. The
rewrite must not silently reinterpret the global adaptive default as PID:
option translation must select a documented supported policy or reject the
request, and startup logging must report the effective policy exactly once.

Policies own their parameters and scalar history. The solver owns iterate
handles and orchestrates calls. A policy may request backend reductions but
must not access device pointers.

## 5. Backend-Specific Requirements

### 5.1 Serial CPU

- Use straightforward single-threaded FP64 loops.
- Favor auditability, deterministic reduction order, and bounds checks in
  debug builds over speed.
- Provide the oracle for backend primitive tests and short iteration traces.
- Avoid external linear algebra dependencies.

### 5.2 Deferred OpenMP placeholder

- Reserve the `openmp` backend name and capability identifier.
- Return a structured unavailable-backend error for explicit requests.
- Keep the serial CPU matrix representation and coarse backend contract usable
  by a future OpenMP implementation without adding threading in this effort.
- Do not claim optimized CPU support or performance.

### 5.3 CUDA

- Encapsulate CUDA stream, cuBLAS handle, cuSPARSE handle, descriptors, graph
  objects, events, buffers, and error translation in RAII wrappers.
- Upload matrix and bounds once.
- Reuse SpMV descriptors and external buffers.
- Keep iterates and scalar state on device between convergence checks.
- Fuse projection, update, averaging, and simple reductions where profiling
  shows launch or bandwidth savings.
- Synchronize only for required host decisions and trustworthy timings.
- Treat CUDA Graph capture as a backend optimization, not algorithm structure.
- Select sparse algorithms by capability and benchmark evidence; record the
  selected algorithm in diagnostics.
- Return structured errors to HiGHS. Library failures must not call `exit`.

### 5.4 Deferred ROCm placeholder

- Reserve the `rocm` backend name and capability identifier.
- Return a structured unavailable-backend error for explicit requests.
- Keep backend semantics and result types suitable for a future
  HIP/rocBLAS/rocSPARSE implementation without adding ROCm source or toolchain
  dependencies in this effort.
- Do not claim ROCm compilation, hardware correctness, or performance support.

## 6. Correctness Strategy

### 6.1 Independent KKT verifier

Implement one backend-independent verifier against the original, unscaled
`HighsLp`. Given primal and dual values, independently recompute:

- row activities;
- variable- and row-bound violations;
- reduced costs;
- dual feasibility under HiGHS sign conventions;
- primal and dual objectives including sense and offset;
- relative primal residual;
- relative dual residual;
- relative objective gap;
- complementarity where defined;
- finiteness of every input and result.

The verifier returns raw and normalized metrics plus the worst offending row or
column. Solver-reported `OPTIMAL` is accepted only if independent metrics pass.

Primary optimality gates:

```text
relative primal residual <= 1e-4
relative dual residual   <= 1e-4
relative objective gap   <= 1e-4
```

Reference-objective relative error should retain the existing benchmark-harness
gate of `1e-3` where a trusted reference objective exists. Near-zero objectives
must use a mixed absolute/relative comparison.

### 6.2 Test layers

1. **Primitive unit tests**
  Compare every backend operation against serial CPU on generated matrices
   covering empty, diagonal, rectangular, skewed, duplicate-entry-normalized,
   and large-index cases.

2. **Transformation tests**
   Round-trip preprocessing, scaling, unscaling, row splitting, objective sense,
  and postsolve maps without invoking PDHG. Curtis-Reid tests must cover fixed
  expected row/column scales, 20-sweep determinism, empty rows/columns,
  extreme coefficient ranges, composed scaling order, and primal/dual/reduced-
  cost round trips. Compare fixtures with both pinned HPR-LP implementations.

3. **Policy tests**
   Feed synthetic metric sequences into restart, step-size, averaging, and
   stopping policies.

4. **Iteration differential tests**
   Run a fixed number of iterations from identical initial state and compare
   iterates, products, step sizes, restart decisions, and metrics between
  serial CPU and CUDA. Use tolerant comparisons for
   parallel reductions; require identical decisions away from thresholds.

5. **HiGHS integration tests**
   Cover `solver = "hipdlp"`, backend selection, relevant PDLP options, model
   statuses, iteration/time limits, solution validity flags, logging disabled,
  and unavailable-backend errors. Include `scrs8` regressions for the default
  step-size translation and the fixed-step false-optimal legacy case.

6. **Problem-class tests**
   Include feasible, infeasible, unbounded, boxed-row, free-variable,
   fixed-variable, maximization, objective-offset, empty-row/column, and
   numerically scaled LPs.

7. **Cross-solver tests**
   Run the same original MPS and tolerances through rewritten HiPDLP, legacy
   HiPDLP, in-tree cuPDLP, frozen external cuPDLPx, cuOpt, and HiGHS
   simplex/IPM references.

Baseline disagreement is not automatically a rewrite failure. The independent
KKT verifier and trusted high-accuracy HiGHS reference run arbitrate objective
or status disagreements, and the discrepancy is retained as an artifact.

### 6.3 External baselines

- Freeze external cuPDLPx at
  `/home/yanyu/github/cuPDLPx/build-baseline/cupdlpx`, recording executable
  hash, source revision, build cache, CUDA version, and command.
- Obtain cuOpt from `https://github.com/NVIDIA/cuopt`, pin a commit or release,
  and record whether the comparison uses its CLI, C API, Python package, or
  container. Do not benchmark an unpinned `main`.
- Pin in-tree cuPDLP and legacy HiPDLP to the same repository commit used to
  begin the rewrite.
- Pin the HPR-LP and HPR-LP-C revisions listed in Section 4.3.1 as
  Curtis-Reid behavioral references, recording source hashes, license, fixture
  generation commands, and expected scaling vectors.
- Run high-accuracy HiGHS simplex/IPM to generate or confirm reference status
  and objective independently of first-order solvers.

## 7. Benchmark and Performance Plan

### 7.1 Existing assets

Use:

- problem data under `/home/yanyu/benchmark-data/lp`;
- environments under `/home/yanyu/benchmark-envs`;
- immutable result directories under `/home/yanyu/benchmark-results/lp`;
- the existing harness in `/home/yanyu/github/lp-benchmarks`;
- existing Netlib feasible and energy continuous manifests and reference data.

Extend the harness with solver configurations and parsers for:

- legacy HiPDLP serial CPU;
- rewritten HiPDLP serial CPU;
- legacy HiPDLP CUDA;
- rewritten HiPDLP CUDA;
- in-tree cuPDLP;
- frozen external cuPDLPx;
- pinned cuOpt;
- HiGHS simplex/IPM reference runs.

All runs use the same input file, presolve policy, tolerances, time limit,
thread count, and objective/status verifier unless a solver cannot express an
option. Every exception is recorded in metadata.

### 7.2 Timing dimensions

Record separately:

- model read/parse;
- preprocessing and scaling;
- backend initialization and allocation;
- host-to-device upload;
- power-method/initialization work;
- steady-state iterations;
- convergence checks and synchronization;
- device-to-host download;
- unscale/postsolve;
- independent verification;
- solver total;
- end-to-end wall time.

The performance gate uses solver time with a documented boundary and reports
end-to-end time alongside it. Never compare a kernel-only time with another
solver's end-to-end time.

### 7.3 Experimental protocol

- Build release binaries with symbols retained for profiling.
- Freeze compiler, flags, dependency versions, GPU clocks/power mode where
  controllable, thread affinity, and environment.
- Record CPU, GPU, driver, CUDA, compiler, executable hash, git revision,
  manifest hash, solver options, and environment variables.
- Run at least one warm-up and five measured repetitions for performance gates.
- Randomize solver/instance order with a recorded seed.
- Use median per-instance time, then compare geometric means of paired ratios.
- Report bootstrap confidence intervals for aggregate ratios.
- Separate cold-start and warm-start results.
- Check GPU completion with events or explicit synchronization before stopping
  timers.
- Preserve stdout/stderr, parsed records, verifier output, and metadata for
  every run.
- Never overwrite a result directory.

### 7.4 Performance gates

For each backend, compare the rewrite only with its frozen legacy HiPDLP
counterpart:

```text
geomean(rewrite solver time / legacy solver time) <= 1.05
```

Also flag any instance regression above 20% for investigation. It may be waived
only with an artifact-backed explanation such as corrected behavior, different
termination iteration, or noise demonstrated by additional repetitions.

cuPDLP, cuPDLPx, cuOpt, and HiGHS reference solver timings are comparative
reports, not the 5% cutover denominator.

Benchmark tiers:

| Tier | Cadence | Scope |
| --- | --- | --- |
| Smoke | Every relevant change | Tiny representative LPs on serial CPU and available CUDA |
| PR | Before merge | Curated Netlib and energy set spanning tiny/large, tall/wide, and sparse/dense |
| Nightly | Scheduled | Full Netlib feasible and energy continuous manifests |
| Release | Before cutover/release | Nightly suites plus curated Mittelmann and custom problems |

Before using Mittelmann/custom in a gate, create versioned manifests with
checksums, dimensions, roles, and trusted reference objectives/statuses.

## 8. Phased Implementation

### Phase 0: Freeze baselines and acceptance artifacts

- Freeze the current git revision and build configurations.
- Build legacy serial CPU and CUDA variants.
- Record current unit/integration outcomes, including known failures.
- Add benchmark-harness solver definitions for legacy HiPDLP and in-tree
  cuPDLP.
- Pin external cuPDLPx and cuOpt.
- Run a small correctness/performance baseline, then full Netlib and energy
  baselines.
- Define the exact solver-time boundary and result schema.

**Exit criteria:** reproducible baseline commands, executable hashes, parsed
results, independent reference objectives/statuses, and immutable result
directories exist.

**Execution record (2026-08-25):** completed on branch
`feature/hipdlp-rewrite` from base commit
`e29a6c04bc28a0b5ebcfa7b42b53f221d4f3d73d`. Immutable artifacts are under
`/home/yanyu/benchmark-results/lp/hipdlp-rewrite/`
`m0-legacy-serial-e29a6c0-20260825`. The tagged PDLP tests pass (9 cases,
35 assertions). In the 13-instance, 60-second legacy HiPDLP baseline, 10
instances report `OPTIMAL`, while `scrs8`, `sctest`, and `stair` reach the time
limit. Five of the 10 optimal claims fail the strict `1e-4` final KKT gate.
High-accuracy simplex references solve all 13 instances. These failures are
frozen characterization results and entry criteria for Phase 1, not baseline
setup failures.

### Phase 1: Add verifier and characterization tests

- Implement the independent KKT verifier before changing algorithm code.
- Expand HiPDLP integration tests around current public behavior.
- Add transformation and edge-case characterization tests.
- Capture short legacy iteration traces for representative models.
- Inventory which current options are honored, ignored, or forced.

**Exit criteria:** current behavior and known defects are documented by tests;
the verifier can accept/reject solutions without using HiPDLP internals.

### Phase 2: Create architecture seams and side-by-side selection

- Add new directories, result/error types, validated option translation, and
  the backend contract.
- Add an internal `legacy|rewrite` selector used only during migration.
- Keep public `solver = "hipdlp"` behavior unchanged by default.
- Separate transformed problem ownership from the legacy `PDLPSolver`.
- Add CMake feature detection and backend registration without enabling new
  implementations.

**Exit criteria:** legacy tests still pass through the new wrapper boundary;
an empty rewrite path can be selected and reports a structured not-ready error.

### Phase 3: Implement serial CPU and algorithm parity

- Implement reference sparse/vector primitives and RAII host buffers.
- Implement deterministic Curtis-Reid scaling and its scale metadata first,
  then port legacy Ruiz/Pock-Chambolle/L2 composition, preprocessing,
  Halpern-PDHG, restart, step-size, convergence, unscale, and postsolve behavior
  into the new modules.
- Remove algorithm decisions from backend calls.
- Differential-test primitives, iteration traces, restart decisions, final
  metrics, status, and solution against legacy CPU.
- Run algorithm-parity differentials with `hipdlp_scaling=legacy`; validate the
  rewrite default separately with `hipdlp_scaling=auto` and the independent KKT
  verifier.
- Fix defects found during porting only when the verifier demonstrates the
  correction; record intentional behavior changes.

**Exit criteria:** serial CPU passes unit/integration tests and correctness
gates on the PR suite, with explained parity differences.

**Bounded Phase 3B execution record (2026-08-26):** deterministic rewrite-side
Ruiz, Pock-Chambolle, and L2 scaling now compose after Curtis-Reid through one
cumulative multiplicative row/column scale convention. `auto`, `off`,
`curtis_reid`, `legacy`, and the existing scaling-disable bit are covered by
fixed-vector, order, determinism, sparse-pattern, extreme-input, transactional,
and full scale/unscale/postsolve tests. Non-iterative postsolve recovers original
primal columns, row activities, row duals, reduced costs, objective sense and
offset, row permutation/sign corrections, and removes synthetic slack columns.
The public solver remains on the frozen legacy route, and the internal rewrite
route still returns structured `kNotReady` without a solution or backend bind.
The fresh CPU-only build is
`/home/yanyu/builds/hipdlp2-m3b-cpu-20260826`: `[phase3b]` passes 558
assertions in 7 cases, `[phase3]` 878 in 17, `[migration]` 122 in 1,
`[hipdlp]` 1271 in 22, and `[pdlp]` 1301 in 30. Its compile commands and linked
libraries contain no CUDA or ROCm dependency. Meson, Bazel/Bazelisk, and
buildifier were unavailable, so their registrations were validated statically.
The PDHG loop, restart, step-size, averaging, convergence, and stopping-policy
work remains deferred to the next bounded Phase 3 slice.

**Bounded Phase 3C execution record (2026-08-26):** an internal fixed-count
trace now builds and scales the solver-owned `Problem` through the Phase 3B
pipeline, binds only the deterministic serial FP64 CPU backend, preallocates
all handles and trace storage, and executes the legacy-equivalent projected,
reflected Halpern-PDHG update with explicit positive finite primal and dual
step sizes. It records deterministic primal, dual, `Ax`, and `A^T y` values
without allocation or transfer calls inside the iteration loop. Exact zero-,
one-, two-, and four-step fixtures, equality/inequality and column-bound
projections, anchor weights, bitwise repeated runs, pre-restart legacy
differentials, all scaling selections, transactional errors, and loop resource
auditing pass in `[phase3c]` (286 assertions in 8 cases). `[phase3b]` passes 643
assertions in 8 cases, `[phase3]` 1164 in 25, `[migration]` 122 in 1,
`[hipdlp]` 1557 in 30, and `[pdlp]` 1587 in 38. Public `solver = "hipdlp"`
still routes to legacy and rewrite `algorithm::solve` still returns structured
`kNotReady` without a valid solution. The fresh CPU-only build is
`/home/yanyu/builds/hipdlp2-m3c-cpu-20260826`; an unusable `CUDACXX` was set,
both GPU options remained off, 304 Ninja commands and 298 compile-database
entries contained no GPU compilation, and neither `unit_tests` nor
`libhighs.so` linked CUDA, HIP, or ROCm libraries. Meson, Bazel/Bazelisk, and
buildifier remained unavailable, so registrations were validated statically.
Restart, adaptive/PID/Malitsky-Pock step sizes, convergence, stopping, result
acceptance, public status mapping, and rewrite routing remain deferred.

### Phase 4: Reserve optimized CPU extension

- Keep `openmp` as an explicit unavailable backend capability.
- Preserve owned matrix layouts and coarse backend operations needed by a
  future OpenMP implementation.
- Do not add OpenMP pragmas, thread controls, BLAS dependencies, or optimized
  CPU performance claims in the current rewrite.

**Exit criteria:** explicit `openmp` requests fail clearly, while the reference
CPU remains deterministic and no interface change is needed to add OpenMP
later.

**Bounded Phase 4 execution record (2026-08-26):** no production changes were
necessary. Focused contract tests confirm that the reserved `openmp`
registration remains disabled and unimplemented, explicit requests repeatedly
return the same structured `kBackendUnavailable` error without a backend,
`auto` selects the serial FP64 CPU backend, and fixed traces reject `openmp`
without trace points or backend resource activity. The existing coarse backend
contract, owned sparse layouts, opaque move-only handles, preallocated
workspace, public legacy route, and rewrite `kNotReady` result remain
unchanged. In the fresh CPU-only build
`/home/yanyu/builds/hipdlp2-m4-cpu-20260826`, `[phase4]` passes 360 assertions
in 8 cases, `[phase3c]` 286 in 8, `[phase3b]` 643 in 8, `[phase3]` 1164 in 25,
`[migration]` 122 in 1, `[hipdlp]` 1592 in 32, and `[pdlp]` 1622 in 40. CMake
used `CUPDLP_GPU=OFF`, `CUPDLP_FIND_CUDA=OFF`,
`HIGHS_HIPDLP_CUDA=OFF`, and the unusable
`CUDACXX=/nonexistent/hipdlp-phase4-nvcc`. All 298 compile-database entries,
304 Ninja commands, and 82,919 resolved dependency lines contain no active
OpenMP or GPU toolchain flags, headers, or commands; `unit_tests` and
`libhighs.so.1.15.1` contain no OpenMP, CUDA, HIP, or ROCm runtime links or
symbols. CMake and Ninja were available; Meson, Bazel, Bazelisk, and buildifier
were unavailable, so their existing source and test registrations were
validated statically. OpenMP implementation, threading, dispatch, tuning,
BLAS integration, and performance work remain deferred.

### Phase 5: Implement CUDA

- Add RAII runtime/library wrappers and error translation.
- Upload problem data and allocate all workspace once.
- Implement SpMV, transpose SpMV, fused updates, reductions, convergence
  metrics, and solution transfer.
- First establish iteration parity without graph capture or aggressive fusion.
- Add graph capture, fusion, and sparse-algorithm tuning incrementally, with an
  unfused diagnostic mode retained for testing.
- Profile representative small, medium, and large instances.

**Exit criteria:** CUDA passes reference-CPU differential tests, full available
correctness suites, and the 5% geomean no-regression gate against frozen legacy
HiPDLP CUDA.

**Bounded Phase 5 execution record (2026-08-26):** the rewrite now has a
deterministic FP64 NVIDIA CUDA backend for the existing coarse primitive
contract and Phase 3C fixed trace. Private move-only RAII objects own the CUDA
stream, runtime allocations, cuBLAS/cuSPARSE handles, sparse/dense descriptors,
and reusable SpMV storage. Binding uploads the owned CSR and transpose layouts,
costs, column/row bounds, and row classifications once; all trace vectors and
workspace are allocated before iteration. Device-side cuSPARSE scalars and
resource auditing confirm zero allocation, descriptor creation, or host/device
transfer in the fixed-trace iteration sequence. Explicit batched vectors and
SpMM remain structured unsupported-capability errors. Explicit CUDA requests
return no object when the build or device is unavailable, `auto` selects CUDA
only after stream, vendor-handle, allocation, and kernel-execution probes pass,
and otherwise selects serial CPU. OpenMP and ROCm remain unavailable; public
`solver = "hipdlp"` remains on legacy and rewrite `algorithm::solve` remains
structured `kNotReady` without a valid solution.

The fresh CUDA build is
`/home/yanyu/builds/hipdlp2-m5-cuda-20260826`, configured with
`HIGHS_HIPDLP_CUDA=ON`, both `CUPDLP_GPU` and `CUPDLP_FIND_CUDA` off,
testing and compile commands enabled. On NVIDIA GB10
`GPU-bac48992-fc11-93f0-a8c1-f604826c83c5` (compute capability 12.1), driver
580.159.03 and CUDA 13.0 / nvcc 13.0.88, `[phase5]` passes 1,058 assertions in
11 cases, `[phase4]` 359/8, `[phase3c]` 892/11, `[phase3b]` 643/8, `[phase3]`
1,160/25, `[migration]` 122/1, `[hipdlp]` 2,626/42, and `[pdlp]` 2,656/50.
Focused primitive, trace, and resource filters pass 288/3, 659/4, and 117/4;
the normal and no-visible-device factory runs pass 34/1 and 28/1, with the
latter repeatedly returning structured `kBackendUnavailable`. CUDA memcheck
passes the repeated ownership test with 29 assertions, zero errors, and zero
leaked bytes.

The fresh CPU-only build is
`/home/yanyu/builds/hipdlp2-m5-cpu-20260826`, configured with
`HIGHS_HIPDLP_CUDA=OFF`, both cuPDLP CUDA options off, and
`CUDACXX=/nonexistent/hipdlp-phase5-nvcc`. `[phase5]` passes 26/1, `[phase4]`
365/8, `[phase3c]` 286/8, `[phase3b]` 643/8, `[phase3]` 1,169/25,
`[migration]` 122/1, `[hipdlp]` 1,597/32, and `[pdlp]` 1,627/40. Its 298
compile-database entries, 304 Ninja commands, resolved headers, dynamic links,
and symbols contain no CUDA, OpenMP, HIP, or ROCm dependency. The CUDA build
has exactly two CUDA translation units under `backend/cuda` among 300 compile
entries and 306 Ninja commands; no legacy HiPDLP/cuPDLP CUDA source is built.
`unit_tests` and `libhighs.so.1.15.1` link only the expected CUDA runtime/vendor
set (`libcudart.so.13`, `libcublas.so.13`, `libcublasLt.so.13`,
`libcusparse.so.12`, and `libnvJitLink.so.13`) and no OpenMP, HIP, or ROCm
library, flag, header, or symbol. No architecture or tuning flag is set in
source; CMake/toolkit detection populated architecture 75 in the build cache.
CMake 3.28.3, Ninja 1.11.1, compute-sanitizer 2025.3.1.0, and Nsight Systems
2025.3.2 are available. Meson, Bazel, Bazelisk, and buildifier are unavailable,
so their CPU-only source, configuration, and test registrations were validated
statically. The full solve loop, convergence/restart/adaptive policies, public
cutover, CUDA Graphs, mixed precision, tuning, benchmarking, and performance
claims/gates remain deferred.

**Bounded Phase 5G follow-up execution record (2026-08-26):** the rewrite CUDA
backend now has a reusable deterministic FP64 execution plan for only the
existing fixed Halpern-PDHG trace. A backend-neutral move-only
`ExecutionPlanHandle` and `kExecutionPlan` capability keep CUDA types private;
the serial CPU backend rejects that capability explicitly, and the existing
direct CUDA primitive path remains the default and the correctness oracle.
The CUDA implementation privately owns move-only `cudaGraph_t` and
`cudaGraphExec_t` resources and uses the existing one non-blocking stream,
uploaded CSR/transpose matrices, cuSPARSE descriptors, dense descriptors,
device vectors, and reusable SpMV external buffer.

Each plan instantiates six fixed executables before replay: reset-bearing first
one- and two-iteration graphs, reset-free continuation one- and two-iteration
graphs, and equivalent first-block graphs used to inject an asynchronous device
failure transactionally in tests. The first-block graphs restore the initial
ping-pong role and point-zero products, reset the sticky error state, and record
the initial snapshot. The two-iteration graphs contain both ping-pong
directions and therefore finish with the same logical current-buffer role.
Preallocated device scalars hold the iteration index, the device-computed
Halpern weight `1 / (iteration + 1)`, fixed primal and dual step sizes, and the
sticky error state. Captured work contains the existing cuSPARSE SpMV path,
cuBLAS copies, custom projected/reflection/Halpern kernels, and device-indexed
trace snapshots. No host node, callback, persistent replacement kernel,
device-side allocation, graph update, or graph node parameter update is used.

Successful replay calls only `cudaGraphLaunch`: zero, one, two, and four
iterations issue exactly 0, 1, 1, and 2 launches, respectively; a five-step
two-block-plus-remainder run issues 3. Replay auditing reports zero graph
capture, graph instantiation, graph update, allocation, descriptor creation,
host-to-device transfer, device-to-host transfer, direct primitive enqueue,
and synchronization calls. Four final bulk trace downloads, one sticky-error
download for nonzero traces, and one final synchronization occur only after
replay. Construction, instantiation, launch, and asynchronous execution
failures return structured errors with no plan or result; captured-resource
lifetime tokens reject a stale plan before launch. Repeated execution of the
same instantiated plan is bitwise deterministic, including recovery after a
sticky device failure and a zero-step run with no stale state.

Real capture support was established before production edits with a standalone
non-blocking-stream probe using device-pointer cuSPARSE FP64 SpMV plus a custom
kernel. On NVIDIA GB10
`GPU-bac48992-fc11-93f0-a8c1-f604826c83c5` (compute capability 12.1), it
captured 4 kernel nodes, instantiated and replayed twice, and returned
`output=[8,7]`, `error=0`. The driver is 580.159.03; the CUDA runtime and driver
API report 13.0; nvcc is 13.0.88. CMake is 3.28.3, Ninja is 1.11.1, and
compute-sanitizer is 2025.3.1.0.

The fresh CUDA build is
`/home/yanyu/builds/hipdlp2-m5g-cuda-20260826`, configured with
`HIGHS_HIPDLP_CUDA=ON`, both cuPDLP CUDA options off, testing and compile
commands enabled. `[phase5g]` passes 3,027 assertions in 12 cases; focused
factory, graph parity, one-step remainder, same-plan reuse, resource, and
failure filters pass 4/1, 2,540/3, 102/1, 61/1, 276/4, and 36/1. `[phase5]`
passes 1,067/11, `[phase4]` 366/8, `[phase3c]` 892/11, `[phase3b]` 643/8,
`[phase3]` 1,167/25, `[migration]` 122/1, `[hipdlp]` 5,667/54, and `[pdlp]`
5,697/62. Compute-sanitizer memcheck over four complete plan
construction/capture/instantiation/replay/download/destruction cycles passes
129 assertions with zero errors, invalid accesses, or leaked bytes.

The fresh CPU-only build is
`/home/yanyu/builds/hipdlp2-m5g-cpu-20260826`, configured with all CUDA options
off and `CUDACXX=/nonexistent/hipdlp-phase5g-nvcc`. `[phase5g]` passes 8/1,
`[phase5]` 28/1, `[phase4]` 372/8, `[phase3c]` 286/8, `[phase3b]` 643/8,
`[phase3]` 1,176/25, `[migration]` 122/1, `[hipdlp]` 1,612/33, and `[pdlp]`
1,642/41. Its 298 compile-database entries, 304 Ninja commands, and 82,918
resolved dependency lines contain zero CUDA, cuBLAS, cuSPARSE, OpenMP, HIP, or
ROCm compiler, command, or header dependencies; its test and shared-library
artifacts contain no corresponding links or runtime symbols.

The CUDA build has 300 compile-database entries (261 `clang++`, 37 `clang`, and
exactly 2 `nvcc`) and 306 Ninja commands. The only CUDA translation units are
`backend/cuda/cuda_backend.cu` and `backend/cuda/cuda_kernels.cu`; no legacy
HiPDLP or cuPDLP CUDA source is built. Its 83,480 resolved dependency lines
contain 117 active CUDA headers and no OpenMP, HIP, or ROCm headers. The test
and shared library link only the expected CUDA set (`libcudart.so.13`,
`libcublas.so.13`, `libcublasLt.so.13`, `libcusparse.so.12`, and
`libnvJitLink.so.13`) and no OpenMP, HIP, or ROCm library. Dynamic symbols
confirm stream begin/end capture, graph instantiation, launch, and RAII
destruction. CMake records architecture 75 from toolkit detection; this slice
adds no architecture-specific flag, tuning, benchmark, or performance claim.
The legacy implementation-local graph calls in `pdhg.cc` predate the rewrite
base and remain unchanged; Phase 5G introduces no graph type or API outside the
rewrite CUDA backend.

Public `solver = "hipdlp"` still routes to legacy, rewrite
`algorithm::solve` still returns structured `kNotReady` without a solution,
OpenMP and ROCm remain unavailable, and the independent cuPDLP CUDA options
remain disabled. The full solve loop, convergence, restart, adaptive step-size
policies, public routing/cutover, multi-stream scheduling, fusion and sparse
algorithm tuning, mixed precision, multistart, benchmarking, and performance
claims remain deferred.

**Bounded Phase 5H architecture-correction execution record (2026-08-26):**
Phase 5G's numerical evidence remains the regression baseline, but its
fixed-iteration ownership model is superseded. The backend-specific
`FixedIterationPlanVectors`, `ExecutionPlanHandle`, and fixed-trace plan
methods duplicated the ordered Halpern-PDHG step in the CUDA backend. They are
removed from `backend/types.hpp`, the universal `Backend` contract, and the
serial CPU backend. CPU has no recording overrides or graph-specific stubs,
and backend contracts no longer expose algorithm-state names such as anchors,
reflected iterates, current products, or Halpern iteration. The paired
primitive needed by the algorithm is the backend-neutral `fusedVectorBlend`.

`algorithm/halpern_pdhg_step.{hpp,cc}` now owns `HalpernPdhgState` and the one
source-level definition of `enqueueHalpernPdhgStep`. That definition alone
orders projected primal update, primal reflection, reflected `Ax`, projected
dual update, dual reflection, anchor weight `1 / (iteration + 1)`, paired
blend, refreshed `Ax` and `A^T y`, ping-pong role exchange, successful
iteration advance, and first-error return. Direct serial CPU, direct CUDA,
and recorded CUDA fixed traces all call that definition. Fixed-trace snapshots
and their device handles remain in the algorithm diagnostic/test layer; CUDA
does not own trace shape, iteration state, step-size state, anchor-weight
state, ping-pong roles, or a standalone ordered PDHG sequence. Fixed-trace
backend selection now uses only FP64 precision and advertised capabilities,
with no backend-kind check.

The optional acceleration boundary is the generic
`backend::CommandBufferRecorder` extension and opaque
`CommandBufferHandle`. Its narrow lifecycle is begin recording, end or cancel
recording, and replay. It contains no LP, PDHG, iteration, or trace contract.
The CUDA implementation translates CUDA errors and privately owns graph and
executable RAII resources, but it only captures primitive calls issued by the
algorithm. The algorithm records at most two calls to the shared step in each
command buffer, so 0, 1, 2, 4, and 5 iterations replay with exactly 0, 1, 1,
2, and 3 CUDA Graph launches. Replay remains launch-only: tests observe zero
direct primitive enqueue, allocation, descriptor creation, graph update,
host-to-device transfer, device-to-host transfer, or synchronization calls.
The existing single non-blocking stream, device-pointer cuSPARSE FP64 SpMV,
dense/sparse descriptors, and reusable external SpMV buffer are unchanged.

The narrow `[phase5h][algorithm][shared-step]` test falsifies the ownership
hypothesis with a five-iteration odd-remainder trace: CPU direct, CUDA direct,
and CUDA recorded points agree within the established `5e-13` CUDA margin and
the audit reports exactly three graph launches. `[phase5h][structure]` reads
the source boundary and asserts that the removed fixed-plan types/methods are
absent, CPU has no recorder overrides, CUDA has no ordered iteration helper,
both direct and recording paths reference the shared definition, and no CUDA
type or graph vocabulary enters algorithm, API, model, legacy, or public
solver options. Existing Phase 5G tests preserve parity for 0, 1, 2, and 4
iterations, all scaling modes, equality/inequality and column projections,
empty/rectangular/skewed fixtures, deterministic replay, remainder handling,
transactionality, stale-resource rejection, and balanced lifetimes.

Fresh CMake builds are
`/home/yanyu/builds/hipdlp2-m5h-cuda-20260826` and
`/home/yanyu/builds/hipdlp2-m5h-cpu-20260826`. Both use Ninja, Release,
Clang, shared libraries, `FAST_BUILD=OFF`, `BUILD_TESTING=ON`, examples off,
`HIGHSINT64=OFF`, and both cuPDLP CUDA options off. The CUDA build alone sets
`HIGHS_HIPDLP_CUDA=ON` with `/usr/local/cuda/bin/nvcc`; the CPU build sets it
off and was configured with `CUDACXX=/definitely/unavailable/nvcc`.

The final CUDA matrix passes `[phase5h]` 252 assertions/2 cases, `[phase5g]`
2,842/12, `[phase5]` 1,068/11, `[phase4]` 367/8, `[phase3c]` 892/11,
`[phase3b]` 643/8, `[phase3]` 1,168/25, `[migration]` 122/1, `[hipdlp]`
5,735/56, and `[pdlp]` 5,765/64. The final CPU matrix passes `[phase5h]`
149/1, `[phase5g]` 8/1, `[phase5]` 29/1, `[phase4]` 373/8, `[phase3c]`
286/8, `[phase3b]` 643/8, `[phase3]` 1,177/25, `[migration]` 122/1,
`[hipdlp]` 1,762/34, and `[pdlp]` 1,792/42. The expected rewrite-not-ready
log appears in routing tests; public `solver = "hipdlp"` still uses legacy and
`algorithm::solve` still returns structured `kNotReady` with no solution.

Compute-sanitizer memcheck over four command-buffer
construction/recording/finalization/replay/destruction cycles passes 57
assertions in one case with zero errors and zero leaked bytes. Failure
injection and `ResourceAudit` state are compiled only under the private
`HIGHS_HIPDLP_TESTING` definition. Normal recording creates no failure-only
command buffer. The separate production build
`/home/yanyu/builds/hipdlp2-m5h-cuda-production-20260826`, configured with
`BUILD_TESTING=OFF`, builds `libhighs.so.1.15.1`; its compile database contains
no test definition, and static symbols, dynamic symbols, and strings contain
zero failure-point, `ResourceAudit`, test-state, or injected-failure matches.

The CPU build has 300 compile-database entries, 306 Ninja commands, and 83,534
resolved dependency lines with zero CUDA, cuBLAS, cuSPARSE, OpenMP, HIP, or
ROCm compiler, command, exact vendor-header, or linked-library dependencies.
The CUDA build has 302 compile entries, 308 Ninja commands, and 84,098
dependency lines. Its only CUDA translation units and vendor-header owners are
`backend/cuda/cuda_backend.cu` and `backend/cuda/cuda_kernels.cu`; no legacy
HiPDLP/cuPDLP CUDA, OpenMP, HIP, or ROCm source is built. `unit_tests` and
`libhighs.so.1.15.1` link only the expected CUDA runtime/vendor set
(`libcudart.so.13`, `libcublas.so.13`, `libcublasLt.so.13`,
`libcusparse.so.12`, and `libnvJitLink.so.13`). Meson, Bazel, Bazelisk, and
buildifier remain unavailable; their CPU-only registrations were checked
statically.

Deferred Tier A scope remains the full solve loop, convergence and stopping,
restart, adaptive/PID/Malitsky-Pock policies, public rewrite routing, Tier-A
iteration blocks, and conditional recorded nodes. Deferred Tier B scope
remains multi-stream scheduling, fusion and sparse-algorithm tuning, mixed
precision, multistart, ROCm implementation, benchmarking, performance gates,
and performance claims. Phase 5H implements none of those items.

### Phase 6: Reserve ROCm extension

- Keep `rocm` as an explicit unavailable backend capability.
- Preserve backend semantics, precision descriptors, and error types required
  by a future HIP/rocBLAS/rocSPARSE implementation.
- Do not add ROCm source files, toolchain discovery, compile-validation claims,
  or hardware-performance claims in the current rewrite.

**Exit criteria:** explicit `rocm` requests fail clearly and CPU/CUDA builds do
not require ROCm headers, tools, libraries, or environment variables.

**Bounded Phase 6 execution record (2026-08-26):** no production change was
necessary. The existing backend-neutral `Kind`, `Descriptor`, and factory
contract already reserves `Kind::kRocm` as the valid registration named
`"rocm"`, with `build_enabled=false`, `implemented=false`, default FP64
storage/compute/reduction precision, and an empty capability set. Repeated
explicit requests return the same structured `kBackendUnavailable` error and
message with no backend object. Automatic selection remains CUDA-then-CPU and
contains no ROCm probe or selection path.

Focused tests were added only to `TestHiPdlpBackend.cpp` and
`TestHiPdlpTrace.cpp`. The first post-edit CPU check,
`[phase6][rocm][factory]`, passed 25 assertions in one case. The final
`[phase6]` filter passes 1,270 assertions in four cases in both CPU and CUDA
builds. It covers the neutral descriptor and deterministic factory error,
CUDA-then-CPU automatic selection, direct and recorded fixed-trace rejection
with no partial trace, resources, transfers, command replay, or
synchronization, and source/build isolation. The structural test confirms
that no `RocmBackend` or `backend/rocm` directory exists; no HIP, hipBLAS,
hipSPARSE, rocBLAS, rocSPARSE, or ROCm type enters algorithm, API, model,
legacy, public options, the CPU backend, or the universal `Backend` contract;
and CMake, Bazel, and Meson contain no ROCm option, language enablement,
compiler detection, source registration, compile definition, include path, or
linked library. The existing Phase 5H structure case is also tagged for Phase
6 and continues to enforce the single algorithm-owned Halpern-PDHG step,
generic optional command-buffer boundary, and absence of CPU recording stubs.

Fresh test builds are `/home/yanyu/builds/hipdlp2-m6-cuda-20260826` and
`/home/yanyu/builds/hipdlp2-m6-cpu-20260826`. Both use Ninja, Release, Clang,
shared libraries, `FAST_BUILD=OFF`, `BUILD_TESTING=ON`, examples off,
`HIGHSINT64=OFF`, and both cuPDLP CUDA options off. The CUDA build alone uses
`HIGHS_HIPDLP_CUDA=ON` and `/usr/local/cuda/bin/nvcc`; the CPU build uses
`HIGHS_HIPDLP_CUDA=OFF` with
`CUDACXX=/definitely/unavailable/nvcc`. Before the first edit, both fresh
`unit_tests` targets built and all requested Phase 3 through Phase 5H,
migration, HiPDLP, and PDLP filters passed with the Phase 5H counts.

The final CUDA matrix passes `[phase6]` 1,270 assertions/4 cases,
`[phase5h]` 252/2, `[phase5g]` 2,842/12, `[phase5]` 1,068/11, `[phase4]`
367/8, `[phase3c]` 892/11, `[phase3b]` 643/8, `[phase3]` 1,168/25,
`[migration]` 122/1, `[hipdlp]` 6,856/59, and `[pdlp]` 6,886/67. The final
CPU matrix passes `[phase6]` 1,270/4, `[phase5h]` 149/1, `[phase5g]` 8/1,
`[phase5]` 29/1, `[phase4]` 373/8, `[phase3c]` 286/8, `[phase3b]` 643/8,
`[phase3]` 1,177/25, `[migration]` 122/1, `[hipdlp]` 2,883/37, and `[pdlp]`
2,913/45. The expected rewrite-not-ready log remains in routing tests: public
`solver = "hipdlp"` still uses legacy and rewrite `algorithm::solve` still
returns structured `kNotReady` without a solution.

The CPU build has 300 compile-database entries, 306 Ninja commands, and
83,551 resolved dependency lines with no CUDA, cuBLAS, cuSPARSE, OpenMP, HIP,
hipBLAS, hipSPARSE, rocBLAS, rocSPARSE, or ROCm compiler, command, header,
link, or vendor-runtime symbol dependency. The CUDA build has 302 compile
entries, 308 Ninja commands, and 84,115 dependency lines. Its only CUDA
translation units and CUDA-header owners are
`backend/cuda/cuda_backend.cu` and `backend/cuda/cuda_kernels.cu`; no legacy
HiPDLP/cuPDLP CUDA, OpenMP, HIP, or ROCm source or toolchain is present.
`unit_tests` and `libhighs.so.1.15.1` link only `libcudart.so.13`,
`libcublas.so.13`, `libcublasLt.so.13`, `libcusparse.so.12`, and
`libnvJitLink.so.13`, with no OpenMP, HIP, or ROCm library or API symbol.

The fresh production build is
`/home/yanyu/builds/hipdlp2-m6-cuda-production-20260826`, configured with the
same CUDA isolation settings and `BUILD_TESTING=OFF`. It builds
`libhighs.so.1.15.1` from 240 compile entries and 244 Ninja commands, including
only the same two rewrite CUDA translation units. Its compile definitions,
full and dynamic symbols, and strings contain zero
`HIGHS_HIPDLP_TESTING`, `FailurePoint`, `ResourceAudit`, failure-injection, or
test-state matches. Compute Sanitizer 2025.3.1.0 memcheck over the existing
four-cycle `hi-pdlp-cuda-graph-resource-lifetimes-are-balanced` regression
passes 57 assertions in one case with zero errors and zero leaked bytes. This
is only a CUDA command-buffer lifecycle regression; no ROCm sanitizer,
runtime, device, compilation, or hardware claim is made.

CMake 3.28.3, Ninja 1.11.1, clang-format 22.1.8, nvcc 13.0.88, Compute
Sanitizer 2025.3.1.0, and Nsight Systems 2025.3.2 are available. Meson, Bazel,
Bazelisk, buildifier, `hipcc`, `rocminfo`, and `rocm-smi` are unavailable, so
Bazel and Meson CPU registrations were validated statically and no ROCm
toolchain or hardware check was attempted.

Deferred Tier A scope remains the full solve loop, convergence and stopping,
restart, adaptive/PID/Malitsky-Pock policies, public rewrite routing, Tier-A
iteration blocks, and conditional recorded nodes. Deferred Tier B scope
remains multi-stream scheduling, fusion and sparse-algorithm tuning, mixed
precision, multistart, ROCm implementation and toolchain support,
benchmarking, performance gates, and performance claims. Phase 6 implements
none of those items.

### Tier A1 prerequisite to Phase 7: convergence and stopping policy

This bounded prerequisite tests the hypothesis that a pure algorithm-owned
policy can normalize backend-neutral scalar metrics and choose continue,
optimal, iteration-limit, time-limit, or numerical error without changing the
`Backend` contract or solver routing. `algorithm/convergence.{hpp,cc}` now owns
that policy. It accepts only scalar metrics, tolerances, progress, and limits;
it contains no backend handle, device type, explicit allocation,
synchronization, or KKT-verifier call. The independent verifier remains the
eventual final-solution acceptance gate.

The policy computes
`primal_residual / (1 + unscaled_rhs_norm)`,
`dual_residual / (1 + unscaled_cost_norm)`, and
`objective_gap / (1 + abs(primal_objective) + abs(dual_objective))` in a fixed
order. All three results must be at or below their separate finite positive
tolerances for an optimal decision. Validation rejects invalid tolerances,
negative residuals, gaps, norms, iterations, elapsed time, or limits,
non-finite scalar inputs other than the exact unlimited time sentinel, and
non-finite normalization intermediates or results. Errors return
`kNumericalError`, an `api::Error`, and invalid zero-initialized normalized
metrics rather than a success-shaped decision. `kHighsIInf` and `kHighsInf`
disable their corresponding limits; finite limits are inclusive. The fixed
precedence is numerical error, optimal, iteration limit, time limit, then
continue.

TDD used `[phase7a][convergence][policy]` as the narrow falsifier. The first
substantive edit added only `TestHiPdlpConvergence.cpp` and test registrations.
The immediate CPU `unit_tests` build failed as expected at the missing
`pdlp/hipdlp/algorithm/convergence.hpp`, before link, so no test binary was
run. After the minimum implementation and source registration, the narrow
filter passed 427 assertions in eight cases in both CPU and CUDA builds. It
covers exact normalization, equality at all three thresholds, each metric
failing independently, zero and subnormal objective scales, every invalid
scalar category, normalization overflow, inclusive and unlimited limits,
decision precedence, bitwise-repeatable synthetic sequences, and the absence
of backend/resource/accelerator/verifier vocabulary in the policy. Only
`TestHiPdlpConvergence.cpp` and the two new policy files were formatted, using
clang-format 22.1.8.

The fresh test builds are
`/home/yanyu/builds/hipdlp2-m7a-cpu-20260826` and
`/home/yanyu/builds/hipdlp2-m7a-cuda-20260826`. Both use Ninja, Release,
Clang 18.1.3, shared libraries, `FAST_BUILD=OFF`, `BUILD_TESTING=ON`, examples
off, `HIGHSINT64=OFF`, exported compile commands, and both cuPDLP CUDA options
off. The CPU build sets `HIGHS_HIPDLP_CUDA=OFF` and was configured with
`CUDACXX=/definitely/unavailable/nvcc`; the CUDA build sets
`HIGHS_HIPDLP_CUDA=ON`, uses `/usr/local/cuda/bin/nvcc` 13.0.88, and retains
architecture 75. Before the first edit, both fresh `unit_tests` targets built
and every inherited filter passed with the Phase 6 counts.

The final CPU matrix passes `[phase7a]` 427 assertions/8 cases, `[phase6]`
1,270/4, `[phase5h]` 149/1, `[phase5g]` 8/1, `[phase5]` 29/1, `[phase4]`
373/8, `[phase3c]` 286/8, `[phase3b]` 643/8, `[phase3]` 1,177/25,
`[migration]` 122/1, `[hipdlp]` 3,310/45, and `[pdlp]` 3,340/53. The final
CUDA matrix passes `[phase7a]` 427/8, `[phase6]` 1,270/4, `[phase5h]` 252/2,
`[phase5g]` 2,842/12, `[phase5]` 1,068/11, `[phase4]` 367/8, `[phase3c]`
892/11, `[phase3b]` 643/8, `[phase3]` 1,168/25, `[migration]` 122/1,
`[hipdlp]` 7,283/67, and `[pdlp]` 7,313/75. The expected
rewrite-not-ready log remains in routing tests.

The CPU audit records 302 compile-database entries, 308 Ninja commands, and
84,166 resolved dependency lines. Exact compiler, command, header, dynamic
link, and dynamic-symbol searches find zero CUDA, cuBLAS, cuSPARSE, OpenMP,
HIP, hipBLAS, hipSPARSE, rocBLAS, rocSPARSE, or ROCm dependency. The CUDA
audit records 304 compile entries (265 `clang++`, 37 `clang`, and exactly two
`nvcc`), 310 Ninja commands, and 84,730 dependency lines. The only CUDA
translation units and owners of all 117 resolved CUDA-header lines remain
`backend/cuda/cuda_backend.cu` and `backend/cuda/cuda_kernels.cu`; no legacy
HiPDLP/cuPDLP CUDA source is built. `unit_tests` and `libhighs.so.1.15.1`
link only `libcudart.so.13`, `libcublas.so.13`, `libcublasLt.so.13`,
`libcusparse.so.12`, and `libnvJitLink.so.13`, with zero OpenMP, HIP, or ROCm
link or symbol match.

The fresh production build is
`/home/yanyu/builds/hipdlp2-m7a-cuda-production-20260826`, configured with the
same CUDA isolation settings and `BUILD_TESTING=OFF`. It builds
`libhighs.so.1.15.1` and `highs` from 241 compile entries, 245 Ninja commands,
and 65,281 resolved dependency lines. It contains no test source or
`unit_tests` target, no `HIGHS_HIPDLP_TESTING` compile or Ninja definition,
and zero full-symbol, dynamic-symbol, or string matches for `FailurePoint`,
`ResourceAudit`, failure injection, or test state. Exactly the same two CUDA
translation units own all 117 CUDA-header lines and the library links the same
five expected CUDA libraries with no OpenMP, HIP, or ROCm dependency.

CMake, Bazel, and Meson registrations include the new source and test. Bazel
also declares the two policy source files as runfiles for the structural test.
CMake 3.28.3, Ninja 1.11.1, Clang 18.1.3, CUDA 13.0.88, Compute Sanitizer, and
Nsight Systems are available. Meson, Bazel, Bazelisk, buildifier, `hipcc`,
`rocminfo`, and `rocm-smi` remain unavailable, so Bazel and Meson registrations
were validated statically and no ROCm build, runtime, or hardware claim is
made. Pre-edit and final SHA-256 checks confirm that the user changes in
`TestHiPdlpBackend.cpp` and `TestHiPdlpTrace.cpp`, the universal backend
contract, every CUDA source/header, `algorithm::solve`, and legacy-default
routing are byte-for-byte unchanged.

### Tier A2 prerequisite to Phase 7: weighted averaging integration

The existing scalar `algorithm/averaging_policy.{hpp,cc}` was retained
unchanged after `[phase7b]` passed in Clang CPU, Clang CUDA, and GCC CPU-only
builds. No additional policy tests were added. The fixed trace now allocates
primal and dual weighted sums before iteration, calls `updateAveraging` after
each completed Halpern step, and passes the accepted sample weight and current
iterates to the existing `Backend::fusedAverageUpdate` primitive. Optional
per-iteration step sequences exercise variable-step weighting without changing
the Halpern or backend contracts.

After the fixed iteration count, the trace applies the final policy
normalization and downloads the normalized primal and dual averages. The
result exposes those values and `AveragingState`; zero iterations have no valid
average, and scalar or backend failures expose neither a partial trace nor a
partial average. Focused CPU/CUDA tests cover exact one-step and variable-step
averages, failure transactionality, and direct CPU/direct CUDA/recorded CUDA
parity. Recorded execution retains 0/1/1/2/3 graph launches for 0/1/2/4/5
iterations, with no iteration-sequence allocation, transfer, or
synchronization.

All requested Clang CPU and CUDA regression filters pass. GCC passes
`[phase7a]` and `[phase7b]`; its broad `[hipdlp]` filter still differs only in
the pre-existing compiler-sensitive legacy `scrs8` frozen trajectory. Restart,
convergence cadence, final iterate selection, public rewrite routing, and the
structured `algorithm::solve` `kNotReady` result remain deferred.

Remaining Tier A work includes the full solve loop, restart, convergence
cadence and final iterate selection, step-size initialization,
fixed/PID/Malitsky-Pock updates, power iteration, infeasibility certificates,
Tier-A iteration blocks, conditional command recording, and public rewrite
routing. Cross-solver validation has not begun. Tier B benchmarking, tuning,
mixed precision, and ROCm implementation also remain deferred.

### Tier A1 execution integration: final convergence observation

The existing `algorithm/convergence.{hpp,cc}` policy was reverified and left
unchanged. `[phase7a]` passed before any edit in the Clang CPU, Clang CUDA, and
GCC CPU-only builds, so no compiler or focused-test failure proved a defect.

The existing surfaces are sufficient, so the feasibility stop was not reached.
Raw scalar metrics are computed in the transformed problem's unscaled units,
the space that `Problem::unscaledNorms()` measures, from `costs`, `rowLower`,
`colLower`, `colUpper`, `equalityRows`, `scaling().row_scale`/`col_scale`,
`objectiveOffset`, and the host trace's completed primal, dual, `Ax`, and
`A^T y`. Residuals reproduce the legacy projected row violation and
`c - A^T y - slack_+ + slack_-` definitions. The averaged iterate reuses the
existing `Backend::multiply` and `Backend::multiplyTranspose` after the fixed
iterations. No `Backend` method, CUDA kernel, or in-loop host download was
added. `model::postsolve` is not used for raw metrics because it returns the
original-space solution, whose rows and bounds no longer match the transformed
`unscaled_rhs_norm` and `unscaled_cost_norm` denominators.

The fixed trace evaluates convergence exactly once, after all fixed iterations.
It observes the completed current iterate and, when `AveragingState::
hasAverage()` holds, the completed normalized averaged iterate, each through
one `StoppingInput` and one `evaluateStopping` call, and exposes both as
`FixedTraceResult::current_convergence` and `averaged_convergence`. Nothing
selects between them. No cadence, periodic host decision, or early stop was
added. `FixedTraceOptions` gains only `convergence_tolerances`, defaulting to
`1e-4`; stopping limits stay unlimited. Zero iterations produce no averaged
observation, and scalar, backend, or metric failures return the existing
transactional failure result with neither observation. The independent KKT
verifier is not called and remains the separate acceptance gate.

Two vectors for the averaged products are created before the iteration
sequence; nothing was added inside it. After-iteration downloads rise from
`4(n+1)+2` to `4(n+1)+4` when an average exists. Recorded CUDA execution still
records generically and retains 0/1/1/2/3 graph launches for 0/1/2/4/5
iterations, and CPU, direct CUDA, and recorded CUDA observations agree.

Focused `[phase7a]` tests in `TestHiPdlpTrace.cpp` cover exact final
current-iterate metrics on the one-step fixture, exact averaged-iterate metrics
on a two-iteration free-column fixture whose average `1.5` differs from the
current iterate `2.0`, zero-iteration and failure transactionality, and
CPU/direct-CUDA/recorded-CUDA parity with unchanged graph launch counts. The
existing bitwise determinism and CUDA parity comparisons were extended to the
two observations.

Clang CPU passes `[phase7a]` 569 assertions/11 cases, `[phase7b]` 344/13,
`[phase5h]` 149/1, `[phase5g]` 8/1, `[phase5]` 29/1, `[phase3c]` 461/8,
`[hipdlp]` 3,971/61, and `[pdlp]` 4,001/69. Clang CUDA passes 1,218/12,
993/14, 320/2, 5,807/12, 1,705/11, 1,704/11, 11,614/83, and 11,644/91. GCC
13.3 CPU-only rebuilds and passes `[phase7a]` 569/11 and `[phase7b]` 344/13;
the known GCC-only legacy `scrs8` frozen-trajectory difference remains a
pre-existing baseline issue outside this slice.

`clang-format` 22.1.8 ran only on the three touched files,
`algorithm/fixed_trace.hpp`, `algorithm/fixed_trace.cc`, and
`check/TestHiPdlpTrace.cpp`. `git diff --check`, the cached check, and editor
diagnostics are clean, and nothing is staged or committed. Restart remains
deferred, public `solver="hipdlp"` still routes to legacy, `algorithm::solve`
still returns the structured `kNotReady` result, and Phase 7 cross-solver
validation has not begun.

### Tier A3 prerequisite to Phase 7: scalar restart policy

`algorithm/restart_policy.{hpp,cc}` now owns only normalized scalar restart
history, candidate reporting, and the deterministic restart-or-continue rule.
The focused test was first confirmed to fail at the missing header; final
`[phase7c]` passes 178 assertions in four cases in Clang CPU, Clang CUDA, and
GCC CPU-only builds. The policy remains unintegrated, public `solver="hipdlp"`
still uses legacy, `algorithm::solve` still returns structured `kNotReady`, and
Phase 7 cross-solver validation has not begun.

### Phase 7: Full cross-solver validation

- Run rewritten CPU and CUDA against legacy HiPDLP, in-tree cuPDLP, external
  cuPDLPx, pinned cuOpt, and HiGHS reference solvers.
- Verify every claimed optimal result independently.
- Triage status, objective, or iteration-count disagreements.
- Test time/iteration limits and failure propagation.
- Produce machine-readable pass/fail summaries and a human-readable report.

**Exit criteria:** no unexplained incorrect status, non-finite result, or KKT
failure remains in the supported suite.

### Phase 8: Cut over and remove legacy code

- Switch public `solver = "hipdlp"` to the rewrite.
- Retain the legacy selector for one stabilization window if maintenance cost
  is acceptable.
- Confirm CMake CPU/CUDA configurations, unavailable-placeholder behavior, and
  Bazel/Meson CPU builds.
- Update user and developer documentation.
- Remove legacy solver, duplicated CUDA state, obsolete macros, and temporary
  differential hooks after the stabilization window.

**Exit criteria:** all Definition of Done items below pass from a clean build,
and the legacy implementation is no longer needed to diagnose routine failures.

## 9. Build and CI Matrix

Minimum matrix:

| Build | Compiler/toolchain | Backends | Validation |
| --- | --- | --- | --- |
| CMake CPU minimal | GCC and Clang | CPU | Unit and integration |
| CMake CUDA | NVCC + host compiler | CPU + CUDA | Unit, differential, GPU integration |
| Bazel CPU | Supported compiler | CPU | Existing build/tests remain green |
| Meson CPU | Supported compiler | CPU | Existing build/tests remain green |

Factory tests must additionally verify that explicit `openmp` and `rocm`
requests return clear unavailable-backend errors.

GPU dependencies must be optional. A normal CPU build must not require CUDA or
ROCm headers, tools, libraries, or environment variables.

## 10. Diagnostics and Error Handling

- Replace process termination and development assertions with typed errors.
- Attach operation, backend, vendor status code, and device context to backend
  failures.
- Preserve the first failure; do not return a success-shaped partial solution.
- Check all dimensions, buffer ownership, and backend compatibility at setup.
- Make unavailable backend, unsupported index width, allocation failure,
  non-finite iterate, and vendor-library failure distinguishable.
- Use one iteration log schema across all backends.
- Include backend, precision, device, algorithm policy, scaling, and effective
  options in the solve summary.
- Keep detailed profiling disabled by default and compile/runtime selectable.

## 11. Research Extension Points

### 11.1 Low-precision GPU SpMV

After FP64 cutover:

- add storage/compute/reduction precision combinations to backend capabilities;
- begin with low-precision matrix storage and SpMV while retaining FP64 master
  iterates and convergence reductions;
- add residual replacement or iterative correction before accepting a mode;
- verify stopping criteria using FP64 independent KKT calculations;
- compare convergence iterations, total time, memory traffic, and final KKT,
  not kernel throughput alone;
- fall back only when `auto` selected the mode; explicit precision requests
  must fail clearly if unsupported.

### 11.2 Multistart warm-up

After the single-start solver is stable:

- extend backend workspace to a batch dimension;
- evaluate SpMM only when effective batch size is greater than one;
- keep candidate generation, deduplication, scoring, and selection outside the
  linear algebra backend;
- define a fixed work or time budget so baseline comparisons are fair;
- transfer only the selected start into the main solve state;
- require end-to-end improvement, not isolated SpMM speedup.

Prior multistart results under
`/home/yanyu/benchmark-results/lp/mlo-spark/cupdlpx/multistart` show that
requested batch size can collapse to effective size one and that SpMM is not
automatically faster. The future design must record requested/effective batch
size and retain a SpMV path.

## 12. Risks and Mitigations

| Risk | Mitigation |
| --- | --- |
| Backend abstraction costs performance | Use coarse fused calls, preallocated workspace, and measured dispatch overhead |
| CPU/GPU algorithm drift | Keep one solve loop and policy implementation; differential-test decisions |
| Vendor APIs force incompatible data layouts | Let each backend own an uploaded representation behind one semantic contract |
| Corrected convergence changes iteration counts | Gate on independent KKT and record intentional differences |
| GPU timing hides synchronization | Define synchronized timing boundaries and retain end-to-end wall time |
| Tiny LPs regress from backend abstraction overhead | Include tiny LPs in the gate and measure coarse dispatch overhead |
| Placeholder backend appears supported | Explicit requests return unavailable; do not add source/toolchain options or support claims |
| External solver versions drift | Pin revisions/releases and hash executables/configuration |
| Benchmark selection bias | Use versioned manifests, randomized order, paired ratios, and immutable raw artifacts |
| Future low precision compromises correctness | Keep FP64 verifier/reductions and make precision an explicit capability |
| Multistart adds work without benefit | Gate on effective batch size and end-to-end time |

## 13. Definition of Done

The rewrite is complete when:

- `solver = "hipdlp"` and relevant existing options retain documented behavior;
- serial CPU and CUDA are selectable at runtime when built;
- rewritten HiPDLP uses deterministic Curtis-Reid scaling by default, with
  tested scale composition and exact primal/dual/reduced-cost recovery;
- explicit unavailable backend requests fail clearly;
- OpenMP and ROCm are reported as unavailable placeholders without support
  claims;
- CPU-only builds have no GPU dependency;
- algorithm code contains no CUDA/ROCm types or conditional branches;
- no steady-state iteration allocates memory or performs unintended transfers;
- all backend primitives pass reference differential tests;
- supported feasible, infeasible, unbounded, limit, and numerical-error cases
  map to correct HiGHS statuses;
- every claimed optimal result passes the independent `1e-4` KKT gates;
- full Netlib and energy comparisons have no unexplained correctness
  disagreement against the baseline set;
- serial CPU and CUDA each satisfy their documented paired legacy
  performance gates;
- every greater-than-20% instance regression is resolved or artifact-backed;
- CMake backend builds pass and Bazel/Meson CPU builds remain working;
- benchmark results include reproducible commands, hashes, environment,
  metadata, raw logs, parsed metrics, and verifier output;
- legacy code is removed after the stabilization window;
- low-precision SpMV and multistart can be added through the defined precision
  and batch capabilities without changing the core solve loop.

