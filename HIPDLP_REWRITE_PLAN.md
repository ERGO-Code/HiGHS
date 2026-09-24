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
- Use the pinned cuOpt 26.08 installation at `/home/yanyu/github/cuopt`
  (tag `v26.08.00`, commit
  `400863c1a92914bd47f5d68245f9a00e2d0d2df4`) through
  `/home/yanyu/github/cuopt-26.08-env`. Phase 7 uses the installed
  `cuopt_cli` with `--method 1 --pdlp-precision double`; it does not use an
  unpinned `main`, a server, or a container.
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

### Tier A4 prerequisite to Phase 7: scalar step-size initialization

`algorithm/step_size_policy.{hpp,cc}` owns only initialization from an
already-computed squared operator norm, reproducing the legacy fixed/PID
formulas: `primal_weight = (1 + unscaled_cost_norm) / (1 + unscaled_rhs_norm)`,
`base_step = safety_factor / sqrt(operator_norm_squared)` with a default
safety factor of `0.998`, `primal_step = base_step / primal_weight`, and
`dual_step = base_step * primal_weight`. It consumes the operator norm as a
scalar, so power iteration was not required and the Backend contract, CUDA
code, and execution path are unchanged. Invalid inputs and every non-finite,
overflowed, underflowed, or non-positive derived step fail transactionally with
zero-initialized step sizes. The focused test was first confirmed to fail at
the missing header; final `[phase7d]` passes 165 assertions in four cases in
Clang CPU, Clang CUDA, and GCC CPU-only builds. Step-size initialization and
restart remain unintegrated, public `solver="hipdlp"` still uses legacy,
`algorithm::solve` still returns structured `kNotReady`, and Phase 7
cross-solver validation has not begun.

### Tier A5 prerequisite to Phase 7: bounded operator-norm estimate

`algorithm/operator_norm.{hpp,cc}` owns only a deterministic, bounded AA'
power iteration over the existing Backend primitives. It starts from an
all-ones row vector, preallocates all three vectors, reports two scalar
readbacks per completed iteration, and returns the squared norm, iteration
count, final relative change, tolerance status, and an explicit empty-matrix
outcome. It remains unintegrated; no Backend contract, backend implementation,
fixed trace, step-size initialization, restart path, solver loop, or public
routing changed. `[phase7e]` passes 208 assertions in Clang and GCC CPU builds
and 220 assertions in the Clang CUDA build, including CPU/CUDA agreement.

Phase 7 triage must reconcile the unchanged legacy `powerMethod()` and
`powerMethodGpu()` divergence: all-ones versus `mt19937(12345)` starts, 20
versus 5000 iterations, and no convergence check versus a `1e-4` residual.
`pdhg.cc` remains unchanged, and Phase 7 cross-solver validation has not begun.

### Tier A6 prerequisite to Phase 7: scalar PID primal-weight update

`algorithm/step_size_policy.{hpp,cc}` gained `updatePrimalWeight()` and
`resetPrimalWeight()`, a backend-neutral scalar model of the active legacy
`updatePrimalWeightAtRestart()`. It consumes only the current step sizes, the
PID state, the primal/dual anchor distances, and the relative primal/dual
residuals; it owns no vector, Backend handle, allocation, reduction,
synchronization, convergence cadence, or restart decision. Defaults are
`k_p=0.99`, `k_i=0.01`, `k_d=0.0`, and `integral_smoothing=0.3`, and the guard
windows are the strict legacy ones: distances in `(1e-16, 1e12)` and a residual
ratio in `(1e-8, 1e8)`, with an undefined ratio represented by the legacy
`1e300` sentinel. Inside the windows the update is
`error = log(dual_distance) - log(primal_distance) - log(primal_weight)`,
`error_sum = integral_smoothing * previous_error_sum + error`,
`delta_error = error - previous_error`, and
`primal_weight *= exp(k_p*error + k_i*error_sum + k_d*delta_error)`. Finite
values outside a window succeed as the legacy fallback, restoring
`best_primal_weight` and zeroing the PID history; the result reports
`kPidUpdate` or `kFallback` explicitly. Best residual balance is tracked with
`abs(log10(relative_dual / relative_primal))` behind an explicit
`best_gap_valid` flag rather than an infinity sentinel, and only a finite gap
from positive residuals can improve it. The geometric mean `eta` of the current
steps is preserved before returning `eta/primal_weight` and
`eta*primal_weight`. Every state, parameter, input, logarithm, exponential,
product, division, weight, and derived step is validated, and any non-finite,
overflowed, underflowed, or non-positive required value fails transactionally
with zero-initialized output.

The update is deliberately unintegrated: `restart_policy`, `fixed_trace`,
`halpern_pdhg_step`, `algorithm::solve`, option translation, and every backend
are unchanged. The TDD falsifier was the missing update API, and `[phase7f]`
then passed 1,088 assertions in five cases in the Clang CPU, Clang CUDA, and
GCC CPU-only builds, with CUDA graph launch counts unchanged at 0/1/1/2/3.

Recorded for Phase 7, not fixed here: `pdhg.cc` maps unsupported adaptive and
Malitsky-Pock `pdlp_step_size_strategy` values onto PID, and its
`performPdhgStep()` PID branch holds a dead `assert(111 == 999)` because the
function is never called and the active Halpern path changes the primal weight
only at restart. `pdhg.cc` and option translation remain unchanged, and Phase 7
cross-solver validation has not begun.

### Tier A7 prerequisite to Phase 7: scalar fixed-point restart trigger

`algorithm/restart_policy.{hpp,cc}` gained `evaluateFixedPointRestart()` and
`resetFixedPointRestart()`, a backend-neutral scalar model of the active legacy
`PDLPSolver::checkRestartCriteria()` in `pdhg.cc`. It consumes only the current
fixed-point error, the reference error for the current restart cycle, the
previous trial error and its validity, the Halpern iterations since the last
restart, and the total iteration count. It owns no vector, Backend handle,
allocation, reduction, synchronization, norm computation, or primal-weight
update. Defaults are `sufficient_decay=0.2`, `necessary_decay=0.8`,
`artificial_threshold=0.36`, and `check_interval=40`.

The rule reproduces the legacy order and comparison strictness:
`total == check_interval` restarts unconditionally;
above the interval, `current <= sufficient_decay * reference`,
then `current <= necessary_decay * reference` with a valid previous trial and
`current > previous_trial`, then
`halpern >= artificial_threshold * total`; otherwise continue. The legacy
inclusive `<=` decay comparisons and strict `>` non-improvement comparison are
preserved, and a necessary-decay window that fails non-improvement still falls
through to the artificial check. The legacy `last_trial_fpe = infinity`
sentinel is replaced by an explicit `previous_trial_valid` flag, so an invalid
previous trial can never satisfy the necessary-decay branch.
`FixedPointRestartCondition` reports which condition fired, including
`kFirstCheckInterval`. A restart returns an invalid reference and previous
trial, so the caller must supply a freshly computed reference through
`resetFixedPointRestart()`; a continue carries the reference and records the
current error as the next previous trial. `artificial_threshold` is
deliberately uncapped so the threshold product stays a reachable guard. Every
parameter, state, input, counter, and product is validated, and any non-finite,
negative, inconsistent, or overflowed value fails transactionally with a
zero-initialized decision, condition, and state.

The trigger is deliberately unintegrated: `fixed_trace`, `halpern_pdhg_step`,
`algorithm::solve`, option translation, `evaluateRestart()`/`RestartState`,
every backend, `pdhg.cc`, and `restart.cc` are unchanged. The TDD falsifier was
the missing trigger API; the first build failed at compile time on the absent
`FixedPointRestart*` names. `[phase7g]` then passed 304 assertions in five
cases in the Clang CPU, Clang CUDA, and GCC CPU-only builds, and the CUDA graph
launch counts remain 0/1/1/2/3 for 0/1/2/4/5 iterations.

Recorded for Phase 7, not fixed here: the Tier A3 `evaluateRestart()` models
`RestartScheme::checkRestart()` in `restart.cc`, which the active Halpern loop
never calls; the two rules differ in metric, comparison strictness, and
reference point. `RestartScheme::updateBeta(primal_weight^2)` ties that
inactive scheme to the legacy primal-weight update. The legacy first restart is
unconditional at `final_iter_count == PDHG_CHECK_INTERVAL`, and
`PDHG_CHECK_INTERVAL` is a compile-time constant rather than an option.

### Tier A8 prerequisite to Phase 7: fixed-point error evaluator

`algorithm/fixed_point_error.{hpp,cc}` now evaluates the active legacy Halpern
fixed-point error through only `copy`, `axpy`, `multiplyTranspose`, `norm2`, and
`dot`. It consumes the current `omega`, the stored `eta`, four immutable
iterate handles, and three distinct caller-owned scratch handles. It allocates
nothing, reports exactly three scalar readbacks on success, and returns a
zero-initialized error value and readback count on any validation or backend
failure. The TDD build first failed on the missing module; final `[phase7h]`
passes 692 assertions in six Clang/GCC CPU cases and 762 assertions in seven
Clang CUDA cases. CPU/CUDA agreement uses the established `5e-13` margin, and
CUDA graph launch counts remain 0/1/1/2/3.

Recorded, not fixed: legacy PID restarts update `params_.omega` without
refreshing `params_.eta`; legacy checks allocate three temporary vectors; and
CPU/GPU reduction order permits tolerant rather than bitwise cross-backend
agreement. The evaluator and Tier A7 trigger remain unintegrated, public
`solver="hipdlp"` still uses legacy, `algorithm::solve` still returns structured
`kNotReady`, and Phase 7 cross-solver validation has not begun.

### Tier A9 prerequisite to Phase 7: final fixed-point-error observation

`algorithm/fixed_trace.{hpp,cc}` now optionally calls the Tier A8 evaluator
exactly once after a positive completed trace and after all recorded replays.
It forwards the supplied `omega` and `eta`, reuses the retired next-iterate
handles and primal candidate as scratch, allocates no vector, and publishes the
three scalar readbacks only on success. The default is disabled with invalid
scalar sentinels; enabled zero-iteration traces validate those scalars but stay
unobserved with no readback. Any evaluator failure discards the trace,
averages, convergence observations, fixed-point observation, and partial
scalar diagnostics. `[phase7i]` passes in Clang CPU/CUDA and GCC CPU, with
CPU/CUDA agreement and unchanged recorded launches 0/1/1/2/3.

Recorded, not fixed: legacy `initial_fpe_` starts at zero and the first
iteration-40 check restarts before a later major step refreshes the reference;
legacy convergence precedes the restart trigger. Tier A9 owns no reference
error or precedence, does not invoke Tier A7, and leaves restart execution,
PID updates, cadence, final selection, public routing, and Phase 7 validation
deferred.

### Tier A10 prerequisite to Phase 7: final fixed-point restart observation

`algorithm/fixed_trace.{hpp,cc}` now composes the Tier A9 fixed-point-error
observation with the Tier A7 scalar trigger exactly once, as an opt-in
observation reported as data. `FixedTraceOptions::fixed_point_restart` carries
an `enabled` flag, the caller-owned previous `FixedPointRestartState`, and the
`FixedPointRestartParameters`; `FixedTraceResult::fixed_point_restart` exposes
the unchanged `FixedPointRestartResult`, including its next state, decision,
and condition. The default is disabled, and enabling the trigger without the
fixed-point-error observation is rejected with `kInvalidOption` before any
backend resource is created.

For a positive trace the trigger runs only after the final convergence
observations succeed. It forwards the Tier A9 error unchanged as
`current_error`, passes the completed iteration count as both
`halpern_iterations` and `total_iterations`, and forwards `previous` and
`parameters` exactly as supplied. The trace neither applies the reported next
state nor synthesizes, refreshes, or retains reference-error state, and it does
not call `resetFixedPointRestart()`. A zero-iteration trace has no
fixed-point-error observation and therefore no trigger observation. A
`kRestart` decision is still a successful trace that publishes the normal
points, averages, and convergence observations; an optimal convergence
decision and a restart decision are exposed independently without choosing
between them. Any trigger failure preserves the first Tier A7 error code and
message and discards the entire trace, including the successful Tier A9
observation.

The trigger is scalar-only: relative to an otherwise identical Tier A9 trace it
adds no allocation, descriptor, backend operation, transfer, scalar readback,
synchronization, recording, or graph launch. No production module, test
executable, backend capability, or build registration was added, and no
restart execution, `updatePrimalWeight()`, KKT verifier, postsolve, logging, or
legacy call was introduced. The TDD falsifier was the absent
`FixedTraceFixedPointRestart*` API; the first build failed at compile time.
`[phase7j]` then passed 569 assertions in seven cases in the Clang CPU and GCC
CPU builds and 1801 assertions in nine cases in the Clang CUDA build, with
serial CPU, direct CUDA, and recorded CUDA agreeing, identical CUDA resource
audits between the Tier A9-only and Tier A10 runs, and recorded launches still
0/1/1/2/3 for 0/1/2/4/5 iterations.

Recorded, not fixed: the two-iteration fixed-point error of the one-step
fixture is exactly `sqrt(289/576) = 17/24` in exact arithmetic but lands one
unit in the last place below it in FP64, because the Halpern anchor weight
`1 - 1/3` is not the correctly rounded `2/3`. Restart execution, reset, PID
updates, reference ownership, periodic cadence, final iterate selection,
operator-norm and step-size initialization, public routing, and Phase 7
cross-solver validation remain unimplemented.

### Tier A11 prerequisite to Phase 7: periodic fixed-point restart observation

`algorithm/fixed_trace.{hpp,cc}` now observes the Tier A7 trigger on a
caller-controlled cadence instead of only once.
`FixedTraceFixedPointRestartOptions::check_cadence` defaults to `0`, which keeps
the Tier A10 final-only check; a positive cadence adds a check at every multiple
of it below the final iteration. `FixedTraceResult` exposes an ordered
`fixed_point_restart_checks` of `FixedTraceFixedPointRestartCheck{iteration,
error, evaluation}`, and `fixed_point_error` and `fixed_point_restart` remain
the last check, so Tier A9 and Tier A10 outputs are unchanged. A negative
cadence is rejected with `kInvalidOption`, as is a cadence that is not a
multiple of the recorded block length while blocks would still be pending,
because recorded blocks are never split, resized, or reordered and no
conditional recorded node was added.

Each check evaluates the Tier A8 error exactly once, reusing only the scratch
handles Tier A9 already reuses: the retired next-iterate slots and the primal
candidate. All three are pure outputs of the Halpern step, so a mid-loop check
cannot perturb a later iteration. Direct traces check at an iteration boundary
after the point is recorded; recorded traces check between block replays using
the handle assignment snapshotted when that block was recorded. Reference state
is chained explicitly and derived from the recorded checks alone: `kContinue`
carries `result.next`, `kRestart` carries `resetFixedPointRestart()` of the
error at that same check, `halpern_iterations` measures the distance to the last
`kRestart` check, and `total_iterations` is the completed iteration count. No
vector, anchor, average, step size, primal weight, counter, or trace result is
altered, and no restart execution, `updatePrimalWeight()`, KKT verifier,
postsolve, logging, or legacy call was introduced. A failure at any check
preserves the first error code and message and discards the whole trace,
including earlier successful checks.

No production module, test executable, backend capability, or build
registration was added. The TDD falsifier was the absent cadence API; the first
build failed at compile time on `check_cadence`,
`FixedTraceFixedPointRestartCheck`, and `fixed_point_restart_checks`.
`[phase7k]` then passed 609 assertions in six cases in the Clang CPU and GCC CPU
builds and 2194 assertions in seven cases in the Clang CUDA build. Cadence `0`
reproduces Tier A10 bitwise, including diagnostics and identical CUDA resource
audits; serial CPU, direct CUDA, and recorded CUDA agree on every decision and
condition; recorded launches remain 0/1/1/2/3 for 0/1/2/4/5 iterations; repeated
runs are bitwise deterministic per backend; and the CUDA audit delta per extra
check is exactly six device-to-host events and three synchronizations, with no
extra allocation, descriptor, host-to-device transfer, graph creation, or graph
launch.

Recorded, not fixed: the cadence is a pure observation sequence with no
consumer. Restart execution, anchor and averaging reset, the full solve loop,
final iterate selection, PID primal-weight invocation, operator-norm and
step-size initialization, adaptive line search, Malitsky-Pock behavior,
certificates, postsolve, status mapping, public rewrite routing, and Phase 7
cross-solver validation all remain unimplemented.

### Tier A12 prerequisite to Phase 7: scalar restart primal-weight update

`algorithm/step_size_policy.{hpp,cc}` now exposes
`updatePrimalWeightAtRestart()` as a pure scalar extension of the existing Tier
A6 policy. The active legacy body reads `params_.i_smooth`, `params_.k_p`,
`params_.k_i`, and `params_.k_d`; the current and best weights, PID sum and last
error, best residual gap, and incoming primal/dual steps. Its admissibility
order is the strict distance window `(1e-16, 1e12)` followed by the strict
relative-infeasibility-ratio window `(1e-8, 1e8)`. A zero relative primal
residual selects the legacy `1e300` sentinel, which always fails the upper
bound and therefore reverts. The controlled PID evaluation order, fallback to
the best weight with both accumulators zeroed, strict post-branch best-gap
tracking, and incoming `sqrt(primal_step * dual_step)` eta invariant are
preserved.

Every input scalar and intermediate is checked for FP64 range; any failure
returns the first error with `kNone`, zero state, and zero step sizes. Unlike
the legacy tail, this routine does not write beta or omega and does not call
the inactive `RestartScheme::updateBeta()`. It owns no vector, Backend,
reduction, allocation, logging, or legacy dependency, and only `[phase7l]`
tests call it. The requested TDD compile failure was observed before adding the
API. Final Tier A aggregates pass 5,656 assertions in 77 cases in Clang CPU
and GCC CPU, and 10,053 assertions in 85 cases in Clang CUDA; all requested
`[phase3c]`, `[phase5g]`, `[hipdlp]`, and `[pdlp]` Clang filters also pass.
Restart execution, reset application, solve-loop integration, public routing,
and Phase 7 cross-solver validation remain unimplemented.

### Tier A13 correction: candidate-based fixed-trace error observation

`algorithm/fixed_trace.cc` now passes the projected primal and dual candidates,
rather than the post-Halpern current iterate, to every Tier A8 fixed-point-error
evaluation. Reflected handles remain immutable inputs, the retired iterate
slots remain delta scratch, and the already allocated `averaged_aty` handle is
reused as distinct transpose-product scratch before its later overwrite. This
adds no resource, transfer, synchronization, readback, capability, or graph
operation; each evaluation still reports exactly three scalar readbacks.

The `[phase7m]` TDD fixture first observed `17/24` instead of the independently
derived candidate value `11/16` at iteration two. It now passes 226 assertions
in four Clang/GCC CPU cases and 904 assertions in five Clang CUDA cases,
covering exact one- and multi-step arithmetic, cadence boundaries, direct and
recorded agreement, audit and graph-launch invariants, determinism, and
transactional failure. Affected Tier A10-A11 expected values now follow the
legacy candidate formula. Restart decisions remain observations only: no
restart, anchor or average reset, primal-weight update, solve loop, final
iterate selection, postsolve/status handling, or public routing is consumed or
changed.

### Tier A14 prerequisite to Phase 7: fixed-trace restart state transition

`algorithm/fixed_trace.cc` now consumes the Tier A11 decisions instead of only
recording them. On `kRestart` the trace copies the projected primal and dual
candidates into both the current iterate and both anchors, recomputes the
cached `Ax` and `A^T y`, resets the Halpern age, resets `AveragingState`, and
clears both weighted sums, exactly as the active legacy
`PDLPSolver::solve()` restart body does. It calls neither
`updatePrimalWeightAtRestart()` nor any step-size update.

Each check is evaluated and pushed before any iterate, anchor, age, average, or
sum is touched, so the record keeps the pre-reset candidate error and decision
while the boundary trace point is rewritten to the post-decision state that
continues into the next iteration. The same transition runs at a final
restart, which leaves no post-restart sample and therefore no valid average or
averaged convergence observation. `HalpernPdhgState::iteration` is now used
only as the resettable Halpern age; the trace owns a separate monotone
`total_iterations` that drives the loop bound, trace indices, cadence
boundaries, step-sequence indices, recorded check iterations, and the trigger's
`total_iterations` input, so no header or shared-step API changed.

Recorded mode now records, finalizes, and replays each existing one- or
two-step block before capturing the next one, because a captured anchor weight
bakes in the Halpern age that the block's own decision may reset. Blocks are
neither split nor reordered, no conditional graph node was added, and the
mid-block-cadence rejection is unchanged; command-buffer, graph-capture,
instantiation, and launch counts stay at 0/1/1/2/3 for 0/1/2/4/5 iterations.
An executed restart adds only device-to-device copies and direct kernel
enqueues: no allocation, host-to-device transfer, download, scalar readback,
explicit synchronization, capability, module, executable, or build
registration. Every fixed-point-error evaluation still reports exactly three
scalar readbacks, and disabled or never-triggering restart traces stay bitwise
identical to fixed-point-error-only traces with identical diagnostics and CUDA
audits.

The `[phase7n]` TDD fixture failed first against the observing implementation:
at a restart on iteration two it saw the Halpern blend `7/6` and the follow-on
`29/16` instead of the candidate `9/8` and the post-restart `57/32`, and the
third error `25/48` instead of `35/64`. `[phase7n]` now passes 423 assertions
in nine cases in Clang CPU and GCC CPU and 2,065 assertions in eleven cases in
Clang CUDA, covering exact restart-at-two arithmetic and the following step,
candidate-to-current and candidate-to-anchor adoption, divergent total and
Halpern counters, averaging reset with only later samples, mixed periodic
decisions, final-boundary semantics, direct CPU, direct CUDA and recorded CUDA
agreement, resource and graph-launch invariants, determinism, and transactional
failure including injected CUDA failures. Affected Tier A9-A11 expectations
now follow the executed trajectory. Final Clang totals are `[hipdlp]`
9,424/127 and `[pdlp]` 9,454/135 on CPU and 24,159/159 and 24,189/167 on CUDA,
with `[phase3c]` and `[phase5g]` unchanged.

Recorded, not fixed: restart primal-weight updates stay unconsumed. Operator
norm and step-size initialization, convergence cadence and early stopping,
final iterate selection, the solve loop, postsolve and status handling, public
rewrite routing, and Phase 7 cross-solver validation all remain unimplemented.

### Tier A15 prerequisite to Phase 7: fixed-trace convergence cadence

`algorithm/fixed_trace.{hpp,cc}` now optionally records Tier A1 convergence
checks from the projected primal and dual candidates at a caller-selected
cadence. The two candidate products are preallocated only when enabled; each
check recomputes both products, downloads exactly the four candidate vectors,
passes the monotone boundary total to `evaluateStopping()`, and is stored
immediately before the restart check. Recorded mode observes only after the
same block has replayed. Decisions are diagnostic only: even `kOptimal` does
not stop or shorten the fixed trace, and Tier A14 restart transitions remain
unchanged.

The one-column equality fixture gives raw
`(primal, dual, gap, primal objective, dual objective)` values
`(3/2, 0, 3/2, -1/2, -2)` at iteration one and
`(7/8, 0, 7/8, -9/8, -2)` at iteration two, normalized to
`(1/2, 0, 3/7)` and `(7/24, 0, 7/33)`. The first Halpern point is
algebraically the candidate; at iteration two its metrics instead use
`(5/6, 0, 5/6, -7/6, -2)` and normalize to `(5/18, 0, 1/5)`. A restart at
the shared boundary leaves the recorded candidate metrics unchanged while
adopting that candidate as the Tier A14 boundary point and resetting the later
trajectory.

TDD first failed to compile because the periodic API and ordered history did
not exist. Final `[phase7o]` passes 771 assertions in five Clang/GCC CPU cases
and 1,988 assertions in seven Clang CUDA cases. Enabled runs add exactly two
pre-loop resources and four diagnostic transfers per check; CUDA adds eight
device-to-host events and four synchronizations per check, with unchanged
command-buffer and graph-launch counts. All 53 requested post-edit filters
pass separately, including `[phase7a]` through `[phase7o]`, `[phase3c]`,
`[phase5g]`, `[hipdlp]`, and `[pdlp]` in both Clang builds and all Tier A
filters in GCC. Early stopping, restart primal-weight updates, initialization,
final selection, the solve loop, postsolve/status handling, public routing,
and Phase 7 remain unconsumed.

### Tier A16 prerequisite to Phase 7: restart primal-weight consumption

`algorithm/fixed_trace.{hpp,cc}` now consumes the Tier A12
`updatePrimalWeightAtRestart()` policy, but only where the Tier A11 decision at
a shared boundary is `kRestart`. The boundary order is the Tier A15 candidate
convergence observation, the Tier A11 evaluation and record, the two
candidate-to-anchor distances taken through `backend::differenceNorm2()` before
either anchor is overwritten, the policy call, installation of the returned
primal and dual steps for later iterations only, and finally the unchanged Tier
A14 transition. No iterate is downloaded and no vector is added; recorded mode
performs the reductions and the policy evaluation after the block has replayed
and before the next block is captured. A disabled-by-default option carries the
initial `PrimalWeightState` and `PrimalWeightParameters`, and an ordered history
records the iteration, both distances, and the `PrimalWeightUpdateResult`. A
`kContinue` boundary updates nothing, and a final-boundary restart records its
update even though no later iteration consumes the returned steps.

Enabled updates are rejected before any backend resource when the fixed-point
restart trigger is off, when a possible restart boundary would lack a Tier A15
observation, or when caller-provided step-size sequences are present. The
boundary convergence and restart records, the Tier A14 candidate adoption and
averaging reset, the fixed-point-error values with exactly three scalar
readbacks per evaluation, trace point count and cadence indexing, disabled-mode
bitwise identity, and command-buffer and graph-launch counts are all preserved.
Each executed update costs exactly two difference-norm reductions: two direct
kernel enqueues, four device-to-host events, and two synchronizations, with no
allocation, vector creation, host-to-device transfer, or captured primitive.

TDD first failed to compile because the option, the history, and the policy
include did not exist. On the restart-at-two fixture the pre-reset distances are
exactly `9/8` and `5/16`; the bounded column has relative residuals
`(7/24, 0)`, whose undefined ratio reverts to the best weight, while a free
column gives `(7/24, 21/32)` and ratio `9/4`, which drives the controlled PID
branch. Both keep `eta = sqrt(tau * sigma)`. An independent oracle that replays
the installed steps as a caller sequence reproduces the enabled trace bitwise.
`[phase7p]` passes 949 assertions in nine Clang and GCC CPU cases and 2,174
assertions in eleven Clang CUDA cases. All requested post-edit filters pass:
Clang CPU `[hipdlp]` 11,209/141 and `[pdlp]` 11,239/149; Clang CUDA
`[hipdlp]` 28,428/177 and `[pdlp]` 28,458/185; `[phase3c]` and `[phase5g]`
unchanged in behavior.

Recorded, not fixed: operator-norm and initial step-size wiring, trace
termination or truncation, final iterate selection, the solve loop,
postsolve and status handling, public rewrite routing, and Phase 7
cross-solver validation all remain unimplemented.

### Cutover readiness: public options, CPU solve completeness, CUDA and auto

Three sequential checkpoints made the rewrite usable from ordinary HiGHS-style
options on CPU, direct CUDA, and auto, while `solver = "hipdlp"` still routes to
legacy by default. Nothing about the public default routing changed and no
legacy file was moved or removed.

Checkpoint 1 added the public options `hipdlp_backend`
(`auto`/`cpu`/`cuda`/`openmp`/`rocm`, default `auto`), `hipdlp_scaling`
(`auto`/`off`/`curtis_reid`/`legacy`, default `auto`), and
`hipdlp_curtis_reid_iterations` (positive integer, default 20), each with a
field, an option record, string validation, and a documentation entry. Option
translation now produces a strongly typed `EffectivePolicy`, and
`SolverOptions` no longer carries raw `HighsInt` restart or step-size
strategies. The requested-to-effective mapping is: `fixed` to fixed; `PID` to
PID; `adaptive`, including the HiGHS default, to fixed with one explicit
warning and never to PID as the legacy solver does; `Malitsky-Pock` to a
structured `kInvalidOption` error; every restart request and the restart-off
feature bit to Halpern with one warning, matching the legacy public message.
The policy resolution is applied only to the rewrite selector, so the legacy
public path keeps its own semantics. A non-`auto` `MigrationOptions` backend
remains the test and migration override; otherwise the public option decides.
TDD RED was observed against the not-yet-rebuilt binaries: `bin/highs` rejected
an options file containing `hipdlp_backend` with `getOptionIndex: Option
"hipdlp_backend" is unknown` and exit code 255, and only one `[phase7][options]`
test existed. GREEN: `[phase7][options]` 165 assertions in three cases and
`[migration]` 193 in two, with `[phase7][solver]` unchanged at 98 in eight.

Checkpoint 2 completed the CPU production loop inside the single
`algorithm::solve()` control loop. An iteration-0 convergence evaluation now
runs on the initialized current iterate, and an initial point that passes it
still has to survive postsolve and the independent KKT verifier before it is
reported optimal. Remaining time is checked immediately before and after the
operator-norm power iteration; no stop predicate was added inside the power
iteration, so a time limit can still be observed only at those boundaries.
`updatePrimalWeightAtRestart()` is consumed only when the fixed-point restart
decision is `kRestart` and the effective policy is PID: both candidate-to-anchor
distances are reduced before any anchor moves, the normalized primal and dual
residuals already computed at that boundary are reused, the returned steps are
installed for later iterations only, and a policy failure returns the first
`api::Error` before any state is mutated. The fixed-point error now uses the
current omega with the unchanged initialized base step. `api::solve()` passes
the rewrite only what is left of the global budget, `time_limit - timer.read()`
clamped at zero, and preserves the infinite sentinel. A small
`algorithm::SolveContext` carries an injectable log sink and self-audit
diagnostics; the HiGHS adapter reports backend, requested and effective
step/restart/scaling, and the final status, iteration count, elapsed time and
relative KKT metrics exactly once, and emits nothing when `output_flag` is
false. Zero-operator LPs keep the explicit `kInvalidModel` precondition, and
tests confirm the public `Highs` path routes no-row and zero-matrix models to
`solveUnconstrainedLp` without entering PDHG.

Checkpoint 3 removed the CPU-only restriction without adding a backend branch.
`algorithm::solve()` validates only that the kind is valid and hands it to the
factory, so explicit CPU always uses CPU, explicit CUDA runs the same direct
loop when built and probed, explicit CUDA returns `kBackendUnavailable` with no
fallback otherwise, `auto` prefers CUDA and otherwise uses CPU, and OpenMP and
ROCm stay structurally unavailable. Backend names for logging come from the
factory registry, so the algorithm layer still contains no CUDA or ROCm symbol.
No CUDA Graph or command-buffer path is used. The Phase 5H and Phase 6
structural source scans were narrowed so that `api/options.cc`, `api/options.hpp`
and `HighsOptions.h` may name the reserved public backend strings while the bans
on CUDA and HIP headers, types, graph symbols and implementation names remain in
force for those files and unchanged everywhere else.

Verification used the existing builds. Clang CPU
`/home/yanyu/builds/hipdlp2-m7b-cpu-20260826`: `[phase7][options]` 185/4,
`[phase7][solver]` 253/17, `[migration]` 193/2, `[phase7]` 406/19,
`[phase7a]` 569/11, `[phase7f]` 1088/5, `[phase7p]` 949/9, `[phase5h]` 180/1,
`[phase6]` 1301/4, `[hipdlp]` 11651/160, `[pdlp]` 11681/168. Clang CUDA
`/home/yanyu/builds/hipdlp2-m7b-cuda-20260826`: `[phase7][solver][cuda]` 61/1,
`[phase7][solver]` 280/17, `[phase7]` 433/19, `[phase7a]` 1348/12,
`[phase7p]` 2174/11, `[phase5h]` 377/2, `[phase6]` 1301/4, `[hipdlp]`
28897/196, `[pdlp]` 28927/204. With `CUDA_VISIBLE_DEVICES` empty the same CUDA
binary passes `[phase7][solver]` 253/17 and `[phase7][solver][cuda]` 34/1,
taking the unavailable and CPU-selected branches. The existing CUDA
ownership/lifecycle regression under Compute Sanitizer passes 29/1 with zero
errors and zero leaked bytes; no CUDA backend resource code was changed, so no
new sanitizer coverage was added. GCC CPU
`/home/yanyu/builds/hipdlp2-tier-a-gcc-cpu-20260827` matches Clang CPU on every
focused filter; its only `[hipdlp]`/`[pdlp]` failure is the known
compiler-sensitive legacy `hi-pdlp-scrs8-legacy-characterization` trajectory
difference, which is reported and not modified.

Not executed: Bazel, Bazelisk, Meson and buildifier are unavailable, so their
registrations were only checked statically; this task added no source or test
file, so no registration changed.

Recorded, not fixed: `defaultMigrationOptions()` still selects legacy and
`solver = "hipdlp"` still routes to legacy; infeasibility and unboundedness
certificates, the adaptive and Malitsky-Pock algorithms, averaged-iterate final
selection, trivial zero-operator solving inside the rewrite, CUDA Graph solving,
OpenMP and ROCm, performance tuning, mixed precision, multistart, Phase 7
cross-solver validation and the performance gates, and legacy deletion all
remain unimplemented. No performance claim is made.

### Phase 7: Full cross-solver validation

Execute this phase through `HIPDLP_PHASE7_IMPLEMENTATION_PLAN.md`.

- Run rewritten CPU and CUDA against legacy HiPDLP, in-tree cuPDLP, external
  cuPDLPx, pinned cuOpt, and HiGHS reference solvers.
- Verify every claimed optimal result independently.
- Triage status, objective, or iteration-count disagreements.
- Test time/iteration limits and failure propagation.
- Produce machine-readable pass/fail summaries and a human-readable report.

**Exit criteria:** no unexplained incorrect status, non-finite result, or KKT
failure remains in the supported suite.

**Completion audit (2026-09-04):** Phase 7 is complete. Provenance-stable smoke
passes rewrite CPU and CUDA 13/13. The eight-solver Netlib regression, full
Netlib, energy regression, and full-energy matrices contain 824, 912, 32, and
72 canonical rows, respectively. Full Netlib accepts rewrite CPU on 111/114
and CUDA on 107/114; full energy accepts each on 6/9. Every rewrite `OPTIMAL`
passes the independent original-LP `1e-4` KKT gate, all rewrite rejections are
time or process limits, and all 626 disagreement rows have non-`UNRESOLVED`
classifications. No new rewrite defect was found. The canonical summary,
directory audit, resource evidence, classifications, and human report are under
the Phase 7 result root and `docs/hipdlp-phase7-cross-solver-report.md`. No
legacy CUDA comparison is claimed. Phase 7 clears the correctness prerequisite
for Phase 8, while the separate performance and stabilization gates still
govern legacy removal.

### Phase 8: Cut over and remove legacy code

- [x] Switch public `solver = "hipdlp"` to the rewrite (already completed before Phase 8).
- Retain the legacy selector for one stabilization window if maintenance cost
  is acceptable.
- Confirm CMake CPU/CUDA configurations, unavailable-placeholder behavior, and
  Bazel/Meson CPU builds.
- Update user and developer documentation.
- Remove legacy solver, duplicated CUDA state, obsolete macros, and temporary
  differential hooks after the stabilization window.

**Exit criteria:** all Definition of Done items below pass from a clean build,
and the legacy implementation is no longer needed to diagnose routine failures.

**Historical four-case checkpoint (2026-09-08; superseded denominator):
BLOCKED, not COMPLETE.** Production
GPU checks now use a backend-neutral 112-byte result packet, conservative
certificate screening, transactional prepared downloads and the unchanged
host policy/independent verifier. The original rewrite and distinct true legacy
CPU/CUDA baselines are frozen outside the source tree. Matched production
builds and twelve original-model correctness preflights pass. The representative
one-warmup/five-measurement experiment executed all 24 jobs and 144 records;
strict revalidation accepts 126 and rejects 18 legacy results. CPU has only
2/4 accepted legacy pairs and CUDA 3/4; their diagnostic accepted-subset ratios
are 3.9903 and 5.2256, not valid substitutes for the full four-pair gate. Both
performance gates are NOT_PASS. Candidate/preoptimization full-solve ratios
are 1.0269 CPU and 1.0213 CUDA, so this batch does not establish full-solve
acceleration. Current CMake Clang/GCC, Bazel/Meson CPU focused checks pass
175 cases each; CUDA passes 236 cases and memcheck is clean. The strict process-wide
allocation audit still sees recurring vendor worker allocations, and the large
energy fixed-workload CPU/CUDA trajectory check fails its unchanged `1e-9`
threshold in both before/after CUDA. Details, commands, timing boundaries and
remaining gates are in `docs/hipdlp-phase8-gpu-performance-report.md`. Legacy
remains present until the independent removal conditions are met.

**Legacy-qualified rerun (2026-09-08): supersedes the four-case denominator.**
The available Phase 7 data contained 30 unique full-Netlib legacy cases that
passed independent verification within 60 seconds and zero eligible energy
cases, so an honest common set could not reach the requested approximate 50.
Fresh frozen legacy CPU/CUDA qualification at uniform `1e-4` retained their
strict 23-case intersection before candidate execution. The formal run used
all 23 cases, PID, FP64, one warmup and five measurements for 552/552 records.
Legacy CPU/CUDA each passed 138/138 records; candidate CPU passed 132/138 and
CUDA 120/138. CPU `fit2d` and CUDA `fit1d`, `fit2d`, `greenbeb` reached the
60-second limit. Thus CPU has 22/23 and CUDA 20/23 complete pairs; formal GM/CI
are N/A by the no-denominator-shrink rule and both gates remain NOT_PASS.
Accepted-subset diagnostic GMs (not gates) are 2.3907 CPU and 2.1477 CUDA, with
34 ratios over 1.20. Evidence is under
`phase8/20260908T205450Z-legacy-qualified-1e4`; the previous four-case result is
historical only. Legacy removal remains blocked.

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


## 14. Phase 4 certificates and Phase 5 cutover evaluation (2026-09-01)

### 14.1 Phase 4 — implemented and verified

Certificate mathematics live in `highs/pdlp/hipdlp/verification/certificate_verifier.{hpp,cc}`,
separate from `kkt_verifier`, which still only validates complete primal-dual
iterates. For the HiGHS form `l_x <= x <= u_x`, `l_r <= A x <= u_r`:

- a primal ray `d` must lie in the recession cone of both bound systems and give
  `sigma * c^T d < 0`, with `sigma = +1` for minimize and `-1` for maximize;
- a Farkas certificate `(y, z)` must satisfy `z = -A^T y`, the row and column
  multiplier sign conditions, and a strictly positive bound objective
  `sum_i y_i b_i(y_i) + sum_j z_j g_j(z_j)`.

The verifier recomputes `A^T y` itself and reports the stationarity residual, so
a caller cannot certify infeasibility by supplying a self-consistent pair.
Violations, residuals and margins are reported raw and normalized by the
direction norm, together with the worst row and column.

`api::SolverStatus` gained `kUnboundedOrInfeasible`, and solver-neutral
`CertificateKind`/`CertificateMetrics`/`Certificate` were added to
`api::SolveResult`. `HighsSolution` was not changed and no rewrite certificate
is written to `HEkk` ray state. `model::Problem` gained explicit
`recoverPrimalDirection`, `recoverDualDirection` and
`recoverReducedCostDirection`; these are not the solution postsolve, because a
direction must not carry the objective offset and must not carry the
maximization sign flip that `recoverOriginalDual` applies.

Detection runs on the single `algorithm::solve()` convergence boundary for both
CPU and CUDA, using only the host candidate already downloaded there. Precedence
is numerical error, then verified optimal, then verified certificate, then
iteration/time limit, then continue. A rejected candidate never terminates the
solve. Certificate statuses still return the last complete postsolved iterate
and log the certificate kind, residual and margin instead of optimality metrics.

`Highs::getDualRay`/`getPrimalRay` now save the caller's `solver`, `presolve`,
`solve_relaxation`, `allow_unbounded_or_infeasible` and objective, force
`solver = simplex` for the auxiliary solve so no first-order solver is re-entered,
and restore everything afterwards; an early return that skipped restoration was
removed.

Verified: Clang CPU full suite 1,271,637 assertions / 545 cases; Clang CUDA
`[hipdlp]` 29,257/228 and `[pdlp]` 29,287/236; GCC CPU focused matrix; CUDA
hidden-device filters; Compute Sanitizer lifecycle 29/1 with zero errors and
zero leaked bytes; Bazel `//:highs` plus `TestHiPdlpCertificates`,
`TestHiPdlpSolver`, `TestRays`; Meson CPU library and app build.

### 14.2 Phase 5 — cutover not performed in the 2026-09-01 checkpoint

`defaultMigrationOptions()` still selects `kLegacy`.

A non-installed migration/benchmark adapter was added
(`check/HiPdlpMigrationAdapter.cpp`, target `hipdlp_migration_adapter`). It reads
one MPS, starts the solver timer after the read, selects legacy or rewrite
through the internal `MigrationOptions`, and emits stable JSON. No public
`hipdlp_implementation` option was added. Legacy only accepts the `auto` backend
selector, so the adapter records both the requested and effective backend.

Smoke correctness on the existing 13-instance set (CPU, `presolve=off`,
`kkt_tolerance=1e-4`, `time_limit=60`) is stored under
`/home/yanyu/benchmark-results/lp/hipdlp-rewrite/m8-smoke-cpu-*`. The rewrite
reaches independently KKT-verified `OPTIMAL` on 11 of 13 instances, including
`scrs8`, where legacy does not terminate at all. Two instances fail:

- `e226`: rewrite hits the 60 s limit after about 1.1e7 iterations while legacy
  returns optimal in 0.067 s and 15,880 iterations;
- `sctest`: rewrite hits the 60 s limit with a relative primal residual of about
  1.0, so it does not converge.

Legacy `adlittle` and `avgas` report `OPTIMAL` but fail the independent `1e-4`
KKT gate, which is recorded as a legacy defect rather than a rewrite one.

The first failing cutover gate is therefore Checkpoint 5B smoke correctness.
Independently, the cutover is also blocked because cuOpt is not installed and
cannot be pinned here, so the full cross-solver suite cannot run; because there
is no frozen legacy CUDA performance baseline; and because no paired
five-repetition timing study was run, so
`geomean(rewrite / legacy) <= 1.05` is not established. Meson test targets could
not be configured because `catch2` is not available.

Still open: legacy deletion, the adaptive and Malitsky-Pock step sizes, CUDA
Graph production solve, OpenMP/ROCm, and the Tier B research items.

### 14.3 Cutover and cuOpt follow-up (2026-09-02)

The current working tree supersedes the cutover conclusion in Section 14.2:

- `defaultMigrationOptions()` now selects the rewrite with backend `auto`, so
  the public `solver = "hipdlp"` path enters the rewrite. The internal legacy
  selector remains available for Phase 7 comparisons and stabilization.
- The migration adapter accepts `--step-size fixed|pid` and forwards the
  selected strategy through the normal option translation path.
- On the CUDA build at
  `/home/yanyu/builds/hipdlp2-m7b-cuda-20260826`, explicit rewrite CUDA PID
  runs with `presolve=off`, `hipdlp_scaling=auto`, a 60-second limit, and
  `kkt_tolerance=1e-4` solve both former blockers. `e226` reports `OPTIMAL` in
  8,440 iterations and 1.18 seconds; its relative primal, dual, and gap metrics
  are `4.34e-7`, `1.16e-9`, and `3.06e-8`. `sctest` reports `OPTIMAL` in 480
  iterations and 0.27 seconds; its metrics are `1.19e-7`, `8.94e-8`, and
  `1.12e-8`. Both pass the independent KKT verifier.
- NVIDIA cuOpt 26.08 is installed without root privileges. The source checkout
  is `/home/yanyu/github/cuopt` at tag `v26.08.00` and commit
  `400863c1a92914bd47f5d68245f9a00e2d0d2df4`; the micromamba environment is
  `/home/yanyu/github/cuopt-26.08-env`. It contains Python 3.12.14,
  `cuopt=26.08.00`, `libcuopt=26.08.00`, and `cuda-version=13.0` for
  `linux-aarch64`. The environment stores `PYTHONNOUSERSITE=1` so packages in
  `~/.local` cannot override its dependencies.
- A forced FP64 GPU PDLP run of cuOpt's `afiro_original.mps` with
  `--method 1`, `--pdlp-precision double`, `presolve=0`, and `1e-4` absolute
  and relative tolerances reports `Optimal` in 460 iterations and 0.192
  seconds on the NVIDIA GB10. The relative primal, dual, and gap metrics are
  `4.27e-6`, `1.60e-5`, and `8.26e-5`. The documented Python LP smoke test also
  reports `Optimal` with objective and primal value 10.

This checkpoint supplied the prerequisites for Phase 7. The validation campaign
was subsequently completed on 2026-09-04 as recorded in the Phase 7 section and
human report. Results remain identified by their fixed executable, manifest,
configuration, and run metadata rather than repository `HEAD` alone.
