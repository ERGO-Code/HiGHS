# HiPDLP Rewrite Implementation Plan

## 1. Objective

Rewrite HiPDLP as a modular, high-performance solver that remains integrated
with HiGHS while separating:

- HiGHS model and option handling;
- problem transformation and scaling;
- the PDHG/Halpern algorithm;
- restart, step-size, and convergence policies;
- linear algebra and device-memory management;
- CPU, NVIDIA CUDA, and AMD ROCm implementations;
- correctness verification, profiling, and benchmarking.

The first complete rewrite will preserve the current Halpern-PDHG behavior and
the public `solver = "hipdlp"` contract. It will provide deterministic reference
CPU, optimized CPU, CUDA, and compile-validated ROCm backends. Low-precision GPU
SpMV and multistart warm-up are follow-on research projects, but the initial
interfaces must accommodate them without another solver rewrite.

## 2. Confirmed Decisions

| Area | Decision |
| --- | --- |
| Core ownership | Keep the solver coupled to HiGHS and its model/types |
| Backend selection | Conditionally compile any subset of backends and select one through a runtime factory |
| Initial backends | Deterministic reference CPU, optimized OpenMP CPU, CUDA, and ROCm |
| Initial algorithm | Reproduce current Halpern-PDHG behavior, then expose research policies |
| Compatibility | Preserve solver name, relevant options, statuses, solution fields, logging, and iteration reporting |
| Migration | Run legacy and rewritten implementations side by side until cutover |
| Correctness | Independently recomputed KKT metrics plus consensus with in-tree cuPDLP, external cuPDLPx, cuOpt, and HiGHS reference solvers |
| Accuracy target | Relative primal residual, dual residual, and objective gap at or below `1e-4` |
| Performance | Rewritten CPU and CUDA geomean solver time must be within 5% of their frozen legacy HiPDLP baselines |
| CPU optimization | Owned CSR/CSC kernels using OpenMP and compiler vectorization; optional BLAS only where it is beneficial |
| Build systems | Full backend support in CMake; keep Bazel and Meson CPU builds working |
| AMD scope | Implement and compile-test ROCm first; hardware correctness and performance validation are deferred |
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
      reference/
        reference_backend.{hpp,cc}
      cpu/
        cpu_backend.{hpp,cc}
        sparse_kernels.{hpp,cc}
        vector_kernels.{hpp,cc}
      cuda/
        cuda_backend.{hpp,cu}
        cuda_kernels.{cuh,cu}
      rocm/
        rocm_backend.{hpp,hip}
        rocm_kernels.{hpp,hip}
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
hipdlp_backend = auto | reference | cpu | cuda | rocm
```

Rules:

- `auto` selects CUDA or ROCm only when that backend was built and a compatible
  device is available; otherwise it selects optimized CPU.
- `reference` is deterministic and intended for tests and diagnosis.
- requesting an unavailable backend returns a clear `HighsStatus::kError`;
  there is no silent fallback for an explicit request;
- backend availability is reported once in the startup log;
- the runtime factory contains selection only, not algorithm logic.

Use independent CMake options instead of `CUPDLP_GPU`, for example:

```text
HIGHS_HIPDLP_REFERENCE
HIGHS_HIPDLP_OPENMP
HIGHS_HIPDLP_CUDA
HIGHS_HIPDLP_ROCM
```

Exact names should follow the prevailing HiGHS option convention when
implemented. In-tree cuPDLP configuration must remain independent.

### 4.6 Algorithm policies

Keep one backend-independent solve loop. Separate these replaceable policies:

- `RestartPolicy`;
- `StepSizePolicy`;
- `ConvergencePolicy`;
- `AveragingPolicy`.

The first policy set reproduces the current Halpern restart, fixed/PID behavior
that is reachable through HiGHS options, scaling sequence, iterate selection,
and convergence cadence. Unsupported legacy option combinations must continue
to produce a clear warning or error consistent with existing behavior.

Policies own their parameters and scalar history. The solver owns iterate
handles and orchestrates calls. A policy may request backend reductions but
must not access device pointers.

## 5. Backend-Specific Requirements

### 5.1 Reference CPU

- Use straightforward single-threaded FP64 loops.
- Favor auditability, deterministic reduction order, and bounds checks in
  debug builds over speed.
- Provide the oracle for backend primitive tests and short iteration traces.
- Avoid external linear algebra dependencies.

### 5.2 Optimized CPU

- Use owned CSR/transpose sparse kernels with static OpenMP scheduling.
- Partition work to avoid false sharing and atomics.
- Use thread-local or row-exclusive reductions as appropriate.
- Use aligned contiguous buffers and compiler-vectorizable vector loops.
- Make thread count explicit and record it in benchmark metadata.
- Keep a serial fast path for small vectors/matrices where OpenMP startup cost
  is larger than the work.
- Add optional architecture-specific tuning only behind measured dispatch;
  retain a portable implementation.
- Do not add a dense BLAS call unless benchmarks demonstrate a benefit.

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

### 5.4 ROCm

- Mirror the CUDA contract using HIP, rocBLAS, and rocSPARSE without exposing
  those types outside `backend/rocm`.
- Keep backend semantics and result normalization identical to CUDA.
- Add a compile-only CI configuration first.
- Mark hardware correctness and performance gates as pending until an AMD
  runner is available.
- Do not claim ROCm completion or performance parity from compile-only results.

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
   Compare every backend operation against reference CPU on generated matrices
   covering empty, diagonal, rectangular, skewed, duplicate-entry-normalized,
   and large-index cases.

2. **Transformation tests**
   Round-trip preprocessing, scaling, unscaling, row splitting, objective sense,
   and postsolve maps without invoking PDHG.

3. **Policy tests**
   Feed synthetic metric sequences into restart, step-size, averaging, and
   stopping policies.

4. **Iteration differential tests**
   Run a fixed number of iterations from identical initial state and compare
   iterates, products, step sizes, restart decisions, and metrics between
   reference CPU and each optimized backend. Use tolerant comparisons for
   parallel reductions; require identical decisions away from thresholds.

5. **HiGHS integration tests**
   Cover `solver = "hipdlp"`, backend selection, relevant PDLP options, model
   statuses, iteration/time limits, solution validity flags, logging disabled,
   and unavailable-backend errors.

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

- legacy HiPDLP reference CPU;
- rewritten HiPDLP reference CPU;
- legacy HiPDLP optimized CPU, if a comparable path exists;
- rewritten HiPDLP optimized CPU;
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
- Record CPU, GPU, driver, CUDA/ROCm, compiler, executable hash, git revision,
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
| Smoke | Every relevant change | Tiny representative LPs on reference and available optimized backends |
| PR | Before merge | Curated Netlib and energy set spanning tiny/large, tall/wide, and sparse/dense |
| Nightly | Scheduled | Full Netlib feasible and energy continuous manifests |
| Release | Before cutover/release | Nightly suites plus curated Mittelmann and custom problems |

Before using Mittelmann/custom in a gate, create versioned manifests with
checksums, dimensions, roles, and trusted reference objectives/statuses.

## 8. Phased Implementation

### Phase 0: Freeze baselines and acceptance artifacts

- Freeze the current git revision and build configurations.
- Build legacy reference CPU, optimized CPU if present, and CUDA variants.
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

### Phase 3: Implement reference CPU and algorithm parity

- Implement reference sparse/vector primitives and RAII host buffers.
- Port preprocessing, scaling, Halpern-PDHG, restart, step-size, convergence,
  unscale, and postsolve behavior into the new modules.
- Remove algorithm decisions from backend calls.
- Differential-test primitives, iteration traces, restart decisions, final
  metrics, status, and solution against legacy CPU.
- Fix defects found during porting only when the verifier demonstrates the
  correction; record intentional behavior changes.

**Exit criteria:** reference CPU passes unit/integration tests and correctness
gates on the PR suite, with explained parity differences.

### Phase 4: Implement optimized CPU

- Add OpenMP sparse and vector kernels plus small-problem serial dispatch.
- Add thread-count and affinity controls consistent with HiGHS.
- Tune layouts and fused kernels using profiles, not assumptions.
- Run race detection where supported and deterministic repeated-run checks.
- Benchmark against the frozen legacy CPU path.

**Exit criteria:** optimized CPU passes cross-backend correctness and the 5%
geomean no-regression gate on full Netlib and energy suites.

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

### Phase 6: Implement ROCm compile-validation path

- Implement the same backend contract with HIP/rocBLAS/rocSPARSE.
- Add CMake discovery and an isolated ROCm build configuration.
- Add compile-only CI and host-side factory/error tests.
- Produce a hardware-validation runbook and expected artifact schema.

**Exit criteria:** ROCm backend and tests compile cleanly with warnings treated
consistently with other backends. The project explicitly reports hardware
validation as pending.

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
- Confirm CMake CPU/CUDA/ROCm configurations and Bazel/Meson CPU builds.
- Update user and developer documentation.
- Remove legacy solver, duplicated CUDA state, obsolete macros, and temporary
  differential hooks after the stabilization window.

**Exit criteria:** all Definition of Done items below pass from a clean build,
and the legacy implementation is no longer needed to diagnose routine failures.

## 9. Build and CI Matrix

Minimum matrix:

| Build | Compiler/toolchain | Backends | Validation |
| --- | --- | --- | --- |
| CMake CPU minimal | GCC and Clang | reference | Unit and integration |
| CMake CPU optimized | GCC and Clang | reference + OpenMP | Unit, integration, PR benchmark |
| CMake CUDA | NVCC + host compiler | reference + CPU + CUDA | Unit, differential, GPU integration |
| CMake ROCm | HIP/Clang | reference + CPU + ROCm | Compile-only initially |
| Bazel CPU | Supported compiler | reference + CPU | Existing build/tests remain green |
| Meson CPU | Supported compiler | reference + CPU | Existing build/tests remain green |

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
| Tiny LPs regress from runtime/OpenMP overhead | Add serial small-problem dispatch and include tiny LPs in the gate |
| CUDA and ROCm cannot be enabled in one practical build | Permit any compiled subset; test separate CUDA and ROCm build configurations |
| ROCm appears supported without hardware evidence | Label it compile-validated only until AMD runs exist |
| External solver versions drift | Pin revisions/releases and hash executables/configuration |
| Benchmark selection bias | Use versioned manifests, randomized order, paired ratios, and immutable raw artifacts |
| Future low precision compromises correctness | Keep FP64 verifier/reductions and make precision an explicit capability |
| Multistart adds work without benefit | Gate on effective batch size and end-to-end time |

## 13. Definition of Done

The rewrite is complete when:

- `solver = "hipdlp"` and relevant existing options retain documented behavior;
- reference CPU, optimized CPU, and CUDA are selectable at runtime when built;
- explicit unavailable backend requests fail clearly;
- ROCm compiles in its supported toolchain and is labeled compile-validated;
- CPU-only builds have no GPU dependency;
- algorithm code contains no CUDA/ROCm types or conditional branches;
- no steady-state iteration allocates memory or performs unintended transfers;
- all backend primitives pass reference differential tests;
- supported feasible, infeasible, unbounded, limit, and numerical-error cases
  map to correct HiGHS statuses;
- every claimed optimal result passes the independent `1e-4` KKT gates;
- full Netlib and energy comparisons have no unexplained correctness
  disagreement against the baseline set;
- optimized CPU and CUDA each satisfy the `<= 1.05` paired geomean legacy
  performance ratio;
- every greater-than-20% instance regression is resolved or artifact-backed;
- CMake backend builds pass and Bazel/Meson CPU builds remain working;
- benchmark results include reproducible commands, hashes, environment,
  metadata, raw logs, parsed metrics, and verifier output;
- legacy code is removed after the stabilization window;
- low-precision SpMV and multistart can be added through the defined precision
  and batch capabilities without changing the core solve loop.

