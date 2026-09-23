# HiPDLP Rewrite Architecture Visualization Design

Date: 2026-09-23

## 1. Purpose

Create a self-contained Chinese HTML architecture explorer at
`docs/hipdlp-architecture.html` for the current rewritten HiPDLP
implementation. The page must help a future reader:

1. understand how one solve moves through the C++ and CUDA code;
2. distinguish CPU control work, CPU numerical work, GPU numerical work,
   host-to-device transfers, and device-to-host transfers;
3. understand the major solver-owned data structures and their lifetimes;
4. find likely performance bottlenecks without confusing measured evidence
   with hypotheses;
5. learn how to read and review this repository's C++ and CUDA code.

The page is a code-oriented companion to, not a replacement for, the existing
`docs/hipdlp-rewrite-progress.html` performance-progress dashboard.

## 2. Snapshot And Evidence Boundary

The initial page describes:

- branch `feature/hipdlp-rewrite`;
- source commit `49b85f4f3`;
- generation date 2026-09-23;
- the current source architecture under `highs/pdlp/hipdlp`;
- the current historical and diagnostic reports under `docs`.

The repository working tree is dirty. The page must state that its architecture
claims are tied to the source commit while report evidence may also come from
currently uncommitted documentation. It must not imply that every repository
artifact is reproducible from the commit alone.

Every claim in the page must be classified as one of:

- **implemented fact**: directly visible in current source;
- **measured evidence**: tied to a dated report and its stated scope;
- **optimization hypothesis**: plausible but not yet measured;
- **future placeholder**: registered but not implemented.

## 3. Technical Form

Use one dependency-free HTML file with:

- semantic HTML;
- embedded CSS;
- embedded JavaScript;
- inline SVG for diagrams;
- no CDN, package, build step, web server, or network request.

Core explanatory content must remain visible when JavaScript is disabled.
JavaScript enhances navigation, filtering, drawers, calculators, exercises,
and highlighting. The document must support keyboard navigation, reduced
motion, narrow screens, and printing.

## 4. Visual Language

Use a stable legend throughout the page:

| Meaning | Color role |
| --- | --- |
| CPU / host execution or storage | Blue |
| GPU / device execution or storage | Green |
| Host-to-device transfer | Orange |
| Device-to-host transfer | Red |
| Backend-neutral control flow | Gray |
| Correctness gate or independent verification | Purple |
| Measured bottleneck | Amber |
| Unavailable/future backend | Muted dashed treatment |

Animations may show movement through the data-flow diagram, but all meaning
must also be encoded by labels, line styles, and icons. Animation must be
disabled under `prefers-reduced-motion`.

## 5. Information Architecture

The page has a persistent snapshot header and five main views.

### 5.1 System Map

Show the end-to-end dependency direction:

```text
Highs / solveLpHiPdlp
  -> API option translation and migration routing
  -> problem construction and scaling
  -> backend-independent algorithm control
  -> CPU or CUDA backend
  -> postsolve and original-model verification
  -> HighsSolution / HighsInfo / model status
```

The map must expose the current legacy/rewrite branch without allowing the
legacy implementation to obscure the rewrite flow.

Primary source anchors:

- [`HiPdlpWrapper.cpp`](../../../highs/pdlp/HiPdlpWrapper.cpp)
- [`api/solver.cc`](../../../highs/pdlp/hipdlp/api/solver.cc)
- [`algorithm/solver.cc`](../../../highs/pdlp/hipdlp/algorithm/solver.cc)
- [`backend/backend_factory.cc`](../../../highs/pdlp/hipdlp/backend/backend_factory.cc)
- [`verification/kkt_verifier.cc`](../../../highs/pdlp/hipdlp/verification/kkt_verifier.cc)

### 5.2 Data And Residency

Show the transformation:

```text
HighsLp
  -> model::OriginalModel (non-owning original-LP view)
  -> model::Problem::Data (owned transformed/scaled host model)
  -> backend handles and resources
  -> HostIterate
  -> postsolved HighsSolution
```

The `Problem::Data` card must include:

- dimensions, objective sense, and offset;
- costs and column bounds;
- normalized row bounds and equality metadata;
- CSR `A` and separately materialized CSR `A^T`;
- scaling metadata;
- transformation and slack-variable maps;
- unscaled norms;
- index-width metadata.

The backend-neutral handle card must explain that `ProblemHandle`,
`VectorHandle`, `WorkspaceHandle`, and `CommandBufferHandle` are move-only
owners of an opaque `Resource`.

The CPU resource cards must show:

- copied host model arrays and matrices;
- `std::vector<double>` numerical vectors;
- row and column scratch vectors.

The CUDA resource cards must show:

- device copies of costs, bounds, row metadata, `A`, and `A^T`;
- cuSPARSE sparse descriptors and SpMV workspace;
- resident device vectors;
- row/column scratch;
- resident check/certificate storage;
- pinned result and vector download buffers;
- stream, cuBLAS/cuSPARSE handles, error flags, timing events, and graph
  resources.

The page must distinguish ownership, non-owning pointers, shared lifetime
tokens, and graph dependencies.

### 5.3 Iteration And Correctness

Provide three selectable lanes:

- CPU;
- CUDA;
- side-by-side comparison.

Show preprocessing and scaling on the host for both paths. Show CUDA problem
binding and vector initialization as explicit H2D/setup work.

The normal Halpern-PDHG step must follow
[`enqueueHalpernPdhgStep`](../../../highs/pdlp/hipdlp/algorithm/halpern_pdhg_step.cc):

1. projected primal update plus reflection;
2. `A * reflected_primal` plus projected dual update and reflection;
3. Halpern blend of anchor and reflected primal/dual;
4. `A^T * next_dual`;
5. swap current/next handles.

The page must state that the production step performs two SpMVs. Candidate
`A*x` and `A^T*y` products are separately prepared at boundaries when needed.

For CUDA, show that steps 2 through 40 are captured once and replayed as a
39-step command buffer when eligible. A 117-double parameter vector carries
per-step primal/dual steps and blend weights; one replay uploads 936 bytes.
Step 1 after a restart and short tails remain direct. A submitted 39-step
block is not interruptible by the host time-limit check.

At the convergence cadence, show:

1. candidate product preparation;
2. optional fixed-point preparation;
3. backend-resident raw convergence/FPE/anchor/certificate screening;
4. the normal 112-byte CUDA `CheckPacket` D2H path;
5. full candidate download and host oracle when arithmetic uncertainty,
   optional metric fallback, a suspected certificate, or a final result
   requires it;
6. restart and PID decisions in the algorithm layer;
7. final postsolve and independent original-model KKT verification on CPU.

The final KKT gate must be visually separate from backend convergence. The page
must never suggest that backend convergence alone is sufficient for `OPTIMAL`.

### 5.4 Bottlenecks And Optimization

Each hotspot card has five sections:

1. symptom;
2. current evidence;
3. relevant counters or profiler signals;
4. a bounded next experiment;
5. correctness invariants that the experiment must preserve.

Cover at least:

- host preprocessing and small-model fixed cost;
- serial CPU SpMV and fused vector loops;
- CUDA setup, uploads, descriptors, and graph construction;
- steady-state CUDA SpMV and fused kernels;
- launch latency and graph replay;
- cadence-boundary resident checks;
- host-oracle fallback, full candidate D2H, and synchronization;
- convergence trajectory changes from scaling, norm estimation, PID, and
  restart;
- postsolve and original-model KKT.

Measured evidence must retain its original scope:

- the 2026-09-08 release gate remains blocked and its incomplete denominator
  must not be converted into a formal geometric mean;
- the 2026-09-14 two-SpMV/CPU fusion comparison is an engineering comparison,
  not release acceptance;
- the 2026-09-22 three-case CUDA fusion-plus-replay result reports a 41.5% to
  52.4% reduction versus its scoped before version and must not be generalized
  to all instances;
- the 2026-09-22 large-instance screen has no fully accepted legacy/rewrite
  performance pairs and therefore cannot establish a speedup.

Source reports:

- [`hipdlp-phase8-cpu-diagnosis.md`](../../hipdlp-phase8-cpu-diagnosis.md)
- [`hipdlp-phase8-gpu-performance-report.md`](../../hipdlp-phase8-gpu-performance-report.md)
- [`hipdlp-phase8-step-fusion.md`](../../hipdlp-phase8-step-fusion.md)
- [`hipdlp-phase8-large20-screen.md`](../../hipdlp-phase8-large20-screen.md)

#### Memory And Transfer Estimator

Accept `rows`, `columns`, `nonzeros`, and index width. Report labeled lower
bounds for:

- core host/device `A + A^T` compressed storage;
- the eight column-sized and eight row-sized iteration vectors created by the
  control loop;
- one full candidate download
  (`primal + dual + Ax + A^T y = 2 * (rows + columns)` doubles);
- one normal CUDA check packet (112 bytes).

Clearly exclude allocator overhead, descriptors, cuSPARSE workspaces, resident
check/certificate duplicates, reduction scratch, pinned staging, graph
objects, and vendor-internal allocations.

#### Diagnostics Map

Map every relevant `SolveDiagnostics` field back to a phase or transfer on the
diagram. Explain which question is best answered by:

- existing solver diagnostics and trace;
- a CPU profiler;
- Nsight Systems;
- Nsight Compute;
- Compute Sanitizer.

### 5.5 Learning And Review

Provide guided routes for a quick architecture pass, a C++ ownership and
interface pass, a systematic CUDA course, and focused correctness/performance
reviews. This view owns the source-reading sequence, glossary, exercises,
collapsed answers, and reusable review checklists specified in Sections 7 and
8.

## 6. CUDA Kernel Drill-Down

The iteration view must expand through four levels:

```text
algorithm operation
  -> Backend virtual method
  -> CUDA launch wrapper or cuSPARSE call
  -> __global__ kernel and buffers
```

Show these mappings:

| Algorithm/backend operation | CUDA implementation |
| --- | --- |
| `multiply` / `multiplyTranspose` | `executeSpmv` -> `cusparseSpMV` |
| primal project/reflect | `launchProjectedPrimalUpdateAndReflect` -> `projectedPrimalKernel<true>` |
| dual project/reflect | `launchProjectedDualUpdateAndReflect` -> `projectedDualKernel<true>` |
| paired Halpern blend | `launchVectorPairBlend` -> `vectorPairBlendKernel` |
| resident partial check | `launchCheckPartials` -> `checkPartialsKernel` |
| check packet reduction | `launchCheckPacket` -> `checkPacketKernel` |
| certificate screen | recovery, norm, normalize, conditions, and flags kernels |
| replay | stream capture -> graph instantiate -> graph launch |

Execution configuration:

- general vector kernels use 256 threads per block;
- their grid is `min(ceil(N / 256), 65535)`;
- pair blend uses `N = max(num_cols, num_rows)`;
- check/certificate first stages use
  `min(ceil(max(rows, cols) / 256), 1024)` blocks;
- packet/norm/flag final reductions use one block of 256 threads;
- launch wrappers request zero bytes of dynamic shared memory;
- reduction kernels use compile-time shared arrays sized for 256 entries;
- cuSPARSE chooses its internal execution configuration.

For every custom kernel, show:

- logical element count;
- grid/block formula;
- grid-stride behavior;
- device buffers read and written;
- branches and atomics;
- error-flag behavior;
- whether it participates in graph capture;
- likely performance dimension;
- review questions.

## 7. C++ And CUDA Learning Design

### 7.1 Learning Routes

Provide:

- a 30-minute end-to-end route;
- a 90-minute C++ architecture route;
- a systematic CUDA route;
- focused review routes for correctness, performance, and ownership.

The C++ route covers:

- wrapper and API boundaries;
- option translation and status mapping;
- immutable original-model boundary;
- owned transformed problem;
- move-only opaque handles;
- virtual backend contract;
- RAII phase timing and CUDA resources;
- structured error propagation;
- postsolve and verification.

### 7.2 Systematic CUDA Course

Organize the CUDA material into eight lessons:

1. grid, block, warp, thread, and grid-stride execution;
2. host, pinned host, global device, register, and shared memory;
3. asynchronous launch, stream ordering, launch errors, asynchronous errors,
   synchronization, and event timing;
4. coalescing, divergence, occupancy, launch latency, bandwidth, atomics, and
   reductions;
5. cuSPARSE descriptors, buffers, and the boundary between vendor and custom
   kernels;
6. stream capture, graph instantiation, dependencies, parameter data, replay,
   and block-level time-limit granularity;
7. FP64 arithmetic, explicit rounding intrinsics, error flags, deterministic
   reductions, and interval screening;
8. selecting Compute Sanitizer, Nsight Systems, or Nsight Compute for a
   concrete question.

Each lesson includes:

- a concept diagram;
- a real HiPDLP code example;
- three self-check questions;
- one review checklist;
- a collapsed answer for every question.

Each answer has:

1. a one-sentence conclusion;
2. direct code evidence;
3. review reasoning, including what static reading cannot prove and what
   profiler evidence would be needed.

Include integrated exercises that trace `primal_candidate`, judge whether a
D2H transfer is necessary, prove grid-stride coverage for a large vector, and
check graph dependency lifetime.

## 8. Review Mode

Review mode applies four passes:

1. interface and ownership;
2. numerical correctness;
3. concurrency and CUDA semantics;
4. performance.

Every component drawer includes:

- responsibility;
- inputs and outputs;
- owner and lifetime;
- CPU/GPU residency;
- error paths;
- invariants;
- recommended source-reading order;
- relevant tests;
- good-code signals;
- danger signals;
- reusable human/AI review questions.

The page must teach reviewers to separate facts available from source from
claims requiring runtime measurement.

## 9. Error And Edge-Case Presentation

The diagrams must show explicit error exits for:

- option validation;
- unavailable or not-ready backends;
- allocation failures;
- malformed models;
- dimension and alias violations;
- non-finite arithmetic;
- CUDA runtime/cuBLAS/cuSPARSE failures;
- graph capture/instantiation/replay failures;
- invalid resident-check results;
- postsolve or KKT input failures.

Do not render errors as silent fallback. The one intended runtime fallback that
must be explained is `auto` backend selection trying CUDA and then choosing CPU
when CUDA creation is unavailable. An explicitly requested CUDA backend does
not silently fall back.

## 10. Verification And Acceptance

The implementation is complete only when:

1. the file opens directly through `file://` without external assets;
2. all five views are reachable by mouse and keyboard;
3. CPU, CUDA, and comparison filters show the correct lanes;
4. every architecture node has its expected source and detail drawer;
5. the estimator returns exact formulas for known synthetic inputs and rejects
   missing, negative, non-integer, or unsafe values;
6. exercise answers are hidden by default and keyboard accessible;
7. reduced-motion mode removes nonessential animation;
8. the narrow/mobile layout preserves labels and does not require horizontal
   page scrolling outside intentionally scrollable diagrams;
9. print mode exposes all core static content;
10. local source links resolve to existing repository files;
11. the browser console has no errors;
12. the page contains no external requests;
13. factual spot checks match the current source for the solve entry, problem
    representation, iteration vector count, two-SpMV step, CUDA launch
    configuration, 112-byte packet, and command-buffer replay;
14. measured performance claims retain dates, scope, and limitations.

Validation should use a lightweight browser smoke test plus manual responsive,
keyboard, print, and reduced-motion checks. No new production dependency is
needed.

## 11. Deliverables

- `docs/hipdlp-architecture.html`
- directly related tests only if an existing HTML/browser test pattern can be
  reused without adding a new framework
- no changes to `docs/hipdlp-rewrite-progress.html` in this task
- no solver behavior changes
