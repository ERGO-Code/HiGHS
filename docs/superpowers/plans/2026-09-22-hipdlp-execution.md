# HiPDLP Execution Optimization Implementation Plan

> Execute the approved CPU cleanup, CUDA fusion, and production command-buffer integration as separately verified stages. Keep existing user changes and frozen experiments intact. Do not commit or create a branch.

**Goal:** Reduce redundant CPU memory work and CUDA submission work without changing the PDHG trajectory, restart/PID policy, or verification tolerances.

**Architecture:** The algorithm owns the shared step and all scalar mathematics. Backends accept optional vector-backed scalar parameters for projection/reflection and blending. A production command buffer records steps 2 through 40 once; its parameter vector is refreshed before replay. The first step remains direct for restart-reference FPE, and short tails remain direct.

**Tech Stack:** C++11, CUDA, cuSPARSE, existing opaque backend handles and command-buffer recorder, Catch2.

**Spec:** User-approved three-stage scope in the 2026-09-22 conversation.

## Constraints

- Preserve scaling, certified norm guard, FP64 arithmetic, restart, PID, certificates, and final original-model KKT acceptance.
- Preserve the direct path as an internal diagnostic control. CPU has no command-buffer dependency.
- Record at most one reusable 39-step buffer per solve; prepare buffers before iterations.
- Respect iteration limits exactly. Check wall-time before submission and at block boundaries; a submitted block cannot be interrupted midway.
- Keep device types and CUDA calls out of algorithm code. Do not bake changing step sizes or blend weights into a reusable buffer.
- Validate resource ownership, aliasing, dimensions, scalar indices, asynchronous errors, and lifetimes.
- Use new build/evidence directories under `20260922-execution`; do not overwrite historical artifacts.

## Stage 1: CPU SpMV

- [x] Add dirty-output, empty-row, all-empty, and zero-row regression coverage in `check/TestHiPdlpBackend.cpp`.
- [x] Remove only the redundant `std::fill` in `backend/cpu/sparse_kernels.cc`.
- [x] Run focused CPU tests: 54 cases, 1814 assertions.
- [x] Measure a small same-core before/after comparison: 216 samples, matching outputs, 1.6-4.7% less SpMV time; no whole-solve claim.

## Stage 2: CUDA Projection/Reflection

- [x] Require exact composed-operation results and one enqueue/no D2D copy in existing backend tests; observe 24 RED assertions.
- [x] Override both fused backend methods and instantiate dedicated fused kernels with explicit intermediate-doubling checks.
- [x] Preload both new kernel instantiations and retain standalone projection operations.
- [x] Verify focused GREEN (6 cases, 591 assertions) and affected regression union (110 cases, 19667 assertions).
- [x] Freeze the fusion-only library and matching source separately from the unfused library.
- [x] Measure fusion separately: 162/162 accepted production records across three stages, three cases, and three balanced blocks.

## Stage 3: Production Command Buffers

- [x] Add default-enabled internal `SolveContext::use_command_buffers` and diagnostic replay/iteration counters, with a direct comparison path.
- [x] Observe RED for a fixed/PID solver comparison over limits 0, 1, 39, 40, 41, 79, 80, 81, 121, and 400: 70 failed replay assertions, 310 passing numerical assertions.
- [x] Add optional `ScalarParameter { value, vector, index }` references to the three existing coarse operations. Direct callers retain literal scalar values.
- [x] Extend the shared step with an optional scalar vector and three consecutive indices, preserving one shared algorithm.
- [x] Prepare 117 scalar slots and record 39 steps using the pointer orientation after the first direct step. Refresh tau, sigma, and `1/(age+offset+2)` before replay.
- [x] Replay only when the next 39 iterations fit before the next convergence/limit boundary. Update host pointer parity and age once after successful submission.
- [x] Verify exact direct/replay solutions, iteration counts, restart/PID counts, graph reuse, zero/positive time limits, dynamic scalars, and injected failures.
- [x] Fix the reviewed capture-allocation exception path with reproduced RED and subsequent-solve GREEN; follow-up review has no important findings.
- [x] Run final CPU 148/21587 and CUDA 220/56407 cases/assertions; memcheck 10/1255 with zero errors/leaks; all GPU production records pass with setup charged.
- [x] Record measured outcomes in `docs/hipdlp-phase8-step-fusion.md`, including the excluded preexisting structure test and unverified full release gate.