# Parallel triangular solve for HiPO: framework sketch

Status: **experimental sketch**, off by default (`hipo_parallel_solve = true` to enable,
together with `parallel = on`). The goal of this code is to demonstrate a workable
architecture and its correctness argument before investing in performance work, and to
serve as a concrete basis for discussion with the HiPO authors. It is *not* expected to
be faster yet — see "What is deliberately not done" below.

## Why

On PyPSA-class LPs the triangular solve phase (forward/diag/backward, ~8 solves per IPM
iteration including correctors and refinement) is 55–65% of HiPO's numerical work at
every problem size, and it is serial: measured on `pypsa-de-elec-10-1h`, going from 2 to
4 threads speeds the factorisation (165s → 126s per 600s window) while solves don't move
(~263s) — an Amdahl ceiling of ~2× however many cores are available. The HiPO paper
(arXiv:2508.04370) lists solve-phase parallelisation as future work.

## Architecture

Three pieces, all keyed off the supernodal elimination tree that `Symbolic` already
stores (`snParent`):

1. **Solve schedule** (`Numeric::buildSolveSchedule`, built once per factorisation):
   supernodes are grouped into *levels* (`level[sn] = 1 + max level of children`), so
   all dependencies of a supernode lie in strictly lower levels. Within a level,
   supernodes are chunked into ascending-order *tasks* (≥8 supernodes per task, ≤64
   tasks). Consecutive singleton levels — the signature of the chain-shaped trees that
   year-long hourly models produce — are merged into one serial task to avoid paying a
   barrier per chain element.

2. **`ParallelSolveHandler`** mirrors `HybridSolveHandler`'s per-supernode arithmetic
   exactly (same BLAS calls, same order within a supernode) and schedules it:
   - `diagSolve`: `parallel::for_each` over supernodes. Each supernode reads and writes
     only its own entries of `x` → race-free and **bit-identical** to the serial handler.
   - `backwardSolve`: levels top-down, tasks within a level in parallel. A supernode
     *gathers* from ancestor rows (final — ancestors live in higher levels, already
     processed) and *writes* only its own block range → race-free and **bit-identical**.
   - `forwardSolve`: levels bottom-up. The scatter `x[row] -= y[i]` is the shared-write
     hazard: rows within the supernode's own columns are needed by its own later blocks
     and touched only by the owning task, so they are applied immediately; rows in the
     clique (ancestors) are recorded in per-task buffers and merged serially at the end
     of the level, in ascending supernode order. Single-task levels (chains) skip the
     deferral entirely and behave exactly like the serial handler.

3. **Buffers** owned by `Numeric`, sized once per factorisation and reused across the
   ~8 solves per iteration: one gemv workspace per task slot plus growable
   deferred-scatter buffers.

## Correctness argument

- Levels encode the dependency partial order of both triangular solves; tasks never
  read or write another same-level task's rows except through the deferred merge.
- The merge applies contributions in ascending supernode order within each level, and
  levels are processed in a fixed order → the schedule fully determines the arithmetic
  → **the parallel solve is deterministic for a fixed factorisation and schedule,
  independent of thread count and timing.**
- `backwardSolve` and `diagSolve` are bit-identical to the serial handler. `forwardSolve`
  is *not* guaranteed bit-identical: the serial handler applies scatters in global
  ascending supernode order, while level scheduling applies them level-by-level, and
  supernode index order and level order need not agree across subtrees; contributions to
  a shared ancestor row can therefore be summed in a different (but fixed) order.
  Expect iterate trajectories to agree to rounding, not to the last bit; validation is
  therefore (a) determinism (two runs with the flag on produce identical logs), and
  (b) optimality (same objective, clean KKT residuals) — not a log diff against serial.

## What is deliberately not done yet (the real Phase 3 work)

- **Block-level DAG** (Hogg–Reid–Scott, SIAM JSC 2010): dependencies between column
  blocks rather than whole supernodes, which is what actually exposes pipeline
  parallelism *inside* the chains that dominate hourly PyPSA trees. The level model
  here leaves chains serial (merged tasks), so the exploitable parallelism on the worst
  instances is limited — this sketch parallelises the bushy part of the tree only.
- **Task-grain tuning / dependency sparsification** (Park et al., ISC 2014): the
  8-supernode grain and 64-task cap are placeholders.
- **Parallel merge**: the deferred-scatter merge is serial; it could be row-partitioned.
- **Multi-RHS composition** (Phase 2): batching the ~8 solves per iteration through a
  BLAS3 path composes with, and is independent of, this scheduling.
- Handler code duplication with `HybridSolveHandler` would be removed by templating the
  per-supernode kernel over a scatter policy.
- The enable flag lives in `HipoTuning` for convenience; an upstreamable version
  belongs in `hipo::Options`.

## Measurement protocol for the go/no-go gate

On `pypsa-de-elec-10-1h` (public instance, see HIPO_PYPSA_NOTES.md), 600s budget,
`log_dev_level=1`: compare `Solve time:` and iterations reached for serial vs parallel
solve at 4 and 8+ threads. Gate for continuing the full Phase 3: **≥1.5× solve-phase
speedup at 4 threads** with determinism and clean KKT residuals. The sketch as-is is
expected to fail that gate on chain-dominated instances (by design); the gate applies
after the block-level DAG lands.
