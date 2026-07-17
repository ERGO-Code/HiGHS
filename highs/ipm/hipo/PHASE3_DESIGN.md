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

# Block-level DAG (parallel_solve_mode = 2)

Mode 2 implements the block-level DAG in `factorhighs/SolveDag.{h,cpp}` and the
`dag*` methods of `ParallelSolveHandler`.

## Structure

- **Run coarsening**: supernode sn joins the current run iff *all* of its incoming
  update sources are already members (a pure chain step) and the run is below
  `hipo_dag_min_task_ops` estimated flops. Runs are contiguous in the postordered
  index, so each run owns a contiguous global row range. Rows of a supernode's panel
  below `intEnd(sn)` belong to its own run and are applied inline by the run's solve
  task; rows above are covered by **delivery tasks** (segments of at least
  `hipo_dag_min_seg_rows` rows, split only at destination-run boundaries).
- **Forward dependencies**: SolveTask(r) waits for every delivery into r; a delivery
  waits for its source run's SolveTask and for the previous delivery in each covered
  destination's inbox (inboxes ordered by ascending source supernode). All edges
  increase the (source, kind) key, so the graph is acyclic by construction.
- **Backward**: one task per run, members processed in descending order;
  BwdTask(r) waits for the runs owning its members' external panel rows (read-only
  dependencies — values read are final, so execution order cannot affect results).
- **Executor**: LIFO ready-stack with a chain fast-path (a task that unlocks exactly
  one dependent hands it straight to the same worker with no locking — depth-first
  execution preserves the serial traversal's memory locality) and idle backoff.

## Determinism, not bit-identity — and why (measured)

The original design aimed for bit-identity with the serial handler via row-split
gemvs. **Measurement killed that assumption**: OpenBLAS `dgemv('T')` is *not* bitwise
row-splittable — per-output dot products change with the output count n (verified
with a direct test: jb in 2..64, splits at many points, diffs of ~1 ulp appear
routinely). Any genuine block decomposition therefore cannot reproduce serial results
bit-for-bit with BLAS kernels.

What mode 2 guarantees instead, both verified on full IPM runs:
- **Bit-determinism across thread counts and repeats**: the application order into
  every row of x is fixed by the inbox chains and split points, independent of
  scheduling, so 1-thread and N-thread runs produce identical iterate sequences.
- **Optimality equivalence**: same optimal objective to all printed digits and clean
  KKT residuals; iterate trajectories may differ from serial after several
  iterations (a different, equally valid FP summation grouping).

## Measured results (4-vCPU cloud VM, post v1.15.1 merge)

`pypsa-de-elec-10-1h`, 600s budget: serial (default) 45 iterations / 159.7s solve;
mode 2 at 4 threads 26 iterations / 340.0s solve (per-solve 0.49s vs 1.62s).
Small instance full solve: serial 5.5s, mode 2 at 4 threads 20.5s.
**Mode 2 loses decisively on this machine.** Single-thread mode 2 is close to
serial (locality preserved by the LIFO/chain-fast-path executor), so the loss is
specific to multi-worker execution: triangular solves are memory-bandwidth-bound
and the shared-bandwidth vCPUs cannot feed multiple workers — parallel execution
mostly migrates cache lines. The go/no-go question (≥1.5× at 4 threads) remains
open for hardware with real per-core bandwidth; if it fails there too, the
conclusion is that solve parallelism requires layout-level changes (blocked
multi-RHS, GPU offload) rather than scheduling alone.

## Remaining work beyond this implementation

- The backward solve is run-grain (block-splitting its gathers would split dot-product
  summations; with bit-identity to serial already off the table, a segmented
  deterministic accumulation is now a reasonable next step).
- Delivery tasks recompute their slice of every column block's sub-gemv; per-block
  flops match serial but small blocks pay per-call BLAS overhead.
- Worker count is not adapted to DAG width (idle workers on narrow DAGs cost a
  little); `hipo_dag_min_task_ops` / `hipo_dag_min_seg_rows` defaults are first
  guesses, exposed for the tuning loop.
