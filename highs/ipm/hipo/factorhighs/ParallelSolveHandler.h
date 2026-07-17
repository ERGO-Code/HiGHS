#ifndef FACTORHIGHS_PARALLEL_SOLVE_HANDLER_H
#define FACTORHIGHS_PARALLEL_SOLVE_HANDLER_H

#include "SolveDag.h"
#include "SolveHandler.h"

namespace hipo {

// Experimental parallel triangular solve (see PHASE3_DESIGN.md).
//
// Executes the same per-supernode arithmetic as HybridSolveHandler, but
// schedules supernodes by levels of the supernodal elimination tree:
// - forwardSolve : levels bottom-up; tasks within a level run in parallel,
//   writing gemv results into per-task buffers. Scatter contributions to
//   rows of the task's own supernode are applied immediately; contributions
//   to ancestor rows are deferred and merged serially at the end of the
//   level (in ascending supernode order), so the result is deterministic
//   for a fixed schedule.
// - backwardSolve: levels top-down; tasks within a level are independent
//   (each supernode reads ancestor entries, which are final, and writes
//   only its own block range) - bit-identical to the serial handler.
// - diagSolve    : embarrassingly parallel - bit-identical.
class ParallelSolveHandler : public SolveHandler {
  const std::vector<std::vector<Int>>& swaps_;
  const std::vector<std::vector<uint8_t>>& swap_flags_;
  const std::vector<std::vector<double>>& pivot_2x2_;

  // levels -> tasks -> ascending supernode lists (mode 1)
  const std::vector<std::vector<std::vector<Int>>>& schedule_;

  // block-level DAG (mode 2); nullptr in mode 1
  const SolveDag* dag_;

  // forward-solve stash of post-dtrsv (pre-inverse-swap) block values,
  // written by solve tasks and read by delivery tasks (mode 2)
  std::vector<double>& fwd_stash_;

  // per-task-slot buffers, owned by Numeric and reused across solves
  std::vector<std::vector<double>>& task_work_;
  std::vector<std::vector<Int>>& task_def_rows_;
  std::vector<std::vector<double>>& task_def_vals_;

  // Forward-solve one task's supernodes. If defer is true, contributions to
  // rows outside the supernode are recorded in the slot's deferred buffers
  // instead of being applied to x.
  void forwardTask(const std::vector<Int>& sns, std::vector<double>& x,
                   bool defer, Int slot) const;

  // Backward-solve one task's supernodes, in reverse task order.
  void backwardTask(const std::vector<Int>& sns, std::vector<double>& x,
                    Int slot) const;

  // Backward-solve a single supernode.
  void backwardOne(Int sn, std::vector<double>& x, Int slot) const;

  // Diagonal-solve a range of supernodes.
  void diagRange(Int sn_begin, Int sn_end, std::vector<double>& x) const;

  // ==== block-level DAG paths (mode 2) ====
  // SolveTask(run): per member supernode, swaps, dtrsv and the gemv over
  // rows owned by the run; stashes post-dtrsv block values for deliveries.
  void dagSolveTask(Int run, std::vector<double>& x, Int slot) const;
  // DeliveryTask: apply one coarsened external segment, block by block.
  void dagDeliveryTask(const SolveDag::Delivery& D, std::vector<double>& x,
                       Int slot) const;
  void forwardSolveDag(std::vector<double>& x) const;
  void backwardSolveDag(std::vector<double>& x) const;

 public:
  void forwardSolve(std::vector<double>& x) const override;
  void backwardSolve(std::vector<double>& x) const override;
  void diagSolve(std::vector<double>& x) const override;

  ParallelSolveHandler(
      const Symbolic& S, const std::vector<std::vector<double>>& sn_columns,
      const std::vector<std::vector<Int>>& swaps,
      const std::vector<std::vector<uint8_t>>& swap_flags,
      const std::vector<std::vector<double>>& pivot_2x2,
      const std::vector<std::vector<std::vector<Int>>>& schedule,
      const SolveDag* dag, std::vector<double>& fwd_stash,
      std::vector<std::vector<double>>& task_work,
      std::vector<std::vector<Int>>& task_def_rows,
      std::vector<std::vector<double>>& task_def_vals, DataCollector& data);
};

}  // namespace hipo

#endif
