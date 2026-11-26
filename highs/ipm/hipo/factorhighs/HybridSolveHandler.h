#ifndef FACTORHIGHS_HYBRID_SOLVE_HANDLER_H
#define FACTORHIGHS_HYBRID_SOLVE_HANDLER_H

#include "ScratchArena.h"
#include "SolveHandler.h"

namespace hipo {

class HybridSolveHandler : public SolveHandler {
  // swaps, 2x2 pivots and any_swaps: one entry per column, from snStart(sn)
  const std::vector<Int>& swaps_;
  const std::vector<char>& any_swaps_;
  const std::vector<double>& pivot_2x2_;

  // gemv workspaces come from the arena; they never live across a taskWait,
  // so each thread holds at most one at a time
  ScratchArena& scratch_;

  std::vector<Int> first_child_, next_child_;
  mutable std::vector<std::vector<Int>> task_rows_;
  mutable std::vector<std::vector<double>> task_vals_;

  void processForwardTask(Int task, double* x) const;
  void processBackwardTask(Int task, double* x) const;

  void processForwardSn(Int sn, double* x, double* work, Int task = -1,
                        Int end_col_in_task = -1) const;
  void processBackwardSn(Int sn, double* x, double* work) const;
  void processDiagSn(Int sn, double* x) const;

 public:
  void forwardSolve(double* x) const override;
  void backwardSolve(double* x) const override;
  void diagSolve(double* x) const override;
  void inertia(Int& pos, Int& neg, Int& zero, double tol) const override;

  HybridSolveHandler(const Symbolic& S, const FactorStorage& factor,
                     const std::vector<Int>& swaps,
                     const std::vector<char>& any_swap,
                     const std::vector<double>& pivot_2x2, DataCollector& data,
                     const FHoptions& options, ScratchArena& scratch);
};

}  // namespace hipo

#endif