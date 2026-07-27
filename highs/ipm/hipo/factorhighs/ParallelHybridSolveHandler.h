#ifndef FACTORHIGHS_PARALLEL_HYBRID_SOLVE_HANDLER_H
#define FACTORHIGHS_PARALLEL_HYBRID_SOLVE_HANDLER_H

#include "SolveHandler.h"

namespace hipo {

class ParallelHybridSolveHandler : public SolveHandler {
  const std::vector<std::vector<Int>>& swaps_;
  const std::vector<std::vector<double>>& pivot_2x2_;

  std::vector<Int> first_child_, next_child_;

  mutable std::vector<std::vector<Int>> task_rows_;
  mutable std::vector<std::vector<double>> task_vals_;

  void processForwardTask(Int task, double* x) const;
  void processBackwardTask(Int task, double* x) const;

 public:
  void forwardSolve(double* x) const override;
  void backwardSolve(double* x) const override;
  void diagSolve(double* x) const override;

  ParallelHybridSolveHandler(const Symbolic& S,
                             const std::vector<std::vector<double>>& sn_columns,
                             const std::vector<std::vector<Int>>& swaps,
                             const std::vector<std::vector<double>>& pivot_2x2,
                             DataCollector& data, const FHoptions& options);

  void inertia(Int& pos, Int& neg, Int& zero, double tol) const override{};
};

}  // namespace hipo

#endif