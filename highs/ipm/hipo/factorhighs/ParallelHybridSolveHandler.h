#ifndef FACTORHIGHS_PARALLEL_HYBRID_SOLVE_HANDLER_H
#define FACTORHIGHS_PARALLEL_HYBRID_SOLVE_HANDLER_H

#include "SolveHandler.h"

namespace hipo {

class ParallelHybridSolveHandler : public SolveHandler {
  const std::vector<std::vector<Int>>& swaps_;
  const std::vector<std::vector<double>>& pivot_2x2_;

 public:
  void forwardSolve(std::vector<double>& x) const override;
  void backwardSolve(std::vector<double>& x) const override;
  void diagSolve(std::vector<double>& x) const override;

  ParallelHybridSolveHandler(const Symbolic& S,
                             const std::vector<std::vector<double>>& sn_columns,
                             const std::vector<std::vector<Int>>& swaps,
                             const std::vector<std::vector<double>>& pivot_2x2,
                             DataCollector& data);
};

}  // namespace hipo

#endif