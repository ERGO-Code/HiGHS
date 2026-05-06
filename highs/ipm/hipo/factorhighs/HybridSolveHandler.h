#ifndef FACTORHIGHS_HYBRID_SOLVE_HANDLER_H
#define FACTORHIGHS_HYBRID_SOLVE_HANDLER_H

#include "SolveHandler.h"

namespace hipo {

class HybridSolveHandler : public SolveHandler {
  const std::vector<std::vector<Int>>& swaps_;
  const std::vector<std::vector<double>>& pivot_2x2_;

 public:
  void forwardSolve(std::vector<double>& x) const override;
  void backwardSolve(std::vector<double>& x) const override;
  void diagSolve(std::vector<double>& x) const override;

  HybridSolveHandler(const Symbolic& S,
                     const std::vector<std::vector<double>>& sn_columns,
                     const std::vector<std::vector<Int>>& swaps,
                     const std::vector<std::vector<double>>& pivot_2x2,
                     DataCollector& data);
};

}  // namespace hipo

#endif