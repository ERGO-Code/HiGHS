#ifndef FACTORHIGHS_SOLVE_HANDLER_H
#define FACTORHIGHS_SOLVE_HANDLER_H

#include <vector>

#include "Symbolic.h"
#include "ipm/hpm/auxiliary/IntConfig.h"

namespace hipo {

// Interface class to handle different formats of dense matrices during the
// solve phase.
// Any implementation of a specific format needs to define:
// - forwardSolve  : to perform the solve with matrix L
// - backwardSolve : to perform the solve with matrix L^T
// - diagSolve     : to perform the solve with matrix D

class SolveHandler {
 protected:
  const Symbolic& S_;
  const std::vector<std::vector<double>>& sn_columns_;

 public:
  SolveHandler(const Symbolic& S,
               const std::vector<std::vector<double>>& sn_columns);

  // avoid copies
  SolveHandler(const SolveHandler&) = delete;
  SolveHandler& operator=(const SolveHandler&) = delete;

  // virtual destructor
  virtual ~SolveHandler() = default;

  // =================================================================
  // Pure virtual functions.
  // These need to be defined by any derived class.
  // =================================================================
  virtual void forwardSolve(std::vector<double>& x) const = 0;
  virtual void backwardSolve(std::vector<double>& x) const = 0;
  virtual void diagSolve(std::vector<double>& x) const = 0;
};

}  // namespace hipo

#endif