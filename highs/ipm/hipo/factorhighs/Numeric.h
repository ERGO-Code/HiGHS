#ifndef FACTORHIGHS_NUMERIC_H
#define FACTORHIGHS_NUMERIC_H

#include <memory>
#include <vector>

#include "DataCollector.h"
#include "FactorHighsOptions.h"
#include "SolveHandler.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

// Numeric allows to perform solve, though a pointer to the numerical factor,
// that is stored in FHsolver. It also holds auxiliary data about swaps,
// pivots...

namespace hipo {

class Numeric {
  // columns of factorisation, stored by supernode
  const std::vector<std::vector<double>>* sn_columns_ = nullptr;

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<std::vector<Int>> swaps_{};

  // information about 2x2 pivots
  std::vector<std::vector<double>> pivot_2x2_{};

  // symbolic object
  const Symbolic* S_;

  DataCollector* data_ = nullptr;

  const FHoptions* options_;

  friend class Factorise;

  // dynamic regularisation applied to the matrix
  std::vector<double> total_reg_{};

 public:
  Int solve(double* x) const;
  Int solve(double* x, Int k) const;

  Int forwardSolve(double* x) const;
  Int diagSolve(double* x) const;
  Int backwardSolve(double* x) const;

  Int forwardSolve(double* x, Int k) const;
  Int diagSolve(double* x, Int k) const;
  Int backwardSolve(double* x, Int k) const;

  void getReg(double* reg);
  void inertia(Int& pos, Int& neg, Int& zero, double tol) const;
};

}  // namespace hipo

#endif
