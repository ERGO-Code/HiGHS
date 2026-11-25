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

class ScratchArena;

class Numeric {
  // columns of factorisation, stored by supernode
  const std::vector<std::vector<double>>* sn_columns_ = nullptr;

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<std::vector<Int>> swaps_{};

  // indicates whether a block of columns in a supernode has any non-identical
  // swap
  std::vector<std::vector<char>> any_swaps_{};

  std::vector<std::vector<double>> pivot_2x2_{};
  const Symbolic* S_;
  DataCollector* data_ = nullptr;
  const FHoptions* options_;
  std::unique_ptr<SolveHandler> SH_;

  ScratchArena* scratch_ = nullptr;

  friend class Factorise;

  // dynamic regularisation applied to the matrix
  std::vector<double> total_reg_{};

  // 1x1 pivots before regularisation, infinite for 2x2 pivots
  std::vector<double> unreg_pivots_{};

  void computeAnySwaps();

 public:
  Int prepare();

  Int solve(double* x) const;
  Int solve(double* x, Int k) const;

  Int forwardSolve(double* x) const;
  Int diagSolve(double* x) const;
  Int backwardSolve(double* x) const;

  Int forwardSolve(double* x, Int k) const;
  Int diagSolve(double* x, Int k) const;
  Int backwardSolve(double* x, Int k) const;

  void getReg(double* reg);
  void getUnregPivots(std::vector<double>& unreg_pivots) const;
  Int num2x2() const;
  void inertia(Int& pos, Int& neg, Int& zero, double tol) const;
};

}  // namespace hipo

#endif
