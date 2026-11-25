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
  const SnColumns* sn_columns_ = nullptr;

  // Swaps and 2x2 pivots have one entry per column: supernode sn owns the
  // entries from S_->snStart(sn).

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<Int> swaps_{};

  // any_swaps_[S_->snStart(sn) + j] indicates whether block j of columns of
  // supernode sn has any non-identical swap
  std::vector<char> any_swaps_{};

  std::vector<double> pivot_2x2_{};
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
