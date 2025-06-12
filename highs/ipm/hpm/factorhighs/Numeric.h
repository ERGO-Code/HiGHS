#ifndef FACTORHIGHS_NUMERIC_H
#define FACTORHIGHS_NUMERIC_H

#include <memory>
#include <vector>

#include "SolveHandler.h"
#include "Symbolic.h"
#include "ipm/hpm/auxiliary/IntConfig.h"

namespace hipo {

class Numeric {
  // columns of factorisation, stored by supernode
  std::vector<std::vector<double>> sn_columns_{};

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<std::vector<Int>> swaps_{};

  // information about 2x2 pivots
  std::vector<std::vector<double>> pivot_2x2_{};

  // symbolic object
  const Symbolic& S_;

  // object to handle solve phase in different formats
  std::unique_ptr<SolveHandler> SH_;

  // lower triangle of original matrix, permuted
  std::vector<Int> rowsA_{};
  std::vector<Int> ptrA_{};
  std::vector<double> valA_{};

  // norms of columns of matrix
  std::vector<double> inf_norm_cols_, one_norm_cols_;

  friend class Factorise;

 public:
  // dynamic regularisation applied to the matrix
  std::vector<double> total_reg_{};

  Numeric(const Symbolic& S);

  // Full solve with refinement
  Int solve(std::vector<double>& x) const;

  // Iterative refinement
  Int refine(const std::vector<double>& rhs, std::vector<double>& x) const;
  std::vector<double> residual(const std::vector<double>& rhs,
                               const std::vector<double>& x) const;
  std::vector<double> residualQuad(const std::vector<double>& rhs,
                                   const std::vector<double>& x) const;
  double computeOmega(const std::vector<double>& b,
                      const std::vector<double>& x,
                      const std::vector<double>& res) const;

  void conditionNumber() const;
};

}  // namespace hipo

#endif
