#ifndef FACTORHIGHS_FACTORISE_H
#define FACTORHIGHS_FACTORISE_H

#include <cmath>

#include "Numeric.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/Log.h"

namespace hipo {

class Factorise {
 public:
  // matrix to factorise
  std::vector<Int> rowsA_{};
  std::vector<Int> ptrA_{};
  std::vector<double> valA_{};
  Int n_{};
  Int nzA_{};

  // symbolic factorisation
  const Symbolic& S_;

  // children in supernodal elimination tree
  std::vector<Int> first_child_{};
  std::vector<Int> next_child_{};

  // reverse linked lists of chidlren
  std::vector<Int> first_child_reverse_{};
  std::vector<Int> next_child_reverse_{};

  // generated elements, aka Schur complements.
  std::vector<std::vector<double>> schur_contribution_{};

  // columns of L, stored as dense supernodes
  std::vector<std::vector<double>> sn_columns_{};

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<std::vector<Int>> swaps_{};

  // Information about 2x2 pivots.
  // If pivot_2x2[sn][i] == 0, 1x1 pivot was used.
  // If pivot_2x2[sn][i] != 0, 2x2 pivot was used and pivot_2x2[sn][i] stores
  //  the off-diagonal pivot entry (of the 2x2 inverse).
  std::vector<std::vector<double>> pivot_2x2_{};

  // largest diagonal element in the original matrix and norms of columns
  double max_diag_{};
  double min_diag_{};
  double A_norm1_{};
  std::vector<double> one_norm_cols_{};
  std::vector<double> inf_norm_cols_{};

  // regularisation
  std::vector<double> total_reg_{};

  // values for static regularisation
  const Regul& regul_;

  // flag to stop computation
  bool flag_stop_ = false;

  const Log* log_;
  DataCollector& data_;

 public:
  void permute(const std::vector<Int>& iperm);
  void processSupernode(Int sn);

 public:
  Factorise(const Symbolic& S, const std::vector<Int>& rowsA,
            const std::vector<Int>& ptrA, const std::vector<double>& valA,
            const Regul& regul, const Log* log, DataCollector& data);

  bool run(Numeric& num);
};

}  // namespace hipo

#endif