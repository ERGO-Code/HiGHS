#ifndef FACTORHIGHS_FACTORISE_H
#define FACTORHIGHS_FACTORISE_H

#include <atomic>
#include <cmath>

#include "FactorHighsOptions.h"
#include "FactoriseMemory.h"
#include "Numeric.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

class Factorise {
 public:
  // matrix to factorise, permuted; managed outside of Factorise, so that it is
  // built once for each symbolic factorisation
  PermutedMatrix& matrix_;
  Int n_{};

  // symbolic factorisation
  const Symbolic& S_;

  // numerical factorisation, which receives the result
  Numeric& num_;

  // generated elements, aka Schur complements, when the tree is processed in
  // parallel. Managed outside of Factorise, so that it can be reused for all
  // ipm iterations.
  std::vector<std::vector<double>>& schur_contribution_;

  // columns of L, stored as dense supernodes
  // This memory is managed outside of Factorise, so that it can be reused for
  // all ipm iterations.
  FactorStorage& factor_;

  // one norm of the original matrix, sets the dynamic regularisation threshold
  double M_norm1_{};

  // The outputs below are the vectors of the Numeric object, written in place,
  // so that their memory is reused by all factorisations.

  // swaps of columns for each supernode, ordered locally within a block, one
  // entry per column
  std::vector<Int>& swaps_;

  // Information about 2x2 pivots, one entry per column.
  // If pivot_2x2[i] == 0, 1x1 pivot was used.
  // If pivot_2x2[i] != 0, 2x2 pivot was used and pivot_2x2[i] stores
  //  the off-diagonal pivot entry (of the 2x2 inverse).
  std::vector<double>& pivot_2x2_;

  // regularisation
  std::vector<double>& total_reg_;

  // 1x1 pivots before regularisation, infinite for 2x2 pivots
  std::vector<double>& unreg_pivots_;

  // static regularisation of each pivot, permuted. Managed outside of
  // Factorise, so that it can be reused for all ipm iterations.
  std::vector<double>& static_reg_;

  // flag to stop computation
  std::atomic<bool> flag_stop_{false};

  // the input is not compatible with the symbolic factorisation
  bool invalid_input_ = false;

  const Logger& logger_;
  DataCollector& data_;
  const FHoptions& FH_opt_;

  CliqueStack& stack_;

  // scratch memory of the dense kernels, managed outside of Factorise
  ScratchArena& scratch_;

  void processTreeSerial();
  void processTreeParallel();
  void processSupernode(Int sn);

 public:
  Factorise(const Symbolic& S, Numeric& num, Int n, Int nz, const Int* rows,
            const Int* ptr, const double* val, FactoriseMemory& memory,
            const FHoptions& FH_opt, const Logger& logger, DataCollector& data);

  bool run();
};

}  // namespace hipo

#endif
