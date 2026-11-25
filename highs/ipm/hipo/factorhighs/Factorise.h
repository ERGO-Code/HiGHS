#ifndef FACTORHIGHS_FACTORISE_H
#define FACTORHIGHS_FACTORISE_H

#include <atomic>
#include <cmath>

#include "CliqueStack.h"
#include "FactorHighsOptions.h"
#include "Numeric.h"
#include "SnColumns.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

class ScratchArena;

class Factorise {
 public:
  // matrix to factorise
  std::vector<Int> rowsM_{};
  std::vector<Int> ptrM_{};
  std::vector<double> valM_{};
  Int n_{};
  Int nzM_{};

  // symbolic factorisation
  const Symbolic& S_;

  // numerical factorisation, which receives the result
  Numeric& num_;

  // children in supernodal elimination tree
  std::vector<Int> first_child_{};
  std::vector<Int> next_child_{};

  // reverse linked lists of chidlren
  std::vector<Int> first_child_reverse_{};
  std::vector<Int> next_child_reverse_{};

  // generated elements, aka Schur complements.
  std::vector<std::vector<double>> schur_contribution_{};

  // columns of L, stored as dense supernodes
  // This memory is managed outside of Factorise, so that it can be reused for
  // all ipm iterations.
  SnColumns& sn_columns_;

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

  // static regularisation of each pivot, permuted
  std::vector<double> static_reg_{};

  // flag to stop computation
  std::atomic<bool> flag_stop_{false};

  const Logger& logger_;
  DataCollector& data_;
  const FHoptions& FH_opt_;

  CliqueStack* stack_;

  // scratch memory of the dense kernels, managed outside of Factorise
  ScratchArena& scratch_;

  void permute(const std::vector<Int>& iperm);
  void processTreeSerial();
  void processTreeParallel();
  void processSupernode(Int sn);

 public:
  Factorise(const Symbolic& S, Numeric& num, Int n, Int nz, const Int* rowsM,
            const Int* ptrM, const double* valM,
            const std::vector<double>& static_reg, const FHoptions& FH_opt,
            const Logger& logger, DataCollector& data,
            SnColumns& sn_columns, CliqueStack* stack,
            ScratchArena& scratch);

  bool run();
};

}  // namespace hipo

#endif
