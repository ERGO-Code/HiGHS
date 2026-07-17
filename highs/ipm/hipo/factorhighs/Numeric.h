#ifndef FACTORHIGHS_NUMERIC_H
#define FACTORHIGHS_NUMERIC_H

#include <memory>
#include <vector>

#include "DataCollector.h"
#include "SolveDag.h"
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

  // per supernode, per block of columns: 1 iff the block contains any
  // non-identity swap, so that solves can skip the swap permutations
  std::vector<std::vector<uint8_t>> swap_flags_{};

  // workspace for the solve phase, sized once per factorisation to the
  // largest supernode leading dimension; avoids per-block allocations
  mutable std::vector<double> solve_work_{};

  // ===========================================================================
  // Schedule for the experimental parallel triangular solve (see
  // PHASE3_DESIGN.md). Built only when hipoTuning().parallel_solve is set.
  //
  // solve_schedule_[level][task] is an ascending list of supernodes. All
  // supernodes of level L depend only on supernodes of levels < L (forward
  // solve) resp. > L (backward solve). Tasks within a level are independent
  // up to the deferred forward-solve scatters, which are merged serially.
  // ===========================================================================
  std::vector<std::vector<std::vector<Int>>> solve_schedule_{};

  // per-task-slot buffers for the parallel solve (gemv workspace and the
  // deferred scatter contributions of the forward solve)
  mutable std::vector<std::vector<double>> task_work_{};
  mutable std::vector<std::vector<Int>> task_def_rows_{};
  mutable std::vector<std::vector<double>> task_def_vals_{};

  // Block-level DAG for parallel_solve_mode == 2 (see PHASE3_DESIGN.md),
  // plus the forward-solve stash: post-dtrsv (pre-inverse-swap) block values
  // of x, written by each supernode's solve task and read by its delivery
  // tasks so that deliveries need no swap handling of their own.
  SolveDag solve_dag_{};
  mutable std::vector<double> fwd_stash_{};

  // Compute swap_flags_ and size solve_work_; called by Factorise once the
  // factor has been handed over
  void finaliseFactor();

  // Build solve_schedule_ from the supernodal elimination tree
  void buildSolveSchedule();

  // information about 2x2 pivots
  std::vector<std::vector<double>> pivot_2x2_{};

  // symbolic object
  const Symbolic* S_;

  DataCollector* data_ = nullptr;

  friend class Factorise;

  // dynamic regularisation applied to the matrix
  std::vector<double> total_reg_{};

 public:
  // Full solve
  Int solve(std::vector<double>& x) const;
  void getReg(std::vector<double>& reg);
};

}  // namespace hipo

#endif
