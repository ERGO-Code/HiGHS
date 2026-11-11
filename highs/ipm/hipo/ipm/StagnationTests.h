#ifndef HIPO_STAGNATION_TESTS_H
#define HIPO_STAGNATION_TESTS_H

#include <sstream>
#include <vector>

#include "Model.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

struct Iterate;

// Most likely cause of stagnation, used to decide how to recover.
enum class StagnationCause {
  kNone,         // no stagnation
  kNumerical,    // residuals not reduced as predicted by the stepsizes
  kNoProgress,   // no quantity made progress in the last iterations
  kOscillation,  // no quantity improved its best value for a long time
};

// Outcome of the stagnation tests at one iteration, with the information
// needed to decide how to recover.
struct StagnationReport {
  bool stagnation = false;
  StagnationCause cause = StagnationCause::kNone;
  Int iter = 0;

  // residual reduction test
  Int bad_res_count = 0;   // failures in the last kBadResWindow iterations
  Int bad_res_blocks = 0;  // bit k set if residual k+1 failed at this iter
  bool large_increase = false;

  // largest backward error of the linear solves in the last iteration, to
  // distinguish failures of the linear solver from other causes
  double omega = 0.0;

  // progress test
  bool no_window_progress = false;
  bool no_best_progress = false;

  // the bounds of free variables have been changed in the last
  // kFreeVarsDivergeIter iterations
  bool free_vars_diverging = false;
};

struct StagnationTests {
  const Model& model;

  // data for residual reduction test
  double res_norm_old[4]{};

  // infinity norms of A, A^T and Q, used to estimate the rounding errors in the
  // computation of the residuals
  double norm_A{}, norm_At{}, norm_Q{};

  // outcome of the residual reduction test in the last kBadResWindow
  // iterations: bit i is set if the test failed i iterations ago
  Int bad_res_history{};

  // number of consecutive iterations in which the bounds of free variables
  // were changed
  Int free_bounds_changed_iter{};

  // data for progress test, for pinf, dinf, pdgap and infeas_drop:
  // - values in the last kProgressMaxIter + 1 iterations in which the quantity
  //   was available, oldest first;
  // - best value and the iteration at which it was set.
  static constexpr Int kNumProgressQuantities = 4;
  std::vector<double> progress_history[kNumProgressQuantities];
  double progress_ref[kNumProgressQuantities]{kHighsInf, kHighsInf, kHighsInf,
                                              kHighsInf};
  Int progress_iter[kNumProgressQuantities]{};

  explicit StagnationTests(const Model& model_input);

  // ===================================================================================
  // Detect stagnation, which happens if any of the following is true:
  // - residualReductionTest failed in at least kBadResMaxIter of the last
  //   kBadResWindow iterations;
  // - residualReductionTest found a residual increased by more than a factor
  //   kResReductionLargeIncrease;
  // - progressTest failed.
  // Returns a report with the outcome of the tests and the most likely cause
  // of stagnation.
  // Must be called after the step, while it.delta still holds the direction.
  // ===================================================================================
  StagnationReport stagnation(const Iterate& it, double alpha_primal,
                              double alpha_dual, Int iter, double feas_tol,
                              double opt_tol, double crossover_tol,
                              std::stringstream& log_stream);

  // ===================================================================================
  // Save the norms of residuals 1,2,3,4 used as rhs of the Newton system, to be
  // used by residualReductionTest. Must be called before the correctors
  // overwrite the residuals.
  // ===================================================================================
  void saveResidualNorms(const Iterate& it);

  // ===================================================================================
  // Check that the step reduced the residuals by at least a fraction of what
  // the stepsizes predict:
  //  ||rk_new|| <= max( (1 - eta * alpha) * ||rk_old||, floor )
  // with alpha_primal for k=1,2,3 and alpha_dual for k=4.
  // In report, sets bad_res_blocks, and large_increase if a residual increased
  // by more than a factor kResReductionLargeIncrease.
  // Must be called after the step, while it.delta still holds the direction.
  // ===================================================================================
  bool residualReductionTest(const Iterate& it, double alpha_primal,
                             double alpha_dual, double feas_tol,
                             StagnationReport& report,
                             std::stringstream& log_stream) const;

  // ===================================================================================
  // Check progress of pinf, dinf, pdgap and infeas_drop. The test fails if:
  // - none of the quantities, at any point in the last kProgressMaxIter
  //   iterations, was smaller by a factor kProgressFactor than
  //   kProgressMaxIter iterations ago, or
  // - none of the quantities improved its best value by a factor
  //   kProgressFactor in the last kProgressBestMaxIter iterations (this
  //   detects quantities that oscillate).
  // Only the quantities above their tolerance are tracked. If all quantities
  // are below their tolerance, but the ipm is still iterating (e.g. because
  // the kkt check failed), all quantities are tracked.
  // History and best values are recorded for all quantities, tracked or not,
  // and are reset only by clear(). infeas_drop is recorded only when it is
  // computed.
  // In report, sets no_window_progress and no_best_progress.
  // ===================================================================================
  bool progressTest(const Iterate& it, Int iter, double feas_tol,
                    double opt_tol, double crossover_tol,
                    StagnationReport& report, std::stringstream& log_stream);

  // ===================================================================================
  // Clear the data used by the tests. The windows of the progress test restart
  // from iteration iter.
  // ===================================================================================
  void clear(Int iter);
};

}  // namespace hipo

#endif
