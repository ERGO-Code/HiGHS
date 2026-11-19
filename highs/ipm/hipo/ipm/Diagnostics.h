#ifndef HIPO_DIAGNOSTICS_H
#define HIPO_DIAGNOSTICS_H

#include <sstream>
#include <vector>

#include "IterationHistory.h"
#include "IterationReport.h"
#include "Model.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

struct Iterate;

struct Diagnostics {
  const Model& model;

  // outcome of the residual reduction test in the last kBadResWindow
  // iterations: bit i is set if the test failed i iterations ago
  Int bad_res_history{};

  // data for progress test, for pinf, dinf, absolute gap and infeas_drop:
  // - values in the last kProgressMaxIter + 1 iterations in which the quantity
  //   was available, oldest first;
  // - whether the quantity was above its tolerance in the same iterations (for
  //   the gap, this uses the relative gap);
  // - best value and the iteration at which it was set.
  static constexpr Int kNumProgressQuantities = 4;
  std::vector<double> progress_history[kNumProgressQuantities];
  std::vector<HighsBool> above_tol_history[kNumProgressQuantities];
  double progress_ref[kNumProgressQuantities]{kHighsInf, kHighsInf, kHighsInf,
                                              kHighsInf};
  Int progress_iter[kNumProgressQuantities]{};

  explicit Diagnostics(const Model& model_input);

  // ===================================================================================
  // Detect stagnation and add to report:
  // - kBadResReduction, if residualReductionTest failed in at least
  //   kBadResMaxIter of the last kBadResWindow iterations;
  // - kLargeResIncrease, if residualReductionTest found a residual increased
  //   by more than a factor kResReductionLargeIncrease;
  // - kNoProgress and kOscillation, if the window and best-value parts of
  //   progressTest failed.
  // Must be called after the step, while it.delta still holds the direction.
  // ===================================================================================
  void stagnation(const Iterate& it, Int iter, double feas_tol, double opt_tol,
                  double crossover_tol, IterationReport& report,
                  std::stringstream& log_stream);

  // ===================================================================================
  // Add kDivergence to report if mu is larger than its best value by a factor
  // kDivergeTol, or, after iteration 5, an objective is larger than the other
  // by a factor kDivergeTol^2.
  // ===================================================================================
  void divergence(const Iterate& it, Int iter, IterationReport& report) const;

  // ===================================================================================
  // Detect numerical trouble in the linear solves and factorisation of the last
  // iteration, from its record, and add to report:
  // - kVeryInaccurateSolve, if cw_error > kVeryInaccurateSolveTol, otherwise
  //   kInaccurateSolve, if cw_error > kInaccurateSolveTol;
  // - kWrongInertia, if a 1x1 pivot has the wrong sign or the inertia differs
  //   from the expected one;
  // - kDynamicReg, if a pivot was perturbed beyond the static regularisation;
  // - kAbsorbedRows, if a row has its residual absorbed by the regularisation;
  // - kManyNearFree, for LPs, if there are more nearly free columns than rows
  //   (near optimality, up to m basic columns are nearly free).
  // ===================================================================================
  void numerical(const Iterate& it, IterationReport& report,
                 std::stringstream& log_stream) const;

  // ===================================================================================
  // Add kPoorCentrality to report if more than a fraction
  // kPoorCentralityFraction of the products xj * zj / mu are outside
  // [kSmallProduct, kLargeProduct] and a stepsize is below kPoorCentralityAlpha.
  // ===================================================================================
  void centrality(const Iterate& it, IterationReport& report,
                  std::stringstream& log_stream) const;

  // ===================================================================================
  // Check that the step reduced the residuals by at least a fraction of what
  // the stepsizes predict:
  //  ||rk_new|| <= max( (1 - eta * alpha) * ||rk_old||, floor )
  // with alpha_primal for k=1,2,3 and alpha_dual for k=4. Stepsizes and norms
  // before and after the step are taken from record; an iteration abandoned
  // before the step has zero stepsizes and is skipped.
  // Sets large_increase if a residual increased by more than a factor
  // kResReductionLargeIncrease.
  // Must be called after the step, while it.delta still holds the direction.
  // ===================================================================================
  bool residualReductionTest(const Iterate& it, const IterationRecord& record,
                             double feas_tol, bool& large_increase,
                             std::stringstream& log_stream) const;

  // ===================================================================================
  // Check progress of pinf, dinf, mu (measuring the gap) and infeas_drop.
  // The relative gap pdgap is used only to decide if the gap is above its
  // tolerance. The test fails if:
  // - none of the quantities, at any point in the last kProgressMaxIter
  //   iterations, was smaller by a factor kProgressFactor than
  //   kProgressMaxIter iterations ago, or
  // - none of the quantities improved its best value by a factor
  //   kProgressFactor in the last kProgressBestMaxIter iterations (this
  //   detects quantities that oscillate).
  // Only the quantities above their tolerance, or that were above their
  // tolerance at the start of the window, are tracked. If all quantities
  // are below their tolerance, but the ipm is still iterating (e.g. because
  // the kkt check failed), all quantities are tracked.
  // History and best values are recorded for all quantities, tracked or not,
  // and are reset only by clear(). infeas_drop is recorded only when it is
  // computed.
  // Sets no_window_progress and no_best_progress if the two parts failed.
  // ===================================================================================
  bool progressTest(const Iterate& it, Int iter, double feas_tol,
                    double opt_tol, double crossover_tol,
                    bool& no_window_progress, bool& no_best_progress,
                    std::stringstream& log_stream);

  // ===================================================================================
  // Clear the data used by the tests. The windows of the progress test restart
  // from iteration iter.
  // ===================================================================================
  void clear(Int iter);
};

}  // namespace hipo

#endif
