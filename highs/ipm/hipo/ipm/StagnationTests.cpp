#include "StagnationTests.h"

#include <algorithm>
#include <limits>

#include "Iterate.h"
#include "Parameters.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"

namespace hipo {

StagnationTests::StagnationTests(const Model& model_input)
    : model{model_input} {}

StagnationReport StagnationTests::stagnation(const Iterate& it,
                                             double alpha_primal,
                                             double alpha_dual, Int iter,
                                             double feas_tol, double opt_tol,
                                             double crossover_tol,
                                             std::stringstream& log_stream) {
  StagnationReport report;
  report.iter = iter;
  report.omega = it.data.record.empty() ? 0.0 : it.data.record.back().omega;
  report.free_vars_diverging = free_bounds_changed_iter >= kFreeVarsDivergeIter;

  const bool reduction_ok = residualReductionTest(it, alpha_primal, alpha_dual,
                                                  feas_tol, report, log_stream);

  // record the outcome in the history of the last kBadResWindow iterations,
  // and count the failures
  const Int window_mask = (Int{1} << kBadResWindow) - 1;
  bad_res_history =
      ((bad_res_history << 1) | (reduction_ok ? 0 : 1)) & window_mask;
  for (Int h = bad_res_history; h != 0; h >>= 1) report.bad_res_count += h & 1;

  const bool bad_residuals = report.bad_res_count >= kBadResMaxIter;
  if (bad_residuals) {
    log_stream << "Bad reduction of residuals in " << report.bad_res_count
               << " of the last " << kBadResWindow << " iterations\n";
  }

  if (report.large_increase) {
    log_stream << "Residuals increased by more than a factor "
               << sci(kResReductionLargeIncrease, 0, 1) << "\n";
  }

  progressTest(it, iter, feas_tol, opt_tol, crossover_tol, report, log_stream);

  if (bad_residuals || report.large_increase)
    report.cause = StagnationCause::kNumerical;
  else if (report.no_window_progress)
    report.cause = StagnationCause::kNoProgress;
  else if (report.no_best_progress)
    report.cause = StagnationCause::kOscillation;

  report.stagnation = report.cause != StagnationCause::kNone;

  return report;
}

void StagnationTests::saveResidualNorms(const Iterate& it) {
  res_norm_old[0] = infNorm(it.res.r1);
  res_norm_old[1] = infNorm(it.res.r2);
  res_norm_old[2] = infNorm(it.res.r3);
  res_norm_old[3] = infNorm(it.res.r4);
}

bool StagnationTests::residualReductionTest(
    const Iterate& it, double alpha_primal, double alpha_dual, double feas_tol,
    StagnationReport& report, std::stringstream& log_stream) const {
  report.large_increase = false;
  report.bad_res_blocks = 0;

  // below these values, the residuals are small enough that their reduction is
  // not needed
  const double tol_p = feas_tol * (1.0 + model.normScaledRhs());
  const double tol_d = feas_tol * (1.0 + model.normScaledObj());

  // Below these values, the residuals cannot be computed accurately, because of
  // rounding errors, which are estimated as eps times the size of the terms
  // that make up the residuals:
  //  r1 = b - A * x
  //  r2 = lower - x + xl
  //  r3 = upper - x - xu
  //  r4 = c - A^T * y - zl + zu + Q * x
  const double eps = std::numeric_limits<double>::epsilon();
  const double norm_x = infNorm(it.x);
  const double noise[4] = {model.normA() * norm_x + infNorm(model.b()),
                           norm_x + infNorm(it.xl) + model.normScaledRhs(),
                           norm_x + infNorm(it.xu) + model.normScaledRhs(),
                           infNorm(model.c()) + model.normAt() * infNorm(it.y) +
                               infNorm(it.zl) + infNorm(it.zu) +
                               model.normQ() * norm_x};

  const double alpha[4] = {alpha_primal, alpha_primal, alpha_primal,
                           alpha_dual};
  const double tol[4] = {tol_p, tol_p, tol_p, tol_d};
  double floor[4];
  for (Int k = 0; k < 4; ++k)
    floor[k] = std::max(tol[k], kResReductionNoiseFactor * eps * noise[k]);
  double res_norm_new[4] = {infNorm(it.res.r1), infNorm(it.res.r2),
                            infNorm(it.res.r3), infNorm(it.res.r4)};

  // For QPs, r4 depends also on x, so after the step
  //  r4_new = (1 - alpha_dual) * r4_old
  //           + (alpha_primal - alpha_dual) * Q * dx
  // The second term is a known change that is not related to the reduction of
  // the residual, so it is removed from r4_new before testing.
  // Stepsizes should be equal for QPs and the code below should never run, so
  // this is just a safety measure.
  if (model.qp() && alpha_primal != alpha_dual) {
    std::vector<double> res4_corrected = it.res.r4;
    model.Q().alphaProductPlusY(-model.sense() * (alpha_primal - alpha_dual),
                                it.delta.x, res4_corrected);
    res_norm_new[3] = infNorm(res4_corrected);
  }

  // Relaxed prediction of the norms of the residuals after the step, i.e.
  // only a fraction eta of the reduction predicted by the stepsizes is
  // required:
  //  rk_relaxed = (1 - eta * alpha) * rk
  const double eta = kResReductionFraction;
  double res_norm_relaxed[4];
  for (Int k = 0; k < 4; ++k)
    res_norm_relaxed[k] = (1.0 - eta * alpha[k]) * res_norm_old[k];

  bool passed = true;
  for (Int k = 0; k < 4; ++k) {
    // test is meaningless if the step is tiny or the residual is already small
    if (alpha[k] < kResReductionMinAlpha || res_norm_old[k] <= floor[k])
      continue;

    const double thresh = std::max(res_norm_relaxed[k], floor[k]);

    if (res_norm_new[k] > thresh) {
      passed = false;
      report.bad_res_blocks |= Int{1} << k;
      log_stream << "Bad reduction of residual " << k + 1 << ", old "
                 << sci(res_norm_old[k], 0, 1) << ", new "
                 << sci(res_norm_new[k], 0, 1) << ", relaxed prediction "
                 << sci(res_norm_relaxed[k], 0, 1) << ", alpha "
                 << sci(alpha[k], 0, 1) << "\n";
    }

    if (res_norm_new[k] > kResReductionLargeIncrease * res_norm_old[k])
      report.large_increase = true;
  }

  return passed;
}

bool StagnationTests::progressTest(const Iterate& it, Int iter, double feas_tol,
                                   double opt_tol, double crossover_tol,
                                   StagnationReport& report,
                                   std::stringstream& log_stream) {
  const Int num = kNumProgressQuantities;
  const double value[num] = {it.pinf, it.dinf, it.pdgap, it.infeas_drop};
  const double tol[num] = {feas_tol, feas_tol, opt_tol, crossover_tol};
  const char* name[num] = {"pinf", "dinf", "pdgap", "drop"};

  // infeas_drop is computed only in some cases, otherwise it is zero
  const bool available[num] = {true, true, true, it.infeas_drop > 0.0};

  // Record the history and best value of each available quantity, whether it
  // is tracked or not. This data is never reset, except by clear(), so that a
  // quantity that moves above and below its tolerance cannot pass the tests
  // just by restarting its window.
  for (Int k = 0; k < num; ++k) {
    if (!available[k]) continue;

    std::vector<double>& history = progress_history[k];
    history.push_back(value[k]);
    if (static_cast<Int>(history.size()) > kProgressMaxIter + 1)
      history.erase(history.begin());

    if (value[k] < kProgressFactor * progress_ref[k]) {
      progress_ref[k] = value[k];
      progress_iter[k] = iter;
    }
  }

  // a quantity is tracked only if it is above its tolerance
  bool tracked[num];
  bool any_tracked = false;
  for (Int k = 0; k < num; ++k) {
    tracked[k] = available[k] && value[k] > tol[k];
    if (tracked[k]) any_tracked = true;
  }

  // If no quantity is above its tolerance, but the ipm is still iterating,
  // track all of them.
  if (!any_tracked)
    for (Int k = 0; k < num; ++k) tracked[k] = available[k];

  // The tests pass if any tracked quantity passes them.
  bool window_passed = false;
  bool best_passed = false;
  for (Int k = 0; k < num; ++k) {
    if (!tracked[k]) continue;
    const std::vector<double>& history = progress_history[k];

    // Window test: check if, at some point in the last kProgressMaxIter
    // iterations, the quantity was smaller by a factor kProgressFactor than
    // kProgressMaxIter iterations ago. If there is not enough history yet, the
    // test passes.
    if (static_cast<Int>(history.size()) <= kProgressMaxIter) {
      window_passed = true;
    } else {
      const double window_min =
          *std::min_element(history.begin() + 1, history.end());
      if (window_min < kProgressFactor * history.front()) window_passed = true;
    }

    // Best test: check if the best value was improved by a factor
    // kProgressFactor in the last kProgressBestMaxIter iterations. This
    // detects quantities that oscillate, which may pass the window test.
    if (iter - progress_iter[k] < kProgressBestMaxIter) best_passed = true;
  }

  if (!window_passed) {
    log_stream << "No progress in last " << kProgressMaxIter << " iterations:";
    for (Int k = 0; k < num; ++k) {
      if (!tracked[k]) continue;
      log_stream << " " << name[k] << " " << sci(value[k], 0, 1) << " (was "
                 << sci(progress_history[k].front(), 0, 1) << ")";
    }
    log_stream << "\n";
  }

  if (!best_passed) {
    log_stream << "No improvement of best values in last "
               << kProgressBestMaxIter << " iterations:";
    for (Int k = 0; k < num; ++k) {
      if (!tracked[k]) continue;
      log_stream << " " << name[k] << " " << sci(value[k], 0, 1) << " (best "
                 << sci(progress_ref[k], 0, 1) << ")";
    }
    log_stream << "\n";
  }

  report.no_window_progress = !window_passed;
  report.no_best_progress = !best_passed;

  return window_passed && best_passed;
}

void StagnationTests::clear(Int iter) {
  bad_res_history = 0;
  free_bounds_changed_iter = 0;

  // with zero old norms, the residual reduction test is skipped until the
  // norms are saved again
  for (Int k = 0; k < 4; ++k) res_norm_old[k] = 0.0;

  for (Int k = 0; k < kNumProgressQuantities; ++k) {
    progress_history[k].clear();
    progress_ref[k] = kHighsInf;
    progress_iter[k] = iter;
  }
}

}  // namespace hipo
