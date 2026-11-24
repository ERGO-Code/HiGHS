#include "Diagnostics.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "Iterate.h"
#include "Parameters.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"

namespace hipo {

Diagnostics::Diagnostics(const Model& model_input) : model{model_input} {}

// Record of the last iteration. The empty record, used before the first one,
// has zero stepsizes and no signals.
static const IterationRecord& lastRecord(const Iterate& it) {
  static const IterationRecord empty_record{};
  return it.history.records.empty() ? empty_record : it.history.records.back();
}

void Diagnostics::stagnation(const Iterate& it, Int iter, double feas_tol,
                             double opt_tol, double crossover_tol,
                             IterationReport& report,
                             std::stringstream& log_stream) {
  const IterationRecord& record = lastRecord(it);

  bool large_increase = false;
  const bool reduction_ok =
      residualReductionTest(it, record, feas_tol, large_increase, log_stream);

  // record the outcome in the history of the last kBadResWindow iterations,
  // and count the failures
  const Int window_mask = (Int{1} << kBadResWindow) - 1;
  bad_res_history =
      ((bad_res_history << 1) | (reduction_ok ? 0 : 1)) & window_mask;
  Int bad_res_count = 0;
  for (Int h = bad_res_history; h != 0; h >>= 1) bad_res_count += h & 1;

  const bool bad_residuals = bad_res_count >= kBadResMaxIter;
  if (bad_residuals) {
    log_stream << "Bad reduction of residuals in " << bad_res_count
               << " of the last " << kBadResWindow << " iterations\n";
  }

  if (large_increase) {
    log_stream << "Residuals increased by more than a factor "
               << sci(kResReductionLargeIncrease, 0, 1) << "\n";
  }

  bool no_window_progress = false;
  bool no_best_progress = false;
  progressTest(it, iter, feas_tol, opt_tol, crossover_tol, no_window_progress,
               no_best_progress, log_stream);

  if (bad_residuals) report.add(IssueType::kBadResReduction);
  if (large_increase) report.add(IssueType::kLargeResIncrease);
  const bool gap_ok = gapTest(it, iter, feas_tol, opt_tol, log_stream);

  if (no_window_progress || !gap_ok) report.add(IssueType::kNoProgress);
  if (no_best_progress) report.add(IssueType::kOscillation);
}

void Diagnostics::divergence(const Iterate& it, Int iter,
                             IterationReport& report) const {
  const bool mu_is_large =
      it.best_mu > 0.0 ? it.mu > it.best_mu * kDivergeTol : false;
  const bool pobj_is_very_large =
      it.pobj < -std::max(std::abs(it.dobj) * kDivergeTol * kDivergeTol, 1.0);
  const bool dobj_is_very_large =
      it.dobj > std::max(std::abs(it.pobj) * kDivergeTol * kDivergeTol, 1.0);
  const bool objective_diverging =
      iter > 5 && (pobj_is_very_large || dobj_is_very_large);

  if (mu_is_large || objective_diverging) report.add(IssueType::kDivergence);
}

void Diagnostics::numerical(const Iterate& it, IterationReport& report,
                            std::stringstream& log_stream) const {
  const IterationRecord& record = lastRecord(it);

  if (record.cw_error > kVeryInaccurateSolveTol) {
    report.add(IssueType::kVeryInaccurateSolve);
    log_stream << "Very inaccurate solve, componentwise error "
               << sci(record.cw_error, 0, 1) << "\n";
  } else if (record.cw_error > kInaccurateSolveTol) {
    report.add(IssueType::kInaccurateSolve);
    log_stream << "Inaccurate solve, componentwise error "
               << sci(record.cw_error, 0, 1) << "\n";
  }

  const double wrong_thresh =
      std::max((double)kWrongInertiaMinPivots,
               kWrongInertiaFraction * record.num_pivots);
  if (std::max(record.num_wrong_sign, record.wrong_inertia) >= wrong_thresh) {
    report.add(IssueType::kWrongInertia);
    log_stream << "Wrong inertia, " << record.num_wrong_sign
               << " 1x1 pivots with wrong sign, inertia off by "
               << record.wrong_inertia << "\n";
  }

  if (record.num_absorbed_rows > 0 || record.num_absorbed_cols > 0) {
    report.add(IssueType::kAbsorbed);
    log_stream << record.num_absorbed_rows << " rows and "
               << record.num_absorbed_cols
               << " columns with residual absorbed by the regularisation\n";
  }

  if (record.num_near_free > model.m()) {
    report.add(IssueType::kManyNearFree);
    log_stream << record.num_near_free << " nearly free columns, " << model.m()
               << " rows\n";
  }
}

void Diagnostics::centrality(const Iterate& it, IterationReport& report,
                             std::stringstream& log_stream) const {
  const IterationRecord& record = lastRecord(it);

  Int num_finite = 0;
  for (Int j = 0; j < model.n(); ++j)
    num_finite += model.hasLb(j) + model.hasUb(j);

  const Int num_outliers = record.num_small_prod + record.num_large_prod;
  const bool many_outliers =
      num_outliers > kPoorCentralityFraction * num_finite;
  const bool short_step =
      std::min(record.alpha_primal, record.alpha_dual) < kPoorCentralityAlpha;

  if (many_outliers && short_step) {
    report.add(IssueType::kPoorCentrality);
    log_stream << "Poor centrality, " << record.num_small_prod << " small and "
               << record.num_large_prod << " large products of " << num_finite
               << ", stepsizes " << sci(record.alpha_primal, 0, 1) << " "
               << sci(record.alpha_dual, 0, 1) << "\n";
  }
}

bool Diagnostics::residualReductionTest(const Iterate& it,
                                        const IterationRecord& record,
                                        double feas_tol, bool& large_increase,
                                        std::stringstream& log_stream) const {
  const double alpha_primal = record.alpha_primal;
  const double alpha_dual = record.alpha_dual;
  const double* res_norm_old = record.res_norm_before;

  large_increase = false;

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
  double res_norm_new[4] = {record.res_norm_after[0], record.res_norm_after[1],
                            record.res_norm_after[2], record.res_norm_after[3]};

  // For QPs, r4 depends also on x, so after the step
  //  r4_new = (1 - alpha_dual) * r4_old
  //           + (alpha_primal - alpha_dual) * Q * dx
  // The second term is a known change that is not related to the reduction of
  // the residual, so it is removed from r4_new before testing.
  // Stepsizes should be equal for QPs and the code below should never run, so
  // this is just a safety measure.
  if (model.qp() && alpha_primal != alpha_dual) {
    std::vector<double> res4_corrected = it.res.r4;
    model.Q().alphaProductPlusY(-(alpha_primal - alpha_dual), it.delta.x,
                                res4_corrected);
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
      log_stream << "Bad reduction of residual " << k + 1 << ", old "
                 << sci(res_norm_old[k], 0, 1) << ", new "
                 << sci(res_norm_new[k], 0, 1) << ", relaxed prediction "
                 << sci(res_norm_relaxed[k], 0, 1) << ", alpha "
                 << sci(alpha[k], 0, 1) << "\n";
    }

    if (res_norm_new[k] > kResReductionLargeIncrease * res_norm_old[k])
      large_increase = true;
  }

  return passed;
}

bool Diagnostics::progressTest(const Iterate& it, Int iter, double feas_tol,
                               double opt_tol, double crossover_tol,
                               bool& no_window_progress, bool& no_best_progress,
                               std::stringstream& log_stream) {
  const Int num = kNumProgressQuantities;
  // The progress of the gap is measured with mu. The relative gap saturates at
  // 2 when one objective is much larger than the other, and then does not show
  // any progress. The absolute gap |pobj - dobj| can pass close to zero when
  // the objectives cross, which happens at infeasible points, and this gives a
  // best value that later iterates cannot improve. mu is non-negative, does not
  // saturate, and its ratios are scale-invariant. Whether the gap is above its
  // tolerance is decided with the relative gap instead.
  const double value[num] = {it.pinf, it.dinf, it.mu, it.infeas_drop};
  const double tol_value[num] = {it.pinf, it.dinf, it.pdgap, it.infeas_drop};
  const double tol[num] = {feas_tol, feas_tol, opt_tol, crossover_tol};
  const char* name[num] = {"pinf", "dinf", "mu", "drop"};

  // infeas_drop is computed only in some cases, otherwise it is zero
  const bool available[num] = {true, true, true, it.infeas_drop > 0.0};

  // Record the history and best value of each available quantity, whether it
  // is tracked or not. This data is never reset, except by clear(), so that a
  // quantity that moves above and below its tolerance cannot pass the tests
  // just by restarting its window.
  for (Int k = 0; k < num; ++k) {
    if (!available[k]) continue;

    std::vector<double>& history = progress_history[k];
    std::vector<HighsBool>& above = above_tol_history[k];
    history.push_back(value[k]);
    above.push_back(tol_value[k] > tol[k]);
    if (static_cast<Int>(history.size()) > kProgressMaxIter + 1) {
      history.erase(history.begin());
      above.erase(above.begin());
    }

    if (value[k] < kProgressFactor * progress_ref[k]) {
      progress_ref[k] = value[k];
      progress_iter[k] = iter;
    }
  }

  // A quantity is tracked only if it is above its tolerance, or if it was above
  // its tolerance at the start of the window. Otherwise, a quantity that has
  // just reached its tolerance would be ignored, even though it made progress
  // in the window.
  bool tracked[num];
  bool any_tracked = false;
  for (Int k = 0; k < num; ++k) {
    const std::vector<HighsBool>& above = above_tol_history[k];
    const bool was_above = !above.empty() && above.front();
    tracked[k] = available[k] && (tol_value[k] > tol[k] || was_above);
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

  no_window_progress = !window_passed;
  no_best_progress = !best_passed;

  return window_passed && best_passed;
}

bool Diagnostics::gapTest(const Iterate& it, Int iter, double feas_tol,
                          double opt_tol, std::stringstream& log_stream) {
  if (it.pinf > feas_tol || it.dinf > feas_tol) {
    gap_ref = kHighsInf;
    return true;
  }
  if (it.pdgap <= opt_tol) return true;

  if (it.pdgap < kProgressFactor * gap_ref) {
    gap_ref = it.pdgap;
    gap_iter = iter;
  }
  if (iter - gap_iter < kGapProgressMaxIter) return true;

  log_stream << "No progress of pdgap at feasible points in last "
             << kGapProgressMaxIter << " iterations: pdgap "
             << sci(it.pdgap, 0, 1) << " (best " << sci(gap_ref, 0, 1) << ")\n";
  return false;
}

void Diagnostics::clear(Int iter) {
  bad_res_history = 0;
  gap_ref = kHighsInf;
  gap_iter = iter;

  for (Int k = 0; k < kNumProgressQuantities; ++k) {
    progress_history[k].clear();
    above_tol_history[k].clear();
    progress_ref[k] = kHighsInf;
    progress_iter[k] = iter;
  }
}

}  // namespace hipo
