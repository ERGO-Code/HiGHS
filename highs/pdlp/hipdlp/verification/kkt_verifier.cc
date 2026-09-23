#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "lp_data/HighsLp.h"

namespace hipdlp {
namespace verification {
namespace {

bool hasLower(const double lower) { return lower > -kHighsInf; }

bool hasUpper(const double upper) { return upper < kHighsInf; }

bool finiteVector(const std::vector<double>& values) {
  for (const double value : values) {
    if (!std::isfinite(value)) return false;
  }
  return true;
}

bool validBounds(const std::vector<double>& lower,
                 const std::vector<double>& upper) {
  if (lower.size() != upper.size()) return false;
  for (size_t index = 0; index < lower.size(); ++index) {
    if (std::isnan(lower[index]) || std::isnan(upper[index]) ||
        lower[index] >= kHighsInf || upper[index] <= -kHighsInf ||
        lower[index] > upper[index]) {
      return false;
    }
  }
  return true;
}

bool validSparseStorage(const HighsSparseMatrix& matrix) {
  if (!matrix.formatOk() || matrix.num_col_ < 0 || matrix.num_row_ < 0)
    return false;

  const HighsInt num_vector =
      matrix.isColwise() ? matrix.num_col_ : matrix.num_row_;
  const HighsInt index_limit =
      matrix.isColwise() ? matrix.num_row_ : matrix.num_col_;
  if (matrix.start_.size() < static_cast<size_t>(num_vector) + 1 ||
      matrix.start_[0] != 0) {
    return false;
  }
  for (HighsInt vector = 0; vector < num_vector; ++vector) {
    if (matrix.start_[vector] < 0 ||
        matrix.start_[vector] > matrix.start_[vector + 1]) {
      return false;
    }
  }
  if (matrix.format_ == MatrixFormat::kRowwisePartitioned) {
    if (matrix.p_end_.size() < static_cast<size_t>(num_vector)) return false;
    for (HighsInt vector = 0; vector < num_vector; ++vector) {
      if (matrix.p_end_[vector] < matrix.start_[vector] ||
          matrix.p_end_[vector] > matrix.start_[vector + 1]) {
        return false;
      }
    }
  }

  const HighsInt num_nonzero = matrix.start_[num_vector];
  if (num_nonzero < 0 ||
      matrix.index_.size() < static_cast<size_t>(num_nonzero) ||
      matrix.value_.size() < static_cast<size_t>(num_nonzero)) {
    return false;
  }
  for (HighsInt nonzero = 0; nonzero < num_nonzero; ++nonzero) {
    if (matrix.index_[nonzero] < 0 || matrix.index_[nonzero] >= index_limit) {
      return false;
    }
  }
  return true;
}

void updateMaximum(const double value, const HighsInt index, double& maximum,
                   HighsInt& worst_index) {
  if (value > maximum) {
    maximum = value;
    worst_index = index;
  }
}

void updatePrimalViolation(const double lower, const double upper,
                           const double value, const HighsInt index,
                           double& raw_maximum, HighsInt& worst_index) {
  double violation = 0.0;
  if (hasLower(lower) && value < lower) {
    violation = lower - value;
  } else if (hasUpper(upper) && value > upper) {
    violation = value - upper;
  }
  updateMaximum(violation, index, raw_maximum, worst_index);
}

void updateActiveBoundNorm(const double lower, const double upper,
                           const double value, const double near_tolerance,
                           double& bound_norm) {
  if (hasLower(lower) && std::fabs(value - lower) <= near_tolerance) {
    bound_norm = std::max(bound_norm, std::fabs(lower));
  } else if (hasUpper(upper) && std::fabs(value - upper) <= near_tolerance) {
    bound_norm = std::max(bound_norm, std::fabs(upper));
  }
}

double relativePrimalViolation(const double lower, const double upper,
                               const double value,
                               const double active_bound_norm) {
  double violation = 0.0;
  double bound_measure = active_bound_norm;
  if (hasLower(lower) && value < lower) {
    violation = lower - value;
    bound_measure = std::max(bound_measure, std::fabs(lower));
  } else if (hasUpper(upper) && value > upper) {
    violation = value - upper;
    bound_measure = std::max(bound_measure, std::fabs(upper));
  }
  return violation / (1.0 + bound_measure);
}

double dualSignViolation(const double lower, const double upper,
                         const double value, const double adjusted_dual,
                         const double near_tolerance) {
  const bool lower_finite = hasLower(lower);
  const bool upper_finite = hasUpper(upper);
  if (!lower_finite && !upper_finite) return std::fabs(adjusted_dual);
  if (lower_finite && !upper_finite) return std::max(0.0, -adjusted_dual);
  if (!lower_finite && upper_finite) return std::max(0.0, adjusted_dual);
  if (upper - lower <= near_tolerance) return 0.0;
  const double midpoint = 0.5 * lower + 0.5 * upper;
  return value < midpoint ? std::max(0.0, -adjusted_dual)
                          : std::max(0.0, adjusted_dual);
}

double complementarityResidual(const double lower, const double upper,
                               const double value, const double adjusted_dual) {
  if (adjusted_dual > 0.0) {
    if (!hasLower(lower)) return std::fabs(adjusted_dual);
    return std::fabs(adjusted_dual * (value - lower));
  }
  if (adjusted_dual < 0.0) {
    if (!hasUpper(upper)) return std::fabs(adjusted_dual);
    return std::fabs(adjusted_dual * (upper - value));
  }
  return 0.0;
}

class ObjectiveSum {
 public:
  bool add(const double term) {
    if (!std::isfinite(term)) return false;
    double combined = sum_ + (exponent_ == 0 ? term : std::scalbn(term, -exponent_));
    if (!std::isfinite(combined)) {
      sum_ *= 0.5;
      ++exponent_;
      combined = sum_ + std::scalbn(term, -exponent_);
    }
    sum_ = combined;
    return std::isfinite(sum_);
  }

  double value() const { return std::scalbn(sum_, exponent_); }

 private:
  double sum_ = 0.0;
  int exponent_ = 0;
};

bool addDualObjectiveContribution(const double lower, const double upper,
                                  const double raw_dual,
                                  const double adjusted_dual,
                                  ObjectiveSum& dual_objective) {
  const bool lower_finite = hasLower(lower);
  const bool upper_finite = hasUpper(upper);
  if (!lower_finite && !upper_finite) return true;
  double bound;
  if (!upper_finite) {
    bound = lower;
  } else if (!lower_finite) {
    bound = upper;
  } else {
    bound = adjusted_dual >= 0.0 ? lower : upper;
  }
  return dual_objective.add(bound * raw_dual);
}

bool finiteMetrics(const KktMetrics& metrics) {
  return std::isfinite(metrics.max_primal_violation) &&
         std::isfinite(metrics.max_dual_violation) &&
         std::isfinite(metrics.max_stationarity_residual) &&
         std::isfinite(metrics.max_row_activity_residual) &&
         std::isfinite(metrics.max_complementarity) &&
         std::isfinite(metrics.primal_objective) &&
         std::isfinite(metrics.dual_objective) &&
         std::isfinite(metrics.absolute_objective_gap) &&
         std::isfinite(metrics.relative_primal_residual) &&
         std::isfinite(metrics.relative_dual_residual) &&
         std::isfinite(metrics.relative_objective_gap) &&
         std::isfinite(metrics.absolute_complementarity) &&
         std::isfinite(metrics.relative_complementarity);
}

}  // namespace

KktVerificationResult verifyKkt(const HighsLp& lp,
                                const HighsSolution& solution,
                                const double tolerance) {
  KktVerificationResult result;
  if (!std::isfinite(tolerance) || tolerance < 0.0) {
    result.error = "KKT tolerance must be finite and nonnegative";
    return result;
  }
  if (lp.num_col_ < 0 || lp.num_row_ < 0 ||
      lp.col_cost_.size() != static_cast<size_t>(lp.num_col_) ||
      lp.col_lower_.size() != static_cast<size_t>(lp.num_col_) ||
      lp.col_upper_.size() != static_cast<size_t>(lp.num_col_) ||
      lp.row_lower_.size() != static_cast<size_t>(lp.num_row_) ||
      lp.row_upper_.size() != static_cast<size_t>(lp.num_row_) ||
      lp.a_matrix_.num_col_ != lp.num_col_ ||
      lp.a_matrix_.num_row_ != lp.num_row_ ||
      !validSparseStorage(lp.a_matrix_)) {
    result.error = "LP dimensions are inconsistent";
    return result;
  }
  if (!solution.value_valid || !solution.dual_valid ||
      solution.col_value.size() != static_cast<size_t>(lp.num_col_) ||
      solution.col_dual.size() != static_cast<size_t>(lp.num_col_) ||
      solution.row_value.size() != static_cast<size_t>(lp.num_row_) ||
      solution.row_dual.size() != static_cast<size_t>(lp.num_row_)) {
    result.error = "A complete primal-dual solution is required";
    return result;
  }
  if (!std::isfinite(lp.offset_) || !finiteVector(lp.col_cost_) ||
      !finiteVector(lp.a_matrix_.value_) ||
      !validBounds(lp.col_lower_, lp.col_upper_) ||
      !validBounds(lp.row_lower_, lp.row_upper_) ||
      !finiteVector(solution.col_value) || !finiteVector(solution.col_dual) ||
      !finiteVector(solution.row_value) || !finiteVector(solution.row_dual)) {
    result.error = "LP or solution contains a non-finite or invalid value";
    return result;
  }

  const double objective_sense = static_cast<HighsInt>(lp.sense_);
  if (objective_sense != 1.0 && objective_sense != -1.0) {
    result.error = "LP objective sense is invalid";
    return result;
  }
  result.input_valid = true;

  std::vector<double> row_activity;
  std::vector<double> transpose_product;
  lp.a_matrix_.product(row_activity, solution.col_value);
  lp.a_matrix_.productTranspose(transpose_product, solution.row_dual);
  if (row_activity.size() != static_cast<size_t>(lp.num_row_) ||
      transpose_product.size() != static_cast<size_t>(lp.num_col_) ||
      !finiteVector(row_activity) || !finiteVector(transpose_product)) {
    result.error = "Sparse products produced non-finite values";
    return result;
  }

  double max_col_primal_violation = 0.0;
  double max_row_primal_violation = 0.0;
  const double near_tolerance = std::sqrt(tolerance);
  double active_bound_norm = 0.0;
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    updateActiveBoundNorm(lp.col_lower_[col], lp.col_upper_[col],
                          solution.col_value[col], near_tolerance,
                          active_bound_norm);
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    updateActiveBoundNorm(lp.row_lower_[row], lp.row_upper_[row],
                          row_activity[row], near_tolerance, active_bound_norm);
  }
  double max_relative_primal_violation = 0.0;
  double max_relative_activity_residual = 0.0;
  HighsInt worst_row_bound_violation = -1;
  HighsInt worst_row_activity_residual = -1;
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    updatePrimalViolation(lp.col_lower_[col], lp.col_upper_[col],
                          solution.col_value[col], col,
                          max_col_primal_violation, result.worst_primal_col);
    max_relative_primal_violation = std::max(
        max_relative_primal_violation,
        relativePrimalViolation(lp.col_lower_[col], lp.col_upper_[col],
                                solution.col_value[col], active_bound_norm));
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    updatePrimalViolation(lp.row_lower_[row], lp.row_upper_[row],
                          row_activity[row], row, max_row_primal_violation,
                          worst_row_bound_violation);
    max_relative_primal_violation =
        std::max(max_relative_primal_violation,
                 relativePrimalViolation(lp.row_lower_[row], lp.row_upper_[row],
                                         row_activity[row], active_bound_norm));
    const double activity_residual =
        std::fabs(row_activity[row] - solution.row_value[row]);
    updateMaximum(activity_residual, row,
                  result.metrics.max_row_activity_residual,
                  worst_row_activity_residual);
    max_relative_activity_residual =
        std::max(max_relative_activity_residual,
                 activity_residual / (1.0 + active_bound_norm));
  }
  result.metrics.max_primal_violation =
      std::max({max_col_primal_violation, max_row_primal_violation,
                result.metrics.max_row_activity_residual});
  result.worst_primal_row =
      max_row_primal_violation >= result.metrics.max_row_activity_residual
          ? worst_row_bound_violation
          : worst_row_activity_residual;
  result.metrics.relative_primal_residual =
      std::max(max_relative_primal_violation, max_relative_activity_residual);

  std::vector<double> reduced_cost(lp.num_col_);
  double near_zero_cost_norm = 0.0;
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    reduced_cost[col] = lp.col_cost_[col] - transpose_product[col];
    if (std::fabs(reduced_cost[col]) < near_tolerance) {
      near_zero_cost_norm =
          std::max(near_zero_cost_norm, std::fabs(lp.col_cost_[col]));
    }
  }
  double max_col_dual_violation = 0.0;
  double max_row_dual_violation = 0.0;
  double max_col_complementarity = 0.0;
  double max_row_complementarity = 0.0;
  HighsInt worst_col_sign_violation = -1;
  HighsInt worst_stationarity_col = -1;
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    const double adjusted_dual = objective_sense * reduced_cost[col];
    const double sign_violation = dualSignViolation(
        lp.col_lower_[col], lp.col_upper_[col], solution.col_value[col],
        adjusted_dual, near_tolerance);
    updateMaximum(sign_violation, col, max_col_dual_violation,
                  worst_col_sign_violation);
    double relative_cost_measure = near_zero_cost_norm;
    if (sign_violation > 0.0 &&
        std::fabs(reduced_cost[col]) >= near_tolerance) {
      relative_cost_measure =
          std::max(relative_cost_measure, std::fabs(lp.col_cost_[col]));
    }
    result.metrics.relative_dual_residual =
        std::max(result.metrics.relative_dual_residual,
                 sign_violation / (1.0 + relative_cost_measure));

    const double stationarity_residual =
        std::fabs(reduced_cost[col] - solution.col_dual[col]);
    if (stationarity_residual > result.metrics.max_stationarity_residual) {
      result.metrics.max_stationarity_residual = stationarity_residual;
      worst_stationarity_col = col;
    }
    result.metrics.relative_dual_residual =
        std::max(result.metrics.relative_dual_residual,
                 stationarity_residual / (1.0 + near_zero_cost_norm));

    const double complementarity =
        complementarityResidual(lp.col_lower_[col], lp.col_upper_[col],
                                solution.col_value[col], adjusted_dual);
    result.metrics.absolute_complementarity += complementarity;
    updateMaximum(complementarity, col, max_col_complementarity,
                  result.worst_complementarity_col);
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    const double adjusted_dual = objective_sense * solution.row_dual[row];
    const double sign_violation =
        dualSignViolation(lp.row_lower_[row], lp.row_upper_[row],
                          row_activity[row], adjusted_dual, near_tolerance);
    updateMaximum(sign_violation, row, max_row_dual_violation,
                  result.worst_dual_row);
    result.metrics.relative_dual_residual =
        std::max(result.metrics.relative_dual_residual,
                 sign_violation / (1.0 + near_zero_cost_norm));

    const double complementarity =
        complementarityResidual(lp.row_lower_[row], lp.row_upper_[row],
                                row_activity[row], adjusted_dual);
    result.metrics.absolute_complementarity += complementarity;
    updateMaximum(complementarity, row, max_row_complementarity,
                  result.worst_complementarity_row);
  }
  if (max_col_complementarity >= max_row_complementarity) {
    result.metrics.max_complementarity = max_col_complementarity;
    result.worst_complementarity_row = -1;
  } else {
    result.metrics.max_complementarity = max_row_complementarity;
    result.worst_complementarity_col = -1;
  }
  result.metrics.max_dual_violation =
      std::max({max_col_dual_violation, max_row_dual_violation,
                result.metrics.max_stationarity_residual});
  result.worst_dual_col =
      max_col_dual_violation >= result.metrics.max_stationarity_residual
          ? worst_col_sign_violation
          : worst_stationarity_col;

  ObjectiveSum primal_sum;
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    if (!primal_sum.add(lp.col_cost_[col] * solution.col_value[col])) {
      result.error = "Primal objective calculation produced a non-finite value";
      return result;
    }
  }

  bool dual_objective_valid = true;
  ObjectiveSum dual_sum;
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    dual_objective_valid &= addDualObjectiveContribution(
        lp.col_lower_[col], lp.col_upper_[col], reduced_cost[col],
        objective_sense * reduced_cost[col], dual_sum);
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    dual_objective_valid &= addDualObjectiveContribution(
        lp.row_lower_[row], lp.row_upper_[row], solution.row_dual[row],
        objective_sense * solution.row_dual[row],
        dual_sum);
  }
  if (!dual_objective_valid) {
    result.error = "Dual objective calculation produced a non-finite value";
    return result;
  }

  const double primal_objective = primal_sum.value();
  const double dual_objective = dual_sum.value();
  if (!std::isfinite(primal_objective) || !std::isfinite(dual_objective)) {
    result.error = "Objective calculation exceeded FP64 range";
    return result;
  }

    result.metrics.primal_objective = lp.offset_ + primal_objective;
    result.metrics.dual_objective = lp.offset_ + dual_objective;
    result.metrics.absolute_objective_gap =
      std::fabs(primal_objective - dual_objective);
    const double objective_scale =
      std::max({1.0, std::fabs(primal_objective), std::fabs(dual_objective)});
    const double objective_denominator =
      1.0 / objective_scale + std::fabs(primal_objective) / objective_scale +
      std::fabs(dual_objective) / objective_scale;
  result.metrics.relative_objective_gap =
      (result.metrics.absolute_objective_gap / objective_scale) /
      objective_denominator;
    result.metrics.relative_complementarity =
      (result.metrics.absolute_complementarity / objective_scale) /
      objective_denominator;
  result.metrics_finite = finiteMetrics(result.metrics);
  if (!result.metrics_finite) {
    result.error = "KKT verification produced a non-finite metric";
    return result;
  }

  result.primal_feasible = result.metrics.relative_primal_residual <= tolerance;
  result.dual_feasible = result.metrics.relative_dual_residual <= tolerance;
  result.objective_gap_ok = result.metrics.relative_objective_gap <= tolerance;
    result.complementarity_ok = result.metrics.relative_complementarity <= tolerance;
  result.pass =
      result.primal_feasible && result.dual_feasible && result.objective_gap_ok &&
      result.complementarity_ok;
  return result;
}

}  // namespace verification
}  // namespace hipdlp