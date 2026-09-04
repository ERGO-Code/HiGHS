#include "pdlp/hipdlp/verification/certificate_verifier.hpp"

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
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index])) return false;
  }
  return true;
}

bool validBounds(const std::vector<double>& lower,
                 const std::vector<double>& upper) {
  if (lower.size() != upper.size()) return false;
  for (std::size_t index = 0; index < lower.size(); ++index) {
    if (std::isnan(lower[index]) || std::isnan(upper[index]) ||
        lower[index] >= kHighsInf || upper[index] <= -kHighsInf ||
        lower[index] > upper[index]) {
      return false;
    }
  }
  return true;
}

bool validSparseStorage(const HighsSparseMatrix& matrix) {
  if (matrix.num_col_ < 0 || matrix.num_row_ < 0) return false;
  if (matrix.format_ != MatrixFormat::kColwise &&
      matrix.format_ != MatrixFormat::kRowwise &&
      matrix.format_ != MatrixFormat::kRowwisePartitioned) {
    return false;
  }
  const bool colwise = matrix.format_ == MatrixFormat::kColwise;
  const HighsInt num_vector = colwise ? matrix.num_col_ : matrix.num_row_;
  const HighsInt index_limit = colwise ? matrix.num_row_ : matrix.num_col_;
  if (matrix.start_.size() < static_cast<std::size_t>(num_vector) + 1 ||
      matrix.start_[0] != 0) {
    return false;
  }
  for (HighsInt vector = 0; vector < num_vector; ++vector) {
    if (matrix.start_[vector] < 0 ||
        matrix.start_[vector] > matrix.start_[vector + 1]) {
      return false;
    }
  }
  const HighsInt num_nonzero = matrix.start_[num_vector];
  if (num_nonzero < 0 ||
      matrix.index_.size() < static_cast<std::size_t>(num_nonzero) ||
      matrix.value_.size() < static_cast<std::size_t>(num_nonzero)) {
    return false;
  }
  for (HighsInt nonzero = 0; nonzero < num_nonzero; ++nonzero) {
    if (matrix.index_[nonzero] < 0 || matrix.index_[nonzero] >= index_limit)
      return false;
    if (!std::isfinite(matrix.value_[nonzero])) return false;
  }
  return true;
}

bool validLp(const HighsLp& lp) {
  if (lp.num_col_ < 0 || lp.num_row_ < 0) return false;
  if (lp.sense_ != ObjSense::kMinimize && lp.sense_ != ObjSense::kMaximize)
    return false;
  if (lp.col_cost_.size() != static_cast<std::size_t>(lp.num_col_) ||
      lp.col_lower_.size() != static_cast<std::size_t>(lp.num_col_) ||
      lp.col_upper_.size() != static_cast<std::size_t>(lp.num_col_) ||
      lp.row_lower_.size() != static_cast<std::size_t>(lp.num_row_) ||
      lp.row_upper_.size() != static_cast<std::size_t>(lp.num_row_)) {
    return false;
  }
  if (!finiteVector(lp.col_cost_)) return false;
  if (!validBounds(lp.col_lower_, lp.col_upper_)) return false;
  if (!validBounds(lp.row_lower_, lp.row_upper_)) return false;
  if (lp.a_matrix_.num_col_ != lp.num_col_ ||
      lp.a_matrix_.num_row_ != lp.num_row_) {
    return false;
  }
  return validSparseStorage(lp.a_matrix_);
}

double euclideanNorm(const std::vector<double>& values) {
  double norm = 0.0;
  for (std::size_t index = 0; index < values.size(); ++index)
    norm = std::hypot(norm, values[index]);
  return norm;
}

void updateWorst(const double violation, const HighsInt index, double& worst,
                 HighsInt& worst_index) {
  if (violation > worst) {
    worst = violation;
    worst_index = index;
  }
}

// Distance from value to the recession cone induced by [lower, upper].
double recessionViolation(const double value, const double lower,
                          const double upper) {
  const bool bounded_below = hasLower(lower);
  const bool bounded_above = hasUpper(upper);
  if (bounded_below && bounded_above) return std::fabs(value);
  if (bounded_below) return value < 0.0 ? -value : 0.0;
  if (bounded_above) return value > 0.0 ? value : 0.0;
  return 0.0;
}

// A multiplier may only be positive where the lower bound is finite and
// negative where the upper bound is finite.
double multiplierSignViolation(const double value, const double lower,
                               const double upper) {
  const bool bounded_below = hasLower(lower);
  const bool bounded_above = hasUpper(upper);
  if (bounded_below && bounded_above) return 0.0;
  if (bounded_below) return value < 0.0 ? -value : 0.0;
  if (bounded_above) return value > 0.0 ? value : 0.0;
  return std::fabs(value);
}

std::vector<double> rowActivity(const HighsLp& lp,
                                const std::vector<double>& column_values) {
  std::vector<double> activity(static_cast<std::size_t>(lp.num_row_), 0.0);
  const HighsSparseMatrix& matrix = lp.a_matrix_;
  if (matrix.format_ == MatrixFormat::kColwise) {
    for (HighsInt col = 0; col < lp.num_col_; ++col) {
      for (HighsInt entry = matrix.start_[col]; entry < matrix.start_[col + 1];
           ++entry) {
        activity[matrix.index_[entry]] +=
            matrix.value_[entry] * column_values[col];
      }
    }
  } else {
    for (HighsInt row = 0; row < lp.num_row_; ++row) {
      double sum = 0.0;
      for (HighsInt entry = matrix.start_[row]; entry < matrix.start_[row + 1];
           ++entry) {
        sum += matrix.value_[entry] * column_values[matrix.index_[entry]];
      }
      activity[row] = sum;
    }
  }
  return activity;
}

std::vector<double> transposeActivity(const HighsLp& lp,
                                      const std::vector<double>& row_values) {
  std::vector<double> activity(static_cast<std::size_t>(lp.num_col_), 0.0);
  const HighsSparseMatrix& matrix = lp.a_matrix_;
  if (matrix.format_ == MatrixFormat::kColwise) {
    for (HighsInt col = 0; col < lp.num_col_; ++col) {
      double sum = 0.0;
      for (HighsInt entry = matrix.start_[col]; entry < matrix.start_[col + 1];
           ++entry) {
        sum += matrix.value_[entry] * row_values[matrix.index_[entry]];
      }
      activity[col] = sum;
    }
  } else {
    for (HighsInt row = 0; row < lp.num_row_; ++row) {
      for (HighsInt entry = matrix.start_[row]; entry < matrix.start_[row + 1];
           ++entry) {
        activity[matrix.index_[entry]] +=
            matrix.value_[entry] * row_values[row];
      }
    }
  }
  return activity;
}

bool addBoundSupport(const double multiplier, const double lower,
                     const double upper, double& support) {
  double contribution = 0.0;
  if (multiplier > 0.0) {
    if (!hasLower(lower)) return false;
    contribution = multiplier * lower;
  } else if (multiplier < 0.0) {
    if (!hasUpper(upper)) return false;
    contribution = multiplier * upper;
  }
  if (!std::isfinite(contribution) ||
      !std::isfinite(support + contribution))
    return false;
  support += contribution;
  return true;
}

}  // namespace

PrimalRayVerificationResult verifyPrimalRay(
    const HighsLp& lp, const std::vector<double>& column_direction,
    const double tolerance) {
  PrimalRayVerificationResult result;
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    result.error = "certificate tolerance must be finite and positive";
    return result;
  }
  if (!validLp(lp)) {
    result.error = "primal ray verification received an invalid LP";
    return result;
  }
  if (column_direction.size() != static_cast<std::size_t>(lp.num_col_)) {
    result.error = "primal ray has the wrong column dimension";
    return result;
  }
  if (!finiteVector(column_direction)) {
    result.error = "primal ray contains a non-finite value";
    return result;
  }
  result.input_valid = true;

  const double direction_norm = euclideanNorm(column_direction);
  result.metrics.direction_norm = direction_norm;
  if (!std::isfinite(direction_norm)) {
    result.error = "primal ray norm is non-finite";
    return result;
  }
  if (direction_norm <= 0.0) {
    result.metrics_finite = true;
    result.error = "primal ray direction is zero";
    return result;
  }
  result.direction_nonzero = true;

  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    const double violation = recessionViolation(
        column_direction[col], lp.col_lower_[col], lp.col_upper_[col]);
    updateWorst(violation, col, result.metrics.max_raw_column_violation,
                result.metrics.worst_col);
  }

  const std::vector<double> activity = rowActivity(lp, column_direction);
  if (!finiteVector(activity)) {
    result.error = "primal ray row activity is non-finite";
    return result;
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    const double violation = recessionViolation(
        activity[row], lp.row_lower_[row], lp.row_upper_[row]);
    updateWorst(violation, row, result.metrics.max_raw_row_violation,
                result.metrics.worst_row);
  }

  double objective_direction = 0.0;
  for (HighsInt col = 0; col < lp.num_col_; ++col)
    objective_direction += lp.col_cost_[col] * column_direction[col];
  result.metrics.raw_objective_direction = objective_direction;

  const double sense_sign = lp.sense_ == ObjSense::kMinimize ? 1.0 : -1.0;
  result.metrics.max_normalized_column_violation =
      result.metrics.max_raw_column_violation / direction_norm;
  result.metrics.max_normalized_row_violation =
      result.metrics.max_raw_row_violation / direction_norm;
  result.metrics.normalized_objective_direction =
      objective_direction / direction_norm;
  result.metrics.objective_margin =
      -sense_sign * result.metrics.normalized_objective_direction;

  result.metrics_finite =
      std::isfinite(result.metrics.max_normalized_column_violation) &&
      std::isfinite(result.metrics.max_normalized_row_violation) &&
      std::isfinite(result.metrics.raw_objective_direction) &&
      std::isfinite(result.metrics.objective_margin);
  if (!result.metrics_finite) {
    result.error = "primal ray metrics are non-finite";
    return result;
  }

  result.column_cone_ok =
      result.metrics.max_normalized_column_violation <= tolerance;
  result.row_cone_ok =
      result.metrics.max_normalized_row_violation <= tolerance;
  result.objective_improving = result.metrics.objective_margin > tolerance;
  result.pass = result.column_cone_ok && result.row_cone_ok &&
                result.objective_improving;
  return result;
}

FarkasVerificationResult verifyFarkasCertificate(
    const HighsLp& lp, const std::vector<double>& row_direction,
    const std::vector<double>& column_direction, const double tolerance) {
  FarkasVerificationResult result;
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    result.error = "certificate tolerance must be finite and positive";
    return result;
  }
  if (!validLp(lp)) {
    result.error = "Farkas verification received an invalid LP";
    return result;
  }
  if (row_direction.size() != static_cast<std::size_t>(lp.num_row_)) {
    result.error = "Farkas certificate has the wrong row dimension";
    return result;
  }
  if (column_direction.size() != static_cast<std::size_t>(lp.num_col_)) {
    result.error = "Farkas certificate has the wrong column dimension";
    return result;
  }
  if (!finiteVector(row_direction) || !finiteVector(column_direction)) {
    result.error = "Farkas certificate contains a non-finite value";
    return result;
  }
  result.input_valid = true;

  const double direction_norm =
      std::hypot(euclideanNorm(row_direction), euclideanNorm(column_direction));
  result.metrics.direction_norm = direction_norm;
  if (!std::isfinite(direction_norm)) {
    result.error = "Farkas certificate norm is non-finite";
    return result;
  }
  if (direction_norm <= 0.0) {
    result.metrics_finite = true;
    result.error = "Farkas certificate direction is zero";
    return result;
  }
  result.direction_nonzero = true;

  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    const double violation = multiplierSignViolation(
        row_direction[row], lp.row_lower_[row], lp.row_upper_[row]);
    updateWorst(violation, row, result.metrics.max_raw_row_sign_violation,
                result.metrics.worst_row);
  }

  // The verifier rebuilds A^T y instead of trusting the supplied direction.
  const std::vector<double> transpose_product =
      transposeActivity(lp, row_direction);
  if (!finiteVector(transpose_product)) {
    result.error = "Farkas certificate transpose product is non-finite";
    return result;
  }
  double worst_column_violation = 0.0;
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    const double sign_violation = multiplierSignViolation(
        column_direction[col], lp.col_lower_[col], lp.col_upper_[col]);
    if (sign_violation > result.metrics.max_raw_column_sign_violation)
      result.metrics.max_raw_column_sign_violation = sign_violation;
    const double residual =
        std::fabs(column_direction[col] + transpose_product[col]);
    if (residual > result.metrics.max_raw_stationarity_residual)
      result.metrics.max_raw_stationarity_residual = residual;
    const double combined =
        sign_violation > residual ? sign_violation : residual;
    updateWorst(combined, col, worst_column_violation, result.metrics.worst_col);
  }

  double bound_objective = 0.0;
  bool bound_objective_finite = true;
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    bound_objective_finite &=
        addBoundSupport(row_direction[row], lp.row_lower_[row],
                        lp.row_upper_[row], bound_objective);
  }
  for (HighsInt col = 0; col < lp.num_col_; ++col) {
    const double exact_column_multiplier = -transpose_product[col];
    bound_objective_finite &=
        addBoundSupport(exact_column_multiplier, lp.col_lower_[col],
                        lp.col_upper_[col], bound_objective);
  }
  result.metrics.raw_bound_objective = bound_objective;

  result.metrics.max_normalized_row_sign_violation =
      result.metrics.max_raw_row_sign_violation / direction_norm;
  result.metrics.max_normalized_column_sign_violation =
      result.metrics.max_raw_column_sign_violation / direction_norm;
  result.metrics.max_normalized_stationarity_residual =
      result.metrics.max_raw_stationarity_residual / direction_norm;
    result.metrics.certificate_margin =
      bound_objective_finite ? bound_objective / direction_norm : 0.0;
    result.bound_objective_finite = bound_objective_finite;

  result.metrics_finite =
      std::isfinite(result.metrics.max_normalized_row_sign_violation) &&
      std::isfinite(result.metrics.max_normalized_column_sign_violation) &&
      std::isfinite(result.metrics.max_normalized_stationarity_residual) &&
      std::isfinite(result.metrics.raw_bound_objective) &&
      std::isfinite(result.metrics.certificate_margin);
  if (!result.metrics_finite) {
    result.error = "Farkas certificate metrics are non-finite";
    return result;
  }

  result.row_signs_ok =
      result.metrics.max_normalized_row_sign_violation <= tolerance;
  result.column_signs_ok =
      result.metrics.max_normalized_column_sign_violation <= tolerance;
  result.stationarity_ok =
      result.metrics.max_normalized_stationarity_residual <= tolerance;
  result.bound_objective_positive =
      result.bound_objective_finite &&
      result.metrics.certificate_margin > tolerance;
  result.pass = result.row_signs_ok && result.column_signs_ok &&
                result.stationarity_ok && result.bound_objective_positive;
  return result;
}

}  // namespace verification
}  // namespace hipdlp
