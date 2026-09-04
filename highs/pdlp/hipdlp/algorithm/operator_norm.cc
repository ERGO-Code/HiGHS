#include "pdlp/hipdlp/algorithm/operator_norm.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include "lp_data/HConst.h"
#include "pdlp/hipdlp/backend/backend.hpp"
#include "pdlp/hipdlp/model/problem.hpp"

namespace hipdlp {
namespace algorithm {
namespace {

const std::size_t kScalarReadbacksPerIteration = 2;

bool finitePositive(const double value) {
  return std::isfinite(value) && value > 0.0;
}

bool finiteNonnegative(const double value) {
  return std::isfinite(value) && value >= 0.0;
}

OperatorNormResult failure(const api::Error& error,
                           const OperatorNormDiagnostics& diagnostics) {
  OperatorNormResult result;
  result.error = error;
  result.diagnostics = diagnostics;
  return result;
}

api::Error createVector(backend::Backend& backend, const HighsInt size,
                        backend::VectorHandle& vector,
                        OperatorNormDiagnostics& diagnostics) {
  ++diagnostics.resource_creation_calls_before_iterations;
  return backend.createVector(size, vector);
}

api::Error norm2(backend::Backend& backend, const backend::VectorHandle& vector,
                 double& result, OperatorNormDiagnostics& diagnostics) {
  ++diagnostics.device_to_host_scalar_readbacks_during_iterations;
  return backend.norm2(vector, result);
}

api::Error dot(backend::Backend& backend, const backend::VectorHandle& left,
               const backend::VectorHandle& right, double& result,
               OperatorNormDiagnostics& diagnostics) {
  ++diagnostics.device_to_host_scalar_readbacks_during_iterations;
  return backend.dot(left, right, result);
}

}  // namespace

OperatorNormParameters::OperatorNormParameters()
    : iteration_budget(20), relative_tolerance(1e-6) {}

OperatorNormDiagnostics::OperatorNormDiagnostics()
    : resource_creation_calls_before_iterations(0),
      resource_creation_calls_during_iterations(0),
      device_to_host_scalar_readbacks_per_iteration(0),
      device_to_host_scalar_readbacks_during_iterations(0) {}

OperatorNormResult::OperatorNormResult()
    : operator_norm_squared(0.0),
      performed_iterations(0),
      final_relative_change(0.0),
      tolerance_met(false),
      empty_matrix(false) {}

bool OperatorNormResult::ok() const {
  if (!error.ok()) return false;
  if (empty_matrix)
    return operator_norm_squared == 0.0 && performed_iterations == 0 &&
           final_relative_change == 0.0 && !tolerance_met;
  return finitePositive(operator_norm_squared) && performed_iterations > 0 &&
         finiteNonnegative(final_relative_change);
}

api::Error computeOperatorNormSquaredUpperBound(
    const model::SparseMatrix& matrix, double& upper_bound_squared) {
  upper_bound_squared = 0.0;
  if (matrix.num_rows < 0 || matrix.num_cols < 0 || matrix.starts.empty() ||
      matrix.starts.size() != static_cast<std::size_t>(matrix.num_rows) + 1 ||
      matrix.starts[0] != 0) {
    return api::Error(api::ErrorCode::kInvalidArgument,
                      "HiPDLP operator-norm bound matrix is malformed");
  }
  std::vector<long double> column_sums(matrix.num_cols, 0.0L);
  long double maximum_row_sum = 0.0L;
  for (HighsInt row = 0; row < matrix.num_rows; ++row) {
    if (matrix.starts[row] < 0 ||
        matrix.starts[row] > matrix.starts[row + 1]) {
      return api::Error(api::ErrorCode::kInvalidArgument,
                        "HiPDLP operator-norm bound starts are malformed");
    }
    long double row_sum = 0.0L;
    for (HighsInt entry = matrix.starts[row];
         entry < matrix.starts[row + 1]; ++entry) {
      if (entry < 0 ||
          static_cast<std::size_t>(entry) >= matrix.indices.size() ||
          static_cast<std::size_t>(entry) >= matrix.values.size()) {
        return api::Error(api::ErrorCode::kInvalidArgument,
                          "HiPDLP operator-norm bound entries are malformed");
      }
      const HighsInt col = matrix.indices[entry];
      const double value = matrix.values[entry];
      if (col < 0 || col >= matrix.num_cols || !std::isfinite(value)) {
        return api::Error(api::ErrorCode::kInvalidArgument,
                          "HiPDLP operator-norm bound entry is invalid");
      }
      const long double magnitude = std::fabs(static_cast<long double>(value));
      row_sum += magnitude;
      column_sums[col] += magnitude;
      if (!std::isfinite(row_sum) || !std::isfinite(column_sums[col])) {
        return api::Error(api::ErrorCode::kInvalidArgument,
                          "HiPDLP operator-norm bound sum overflowed");
      }
    }
    if (row_sum > maximum_row_sum) maximum_row_sum = row_sum;
  }
  if (matrix.starts[matrix.num_rows] < 0 ||
      static_cast<std::size_t>(matrix.starts[matrix.num_rows]) !=
          matrix.indices.size() ||
      matrix.indices.size() != matrix.values.size()) {
    return api::Error(api::ErrorCode::kInvalidArgument,
                      "HiPDLP operator-norm bound arrays are malformed");
  }
  long double maximum_column_sum = 0.0L;
  for (HighsInt col = 0; col < matrix.num_cols; ++col) {
    if (column_sums[col] > maximum_column_sum)
      maximum_column_sum = column_sums[col];
  }
  const long double bound = maximum_row_sum * maximum_column_sum;
  if (!std::isfinite(bound) ||
      bound > static_cast<long double>(std::numeric_limits<double>::max())) {
    return api::Error(api::ErrorCode::kInvalidArgument,
                      "HiPDLP operator-norm bound is outside FP64 range");
  }
  if (bound == 0.0L) return api::Error();
  upper_bound_squared = std::nextafter(
      static_cast<double>(bound), std::numeric_limits<double>::infinity());
  return api::Error();
}

OperatorNormResult estimateOperatorNormSquared(
    backend::Backend& backend, const backend::ProblemHandle& problem,
    const HighsInt num_rows, const HighsInt num_cols,
    const OperatorNormParameters& parameters) {
  OperatorNormDiagnostics diagnostics;
  if (parameters.iteration_budget <= 0 ||
      parameters.iteration_budget == kHighsIInf ||
      !finitePositive(parameters.relative_tolerance)) {
    return failure(api::Error(api::ErrorCode::kInvalidOption,
                              "HiPDLP operator-norm iteration budget and "
                              "relative tolerance must be finite and positive"),
                   diagnostics);
  }
  if (num_rows < 0 || num_cols < 0 || !problem.valid()) {
    return failure(api::Error(api::ErrorCode::kInvalidArgument,
                              "HiPDLP operator-norm problem is invalid"),
                   diagnostics);
  }
  if (num_rows == 0 || num_cols == 0) {
    OperatorNormResult result;
    result.empty_matrix = true;
    return result;
  }

  backend::VectorHandle row_vector;
  backend::VectorHandle transpose_product;
  backend::VectorHandle next_row_vector;
  api::Error error = createVector(backend, num_rows, row_vector, diagnostics);
  if (!error.ok()) return failure(error, diagnostics);
  error = createVector(backend, num_cols, transpose_product, diagnostics);
  if (!error.ok()) return failure(error, diagnostics);
  error = createVector(backend, num_rows, next_row_vector, diagnostics);
  if (!error.ok()) return failure(error, diagnostics);
  diagnostics.device_to_host_scalar_readbacks_per_iteration =
      kScalarReadbacksPerIteration;

  error = backend.fill(row_vector, 1.0);
  if (!error.ok()) return failure(error, diagnostics);

  double previous_operator_norm_squared = 0.0;
  double operator_norm_squared = 0.0;
  double final_relative_change = 0.0;
  HighsInt performed_iterations = 0;
  bool tolerance_met = false;
  for (HighsInt iteration = 0; iteration < parameters.iteration_budget;
       ++iteration) {
    error = backend.multiplyTranspose(problem, row_vector, transpose_product);
    if (!error.ok()) return failure(error, diagnostics);
    error = backend.multiply(problem, transpose_product, next_row_vector);
    if (!error.ok()) return failure(error, diagnostics);

    double next_row_norm = 0.0;
    error = norm2(backend, next_row_vector, next_row_norm, diagnostics);
    if (!error.ok()) return failure(error, diagnostics);
    if (!finitePositive(next_row_norm)) {
      return failure(api::Error(api::ErrorCode::kInvalidArgument,
                                "HiPDLP operator-norm AA' iterate norm is not "
                                "finite and positive"),
                     diagnostics);
    }
    const double inverse_next_row_norm = 1.0 / next_row_norm;
    if (!finitePositive(inverse_next_row_norm)) {
      return failure(api::Error(api::ErrorCode::kInvalidArgument,
                                "HiPDLP operator-norm normalization is "
                                "outside FP64 range"),
                     diagnostics);
    }
    error = backend.scale(inverse_next_row_norm, next_row_vector);
    if (!error.ok()) return failure(error, diagnostics);
    error = backend.copy(next_row_vector, row_vector);
    if (!error.ok()) return failure(error, diagnostics);

    error = backend.multiplyTranspose(problem, row_vector, transpose_product);
    if (!error.ok()) return failure(error, diagnostics);
    double next_operator_norm_squared = 0.0;
    error = dot(backend, transpose_product, transpose_product,
                next_operator_norm_squared, diagnostics);
    if (!error.ok()) return failure(error, diagnostics);
    if (!finitePositive(next_operator_norm_squared)) {
      return failure(api::Error(api::ErrorCode::kInvalidArgument,
                                "HiPDLP operator-norm Rayleigh quotient is "
                                "not finite and positive"),
                     diagnostics);
    }

    const double absolute_change =
        std::fabs(next_operator_norm_squared - previous_operator_norm_squared);
    const double relative_change = absolute_change / next_operator_norm_squared;
    if (!finiteNonnegative(absolute_change) ||
        !finiteNonnegative(relative_change)) {
      return failure(api::Error(api::ErrorCode::kInvalidArgument,
                                "HiPDLP operator-norm relative change is "
                                "outside FP64 range"),
                     diagnostics);
    }

    operator_norm_squared = next_operator_norm_squared;
    previous_operator_norm_squared = next_operator_norm_squared;
    final_relative_change = relative_change;
    performed_iterations = iteration + 1;
    tolerance_met = relative_change <= parameters.relative_tolerance;
    if (tolerance_met) break;
  }

  OperatorNormResult result;
  result.operator_norm_squared = operator_norm_squared;
  result.performed_iterations = performed_iterations;
  result.final_relative_change = final_relative_change;
  result.tolerance_met = tolerance_met;
  result.diagnostics = diagnostics;
  return result;
}

}  // namespace algorithm
}  // namespace hipdlp