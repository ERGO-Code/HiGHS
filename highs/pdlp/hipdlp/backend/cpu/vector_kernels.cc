#include "pdlp/hipdlp/backend/cpu/vector_kernels.hpp"

#include <algorithm>
#include <cmath>

namespace hipdlp {
namespace backend {
namespace cpu {

namespace {

void updateScaledSumOfSquares(const double value, double& scale,
                              double& sum_of_squares) {
  const double magnitude = std::abs(value);
  if (magnitude == 0.0) return;
  if (scale < magnitude) {
    const double ratio = scale / magnitude;
    sum_of_squares = 1.0 + sum_of_squares * ratio * ratio;
    scale = magnitude;
  } else {
    const double ratio = magnitude / scale;
    sum_of_squares += ratio * ratio;
  }
}

}  // namespace

void fillVector(std::vector<double>& vector, const double value) {
  std::fill(vector.begin(), vector.end(), value);
}

void copyVector(const std::vector<double>& source,
                std::vector<double>& destination) {
  for (std::size_t index = 0; index < source.size(); ++index)
    destination[index] = source[index];
}

void scaleVector(const double factor, std::vector<double>& vector) {
  for (std::size_t index = 0; index < vector.size(); ++index)
    vector[index] *= factor;
}

void axpyVector(const double factor, const std::vector<double>& input,
                std::vector<double>& output) {
  for (std::size_t index = 0; index < input.size(); ++index)
    output[index] += factor * input[index];
}

double dotVector(const std::vector<double>& left,
                 const std::vector<double>& right) {
  double result = 0.0;
  for (std::size_t index = 0; index < left.size(); ++index)
    result += left[index] * right[index];
  return result;
}

double norm2Vector(const std::vector<double>& vector) {
  double scale = 0.0;
  double sum_of_squares = 1.0;
  for (std::size_t index = 0; index < vector.size(); ++index)
    updateScaledSumOfSquares(vector[index], scale, sum_of_squares);
  return scale == 0.0 ? 0.0 : scale * std::sqrt(sum_of_squares);
}

double differenceNorm2Vector(const std::vector<double>& left,
                             const std::vector<double>& right) {
  double scale = 0.0;
  double sum_of_squares = 1.0;
  for (std::size_t index = 0; index < left.size(); ++index) {
    const double difference = left[index] - right[index];
    updateScaledSumOfSquares(difference, scale, sum_of_squares);
  }
  return scale == 0.0 ? 0.0 : scale * std::sqrt(sum_of_squares);
}

void projectedPrimalUpdate(const std::vector<double>& current,
                           const std::vector<double>& transpose_product,
                           const std::vector<double>& costs,
                           const std::vector<double>& lower,
                           const std::vector<double>& upper,
                           const double step_size,
                           std::vector<double>& output) {
  for (std::size_t col = 0; col < current.size(); ++col) {
    const double candidate =
        current[col] - step_size * (costs[col] - transpose_product[col]);
    output[col] = std::max(lower[col], std::min(candidate, upper[col]));
  }
}

void projectedDualUpdate(const std::vector<double>& current,
                         const std::vector<double>& extrapolated_product,
                         const std::vector<double>& row_lower,
                         const std::vector<HighsBool>& equality_rows,
                         const double step_size, std::vector<double>& output) {
  for (std::size_t row = 0; row < current.size(); ++row) {
    const double candidate =
        current[row] + step_size * (row_lower[row] - extrapolated_product[row]);
    output[row] = equality_rows[row] ? candidate : std::max(0.0, candidate);
  }
}

void blendVector(const std::vector<double>& left,
                 const std::vector<double>& right, const double left_weight,
                 std::vector<double>& output) {
  const double right_weight = 1.0 - left_weight;
  for (std::size_t index = 0; index < left.size(); ++index) {
    output[index] = left_weight * left[index] + right_weight * right[index];
  }
}

void averageUpdate(const std::vector<double>& current, const double weight,
                   std::vector<double>& sum) {
  for (std::size_t index = 0; index < current.size(); ++index)
    sum[index] += weight * current[index];
}

}  // namespace cpu
}  // namespace backend
}  // namespace hipdlp