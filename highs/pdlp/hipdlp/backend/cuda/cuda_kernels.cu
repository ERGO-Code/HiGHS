#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>

#include "pdlp/hipdlp/backend/cuda/cuda_kernels.cuh"

namespace hipdlp {
namespace backend {
namespace cuda {

namespace {

const int kBlockSize = 256;
const int kMaximumBlocks = 65535;
const int kWarpSize = 32;
using GridIndex = std::uint64_t;
static_assert(std::numeric_limits<GridIndex>::max() -
          static_cast<std::uint64_t>(std::numeric_limits<HighsInt>::max()) >=
        static_cast<std::uint64_t>(kMaximumBlocks) * kBlockSize,
        "grid stride must not overflow at the vector size limit");

int blockCount(const HighsInt size) {
  if (size <= 0) return 0;
  const HighsInt count = 1 + (size - 1) / kBlockSize;
  return static_cast<int>(count < kMaximumBlocks ? count : kMaximumBlocks);
}

__device__ void recordNonFinite(const double value, int* error_flag) {
  if (!isfinite(value)) atomicExch(error_flag, 1);
}

__global__ void fillKernel(const HighsInt size, const double value,
                           double* output) {
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride)
    output[index] = value;
}

__global__ void scaleKernel(const HighsInt size, const double factor,
                            double* values, int* error_flag) {
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    const double value = values[index] * factor;
    values[index] = value;
    recordNonFinite(value, error_flag);
  }
}

__global__ void axpyKernel(const HighsInt size, const double factor,
                           const double* input, double* output,
                           int* error_flag) {
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    const double value = output[index] + factor * input[index];
    output[index] = value;
    recordNonFinite(value, error_flag);
  }
}

__global__ void differenceKernel(const HighsInt size, const double* left,
                                 const double* right, double* output,
                                 int* error_flag) {
  const GridIndex stride = static_cast<GridIndex>(blockDim.x) * gridDim.x;
  for (GridIndex index = static_cast<GridIndex>(blockIdx.x) * blockDim.x + threadIdx.x;
    index < static_cast<GridIndex>(size); index += stride) {
    const double value = __dsub_rn(left[index], right[index]);
    output[index] = value;
    recordNonFinite(value, error_flag);
  }
}

template <bool Reflect>
__global__ void projectedPrimalKernel(
    const HighsInt size, const double* current, const double* transpose_product,
    const double* costs, const double* lower, const double* upper,
    const double literal_step, double* output, double* reflected, int* error_flag,
    const double* step_parameter) {
  const double step_size = step_parameter == nullptr ? literal_step : *step_parameter;
  if (step_parameter != nullptr && (!isfinite(step_size) || step_size < 0.0)) {
    atomicExch(error_flag, 1);
    return;
  }
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    double value =
        current[index] - step_size * (costs[index] - transpose_product[index]);
    if (value < lower[index]) value = lower[index];
    if (value > upper[index]) value = upper[index];
    output[index] = value;
    recordNonFinite(value, error_flag);
    if (Reflect) {
      const double doubled = __dmul_rn(2.0, value);
      const double reflection = __dsub_rn(doubled, current[index]);
      reflected[index] = reflection;
      recordNonFinite(doubled, error_flag);
      recordNonFinite(reflection, error_flag);
    }
  }
}

template <bool Reflect>
__global__ void projectedDualKernel(const HighsInt size, const double* current,
                                    const double* extrapolated_product,
                                    const double* row_lower,
                                    const HighsBool* equality_rows,
                                    const double literal_step, double* output,
                                    double* reflected, int* error_flag,
                                    const double* step_parameter) {
  const double step_size = step_parameter == nullptr ? literal_step : *step_parameter;
  if (step_parameter != nullptr && (!isfinite(step_size) || step_size < 0.0)) {
    atomicExch(error_flag, 1);
    return;
  }
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
     recordNonFinite(extrapolated_product[index], error_flag);
    double value = current[index] +
                   step_size * (row_lower[index] - extrapolated_product[index]);
    if (!equality_rows[index] && value < 0.0) value = 0.0;
    output[index] = value;
    recordNonFinite(value, error_flag);
    if (Reflect) {
      const double doubled = __dmul_rn(2.0, value);
      const double reflection = __dsub_rn(doubled, current[index]);
      reflected[index] = reflection;
      recordNonFinite(doubled, error_flag);
      recordNonFinite(reflection, error_flag);
    }
  }
}

__global__ void vectorBlendKernel(const HighsInt size, const double* left,
                                  const double* right, const double literal_weight,
                                  double* output, int* error_flag,
                                  const double* weight_parameter) {
  const double left_weight = weight_parameter == nullptr
      ? literal_weight : *weight_parameter;
  if (weight_parameter != nullptr &&
      (!isfinite(left_weight) || left_weight < 0.0 || left_weight > 1.0)) {
    atomicExch(error_flag, 1);
    return;
  }
  const double right_weight = 1.0 - left_weight;
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    const double value =
        left_weight * left[index] + right_weight * right[index];
    output[index] = value;
    recordNonFinite(value, error_flag);
  }
}

__global__ void vectorPairBlendKernel(const HighsInt primal_size,
    const HighsInt dual_size, const double* primal_left, const double* primal_right,
    const double* dual_left, const double* dual_right, const double literal_weight,
    double* primal_output, double* dual_output, int* error_flag,
    const double* weight_parameter) {
  const double left_weight = weight_parameter == nullptr ? literal_weight : *weight_parameter;
  if (!isfinite(left_weight) || left_weight < 0.0 || left_weight > 1.0) {
    atomicExch(error_flag, 1);
    return;
  }
  const double right_weight = 1.0 - left_weight;
  const HighsInt size = primal_size > dual_size ? primal_size : dual_size;
  const GridIndex stride = static_cast<GridIndex>(blockDim.x) * gridDim.x;
  for (GridIndex index = static_cast<GridIndex>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < static_cast<GridIndex>(size); index += stride) {
    if (index < static_cast<GridIndex>(primal_size)) {
      const double value = left_weight * primal_left[index] + right_weight * primal_right[index];
      primal_output[index] = value;
      recordNonFinite(value, error_flag);
    }
    if (index < static_cast<GridIndex>(dual_size)) {
      const double value = left_weight * dual_left[index] + right_weight * dual_right[index];
      dual_output[index] = value;
      recordNonFinite(value, error_flag);
    }
  }
}

__global__ void averageKernel(const HighsInt size, const double* current,
                              const double weight, double* sum,
                              int* error_flag) {
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    const double value = sum[index] + weight * current[index];
    sum[index] = value;
    recordNonFinite(value, error_flag);
  }
}

__global__ void finiteCheckKernel(const HighsInt size, const double* values,
                                  int* error_flag) {
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride)
    recordNonFinite(values[index], error_flag);
}

__global__ void dotKernel(const HighsInt size, const double* left,
                          const double* right, double* result,
                          int* error_flag) {
  if (blockIdx.x != 0 || threadIdx.x != 0) return;
  double value = 0.0;
  for (HighsInt index = 0; index < size; ++index)
    value += left[index] * right[index];
  *result = value;
  recordNonFinite(value, error_flag);
}

__device__ double scaledNorm(const HighsInt size, const double* values,
                             const double* right, const bool difference) {
  double scale = 0.0;
  double sum_of_squares = 1.0;
  for (HighsInt index = 0; index < size; ++index) {
    const double value =
        difference ? values[index] - right[index] : values[index];
    const double magnitude = fabs(value);
    if (magnitude == 0.0) continue;
    if (scale < magnitude) {
      const double ratio = scale / magnitude;
      sum_of_squares = 1.0 + sum_of_squares * ratio * ratio;
      scale = magnitude;
    } else {
      const double ratio = magnitude / scale;
      sum_of_squares += ratio * ratio;
    }
  }
  return scale == 0.0 ? 0.0 : scale * sqrt(sum_of_squares);
}

__global__ void normKernel(const HighsInt size, const double* values,
                           double* result, int* error_flag) {
  if (blockIdx.x != 0 || threadIdx.x != 0) return;
  const double value = scaledNorm(size, values, nullptr, false);
  *result = value;
  recordNonFinite(value, error_flag);
}

__global__ void differenceNormKernel(const HighsInt size, const double* left,
                                     const double* right, double* result,
                                     int* error_flag) {
  if (blockIdx.x != 0 || threadIdx.x != 0) return;
  const double value = scaledNorm(size, left, right, true);
  *result = value;
  recordNonFinite(value, error_flag);
}

const int kMaximumCheckBlocks = 1024;

enum CheckNormIndex {
  kCheckPrimalResidual,
  kCheckDualResidual,
  kCheckPrimalDelta,
  kCheckDualDelta,
  kCheckPrimalAnchor,
  kCheckDualAnchor
};

enum CheckSumIndex { kCheckPrimalObjective, kCheckDualObjective, kCheckDeltaDot };

__device__ void requestCheckFallback(CheckPartial& partial, const int reason) {
  if (reason == kCheckRawFallback) partial.requires_host_oracle = 1;
  else partial.optional_invalid |= reason;
}

__device__ void addCheckNorm(CheckNormState& norm, const double value,
                             CheckPartial& partial, const int reason = kCheckRawFallback) {
  if (!isfinite(value)) {
    requestCheckFallback(partial, reason);
    if (reason == kCheckFixedPointFallback)
      partial.optional_invalid |= kCheckUnpreparedProducts;
    return;
  }
  const double magnitude = fabs(value);
  if (magnitude == 0.0) return;
  if (norm.scale < magnitude) {
    const double ratio = __ddiv_rn(norm.scale, magnitude);
    norm.sum_of_squares = __dadd_rn(
        1.0, __dmul_rn(__dmul_rn(norm.sum_of_squares, ratio), ratio));
    norm.scale = magnitude;
  } else {
    const double ratio = __ddiv_rn(magnitude, norm.scale);
    norm.sum_of_squares =
        __dadd_rn(norm.sum_of_squares, __dmul_rn(ratio, ratio));
  }
}

__device__ void mergeCheckNorm(CheckNormState& left,
                               const CheckNormState& right) {
  if (right.scale == 0.0) return;
  if (left.scale < right.scale) {
    const double ratio = __ddiv_rn(left.scale, right.scale);
    left.sum_of_squares = __dadd_rn(
        right.sum_of_squares,
        __dmul_rn(__dmul_rn(left.sum_of_squares, ratio), ratio));
    left.scale = right.scale;
  } else {
    const double ratio = __ddiv_rn(right.scale, left.scale);
    left.sum_of_squares = __dadd_rn(
        left.sum_of_squares,
        __dmul_rn(__dmul_rn(right.sum_of_squares, ratio), ratio));
  }
}

__device__ void mergeFixedPointPartial(FixedPointPartial& left,
                                        const FixedPointPartial& right) {
  mergeCheckNorm(left.primal, right.primal);
  mergeCheckNorm(left.dual, right.dual);
  left.cross = __dadd_rn(left.cross, right.cross);
  left.invalid |= right.invalid || !isfinite(left.cross);
}

__device__ FixedPointPartial reduceFixedPointBlock(FixedPointPartial partial) {
  __shared__ FixedPointPartial shared[kBlockSize];
  shared[threadIdx.x] = partial;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      mergeFixedPointPartial(shared[threadIdx.x], shared[threadIdx.x + distance]);
    __syncthreads();
  }
  return shared[0];
}

__global__ void fixedPointPartialsKernel(
    const HighsInt primal_size, const HighsInt dual_size,
    const double* primal, const double* dual, const double* product,
    FixedPointPartial* partials) {
  FixedPointPartial partial = {};
  CheckPartial validity = {};
  const GridIndex size = primal_size > dual_size ? primal_size : dual_size;
  const GridIndex stride = static_cast<GridIndex>(blockDim.x) * gridDim.x;
  for (GridIndex index = static_cast<GridIndex>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    if (index < static_cast<GridIndex>(primal_size)) {
      const double value = primal[index];
      addCheckNorm(partial.primal, value, validity, kCheckFixedPointFallback);
      partial.cross = __dadd_rn(partial.cross, __dmul_rn(value, product[index]));
      partial.invalid |= !isfinite(partial.cross);
    }
    if (index < static_cast<GridIndex>(dual_size))
      addCheckNorm(partial.dual, dual[index], validity, kCheckFixedPointFallback);
  }
  partial.invalid |= validity.optional_invalid != 0;
  const FixedPointPartial total = reduceFixedPointBlock(partial);
  if (threadIdx.x == 0) partials[blockIdx.x] = total;
}

__global__ void fixedPointResultKernel(const int count,
                                       const FixedPointPartial* partials,
                                       double* result, int* error_flag) {
  FixedPointPartial partial = {};
  for (int index = threadIdx.x; index < count; index += blockDim.x)
    mergeFixedPointPartial(partial, partials[index]);
  const FixedPointPartial total = reduceFixedPointBlock(partial);
  if (threadIdx.x != 0) return;
  result[0] = __dmul_rn(total.primal.scale, __dsqrt_rn(total.primal.sum_of_squares));
  result[1] = __dmul_rn(total.dual.scale, __dsqrt_rn(total.dual.sum_of_squares));
  result[2] = total.cross;
  if (total.invalid) atomicExch(error_flag, 1);
  for (int index = 0; index < 3; ++index) recordNonFinite(result[index], error_flag);
}

__device__ void addCheckSum(CheckSumState& sum, const double value,
                            CheckPartial& partial, const int reason = kCheckRawFallback) {
  if (!isfinite(value)) {
    requestCheckFallback(partial, reason);
    return;
  }
  const double magnitude = fabs(value);
  if (magnitude > 0.0 &&
      (sum.minimum_nonzero == 0.0 || magnitude < sum.minimum_nonzero))
    sum.minimum_nonzero = magnitude;
  const double combined = __dadd_rn(sum.value, value);
  const double absolute = __dadd_rn(sum.absolute_sum, fabs(value));
  if (!isfinite(combined) || !isfinite(absolute))
    requestCheckFallback(partial, reason);
  sum.value = isfinite(combined) ? combined : 0.0;
  sum.absolute_sum = isfinite(absolute) ? absolute : 0.0;
}

__device__ void mergeCheckPartial(CheckPartial& left,
                                  const CheckPartial& right) {
  left.numerical_error |= right.numerical_error;
  left.optional_invalid |= right.optional_invalid;
  left.requires_host_oracle |= right.requires_host_oracle;
  for (int index = 0; index < 6; ++index)
    mergeCheckNorm(left.norms[index], right.norms[index]);
  for (int index = 0; index < 3; ++index) {
    const double minimum = right.sums[index].minimum_nonzero;
    if (minimum > 0.0 && (left.sums[index].minimum_nonzero == 0.0 ||
                          minimum < left.sums[index].minimum_nonzero))
      left.sums[index].minimum_nonzero = minimum;
    const double combined = __dadd_rn(left.sums[index].value, right.sums[index].value);
    const double absolute = __dadd_rn(left.sums[index].absolute_sum,
                                     right.sums[index].absolute_sum);
    if (!isfinite(combined) || !isfinite(absolute))
      requestCheckFallback(left, index == kCheckDeltaDot
                                     ? kCheckFixedPointFallback : kCheckRawFallback);
    left.sums[index].value = isfinite(combined) ? combined : 0.0;
    left.sums[index].absolute_sum = isfinite(absolute) ? absolute : 0.0;
  }
}

__device__ double finishCheckNorm(const CheckNormState& norm,
                                 CheckPartial& partial,
                                 const int reason = kCheckRawFallback) {
  const double value = norm.scale == 0.0
                           ? 0.0
                           : __dmul_rn(norm.scale, __dsqrt_rn(norm.sum_of_squares));
  if (!isfinite(value))
    requestCheckFallback(partial, reason);
  return value;
}

// interval_guarded: guardCheckStopping will judge this sum against a strict
// a-priori interval, so the heuristics below only cover what that interval
// does not.
//
// Overflow and cancellation always fall back. Absorption of a small term into
// a large partial sum (minimum_nonzero <= eps * sum|a_i|) falls back only for
// the objective sums, and only when no interval will judge them. When one will,
// absorption is redundant:
// every summation order, serial or tree, satisfies
//   |fl(sum a_i) - sum a_i| <= gamma_k * sum|a_i|,
// and a rounding that swallows a small term is one of the k roundings that
// bound already counts. The host oracle and the GPU each lie within that
// distance of the exact sum, so guardCheckStopping's interval
// (relative_error * sum|a_i| + a subnormal floor) contains both, and it falls
// back whenever the interval straddles a tolerance. Absorbed terms therefore
// cannot flip a stopping decision unnoticed; at worst they widen the interval.
// The bound only breaks in the subnormal range, where a rounding error is
// absolute rather than relative; the interval carries a floor for it, and a
// subnormal term still falls back here.
//
// The fixed-point sum delta_dot has no interval, but it needs none for
// absorption: its fallback recomputes the same dot product with backend
// primitives, another rounding of the same sum, so it can differ from the
// packet by no more than the bound above. It falls back on overflow and
// cancellation only.
__device__ double finishCheckSum(CheckSumState& sum, const double offset,
                                const double count, CheckPartial& partial,
                                const int reason = kCheckRawFallback,
                                const bool interval_guarded = false) {
  addCheckSum(sum, offset, partial, reason);
  const double cancellation_limit =
      __dmul_rn(__dmul_rn(64.0 * DBL_EPSILON, count + 1.0), sum.absolute_sum);
  const bool has_term = reason == kCheckRawFallback && sum.minimum_nonzero > 0.0;
  const bool subnormal_term = has_term && sum.minimum_nonzero < DBL_MIN;
  const bool absorbed_term =
      has_term &&
      sum.minimum_nonzero <= __dmul_rn(DBL_EPSILON, sum.absolute_sum);
  if (sum.absolute_sum > DBL_MAX / 4.0 ||
      (sum.absolute_sum > 0.0 && fabs(sum.value) <= cancellation_limit) ||
      (interval_guarded ? subnormal_term : absorbed_term))
    requestCheckFallback(partial, reason);
  return sum.value;
}

  __device__ certificate::Interval checkRoundingInterval(
    const double value, const double magnitude, const double relative_error,
    const double absolute_error) {
    const double error = certificate::addUp(
      certificate::multiplyUp(relative_error, magnitude), absolute_error);
    return {certificate::addDown(value, -error), certificate::addUp(value, error)};
  }

  __device__ certificate::Interval checkNormInterval(
    const CheckNormState& norm, const double value, const double relative_error,
    const double scaled_floor, const double absolute_floor) {
    if (norm.scale == 0.0) return {0.0, 0.0};
    const double absolute_error = certificate::addUp(
      certificate::multiplyUp(norm.scale, scaled_floor), absolute_floor);
    certificate::Interval interval =
      checkRoundingInterval(value, value, relative_error, absolute_error);
    interval.lower = certificate::maximum(0.0, interval.lower);
    return interval;
  }

  __device__ bool checkThresholdUncertain(const certificate::Interval value,
                       const double tolerance) {
    return !certificate::finite(value) ||
       (value.lower <= tolerance && value.upper > tolerance);
  }

  __device__ void guardCheckStopping(const DeviceCheckProblem& problem,
                     const DeviceCheckInput& input,
                     const CheckPacket& packet,
                     CheckPartial& total) {
    using namespace certificate;
    if (total.requires_host_oracle || total.numerical_error ||
      (input.primal_tolerance == 0.0 && input.dual_tolerance == 0.0 &&
       input.gap_tolerance == 0.0))
    return;
    const double dimension = addUp(
      addUp(static_cast<double>(problem.num_rows),
        static_cast<double>(problem.num_cols)), 1.0);
    // Five rounding factors per scaled merge, plus sqrt/product; eight per
    // component also bounds products, offset and serial/tree objective additions.
    const double operations = multiplyUp(8.0, dimension);
    const double gamma = errorGamma(operations);
    const double relative_error = errorGamma(multiplyUp(2.0, operations));
    if (!finiteValue(relative_error) || relative_error >= 1.0) {
    total.requires_host_oracle = 1;
    return;
    }
    const double floor_amplification = divideUp(2.0, addDown(1.0, -gamma));
    // At most 12*dimension denorm_min absolute error per rescaled merge;
    // there are fewer than dimension nontrivial merges, in either ordering.
    const double scaled_floor = multiplyUp(
      floor_amplification,
      subnormalError(multiplyUp(12.0, multiplyUp(dimension, dimension)), gamma));
    const double norm_floor =
      multiplyUp(floor_amplification, subnormalError(2.0, gamma));
    if (input.primal_tolerance > 0.0) {
    const Interval residual = checkNormInterval(
      total.norms[kCheckPrimalResidual], packet.primal_residual,
      relative_error, scaled_floor, norm_floor);
    const Interval denominator =
      add({1.0, 1.0}, {packet.unscaled_rhs_norm, packet.unscaled_rhs_norm});
    if (checkThresholdUncertain(divide(residual, denominator),
                  input.primal_tolerance))
      total.requires_host_oracle = 1;
    }
    if (input.dual_tolerance > 0.0) {
    const Interval residual = checkNormInterval(
      total.norms[kCheckDualResidual], packet.dual_residual,
      relative_error, scaled_floor, norm_floor);
    const Interval denominator =
      add({1.0, 1.0}, {packet.unscaled_cost_norm, packet.unscaled_cost_norm});
    if (checkThresholdUncertain(divide(residual, denominator),
                  input.dual_tolerance))
      total.requires_host_oracle = 1;
    }
    if (input.gap_tolerance > 0.0) {
    const double objective_floor =
      multiplyUp(floor_amplification, subnormalError(operations, gamma));
    const Interval primal = checkRoundingInterval(
      packet.primal_objective,
      total.sums[kCheckPrimalObjective].absolute_sum,
      relative_error, objective_floor);
    const Interval dual = checkRoundingInterval(
      packet.dual_objective, total.sums[kCheckDualObjective].absolute_sum,
      relative_error, objective_floor);
    const Interval difference = add(primal, multiply(dual, -1.0));
    const Interval gap = {absoluteLower(difference), absoluteUpper(difference)};
    const Interval denominator =
      add(add({1.0, 1.0}, {absoluteLower(primal), absoluteUpper(primal)}),
        {absoluteLower(dual), absoluteUpper(dual)});
    if (checkThresholdUncertain(divide(gap, denominator), input.gap_tolerance))
      total.requires_host_oracle = 1;
    }
  }

__global__ void checkPartialsKernel(const DeviceCheckProblem problem,
                                    const DeviceCheckInput input,
                                    CheckPartial* partials) {
  CheckPartial partial = {};
  const HighsInt size = problem.num_rows > problem.num_cols
                            ? problem.num_rows : problem.num_cols;
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index = static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    if (index < problem.num_rows) {
      const double dual = input.dual[index];
      const double ax = input.ax[index];
      if (!isfinite(dual) || !isfinite(ax)) partial.numerical_error = 1;
      double violation = __dsub_rn(ax, problem.row_lower[index]);
      if (!problem.equality_rows[index] && violation > 0.0) violation = 0.0;
      addCheckNorm(partial.norms[kCheckPrimalResidual],
                   __ddiv_rn(violation, problem.row_scale[index]), partial);
      addCheckSum(partial.sums[kCheckDualObjective],
                  __dmul_rn(problem.row_lower[index], dual), partial);
      if (input.fixed_point)
        addCheckNorm(partial.norms[kCheckDualDelta], input.delta_dual[index],
                     partial, kCheckFixedPointFallback);
      if (input.anchor_distances) {
        const double anchor = input.dual_anchor[index];
        if (!isfinite(anchor)) partial.optional_invalid |= kCheckAnchorFallback;
        addCheckNorm(partial.norms[kCheckDualAnchor], __dsub_rn(dual, anchor),
               partial, kCheckAnchorFallback);
      }
    }
    if (index < problem.num_cols) {
      const double primal = input.primal[index];
      const double aty = input.aty[index];
      if (!isfinite(primal) || !isfinite(aty)) partial.numerical_error = 1;
      const double reduced_cost = __dsub_rn(problem.costs[index], aty);
      if (!isfinite(reduced_cost)) partial.requires_host_oracle = 1;
      const bool has_lower = isfinite(problem.col_lower[index]);
      const bool has_upper = isfinite(problem.col_upper[index]);
      const double lower_slack = has_lower && reduced_cost > 0.0 ? reduced_cost : 0.0;
      const double upper_slack = has_upper && reduced_cost < 0.0 ? -reduced_cost : 0.0;
      const double residual =
          __dadd_rn(__dsub_rn(reduced_cost, lower_slack), upper_slack);
      addCheckNorm(partial.norms[kCheckDualResidual],
                   __ddiv_rn(residual, problem.col_scale[index]), partial);
      addCheckSum(partial.sums[kCheckPrimalObjective],
                  __dmul_rn(problem.costs[index], primal), partial);
      if (has_lower)
        addCheckSum(partial.sums[kCheckDualObjective],
                    __dmul_rn(problem.col_lower[index], lower_slack), partial);
      if (has_upper)
        addCheckSum(partial.sums[kCheckDualObjective],
                    -__dmul_rn(problem.col_upper[index], upper_slack), partial);
      if (input.fixed_point) {
        const double delta = input.delta_primal[index];
        const double product = input.aty_delta[index];
        if (!isfinite(product))
          partial.optional_invalid |= kCheckFixedPointFallback | kCheckUnpreparedProducts;
        addCheckNorm(partial.norms[kCheckPrimalDelta], delta, partial, kCheckFixedPointFallback);
        addCheckSum(partial.sums[kCheckDeltaDot], __dmul_rn(delta, product),
                    partial, kCheckFixedPointFallback);
      }
      if (input.anchor_distances) {
        const double anchor = input.primal_anchor[index];
        if (!isfinite(anchor)) partial.optional_invalid |= kCheckAnchorFallback;
        addCheckNorm(partial.norms[kCheckPrimalAnchor], __dsub_rn(primal, anchor),
               partial, kCheckAnchorFallback);
      }
    }
  }

  __shared__ CheckPartial shared[kBlockSize];
  shared[threadIdx.x] = partial;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      mergeCheckPartial(shared[threadIdx.x], shared[threadIdx.x + distance]);
    __syncthreads();
  }
  if (threadIdx.x == 0) partials[blockIdx.x] = shared[0];
}

__global__ void checkPacketKernel(const DeviceCheckProblem problem,
                                  const DeviceCheckInput input,
                                  const int count, const CheckPartial* partials,
                                  CheckPacket* output, int* error_flag) {
  CheckPartial partial = {};
  for (int index = threadIdx.x; index < count; index += blockDim.x)
    mergeCheckPartial(partial, partials[index]);
  __shared__ CheckPartial shared[kBlockSize];
  shared[threadIdx.x] = partial;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      mergeCheckPartial(shared[threadIdx.x], shared[threadIdx.x + distance]);
    __syncthreads();
  }
  if (threadIdx.x != 0) return;
  CheckPartial& total = shared[0];
  CheckPacket packet = {};
  packet.primal_residual = finishCheckNorm(total.norms[kCheckPrimalResidual], total);
  packet.dual_residual = finishCheckNorm(total.norms[kCheckDualResidual], total);
  // The gap interval in guardCheckStopping judges both objective sums.
  const bool objective_guarded = input.gap_tolerance > 0.0;
  packet.primal_objective = finishCheckSum(
      total.sums[kCheckPrimalObjective], problem.objective_offset,
      static_cast<double>(problem.num_cols), total, kCheckRawFallback,
      objective_guarded);
  packet.dual_objective = finishCheckSum(
      total.sums[kCheckDualObjective], problem.objective_offset,
      static_cast<double>(problem.num_rows) + 2.0 * problem.num_cols, total,
      kCheckRawFallback, objective_guarded);
  packet.objective_gap = fabs(__dsub_rn(packet.primal_objective, packet.dual_objective));
  if (!isfinite(packet.objective_gap)) total.requires_host_oracle = 1;
  packet.unscaled_rhs_norm = problem.unscaled_rhs_norm;
  packet.unscaled_cost_norm = problem.unscaled_cost_norm;
  guardCheckStopping(problem, input, packet, total);
  if (input.fixed_point) {
    packet.primal_delta_norm =
      finishCheckNorm(total.norms[kCheckPrimalDelta], total, kCheckFixedPointFallback);
    packet.dual_delta_norm =
      finishCheckNorm(total.norms[kCheckDualDelta], total, kCheckFixedPointFallback);
    packet.delta_dot = finishCheckSum(total.sums[kCheckDeltaDot], 0.0,
                      static_cast<double>(problem.num_cols), total,
                      kCheckFixedPointFallback);
  }
  if (input.anchor_distances) {
    packet.primal_anchor_distance =
      finishCheckNorm(total.norms[kCheckPrimalAnchor], total, kCheckAnchorFallback);
    packet.dual_anchor_distance =
      finishCheckNorm(total.norms[kCheckDualAnchor], total, kCheckAnchorFallback);
  }
  if (*error_flag == 0 && total.numerical_error != 0) *error_flag = 1;
  packet.error_flag = *error_flag;
  packet.fallback_flags = (total.requires_host_oracle ? kCheckRawFallback : 0) |
                           total.optional_invalid;
  if (input.fixed_point && *input.optional_error != 0)
    packet.fallback_flags |= kCheckFixedPointFallback | kCheckUnpreparedProducts;
  packet.primal_ray_possible = 1;
  packet.farkas_possible = 1;
  *output = packet;
}

__device__ void mergeCertificateNorm(CertificateNormPartial& left,
                                     const CertificateNormPartial& right) {
  certificate::mergeNorm(left.primal, right.primal);
  certificate::mergeNorm(left.dual, right.dual);
  left.unsafe |= right.unsafe;
}

__global__ void certificateRecoveryKernel(
    const DeviceCertificateProblem problem, const DeviceCheckInput input,
    certificate::Interval* primal, certificate::Interval* dual,
    certificate::Interval* column, CertificateNormPartial* partials) {
  CertificateNormPartial partial = {};
  const HighsInt size =
      problem.num_rows > problem.num_cols ? problem.num_rows : problem.num_cols;
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    if (index < problem.num_cols) {
      primal[index] = certificate::recoverProduct(input.primal[index],
                                                  problem.col_scale[index]);
      column[index] = certificate::recoverQuotient(-input.aty[index],
                                                   problem.col_scale[index]);
      partial.unsafe |=
          !certificate::safeRecovery(primal[index], input.primal[index]) ||
          !certificate::safeRecovery(column[index], -input.aty[index]);
      certificate::addNorm(partial.primal, primal[index]);
      certificate::addNorm(partial.dual, column[index]);
    }
    if (index < problem.num_rows) {
      const HighsInt transformed = problem.row_map[index];
      dual[index] = certificate::multiply(
          certificate::recoverProduct(input.dual[transformed],
                                      problem.row_scale[transformed]),
          problem.row_sign[index]);
      partial.unsafe |=
          !certificate::safeRecovery(dual[index], input.dual[transformed]);
      certificate::addNorm(partial.dual, dual[index]);
    }
  }
  __shared__ CertificateNormPartial shared[kBlockSize];
  shared[threadIdx.x] = partial;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      mergeCertificateNorm(shared[threadIdx.x], shared[threadIdx.x + distance]);
    __syncthreads();
  }
  if (threadIdx.x == 0) partials[blockIdx.x] = shared[0];
}

__global__ void certificateNormKernel(const DeviceCertificateProblem problem,
                                      const int count,
                                      const CertificateNormPartial* partials,
                                      const bool verifier,
                                      CertificateState* state) {
  CertificateNormPartial partial = {};
  for (int index = threadIdx.x; index < count; index += blockDim.x)
    mergeCertificateNorm(partial, partials[index]);
  __shared__ CertificateNormPartial shared[kBlockSize];
  shared[threadIdx.x] = partial;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      mergeCertificateNorm(shared[threadIdx.x], shared[threadIdx.x + distance]);
    __syncthreads();
  }
  if (threadIdx.x != 0) return;
  const certificate::Interval primal = certificate::hostNormBounds(
      shared[0].primal, static_cast<double>(problem.num_cols));
  const certificate::Interval dual = certificate::hostNormBounds(
      shared[0].dual,
      static_cast<double>(problem.num_rows) + problem.num_cols + 1.0);
  const int unsafe = shared[0].unsafe || !certificate::positiveFinite(primal) ||
                     !certificate::positiveFinite(dual);
  if (verifier) {
    state->primal_verifier_norm = primal;
    state->dual_verifier_norm = dual;
    state->unsafe |= unsafe;
  } else {
    state->primal_normalization = primal;
    state->dual_normalization = dual;
    state->unsafe = unsafe;
  }
}

__global__ void certificateNormalizeKernel(
    const DeviceCertificateProblem problem, const CertificateState* state,
    certificate::Interval* primal, certificate::Interval* dual,
    certificate::Interval* column, CertificateNormPartial* partials) {
  CertificateNormPartial partial = {};
  partial.unsafe = state->unsafe;
  const HighsInt size =
      problem.num_rows > problem.num_cols ? problem.num_rows : problem.num_cols;
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  if (!state->unsafe) {
    for (HighsInt index =
             static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
         index < size; index += stride) {
      if (index < problem.num_cols) {
        primal[index] =
            certificate::divide(primal[index], state->primal_normalization);
        column[index] =
            certificate::divide(column[index], state->dual_normalization);
        partial.unsafe |= !certificate::finite(primal[index]) ||
                          !certificate::finite(column[index]);
        certificate::addNorm(partial.primal, primal[index]);
        certificate::addNorm(partial.dual, column[index]);
      }
      if (index < problem.num_rows) {
        dual[index] =
            certificate::divide(dual[index], state->dual_normalization);
        partial.unsafe |= !certificate::finite(dual[index]);
        certificate::addNorm(partial.dual, dual[index]);
      }
    }
  }
  __shared__ CertificateNormPartial shared[kBlockSize];
  shared[threadIdx.x] = partial;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      mergeCertificateNorm(shared[threadIdx.x], shared[threadIdx.x + distance]);
    __syncthreads();
  }
  if (threadIdx.x == 0) partials[blockIdx.x] = shared[0];
}

// One warp per row segment: the lanes stride over the segment and their
// outward-rounded partial sums are combined with outward-rounded additions,
// which encloses the activity for any summation order. Splitting a long row
// into segments, and adding the segment sums afterwards, is one more order.
__global__ void certificateRowActivityKernel(
    const DeviceCertificateProblem problem, const CertificateState* state,
    const certificate::Interval* primal, certificate::Interval* activity) {
  using namespace certificate;
  if (state->unsafe || state->rows_skipped) return;
  const HighsInt lane = static_cast<HighsInt>(threadIdx.x) % kWarpSize;
  const HighsInt warps =
      static_cast<HighsInt>(gridDim.x) * (blockDim.x / kWarpSize);
  for (HighsInt segment = (static_cast<HighsInt>(blockIdx.x) * blockDim.x +
                           threadIdx.x) / kWarpSize;
       segment < problem.num_segments; segment += warps) {
    const HighsInt row = problem.segment_row[segment];
    const HighsInt begin = problem.segment_begin[segment];
    const HighsInt row_end = problem.row_starts[row + 1];
    const HighsInt end = row_end - begin < kCertificateSegmentLength
                             ? row_end : begin + kCertificateSegmentLength;
    Interval sum = {0.0, 0.0};
    for (HighsInt entry = begin + lane; entry < end; entry += kWarpSize)
      sum = add(sum, multiply(primal[problem.column_indices[entry]],
                              problem.matrix_values[entry]));
    for (int offset = kWarpSize / 2; offset > 0; offset /= 2) {
      const Interval other = {
          __shfl_down_sync(0xffffffffu, sum.lower, offset),
          __shfl_down_sync(0xffffffffu, sum.upper, offset)};
      sum = add(sum, other);
    }
    if (lane == 0) {
      if (row_end - problem.row_starts[row] <= kCertificateSegmentLength)
        activity[row] = sum;
      else
        problem.segment_partial[segment] = sum;
    }
  }
}

// Adds the segment sums of each row that was split.
__global__ void certificateLongRowKernel(
    const DeviceCertificateProblem problem, const CertificateState* state,
    certificate::Interval* activity) {
  using namespace certificate;
  if (state->unsafe || state->rows_skipped) return;
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index = static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < problem.num_long_rows; index += stride) {
    Interval sum = {0.0, 0.0};
    const HighsInt first = problem.long_first[index];
    for (HighsInt part = 0; part < problem.long_count[index]; ++part)
      sum = add(sum, problem.segment_partial[first + part]);
    activity[problem.long_row[index]] = sum;
  }
}

// Merges the per-block screen partials; the result is valid in thread 0.
__device__ certificate::ScreenPartial mergeScreenPartials(
    const certificate::ScreenPartial* partials, const int count) {
  certificate::ScreenPartial partial = {};
  for (int index = threadIdx.x; index < count; index += blockDim.x)
    certificate::mergeScreen(partial, partials[index]);
  __shared__ certificate::ScreenPartial shared[kBlockSize];
  shared[threadIdx.x] = partial;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      certificate::mergeScreen(shared[threadIdx.x],
                               shared[threadIdx.x + distance]);
    __syncthreads();
  }
  return shared[0];
}

__global__ void certificateConditionsKernel(
    const DeviceCertificateProblem problem, const CertificateState* state,
    const certificate::Interval* primal, const certificate::Interval* dual,
    const certificate::Interval* column,
    certificate::ScreenPartial* partials) {
  using namespace certificate;
  ScreenPartial partial = {};
  partial.unsafe = state->unsafe;
  const HighsInt size =
      problem.num_rows > problem.num_cols ? problem.num_rows : problem.num_cols;
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  if (!state->unsafe) {
    for (HighsInt index =
             static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
         index < size; index += stride) {
      if (index < problem.num_cols) {
        partial.primal_violation =
            maximum(partial.primal_violation,
                    violationLower(primal[index], problem.col_lower[index],
                                   problem.col_upper[index], false));
        partial.farkas_violation =
            maximum(partial.farkas_violation,
                    violationLower(column[index], problem.col_lower[index],
                                   problem.col_upper[index], true));
        addSum(partial.objective,
               multiply(primal[index], problem.costs[index]));
      }
      if (index < problem.num_rows) {
        partial.farkas_violation =
            maximum(partial.farkas_violation,
                    violationLower(dual[index], problem.row_lower[index],
                                   problem.row_upper[index], true));
        if (problem.column_bounds_straddle_zero)
          addSum(partial.row_support,
                 rowSupport(dual[index], problem.row_lower[index],
                            problem.row_upper[index]));
      }
    }
  }
  __shared__ ScreenPartial shared[kBlockSize];
  shared[threadIdx.x] = partial;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      mergeScreen(shared[threadIdx.x], shared[threadIdx.x + distance]);
    __syncthreads();
  }
  if (threadIdx.x == 0) partials[blockIdx.x] = shared[0];
}

// Decides from the column part and the objective alone whether a primal ray is
// already ruled out. finishScreen rejects a ray when either the violation or
// the objective test rejects, and the violation is a maximum over columns and
// rows, so a rejection that holds without the rows holds with them and the
// final flags cannot change.
__global__ void certificateRaySkipKernel(
    const DeviceCertificateProblem problem, const DeviceCheckInput input,
    const int count, const certificate::ScreenPartial* partials,
    CertificateState* state) {
  const certificate::ScreenPartial partial =
      mergeScreenPartials(partials, count);
  if (threadIdx.x != 0) return;
  state->row_violation_bits = 0;
  state->rows_skipped =
      !state->unsafe && !partial.unsafe &&
      (certificate::rejectViolation(partial.primal_violation,
                                    state->primal_verifier_norm,
                                    input.certificate_tolerance) ||
       certificate::rejectObjective(
           certificate::serialSumBounds(partial.objective,
                                        static_cast<double>(problem.num_cols)),
           state->primal_verifier_norm, problem.sense_sign,
           input.certificate_tolerance));
}

__global__ void certificateRowViolationKernel(
    const DeviceCertificateProblem problem, const CertificateState* state,
    const certificate::Interval* activities, CertificateState* output) {
  using namespace certificate;
  if (state->unsafe || state->rows_skipped) return;
  double violation = 0.0;
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt row = static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       row < problem.num_rows; row += stride)
    violation = maximum(violation,
                        violationLower(activities[row], problem.row_lower[row],
                                       problem.row_upper[row], false));
  __shared__ double shared[kBlockSize];
  shared[threadIdx.x] = violation;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      shared[threadIdx.x] = maximum(shared[threadIdx.x],
                                    shared[threadIdx.x + distance]);
    __syncthreads();
  }
  if (threadIdx.x == 0 && shared[0] > 0.0)
    atomicMax(&output->row_violation_bits,
              static_cast<unsigned long long>(__double_as_longlong(shared[0])));
}

__global__ void certificateFlagsKernel(
    const DeviceCertificateProblem problem, const DeviceCheckInput input,
    const int count, const CertificateState* state,
    const certificate::ScreenPartial* partials, CheckPacket* packet) {
  certificate::ScreenPartial partial = mergeScreenPartials(partials, count);
  if (threadIdx.x != 0) return;
  partial.primal_violation = certificate::maximum(
      partial.primal_violation, __longlong_as_double(
                                    static_cast<long long>(state->row_violation_bits)));
  const certificate::ScreenFlags flags = certificate::finishScreen(
      partial, state->primal_verifier_norm, state->dual_verifier_norm,
      static_cast<double>(problem.num_cols),
      static_cast<double>(problem.num_rows), problem.sense_sign,
      problem.column_bounds_straddle_zero,
      state->unsafe || packet->error_flag != 0, input.certificate_tolerance);
  packet->primal_ray_possible = flags.primal_ray_possible;
  packet->farkas_possible = flags.farkas_possible;
}

#ifdef HIGHS_HIPDLP_TESTING
__global__ void injectedFailureKernel(int* error_flag) {
  if (blockIdx.x != 0 || threadIdx.x != 0) return;
  atomicExch(error_flag, 1);
}
#endif


// Original-model KKT screen. Each quantity carries an error bound covering
// the host verifier's different summation order and long-double unscaling, so
// the reported values are lower bounds of the verifier's metrics.
const int kScreenWarps = kBlockSize / 32;
const int kMaximumScreenBlocks = 1024;

__device__ void mergeScreenPartial(KktScreenPartial& left,
                                   const KktScreenPartial& right) {
  left.active_bound = fmax(left.active_bound, right.active_bound);
  left.near_zero_cost = fmax(left.near_zero_cost, right.near_zero_cost);
  left.primal_objective += right.primal_objective;
  left.primal_absolute += right.primal_absolute;
  left.dual_objective += right.dual_objective;
  left.dual_absolute += right.dual_absolute;
  left.dual_error += right.dual_error;
  left.complementarity += right.complementarity;
  left.invalid |= right.invalid;
}

// Reduces every thread's partial in a fixed order; thread 0 writes the block's.
__device__ void writeScreenBlock(KktScreenPartial partial,
                                 KktScreenPartial* partials) {
  for (int offset = 16; offset > 0; offset >>= 1) {
    partial.active_bound = fmax(partial.active_bound,
        __shfl_down_sync(0xffffffffu, partial.active_bound, offset));
    partial.near_zero_cost = fmax(partial.near_zero_cost,
        __shfl_down_sync(0xffffffffu, partial.near_zero_cost, offset));
    partial.primal_objective +=
        __shfl_down_sync(0xffffffffu, partial.primal_objective, offset);
    partial.primal_absolute +=
        __shfl_down_sync(0xffffffffu, partial.primal_absolute, offset);
    partial.dual_objective +=
        __shfl_down_sync(0xffffffffu, partial.dual_objective, offset);
    partial.dual_absolute +=
        __shfl_down_sync(0xffffffffu, partial.dual_absolute, offset);
    partial.dual_error += __shfl_down_sync(0xffffffffu, partial.dual_error, offset);
    partial.complementarity +=
        __shfl_down_sync(0xffffffffu, partial.complementarity, offset);
    partial.invalid |= __shfl_down_sync(0xffffffffu, partial.invalid, offset);
  }
  __shared__ KktScreenPartial shared[kScreenWarps];
  if ((threadIdx.x & 31) == 0) shared[threadIdx.x >> 5] = partial;
  __syncthreads();
  if (threadIdx.x != 0) return;
  KktScreenPartial total = shared[0];
  for (int warp = 1; warp < kScreenWarps; ++warp)
    mergeScreenPartial(total, shared[warp]);
  partials[blockIdx.x] = total;
}

// One thread per row segment: partial activity of the original primal.
__global__ void kktScreenSegmentsKernel(const DeviceKktScreenProblem p,
                                        const double* primal) {
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt segment = static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       segment < p.num_segments; segment += stride) {
    const HighsInt row_end = p.row_starts[p.segment_row[segment] + 1];
    const HighsInt begin = p.segment_begin[segment];
    const HighsInt end = begin + kKktScreenSegmentLength < row_end
                             ? begin + kKktScreenSegmentLength : row_end;
    double sum = 0.0;
    double absolute = 0.0;
    for (HighsInt entry = begin; entry < end; ++entry) {
      const HighsInt col = p.column_indices[entry];
      const double term = p.matrix_values[entry] * (primal[col] * p.col_scale[col]);
      sum += term;
      absolute += fabs(term);
    }
    p.segment_sum[segment] = sum;
    p.segment_absolute[segment] = absolute;
  }
}

// One thread per original row: its activity from the segment sums, possibly
// active bounds, dual-objective and complementarity terms.
__global__ void kktScreenRowsKernel(const DeviceKktScreenProblem p,
                                    const double* primal, const double* dual,
                                    double* activity, double* activity_error,
                                    KktScreenPartial* partials) {
  KktScreenPartial acc = {};
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt row = static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       row < p.num_rows; row += stride) {
    const HighsInt begin = p.row_starts[row];
    const HighsInt end = p.row_starts[row + 1];
    double sum = 0.0;
    double absolute = 0.0;
    for (HighsInt segment = p.row_first_segment[row];
         segment < p.row_first_segment[row + 1]; ++segment) {
      sum += p.segment_sum[segment];
      absolute += p.segment_absolute[segment];
    }
    // Both summation orders, the products and the unscaled primal.
    const double error = (2.0 * static_cast<double>(end - begin) + 8.0) *
                         DBL_EPSILON * absolute * 1.01;
    activity[row] = sum;
    activity_error[row] = error;
    const HighsInt transformed = p.row_map[row];
    const double unscaled = dual[transformed] * p.row_scale[transformed];
    if (!isfinite(sum) || !isfinite(error) || !isfinite(unscaled)) {
      acc.invalid = 1;
      continue;
    }
    const double lower = p.row_lower[row];
    const double upper = p.row_upper[row];
    const bool has_lower = isfinite(lower);
    const bool has_upper = isfinite(upper);
    if (has_lower && fabs(sum - lower) <= p.near_tolerance + error)
      acc.active_bound = fmax(acc.active_bound, fabs(lower));
    if (has_upper && fabs(sum - upper) <= p.near_tolerance + error)
      acc.active_bound = fmax(acc.active_bound, fabs(upper));
    // The verifier's row dual and its sense-adjusted value.
    const double adjusted = p.row_sign[row] * unscaled;
    const double row_dual = p.sense * adjusted;
    const double dual_error = 4.0 * DBL_EPSILON * fabs(unscaled);
    if (has_lower || has_upper) {
      const double bound = !has_upper   ? lower
                           : !has_lower ? upper
                           : adjusted >= 0.0 ? lower : upper;
      const double term = bound * row_dual;
      acc.dual_objective += term;
      acc.dual_absolute += fabs(term);
      acc.dual_error += fabs(bound) * dual_error;
    }
    const double magnitude = fmax(0.0, fabs(adjusted) - dual_error);
    if (adjusted > 0.0)
      acc.complementarity +=
          has_lower ? magnitude * fmax(0.0, fabs(sum - lower) - error) : magnitude;
    else if (adjusted < 0.0)
      acc.complementarity +=
          has_upper ? magnitude * fmax(0.0, fabs(upper - sum) - error) : magnitude;
  }
  writeScreenBlock(acc, partials);
}

// One thread per original column, reading the scaled transpose: the reduced
// cost as postsolve forms it, plus the column's objective, possibly active
// bounds and complementarity terms.
__global__ void kktScreenColumnsKernel(const DeviceKktScreenProblem p,
                                       const double* primal, const double* dual,
                                       double* reduced_cost,
                                       double* reduced_cost_error,
                                       KktScreenPartial* partials) {
  KktScreenPartial acc = {};
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt col = static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       col < p.num_cols; col += stride) {
    const HighsInt begin = p.transpose_starts[col];
    const HighsInt end = p.transpose_starts[col + 1];
    double sum = 0.0;
    double absolute = 0.0;
    for (HighsInt entry = begin; entry < end; ++entry) {
      const double term = p.transpose_values[entry] * dual[p.transpose_indices[entry]];
      sum += term;
      absolute += fabs(term);
    }
    const double scale = p.col_scale[col];
    const double cost = p.costs[col];
    const double rc = p.sense * ((p.scaled_costs[col] - sum) / scale);
    // Summation orders, scaled data that carries the scaling pipeline's
    // roundings, and the verifier's own original-space products.
    const double rc_error = 8.0 * (static_cast<double>(end - begin) + 64.0) *
                            DBL_EPSILON * (fabs(cost) + absolute / scale) * 1.01;
    reduced_cost[col] = rc;
    reduced_cost_error[col] = rc_error;
    const double value = primal[col] * scale;
    if (!isfinite(rc) || !isfinite(rc_error) || !isfinite(value)) {
      acc.invalid = 1;
      continue;
    }
    const double value_error = 4.0 * DBL_EPSILON * fabs(value);
    const double lower = p.col_lower[col];
    const double upper = p.col_upper[col];
    const bool has_lower = isfinite(lower);
    const bool has_upper = isfinite(upper);
    if (has_lower && fabs(value - lower) <= p.near_tolerance + value_error)
      acc.active_bound = fmax(acc.active_bound, fabs(lower));
    if (has_upper && fabs(value - upper) <= p.near_tolerance + value_error)
      acc.active_bound = fmax(acc.active_bound, fabs(upper));
    if (fabs(rc) < p.near_tolerance + rc_error)
      acc.near_zero_cost = fmax(acc.near_zero_cost, fabs(cost));
    const double objective = cost * value;
    acc.primal_objective += objective;
    acc.primal_absolute += fabs(objective);
    const double adjusted = p.sense * rc;
    if (has_lower || has_upper) {
      double term = 0.0;
      double term_error = 0.0;
      if (!has_upper) {
        term = lower * rc;
        term_error = fabs(lower) * rc_error;
      } else if (!has_lower) {
        term = upper * rc;
        term_error = fabs(upper) * rc_error;
      } else if (fabs(adjusted) <= rc_error) {
        // Either bound may be chosen, and the reduced cost is tiny.
        term_error = 2.0 * fmax(fabs(lower), fabs(upper)) * rc_error;
      } else {
        const double bound = adjusted >= 0.0 ? lower : upper;
        term = bound * rc;
        term_error = fabs(bound) * rc_error;
      }
      acc.dual_objective += term;
      acc.dual_absolute += fabs(term);
      acc.dual_error += term_error;
    }
    if (adjusted > rc_error)
      acc.complementarity +=
          has_lower ? (adjusted - rc_error) * fmax(0.0, fabs(value - lower) - value_error)
                    : adjusted - rc_error;
    else if (adjusted < -rc_error)
      acc.complementarity +=
          has_upper ? (-adjusted - rc_error) * fmax(0.0, fabs(upper - value) - value_error)
                    : -adjusted - rc_error;
  }
  writeScreenBlock(acc, partials);
}

__global__ void kktScreenTotalsKernel(const int row_count,
                                      const KktScreenPartial* row_partials,
                                      const int column_count,
                                      const KktScreenPartial* column_partials,
                                      KktScreenPacket* packet) {
  KktScreenPartial total = {};
  for (int index = threadIdx.x; index < row_count; index += blockDim.x)
    mergeScreenPartial(total, row_partials[index]);
  for (int index = threadIdx.x; index < column_count; index += blockDim.x)
    mergeScreenPartial(total, column_partials[index]);
  __shared__ KktScreenPartial shared[kBlockSize];
  shared[threadIdx.x] = total;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance)
      mergeScreenPartial(shared[threadIdx.x], shared[threadIdx.x + distance]);
    __syncthreads();
  }
  if (threadIdx.x != 0) return;
  packet->totals = shared[0];
  packet->invalid = shared[0].invalid;
}

// Lower bounds of the relative primal and dual violations, given the
// largest possible normalizations from the totals.
__global__ void kktScreenViolationsKernel(
    const DeviceKktScreenProblem p, const double* primal, const double* dual,
    const double* activity, const double* activity_error,
    const double* reduced_cost, const double* reduced_cost_error,
    const KktScreenPacket* packet, double* partials) {
  const double active_bound = packet->totals.active_bound;
  const double near_zero_cost = packet->totals.near_zero_cost;
  double primal_bound = 0.0;
  double dual_bound = 0.0;
  const HighsInt size = p.num_rows > p.num_cols ? p.num_rows : p.num_cols;
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index = static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    if (index < p.num_rows) {
      const double value = activity[index];
      const double error = activity_error[index];
      const double lower = p.row_lower[index];
      const double upper = p.row_upper[index];
      const bool has_lower = isfinite(lower);
      const bool has_upper = isfinite(upper);
      if (has_lower && lower - value - error > 0.0)
        primal_bound = fmax(primal_bound, (lower - value - error) /
                                              (1.0 + fmax(active_bound, fabs(lower))));
      if (has_upper && value - upper - error > 0.0)
        primal_bound = fmax(primal_bound, (value - upper - error) /
                                              (1.0 + fmax(active_bound, fabs(upper))));
      const HighsInt transformed = p.row_map[index];
      const double unscaled = dual[transformed] * p.row_scale[transformed];
      const double adjusted = p.row_sign[index] * unscaled;
      double violation = 0.0;
      if (!has_lower && !has_upper) {
        violation = fabs(adjusted);
      } else if (!has_upper) {
        violation = -adjusted;
      } else if (!has_lower) {
        violation = adjusted;
      } else if (upper - lower > p.near_tolerance) {
        const double midpoint = 0.5 * lower + 0.5 * upper;
        if (fabs(value - midpoint) > error)
          violation = value < midpoint ? -adjusted : adjusted;
      }
      violation -= 4.0 * DBL_EPSILON * fabs(unscaled);
      if (violation > 0.0)
        dual_bound = fmax(dual_bound, violation / (1.0 + near_zero_cost));
    }
    if (index < p.num_cols) {
      const double value = primal[index] * p.col_scale[index];
      const double value_error = 4.0 * DBL_EPSILON * fabs(value);
      const double lower = p.col_lower[index];
      const double upper = p.col_upper[index];
      const bool has_lower = isfinite(lower);
      const bool has_upper = isfinite(upper);
      if (has_lower && lower - value - value_error > 0.0)
        primal_bound = fmax(primal_bound, (lower - value - value_error) /
                                              (1.0 + fmax(active_bound, fabs(lower))));
      if (has_upper && value - upper - value_error > 0.0)
        primal_bound = fmax(primal_bound, (value - upper - value_error) /
                                              (1.0 + fmax(active_bound, fabs(upper))));
      const double adjusted = p.sense * reduced_cost[index];
      double violation = 0.0;
      if (!has_lower && !has_upper) {
        violation = fabs(adjusted);
      } else if (!has_upper) {
        violation = -adjusted;
      } else if (!has_lower) {
        violation = adjusted;
      } else if (upper - lower > p.near_tolerance) {
        const double midpoint = 0.5 * lower + 0.5 * upper;
        if (fabs(value - midpoint) > value_error)
          violation = value < midpoint ? -adjusted : adjusted;
      }
      violation -= reduced_cost_error[index];
      if (violation > 0.0)
        dual_bound = fmax(dual_bound,
                          violation / (1.0 + fmax(near_zero_cost, fabs(p.costs[index]))));
    }
  }
  __shared__ double shared_primal[kBlockSize];
  __shared__ double shared_dual[kBlockSize];
  shared_primal[threadIdx.x] = primal_bound;
  shared_dual[threadIdx.x] = dual_bound;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance) {
      shared_primal[threadIdx.x] =
          fmax(shared_primal[threadIdx.x], shared_primal[threadIdx.x + distance]);
      shared_dual[threadIdx.x] =
          fmax(shared_dual[threadIdx.x], shared_dual[threadIdx.x + distance]);
    }
    __syncthreads();
  }
  if (threadIdx.x != 0) return;
  partials[2 * blockIdx.x] = shared_primal[0];
  partials[2 * blockIdx.x + 1] = shared_dual[0];
}

__global__ void kktScreenMaximaKernel(const int count, const double* partials,
                                      KktScreenPacket* packet) {
  double primal_bound = 0.0;
  double dual_bound = 0.0;
  for (int index = threadIdx.x; index < count; index += blockDim.x) {
    primal_bound = fmax(primal_bound, partials[2 * index]);
    dual_bound = fmax(dual_bound, partials[2 * index + 1]);
  }
  __shared__ double shared_primal[kBlockSize];
  __shared__ double shared_dual[kBlockSize];
  shared_primal[threadIdx.x] = primal_bound;
  shared_dual[threadIdx.x] = dual_bound;
  __syncthreads();
  for (int distance = kBlockSize / 2; distance > 0; distance /= 2) {
    if (threadIdx.x < distance) {
      shared_primal[threadIdx.x] =
          fmax(shared_primal[threadIdx.x], shared_primal[threadIdx.x + distance]);
      shared_dual[threadIdx.x] =
          fmax(shared_dual[threadIdx.x], shared_dual[threadIdx.x + distance]);
    }
    __syncthreads();
  }
  if (threadIdx.x != 0) return;
  packet->primal = shared_primal[0];
  packet->dual = shared_dual[0];
}

}  // namespace

cudaError_t preloadKernels() {
  const void* const kernels[] = {
      reinterpret_cast<const void*>(fillKernel),
      reinterpret_cast<const void*>(scaleKernel),
      reinterpret_cast<const void*>(axpyKernel),
      reinterpret_cast<const void*>(differenceKernel),
      reinterpret_cast<const void*>(projectedPrimalKernel<false>),
      reinterpret_cast<const void*>(projectedPrimalKernel<true>),
      reinterpret_cast<const void*>(projectedDualKernel<false>),
      reinterpret_cast<const void*>(projectedDualKernel<true>),
      reinterpret_cast<const void*>(vectorBlendKernel),
      reinterpret_cast<const void*>(vectorPairBlendKernel),
      reinterpret_cast<const void*>(averageKernel),
      reinterpret_cast<const void*>(finiteCheckKernel),
      reinterpret_cast<const void*>(dotKernel),
      reinterpret_cast<const void*>(normKernel),
      reinterpret_cast<const void*>(differenceNormKernel),
      reinterpret_cast<const void*>(fixedPointPartialsKernel),
      reinterpret_cast<const void*>(fixedPointResultKernel),
      reinterpret_cast<const void*>(checkPartialsKernel),
      reinterpret_cast<const void*>(checkPacketKernel),
      reinterpret_cast<const void*>(certificateRecoveryKernel),
      reinterpret_cast<const void*>(certificateNormKernel),
      reinterpret_cast<const void*>(certificateNormalizeKernel),
      reinterpret_cast<const void*>(certificateRowActivityKernel),
      reinterpret_cast<const void*>(certificateLongRowKernel),
      reinterpret_cast<const void*>(certificateConditionsKernel),
      reinterpret_cast<const void*>(certificateRaySkipKernel),
      reinterpret_cast<const void*>(certificateRowViolationKernel),
      reinterpret_cast<const void*>(certificateFlagsKernel),
#ifdef HIGHS_HIPDLP_TESTING
      reinterpret_cast<const void*>(injectedFailureKernel),
#endif
  };
  for (const void* kernel : kernels) {
    cudaFuncAttributes attributes;
    const cudaError_t status = cudaFuncGetAttributes(&attributes, kernel);
    if (status != cudaSuccess) return status;
  }
  return cudaSuccess;
}

std::size_t checkPartialCount(const HighsInt rows, const HighsInt cols) {
  const HighsInt size = rows > cols ? rows : cols;
  if (size == 0) return 1;
  const HighsInt count = 1 + (size - 1) / kBlockSize;
  return static_cast<std::size_t>(count < kMaximumCheckBlocks
                                      ? count : kMaximumCheckBlocks);
}

static int kktScreenBlocks(const HighsInt size) {
  const HighsInt blocks = size <= 0 ? 1 : 1 + (size - 1) / kBlockSize;
  return static_cast<int>(blocks < kMaximumScreenBlocks ? blocks
                                                        : kMaximumScreenBlocks);
}

std::size_t kktScreenPartialCount(const HighsInt, const HighsInt) {
  return static_cast<std::size_t>(kMaximumScreenBlocks);
}

cudaError_t launchKktScreen(const cudaStream_t stream,
                            const DeviceKktScreenProblem& problem,
                            const double* primal, const double* dual,
                            double* activity, double* activity_error,
                            double* reduced_cost, double* reduced_cost_error,
                            KktScreenPartial* row_partials,
                            KktScreenPartial* column_partials,
                            double* violation_partials, KktScreenPacket* packet) {
  const int row_blocks = kktScreenBlocks(problem.num_rows);
  const int column_blocks = kktScreenBlocks(problem.num_cols);
  kktScreenSegmentsKernel<<<kktScreenBlocks(problem.num_segments), kBlockSize, 0,
                            stream>>>(problem, primal);
  cudaError_t segment_status = cudaPeekAtLastError();
  if (segment_status != cudaSuccess) return segment_status;
  kktScreenRowsKernel<<<row_blocks, kBlockSize, 0, stream>>>(
      problem, primal, dual, activity, activity_error, row_partials);
  cudaError_t status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  kktScreenColumnsKernel<<<column_blocks, kBlockSize, 0, stream>>>(
      problem, primal, dual, reduced_cost, reduced_cost_error, column_partials);
  status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  kktScreenTotalsKernel<<<1, kBlockSize, 0, stream>>>(
      row_blocks, row_partials, column_blocks, column_partials, packet);
  status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  const int violation_blocks =
      static_cast<int>(checkPartialCount(problem.num_rows, problem.num_cols));
  kktScreenViolationsKernel<<<violation_blocks, kBlockSize, 0, stream>>>(
      problem, primal, dual, activity, activity_error, reduced_cost,
      reduced_cost_error, packet, violation_partials);
  status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  kktScreenMaximaKernel<<<1, kBlockSize, 0, stream>>>(violation_blocks,
                                                      violation_partials, packet);
  return cudaPeekAtLastError();
}

cudaError_t launchFixedPointReduction(
    const cudaStream_t stream, const HighsInt primal_size, const HighsInt dual_size,
    const double* primal, const double* dual, const double* product,
    FixedPointPartial* partials, double* result, int* error_flag) {
  const int blocks = static_cast<int>(checkPartialCount(primal_size, dual_size));
  fixedPointPartialsKernel<<<blocks, kBlockSize, 0, stream>>>(
      primal_size, dual_size, primal, dual, product, partials);
  const cudaError_t status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  fixedPointResultKernel<<<1, kBlockSize, 0, stream>>>(
      blocks, partials, result, error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchCheckPartials(const cudaStream_t stream,
                                const DeviceCheckProblem& problem,
                                const DeviceCheckInput& input,
                                CheckPartial* partials) {
  const int blocks = static_cast<int>(checkPartialCount(problem.num_rows, problem.num_cols));
  checkPartialsKernel<<<blocks, kBlockSize, 0, stream>>>(problem, input, partials);
  return cudaPeekAtLastError();
}

cudaError_t launchCheckPacket(const cudaStream_t stream,
                              const DeviceCheckProblem& problem,
                              const DeviceCheckInput& input,
                              const CheckPartial* partials, CheckPacket* packet,
                              int* error_flag) {
  const int count = static_cast<int>(checkPartialCount(problem.num_rows, problem.num_cols));
  checkPacketKernel<<<1, kBlockSize, 0, stream>>>(problem, input, count, partials,
                                               packet, error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchCertificateScreen(
    const cudaStream_t stream, const DeviceCertificateProblem& problem,
    const DeviceCheckInput& input, certificate::Interval* primal,
    certificate::Interval* dual, certificate::Interval* column,
    certificate::Interval* activity, CertificateNormPartial* norm_partials,
    CertificateState* state, certificate::ScreenPartial* partials,
    CheckPacket* packet) {
  const int blocks =
      static_cast<int>(checkPartialCount(problem.num_rows, problem.num_cols));
  certificateRecoveryKernel<<<blocks, kBlockSize, 0, stream>>>(
      problem, input, primal, dual, column, norm_partials);
  cudaError_t status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  certificateNormKernel<<<1, kBlockSize, 0, stream>>>(
      problem, blocks, norm_partials, false, state);
  status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  certificateNormalizeKernel<<<blocks, kBlockSize, 0, stream>>>(
      problem, state, primal, dual, column, norm_partials);
  status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  certificateNormKernel<<<1, kBlockSize, 0, stream>>>(
      problem, blocks, norm_partials, true, state);
  status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  certificateConditionsKernel<<<blocks, kBlockSize, 0, stream>>>(
      problem, state, primal, dual, column, partials);
  status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  certificateRaySkipKernel<<<1, kBlockSize, 0, stream>>>(problem, input, blocks,
                                                         partials, state);
  status = cudaPeekAtLastError();
  if (status != cudaSuccess) return status;
  const HighsInt segments_per_block = kBlockSize / kWarpSize;
  const HighsInt activity_blocks =
      problem.num_segments == 0
          ? 0 : 1 + (problem.num_segments - 1) / segments_per_block;
  if (activity_blocks > 0) {
    certificateRowActivityKernel<<<
        static_cast<int>(activity_blocks < kMaximumBlocks ? activity_blocks
                                                          : kMaximumBlocks),
        kBlockSize, 0, stream>>>(problem, state, primal, activity);
    status = cudaPeekAtLastError();
    if (status != cudaSuccess) return status;
    if (problem.num_long_rows > 0) {
      certificateLongRowKernel<<<blockCount(problem.num_long_rows), kBlockSize,
                                 0, stream>>>(problem, state, activity);
      status = cudaPeekAtLastError();
      if (status != cudaSuccess) return status;
    }
    certificateRowViolationKernel<<<blocks, kBlockSize, 0, stream>>>(
        problem, state, activity, state);
    status = cudaPeekAtLastError();
    if (status != cudaSuccess) return status;
  }
  certificateFlagsKernel<<<1, kBlockSize, 0, stream>>>(problem, input, blocks,
                                                       state, partials, packet);
  return cudaPeekAtLastError();
}

cudaError_t launchFill(const cudaStream_t stream, const HighsInt size,
                       const double value, double* output) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  fillKernel<<<blocks, kBlockSize, 0, stream>>>(size, value, output);
  return cudaPeekAtLastError();
}

cudaError_t launchScale(const cudaStream_t stream, const HighsInt size,
                        const double factor, double* values, int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  scaleKernel<<<blocks, kBlockSize, 0, stream>>>(size, factor, values,
                                                 error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchAxpy(const cudaStream_t stream, const HighsInt size,
                       const double factor, const double* input, double* output,
                       int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  axpyKernel<<<blocks, kBlockSize, 0, stream>>>(size, factor, input, output,
                                                error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchDifference(const cudaStream_t stream, const HighsInt size,
                             const double* left, const double* right,
                             double* output, int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  differenceKernel<<<blocks, kBlockSize, 0, stream>>>(size, left, right, output,
                                                    error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchProjectedPrimalUpdate(
    const cudaStream_t stream, const HighsInt size, const double* current,
    const double* transpose_product, const double* costs, const double* lower,
    const double* upper, const double step_size, double* output,
    int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
    projectedPrimalKernel<false><<<blocks, kBlockSize, 0, stream>>>(
      size, current, transpose_product, costs, lower, upper, step_size, output,
      nullptr, error_flag, nullptr);
  return cudaPeekAtLastError();
}

cudaError_t launchProjectedDualUpdate(
    const cudaStream_t stream, const HighsInt size, const double* current,
    const double* extrapolated_product, const double* row_lower,
    const HighsBool* equality_rows, const double step_size, double* output,
    int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  projectedDualKernel<false><<<blocks, kBlockSize, 0, stream>>>(
      size, current, extrapolated_product, row_lower, equality_rows, step_size,
      output, nullptr, error_flag, nullptr);
  return cudaPeekAtLastError();
}

cudaError_t launchProjectedPrimalUpdateAndReflect(
    const cudaStream_t stream, const HighsInt size, const double* current,
    const double* transpose_product, const double* costs, const double* lower,
    const double* upper, const double step_size, double* output,
    double* reflected, int* error_flag, const double* step_parameter) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  projectedPrimalKernel<true><<<blocks, kBlockSize, 0, stream>>>(
      size, current, transpose_product, costs, lower, upper, step_size, output,
      reflected, error_flag, step_parameter);
  return cudaPeekAtLastError();
}

cudaError_t launchProjectedDualUpdateAndReflect(
    const cudaStream_t stream, const HighsInt size, const double* current,
    const double* extrapolated_product, const double* row_lower,
    const HighsBool* equality_rows, const double step_size, double* output,
    double* reflected, int* error_flag, const double* step_parameter) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  projectedDualKernel<true><<<blocks, kBlockSize, 0, stream>>>(
      size, current, extrapolated_product, row_lower, equality_rows, step_size,
      output, reflected, error_flag, step_parameter);
  return cudaPeekAtLastError();
}

cudaError_t launchVectorBlend(const cudaStream_t stream, const HighsInt size,
                              const double* left, const double* right,
                              const double left_weight, double* output,
                              int* error_flag, const double* weight_parameter) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  vectorBlendKernel<<<blocks, kBlockSize, 0, stream>>>(
      size, left, right, left_weight, output, error_flag, weight_parameter);
  return cudaPeekAtLastError();
}

cudaError_t launchVectorPairBlend(const cudaStream_t stream, const HighsInt primal_size,
    const HighsInt dual_size, const double* primal_left, const double* primal_right,
    const double* dual_left, const double* dual_right, const double left_weight,
    double* primal_output, double* dual_output, int* error_flag,
    const double* weight_parameter) {
  const int blocks = blockCount(primal_size > dual_size ? primal_size : dual_size);
  if (blocks == 0) return cudaSuccess;
  vectorPairBlendKernel<<<blocks, kBlockSize, 0, stream>>>(primal_size, dual_size,
      primal_left, primal_right, dual_left, dual_right, left_weight,
      primal_output, dual_output, error_flag, weight_parameter);
  return cudaPeekAtLastError();
}

cudaError_t launchAverageUpdate(const cudaStream_t stream, const HighsInt size,
                                const double* current, const double weight,
                                double* sum, int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  averageKernel<<<blocks, kBlockSize, 0, stream>>>(size, current, weight, sum,
                                                   error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchFiniteCheck(const cudaStream_t stream, const HighsInt size,
                              const double* values, int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  finiteCheckKernel<<<blocks, kBlockSize, 0, stream>>>(size, values,
                                                       error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchDot(const cudaStream_t stream, const HighsInt size,
                      const double* left, const double* right, double* result,
                      int* error_flag) {
  dotKernel<<<1, 1, 0, stream>>>(size, left, right, result, error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchNorm2(const cudaStream_t stream, const HighsInt size,
                        const double* values, double* result, int* error_flag) {
  normKernel<<<1, 1, 0, stream>>>(size, values, result, error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchDifferenceNorm2(const cudaStream_t stream,
                                  const HighsInt size, const double* left,
                                  const double* right, double* result,
                                  int* error_flag) {
  differenceNormKernel<<<1, 1, 0, stream>>>(size, left, right, result,
                                            error_flag);
  return cudaPeekAtLastError();
}

#ifdef HIGHS_HIPDLP_TESTING
cudaError_t launchInjectedFailure(const cudaStream_t stream, int* error_flag) {
  injectedFailureKernel<<<1, 1, 0, stream>>>(error_flag);
  return cudaPeekAtLastError();
}
#endif

}  // namespace cuda
}  // namespace backend
}  // namespace hipdlp