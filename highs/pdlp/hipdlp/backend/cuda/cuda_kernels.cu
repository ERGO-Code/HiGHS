#include <cmath>

#include "pdlp/hipdlp/backend/cuda/cuda_kernels.cuh"

namespace hipdlp {
namespace backend {
namespace cuda {

namespace {

const int kBlockSize = 256;
const int kMaximumBlocks = 65535;

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

__global__ void projectedPrimalKernel(
    const HighsInt size, const double* current, const double* transpose_product,
    const double* costs, const double* lower, const double* upper,
    const double step_size, double* output, int* error_flag) {
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
  }
}

__global__ void projectedDualKernel(const HighsInt size, const double* current,
                                    const double* extrapolated_product,
                                    const double* row_lower,
                                    const HighsBool* equality_rows,
                                    const double step_size, double* output,
                                    int* error_flag) {
  const HighsInt stride = static_cast<HighsInt>(blockDim.x) * gridDim.x;
  for (HighsInt index =
           static_cast<HighsInt>(blockIdx.x) * blockDim.x + threadIdx.x;
       index < size; index += stride) {
    double value = current[index] +
                   step_size * (row_lower[index] - extrapolated_product[index]);
    if (!equality_rows[index] && value < 0.0) value = 0.0;
    output[index] = value;
    recordNonFinite(value, error_flag);
  }
}

__global__ void vectorBlendKernel(const HighsInt size, const double* left,
                                  const double* right, const double left_weight,
                                  double* output, int* error_flag) {
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

#ifdef HIGHS_HIPDLP_TESTING
__global__ void injectedFailureKernel(int* error_flag) {
  if (blockIdx.x != 0 || threadIdx.x != 0) return;
  atomicExch(error_flag, 1);
}
#endif

}  // namespace

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

cudaError_t launchProjectedPrimalUpdate(
    const cudaStream_t stream, const HighsInt size, const double* current,
    const double* transpose_product, const double* costs, const double* lower,
    const double* upper, const double step_size, double* output,
    int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  projectedPrimalKernel<<<blocks, kBlockSize, 0, stream>>>(
      size, current, transpose_product, costs, lower, upper, step_size, output,
      error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchProjectedDualUpdate(
    const cudaStream_t stream, const HighsInt size, const double* current,
    const double* extrapolated_product, const double* row_lower,
    const HighsBool* equality_rows, const double step_size, double* output,
    int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  projectedDualKernel<<<blocks, kBlockSize, 0, stream>>>(
      size, current, extrapolated_product, row_lower, equality_rows, step_size,
      output, error_flag);
  return cudaPeekAtLastError();
}

cudaError_t launchVectorBlend(const cudaStream_t stream, const HighsInt size,
                              const double* left, const double* right,
                              const double left_weight, double* output,
                              int* error_flag) {
  const int blocks = blockCount(size);
  if (blocks == 0) return cudaSuccess;
  vectorBlendKernel<<<blocks, kBlockSize, 0, stream>>>(
      size, left, right, left_weight, output, error_flag);
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