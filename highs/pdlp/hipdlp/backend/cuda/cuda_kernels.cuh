#ifndef PDLP_HIPDLP_BACKEND_CUDA_CUDA_KERNELS_CUH_
#define PDLP_HIPDLP_BACKEND_CUDA_CUDA_KERNELS_CUH_

#include <cuda_runtime_api.h>

#include "util/HighsInt.h"
#include "util/HighsType.h"

namespace hipdlp {
namespace backend {
namespace cuda {

cudaError_t launchFill(cudaStream_t stream, HighsInt size, double value,
                       double* output);
cudaError_t launchScale(cudaStream_t stream, HighsInt size, double factor,
                        double* values, int* error_flag);
cudaError_t launchAxpy(cudaStream_t stream, HighsInt size, double factor,
                       const double* input, double* output, int* error_flag);
cudaError_t launchProjectedPrimalUpdate(
    cudaStream_t stream, HighsInt size, const double* current,
    const double* transpose_product, const double* costs, const double* lower,
    const double* upper, double step_size, double* output, int* error_flag);
cudaError_t launchProjectedDualUpdate(cudaStream_t stream, HighsInt size,
                                      const double* current,
                                      const double* extrapolated_product,
                                      const double* row_lower,
                                      const HighsBool* equality_rows,
                                      double step_size, double* output,
                                      int* error_flag);
cudaError_t launchVectorBlend(cudaStream_t stream, HighsInt size,
                              const double* left, const double* right,
                              double left_weight, double* output,
                              int* error_flag);
cudaError_t launchAverageUpdate(cudaStream_t stream, HighsInt size,
                                const double* current, double weight,
                                double* sum, int* error_flag);
cudaError_t launchFiniteCheck(cudaStream_t stream, HighsInt size,
                              const double* values, int* error_flag);
cudaError_t launchDot(cudaStream_t stream, HighsInt size, const double* left,
                      const double* right, double* result, int* error_flag);
cudaError_t launchNorm2(cudaStream_t stream, HighsInt size,
                        const double* values, double* result, int* error_flag);
cudaError_t launchDifferenceNorm2(cudaStream_t stream, HighsInt size,
                                  const double* left, const double* right,
                                  double* result, int* error_flag);
#ifdef HIGHS_HIPDLP_TESTING
cudaError_t launchInjectedFailure(cudaStream_t stream, int* error_flag);
#endif

}  // namespace cuda
}  // namespace backend
}  // namespace hipdlp

#endif