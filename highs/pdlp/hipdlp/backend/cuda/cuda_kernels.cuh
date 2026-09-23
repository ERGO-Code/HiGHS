#ifndef PDLP_HIPDLP_BACKEND_CUDA_CUDA_KERNELS_CUH_
#define PDLP_HIPDLP_BACKEND_CUDA_CUDA_KERNELS_CUH_

#include <cuda_runtime_api.h>

#include <cstddef>

#include "pdlp/hipdlp/backend/cuda/certificate_screen.hpp"
#include "util/HighsInt.h"
#include "util/HighsType.h"

namespace hipdlp {
namespace backend {
namespace cuda {

struct DeviceCheckProblem {
  HighsInt num_rows;
  HighsInt num_cols;
  const double* costs;
  const double* col_lower;
  const double* col_upper;
  const double* row_lower;
  const HighsBool* equality_rows;
  const double* row_scale;
  const double* col_scale;
  double objective_offset;
  double unscaled_rhs_norm;
  double unscaled_cost_norm;
};

struct DeviceCheckInput {
  const double* primal;
  const double* dual;
  const double* ax;
  const double* aty;
  const double* delta_primal;
  const double* delta_dual;
  const double* aty_delta;
  const double* primal_anchor;
  const double* dual_anchor;
  const int* optional_error;
  bool fixed_point;
  bool anchor_distances;
  double certificate_tolerance;
  double primal_tolerance;
  double dual_tolerance;
  double gap_tolerance;
};

struct DeviceCertificateProblem {
  HighsInt num_rows;
  HighsInt num_cols;
  const double* costs;
  const double* col_lower;
  const double* col_upper;
  const double* row_lower;
  const double* row_upper;
  const HighsInt* row_starts;
  const HighsInt* column_indices;
  const double* matrix_values;
  const HighsInt* row_map;
  const double* row_sign;
  const double* row_scale;
  const double* col_scale;
  double sense_sign;
  bool column_bounds_straddle_zero;
};

struct CertificateNormPartial {
  certificate::NormBounds primal;
  certificate::NormBounds dual;
  int unsafe;
};

struct CertificateState {
  certificate::Interval primal_normalization;
  certificate::Interval dual_normalization;
  certificate::Interval primal_verifier_norm;
  certificate::Interval dual_verifier_norm;
  int unsafe;
};

struct CheckNormState {
  double scale;
  double sum_of_squares;
};

struct CheckSumState {
  double value;
  double absolute_sum;
  double minimum_nonzero;
};

enum CheckFallbackFlags {
  kCheckRawFallback = 1,
  kCheckFixedPointFallback = 2,
  kCheckAnchorFallback = 4,
  kCheckUnpreparedProducts = 8
};

struct CheckPartial {
  CheckNormState norms[6];
  CheckSumState sums[3];
  int numerical_error;
  int optional_invalid;
  int requires_host_oracle;
};

struct CheckPacket {
  double primal_residual;
  double dual_residual;
  double objective_gap;
  double primal_objective;
  double dual_objective;
  double unscaled_rhs_norm;
  double unscaled_cost_norm;
  double primal_delta_norm;
  double dual_delta_norm;
  double delta_dot;
  double primal_anchor_distance;
  double dual_anchor_distance;
  int error_flag;
  int fallback_flags;
  int primal_ray_possible;
  int farkas_possible;
};

static_assert(sizeof(CheckPacket) == 112,
              "unexpected CUDA check packet layout");

cudaError_t preloadKernels();

cudaError_t launchCertificateScreen(
    cudaStream_t stream, const DeviceCertificateProblem& problem,
    const DeviceCheckInput& input, certificate::Interval* primal,
    certificate::Interval* dual, certificate::Interval* column,
    CertificateNormPartial* norm_partials, CertificateState* state,
    certificate::ScreenPartial* partials, CheckPacket* packet);

std::size_t checkPartialCount(HighsInt rows, HighsInt cols);
cudaError_t launchCheckPartials(cudaStream_t stream,
                                const DeviceCheckProblem& problem,
                                const DeviceCheckInput& input,
                                CheckPartial* partials);
cudaError_t launchCheckPacket(cudaStream_t stream,
                              const DeviceCheckProblem& problem,
                              const DeviceCheckInput& input,
                              const CheckPartial* partials, CheckPacket* packet,
                              int* error_flag);

cudaError_t launchFill(cudaStream_t stream, HighsInt size, double value,
                       double* output);
cudaError_t launchScale(cudaStream_t stream, HighsInt size, double factor,
                        double* values, int* error_flag);
cudaError_t launchAxpy(cudaStream_t stream, HighsInt size, double factor,
                       const double* input, double* output, int* error_flag);
cudaError_t launchDifference(cudaStream_t stream, HighsInt size,
                             const double* left, const double* right,
                             double* output, int* error_flag);
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
                    cudaError_t launchProjectedPrimalUpdateAndReflect(
                      cudaStream_t stream, HighsInt size, const double* current,
                      const double* transpose_product, const double* costs, const double* lower,
                      const double* upper, double step_size, double* output, double* reflected,
                      int* error_flag, const double* step_parameter = nullptr);
                    cudaError_t launchProjectedDualUpdateAndReflect(
                      cudaStream_t stream, HighsInt size, const double* current,
                      const double* extrapolated_product, const double* row_lower,
                      const HighsBool* equality_rows, double step_size, double* output,
                      double* reflected, int* error_flag, const double* step_parameter = nullptr);
cudaError_t launchVectorBlend(cudaStream_t stream, HighsInt size,
                              const double* left, const double* right,
                              double left_weight, double* output,
                              int* error_flag,
                              const double* weight_parameter = nullptr);
                cudaError_t launchVectorPairBlend(cudaStream_t stream, HighsInt primal_size,
                  HighsInt dual_size, const double* primal_left, const double* primal_right,
                  const double* dual_left, const double* dual_right, double left_weight,
                  double* primal_output, double* dual_output, int* error_flag,
                  const double* weight_parameter);
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