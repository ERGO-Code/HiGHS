#ifndef PDLP_HIPDLP_MODEL_PREPROCESSING_HPP_
#define PDLP_HIPDLP_MODEL_PREPROCESSING_HPP_

#include <chrono>

namespace hipdlp {
namespace model {

// Optional wall timings. Child phases are nested in model_build/scaling and
// accumulate across both candidate and fallback scaling pipelines.
struct PreprocessingDiagnostics {
  double model_build_seconds = 0.0;
  double validation_seconds = 0.0;
  double input_matrix_seconds = 0.0;
  double row_transform_seconds = 0.0;
  double transformed_matrix_seconds = 0.0;
  double transpose_seconds = 0.0;
  double unscaled_norms_seconds = 0.0;
  double scaling_seconds = 0.0;
  double scaling_candidate_seconds = 0.0;
  double scaling_reference_seconds = 0.0;
  double scaling_guard_seconds = 0.0;
  double scaling_apply_seconds = 0.0;
  double scaling_copy_seconds = 0.0;
  double curtis_reid_seconds = 0.0;
  double ruiz_seconds = 0.0;
  double pock_chambolle_seconds = 0.0;
  double l2_seconds = 0.0;
};

namespace detail {
class PreprocessingTimer {
 public:
  explicit PreprocessingTimer(double* seconds) { reset(seconds); }
  ~PreprocessingTimer() { reset(nullptr); }

  void reset(double* seconds) {
    if (seconds_ != nullptr)
      *seconds_ += std::chrono::duration<double>(Clock::now() - start_).count();
    seconds_ = seconds;
    if (seconds_ != nullptr) start_ = Clock::now();
  }

  PreprocessingTimer(const PreprocessingTimer&) = delete;
  PreprocessingTimer& operator=(const PreprocessingTimer&) = delete;

 private:
  using Clock = std::chrono::steady_clock;
  double* seconds_ = nullptr;
  Clock::time_point start_;
};
}  // namespace detail
}  // namespace model
}  // namespace hipdlp
#endif
