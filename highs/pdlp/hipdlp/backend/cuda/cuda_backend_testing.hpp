#ifndef PDLP_HIPDLP_BACKEND_CUDA_CUDA_BACKEND_TESTING_HPP_
#define PDLP_HIPDLP_BACKEND_CUDA_CUDA_BACKEND_TESTING_HPP_

#ifndef HIGHS_HIPDLP_TESTING
#error "HiPDLP CUDA test hooks require HIGHS_HIPDLP_TESTING"
#endif

#include <cstddef>

namespace hipdlp {
namespace backend {
namespace cuda {
namespace testing {

enum class FailurePoint {
  kNone,
  kBackendInitialization,
  kDeviceAllocation,
  kGraphCapture,
  kGraphInstantiation,
  kGraphLaunch,
  kGraphAsynchronousExecution
};

struct ResourceAudit {
  std::size_t stream_creations;
  std::size_t stream_destructions;
  std::size_t cublas_creations;
  std::size_t cublas_destructions;
  std::size_t cusparse_creations;
  std::size_t cusparse_destructions;
  std::size_t allocation_calls;
  std::size_t free_calls;
  std::size_t descriptor_creations;
  std::size_t descriptor_destructions;
  std::size_t host_to_device_transfers;
  std::size_t device_to_host_transfers;
  std::size_t device_to_device_transfers;
  std::size_t graph_capture_calls;
  std::size_t graph_creations;
  std::size_t graph_destructions;
  std::size_t graph_instantiation_calls;
  std::size_t graph_exec_creations;
  std::size_t graph_exec_destructions;
  std::size_t graph_launch_calls;
  std::size_t graph_update_calls;
  std::size_t graph_capture_primitive_enqueue_calls;
  std::size_t direct_primitive_enqueue_calls;
  std::size_t synchronization_calls;

  ResourceAudit();
};

void setFailurePoint(FailurePoint failure_point);
void resetResourceAudit();
ResourceAudit resourceAudit();

}  // namespace testing
}  // namespace cuda
}  // namespace backend
}  // namespace hipdlp

#endif