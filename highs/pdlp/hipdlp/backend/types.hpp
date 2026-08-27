#ifndef PDLP_HIPDLP_BACKEND_TYPES_HPP_
#define PDLP_HIPDLP_BACKEND_TYPES_HPP_

#include <memory>
#include <string>
#include <utility>

#include "util/HighsInt.h"

namespace hipdlp {
namespace backend {

namespace cpu {
class CpuBackend;
}
namespace cuda {
class CudaBackend;
}

enum class Kind { kAuto, kCpu, kCuda, kOpenmp, kRocm };

bool isValidKind(Kind kind);

enum class Precision { kFloat32, kFloat64 };

struct PrecisionDescriptor {
  Precision storage;
  Precision compute;
  Precision reduction;

  PrecisionDescriptor();
  PrecisionDescriptor(Precision storage, Precision compute,
                      Precision reduction);

  bool valid() const;
  bool isFp64() const;
};

enum class Capability : unsigned long long {
  kSparseMatrixVector = 1ULL << 0,
  kFusedPrimalDualUpdate = 1ULL << 1,
  kFusedAverageUpdate = 1ULL << 2,
  kConvergenceReduction = 1ULL << 3,
  kBatchedVectors = 1ULL << 4,
  kSparseMatrixMatrix = 1ULL << 5,
  kVectorOperations = 1ULL << 6,
  kReductions = 1ULL << 7,
  kHostDownload = 1ULL << 8,
  kCommandBuffer = 1ULL << 9
};

struct Capabilities {
  unsigned long long mask;

  Capabilities();
  explicit Capabilities(unsigned long long mask);

  bool supports(Capability capability) const;
};

struct Resource {
  virtual ~Resource() {}
};

struct DeviceInfo {
  bool available;
  std::string name;

  DeviceInfo();
};

template <typename Tag>
class Handle {
 public:
  Handle() {}
  explicit Handle(std::unique_ptr<Resource> resource)
      : resource_(std::move(resource)) {}

  Handle(Handle&& other) noexcept : resource_(std::move(other.resource_)) {}
  Handle& operator=(Handle&& other) noexcept {
    resource_ = std::move(other.resource_);
    return *this;
  }

  Handle(const Handle&) = delete;
  Handle& operator=(const Handle&) = delete;

  bool valid() const { return static_cast<bool>(resource_); }
  void reset() { resource_.reset(); }

 private:
  std::unique_ptr<Resource> resource_;

  friend class cpu::CpuBackend;
  friend class cuda::CudaBackend;
};

struct ProblemHandleTag {};
struct VectorHandleTag {};
struct WorkspaceHandleTag {};
struct CommandBufferHandleTag {};

using ProblemHandle = Handle<ProblemHandleTag>;
using VectorHandle = Handle<VectorHandleTag>;
using WorkspaceHandle = Handle<WorkspaceHandleTag>;
using CommandBufferHandle = Handle<CommandBufferHandleTag>;

struct ConvergenceMetrics {
  double primal_residual;
  double dual_residual;
  double objective_gap;

  ConvergenceMetrics();
};

}  // namespace backend
}  // namespace hipdlp

#endif