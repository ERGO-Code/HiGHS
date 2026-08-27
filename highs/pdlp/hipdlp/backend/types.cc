#include "pdlp/hipdlp/backend/types.hpp"

namespace hipdlp {
namespace backend {

bool isValidKind(const Kind kind) {
  switch (kind) {
    case Kind::kAuto:
    case Kind::kCpu:
    case Kind::kCuda:
    case Kind::kOpenmp:
    case Kind::kRocm:
      return true;
  }
  return false;
}

PrecisionDescriptor::PrecisionDescriptor()
    : storage(Precision::kFloat64),
      compute(Precision::kFloat64),
      reduction(Precision::kFloat64) {}

PrecisionDescriptor::PrecisionDescriptor(const Precision storage,
                                         const Precision compute,
                                         const Precision reduction)
    : storage(storage), compute(compute), reduction(reduction) {}

bool PrecisionDescriptor::valid() const {
  const bool storage_valid =
      storage == Precision::kFloat32 || storage == Precision::kFloat64;
  const bool compute_valid =
      compute == Precision::kFloat32 || compute == Precision::kFloat64;
  const bool reduction_valid =
      reduction == Precision::kFloat32 || reduction == Precision::kFloat64;
  return storage_valid && compute_valid && reduction_valid;
}

bool PrecisionDescriptor::isFp64() const {
  return storage == Precision::kFloat64 && compute == Precision::kFloat64 &&
         reduction == Precision::kFloat64;
}

Capabilities::Capabilities() : mask(0) {}

Capabilities::Capabilities(const unsigned long long mask) : mask(mask) {}

bool Capabilities::supports(const Capability capability) const {
  return (mask & static_cast<unsigned long long>(capability)) != 0;
}

DeviceInfo::DeviceInfo() : available(false) {}

ConvergenceMetrics::ConvergenceMetrics()
    : primal_residual(0.0), dual_residual(0.0), objective_gap(0.0) {}

}  // namespace backend
}  // namespace hipdlp