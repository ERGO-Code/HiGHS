#include "pdlp/hipdlp/backend/backend_factory.hpp"

#include <new>
#include <utility>

#include "HConfig.h"
#include "pdlp/hipdlp/backend/cpu/cpu_backend.hpp"
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
#include "pdlp/hipdlp/backend/cuda/cuda_backend.hpp"
#endif

namespace hipdlp {
namespace backend {

namespace {

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
const bool kCudaBuildEnabled = true;
#else
const bool kCudaBuildEnabled = false;
#endif

const Descriptor kRegistrations[] = {
    {Kind::kCpu, "cpu", true, true, PrecisionDescriptor(),
     cpu::CpuBackend().descriptor().capabilities},
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    cuda::cudaDescriptor(),
#else
    {Kind::kCuda, "cuda", kCudaBuildEnabled, false, PrecisionDescriptor(),
     Capabilities()},
#endif
    {Kind::kOpenmp, "openmp", false, false, PrecisionDescriptor(),
     Capabilities()},
    {Kind::kRocm, "rocm", false, false, PrecisionDescriptor(), Capabilities()}};

}  // namespace

FactoryResult::FactoryResult() {}

FactoryResult::FactoryResult(FactoryResult&& other)
    : error(std::move(other.error)), backend(std::move(other.backend)) {}

FactoryResult& FactoryResult::operator=(FactoryResult&& other) {
  error = std::move(other.error);
  backend = std::move(other.backend);
  return *this;
}

bool FactoryResult::ok() const {
  return error.ok() && backend.get() != nullptr;
}

const Descriptor* registrations(std::size_t& count) {
  count = sizeof(kRegistrations) / sizeof(kRegistrations[0]);
  return kRegistrations;
}

const Descriptor* findRegistration(const Kind kind) {
  std::size_t count = 0;
  const Descriptor* entries = registrations(count);
  for (std::size_t index = 0; index < count; ++index) {
    if (entries[index].kind == kind) return &entries[index];
  }
  return nullptr;
}

FactoryResult createBackend(const Kind kind) {
  FactoryResult result;
  if (!isValidKind(kind)) {
    result.error =
        api::Error(api::ErrorCode::kInvalidOption, "invalid HiPDLP backend");
    return result;
  }
  Kind selected_kind = kind;
  if (kind == Kind::kAuto) {
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    api::Error cuda_error = cuda::createCudaBackend(result.backend);
    if (cuda_error.ok()) return result;
#endif
    selected_kind = Kind::kCpu;
  }

  const Descriptor* registration = findRegistration(selected_kind);
  if (registration == nullptr || !registration->build_enabled) {
    result.error = api::Error(api::ErrorCode::kBackendUnavailable,
                              "requested HiPDLP backend is unavailable");
    return result;
  }
  if (!registration->implemented) {
    result.error = api::Error(api::ErrorCode::kNotReady,
                              "requested HiPDLP backend is not ready");
    return result;
  }

  if (selected_kind == Kind::kCpu) {
    try {
      result.backend = cpu::createCpuBackend();
    } catch (const std::bad_alloc&) {
      result.error = api::Error(api::ErrorCode::kAllocationFailure,
                                "failed to allocate serial CPU backend");
      return result;
    }
    return result;
  }

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  if (selected_kind == Kind::kCuda) {
    result.error = cuda::createCudaBackend(result.backend);
    return result;
  }
#endif

  result.error = api::Error(api::ErrorCode::kInternal,
                            "HiPDLP backend registration has no factory");
  return result;
}

}  // namespace backend
}  // namespace hipdlp