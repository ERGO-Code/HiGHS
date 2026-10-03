#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>

#include <cfenv>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef __GLIBC__
#include <gnu/libc-version.h>
#endif
#if defined(__x86_64__) && defined(__SSE2_MATH__)
#include <xmmintrin.h>
#endif

#include "pdlp/hipdlp/backend/cuda/cuda_backend.hpp"
#ifdef HIGHS_HIPDLP_TESTING
#include "pdlp/hipdlp/backend/cuda/cuda_backend_testing.hpp"
#endif
#include "pdlp/hipdlp/backend/cuda/cuda_kernels.cuh"

namespace hipdlp {
namespace backend {
namespace cuda {

namespace testing {

#ifdef HIGHS_HIPDLP_TESTING
namespace {
thread_local FailurePoint current_failure_point = FailurePoint::kNone;
thread_local ResourceAudit current_resource_audit;
thread_local bool graph_capture_active = false;
}  // namespace

ResourceAudit::ResourceAudit()
    : stream_creations(0),
      stream_destructions(0),
      cublas_creations(0),
      cublas_destructions(0),
      cusparse_creations(0),
      cusparse_destructions(0),
      allocation_calls(0),
      allocation_bytes(0),
      max_allocation_bytes(0),
      free_calls(0),
      free_bytes(0),
      pinned_allocation_calls(0),
      pinned_allocation_bytes(0),
      max_pinned_allocation_bytes(0),
      pinned_free_calls(0),
      pinned_free_bytes(0),
      descriptor_creations(0),
      descriptor_destructions(0),
      host_to_device_transfers(0),
      host_to_device_bytes(0),
      device_to_host_transfers(0),
      device_to_host_bytes(0),
      device_to_device_transfers(0),
      device_to_device_bytes(0),
      graph_capture_calls(0),
      graph_creations(0),
      graph_destructions(0),
      graph_instantiation_calls(0),
      graph_exec_creations(0),
      graph_exec_destructions(0),
      graph_launch_calls(0),
      graph_update_calls(0),
      graph_capture_primitive_enqueue_calls(0),
      direct_primitive_enqueue_calls(0),
      synchronization_calls(0) {}

void setFailurePoint(const FailurePoint failure_point) {
  current_failure_point = failure_point;
}

void resetResourceAudit() { current_resource_audit = ResourceAudit(); }

ResourceAudit resourceAudit() { return current_resource_audit; }

static FailurePoint failurePoint() { return current_failure_point; }

static ResourceAudit& mutableResourceAudit() { return current_resource_audit; }

static void setGraphCaptureActive(const bool active) {
  graph_capture_active = active;
}

static void recordPrimitiveEnqueue() {
  if (graph_capture_active)
    ++current_resource_audit.graph_capture_primitive_enqueue_calls;
  else
    ++current_resource_audit.direct_primitive_enqueue_calls;
}
#endif

}  // namespace testing

#ifdef HIGHS_HIPDLP_TESTING
#define HIPDLP_CUDA_AUDIT(field) (++testing::mutableResourceAudit().field)
#define HIPDLP_CUDA_AUDIT_ADD(field, amount) \
  (testing::mutableResourceAudit().field += (amount))
#define HIPDLP_CUDA_AUDIT_MAX(field, amount)                   \
  do {                                                      \
    if (testing::mutableResourceAudit().field < (amount))     \
      testing::mutableResourceAudit().field = (amount);       \
  } while (false)
#define HIPDLP_CUDA_FAILURE(point) \
  (testing::failurePoint() == testing::FailurePoint::point)
#define HIPDLP_CUDA_CAPTURE_ACTIVE(active) \
  testing::setGraphCaptureActive(active)
#define HIPDLP_CUDA_RECORD_PRIMITIVE() testing::recordPrimitiveEnqueue()
#else
#define HIPDLP_CUDA_AUDIT(field) ((void)0)
#define HIPDLP_CUDA_AUDIT_ADD(field, amount) ((void)0)
#define HIPDLP_CUDA_AUDIT_MAX(field, amount) ((void)0)
#define HIPDLP_CUDA_FAILURE(point) false
#define HIPDLP_CUDA_CAPTURE_ACTIVE(active) ((void)0)
#define HIPDLP_CUDA_RECORD_PRIMITIVE() ((void)0)
#endif

namespace {

const unsigned long long kCudaCapabilityMask =
    static_cast<unsigned long long>(Capability::kSparseMatrixVector) |
    static_cast<unsigned long long>(Capability::kFusedPrimalDualUpdate) |
    static_cast<unsigned long long>(Capability::kFusedAverageUpdate) |
    static_cast<unsigned long long>(Capability::kConvergenceReduction) |
    static_cast<unsigned long long>(Capability::kVectorOperations) |
    static_cast<unsigned long long>(Capability::kReductions) |
    static_cast<unsigned long long>(Capability::kHostDownload) |
    static_cast<unsigned long long>(Capability::kCommandBuffer);

api::Error unavailable(const std::string& message) {
  return api::Error(api::ErrorCode::kBackendUnavailable,
                    "HiPDLP CUDA backend is unavailable: " + message);
}

bool unavailableStatus(const cudaError_t status) {
  return status == cudaErrorNoDevice || status == cudaErrorInsufficientDriver ||
         status == cudaErrorCallRequiresNewerDriver ||
         status == cudaErrorIncompatibleDriverContext ||
         status == cudaErrorDevicesUnavailable ||
         status == cudaErrorNoKernelImageForDevice ||
         status == cudaErrorInvalidDeviceFunction ||
         status == cudaErrorInvalidPtx ||
         status == cudaErrorUnsupportedPtxVersion ||
         status == cudaErrorSystemDriverMismatch;
}

api::Error cudaUnavailable(const char* operation, const cudaError_t status) {
  std::ostringstream message;
  message << operation << " failed with " << cudaGetErrorName(status) << " ("
          << static_cast<int>(status) << "): " << cudaGetErrorString(status);
  return unavailable(message.str());
}

api::Error cudaFailure(const char* operation, const cudaError_t status) {
  std::ostringstream message;
  message << "HiPDLP CUDA " << operation << " failed with "
          << cudaGetErrorName(status) << " (" << static_cast<int>(status)
          << "): " << cudaGetErrorString(status);
  return api::Error(status == cudaErrorMemoryAllocation
                        ? api::ErrorCode::kAllocationFailure
                        : api::ErrorCode::kBackendFailure,
                    message.str());
}

api::Error cublasFailure(const char* operation, const cublasStatus_t status) {
  std::ostringstream message;
  message << "HiPDLP CUDA cuBLAS " << operation << " failed with status "
          << static_cast<int>(status);
  return api::Error(status == CUBLAS_STATUS_ALLOC_FAILED
                        ? api::ErrorCode::kAllocationFailure
                    : status == CUBLAS_STATUS_ARCH_MISMATCH
                        ? api::ErrorCode::kBackendUnavailable
                        : api::ErrorCode::kBackendFailure,
                    message.str());
}

api::Error cusparseFailure(const char* operation,
                           const cusparseStatus_t status) {
  std::ostringstream message;
  message << "HiPDLP CUDA cuSPARSE " << operation << " failed with status "
          << static_cast<int>(status);
  return api::Error(status == CUSPARSE_STATUS_ALLOC_FAILED
                        ? api::ErrorCode::kAllocationFailure
                    : status == CUSPARSE_STATUS_ARCH_MISMATCH
                        ? api::ErrorCode::kBackendUnavailable
                        : api::ErrorCode::kBackendFailure,
                    message.str());
}

api::Error invalidArgument(const std::string& message) {
  return api::Error(api::ErrorCode::kInvalidArgument, message);
}

api::Error backendFailure(const std::string& message) {
  return api::Error(api::ErrorCode::kBackendFailure, message);
}

api::Error allocationFailure(const std::string& operation) {
  return api::Error(api::ErrorCode::kAllocationFailure,
                    "HiPDLP CUDA allocation failed while " + operation);
}

api::Error unsupported(const char* capability) {
  return api::Error(
      api::ErrorCode::kUnsupportedCapability,
      std::string("HiPDLP CUDA backend does not support ") + capability);
}

bool finiteVector(const std::vector<double>& values) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index])) return false;
  }
  return true;
}

class TimingEvents {
 public:
  TimingEvents() : start_(nullptr), stop_(nullptr) {}
  ~TimingEvents() {
    if (stop_ != nullptr) cudaEventDestroy(stop_);
    if (start_ != nullptr) cudaEventDestroy(start_);
  }

  TimingEvents(const TimingEvents&) = delete;
  TimingEvents& operator=(const TimingEvents&) = delete;

  api::Error prepare() {
    if (start_ != nullptr) return api::Error();
    cudaEvent_t start = nullptr;
    cudaError_t status = cudaEventCreate(&start);
    if (status != cudaSuccess) return cudaFailure("timer start creation", status);
    cudaEvent_t stop = nullptr;
    status = cudaEventCreate(&stop);
    if (status != cudaSuccess) {
      cudaEventDestroy(start);
      return cudaFailure("timer stop creation", status);
    }
    start_ = start;
    stop_ = stop;
    return api::Error();
  }

  cudaEvent_t start() const { return start_; }
  cudaEvent_t stop() const { return stop_; }

 private:
  cudaEvent_t start_;
  cudaEvent_t stop_;
};

class Stream {
 public:
  Stream() : value_(nullptr) {}
  ~Stream() { reset(); }

  Stream(Stream&& other) noexcept : value_(other.value_) {
    other.value_ = nullptr;
  }
  Stream& operator=(Stream&& other) noexcept {
    if (this != &other) {
      reset();
      value_ = other.value_;
      other.value_ = nullptr;
    }
    return *this;
  }

  Stream(const Stream&) = delete;
  Stream& operator=(const Stream&) = delete;

  api::Error create() {
    const cudaError_t status =
        cudaStreamCreateWithFlags(&value_, cudaStreamNonBlocking);
    if (status != cudaSuccess) return cudaFailure("stream creation", status);
    HIPDLP_CUDA_AUDIT(stream_creations);
    return api::Error();
  }

  cudaStream_t get() const { return value_; }

 private:
  void reset() {
    if (value_ != nullptr) {
      cudaStreamDestroy(value_);
      HIPDLP_CUDA_AUDIT(stream_destructions);
    }
    value_ = nullptr;
  }

  cudaStream_t value_;
};

class Graph {
 public:
  Graph() : value_(nullptr) {}
  ~Graph() { reset(); }

  Graph(Graph&& other) noexcept : value_(other.value_) {
    other.value_ = nullptr;
  }
  Graph& operator=(Graph&& other) noexcept {
    if (this != &other) {
      reset();
      value_ = other.value_;
      other.value_ = nullptr;
    }
    return *this;
  }

  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  void adopt(const cudaGraph_t graph) {
    reset();
    value_ = graph;
    if (value_ != nullptr) HIPDLP_CUDA_AUDIT(graph_creations);
  }

  cudaGraph_t get() const { return value_; }

 private:
  void reset() {
    if (value_ != nullptr) {
      cudaGraphDestroy(value_);
      HIPDLP_CUDA_AUDIT(graph_destructions);
    }
    value_ = nullptr;
  }

  cudaGraph_t value_;
};

class GraphExecutable {
 public:
  GraphExecutable() : value_(nullptr) {}
  ~GraphExecutable() { reset(); }

  GraphExecutable(GraphExecutable&& other) noexcept : value_(other.value_) {
    other.value_ = nullptr;
  }
  GraphExecutable& operator=(GraphExecutable&& other) noexcept {
    if (this != &other) {
      reset();
      value_ = other.value_;
      other.value_ = nullptr;
    }
    return *this;
  }

  GraphExecutable(const GraphExecutable&) = delete;
  GraphExecutable& operator=(const GraphExecutable&) = delete;

  void adopt(const cudaGraphExec_t executable) {
    reset();
    value_ = executable;
    if (value_ != nullptr) HIPDLP_CUDA_AUDIT(graph_exec_creations);
  }

  cudaGraphExec_t get() const { return value_; }

 private:
  void reset() {
    if (value_ != nullptr) {
      cudaGraphExecDestroy(value_);
      HIPDLP_CUDA_AUDIT(graph_exec_destructions);
    }
    value_ = nullptr;
  }

  cudaGraphExec_t value_;
};

template <typename T>
class DeviceBuffer {
 public:
  DeviceBuffer() : value_(nullptr), size_(0) {}
  ~DeviceBuffer() { reset(); }

  DeviceBuffer(DeviceBuffer&& other) noexcept
      : value_(other.value_), size_(other.size_) {
    other.value_ = nullptr;
    other.size_ = 0;
  }
  DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
    if (this != &other) {
      reset();
      value_ = other.value_;
      size_ = other.size_;
      other.value_ = nullptr;
      other.size_ = 0;
    }
    return *this;
  }

  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;

  api::Error allocate(const std::size_t size, const std::string& operation) {
    reset();
    if (size == 0) return api::Error();
    if (size > std::numeric_limits<std::size_t>::max() / sizeof(T))
      return allocationFailure(operation);
    if (HIPDLP_CUDA_FAILURE(kDeviceAllocation))
      return allocationFailure(operation + " (injected failure)");
    const cudaError_t status =
        cudaMalloc(reinterpret_cast<void**>(&value_), size * sizeof(T));
    if (status != cudaSuccess) return cudaFailure(operation.c_str(), status);
    size_ = size;
    HIPDLP_CUDA_AUDIT(allocation_calls);
    HIPDLP_CUDA_AUDIT_ADD(allocation_bytes, bytes());
    HIPDLP_CUDA_AUDIT_MAX(max_allocation_bytes, bytes());
    return api::Error();
  }

  T* get() { return value_; }
  const T* get() const { return value_; }
  std::size_t size() const { return size_; }
  std::size_t bytes() const { return size_ * sizeof(T); }

 private:
  void reset() {
    if (value_ != nullptr) {
      cudaFree(value_);
      HIPDLP_CUDA_AUDIT(free_calls);
      HIPDLP_CUDA_AUDIT_ADD(free_bytes, bytes());
    }
    value_ = nullptr;
    size_ = 0;
  }

  T* value_;
  std::size_t size_;
};

template <typename T>
class PinnedDownloadBuffer {
 public:
  PinnedDownloadBuffer() : value_(nullptr), size_(0) {}
  ~PinnedDownloadBuffer() {
    if (value_ != nullptr) {
      cudaFreeHost(value_);
      HIPDLP_CUDA_AUDIT(pinned_free_calls);
      HIPDLP_CUDA_AUDIT_ADD(pinned_free_bytes, size_ * sizeof(T));
    }
  }

  PinnedDownloadBuffer(const PinnedDownloadBuffer&) = delete;
  PinnedDownloadBuffer& operator=(const PinnedDownloadBuffer&) = delete;

  api::Error allocate(const std::size_t size) {
    if (size <= size_) return api::Error();
    if (size > std::numeric_limits<std::size_t>::max() / sizeof(T))
      return allocationFailure("allocating pinned vector download storage");
    if (HIPDLP_CUDA_FAILURE(kPinnedAllocation))
      return allocationFailure("allocating pinned download storage (injected failure)");
    T* allocated = nullptr;
    const cudaError_t status =
        cudaMallocHost(reinterpret_cast<void**>(&allocated), size * sizeof(T));
    if (status != cudaSuccess)
      return cudaFailure("pinned download allocation", status);
    HIPDLP_CUDA_AUDIT(pinned_allocation_calls);
    HIPDLP_CUDA_AUDIT_ADD(pinned_allocation_bytes, size * sizeof(T));
    HIPDLP_CUDA_AUDIT_MAX(max_pinned_allocation_bytes, size * sizeof(T));
    if (value_ != nullptr) {
      cudaFreeHost(value_);
      HIPDLP_CUDA_AUDIT(pinned_free_calls);
      HIPDLP_CUDA_AUDIT_ADD(pinned_free_bytes, size_ * sizeof(T));
    }
    value_ = allocated;
    size_ = size;
    return api::Error();
  }

  T* get() const { return value_; }

 private:
  T* value_;
  std::size_t size_;
};

class PinnedCheckPacket {
 public:
  PinnedCheckPacket() : value_(nullptr) {}
  ~PinnedCheckPacket() {
    if (value_ != nullptr) {
      cudaFreeHost(value_);
      HIPDLP_CUDA_AUDIT(pinned_free_calls);
      HIPDLP_CUDA_AUDIT_ADD(pinned_free_bytes, sizeof(CheckPacket));
    }
  }

  PinnedCheckPacket(const PinnedCheckPacket&) = delete;
  PinnedCheckPacket& operator=(const PinnedCheckPacket&) = delete;

  api::Error allocate() {
    if (HIPDLP_CUDA_FAILURE(kPinnedAllocation))
      return allocationFailure("allocating pinned check packet (injected failure)");
    const cudaError_t status =
        cudaMallocHost(reinterpret_cast<void**>(&value_), sizeof(CheckPacket));
    if (status != cudaSuccess)
      return cudaFailure("pinned check packet allocation", status);
    HIPDLP_CUDA_AUDIT(pinned_allocation_calls);
    HIPDLP_CUDA_AUDIT_ADD(pinned_allocation_bytes, sizeof(CheckPacket));
    HIPDLP_CUDA_AUDIT_MAX(max_pinned_allocation_bytes, sizeof(CheckPacket));
    return api::Error();
  }

  CheckPacket* get() const { return value_; }

 private:
  CheckPacket* value_;
};

template <typename T>
api::Error upload(const std::vector<T>& values, DeviceBuffer<T>& buffer,
                  const cudaStream_t stream, const char* operation) {
  if (values.empty()) return api::Error();
  const cudaError_t status =
      cudaMemcpyAsync(buffer.get(), values.data(), values.size() * sizeof(T),
                      cudaMemcpyHostToDevice, stream);
  if (status != cudaSuccess) return cudaFailure(operation, status);
  HIPDLP_CUDA_AUDIT(host_to_device_transfers);
  HIPDLP_CUDA_AUDIT_ADD(host_to_device_bytes, values.size() * sizeof(T));
  return api::Error();
}

class CublasHandle {
 public:
  CublasHandle() : value_(nullptr) {}
  ~CublasHandle() { reset(); }

  CublasHandle(CublasHandle&& other) noexcept : value_(other.value_) {
    other.value_ = nullptr;
  }
  CublasHandle& operator=(CublasHandle&& other) noexcept {
    if (this != &other) {
      reset();
      value_ = other.value_;
      other.value_ = nullptr;
    }
    return *this;
  }

  CublasHandle(const CublasHandle&) = delete;
  CublasHandle& operator=(const CublasHandle&) = delete;

  api::Error create(const cudaStream_t stream) {
    cublasStatus_t status = cublasCreate(&value_);
    if (status != CUBLAS_STATUS_SUCCESS)
      return cublasFailure("handle creation", status);
    HIPDLP_CUDA_AUDIT(cublas_creations);
    status = cublasSetStream(value_, stream);
    if (status != CUBLAS_STATUS_SUCCESS)
      return cublasFailure("stream binding", status);
    status = cublasSetPointerMode(value_, CUBLAS_POINTER_MODE_DEVICE);
    return status == CUBLAS_STATUS_SUCCESS
               ? api::Error()
               : cublasFailure("pointer mode selection", status);
  }

  cublasHandle_t get() const { return value_; }

 private:
  void reset() {
    if (value_ != nullptr) {
      cublasDestroy(value_);
      HIPDLP_CUDA_AUDIT(cublas_destructions);
    }
    value_ = nullptr;
  }

  cublasHandle_t value_;
};

class CusparseHandle {
 public:
  CusparseHandle() : value_(nullptr) {}
  ~CusparseHandle() { reset(); }

  CusparseHandle(CusparseHandle&& other) noexcept : value_(other.value_) {
    other.value_ = nullptr;
  }
  CusparseHandle& operator=(CusparseHandle&& other) noexcept {
    if (this != &other) {
      reset();
      value_ = other.value_;
      other.value_ = nullptr;
    }
    return *this;
  }

  CusparseHandle(const CusparseHandle&) = delete;
  CusparseHandle& operator=(const CusparseHandle&) = delete;

  api::Error create(const cudaStream_t stream) {
    cusparseStatus_t status = cusparseCreate(&value_);
    if (status != CUSPARSE_STATUS_SUCCESS)
      return cusparseFailure("handle creation", status);
    HIPDLP_CUDA_AUDIT(cusparse_creations);
    status = cusparseSetStream(value_, stream);
    if (status != CUSPARSE_STATUS_SUCCESS)
      return cusparseFailure("stream binding", status);
    status = cusparseSetPointerMode(value_, CUSPARSE_POINTER_MODE_DEVICE);
    return status == CUSPARSE_STATUS_SUCCESS
               ? api::Error()
               : cusparseFailure("pointer mode selection", status);
  }

  cusparseHandle_t get() const { return value_; }

 private:
  void reset() {
    if (value_ != nullptr) {
      cusparseDestroy(value_);
      HIPDLP_CUDA_AUDIT(cusparse_destructions);
    }
    value_ = nullptr;
  }

  cusparseHandle_t value_;
};

class DenseVectorDescriptor {
 public:
  DenseVectorDescriptor() : value_(nullptr) {}
  ~DenseVectorDescriptor() { reset(); }

  DenseVectorDescriptor(DenseVectorDescriptor&& other) noexcept
      : value_(other.value_) {
    other.value_ = nullptr;
  }
  DenseVectorDescriptor& operator=(DenseVectorDescriptor&& other) noexcept {
    if (this != &other) {
      reset();
      value_ = other.value_;
      other.value_ = nullptr;
    }
    return *this;
  }

  DenseVectorDescriptor(const DenseVectorDescriptor&) = delete;
  DenseVectorDescriptor& operator=(const DenseVectorDescriptor&) = delete;

  api::Error create(const HighsInt size, double* values) {
    if (size == 0) return api::Error();
    const cusparseStatus_t status = cusparseCreateDnVec(
        &value_, static_cast<int64_t>(size), values, CUDA_R_64F);
    if (status != CUSPARSE_STATUS_SUCCESS)
      return cusparseFailure("dense vector descriptor creation", status);
    HIPDLP_CUDA_AUDIT(descriptor_creations);
    return api::Error();
  }

  cusparseDnVecDescr_t get() const { return value_; }

 private:
  void reset() {
    if (value_ != nullptr) {
      cusparseDestroyDnVec(value_);
      HIPDLP_CUDA_AUDIT(descriptor_destructions);
    }
    value_ = nullptr;
  }

  cusparseDnVecDescr_t value_;
};

class SparseMatrixDescriptor {
 public:
  SparseMatrixDescriptor() : value_(nullptr) {}
  ~SparseMatrixDescriptor() { reset(); }

  SparseMatrixDescriptor(SparseMatrixDescriptor&& other) noexcept
      : value_(other.value_) {
    other.value_ = nullptr;
  }
  SparseMatrixDescriptor& operator=(SparseMatrixDescriptor&& other) noexcept {
    if (this != &other) {
      reset();
      value_ = other.value_;
      other.value_ = nullptr;
    }
    return *this;
  }

  SparseMatrixDescriptor(const SparseMatrixDescriptor&) = delete;
  SparseMatrixDescriptor& operator=(const SparseMatrixDescriptor&) = delete;

  api::Error create(const model::SparseMatrix& matrix, HighsInt* starts,
                    HighsInt* indices, double* values) {
    if (matrix.numNonzeros() == 0) return api::Error();
    const cusparseIndexType_t index_type = sizeof(HighsInt) == sizeof(int32_t)
                                               ? CUSPARSE_INDEX_32I
                                               : CUSPARSE_INDEX_64I;
    const cusparseStatus_t status = cusparseCreateCsr(
        &value_, static_cast<int64_t>(matrix.num_rows),
        static_cast<int64_t>(matrix.num_cols),
        static_cast<int64_t>(matrix.numNonzeros()), starts, indices, values,
        index_type, index_type, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
    if (status != CUSPARSE_STATUS_SUCCESS)
      return cusparseFailure("sparse matrix descriptor creation", status);
    HIPDLP_CUDA_AUDIT(descriptor_creations);
    return api::Error();
  }

  cusparseSpMatDescr_t get() const { return value_; }

 private:
  void reset() {
    if (value_ != nullptr) {
      cusparseDestroySpMat(value_);
      HIPDLP_CUDA_AUDIT(descriptor_destructions);
    }
    value_ = nullptr;
  }

  cusparseSpMatDescr_t value_;
};

// CSR_ALG2 is bitwise deterministic across runs, unlike the default ALG1.
const cusparseSpMVAlg_t kSpmvAlgorithm = CUSPARSE_SPMV_CSR_ALG2;

struct SparseStorage {
  HighsInt num_rows;
  HighsInt num_cols;
  HighsInt num_nonzeros;
  std::size_t spmv_buffer_size;
  DeviceBuffer<HighsInt> starts;
  DeviceBuffer<HighsInt> indices;
  DeviceBuffer<double> values;
  // Holds this matrix's preprocessed SpMV plan, so it is never shared.
  DeviceBuffer<unsigned char> spmv_buffer;
  SparseMatrixDescriptor descriptor;

  SparseStorage()
      : num_rows(0), num_cols(0), num_nonzeros(0), spmv_buffer_size(0) {}

  api::Error initialize(const model::SparseMatrix& matrix,
                        const cudaStream_t stream, const char* name) {
    num_rows = matrix.num_rows;
    num_cols = matrix.num_cols;
    num_nonzeros = matrix.numNonzeros();
    api::Error error =
        starts.allocate(matrix.starts.size(),
                        std::string("allocating ") + name + " row starts");
    if (!error.ok()) return error;
    error =
        indices.allocate(matrix.indices.size(),
                         std::string("allocating ") + name + " column indices");
    if (!error.ok()) return error;
    error = values.allocate(matrix.values.size(),
                            std::string("allocating ") + name + " values");
    if (!error.ok()) return error;
    error = upload(matrix.starts, starts, stream, "matrix row-start upload");
    if (!error.ok()) return error;
    error =
        upload(matrix.indices, indices, stream, "matrix column-index upload");
    if (!error.ok()) return error;
    error = upload(matrix.values, values, stream, "matrix value upload");
    if (!error.ok()) return error;
    return descriptor.create(matrix, starts.get(), indices.get(), values.get());
  }
};

}  // namespace

struct CudaProblemResource : Resource {
  std::shared_ptr<int> lifetime;
  const CudaBackend* owner;
  HighsInt num_rows;
  HighsInt num_cols;
  model::IndexWidth index_width;
  DeviceBuffer<double> costs;
  DeviceBuffer<double> col_lower;
  DeviceBuffer<double> col_upper;
  DeviceBuffer<double> row_lower;
  DeviceBuffer<double> row_upper;
  DeviceBuffer<HighsBool> equality_rows;
  SparseStorage matrix;
  SparseStorage transpose;

  explicit CudaProblemResource(const CudaBackend* backend)
      : lifetime(new int(0)),
        owner(backend),
        num_rows(0),
        num_cols(0),
        index_width(model::IndexWidth::kInt32) {}
};

struct CudaVectorResource : Resource {
  std::shared_ptr<int> lifetime;
  const CudaBackend* owner;
  HighsInt size;
  DeviceBuffer<double> values;
  DenseVectorDescriptor descriptor;
  mutable PinnedDownloadBuffer<double> downloaded;

  explicit CudaVectorResource(const CudaBackend* backend)
      : lifetime(new int(0)), owner(backend), size(0) {}
};

struct CudaCertificateStorage {
  DeviceCertificateProblem data;
  DeviceBuffer<double> costs;
  DeviceBuffer<double> col_lower;
  DeviceBuffer<double> col_upper;
  DeviceBuffer<double> row_lower;
  DeviceBuffer<double> row_upper;
  DeviceBuffer<HighsInt> row_starts;
  DeviceBuffer<HighsInt> column_indices;
  DeviceBuffer<double> matrix_values;
  DeviceBuffer<HighsInt> row_map;
  DeviceBuffer<double> row_sign;
  DeviceBuffer<certificate::Interval> primal;
  DeviceBuffer<certificate::Interval> dual;
  DeviceBuffer<certificate::Interval> column;
  DeviceBuffer<certificate::Interval> activity;
  DeviceBuffer<HighsInt> segment_row;
  DeviceBuffer<HighsInt> segment_begin;
  DeviceBuffer<certificate::Interval> segment_partial;
  DeviceBuffer<HighsInt> long_row;
  DeviceBuffer<HighsInt> long_first;
  DeviceBuffer<HighsInt> long_count;
  DeviceBuffer<CertificateNormPartial> norm_partials;
  DeviceBuffer<CertificateState> state;
  DeviceBuffer<certificate::ScreenPartial> partials;

  CudaCertificateStorage() : data() {}
};

struct CudaCheckStorage {
  const model::Problem* source;
  DeviceCheckProblem data;
  DeviceBuffer<double> costs;
  DeviceBuffer<double> col_lower;
  DeviceBuffer<double> col_upper;
  DeviceBuffer<double> row_lower;
  DeviceBuffer<HighsBool> equality_rows;
  DeviceBuffer<double> row_scale;
  DeviceBuffer<double> col_scale;
  DeviceBuffer<CheckPartial> partials;
  DeviceBuffer<CheckPacket> packet;
  DeviceBuffer<int> optional_error;
  PinnedCheckPacket received;
  CudaCertificateStorage certificate;
  api::Error failure;

  CudaCheckStorage() : source(nullptr), data() {}

  api::Error fail(const api::Error& error) {
    if (failure.ok()) failure = error;
    return failure;
  }
};

struct CudaWorkspaceResource : Resource {
  std::shared_ptr<int> lifetime;
  const CudaBackend* owner;
  HighsInt num_rows;
  HighsInt num_cols;
  DeviceBuffer<double> row_scratch;
  DeviceBuffer<double> col_scratch;
  std::unique_ptr<CudaCheckStorage> check;

  explicit CudaWorkspaceResource(const CudaBackend* backend)
      : lifetime(new int(0)), owner(backend), num_rows(0), num_cols(0) {}
};

struct CudaCommandBufferResource : Resource {
  std::vector<std::weak_ptr<int>> dependencies;
  Graph graph;
  GraphExecutable executable;
};

struct CudaBackend::Impl {
  int device;
  std::string device_name;
  Stream stream;
  TimingEvents timing_events;
  bool timing_active = false;
  CublasHandle cublas;
  CusparseHandle cusparse;
  DeviceBuffer<int> error_flag;
  PinnedDownloadBuffer<int> download_error;
  DeviceBuffer<double> reduction_result;
  DeviceBuffer<FixedPointPartial> fixed_point_partials;
  DeviceBuffer<double> difference_scratch;
  DeviceBuffer<double> spmv_alpha;
  DeviceBuffer<double> spmv_beta;
  bool recording;
  std::vector<std::weak_ptr<int>> recording_dependencies;

  Impl() : device(-1), recording(false) {}

  void recordDependency(const std::shared_ptr<int>& dependency) {
    if (recording) recording_dependencies.push_back(dependency);
  }

  api::Error synchronizeAndCheck(const int* device_error_flag) {
    int host_error = 0;
    cudaError_t status =
        cudaMemcpyAsync(&host_error, device_error_flag, sizeof(host_error),
                        cudaMemcpyDeviceToHost, stream.get());
    if (status != cudaSuccess)
      return cudaFailure("error-state download", status);
    HIPDLP_CUDA_AUDIT(device_to_host_transfers);
    HIPDLP_CUDA_AUDIT_ADD(device_to_host_bytes, sizeof(host_error));
    HIPDLP_CUDA_AUDIT(synchronization_calls);
    status = cudaStreamSynchronize(stream.get());
    if (status != cudaSuccess)
      return cudaFailure("stream synchronization", status);
    if (host_error != 0)
      return backendFailure(
          "HiPDLP CUDA operation produced a non-finite value");
    return api::Error();
  }

  api::Error synchronizeAndCheck() {
    return synchronizeAndCheck(error_flag.get());
  }
};

namespace {

template <typename T>
api::Error allocateAndUpload(const std::vector<T>& values,
                             DeviceBuffer<T>& buffer, const cudaStream_t stream,
                             const std::string& allocation_operation,
                             const char* upload_operation) {
  api::Error error = buffer.allocate(values.size(), allocation_operation);
  if (!error.ok()) return error;
  return upload(values, buffer, stream, upload_operation);
}

bool certificateNormEnvironmentSupported() {
#if defined(__aarch64__) && defined(__GLIBC__) && __GLIBC__ == 2 &&       \
    __GLIBC_MINOR__ == 39 && !defined(__FAST_MATH__) && FLT_RADIX == 2 && \
    DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 && LDBL_MANT_DIG >= 53 &&   \
    LDBL_MAX_EXP >= DBL_MAX_EXP && FLT_EVAL_METHOD == 0
  if (std::strcmp(gnu_get_libc_version(), "2.39") != 0 ||
      std::fegetround() != FE_TONEAREST)
    return false;
  uint64_t control = 0;
  asm volatile("mrs %0, fpcr" : "=r"(control));
  return (control & ((uint64_t(3) << 22) | (uint64_t(1) << 24) | 3)) == 0;
#elif defined(__x86_64__) && defined(__SSE2_MATH__) && defined(__GLIBC__) && \
    __GLIBC__ == 2 && __GLIBC_MINOR__ == 39 && !defined(__FAST_MATH__) &&   \
    FLT_RADIX == 2 && DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 &&          \
    LDBL_MANT_DIG >= 53 && LDBL_MAX_EXP >= DBL_MAX_EXP && FLT_EVAL_METHOD == 0
  if (std::strcmp(gnu_get_libc_version(), "2.39") != 0 ||
      std::fegetround() != FE_TONEAREST)
    return false;
  // SSE rounding control must be nearest (bits 13-14), with neither
  // flush-to-zero (bit 15) nor denormals-are-zero (bit 6).
  return (_mm_getcsr() & ((3u << 13) | (1u << 15) | (1u << 6))) == 0;
#else
  return false;
#endif
}

bool certificateBoundsValid(const std::vector<double>& lower,
                            const std::vector<double>& upper) {
  for (std::size_t index = 0; index < lower.size(); ++index) {
    if (std::isnan(lower[index]) || std::isnan(upper[index]) ||
        lower[index] >= kHighsInf || upper[index] <= -kHighsInf ||
        lower[index] > upper[index])
      return false;
  }
  return true;
}

api::Error originalCertificateMatrix(const model::Problem& problem,
                                     const HighsLp& lp,
                                     model::SparseMatrix& original) {
  const model::TransformationMetadata& transform = problem.transformation();
  if (lp.num_row_ < 0 || lp.num_col_ < 0 ||
      lp.num_row_ != transform.original_num_rows ||
      lp.num_col_ != transform.original_num_cols ||
      lp.num_col_ > problem.numCols() ||
      lp.sense_ != transform.original_sense ||
      (lp.sense_ != ObjSense::kMinimize && lp.sense_ != ObjSense::kMaximize) ||
      lp.col_cost_.size() != static_cast<std::size_t>(lp.num_col_) ||
      lp.col_lower_.size() != static_cast<std::size_t>(lp.num_col_) ||
      lp.col_upper_.size() != static_cast<std::size_t>(lp.num_col_) ||
      lp.row_lower_.size() != static_cast<std::size_t>(lp.num_row_) ||
      lp.row_upper_.size() != static_cast<std::size_t>(lp.num_row_) ||
      transform.original_to_transformed_row.size() !=
          static_cast<std::size_t>(lp.num_row_) ||
      transform.original_row_sign.size() !=
          static_cast<std::size_t>(lp.num_row_))
    return invalidArgument(
        "HiPDLP CUDA original certificate dimensions are invalid");
  if (!certificateBoundsValid(lp.col_lower_, lp.col_upper_) ||
      !certificateBoundsValid(lp.row_lower_, lp.row_upper_))
    return invalidArgument(
        "HiPDLP CUDA original certificate bounds are invalid");
  for (const double cost : lp.col_cost_) {
    if (!std::isfinite(cost))
      return invalidArgument(
          "HiPDLP CUDA original certificate cost is invalid");
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    if (transform.original_to_transformed_row[row] < 0 ||
        transform.original_to_transformed_row[row] >= problem.numRows() ||
        (transform.original_row_sign[row] != 1.0 &&
         transform.original_row_sign[row] != -1.0))
      return invalidArgument(
          "HiPDLP CUDA original certificate row mapping is invalid");
  }
  const HighsSparseMatrix& matrix = lp.a_matrix_;
  const bool colwise = matrix.format_ == MatrixFormat::kColwise;
  const HighsInt vectors = colwise ? lp.num_col_ : lp.num_row_;
  if (matrix.num_col_ != lp.num_col_ || matrix.num_row_ != lp.num_row_ ||
      (!colwise && matrix.format_ != MatrixFormat::kRowwise &&
       matrix.format_ != MatrixFormat::kRowwisePartitioned) ||
      matrix.start_.size() < static_cast<std::size_t>(vectors) + 1 ||
      matrix.start_[0] != 0)
    return invalidArgument(
        "HiPDLP CUDA original certificate matrix is invalid");
  for (HighsInt vector = 0; vector < vectors; ++vector) {
    if (matrix.start_[vector] < 0 ||
        matrix.start_[vector] > matrix.start_[vector + 1])
      return invalidArgument(
          "HiPDLP CUDA original certificate matrix starts are invalid");
  }
  const HighsInt entries = matrix.start_[vectors];
  if (entries < 0 || matrix.index_.size() < static_cast<std::size_t>(entries) ||
      matrix.value_.size() < static_cast<std::size_t>(entries))
    return invalidArgument(
        "HiPDLP CUDA original certificate matrix storage is invalid");
  original.num_rows = lp.num_row_;
  original.num_cols = lp.num_col_;
  original.starts.assign(static_cast<std::size_t>(lp.num_row_) + 1, 0);
  original.indices.resize(entries);
  original.values.resize(entries);
  for (HighsInt vector = 0; vector < vectors; ++vector) {
    for (HighsInt entry = matrix.start_[vector];
         entry < matrix.start_[vector + 1]; ++entry) {
      if (matrix.index_[entry] < 0 ||
          matrix.index_[entry] >= (colwise ? lp.num_row_ : lp.num_col_) ||
          !std::isfinite(matrix.value_[entry]))
        return invalidArgument(
            "HiPDLP CUDA original certificate matrix entry is invalid");
      const HighsInt row = colwise ? matrix.index_[entry] : vector;
      ++original.starts[row + 1];
    }
  }
  for (HighsInt row = 0; row < lp.num_row_; ++row)
    original.starts[row + 1] += original.starts[row];
  std::vector<HighsInt> next = original.starts;
  for (HighsInt vector = 0; vector < vectors; ++vector) {
    for (HighsInt entry = matrix.start_[vector];
         entry < matrix.start_[vector + 1]; ++entry) {
      const HighsInt row = colwise ? matrix.index_[entry] : vector;
      const HighsInt target = next[row]++;
      original.indices[target] = colwise ? vector : matrix.index_[entry];
      original.values[target] = matrix.value_[entry];
    }
  }
  return api::Error();
}

api::Error prepareCertificateStorage(const model::Problem& problem,
                                     const HighsLp& lp,
                                     const model::SparseMatrix& original,
                                     const cudaStream_t stream,
                                     CudaCheckStorage& check) {
  CudaCertificateStorage& storage = check.certificate;
  const std::vector<double>* values[] = {&lp.col_cost_, &lp.col_lower_,
                                         &lp.col_upper_, &lp.row_lower_,
                                         &lp.row_upper_};
  DeviceBuffer<double>* buffers[] = {&storage.costs, &storage.col_lower,
                                     &storage.col_upper, &storage.row_lower,
                                     &storage.row_upper};
  for (std::size_t index = 0; index < 5; ++index) {
    api::Error error =
        allocateAndUpload(*values[index], *buffers[index], stream,
                          "allocating original certificate data",
                          "original certificate data upload");
    if (!error.ok()) return error;
  }
  api::Error error =
      allocateAndUpload(original.starts, storage.row_starts, stream,
                        "allocating original certificate row starts",
                        "original certificate row-start upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(original.indices, storage.column_indices, stream,
                            "allocating original certificate indices",
                            "original certificate index upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(original.values, storage.matrix_values, stream,
                            "allocating original certificate matrix values",
                            "original certificate matrix upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(
      problem.transformation().original_to_transformed_row, storage.row_map,
      stream, "allocating certificate row map", "certificate row-map upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(
      problem.transformation().original_row_sign, storage.row_sign, stream,
      "allocating certificate row signs", "certificate row-sign upload");
  if (!error.ok()) return error;
  error = storage.primal.allocate(lp.num_col_,
                                  "allocating certificate primal intervals");
  if (!error.ok()) return error;
  error = storage.column.allocate(lp.num_col_,
                                  "allocating certificate column intervals");
  if (!error.ok()) return error;
  error = storage.dual.allocate(lp.num_row_,
                                "allocating certificate dual intervals");
  if (!error.ok()) return error;
  error = storage.activity.allocate(
      lp.num_row_, "allocating certificate row activity intervals");
  if (!error.ok()) return error;
  // Row segments: every row gets at least one, and a row longer than
  // kCertificateSegmentLength is split so several warps share it.
  std::vector<HighsInt> segment_row;
  std::vector<HighsInt> segment_begin;
  std::vector<HighsInt> long_row;
  std::vector<HighsInt> long_first;
  std::vector<HighsInt> long_count;
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    const HighsInt start = original.starts[row];
    const HighsInt length = original.starts[row + 1] - start;
    const HighsInt parts =
        length <= kCertificateSegmentLength
            ? 1 : 1 + (length - 1) / kCertificateSegmentLength;
    if (parts > 1) {
      long_row.push_back(row);
      long_first.push_back(static_cast<HighsInt>(segment_row.size()));
      long_count.push_back(parts);
    }
    for (HighsInt part = 0; part < parts; ++part) {
      segment_row.push_back(row);
      segment_begin.push_back(start + part * kCertificateSegmentLength);
    }
  }
  error = allocateAndUpload(segment_row, storage.segment_row, stream,
                            "allocating certificate segment rows",
                            "certificate segment-row upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(segment_begin, storage.segment_begin, stream,
                            "allocating certificate segment starts",
                            "certificate segment-start upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(long_row, storage.long_row, stream,
                            "allocating certificate long rows",
                            "certificate long-row upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(long_first, storage.long_first, stream,
                            "allocating certificate long-row starts",
                            "certificate long-row-start upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(long_count, storage.long_count, stream,
                            "allocating certificate long-row counts",
                            "certificate long-row-count upload");
  if (!error.ok()) return error;
  error = storage.segment_partial.allocate(
      segment_row.size(), "allocating certificate segment partial sums");
  if (!error.ok()) return error;
  const std::size_t count = checkPartialCount(lp.num_row_, lp.num_col_);
  error = storage.norm_partials.allocate(
      count, "allocating certificate norm partials");
  if (!error.ok()) return error;
  error = storage.partials.allocate(
      count, "allocating certificate condition partials");
  if (!error.ok()) return error;
  error = storage.state.allocate(1, "allocating certificate norm bounds");
  if (!error.ok()) return error;
  storage.data.num_rows = lp.num_row_;
  storage.data.num_cols = lp.num_col_;
  storage.data.costs = storage.costs.get();
  storage.data.col_lower = storage.col_lower.get();
  storage.data.col_upper = storage.col_upper.get();
  storage.data.row_lower = storage.row_lower.get();
  storage.data.row_upper = storage.row_upper.get();
  storage.data.row_starts = storage.row_starts.get();
  storage.data.column_indices = storage.column_indices.get();
  storage.data.matrix_values = storage.matrix_values.get();
  storage.data.row_map = storage.row_map.get();
  storage.data.row_sign = storage.row_sign.get();
  storage.data.row_scale = check.row_scale.get();
  storage.data.col_scale = check.col_scale.get();
  storage.data.num_segments = static_cast<HighsInt>(segment_row.size());
  storage.data.segment_row = storage.segment_row.get();
  storage.data.segment_begin = storage.segment_begin.get();
  storage.data.segment_partial = storage.segment_partial.get();
  storage.data.num_long_rows = static_cast<HighsInt>(long_row.size());
  storage.data.long_row = storage.long_row.get();
  storage.data.long_first = storage.long_first.get();
  storage.data.long_count = storage.long_count.get();
  storage.data.sense_sign = lp.sense_ == ObjSense::kMinimize ? 1.0 : -1.0;
  storage.data.column_bounds_straddle_zero = true;
  for (HighsInt col = 0; col < lp.num_col_; ++col)
    storage.data.column_bounds_straddle_zero &=
        lp.col_lower_[col] <= 0.0 && lp.col_upper_[col] >= 0.0;
  return api::Error();
}

template <typename ImplType>
api::Error initializeProblemResource(const model::Problem& problem,
                                     ImplType& impl,
                                     CudaProblemResource& resource) {
  resource.num_rows = problem.numRows();
  resource.num_cols = problem.numCols();
  resource.index_width = problem.indexWidth();
  api::Error error =
      allocateAndUpload(problem.costs(), resource.costs, impl.stream.get(),
                        "allocating problem costs", "problem cost upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(problem.colLower(), resource.col_lower,
                            impl.stream.get(), "allocating column lower bounds",
                            "column lower-bound upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(problem.colUpper(), resource.col_upper,
                            impl.stream.get(), "allocating column upper bounds",
                            "column upper-bound upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(problem.rowLower(), resource.row_lower,
                            impl.stream.get(), "allocating row lower bounds",
                            "row lower-bound upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(problem.rowUpper(), resource.row_upper,
                            impl.stream.get(), "allocating row upper bounds",
                            "row upper-bound upload");
  if (!error.ok()) return error;
  error = allocateAndUpload(
      problem.equalityRows(), resource.equality_rows, impl.stream.get(),
      "allocating row classification metadata", "row classification upload");
  if (!error.ok()) return error;
  error = resource.matrix.initialize(problem.matrix(), impl.stream.get(),
                                     "CSR matrix");
  if (!error.ok()) return error;
  return resource.transpose.initialize(problem.transpose(), impl.stream.get(),
                                       "transpose matrix");
}

template <typename ImplType>
api::Error prepareSpmvPlan(ImplType& impl, SparseStorage& matrix) {
  if (matrix.num_nonzeros == 0) {
    matrix.spmv_buffer_size = 0;
    return api::Error();
  }
  DeviceBuffer<double> input;
  DeviceBuffer<double> output;
  api::Error error = input.allocate(static_cast<std::size_t>(matrix.num_cols),
                                    "sizing a cuSPARSE input vector");
  if (!error.ok()) return error;
  error = output.allocate(static_cast<std::size_t>(matrix.num_rows),
                          "sizing a cuSPARSE output vector");
  if (!error.ok()) return error;
  DenseVectorDescriptor input_descriptor;
  DenseVectorDescriptor output_descriptor;
  error = input_descriptor.create(matrix.num_cols, input.get());
  if (!error.ok()) return error;
  error = output_descriptor.create(matrix.num_rows, output.get());
  if (!error.ok()) return error;
  cusparseStatus_t status = cusparseSpMV_bufferSize(
      impl.cusparse.get(), CUSPARSE_OPERATION_NON_TRANSPOSE,
      impl.spmv_alpha.get(), matrix.descriptor.get(), input_descriptor.get(),
      impl.spmv_beta.get(), output_descriptor.get(), CUDA_R_64F,
      kSpmvAlgorithm, &matrix.spmv_buffer_size);
  if (status != CUSPARSE_STATUS_SUCCESS)
    return cusparseFailure("SpMV buffer sizing", status);
  error = matrix.spmv_buffer.allocate(matrix.spmv_buffer_size,
                                      "allocating cuSPARSE SpMV workspace");
  if (!error.ok()) return error;
  // Partition the rows once instead of on every product.
  status = cusparseSpMV_preprocess(
      impl.cusparse.get(), CUSPARSE_OPERATION_NON_TRANSPOSE,
      impl.spmv_alpha.get(), matrix.descriptor.get(), input_descriptor.get(),
      impl.spmv_beta.get(), output_descriptor.get(), CUDA_R_64F,
      kSpmvAlgorithm, matrix.spmv_buffer.get());
  return status == CUSPARSE_STATUS_SUCCESS
             ? api::Error()
             : cusparseFailure("SpMV preprocessing", status);
}

template <typename ImplType>
api::Error prepareSpmvStorage(ImplType& impl, CudaProblemResource& resource) {
  const api::Error error = prepareSpmvPlan(impl, resource.matrix);
  if (!error.ok()) return error;
  return prepareSpmvPlan(impl, resource.transpose);
}

api::Error launchResult(const char* operation, const cudaError_t status) {
  HIPDLP_CUDA_RECORD_PRIMITIVE();
  return status == cudaSuccess ? api::Error() : cudaFailure(operation, status);
}

template <typename ImplType>
api::Error executeSpmv(ImplType& impl, const SparseStorage& matrix,
                       const CudaVectorResource& input,
                       CudaVectorResource& output, int* error_flag,
                       const bool check_output = true) {
  if (matrix.num_nonzeros == 0)
    return launchResult(
        "empty SpMV fill kernel launch",
        launchFill(impl.stream.get(), output.size, 0.0, output.values.get()));
  HIPDLP_CUDA_RECORD_PRIMITIVE();
  const cusparseStatus_t status = cusparseSpMV(
      impl.cusparse.get(), CUSPARSE_OPERATION_NON_TRANSPOSE,
      impl.spmv_alpha.get(), matrix.descriptor.get(), input.descriptor.get(),
      impl.spmv_beta.get(), output.descriptor.get(), CUDA_R_64F,
      kSpmvAlgorithm, const_cast<unsigned char*>(matrix.spmv_buffer.get()));
  if (status != CUSPARSE_STATUS_SUCCESS)
    return cusparseFailure("SpMV execution", status);
  if (!check_output) return api::Error();
  return launchResult("SpMV finite-check kernel launch",
                      launchFiniteCheck(impl.stream.get(), output.size,
                                        output.values.get(), error_flag));
}

template <typename ImplType>
api::Error enqueueNorm(ImplType& impl, const HighsInt size, const double* input,
                        double* output) {
  if (size == 0)
    return launchResult("empty norm reduction",
        launchFill(impl.stream.get(), 1, 0.0, output));
#if CUBLAS_VER_MAJOR >= 12
  HIPDLP_CUDA_RECORD_PRIMITIVE();
  const cublasStatus_t status = cublasDnrm2_64(
      impl.cublas.get(), static_cast<int64_t>(size), input, 1, output);
#else
  if (size > std::numeric_limits<int>::max())
    return launchResult("wide norm reduction",
        launchNorm2(impl.stream.get(), size, input, output, impl.error_flag.get()));
  HIPDLP_CUDA_RECORD_PRIMITIVE();
  const cublasStatus_t status = cublasDnrm2(
      impl.cublas.get(), static_cast<int>(size), input, 1, output);
#endif
  return status == CUBLAS_STATUS_SUCCESS ? api::Error()
                                         : cublasFailure("norm reduction", status);
}

template <typename ImplType>
api::Error enqueueDot(ImplType& impl, const HighsInt size, const double* left,
                       const double* right, double* output) {
  if (size == 0)
    return launchResult("empty dot reduction",
        launchFill(impl.stream.get(), 1, 0.0, output));
#if CUBLAS_VER_MAJOR >= 12
  HIPDLP_CUDA_RECORD_PRIMITIVE();
  const cublasStatus_t status = cublasDdot_64(
      impl.cublas.get(), static_cast<int64_t>(size), left, 1, right, 1, output);
#else
  if (size > std::numeric_limits<int>::max())
    return launchResult("wide dot reduction",
        launchDot(impl.stream.get(), size, left, right, output, impl.error_flag.get()));
  HIPDLP_CUDA_RECORD_PRIMITIVE();
  const cublasStatus_t status = cublasDdot(
      impl.cublas.get(), static_cast<int>(size), left, 1, right, 1, output);
#endif
  return status == CUBLAS_STATUS_SUCCESS ? api::Error()
                                         : cublasFailure("dot reduction", status);
}

bool dependenciesAlive(const CudaCommandBufferResource& command_buffer) {
  for (std::size_t index = 0; index < command_buffer.dependencies.size();
       ++index) {
    if (command_buffer.dependencies[index].expired()) return false;
  }
  return true;
}

}  // namespace

CudaBackend::CudaBackend(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}

CudaBackend::~CudaBackend() {}

Descriptor cudaDescriptor() {
  return {Kind::kCuda,
          "cuda",
          true,
          true,
          PrecisionDescriptor(),
          Capabilities(kCudaCapabilityMask)};
}

Descriptor CudaBackend::descriptor() const { return cudaDescriptor(); }

api::Error CudaBackend::queryDevice(DeviceInfo& device) const {
  device.available = impl_.get() != nullptr && impl_->device >= 0;
  device.name = device.available ? impl_->device_name : std::string();
  return device.available ? api::Error()
                          : unavailable("no compatible NVIDIA device");
}

api::Error CudaBackend::requireCapability(const Capability capability) const {
  if (!descriptor().capabilities.supports(capability))
    return unsupported(capability == Capability::kBatchedVectors
                           ? "batched vectors"
                       : capability == Capability::kSparseMatrixMatrix
                           ? "SpMM"
                           : "the requested capability");
  return api::Error();
}

CudaProblemResource* CudaBackend::problemResource(ProblemHandle& handle) const {
  CudaProblemResource* resource =
      dynamic_cast<CudaProblemResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

const CudaProblemResource* CudaBackend::problemResource(
    const ProblemHandle& handle) const {
  const CudaProblemResource* resource =
      dynamic_cast<const CudaProblemResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

CudaVectorResource* CudaBackend::vectorResource(VectorHandle& handle) const {
  CudaVectorResource* resource =
      dynamic_cast<CudaVectorResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

const CudaVectorResource* CudaBackend::vectorResource(
    const VectorHandle& handle) const {
  const CudaVectorResource* resource =
      dynamic_cast<const CudaVectorResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

CudaWorkspaceResource* CudaBackend::workspaceResource(
    WorkspaceHandle& handle) const {
  CudaWorkspaceResource* resource =
      dynamic_cast<CudaWorkspaceResource*>(handle.resource_.get());
  return resource != nullptr && resource->owner == this ? resource : nullptr;
}

const CudaCommandBufferResource* CudaBackend::commandBufferResource(
    const CommandBufferHandle& handle) const {
  return dynamic_cast<const CudaCommandBufferResource*>(handle.resource_.get());
}

api::Error CudaBackend::bindProblem(const model::Problem& problem,
                                    ProblemHandle& handle) {
  if (!problem.initialized())
    return invalidArgument("cannot bind an uninitialized HiPDLP problem");
  try {
    std::unique_ptr<CudaProblemResource> resource(new CudaProblemResource(this));
    api::Error error = initializeProblemResource(problem, *impl_, *resource);
    if (!error.ok()) return error;
    error = prepareSpmvStorage(*impl_, *resource);
    if (!error.ok()) return error;
    error = impl_->synchronizeAndCheck();
    if (!error.ok()) return error;
    handle = ProblemHandle(std::unique_ptr<Resource>(resource.release()));
  } catch (const std::bad_alloc&) {
    return allocationFailure("binding a problem");
  }
  return api::Error();
}

api::Error CudaBackend::createWorkspace(const ProblemHandle& problem,
                                        WorkspaceHandle& workspace) {
  const CudaProblemResource* cuda_problem = problemResource(problem);
  if (cuda_problem == nullptr)
    return invalidArgument("HiPDLP CUDA problem handle is invalid");
  try {
    std::unique_ptr<CudaWorkspaceResource> resource(
      new CudaWorkspaceResource(this));
    resource->num_rows = cuda_problem->num_rows;
    resource->num_cols = cuda_problem->num_cols;
    api::Error error = resource->row_scratch.allocate(
        static_cast<std::size_t>(resource->num_rows),
        "creating row workspace storage");
    if (!error.ok()) return error;
    error = resource->col_scratch.allocate(
        static_cast<std::size_t>(resource->num_cols),
        "creating column workspace storage");
    if (!error.ok()) return error;
    if (resource->row_scratch.bytes() > 0) {
      const cudaError_t status =
          cudaMemsetAsync(resource->row_scratch.get(), 0,
                          resource->row_scratch.bytes(), impl_->stream.get());
      if (status != cudaSuccess)
        return cudaFailure("row workspace initialization", status);
    }
    if (resource->col_scratch.bytes() > 0) {
      const cudaError_t status =
          cudaMemsetAsync(resource->col_scratch.get(), 0,
                          resource->col_scratch.bytes(), impl_->stream.get());
      if (status != cudaSuccess)
        return cudaFailure("column workspace initialization", status);
    }
    workspace = WorkspaceHandle(std::unique_ptr<Resource>(resource.release()));
  } catch (const std::bad_alloc&) {
    return allocationFailure("creating a workspace");
  }
  return api::Error();
}

api::Error CudaBackend::prepareCheckWorkspace(const model::Problem& problem,
                                              const HighsLp& original_lp,
                                              WorkspaceHandle& workspace) {
  if (impl_->recording) return unsupported("resident checks during recording");
  CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  if (cuda_workspace == nullptr)
    return invalidArgument("HiPDLP CUDA check workspace is invalid");
  cuda_workspace->check.reset();
  if (!problem.initialized() || cuda_workspace->num_rows != problem.numRows() ||
      cuda_workspace->num_cols != problem.numCols())
    return invalidArgument("HiPDLP CUDA check problem dimensions are invalid");
  if (!std::isfinite(problem.objectiveOffset()) ||
      !std::isfinite(problem.unscaledNorms().row_bounds_l2) ||
      !std::isfinite(problem.unscaledNorms().objective_l2))
    return invalidArgument("HiPDLP CUDA check model constants are non-finite");
  const model::ScalingMetadata& scaling = problem.scaling();
  if (scaling.row_scale.size() != static_cast<std::size_t>(problem.numRows()) ||
      scaling.col_scale.size() != static_cast<std::size_t>(problem.numCols()))
    return invalidArgument("HiPDLP CUDA check scaling dimensions are invalid");
  for (const double scale : scaling.row_scale) {
    if (!std::isfinite(scale) || scale <= 0.0)
      return invalidArgument("HiPDLP CUDA check row scale is invalid");
  }
  for (const double scale : scaling.col_scale) {
    if (!std::isfinite(scale) || scale <= 0.0)
      return invalidArgument("HiPDLP CUDA check column scale is invalid");
  }
  try {
    model::SparseMatrix original_matrix;
    api::Error error =
        originalCertificateMatrix(problem, original_lp, original_matrix);
    if (!error.ok()) return error;
    std::unique_ptr<CudaCheckStorage> check(new CudaCheckStorage());
    error = check->received.allocate();
    if (!error.ok()) return error;
    error = check->packet.allocate(1, "allocating check result packet");
    if (!error.ok()) return error;
    error = check->optional_error.allocate(1, "allocating check optional error");
    if (!error.ok()) return error;
    const cudaError_t status = cudaMemsetAsync(
        check->optional_error.get(), 0, check->optional_error.bytes(),
        impl_->stream.get());
    if (status != cudaSuccess)
      return cudaFailure("check optional error initialization", status);
    error = check->partials.allocate(
        checkPartialCount(problem.numRows(), problem.numCols()),
        "allocating check reduction workspace");
    if (!error.ok()) return error;
    error =
        allocateAndUpload(problem.costs(), check->costs, impl_->stream.get(),
                          "allocating check costs", "check cost upload");
    if (!error.ok()) return error;
    error = allocateAndUpload(
        problem.colLower(), check->col_lower, impl_->stream.get(),
        "allocating check lower bounds", "check lower-bound upload");
    if (!error.ok()) return error;
    error = allocateAndUpload(
        problem.colUpper(), check->col_upper, impl_->stream.get(),
        "allocating check upper bounds", "check upper-bound upload");
    if (!error.ok()) return error;
    error = allocateAndUpload(
        problem.rowLower(), check->row_lower, impl_->stream.get(),
        "allocating check row bounds", "check row-bound upload");
    if (!error.ok()) return error;
    error = allocateAndUpload(problem.equalityRows(), check->equality_rows,
                              impl_->stream.get(), "allocating check row types",
                              "check row-type upload");
    if (!error.ok()) return error;
    error = allocateAndUpload(
        scaling.row_scale, check->row_scale, impl_->stream.get(),
        "allocating check row scales", "check row-scale upload");
    if (!error.ok()) return error;
    error = allocateAndUpload(
        scaling.col_scale, check->col_scale, impl_->stream.get(),
        "allocating check column scales", "check column-scale upload");
    if (!error.ok()) return error;
    error = prepareCertificateStorage(problem, original_lp, original_matrix,
                                      impl_->stream.get(), *check);
    if (!error.ok()) return error;
    error = impl_->synchronizeAndCheck();
    if (!error.ok()) return error;
    check->source = &problem;
    check->data.num_rows = problem.numRows();
    check->data.num_cols = problem.numCols();
    check->data.costs = check->costs.get();
    check->data.col_lower = check->col_lower.get();
    check->data.col_upper = check->col_upper.get();
    check->data.row_lower = check->row_lower.get();
    check->data.equality_rows = check->equality_rows.get();
    check->data.row_scale = check->row_scale.get();
    check->data.col_scale = check->col_scale.get();
    check->data.objective_offset = problem.objectiveOffset();
    check->data.unscaled_rhs_norm = problem.unscaledNorms().row_bounds_l2;
    check->data.unscaled_cost_norm = problem.unscaledNorms().objective_l2;
    cuda_workspace->check = std::move(check);
  } catch (const std::bad_alloc&) {
    return allocationFailure("preparing resident check storage");
  }
  return api::Error();
}

api::Error CudaBackend::prepareCheckFixedPoint(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& primal, const VectorHandle& dual,
    const VectorHandle& reflected_primal, const VectorHandle& reflected_dual,
    VectorHandle& delta_primal, VectorHandle& delta_dual,
    VectorHandle& aty_delta) {
  if (impl_->recording) return unsupported("resident checks during recording");
  const CudaProblemResource* cuda_problem = problemResource(problem);
  CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  if (cuda_problem == nullptr || cuda_workspace == nullptr ||
      !cuda_workspace->check ||
      cuda_workspace->num_rows != cuda_problem->num_rows ||
      cuda_workspace->num_cols != cuda_problem->num_cols)
    return invalidArgument("HiPDLP CUDA check preparation workspace is invalid");
  CudaCheckStorage& check = *cuda_workspace->check;
  if (!check.failure.ok()) return check.failure;
  const VectorHandle* handles[] = {
      &primal, &dual, &reflected_primal, &reflected_dual,
      &delta_primal, &delta_dual, &aty_delta};
  const HighsInt sizes[] = {
      cuda_problem->num_cols, cuda_problem->num_rows, cuda_problem->num_cols,
      cuda_problem->num_rows, cuda_problem->num_cols, cuda_problem->num_rows,
      cuda_problem->num_cols};
  const CudaVectorResource* vectors[7] = {};
  for (std::size_t index = 0; index < 7; ++index) {
    vectors[index] = vectorResource(*handles[index]);
    if (vectors[index] == nullptr || vectors[index]->size != sizes[index])
      return invalidArgument(
          "HiPDLP CUDA check preparation vector dimensions or owner are invalid");
    for (std::size_t previous = 0; previous < index; ++previous) {
      if (vectors[index] == vectors[previous])
        return invalidArgument(
            "HiPDLP CUDA check preparation vectors must be distinct");
    }
  }
  const cudaError_t status = cudaMemsetAsync(
      check.optional_error.get(), 0, check.optional_error.bytes(),
      impl_->stream.get());
  if (status != cudaSuccess)
    return check.fail(cudaFailure("check optional error reset", status));
    api::Error error = launchResult(
      "check primal delta kernel launch",
      launchDifference(impl_->stream.get(), cuda_problem->num_cols,
               vectors[0]->values.get(), vectors[2]->values.get(),
               vectorResource(delta_primal)->values.get(), check.optional_error.get()));
  if (!error.ok()) return check.fail(error);
  error = launchResult(
      "check dual delta kernel launch",
      launchDifference(impl_->stream.get(), cuda_problem->num_rows,
               vectors[1]->values.get(), vectors[3]->values.get(),
               vectorResource(delta_dual)->values.get(), check.optional_error.get()));
  if (!error.ok()) return check.fail(error);
  error = executeSpmv(*impl_, cuda_problem->transpose,
                      *vectors[5], *vectorResource(aty_delta),
                      check.optional_error.get());
  return error.ok() ? error : check.fail(error);
}

api::Error CudaBackend::evaluateCheck(const model::Problem& problem,
                                      WorkspaceHandle& workspace,
                                      const CheckInput& input,
                                      CheckResult& result) {
  result = CheckResult();
  if (!input.validTolerances())
    return invalidArgument("HiPDLP check tolerances must be finite and nonnegative");
  if (impl_->recording) return unsupported("resident checks during recording");
  CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  if (cuda_workspace == nullptr || !cuda_workspace->check ||
      cuda_workspace->check->source != &problem || !problem.initialized() ||
      cuda_workspace->num_rows != problem.numRows() ||
      cuda_workspace->num_cols != problem.numCols())
    return invalidArgument("HiPDLP CUDA check workspace is not prepared");
  CudaCheckStorage& check = *cuda_workspace->check;
  if (!check.failure.ok()) return check.failure;
  const VectorHandle* handles[] = {
      input.primal,    input.dual,          input.ax,
      input.aty,       input.delta_primal,  input.delta_dual,
      input.aty_delta, input.primal_anchor, input.dual_anchor};
  const HighsInt sizes[] = {
      problem.numCols(), problem.numRows(), problem.numRows(),
      problem.numCols(), problem.numCols(), problem.numRows(),
      problem.numCols(), problem.numCols(), problem.numRows()};
  const CudaVectorResource* vectors[9] = {};
  for (std::size_t index = 0; index < 9; ++index) {
    const bool required =
        index < 4 || (index < 7 ? input.fixed_point : input.anchor_distances);
    if (!required) continue;
    vectors[index] =
        handles[index] == nullptr ? nullptr : vectorResource(*handles[index]);
    if (vectors[index] == nullptr || vectors[index]->size != sizes[index])
      return invalidArgument(
          "HiPDLP CUDA check vector dimensions or owner are invalid");
  }
  DeviceCheckInput device_input = {};
  device_input.primal = vectors[0]->values.get();
  device_input.dual = vectors[1]->values.get();
  device_input.ax = vectors[2]->values.get();
  device_input.aty = vectors[3]->values.get();
  device_input.optional_error = check.optional_error.get();
  device_input.fixed_point = input.fixed_point;
  device_input.anchor_distances = input.anchor_distances;
  device_input.certificate_tolerance = input.certificate_tolerance;
  device_input.primal_tolerance = input.primal_tolerance;
  device_input.dual_tolerance = input.dual_tolerance;
  device_input.gap_tolerance = input.gap_tolerance;
  if (input.fixed_point) {
    device_input.delta_primal = vectors[4]->values.get();
    device_input.delta_dual = vectors[5]->values.get();
    device_input.aty_delta = vectors[6]->values.get();
  }
  if (input.anchor_distances) {
    device_input.primal_anchor = vectors[7]->values.get();
    device_input.dual_anchor = vectors[8]->values.get();
  }
  api::Error error =
      launchResult("check partial reduction kernel launch",
                   launchCheckPartials(impl_->stream.get(), check.data,
                                       device_input, check.partials.get()));
  if (!error.ok()) return check.fail(error);
  error = launchResult(
      "check packet kernel launch",
      launchCheckPacket(impl_->stream.get(), check.data, device_input,
                        check.partials.get(), check.packet.get(),
                        impl_->error_flag.get()));
  if (!error.ok()) return check.fail(error);

  CudaCertificateStorage& certificate = check.certificate;
  const bool oracle_rounding_supported = certificateNormEnvironmentSupported();
  const bool screen_supported =
      input.certificates && oracle_rounding_supported &&
      static_cast<double>(certificate.data.num_rows) + certificate.data.num_cols <
          4503599627370496.0;
  if (screen_supported) {
    error = launchResult(
        "certificate screen kernel launch",
        launchCertificateScreen(
            impl_->stream.get(), certificate.data, device_input,
            certificate.primal.get(), certificate.dual.get(),
            certificate.column.get(), certificate.activity.get(),
            certificate.norm_partials.get(),
            certificate.state.get(), certificate.partials.get(),
            check.packet.get()));
    if (!error.ok()) return check.fail(error);
  }

  CheckResult computed;
  cudaError_t status = cudaMemcpyAsync(
      check.received.get(), check.packet.get(), check.packet.bytes(),
      cudaMemcpyDeviceToHost, impl_->stream.get());
  if (status != cudaSuccess)
    return check.fail(cudaFailure("check packet download", status));
  ++computed.device_to_host_transfers;
  computed.device_to_host_bytes += check.packet.bytes();
  HIPDLP_CUDA_AUDIT(device_to_host_transfers);
  HIPDLP_CUDA_AUDIT_ADD(device_to_host_bytes, check.packet.bytes());
  ++computed.synchronizations;
  HIPDLP_CUDA_AUDIT(synchronization_calls);
  status = cudaStreamSynchronize(impl_->stream.get());
  if (status != cudaSuccess)
    return check.fail(cudaFailure("check packet synchronization", status));
  const CheckPacket& packet = *check.received.get();
  if (packet.error_flag != 0)
    return check.fail(
        backendFailure("HiPDLP CUDA operation produced a non-finite value"));
  const bool stopping_enabled = input.primal_tolerance > 0.0 ||
                                input.dual_tolerance > 0.0 ||
                                input.gap_tolerance > 0.0;
  computed.requires_host_oracle = (packet.fallback_flags & kCheckRawFallback) != 0 ||
                                 (stopping_enabled && !oracle_rounding_supported);
  computed.fixed_point_prepared = input.fixed_point &&
      (packet.fallback_flags & kCheckUnpreparedProducts) == 0;
  computed.fixed_point_valid = computed.fixed_point_prepared &&
      (packet.fallback_flags & kCheckFixedPointFallback) == 0;
  computed.anchor_distances_valid = input.anchor_distances &&
      (packet.fallback_flags & kCheckAnchorFallback) == 0;
  computed.certificate_screen_valid = screen_supported;
  if (!input.certificates) {
    computed.primal_ray_possible = false;
    computed.farkas_possible = false;
  } else if (screen_supported) {
    computed.primal_ray_possible = packet.primal_ray_possible != 0;
    computed.farkas_possible = packet.farkas_possible != 0;
  }
  if (!computed.requires_host_oracle) {
    computed.raw.primal_residual = packet.primal_residual;
    computed.raw.dual_residual = packet.dual_residual;
    computed.raw.objective_gap = packet.objective_gap;
    computed.raw.primal_objective = packet.primal_objective;
    computed.raw.dual_objective = packet.dual_objective;
    computed.raw.unscaled_rhs_norm = packet.unscaled_rhs_norm;
    computed.raw.unscaled_cost_norm = packet.unscaled_cost_norm;
    computed.scalar_count = 7;
    computed.valid = true;
  }
  if (computed.fixed_point_valid) {
    computed.primal_delta_norm = packet.primal_delta_norm;
    computed.dual_delta_norm = packet.dual_delta_norm;
    computed.delta_dot = packet.delta_dot;
    computed.scalar_count += 3;
  }
  if (computed.anchor_distances_valid) {
    computed.primal_anchor_distance = packet.primal_anchor_distance;
    computed.dual_anchor_distance = packet.dual_anchor_distance;
    computed.scalar_count += 2;
  }
  result = computed;
  return api::Error();
}

api::Error CudaBackend::createVector(const HighsInt size,
                                     VectorHandle& vector) {
  if (size < 0) return invalidArgument("vector size must be nonnegative");
  try {
    std::unique_ptr<CudaVectorResource> resource(new CudaVectorResource(this));
    resource->size = size;
    api::Error error = resource->values.allocate(static_cast<std::size_t>(size),
                                                 "creating vector storage");
    if (!error.ok()) return error;
    if (resource->values.bytes() > 0) {
      const cudaError_t status =
          cudaMemsetAsync(resource->values.get(), 0, resource->values.bytes(),
                          impl_->stream.get());
      if (status != cudaSuccess)
        return cudaFailure("vector initialization", status);
    }
    error = resource->descriptor.create(size, resource->values.get());
    if (!error.ok()) return error;
    if (static_cast<std::size_t>(size) > impl_->difference_scratch.size()) {
      DeviceBuffer<double> scratch;
      error = scratch.allocate(static_cast<std::size_t>(size),
                                "preparing difference reduction scratch");
      if (!error.ok()) return error;
      impl_->difference_scratch = std::move(scratch);
    }
    vector = VectorHandle(std::unique_ptr<Resource>(resource.release()));
  } catch (const std::bad_alloc&) {
    return allocationFailure("creating a vector");
  }
  return api::Error();
}

api::Error CudaBackend::uploadVector(const std::vector<double>& values,
                                     VectorHandle& vector) {
  CudaVectorResource* cuda_vector = vectorResource(vector);
  if (cuda_vector == nullptr ||
      static_cast<std::size_t>(cuda_vector->size) != values.size())
    return invalidArgument("HiPDLP CUDA upload vector dimension mismatch");
  if (!finiteVector(values))
    return invalidArgument("HiPDLP CUDA upload contains a non-finite value");
  return upload(values, cuda_vector->values, impl_->stream.get(),
                "vector upload");
}

api::Error CudaBackend::initialize(VectorHandle& primal, VectorHandle& dual) {
  api::Error error = fill(primal, 0.0);
  if (!error.ok()) return error;
  return fill(dual, 0.0);
}

api::Error CudaBackend::fill(VectorHandle& vector, const double value) {
  CudaVectorResource* cuda_vector = vectorResource(vector);
  if (cuda_vector == nullptr)
    return invalidArgument("HiPDLP CUDA vector handle is invalid");
  if (!std::isfinite(value))
    return invalidArgument("HiPDLP CUDA fill value must be finite");
  impl_->recordDependency(cuda_vector->lifetime);
  return launchResult("fill kernel launch",
                      launchFill(impl_->stream.get(), cuda_vector->size, value,
                                 cuda_vector->values.get()));
}

api::Error CudaBackend::copy(const VectorHandle& source,
                             VectorHandle& destination) {
  const CudaVectorResource* cuda_source = vectorResource(source);
  CudaVectorResource* cuda_destination = vectorResource(destination);
  if (cuda_source == nullptr || cuda_destination == nullptr ||
      cuda_source->size != cuda_destination->size)
    return invalidArgument("HiPDLP CUDA copy vector dimension mismatch");
  impl_->recordDependency(cuda_source->lifetime);
  impl_->recordDependency(cuda_destination->lifetime);
  if (cuda_source == cuda_destination || cuda_source->size == 0)
    return api::Error();
  if (cuda_source->size <= std::numeric_limits<int>::max()) {
    HIPDLP_CUDA_RECORD_PRIMITIVE();
    const cublasStatus_t status = cublasDcopy(
        impl_->cublas.get(), static_cast<int>(cuda_source->size),
        cuda_source->values.get(), 1, cuda_destination->values.get(), 1);
    if (status != CUBLAS_STATUS_SUCCESS)
      return cublasFailure("vector copy", status);
    HIPDLP_CUDA_AUDIT(device_to_device_transfers);
    HIPDLP_CUDA_AUDIT_ADD(device_to_device_bytes, cuda_source->values.bytes());
    return api::Error();
  }
  const cudaError_t status =
      cudaMemcpyAsync(cuda_destination->values.get(), cuda_source->values.get(),
                      cuda_source->values.bytes(), cudaMemcpyDeviceToDevice,
                      impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("vector copy", status);
  HIPDLP_CUDA_AUDIT(device_to_device_transfers);
  HIPDLP_CUDA_AUDIT_ADD(device_to_device_bytes, cuda_source->values.bytes());
  return api::Error();
}

api::Error CudaBackend::scale(const double factor, VectorHandle& vector) {
  CudaVectorResource* cuda_vector = vectorResource(vector);
  if (cuda_vector == nullptr || !std::isfinite(factor))
    return invalidArgument("HiPDLP CUDA scale arguments are invalid");
  impl_->recordDependency(cuda_vector->lifetime);
  return launchResult(
      "scale kernel launch",
      launchScale(impl_->stream.get(), cuda_vector->size, factor,
                  cuda_vector->values.get(), impl_->error_flag.get()));
}

api::Error CudaBackend::axpy(const double factor, const VectorHandle& input,
                             VectorHandle& output) {
  const CudaVectorResource* cuda_input = vectorResource(input);
  CudaVectorResource* cuda_output = vectorResource(output);
  if (cuda_input == nullptr || cuda_output == nullptr ||
      cuda_input->size != cuda_output->size || !std::isfinite(factor))
    return invalidArgument("HiPDLP CUDA axpy arguments are invalid");
  impl_->recordDependency(cuda_input->lifetime);
  impl_->recordDependency(cuda_output->lifetime);
  return launchResult(
      "axpy kernel launch",
      launchAxpy(impl_->stream.get(), cuda_input->size, factor,
                 cuda_input->values.get(), cuda_output->values.get(),
                 impl_->error_flag.get()));
}

api::Error CudaBackend::difference(const VectorHandle& left,
                                   const VectorHandle& right,
                                   VectorHandle& output) {
  const CudaVectorResource* first = vectorResource(left);
  const CudaVectorResource* second = vectorResource(right);
  CudaVectorResource* destination = vectorResource(output);
  if (first == nullptr || second == nullptr || destination == nullptr ||
      first->size != second->size || first->size != destination->size ||
      destination == first || destination == second)
    return invalidArgument("HiPDLP CUDA difference arguments are invalid");
  impl_->recordDependency(first->lifetime);
  impl_->recordDependency(second->lifetime);
  impl_->recordDependency(destination->lifetime);
  return launchResult("difference kernel launch",
      launchDifference(impl_->stream.get(), first->size, first->values.get(),
                       second->values.get(), destination->values.get(),
                       impl_->error_flag.get()));
}

api::Error CudaBackend::multiply(const ProblemHandle& problem,
                                 const VectorHandle& input,
                                 VectorHandle& output) {
  const CudaProblemResource* cuda_problem = problemResource(problem);
  const CudaVectorResource* cuda_input = vectorResource(input);
  CudaVectorResource* cuda_output = vectorResource(output);
  if (cuda_problem == nullptr || cuda_input == nullptr ||
      cuda_output == nullptr || cuda_input->size != cuda_problem->num_cols ||
      cuda_output->size != cuda_problem->num_rows)
    return invalidArgument("HiPDLP CUDA Ax dimensions are inconsistent");
  if (cuda_input == cuda_output)
    return invalidArgument("HiPDLP CUDA Ax does not support aliased vectors");
  impl_->recordDependency(cuda_problem->lifetime);
  impl_->recordDependency(cuda_input->lifetime);
  impl_->recordDependency(cuda_output->lifetime);
  return executeSpmv(*impl_, cuda_problem->matrix, *cuda_input,
                     *cuda_output, impl_->error_flag.get());
}

api::Error CudaBackend::multiplyTranspose(const ProblemHandle& problem,
                                          const VectorHandle& input,
                                          VectorHandle& output) {
  const CudaProblemResource* cuda_problem = problemResource(problem);
  const CudaVectorResource* cuda_input = vectorResource(input);
  CudaVectorResource* cuda_output = vectorResource(output);
  if (cuda_problem == nullptr || cuda_input == nullptr ||
      cuda_output == nullptr || cuda_input->size != cuda_problem->num_rows ||
      cuda_output->size != cuda_problem->num_cols)
    return invalidArgument("HiPDLP CUDA A^T y dimensions are inconsistent");
  if (cuda_input == cuda_output)
    return invalidArgument(
        "HiPDLP CUDA A^T y does not support aliased vectors");
  impl_->recordDependency(cuda_problem->lifetime);
  impl_->recordDependency(cuda_input->lifetime);
  impl_->recordDependency(cuda_output->lifetime);
  return executeSpmv(*impl_, cuda_problem->transpose,
                     *cuda_input, *cuda_output, impl_->error_flag.get());
}

api::Error CudaBackend::fusedProjectedPrimalUpdate(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current_primal, const VectorHandle& transpose_product,
    const double step_size, VectorHandle& next_primal) {
  const CudaProblemResource* cuda_problem = problemResource(problem);
  CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  const CudaVectorResource* current = vectorResource(current_primal);
  const CudaVectorResource* product = vectorResource(transpose_product);
  CudaVectorResource* next = vectorResource(next_primal);
  if (cuda_problem == nullptr || cuda_workspace == nullptr ||
      current == nullptr || product == nullptr || next == nullptr ||
      !std::isfinite(step_size) || step_size < 0.0 ||
      cuda_workspace->num_cols != cuda_problem->num_cols ||
      current->size != cuda_problem->num_cols ||
      product->size != current->size || next->size != current->size)
    return invalidArgument(
        "HiPDLP CUDA fused primal update arguments are invalid");
  impl_->recordDependency(cuda_problem->lifetime);
  impl_->recordDependency(cuda_workspace->lifetime);
  impl_->recordDependency(current->lifetime);
  impl_->recordDependency(product->lifetime);
  impl_->recordDependency(next->lifetime);
  return launchResult(
      "projected primal update kernel launch",
      launchProjectedPrimalUpdate(
          impl_->stream.get(), current->size, current->values.get(),
          product->values.get(), cuda_problem->costs.get(),
          cuda_problem->col_lower.get(), cuda_problem->col_upper.get(),
          step_size, next->values.get(), impl_->error_flag.get()));
}

api::Error CudaBackend::fusedProjectedDualUpdate(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current_dual, const VectorHandle& extrapolated_product,
    const double step_size, VectorHandle& next_dual) {
  const CudaProblemResource* cuda_problem = problemResource(problem);
  CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  const CudaVectorResource* current = vectorResource(current_dual);
  const CudaVectorResource* product = vectorResource(extrapolated_product);
  CudaVectorResource* next = vectorResource(next_dual);
  if (cuda_problem == nullptr || cuda_workspace == nullptr ||
      current == nullptr || product == nullptr || next == nullptr ||
      !std::isfinite(step_size) || step_size < 0.0 ||
      cuda_workspace->num_rows != cuda_problem->num_rows ||
      current->size != cuda_problem->num_rows ||
      product->size != current->size || next->size != current->size)
    return invalidArgument(
        "HiPDLP CUDA fused dual update arguments are invalid");
  impl_->recordDependency(cuda_problem->lifetime);
  impl_->recordDependency(cuda_workspace->lifetime);
  impl_->recordDependency(current->lifetime);
  impl_->recordDependency(product->lifetime);
  impl_->recordDependency(next->lifetime);
  return launchResult(
      "projected dual update kernel launch",
      launchProjectedDualUpdate(impl_->stream.get(), current->size,
                                current->values.get(), product->values.get(),
                                cuda_problem->row_lower.get(),
                                cuda_problem->equality_rows.get(), step_size,
                                next->values.get(), impl_->error_flag.get()));
}

api::Error CudaBackend::fusedProjectedPrimalUpdateAndReflect(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current, const VectorHandle& product,
    const double step_size, VectorHandle& candidate, VectorHandle& reflected) {
  return fusedProjectedPrimalUpdateAndReflect(problem, workspace, current,
      product, ScalarParameter(step_size), candidate, reflected);
}

api::Error CudaBackend::scalarParameter(const ScalarParameter& parameter,
                                       const double*& device_value) {
  device_value = nullptr;
  if (parameter.vector == nullptr) return api::Error();
  const CudaVectorResource* values = vectorResource(*parameter.vector);
  if (values == nullptr || parameter.index < 0 || parameter.index >= values->size)
    return invalidArgument("HiPDLP CUDA scalar parameter owner or index is invalid");
  impl_->recordDependency(values->lifetime);
  device_value = values->values.get() + parameter.index;
  return api::Error();
}

api::Error CudaBackend::fusedProjectedPrimalUpdateAndReflect(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current, const VectorHandle& product,
    const ScalarParameter& step, VectorHandle& candidate, VectorHandle& reflected) {
  const double step_size = step.value;
  const api::Error handles = validateReflectionHandles(current, product,
                                                       candidate, reflected);
  if (!handles.ok()) return handles;
  const CudaProblemResource* cuda_problem = problemResource(problem);
  const CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  const CudaVectorResource* cuda_current = vectorResource(current);
  const CudaVectorResource* cuda_product = vectorResource(product);
  CudaVectorResource* cuda_candidate = vectorResource(candidate);
  CudaVectorResource* cuda_reflected = vectorResource(reflected);
  if (!cuda_problem || !cuda_workspace || !cuda_current || !cuda_product ||
      !cuda_candidate || !cuda_reflected || !std::isfinite(step_size) ||
      step_size < 0.0 || cuda_workspace->num_cols != cuda_problem->num_cols ||
      cuda_current->size != cuda_problem->num_cols ||
      cuda_product->size != cuda_current->size ||
      cuda_candidate->size != cuda_current->size ||
      cuda_reflected->size != cuda_current->size)
    return invalidArgument("HiPDLP CUDA primal projection/reflection arguments are invalid");
  if (step.vector == &candidate || step.vector == &reflected)
    return invalidArgument("HiPDLP CUDA scalar parameter aliases an output");
  const double* device_step = nullptr;
  const api::Error parameter_error = scalarParameter(step, device_step);
  if (!parameter_error.ok()) return parameter_error;
  impl_->recordDependency(cuda_problem->lifetime);
  impl_->recordDependency(cuda_workspace->lifetime);
  impl_->recordDependency(cuda_current->lifetime);
  impl_->recordDependency(cuda_product->lifetime);
  impl_->recordDependency(cuda_candidate->lifetime);
  impl_->recordDependency(cuda_reflected->lifetime);
  return launchResult(
      "primal projection/reflection kernel launch",
      launchProjectedPrimalUpdateAndReflect(
          impl_->stream.get(), cuda_current->size, cuda_current->values.get(),
          cuda_product->values.get(), cuda_problem->costs.get(),
          cuda_problem->col_lower.get(), cuda_problem->col_upper.get(),
          step_size, cuda_candidate->values.get(), cuda_reflected->values.get(),
          impl_->error_flag.get(), device_step));
}

api::Error CudaBackend::fusedProjectedDualUpdateAndReflect(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current, const VectorHandle& product,
    const double step_size, VectorHandle& candidate, VectorHandle& reflected) {
  return fusedProjectedDualUpdateAndReflect(problem, workspace, current,
      product, ScalarParameter(step_size), candidate, reflected);
}

api::Error CudaBackend::fusedProjectedDualUpdateAndReflect(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& current, const VectorHandle& product,
    const ScalarParameter& step, VectorHandle& candidate, VectorHandle& reflected) {
  const double step_size = step.value;
  const api::Error handles = validateReflectionHandles(current, product,
                                                       candidate, reflected);
  if (!handles.ok()) return handles;
  const CudaProblemResource* cuda_problem = problemResource(problem);
  const CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  const CudaVectorResource* cuda_current = vectorResource(current);
  const CudaVectorResource* cuda_product = vectorResource(product);
  CudaVectorResource* cuda_candidate = vectorResource(candidate);
  CudaVectorResource* cuda_reflected = vectorResource(reflected);
  if (!cuda_problem || !cuda_workspace || !cuda_current || !cuda_product ||
      !cuda_candidate || !cuda_reflected || !std::isfinite(step_size) ||
      step_size < 0.0 || cuda_workspace->num_rows != cuda_problem->num_rows ||
      cuda_current->size != cuda_problem->num_rows ||
      cuda_product->size != cuda_current->size ||
      cuda_candidate->size != cuda_current->size ||
      cuda_reflected->size != cuda_current->size)
    return invalidArgument("HiPDLP CUDA dual projection/reflection arguments are invalid");
  if (step.vector == &candidate || step.vector == &reflected)
    return invalidArgument("HiPDLP CUDA scalar parameter aliases an output");
  const double* device_step = nullptr;
  const api::Error parameter_error = scalarParameter(step, device_step);
  if (!parameter_error.ok()) return parameter_error;
  impl_->recordDependency(cuda_problem->lifetime);
  impl_->recordDependency(cuda_workspace->lifetime);
  impl_->recordDependency(cuda_current->lifetime);
  impl_->recordDependency(cuda_product->lifetime);
  impl_->recordDependency(cuda_candidate->lifetime);
  impl_->recordDependency(cuda_reflected->lifetime);
  return launchResult(
      "dual projection/reflection kernel launch",
      launchProjectedDualUpdateAndReflect(
          impl_->stream.get(), cuda_current->size, cuda_current->values.get(),
          cuda_product->values.get(), cuda_problem->row_lower.get(),
          cuda_problem->equality_rows.get(), step_size,
          cuda_candidate->values.get(), cuda_reflected->values.get(),
          impl_->error_flag.get(), device_step));
}

api::Error CudaBackend::multiplyAndProjectedDualUpdateAndReflect(
    const ProblemHandle& problem, WorkspaceHandle& workspace,
    const VectorHandle& input, const VectorHandle& current,
    const ScalarParameter& step, VectorHandle& product,
    VectorHandle& candidate, VectorHandle& reflected) {
  const api::Error handles = validateReflectionHandles(current, product,
                                                       candidate, reflected);
  if (!handles.ok()) return handles;
  const CudaProblemResource* cuda_problem = problemResource(problem);
  const CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  const CudaVectorResource* cuda_input = vectorResource(input);
  const CudaVectorResource* cuda_current = vectorResource(current);
  CudaVectorResource* cuda_product = vectorResource(product);
  CudaVectorResource* cuda_candidate = vectorResource(candidate);
  CudaVectorResource* cuda_reflected = vectorResource(reflected);
  if (!cuda_problem || !cuda_workspace || !cuda_input || !cuda_current ||
      !cuda_product || !cuda_candidate || !cuda_reflected ||
      cuda_input == cuda_product || !std::isfinite(step.value) || step.value < 0.0 ||
      cuda_input->size != cuda_problem->num_cols ||
      cuda_workspace->num_rows != cuda_problem->num_rows ||
      cuda_current->size != cuda_problem->num_rows ||
      cuda_product->size != cuda_current->size ||
      cuda_candidate->size != cuda_current->size ||
      cuda_reflected->size != cuda_current->size ||
      step.vector == &product || step.vector == &candidate || step.vector == &reflected)
    return invalidArgument("HiPDLP CUDA product/projection arguments are invalid");
  const double* device_step = nullptr;
  api::Error error = scalarParameter(step, device_step);
  if (!error.ok()) return error;
  impl_->recordDependency(cuda_problem->lifetime);
  impl_->recordDependency(cuda_workspace->lifetime);
  impl_->recordDependency(cuda_input->lifetime);
  impl_->recordDependency(cuda_current->lifetime);
  impl_->recordDependency(cuda_product->lifetime);
  impl_->recordDependency(cuda_candidate->lifetime);
  impl_->recordDependency(cuda_reflected->lifetime);
  error = executeSpmv(*impl_, cuda_problem->matrix, *cuda_input,
                      *cuda_product, impl_->error_flag.get(), false);
  if (!error.ok()) return error;
  return launchResult("product/projection/reflection kernel launch",
      launchProjectedDualUpdateAndReflect(impl_->stream.get(), cuda_current->size,
          cuda_current->values.get(), cuda_product->values.get(),
          cuda_problem->row_lower.get(), cuda_problem->equality_rows.get(),
          step.value, cuda_candidate->values.get(), cuda_reflected->values.get(),
          impl_->error_flag.get(), device_step));
}

api::Error CudaBackend::fusedVectorBlend(
    WorkspaceHandle& workspace, const VectorHandle& primal_left,
    const VectorHandle& dual_left, const VectorHandle& primal_right,
    const VectorHandle& dual_right, const double left_weight,
    VectorHandle& primal_output, VectorHandle& dual_output) {
  return fusedVectorBlend(workspace, primal_left, dual_left, primal_right,
                           dual_right, ScalarParameter(left_weight),
                           primal_output, dual_output);
}

api::Error CudaBackend::fusedVectorBlend(
    WorkspaceHandle& workspace, const VectorHandle& primal_left,
    const VectorHandle& dual_left, const VectorHandle& primal_right,
    const VectorHandle& dual_right, const ScalarParameter& weight,
    VectorHandle& primal_output, VectorHandle& dual_output) {
  const double left_weight = weight.value;
  CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  const CudaVectorResource* p_left = vectorResource(primal_left);
  const CudaVectorResource* d_left = vectorResource(dual_left);
  const CudaVectorResource* p_right = vectorResource(primal_right);
  const CudaVectorResource* d_right = vectorResource(dual_right);
  CudaVectorResource* p_output = vectorResource(primal_output);
  CudaVectorResource* d_output = vectorResource(dual_output);
  if (cuda_workspace == nullptr || p_left == nullptr || d_left == nullptr ||
      p_right == nullptr || d_right == nullptr || p_output == nullptr ||
      d_output == nullptr || !std::isfinite(left_weight) || left_weight < 0.0 ||
      left_weight > 1.0 || p_left->size != p_right->size ||
      p_left->size != p_output->size || d_left->size != d_right->size ||
      d_left->size != d_output->size ||
      p_left->size != cuda_workspace->num_cols ||
      d_left->size != cuda_workspace->num_rows)
    return invalidArgument(
        "HiPDLP CUDA fused vector blend arguments are invalid");
  if (weight.vector == &primal_output || weight.vector == &dual_output)
    return invalidArgument("HiPDLP CUDA scalar parameter aliases an output");
  const double* device_weight = nullptr;
  const api::Error parameter_error = scalarParameter(weight, device_weight);
  if (!parameter_error.ok()) return parameter_error;
  impl_->recordDependency(cuda_workspace->lifetime);
  impl_->recordDependency(p_left->lifetime);
  impl_->recordDependency(d_left->lifetime);
  impl_->recordDependency(p_right->lifetime);
  impl_->recordDependency(d_right->lifetime);
  impl_->recordDependency(p_output->lifetime);
  impl_->recordDependency(d_output->lifetime);
  return launchResult(
      "paired vector blend kernel launch",
      launchVectorPairBlend(impl_->stream.get(), p_left->size, d_left->size,
        p_left->values.get(), p_right->values.get(), d_left->values.get(),
        d_right->values.get(), left_weight, p_output->values.get(),
        d_output->values.get(), impl_->error_flag.get(), device_weight));
}

api::Error CudaBackend::fusedAverageUpdate(WorkspaceHandle& workspace,
                                           const VectorHandle& primal,
                                           const VectorHandle& dual,
                                           const double weight,
                                           VectorHandle& primal_sum,
                                           VectorHandle& dual_sum) {
  CudaWorkspaceResource* cuda_workspace = workspaceResource(workspace);
  const CudaVectorResource* p = vectorResource(primal);
  const CudaVectorResource* d = vectorResource(dual);
  CudaVectorResource* p_sum = vectorResource(primal_sum);
  CudaVectorResource* d_sum = vectorResource(dual_sum);
  if (cuda_workspace == nullptr || p == nullptr || d == nullptr ||
      p_sum == nullptr || d_sum == nullptr || !std::isfinite(weight) ||
      p->size != p_sum->size || d->size != d_sum->size ||
      p->size != cuda_workspace->num_cols ||
      d->size != cuda_workspace->num_rows)
    return invalidArgument(
        "HiPDLP CUDA fused average update arguments are invalid");
  api::Error error = launchResult(
      "primal average update kernel launch",
      launchAverageUpdate(impl_->stream.get(), p->size, p->values.get(), weight,
                          p_sum->values.get(), impl_->error_flag.get()));
  if (!error.ok()) return error;
  return launchResult(
      "dual average update kernel launch",
      launchAverageUpdate(impl_->stream.get(), d->size, d->values.get(), weight,
                          d_sum->values.get(), impl_->error_flag.get()));
}

api::Error CudaBackend::dot(const VectorHandle& left, const VectorHandle& right,
                            double& result) {
  const CudaVectorResource* cuda_left = vectorResource(left);
  const CudaVectorResource* cuda_right = vectorResource(right);
  if (cuda_left == nullptr || cuda_right == nullptr ||
      cuda_left->size != cuda_right->size)
    return invalidArgument("HiPDLP CUDA dot dimensions are inconsistent");
  api::Error error = enqueueDot(*impl_, cuda_left->size, cuda_left->values.get(),
                                cuda_right->values.get(), impl_->reduction_result.get());
  if (!error.ok()) return error;
  error = launchResult("dot finite-check kernel launch",
        launchFiniteCheck(impl_->stream.get(), 1, impl_->reduction_result.get(),
                           impl_->error_flag.get()));
  if (!error.ok()) return error;
  double reduced = 0.0;
  const cudaError_t status =
      cudaMemcpyAsync(&reduced, impl_->reduction_result.get(), sizeof(reduced),
                      cudaMemcpyDeviceToHost, impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("dot result download", status);
  HIPDLP_CUDA_AUDIT(device_to_host_transfers);
  HIPDLP_CUDA_AUDIT_ADD(device_to_host_bytes, sizeof(reduced));
  error = impl_->synchronizeAndCheck();
  if (!error.ok()) return error;
  result = reduced;
  return api::Error();
}

api::Error CudaBackend::norm2(const VectorHandle& vector, double& result) {
  const CudaVectorResource* cuda_vector = vectorResource(vector);
  if (cuda_vector == nullptr)
    return invalidArgument("HiPDLP CUDA norm vector is invalid");
  api::Error error = enqueueNorm(*impl_, cuda_vector->size, cuda_vector->values.get(),
                                 impl_->reduction_result.get());
  if (!error.ok()) return error;
  error = launchResult("norm finite-check kernel launch",
        launchFiniteCheck(impl_->stream.get(), 1, impl_->reduction_result.get(),
                           impl_->error_flag.get()));
  if (!error.ok()) return error;
  double reduced = 0.0;
  const cudaError_t status =
      cudaMemcpyAsync(&reduced, impl_->reduction_result.get(), sizeof(reduced),
                      cudaMemcpyDeviceToHost, impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("norm result download", status);
  HIPDLP_CUDA_AUDIT(device_to_host_transfers);
  HIPDLP_CUDA_AUDIT_ADD(device_to_host_bytes, sizeof(reduced));
  error = impl_->synchronizeAndCheck();
  if (!error.ok()) return error;
  result = reduced;
  return api::Error();
}

api::Error CudaBackend::differenceNorm2(const VectorHandle& left,
                                        const VectorHandle& right,
                                        double& result) {
  const CudaVectorResource* cuda_left = vectorResource(left);
  const CudaVectorResource* cuda_right = vectorResource(right);
  if (cuda_left == nullptr || cuda_right == nullptr ||
      cuda_left->size != cuda_right->size)
    return invalidArgument(
        "HiPDLP CUDA difference norm dimensions are inconsistent");
  api::Error error;
  if (cuda_left->size > 0) {
    error = launchResult("difference kernel launch",
        launchDifference(impl_->stream.get(), cuda_left->size,
                         cuda_left->values.get(), cuda_right->values.get(),
                         impl_->difference_scratch.get(), impl_->error_flag.get()));
    if (!error.ok()) return error;
    error = enqueueNorm(*impl_, cuda_left->size, impl_->difference_scratch.get(),
               impl_->reduction_result.get());
    if (!error.ok()) return error;
    error = launchResult("difference norm finite-check kernel launch",
        launchFiniteCheck(impl_->stream.get(), 1, impl_->reduction_result.get(),
                           impl_->error_flag.get()));
  } else {
    error = launchResult("difference norm reduction kernel launch",
        launchDifferenceNorm2(impl_->stream.get(), cuda_left->size,
                              cuda_left->values.get(), cuda_right->values.get(),
                              impl_->reduction_result.get(), impl_->error_flag.get()));
  }
  if (!error.ok()) return error;
  double reduced = 0.0;
  const cudaError_t status =
      cudaMemcpyAsync(&reduced, impl_->reduction_result.get(), sizeof(reduced),
                      cudaMemcpyDeviceToHost, impl_->stream.get());
  if (status != cudaSuccess)
    return cudaFailure("difference norm result download", status);
  HIPDLP_CUDA_AUDIT(device_to_host_transfers);
  HIPDLP_CUDA_AUDIT_ADD(device_to_host_bytes, sizeof(reduced));
  error = impl_->synchronizeAndCheck();
  if (!error.ok()) return error;
  result = reduced;
  return api::Error();
}

api::Error CudaBackend::reduceFixedPointMetrics(
    const VectorHandle& delta_primal, const VectorHandle& delta_dual,
    const VectorHandle& aty_delta, double& primal_norm, double& dual_norm,
    double& cross_term) {
  if (impl_->recording) return unsupported("host reductions during recording");
  const CudaVectorResource* primal = vectorResource(delta_primal);
  const CudaVectorResource* dual = vectorResource(delta_dual);
  const CudaVectorResource* product = vectorResource(aty_delta);
  if (primal == nullptr || dual == nullptr || product == nullptr ||
      primal->size != product->size || &primal_norm == &dual_norm ||
      &primal_norm == &cross_term || &dual_norm == &cross_term)
    return invalidArgument("HiPDLP CUDA fixed-point reduction arguments are invalid");
  double* values = impl_->reduction_result.get();
  api::Error error = launchResult("fused fixed-point reduction",
      launchFixedPointReduction(impl_->stream.get(), primal->size, dual->size,
          primal->values.get(), dual->values.get(), product->values.get(),
          impl_->fixed_point_partials.get(), values, impl_->error_flag.get()));
  if (!error.ok()) return error;
  double received[3] = {};
  const cudaError_t status = cudaMemcpyAsync(received, values, sizeof(received),
                                             cudaMemcpyDeviceToHost, impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("fixed-point result download", status);
  HIPDLP_CUDA_AUDIT(device_to_host_transfers);
  HIPDLP_CUDA_AUDIT_ADD(device_to_host_bytes, sizeof(received));
  error = impl_->synchronizeAndCheck();
  if (!error.ok()) return error;
  primal_norm = received[0];
  dual_norm = received[1];
  cross_term = received[2];
  return api::Error();
}

api::Error CudaBackend::reduceConvergenceMetrics(
    const VectorHandle& primal_residual, const VectorHandle& dual_residual,
    const double primal_objective, const double dual_objective,
    ConvergenceMetrics& metrics) {
  if (!std::isfinite(primal_objective) || !std::isfinite(dual_objective))
    return invalidArgument("convergence objectives must be finite");
  if (impl_->recording) return unsupported("host reductions during recording");
  const CudaVectorResource* primal = vectorResource(primal_residual);
  const CudaVectorResource* dual = vectorResource(dual_residual);
  if (primal == nullptr || dual == nullptr)
    return invalidArgument("HiPDLP CUDA convergence reduction arguments are invalid");
  double* values = impl_->reduction_result.get();
  api::Error error = enqueueNorm(*impl_, primal->size, primal->values.get(), values);
  if (!error.ok()) return error;
  error = enqueueNorm(*impl_, dual->size, dual->values.get(), values + 1);
  if (!error.ok()) return error;
  error = launchResult("convergence result finite check",
      launchFiniteCheck(impl_->stream.get(), 2, values, impl_->error_flag.get()));
  if (!error.ok()) return error;
  double received[2] = {};
  const cudaError_t status = cudaMemcpyAsync(received, values, sizeof(received),
                                             cudaMemcpyDeviceToHost, impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("convergence result download", status);
  HIPDLP_CUDA_AUDIT(device_to_host_transfers);
  HIPDLP_CUDA_AUDIT_ADD(device_to_host_bytes, sizeof(received));
  error = impl_->synchronizeAndCheck();
  if (!error.ok()) return error;
  const double gap = std::abs(primal_objective - dual_objective);
  if (!std::isfinite(gap))
    return backendFailure("HiPDLP CUDA convergence gap is non-finite");
  metrics.primal_residual = received[0];
  metrics.dual_residual = received[1];
  metrics.objective_gap = gap;
  return api::Error();
}

api::Error CudaBackend::beginRecording() {
  if (impl_->timing_active)
    return invalidArgument("HiPDLP CUDA backend timing is active");
  if (impl_->recording)
    return invalidArgument("HiPDLP CUDA command recording is already active");
  if (HIPDLP_CUDA_FAILURE(kGraphCapture))
    return backendFailure(
        "HiPDLP CUDA command recording failed (injected failure)");
  api::Error error = impl_->synchronizeAndCheck();
  if (!error.ok()) return error;
  try {
    impl_->recording_dependencies.clear();
    impl_->recording_dependencies.reserve(128);
  } catch (const std::bad_alloc&) {
    return allocationFailure("preparing command recording dependencies");
  }
  HIPDLP_CUDA_AUDIT(graph_capture_calls);
  const cudaError_t status =
      cudaStreamBeginCapture(impl_->stream.get(), cudaStreamCaptureModeGlobal);
  if (status != cudaSuccess)
    return cudaFailure("command recording begin", status);
  impl_->recording = true;
  HIPDLP_CUDA_CAPTURE_ACTIVE(true);
#ifdef HIGHS_HIPDLP_TESTING
  if (HIPDLP_CUDA_FAILURE(kGraphAsynchronousExecution)) {
    error = launchResult(
        "injected asynchronous failure kernel launch",
        launchInjectedFailure(impl_->stream.get(), impl_->error_flag.get()));
    if (!error.ok()) {
      cancelRecording();
      return error;
    }
  }
#endif
  return api::Error();
}

api::Error CudaBackend::endRecording(CommandBufferHandle& command_buffer) {
  command_buffer.reset();
  if (!impl_->recording)
    return invalidArgument("HiPDLP CUDA command recording is not active");
  std::unique_ptr<CudaCommandBufferResource> resource;
  try {
    resource.reset(new CudaCommandBufferResource());
  } catch (const std::bad_alloc&) {
    cancelRecording();
    return allocationFailure("creating a command buffer");
  }

  cudaGraph_t captured_graph = nullptr;
  cudaError_t status =
      cudaStreamEndCapture(impl_->stream.get(), &captured_graph);
  impl_->recording = false;
  HIPDLP_CUDA_CAPTURE_ACTIVE(false);
  if (captured_graph != nullptr) resource->graph.adopt(captured_graph);
  if (status != cudaSuccess)
    return cudaFailure("command recording end", status);
  if (resource->graph.get() == nullptr)
    return backendFailure("HiPDLP CUDA command recording produced no work");

  HIPDLP_CUDA_AUDIT(graph_instantiation_calls);
  if (HIPDLP_CUDA_FAILURE(kGraphInstantiation))
    return backendFailure(
        "HiPDLP CUDA command finalization failed (injected failure)");
  cudaGraphExec_t executable = nullptr;
  status = cudaGraphInstantiate(&executable, resource->graph.get(), nullptr,
                                nullptr, 0);
  if (status != cudaSuccess)
    return cudaFailure("command buffer finalization", status);
  if (executable == nullptr)
    return backendFailure(
        "HiPDLP CUDA command finalization produced no executable");
  resource->executable.adopt(executable);
  resource->dependencies.swap(impl_->recording_dependencies);
  command_buffer =
      CommandBufferHandle(std::unique_ptr<Resource>(resource.release()));
  return api::Error();
}

api::Error CudaBackend::cancelRecording() {
  if (!impl_->recording)
    return invalidArgument("HiPDLP CUDA command recording is not active");
  cudaGraph_t captured_graph = nullptr;
  const cudaError_t status =
      cudaStreamEndCapture(impl_->stream.get(), &captured_graph);
  impl_->recording = false;
  impl_->recording_dependencies.clear();
  HIPDLP_CUDA_CAPTURE_ACTIVE(false);
  if (captured_graph != nullptr) cudaGraphDestroy(captured_graph);
  return status == cudaSuccess
             ? api::Error()
             : cudaFailure("command recording cancel", status);
}

api::Error CudaBackend::replay(const CommandBufferHandle& command_buffer,
                               CommandBufferStatistics& statistics) {
  statistics = CommandBufferStatistics();
  const CudaCommandBufferResource* resource =
      commandBufferResource(command_buffer);
  if (resource == nullptr)
    return invalidArgument("HiPDLP CUDA command buffer is invalid");
  if (!dependenciesAlive(*resource))
    return invalidArgument("HiPDLP CUDA command buffer dependency has expired");
  if (HIPDLP_CUDA_FAILURE(kGraphLaunch))
    return backendFailure(
        "HiPDLP CUDA command-buffer replay failed (injected failure)");
  HIPDLP_CUDA_AUDIT(graph_launch_calls);
  const cudaError_t status =
      cudaGraphLaunch(resource->executable.get(), impl_->stream.get());
  if (status != cudaSuccess)
    return cudaFailure("command-buffer replay", status);
  statistics.replay_calls = 1;
  return api::Error();
}

api::Error CudaBackend::prepareDownloadBuffers(const VectorDownload* downloads,
                                               const std::size_t count) {
  api::Error error = validateVectorDownloads(downloads, count);
  if (!error.ok() || count == 0) return error;
  if (impl_->recording)
    return invalidArgument("HiPDLP CUDA downloads cannot run during recording");
  for (std::size_t index = 0; index < count; ++index) {
    const CudaVectorResource* resource = vectorResource(*downloads[index].vector);
    if (resource == nullptr || resource->size < 0 ||
        resource->values.size() != static_cast<std::size_t>(resource->size))
      return invalidArgument("HiPDLP CUDA download vector or dimensions are invalid");
  }
  error = impl_->download_error.allocate(1);
  if (!error.ok()) return error;
  for (std::size_t index = 0; index < count; ++index) {
    const CudaVectorResource* resource = vectorResource(*downloads[index].vector);
    error = resource->downloaded.allocate(resource->values.size());
    if (!error.ok()) return error;
  }
  return api::Error();
}

api::Error CudaBackend::prepareVectorDownloads(const VectorDownload* downloads,
                                               const std::size_t count) {
  api::Error error = prepareDownloadBuffers(downloads, count);
  if (!error.ok()) return error;
  DownloadOutputs outputs(downloads, count);
  for (std::size_t index = 0; index < count; ++index) {
    error = outputs.reserve(
        index, vectorResource(*downloads[index].vector)->values.size());
    if (!error.ok()) return error;
  }
  outputs.publishReservations();
  return api::Error();
}

api::Error CudaBackend::downloadVectors(const VectorDownload* downloads,
                                        const std::size_t count) {
  api::Error error = prepareDownloadBuffers(downloads, count);
  if (!error.ok() || count == 0) return error;
  DownloadOutputs outputs(downloads, count);
  for (std::size_t index = 0; index < count; ++index) {
    error = outputs.reserve(
        index, vectorResource(*downloads[index].vector)->values.size());
    if (!error.ok()) return error;
  }
  cudaError_t status = cudaSuccess;
  const char* operation = "vector batch download";
  for (std::size_t index = 0; index < count; ++index) {
    const CudaVectorResource* resource = vectorResource(*downloads[index].vector);
    if (resource->values.size() == 0) continue;
    status = cudaMemcpyAsync(resource->downloaded.get(), resource->values.get(),
                             resource->values.bytes(), cudaMemcpyDeviceToHost,
                             impl_->stream.get());
    if (status != cudaSuccess) break;
    HIPDLP_CUDA_AUDIT(device_to_host_transfers);
    HIPDLP_CUDA_AUDIT_ADD(device_to_host_bytes, resource->values.bytes());
  }
  if (status == cudaSuccess) {
    operation = "batch error-state download";
    status = cudaMemcpyAsync(impl_->download_error.get(), impl_->error_flag.get(),
                             sizeof(int), cudaMemcpyDeviceToHost,
                             impl_->stream.get());
    if (status == cudaSuccess) {
      HIPDLP_CUDA_AUDIT(device_to_host_transfers);
      HIPDLP_CUDA_AUDIT_ADD(device_to_host_bytes, sizeof(int));
    }
  }
  HIPDLP_CUDA_AUDIT(synchronization_calls);
  const cudaError_t waited = cudaStreamSynchronize(impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure(operation, status);
  if (waited != cudaSuccess)
    return cudaFailure("batch download synchronization", waited);
  if (*impl_->download_error.get() != 0)
    return backendFailure("HiPDLP CUDA operation produced a non-finite value");
  for (std::size_t index = 0; index < count; ++index) {
    const CudaVectorResource* resource = vectorResource(*downloads[index].vector);
    outputs.publish(index, resource->downloaded.get(), resource->values.size());
  }
  return api::Error();
}

api::Error CudaBackend::downloadVector(const VectorHandle& vector,
                                       std::vector<double>& values) {
  const VectorDownload download = {&vector, &values};
  return downloadVectors(&download, 1);
}

api::Error CudaBackend::downloadSolution(const VectorHandle& primal,
                                         const VectorHandle& dual,
                                         std::vector<double>& primal_values,
                                         std::vector<double>& dual_values) {
  const VectorDownload downloads[] = {{&primal, &primal_values},
                                       {&dual, &dual_values}};
  return downloadVectors(downloads, 2);
}

api::Error CudaBackend::createBatchedVector(const HighsInt size,
                                            const HighsInt batch_size,
                                            VectorHandle& vector) {
  (void)size;
  (void)batch_size;
  (void)vector;
  return unsupported("batched vectors");
}

api::Error CudaBackend::multiplyBatched(const ProblemHandle& problem,
                                        const VectorHandle& input,
                                        VectorHandle& output) {
  (void)problem;
  (void)input;
  (void)output;
  return unsupported("SpMM");
}

api::Error CudaBackend::synchronize() { return impl_->synchronizeAndCheck(); }

bool CudaBackend::timingUsesEvents() const { return true; }

api::Error CudaBackend::prepareTiming() {
  if (impl_->recording || impl_->timing_active)
    return invalidArgument(
        "HiPDLP CUDA timer preparation requires inactive timing and recording");
  return impl_->timing_events.prepare();
}

api::Error CudaBackend::beginTiming() {
  if (impl_->recording || impl_->timing_active ||
      impl_->timing_events.start() == nullptr)
    return invalidArgument("HiPDLP CUDA timer must be prepared and inactive");
  const cudaError_t status =
      cudaEventRecord(impl_->timing_events.start(), impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("timer start recording", status);
  impl_->timing_active = true;
  return api::Error();
}

api::Error CudaBackend::endTiming(double& seconds) {
  if (impl_->recording || !impl_->timing_active)
    return invalidArgument("HiPDLP CUDA timer has not begun or is recording");
  impl_->timing_active = false;
  cudaError_t status =
      cudaEventRecord(impl_->timing_events.stop(), impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("timer stop recording", status);
  status = cudaEventSynchronize(impl_->timing_events.stop());
  if (status != cudaSuccess) return cudaFailure("timer completion", status);
  float milliseconds = 0.0f;
  status = cudaEventElapsedTime(&milliseconds, impl_->timing_events.start(),
                               impl_->timing_events.stop());
  if (status != cudaSuccess) return cudaFailure("timer elapsed time", status);
  const double elapsed = static_cast<double>(milliseconds) * 0.001;
  if (!std::isfinite(elapsed) || elapsed < 0.0)
    return backendFailure("HiPDLP CUDA timer returned invalid elapsed time");
  seconds = elapsed;
  return api::Error();
}

api::Error createCudaBackend(std::unique_ptr<Backend>& backend) {
  backend.reset();
  int device_count = 0;
  cudaError_t cuda_status = cudaGetDeviceCount(&device_count);
  if (cuda_status == cudaSuccess && device_count == 0)
    return unavailable("no compatible NVIDIA device");
  if (unavailableStatus(cuda_status))
    return cudaUnavailable("device query", cuda_status);
  if (cuda_status != cudaSuccess)
    return cudaFailure("device query", cuda_status);

  int selected_device = -1;
  cudaDeviceProp selected_properties;
  for (int device = 0; device < device_count; ++device) {
    cudaDeviceProp properties;
    cuda_status = cudaGetDeviceProperties(&properties, device);
    if (unavailableStatus(cuda_status))
      return cudaUnavailable("device property query", cuda_status);
    if (cuda_status != cudaSuccess)
      return cudaFailure("device property query", cuda_status);
    if (properties.major > 0) {
      selected_device = device;
      selected_properties = properties;
      break;
    }
  }
  if (selected_device < 0) return unavailable("no compatible NVIDIA device");

  cuda_status = cudaSetDevice(selected_device);
  if (unavailableStatus(cuda_status))
    return cudaUnavailable("device selection", cuda_status);
  if (cuda_status != cudaSuccess)
    return cudaFailure("device selection", cuda_status);

  try {
    std::unique_ptr<CudaBackend::Impl> impl(new CudaBackend::Impl());
    impl->device = selected_device;
    impl->device_name = selected_properties.name;
    api::Error error = impl->stream.create();
    if (!error.ok()) return error;
    error = impl->cublas.create(impl->stream.get());
    if (!error.ok()) return error;
    error = impl->cusparse.create(impl->stream.get());
    if (!error.ok()) return error;
    if (HIPDLP_CUDA_FAILURE(kBackendInitialization))
      return backendFailure(
          "HiPDLP CUDA backend initialization failed (injected failure)");
    error = impl->error_flag.allocate(1, "creating error-state storage");
    if (!error.ok()) return error;
    error = impl->reduction_result.allocate(3, "creating reduction storage");
    if (!error.ok()) return error;
    error = impl->fixed_point_partials.allocate(
        checkPartialCount(kHighsIInf, kHighsIInf), "creating fixed-point reduction storage");
    if (!error.ok()) return error;
    error = impl->spmv_alpha.allocate(1, "creating SpMV alpha storage");
    if (!error.ok()) return error;
    error = impl->spmv_beta.allocate(1, "creating SpMV beta storage");
    if (!error.ok()) return error;
    cuda_status = cudaMemsetAsync(impl->error_flag.get(), 0, sizeof(int),
                                  impl->stream.get());
    if (cuda_status != cudaSuccess)
      return cudaFailure("error-state initialization", cuda_status);
    cuda_status = cudaMemsetAsync(impl->reduction_result.get(), 0,
                                  sizeof(double), impl->stream.get());
    if (cuda_status != cudaSuccess)
      return cudaFailure("reduction storage initialization", cuda_status);
    const double alpha = 1.0;
    cuda_status = cudaMemcpyAsync(impl->spmv_alpha.get(), &alpha, sizeof(alpha),
                                  cudaMemcpyHostToDevice, impl->stream.get());
    if (cuda_status != cudaSuccess)
      return cudaFailure("SpMV alpha initialization", cuda_status);
    HIPDLP_CUDA_AUDIT(host_to_device_transfers);
    HIPDLP_CUDA_AUDIT_ADD(host_to_device_bytes, sizeof(alpha));
    const double beta = 0.0;
    cuda_status = cudaMemcpyAsync(impl->spmv_beta.get(), &beta, sizeof(beta),
                                  cudaMemcpyHostToDevice, impl->stream.get());
    if (cuda_status != cudaSuccess)
      return cudaFailure("SpMV beta initialization", cuda_status);
    HIPDLP_CUDA_AUDIT(host_to_device_transfers);
    HIPDLP_CUDA_AUDIT_ADD(host_to_device_bytes, sizeof(beta));
    cuda_status = preloadKernels();
    if (unavailableStatus(cuda_status))
      return cudaUnavailable("kernel preload", cuda_status);
    if (cuda_status != cudaSuccess)
      return cudaFailure("kernel preload", cuda_status);
    cuda_status =
        launchFill(impl->stream.get(), 1, 0.0, impl->reduction_result.get());
    HIPDLP_CUDA_RECORD_PRIMITIVE();
    if (unavailableStatus(cuda_status))
      return cudaUnavailable("kernel compatibility probe", cuda_status);
    if (cuda_status != cudaSuccess)
      return cudaFailure("kernel compatibility probe", cuda_status);
    HIPDLP_CUDA_AUDIT(synchronization_calls);
    cuda_status = cudaStreamSynchronize(impl->stream.get());
    if (unavailableStatus(cuda_status))
      return cudaUnavailable("kernel compatibility probe", cuda_status);
    if (cuda_status != cudaSuccess)
      return cudaFailure("kernel compatibility probe", cuda_status);
    backend.reset(new CudaBackend(std::move(impl)));
  } catch (const std::bad_alloc&) {
    return api::Error(api::ErrorCode::kAllocationFailure,
                      "failed to allocate HiPDLP CUDA backend");
  }
  return api::Error();
}

}  // namespace cuda
}  // namespace backend
}  // namespace hipdlp