#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>

#include <cmath>
#include <limits>
#include <new>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

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
      free_calls(0),
      descriptor_creations(0),
      descriptor_destructions(0),
      host_to_device_transfers(0),
      device_to_host_transfers(0),
      device_to_device_transfers(0),
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
#define HIPDLP_CUDA_FAILURE(point) \
  (testing::failurePoint() == testing::FailurePoint::point)
#define HIPDLP_CUDA_CAPTURE_ACTIVE(active) \
  testing::setGraphCaptureActive(active)
#define HIPDLP_CUDA_RECORD_PRIMITIVE() testing::recordPrimitiveEnqueue()
#else
#define HIPDLP_CUDA_AUDIT(field) ((void)0)
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
    }
    value_ = nullptr;
    size_ = 0;
  }

  T* value_;
  std::size_t size_;
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
    status = cublasSetPointerMode(value_, CUBLAS_POINTER_MODE_HOST);
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

struct SparseStorage {
  HighsInt num_rows;
  HighsInt num_cols;
  HighsInt num_nonzeros;
  std::size_t spmv_buffer_size;
  DeviceBuffer<HighsInt> starts;
  DeviceBuffer<HighsInt> indices;
  DeviceBuffer<double> values;
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
  DeviceBuffer<unsigned char> spmv_buffer;

  CudaProblemResource()
      : lifetime(new int(0)),
        num_rows(0),
        num_cols(0),
        index_width(model::IndexWidth::kInt32) {}
};

struct CudaVectorResource : Resource {
  std::shared_ptr<int> lifetime;
  HighsInt size;
  DeviceBuffer<double> values;
  DenseVectorDescriptor descriptor;

  CudaVectorResource() : lifetime(new int(0)), size(0) {}
};

struct CudaWorkspaceResource : Resource {
  std::shared_ptr<int> lifetime;
  HighsInt num_rows;
  HighsInt num_cols;
  DeviceBuffer<double> row_scratch;
  DeviceBuffer<double> col_scratch;

  CudaWorkspaceResource() : lifetime(new int(0)), num_rows(0), num_cols(0) {}
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
  CublasHandle cublas;
  CusparseHandle cusparse;
  DeviceBuffer<int> error_flag;
  DeviceBuffer<double> reduction_result;
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
api::Error querySpmvBufferSize(ImplType& impl, SparseStorage& matrix) {
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
  const cusparseStatus_t status = cusparseSpMV_bufferSize(
      impl.cusparse.get(), CUSPARSE_OPERATION_NON_TRANSPOSE,
      impl.spmv_alpha.get(), matrix.descriptor.get(), input_descriptor.get(),
      impl.spmv_beta.get(), output_descriptor.get(), CUDA_R_64F,
      CUSPARSE_SPMV_ALG_DEFAULT, &matrix.spmv_buffer_size);
  return status == CUSPARSE_STATUS_SUCCESS
             ? api::Error()
             : cusparseFailure("SpMV buffer sizing", status);
}

template <typename ImplType>
api::Error prepareSpmvStorage(ImplType& impl, CudaProblemResource& resource) {
  api::Error error = querySpmvBufferSize(impl, resource.matrix);
  if (!error.ok()) return error;
  error = querySpmvBufferSize(impl, resource.transpose);
  if (!error.ok()) return error;
  const std::size_t buffer_size =
      resource.matrix.spmv_buffer_size > resource.transpose.spmv_buffer_size
          ? resource.matrix.spmv_buffer_size
          : resource.transpose.spmv_buffer_size;
  return resource.spmv_buffer.allocate(buffer_size,
                                       "allocating cuSPARSE SpMV workspace");
}

api::Error launchResult(const char* operation, const cudaError_t status) {
  HIPDLP_CUDA_RECORD_PRIMITIVE();
  return status == cudaSuccess ? api::Error() : cudaFailure(operation, status);
}

template <typename ImplType>
api::Error executeSpmv(ImplType& impl, const CudaProblemResource& problem,
                       const SparseStorage& matrix,
                       const CudaVectorResource& input,
                       CudaVectorResource& output, int* error_flag) {
  if (matrix.num_nonzeros == 0)
    return launchResult(
        "empty SpMV fill kernel launch",
        launchFill(impl.stream.get(), output.size, 0.0, output.values.get()));
  HIPDLP_CUDA_RECORD_PRIMITIVE();
  const cusparseStatus_t status = cusparseSpMV(
      impl.cusparse.get(), CUSPARSE_OPERATION_NON_TRANSPOSE,
      impl.spmv_alpha.get(), matrix.descriptor.get(), input.descriptor.get(),
      impl.spmv_beta.get(), output.descriptor.get(), CUDA_R_64F,
      CUSPARSE_SPMV_ALG_DEFAULT,
      const_cast<unsigned char*>(problem.spmv_buffer.get()));
  if (status != CUSPARSE_STATUS_SUCCESS)
    return cusparseFailure("SpMV execution", status);
  return launchResult("SpMV finite-check kernel launch",
                      launchFiniteCheck(impl.stream.get(), output.size,
                                        output.values.get(), error_flag));
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
  return dynamic_cast<CudaProblemResource*>(handle.resource_.get());
}

const CudaProblemResource* CudaBackend::problemResource(
    const ProblemHandle& handle) const {
  return dynamic_cast<const CudaProblemResource*>(handle.resource_.get());
}

CudaVectorResource* CudaBackend::vectorResource(VectorHandle& handle) const {
  return dynamic_cast<CudaVectorResource*>(handle.resource_.get());
}

const CudaVectorResource* CudaBackend::vectorResource(
    const VectorHandle& handle) const {
  return dynamic_cast<const CudaVectorResource*>(handle.resource_.get());
}

CudaWorkspaceResource* CudaBackend::workspaceResource(
    WorkspaceHandle& handle) const {
  return dynamic_cast<CudaWorkspaceResource*>(handle.resource_.get());
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
    std::unique_ptr<CudaProblemResource> resource(new CudaProblemResource());
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
        new CudaWorkspaceResource());
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

api::Error CudaBackend::createVector(const HighsInt size,
                                     VectorHandle& vector) {
  if (size < 0) return invalidArgument("vector size must be nonnegative");
  try {
    std::unique_ptr<CudaVectorResource> resource(new CudaVectorResource());
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
    return api::Error();
  }
  const cudaError_t status =
      cudaMemcpyAsync(cuda_destination->values.get(), cuda_source->values.get(),
                      cuda_source->values.bytes(), cudaMemcpyDeviceToDevice,
                      impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("vector copy", status);
  HIPDLP_CUDA_AUDIT(device_to_device_transfers);
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
  return executeSpmv(*impl_, *cuda_problem, cuda_problem->matrix, *cuda_input,
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
  return executeSpmv(*impl_, *cuda_problem, cuda_problem->transpose,
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

api::Error CudaBackend::fusedVectorBlend(
    WorkspaceHandle& workspace, const VectorHandle& primal_left,
    const VectorHandle& dual_left, const VectorHandle& primal_right,
    const VectorHandle& dual_right, const double left_weight,
    VectorHandle& primal_output, VectorHandle& dual_output) {
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
  impl_->recordDependency(cuda_workspace->lifetime);
  impl_->recordDependency(p_left->lifetime);
  impl_->recordDependency(d_left->lifetime);
  impl_->recordDependency(p_right->lifetime);
  impl_->recordDependency(d_right->lifetime);
  impl_->recordDependency(p_output->lifetime);
  impl_->recordDependency(d_output->lifetime);
  api::Error error = launchResult(
      "primal vector blend kernel launch",
      launchVectorBlend(impl_->stream.get(), p_left->size, p_left->values.get(),
                        p_right->values.get(), left_weight,
                        p_output->values.get(), impl_->error_flag.get()));
  if (!error.ok()) return error;
  return launchResult(
      "dual vector blend kernel launch",
      launchVectorBlend(impl_->stream.get(), d_left->size, d_left->values.get(),
                        d_right->values.get(), left_weight,
                        d_output->values.get(), impl_->error_flag.get()));
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
  api::Error error = launchResult(
      "dot reduction kernel launch",
      launchDot(impl_->stream.get(), cuda_left->size, cuda_left->values.get(),
                cuda_right->values.get(), impl_->reduction_result.get(),
                impl_->error_flag.get()));
  if (!error.ok()) return error;
  double reduced = 0.0;
  const cudaError_t status =
      cudaMemcpyAsync(&reduced, impl_->reduction_result.get(), sizeof(reduced),
                      cudaMemcpyDeviceToHost, impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("dot result download", status);
  HIPDLP_CUDA_AUDIT(device_to_host_transfers);
  error = impl_->synchronizeAndCheck();
  if (!error.ok()) return error;
  result = reduced;
  return api::Error();
}

api::Error CudaBackend::norm2(const VectorHandle& vector, double& result) {
  const CudaVectorResource* cuda_vector = vectorResource(vector);
  if (cuda_vector == nullptr)
    return invalidArgument("HiPDLP CUDA norm vector is invalid");
  api::Error error = launchResult(
      "norm reduction kernel launch",
      launchNorm2(impl_->stream.get(), cuda_vector->size,
                  cuda_vector->values.get(), impl_->reduction_result.get(),
                  impl_->error_flag.get()));
  if (!error.ok()) return error;
  double reduced = 0.0;
  const cudaError_t status =
      cudaMemcpyAsync(&reduced, impl_->reduction_result.get(), sizeof(reduced),
                      cudaMemcpyDeviceToHost, impl_->stream.get());
  if (status != cudaSuccess) return cudaFailure("norm result download", status);
  HIPDLP_CUDA_AUDIT(device_to_host_transfers);
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
  api::Error error =
      launchResult("difference norm reduction kernel launch",
                   launchDifferenceNorm2(
                       impl_->stream.get(), cuda_left->size,
                       cuda_left->values.get(), cuda_right->values.get(),
                       impl_->reduction_result.get(), impl_->error_flag.get()));
  if (!error.ok()) return error;
  double reduced = 0.0;
  const cudaError_t status =
      cudaMemcpyAsync(&reduced, impl_->reduction_result.get(), sizeof(reduced),
                      cudaMemcpyDeviceToHost, impl_->stream.get());
  if (status != cudaSuccess)
    return cudaFailure("difference norm result download", status);
  HIPDLP_CUDA_AUDIT(device_to_host_transfers);
  error = impl_->synchronizeAndCheck();
  if (!error.ok()) return error;
  result = reduced;
  return api::Error();
}

api::Error CudaBackend::reduceConvergenceMetrics(
    const VectorHandle& primal_residual, const VectorHandle& dual_residual,
    const double primal_objective, const double dual_objective,
    ConvergenceMetrics& metrics) {
  if (!std::isfinite(primal_objective) || !std::isfinite(dual_objective))
    return invalidArgument("convergence objectives must be finite");
  double primal_norm = 0.0;
  api::Error error = norm2(primal_residual, primal_norm);
  if (!error.ok()) return error;
  double dual_norm = 0.0;
  error = norm2(dual_residual, dual_norm);
  if (!error.ok()) return error;
  const double gap = std::abs(primal_objective - dual_objective);
  if (!std::isfinite(gap))
    return backendFailure("HiPDLP CUDA convergence gap is non-finite");
  metrics.primal_residual = primal_norm;
  metrics.dual_residual = dual_norm;
  metrics.objective_gap = gap;
  return api::Error();
}

api::Error CudaBackend::beginRecording() {
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

api::Error CudaBackend::downloadVector(const VectorHandle& vector,
                                       std::vector<double>& values) {
  const CudaVectorResource* cuda_vector = vectorResource(vector);
  if (cuda_vector == nullptr)
    return invalidArgument("HiPDLP CUDA download vector is invalid");
  try {
    std::vector<double> downloaded(static_cast<std::size_t>(cuda_vector->size));
    if (!downloaded.empty()) {
      const cudaError_t status =
          cudaMemcpyAsync(downloaded.data(), cuda_vector->values.get(),
                          downloaded.size() * sizeof(double),
                          cudaMemcpyDeviceToHost, impl_->stream.get());
      if (status != cudaSuccess) return cudaFailure("vector download", status);
      HIPDLP_CUDA_AUDIT(device_to_host_transfers);
    }
    api::Error error = impl_->synchronizeAndCheck();
    if (!error.ok()) return error;
    values.swap(downloaded);
  } catch (const std::bad_alloc&) {
    return allocationFailure("downloading a vector");
  }
  return api::Error();
}

api::Error CudaBackend::downloadSolution(const VectorHandle& primal,
                                         const VectorHandle& dual,
                                         std::vector<double>& primal_values,
                                         std::vector<double>& dual_values) {
  std::vector<double> downloaded_primal;
  std::vector<double> downloaded_dual;
  api::Error error = downloadVector(primal, downloaded_primal);
  if (!error.ok()) return error;
  error = downloadVector(dual, downloaded_dual);
  if (!error.ok()) return error;
  primal_values.swap(downloaded_primal);
  dual_values.swap(downloaded_dual);
  return api::Error();
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
    error = impl->reduction_result.allocate(1, "creating reduction storage");
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
    const double beta = 0.0;
    cuda_status = cudaMemcpyAsync(impl->spmv_beta.get(), &beta, sizeof(beta),
                                  cudaMemcpyHostToDevice, impl->stream.get());
    if (cuda_status != cudaSuccess)
      return cudaFailure("SpMV beta initialization", cuda_status);
    HIPDLP_CUDA_AUDIT(host_to_device_transfers);
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