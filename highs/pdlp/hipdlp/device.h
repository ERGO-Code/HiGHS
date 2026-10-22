#pragma once

#include <cuda_runtime.h>
#include <cusparse.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

#include "solver.h"

namespace hipdlp {
    
struct DeviceFailure : std::runtime_error {
    Status status;
    DeviceFailure(Status code, const std::string& text) : std::runtime_error(text), status(code) {}
};

inline void Cuda(cudaError_t error) {
    if (error != cudaSuccess) {
        cudaGetLastError(); // Clear the error state
        throw DeviceFailure(
            error == cudaErrorMemoryAllocation ? Status::kOutOfMemory : Status::kDeviceError,
            cudaGetErrorString(error));
    }
}

inline void Sparse(cusparseStatus_t error) {
  if (error != CUSPARSE_STATUS_SUCCESS) {
    cudaGetLastError();
    throw DeviceFailure(
        error == CUSPARSE_STATUS_ALLOC_FAILED ? Status::kOutOfMemory : Status::kDeviceError,
        cusparseGetErrorString(error));
  }
}

template <typename T>
class DeviceBuffer {
 public:
  DeviceBuffer() = default;
  explicit DeviceBuffer(size_t count) : count_(count) {
    void* p = nullptr;
    Cuda(cudaMalloc(&p, std::max(size_t(1), count) * sizeof(T)));
    ptr_ = static_cast<T*>(p);
  }
  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;
  DeviceBuffer(DeviceBuffer&& other) noexcept
      : ptr_(std::exchange(other.ptr_, nullptr)), count_(std::exchange(other.count_, 0)) {}
  DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
    if (this != &other) {
      Reset();
      ptr_ = std::exchange(other.ptr_, nullptr);
      count_ = std::exchange(other.count_, 0);
    }
    return *this;
  }
  ~DeviceBuffer() { Reset(); }
  void Reset() {
    if (ptr_) cudaFree(ptr_);
    ptr_ = nullptr;
    count_ = 0;
  }
  T* Get() const { return ptr_; }
  size_t Size() const { return count_; }
  size_t Bytes() const { return ptr_ ? std::max(size_t(1), count_) * sizeof(T) : 0; }
  explicit operator bool() const { return ptr_ != nullptr; }

 private:
  T* ptr_ = nullptr;
  size_t count_ = 0;
};

struct DeviceCsr {
    int rows = 0, columns = 0, nnz = 0;
    DeviceBuffer<int> row_start, column_index;
    DeviceBuffer<double> value;
    size_t Bytes() const { 
        return row_start.Bytes() + column_index.Bytes() + value.Bytes();
    }
    explicit operator bool() const { return bool(row_start);}
};

struct DeviceTranspose {
  DeviceBuffer<int> start, row, source_slot;
  DeviceBuffer<double> value;
  size_t Bytes() const { return start.Bytes() + row.Bytes() + source_slot.Bytes() + value.Bytes(); }
};

// The unscaled original matrix of one model kept on the device between solves (the daily
// reuse): A, A^T and, for numeric patches, the CSR slot -> transposed position map. Identity,
// revision and lifetime belong to the owner (ResidentModel); the solver only borrows it, and may
// release it (Reset the members) when the device memory is needed for a reduced solve.
struct ResidentMatrix {
  DeviceCsr a;
  DeviceTranspose at;
  DeviceBuffer<int> position_of_slot;
  size_t Bytes() const { return a.Bytes() + at.Bytes() + position_of_slot.Bytes(); }
  explicit operator bool() const { return bool(a); }
  void Reset() {
    a = DeviceCsr();
    at = DeviceTranspose();
    position_of_slot.Reset();
  }
};

// Process-wide device memory in use above the level at construction (the CUDA context excluded).
// Sampling after each allocation step gives the true peak of everything alive together, including
// cuSPARSE/CUB workspaces and buffers owned by other stages; per-owner byte counts do not.
class MemoryMeter {
public:
    MemoryMeter() : baseline_(Used()) {}
    void Sample() {
        size_t used = Used();
        if (used > basedline_) peak_ = std::max(peak_, used - baseline_);
    }
    size_t Peak() const { return peak_; }
private:
    static size_t Used() {
        size_t free_bytes = 0, total_bytes = 0;
        Cuda(cudaMemGetInfo(&free_bytes, &total_bytes));
        return total_bytes - free_bytes;
    }
    size_t baseline_ = 0, peak_ = 0;
};

// cuSPARSE handle owner for the one-off conversions outside a solve session
struct SparseHandle {
  cusparseHandle_t handle = nullptr;
  SparseHandle() { Sparse(cusparseCreate(&handle)); }
  SparseHandle(const SparseHandle&) = delete;
  SparseHandle& operator=(const SparseHandle&) = delete;
  ~SparseHandle() { cusparseDestroy(handle); }
};

DeviceCsr UploadCsr(const HostLp& lp, cudaStream_t stream, Profile* profile = nullptr);

void BuildTranspose(const DeviceCsr& a, cusparseHandle_t sparse, cudaStream_t stream,
                    bool want_slot_map, DeviceTranspose* out, Profile* profile = nullptr,
                    MemoryMeter* meter = nullptr);
// Test helper: copies a device matrix back (synchronizing).
HostLp DownloadCsr(const DeviceCsr& a, cudaStream_t stream);

}  // namespace hipdlp