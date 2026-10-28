#include <algorithm>

#include "device.h"

namespace hipdlp {

DeviceCsr UploadCsr(const HostLp& lp, cudaStream_t stream, Profile* profile) {
  DeviceCsr a;
  a.rows = lp.rows;
  a.columns = lp.columns;
  a.nnz = int(lp.value.size());
  a.row_start = DeviceBuffer<int>(size_t(lp.rows) + 1);
  a.column_index = DeviceBuffer<int>(lp.value.size());
  a.value = DeviceBuffer<double>(lp.value.size());
  Cuda(cudaMemcpyAsync(a.row_start.Get(), lp.row_start.data(), (size_t(lp.rows) + 1) * sizeof(int),
                       cudaMemcpyHostToDevice, stream));
  if (a.nnz) {
    Cuda(cudaMemcpyAsync(a.column_index.Get(), lp.column_index.data(), a.nnz * sizeof(int),
                         cudaMemcpyHostToDevice, stream));
    Cuda(cudaMemcpyAsync(a.value.Get(), lp.value.data(), a.nnz * sizeof(double),
                         cudaMemcpyHostToDevice, stream));
  }
  if (profile) profile->h2d_bytes += a.Bytes();
  return a;
}

static __global__ void FillSlots(double* values, int size) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) values[i] = i;
}
static __global__ void SlotsToIndex(const double* values, int* index, int size) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) index[i] = int(values[i]);
}

void BuildTranspose(const DeviceCsr& a, cusparseHandle_t sparse, cudaStream_t stream,
                    bool want_slot_map, DeviceTranspose* out, Profile* profile,
                    MemoryMeter* meter) {
  int m = a.rows, n = a.columns, nnz = a.nnz;
  out->start = DeviceBuffer<int>(size_t(n) + 1);
  out->row = DeviceBuffer<int>(nnz);
  out->value = DeviceBuffer<double>(nnz);
  if (want_slot_map) out->source_slot = DeviceBuffer<int>(nnz);
  if (nnz == 0 || m == 0 || n == 0) {  // an empty matrix has nothing to convert
    Cuda(cudaMemsetAsync(out->start.Get(), 0, (size_t(n) + 1) * sizeof(int), stream));
    return;
  }
  Sparse(cusparseSetStream(sparse, stream));
  auto convert = [&](const double* from, double* to, DeviceBuffer<char>& workspace) {
    Sparse(cusparseCsr2cscEx2(sparse, m, n, nnz, from, a.row_start.Get(), a.column_index.Get(), to,
                              out->start.Get(), out->row.Get(), CUDA_R_64F, CUSPARSE_ACTION_NUMERIC,
                              CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1, workspace.Get()));
  };
  size_t bytes = 0;
  Sparse(cusparseCsr2cscEx2_bufferSize(sparse, m, n, nnz, a.value.Get(), a.row_start.Get(),
                                       a.column_index.Get(), out->value.Get(), out->start.Get(),
                                       out->row.Get(), CUDA_R_64F, CUSPARSE_ACTION_NUMERIC,
                                       CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1, &bytes));
  if (profile)
    profile->transpose_workspace_bytes = std::max(profile->transpose_workspace_bytes, bytes);
  DeviceBuffer<char> workspace(bytes);
  if (meter) meter->Sample();  // the transient peak: arrays, outputs and the workspace together
  convert(a.value.Get(), out->value.Get(), workspace);
  if (want_slot_map) {
    // The conversion of the slot numbers (exact in FP64) is the position permutation.
    DeviceBuffer<double> slots(nnz), moved(nnz);
    if (meter) meter->Sample();
    int blocks = int((int64_t(nnz) + 255) / 256);
    FillSlots<<<blocks, 256, 0, stream>>>(slots.Get(), nnz);
    convert(slots.Get(), moved.Get(), workspace);
    SlotsToIndex<<<blocks, 256, 0, stream>>>(moved.Get(), out->source_slot.Get(), nnz);
    Cuda(cudaGetLastError());
    Cuda(cudaStreamSynchronize(stream));  // slots/moved are released at scope exit
  }
}

HostLp DownloadCsr(const DeviceCsr& a, cudaStream_t stream) {
  HostLp lp;
  lp.rows = a.rows;
  lp.columns = a.columns;
  lp.row_start.resize(size_t(a.rows) + 1);
  lp.column_index.resize(a.nnz);
  lp.value.resize(a.nnz);
  Cuda(cudaMemcpyAsync(lp.row_start.data(), a.row_start.Get(), lp.row_start.size() * sizeof(int),
                       cudaMemcpyDeviceToHost, stream));
  if (a.nnz) {
    Cuda(cudaMemcpyAsync(lp.column_index.data(), a.column_index.Get(), a.nnz * sizeof(int),
                         cudaMemcpyDeviceToHost, stream));
    Cuda(cudaMemcpyAsync(lp.value.data(), a.value.Get(), a.nnz * sizeof(double),
                         cudaMemcpyDeviceToHost, stream));
  }
  Cuda(cudaStreamSynchronize(stream));
  return lp;
}

}  // namespace hipdlp
