#pragma once

// The solver is written against CUDA. With USE_HIP the CUDA names below are mapped to their
// HIP / ROCm equivalents (same arguments, same order), so the sources keep the CUDA spelling.
// The conditional graph API (cudaGraphConditionalHandle, ...) has no HIP equivalent and is only
// used without HIPDLP_HOST_LOOP, which USE_HIP defines (solver.cu).

#ifndef USE_HIP

#include <cub/block/block_reduce.cuh>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <nvtx3/nvToolsExt.h>

#else

#include <hip/hip_runtime.h>
#include <hipcub/block/block_reduce.hpp>
#include <hipsparse/hipsparse.h>
#include <rocprofiler-sdk-roctx/roctx.h>

// Runtime: types and constants
#define cudaDeviceProp hipDeviceProp_t
#define cudaError_t hipError_t
#define cudaEvent_t hipEvent_t
#define cudaFuncAttributes hipFuncAttributes
#define cudaGraphExec_t hipGraphExec_t
#define cudaGraphNode_t hipGraphNode_t
#define cudaGraph_t hipGraph_t
#define cudaStreamCaptureStatus hipStreamCaptureStatus
#define cudaStream_t hipStream_t
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaStreamCaptureModeGlobal hipStreamCaptureModeGlobal
#define cudaStreamCaptureStatusNone hipStreamCaptureStatusNone
#define cudaStreamNonBlocking hipStreamNonBlocking
#define cudaSuccess hipSuccess

// Runtime: errors (StatusOf in device.h)
#define cudaErrorInsufficientDriver hipErrorInsufficientDriver
#define cudaErrorInvalidDeviceFunction hipErrorInvalidDeviceFunction
#define cudaErrorMemoryAllocation hipErrorOutOfMemory
#define cudaErrorNoDevice hipErrorNoDevice
#define cudaErrorNoKernelImageForDevice hipErrorNoBinaryForGpu
#define cudaErrorNotSupported hipErrorNotSupported
// No HIP counterpart. The closest "this device cannot run it" code, distinct from the above.
#define cudaErrorCompatNotSupportedOnDevice hipErrorInvalidDevice

// Runtime: functions
#define cudaEventCreate hipEventCreate
#define cudaEventDestroy hipEventDestroy
#define cudaEventElapsedTime hipEventElapsedTime
#define cudaEventRecord hipEventRecord
#define cudaFree hipFree
#define cudaGetDevice hipGetDevice
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError
#define cudaGraphDestroy hipGraphDestroy
#define cudaGraphExecDestroy hipGraphExecDestroy
// hipGraphInstantiate takes an error node and a log buffer, the CUDA 12 form takes flags.
#define cudaGraphInstantiate hipGraphInstantiateWithFlags
#define cudaGraphLaunch hipGraphLaunch
#define cudaMalloc hipMalloc
#define cudaMemGetInfo hipMemGetInfo
#define cudaMemcpy hipMemcpy
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemsetAsync hipMemsetAsync
#define cudaStreamBeginCapture hipStreamBeginCapture
#define cudaStreamCreateWithFlags hipStreamCreateWithFlags
#define cudaStreamDestroy hipStreamDestroy
#define cudaStreamEndCapture hipStreamEndCapture
#define cudaStreamIsCapturing hipStreamIsCapturing
#define cudaStreamSynchronize hipStreamSynchronize
// hipFuncGetAttributes takes the kernel as const void*, cudaFuncGetAttributes as a function.
template <typename Kernel>
inline hipError_t HipFuncGetAttributes(hipFuncAttributes* attributes, Kernel kernel) {
  return hipFuncGetAttributes(attributes, reinterpret_cast<const void*>(kernel));
}
#define cudaFuncGetAttributes HipFuncGetAttributes

// cuSPARSE
#define CUDA_R_64F HIP_R_64F
#define CUSPARSE_ACTION_NUMERIC HIPSPARSE_ACTION_NUMERIC
#define CUSPARSE_CSR2CSC_ALG1 HIPSPARSE_CSR2CSC_ALG1
#define CUSPARSE_INDEX_32I HIPSPARSE_INDEX_32I
#define CUSPARSE_INDEX_BASE_ZERO HIPSPARSE_INDEX_BASE_ZERO
#define CUSPARSE_OPERATION_NON_TRANSPOSE HIPSPARSE_OPERATION_NON_TRANSPOSE
#define CUSPARSE_POINTER_MODE_DEVICE HIPSPARSE_POINTER_MODE_DEVICE
#define CUSPARSE_SPMV_CSR_ALG2 HIPSPARSE_SPMV_CSR_ALG2
#define CUSPARSE_STATUS_ALLOC_FAILED HIPSPARSE_STATUS_ALLOC_FAILED
#define CUSPARSE_STATUS_SUCCESS HIPSPARSE_STATUS_SUCCESS
#define cusparseDnVecDescr_t hipsparseDnVecDescr_t
#define cusparseHandle_t hipsparseHandle_t
#define cusparseSpMatDescr_t hipsparseSpMatDescr_t
#define cusparseStatus_t hipsparseStatus_t
#define cusparseCreate hipsparseCreate
#define cusparseCreateCsr hipsparseCreateCsr
#define cusparseCreateDnVec hipsparseCreateDnVec
#define cusparseCsr2cscEx2 hipsparseCsr2cscEx2
#define cusparseCsr2cscEx2_bufferSize hipsparseCsr2cscEx2_bufferSize
#define cusparseDestroy hipsparseDestroy
#define cusparseDestroyDnVec hipsparseDestroyDnVec
#define cusparseDestroySpMat hipsparseDestroySpMat
#define cusparseGetErrorString hipsparseGetErrorString
#define cusparseSetPointerMode hipsparseSetPointerMode
#define cusparseSetStream hipsparseSetStream
#define cusparseSpMV hipsparseSpMV
#define cusparseSpMV_bufferSize hipsparseSpMV_bufferSize
#define cusparseSpMV_preprocess hipsparseSpMV_preprocess

// CUB and NVTX
#define cub hipcub
#define nvtxRangePop roctxRangePop
#define nvtxRangePushA roctxRangePushA

#endif  // USE_HIP
