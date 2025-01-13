#ifndef CUPDLP_CUDA_KERNALS_H
#define CUPDLP_CUDA_KERNALS_H

#include "cuda_runtime.h"
#define CUPDLP_BLOCK_SIZE 512

#ifndef SFLOAT
#ifdef DLONG
typedef long long cupdlp_int;
#else
typedef int cupdlp_int;
#endif
typedef double cupdlp_float;
#define CudaComputeType CUDA_R_64F
#else
#define CudaComputeType CUDA_R_32F
#endif

#define CHECK_CUDA(func)                                               \
  {                                                                    \
    cudaError_t status = (func);                                       \
    if (status != cudaSuccess) {                                       \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n", \
             __LINE__, __FILE__, cudaGetErrorString(status), status);  \
      return EXIT_FAILURE;                                             \
    }                                                                  \
  }

#define CHECK_CUSPARSE(func)                                               \
  {                                                                        \
    cusparseStatus_t status = (func);                                      \
    if (status != CUSPARSE_STATUS_SUCCESS) {                               \
      printf("CUSPARSE API failed at line %d of %s with error: %s (%d)\n", \
             __LINE__, __FILE__, cusparseGetErrorString(status), status);  \
      return EXIT_FAILURE;                                                 \
    }                                                                      \
  }

#define CHECK_CUBLAS(func)                                               \
  {                                                                      \
    cublasStatus_t status = (func);                                      \
    if (status != CUBLAS_STATUS_SUCCESS) {                               \
      printf("CUBLAS API failed at line %d of %s with error: %s (%d)\n", \
             __LINE__, __FILE__, cublasGetStatusString(status), status); \
      return EXIT_FAILURE;                                               \
    }                                                                    \
  }

#define CUPDLP_FREE_VEC(x) \
  {                        \
    cudaFree(x);           \
    x = cupdlp_NULL;       \
  }

#define CUPDLP_COPY_VEC(dst, src, type, size) \
  cudaMemcpy(dst, src, sizeof(type) * (size), cudaMemcpyDefault)

// #define CUPDLP_INIT_VEC(var, size)                                             \
//   {                                                                            \
//     cusparseStatus_t status =                                                  \
//         cudaMalloc((void **)&var, (size) * sizeof(typeof(*var)));              \
//     if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
//       printf("CUSPARSE API failed at line %d with error: %s (%d)\n", __LINE__, \
//              cusparseGetErrorString(status), status);                          \
//       goto exit_cleanup;                                                       \
//     }                                                                          \
//   }

#define cupdlp_init_vec_int(var, size)                                   \
  {                                                                      \
    cudaError_t status = cudaMalloc((void**)&var, (size) * sizeof(int)); \
    if (status != cudaSuccess) {                                         \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",   \
             __LINE__, __FILE__, cudaGetErrorString(status), status);    \
      goto exit_cleanup;                                                 \
    }                                                                    \
  }

#define cupdlp_init_vec_double(var, size)                                   \
  {                                                                         \
    cudaError_t status = cudaMalloc((void**)&var, (size) * sizeof(double)); \
    if (status != cudaSuccess) {                                            \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",      \
             __LINE__, __FILE__, cudaGetErrorString(status), status);       \
      goto exit_cleanup;                                                    \
    }                                                                       \
  }

// #define CUPDLP_INIT_ZERO_VEC(var, size)                                        \
//   {                                                                            \
//     cusparseStatus_t status =                                                  \
//         cudaMalloc((void**)&var, (size) * sizeof(typeof(*var)));               \
//     if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
//       printf("CUSPARSE API failed at line %d with error: %s (%d)\n", __LINE__, \
//              cusparseGetErrorString(status), status);                          \
//       goto exit_cleanup;                                                       \
//     }                                                                          \
//     status = cudaMemset(var, 0, (size) * sizeof(typeof(*var)));                \
//     if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
//       printf("CUSPARSE API failed at line %d with error: %s (%d)\n", __LINE__, \
//              cusparseGetErrorString(status), status);                          \
//       goto exit_cleanup;                                                       \
//     }                                                                          \
//   }

#define cupdlp_init_zero_vec_int(var, size)                              \
  {                                                                      \
    cudaError_t status = cudaMalloc((void**)&var, (size) * sizeof(int)); \
    if (status != cudaSuccess) {                                         \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",   \
             __LINE__, __FILE__, cudaGetErrorString(status), status);    \
      goto exit_cleanup;                                                 \
    }                                                                    \
    status = cudaMemset(var, 0, (size) * sizeof(int));                   \
    if (status != cudaSuccess) {                                         \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",   \
             __LINE__, __FILE__, cudaGetErrorString(status), status);    \
      goto exit_cleanup;                                                 \
    }                                                                    \
  }
#define cupdlp_init_zero_vec_double(var, size)                              \
  {                                                                         \
    cudaError_t status = cudaMalloc((void**)&var, (size) * sizeof(double)); \
    if (status != cudaSuccess) {                                            \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",      \
             __LINE__, __FILE__, cudaGetErrorString(status), status);       \
      goto exit_cleanup;                                                    \
    }                                                                       \
    status = cudaMemset(var, 0, (size) * sizeof(double));                   \
    if (status != cudaSuccess) {                                            \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",      \
             __LINE__, __FILE__, cudaGetErrorString(status), status);       \
      goto exit_cleanup;                                                    \
    }                                                                       \
  }

#define CUPDLP_ZERO_VEC(var, type, size) \
  cudaMemset(var, 0, sizeof(type) * (size))

dim3 cuda_gridsize(cupdlp_int n);

#define get_gpu_vec(vec_ptr, index, result, success)       \
  {                                                                \
    cublasStatus_t stat =                                          \
      cublasGetVector(N, sizeof(*vec_ptr), vec_ptr, 1, result, 1); \
    if (stat != CUBLAS_STATUS_SUCCESS) {                           \
      printf("data upload failed");                                \
      *success = 0;                                                 \
    } else {                                                       \
      *success = 1;                                                 \
    }                                                              \
  }

#define get_gpu_vec_element(vec_ptr, index, result, success)       \
  {                                                                \
    cublasStatus_t stat =                                          \
        cublasGetVector(1, sizeof(double), vec_ptr, 1, result, 1); \
    if (stat != CUBLAS_STATUS_SUCCESS) {                           \
      printf("data upload failed");                                \
      *success = 0;                                                 \
    } else {                                                       \
      *success = 1;                                                 \
    }                                                              \
  }

__global__ void element_wise_dot_kernel(cupdlp_float* x, const cupdlp_float* y,
                                        const cupdlp_int len);

__global__ void element_wise_div_kernel(cupdlp_float* x, const cupdlp_float* y,
                                        const cupdlp_int len);

__global__ void element_wise_projlb_kernel(cupdlp_float* x,
                                           const cupdlp_float* lb,
                                           const cupdlp_int len);

__global__ void element_wise_projub_kernel(cupdlp_float* x,
                                           const cupdlp_float* ub,
                                           const cupdlp_int len);

__global__ void element_wise_projSamelb_kernel(cupdlp_float* x,
                                               const cupdlp_float lb,
                                               const cupdlp_int len);

__global__ void element_wise_projSameub_kernel(cupdlp_float* x,
                                               const cupdlp_float ub,
                                               const cupdlp_int len);

__global__ void element_wise_initHaslb_kernal(cupdlp_float* haslb,
                                              const cupdlp_float* lb,
                                              const cupdlp_float bound,
                                              const cupdlp_int len);

__global__ void element_wise_initHasub_kernal(cupdlp_float* hasub,
                                              const cupdlp_float* ub,
                                              const cupdlp_float bound,
                                              const cupdlp_int len);

__global__ void element_wise_filterlb_kernal(cupdlp_float* x,
                                             const cupdlp_float* lb,
                                             const cupdlp_float bound,
                                             const cupdlp_int len);

__global__ void element_wise_filterub_kernal(cupdlp_float* x,
                                             const cupdlp_float* ub,
                                             const cupdlp_float bound,
                                             const cupdlp_int len);

__global__ void init_cuda_vec_kernal(cupdlp_float* x, const cupdlp_float val,
                                     const cupdlp_int len);

__global__ void primal_grad_step_kernal(cupdlp_float* xUpdate,
                                        const cupdlp_float* x,
                                        const cupdlp_float* cost,
                                        const cupdlp_float* ATy,
                                        const cupdlp_float dPrimalStep,
                                        const cupdlp_int len);

__global__ void dual_grad_step_kernal(
    cupdlp_float* yUpdate, const cupdlp_float* y, const cupdlp_float* b,
    const cupdlp_float* Ax, const cupdlp_float* AxUpdate,
    const cupdlp_float dDualStep, const cupdlp_int len);

__global__ void naive_sub_kernal(cupdlp_float* z, const cupdlp_float* x,
                                 const cupdlp_float* y, const cupdlp_int len);
#endif