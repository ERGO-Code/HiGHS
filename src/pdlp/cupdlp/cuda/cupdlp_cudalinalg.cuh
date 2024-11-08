#ifndef CUPDLP_CUDA_LINALG_H
#define CUPDLP_CUDA_LINALG_H

#include <cublas_v2.h>         // cublas
#include <cuda_runtime_api.h>  // cudaMalloc, cudaMemcpy, etc.
#include <cusparse.h>          // cusparseSpMV

#include "cupdlp_cuda_kernels.cuh"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>   // printf
#include <stdlib.h>  // EXIT_FAILURE

// #include "../cupdlp_defs.h"
// #include "../glbopts.h"
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

extern "C" cupdlp_int cuda_alloc_MVbuffer(
    cusparseHandle_t handle, cusparseSpMatDescr_t cuda_csc,
    cusparseDnVecDescr_t vecX, cusparseDnVecDescr_t vecAx,
    cusparseSpMatDescr_t cuda_csr, cusparseDnVecDescr_t vecY,
    cusparseDnVecDescr_t vecATy, void **dBuffer);

extern "C" cupdlp_int cuda_csc_Ax(cusparseHandle_t handle,
                                  cusparseSpMatDescr_t cuda_csc,
                                  cusparseDnVecDescr_t vecX,
                                  cusparseDnVecDescr_t vecAx, void *dBuffer,
                                  const cupdlp_float alpha,
                                  const cupdlp_float beta);
extern "C" cupdlp_int cuda_csr_Ax(cusparseHandle_t handle,
                                  cusparseSpMatDescr_t cuda_csr,
                                  cusparseDnVecDescr_t vecX,
                                  cusparseDnVecDescr_t vecAx, void *dBuffer,
                                  const cupdlp_float alpha,
                                  const cupdlp_float beta);
extern "C" cupdlp_int cuda_csc_ATy(cusparseHandle_t handle,
                                   cusparseSpMatDescr_t cuda_csc,
                                   cusparseDnVecDescr_t vecY,
                                   cusparseDnVecDescr_t vecATy, void *dBuffer,
                                   const cupdlp_float alpha,
                                   const cupdlp_float beta);
extern "C" cupdlp_int cuda_csr_ATy(cusparseHandle_t handle,
                                   cusparseSpMatDescr_t cuda_csr,
                                   cusparseDnVecDescr_t vecY,
                                   cusparseDnVecDescr_t vecATy, void *dBuffer,
                                   const cupdlp_float alpha,
                                   const cupdlp_float beta);

extern "C" void cupdlp_projSameub_cuda(cupdlp_float *x, const cupdlp_float ub,
                                       const cupdlp_int len);
extern "C" void cupdlp_projSamelb_cuda(cupdlp_float *x, const cupdlp_float lb,
                                       const cupdlp_int len);
extern "C" void cupdlp_projub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                                   const cupdlp_int len);
extern "C" void cupdlp_projlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                                   const cupdlp_int len);
extern "C" void cupdlp_ediv_cuda(cupdlp_float *x, const cupdlp_float *y,
                                 const cupdlp_int len);
extern "C" void cupdlp_edot_cuda(cupdlp_float *x, const cupdlp_float *y,
                                 const cupdlp_int len);
extern "C" void cupdlp_haslb_cuda(cupdlp_float *haslb, const cupdlp_float *lb,
                                  const cupdlp_float bound,
                                  const cupdlp_int len);
extern "C" void cupdlp_hasub_cuda(cupdlp_float *hasub, const cupdlp_float *ub,
                                  const cupdlp_float bound,
                                  const cupdlp_int len);
extern "C" void cupdlp_filterlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                                     const cupdlp_float bound,
                                     const cupdlp_int len);
extern "C" void cupdlp_filterub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                                     const cupdlp_float bound,
                                     const cupdlp_int len);
extern "C" void cupdlp_initvec_cuda(cupdlp_float *x, const cupdlp_float val,
                                    const cupdlp_int len);

extern "C" void cupdlp_pgrad_cuda(cupdlp_float *xUpdate, const cupdlp_float *x,
                                  const cupdlp_float *cost,
                                  const cupdlp_float *ATy,
                                  const cupdlp_float dPrimalStep,
                                  const cupdlp_int len);

extern "C" void cupdlp_dgrad_cuda(cupdlp_float *yUpdate, const cupdlp_float *y,
                                  const cupdlp_float *b, const cupdlp_float *Ax,
                                  const cupdlp_float *AxUpdate,
                                  const cupdlp_float dDualStep,
                                  const cupdlp_int len);

extern "C" void cupdlp_sub_cuda(cupdlp_float *z, const cupdlp_float *x,
                                const cupdlp_float *y, const cupdlp_int len);
#endif
#endif