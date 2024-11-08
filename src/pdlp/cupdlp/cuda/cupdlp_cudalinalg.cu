#include "cupdlp_cudalinalg.cuh"

extern "C" cupdlp_int cuda_alloc_MVbuffer(
    cusparseHandle_t handle, cusparseSpMatDescr_t cuda_csc,
    cusparseDnVecDescr_t vecX, cusparseDnVecDescr_t vecAx,
    cusparseSpMatDescr_t cuda_csr, cusparseDnVecDescr_t vecY,
    cusparseDnVecDescr_t vecATy, void **dBuffer) {
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  size_t AxBufferSize = 0;
  size_t ATyBufferSize = 0;
  cupdlp_float alpha = 1.0;
  cupdlp_float beta = 0.0;
  // cusparseSpSVAlg_t alg = CUSPARSE_SPSV_ALG_DEFAULT;
  cusparseSpMVAlg_t alg = CUSPARSE_SPMV_CSR_ALG2; //deterministic

  // get the buffer size needed by csr Ax
  CHECK_CUSPARSE(cusparseSpMV_bufferSize(
      handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, cuda_csr, vecX, &beta,
      vecAx, computeType, alg, &AxBufferSize))

  // get the buffer size needed by csc ATy
  CHECK_CUSPARSE(cusparseSpMV_bufferSize(
      handle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, cuda_csc, vecY, &beta,
      vecATy, computeType, alg, &ATyBufferSize))

  size_t bufferSize =
      (AxBufferSize > ATyBufferSize) ? AxBufferSize : ATyBufferSize;

  // allocate an external buffer if needed
  CHECK_CUDA(cudaMalloc(dBuffer, bufferSize))

  return EXIT_SUCCESS;
}

extern "C" cupdlp_int cuda_csc_Ax(cusparseHandle_t handle,
                                  cusparseSpMatDescr_t cuda_csc,
                                  cusparseDnVecDescr_t vecX,
                                  cusparseDnVecDescr_t vecAx, void *dBuffer,
                                  const cupdlp_float alpha,
                                  const cupdlp_float beta) {
  // hAx = alpha * Acsc * hX + beta * hAx

  cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csc, vecX, &beta, vecAx,
                              // computeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              computeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

extern "C" cupdlp_int cuda_csr_Ax(cusparseHandle_t handle,
                                  cusparseSpMatDescr_t cuda_csr,
                                  cusparseDnVecDescr_t vecX,
                                  cusparseDnVecDescr_t vecAx, void *dBuffer,
                                  const cupdlp_float alpha,
                                  const cupdlp_float beta) {
  // hAx = alpha * Acsc * hX + beta * hAx

  cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csr, vecX, &beta, vecAx,
                              // computeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              computeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

extern "C" cupdlp_int cuda_csc_ATy(cusparseHandle_t handle,
                                   cusparseSpMatDescr_t cuda_csc,
                                   cusparseDnVecDescr_t vecY,
                                   cusparseDnVecDescr_t vecATy, void *dBuffer,
                                   const cupdlp_float alpha,
                                   const cupdlp_float beta) {
  // hATy = alpha * Acsr^T * hY + beta * hATy
  cusparseOperation_t op = CUSPARSE_OPERATION_TRANSPOSE;
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csc, vecY, &beta, vecATy,
                              // computeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              computeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

extern "C" cupdlp_int cuda_csr_ATy(cusparseHandle_t handle,
                                   cusparseSpMatDescr_t cuda_csr,
                                   cusparseDnVecDescr_t vecY,
                                   cusparseDnVecDescr_t vecATy, void *dBuffer,
                                   const cupdlp_float alpha,
                                   const cupdlp_float beta) {
  // hATy = alpha * Acsr^T * hY + beta * hATy
  cusparseOperation_t op = CUSPARSE_OPERATION_TRANSPOSE;
  cudaDataType computeType = CUDA_R_32F;
#ifndef SFLOAT
  computeType = CUDA_R_64F;
#endif

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csr, vecY, &beta, vecATy,
                              // computeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              computeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

extern "C" void cupdlp_projSameub_cuda(cupdlp_float *x, const cupdlp_float ub,
                                       const cupdlp_int len) {
  element_wise_projSameub_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      x, ub, len);
}

extern "C" void cupdlp_projSamelb_cuda(cupdlp_float *x, const cupdlp_float lb,
                                       const cupdlp_int len) {
  element_wise_projSamelb_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      x, lb, len);
}

extern "C" void cupdlp_projub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                                   const cupdlp_int len) {
  element_wise_projub_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, ub,
                                                                        len);
}

extern "C" void cupdlp_projlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                                   const cupdlp_int len) {
  element_wise_projlb_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, lb,
                                                                        len);
}

extern "C" void cupdlp_ediv_cuda(cupdlp_float *x, const cupdlp_float *y,
                                 const cupdlp_int len) {
  element_wise_div_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, y, len);
}

extern "C" void cupdlp_edot_cuda(cupdlp_float *x, const cupdlp_float *y,
                                 const cupdlp_int len) {
  element_wise_dot_kernel<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, y, len);
}

extern "C" void cupdlp_haslb_cuda(cupdlp_float *haslb, const cupdlp_float *lb,
                                  const cupdlp_float bound,
                                  const cupdlp_int len) {
  element_wise_initHaslb_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      haslb, lb, bound, len);
}

extern "C" void cupdlp_hasub_cuda(cupdlp_float *hasub, const cupdlp_float *ub,
                                  const cupdlp_float bound,
                                  const cupdlp_int len) {
  element_wise_initHasub_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      hasub, ub, bound, len);
}

extern "C" void cupdlp_filterlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                                     const cupdlp_float bound,
                                     const cupdlp_int len) {
  element_wise_filterlb_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      x, lb, bound, len);
}

extern "C" void cupdlp_filterub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                                     const cupdlp_float bound,
                                     const cupdlp_int len) {
  element_wise_filterub_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
      x, ub, bound, len);
}

extern "C" void cupdlp_initvec_cuda(cupdlp_float *x, const cupdlp_float val,
                                    const cupdlp_int len) {
  init_cuda_vec_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(x, val, len);
}

extern "C" void cupdlp_pgrad_cuda(cupdlp_float *xUpdate,
                                        const cupdlp_float *x,
                                        const cupdlp_float *cost,
                                        const cupdlp_float *ATy,
                                        const cupdlp_float dPrimalStep,
                                        const cupdlp_int len) {
    primal_grad_step_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
        xUpdate, x, cost, ATy, dPrimalStep, len);
}

extern "C" void cupdlp_dgrad_cuda(cupdlp_float *yUpdate, const cupdlp_float *y, const cupdlp_float *b,
    const cupdlp_float *Ax, const cupdlp_float *AxUpdate,
    const cupdlp_float dDualStep, const cupdlp_int len) {
      dual_grad_step_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(
          yUpdate, y, b, Ax, AxUpdate, dDualStep, len);
}

extern "C" void cupdlp_sub_cuda(cupdlp_float *z, const cupdlp_float *x,
                                  const cupdlp_float *y, const cupdlp_int len)
{
   naive_sub_kernal<<<cuda_gridsize(len), CUPDLP_BLOCK_SIZE>>>(z, x, y, len);
}