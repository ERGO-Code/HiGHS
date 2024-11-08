#include "cupdlp_cuda_kernels.cuh"

dim3 cuda_gridsize(cupdlp_int n) {
  cupdlp_int k = (n - 1) / CUPDLP_BLOCK_SIZE + 1;
  cupdlp_int x = k;
  cupdlp_int y = 1;
  if (x > 65535) {
    x = ceil(sqrt(k));
    y = (n - 1) / (x * CUPDLP_BLOCK_SIZE) + 1;
  }
  dim3 d = {x, y, 1};
  return d;
}

__global__ void element_wise_dot_kernel(cupdlp_float *x, const cupdlp_float *y,
                                        const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] *= y[i];
}

__global__ void element_wise_div_kernel(cupdlp_float *x, const cupdlp_float *y,
                                        const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] /= y[i];
}

__global__ void element_wise_projlb_kernel(cupdlp_float *x,
                                           const cupdlp_float *lb,
                                           const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = x[i] < lb[i] ? lb[i] : x[i];
}

__global__ void element_wise_projub_kernel(cupdlp_float *x,
                                           const cupdlp_float *ub,
                                           const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = x[i] > ub[i] ? ub[i] : x[i];
}

__global__ void element_wise_projSamelb_kernel(cupdlp_float *x,
                                               const cupdlp_float lb,
                                               const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = x[i] <= lb ? lb : x[i];
}

__global__ void element_wise_projSameub_kernel(cupdlp_float *x,
                                               const cupdlp_float ub,
                                               const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = x[i] >= ub ? ub : x[i];
}

__global__ void element_wise_initHaslb_kernal(cupdlp_float *haslb,
                                              const cupdlp_float *lb,
                                              const cupdlp_float bound,
                                              const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) haslb[i] = lb[i] > bound ? 1.0 : 0.0;
}

__global__ void element_wise_initHasub_kernal(cupdlp_float *hasub,
                                              const cupdlp_float *ub,
                                              const cupdlp_float bound,
                                              const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) hasub[i] = ub[i] < bound ? 1.0 : 0.0;
}

__global__ void element_wise_filterlb_kernal(cupdlp_float *x,
                                             const cupdlp_float *lb,
                                             const cupdlp_float bound,
                                             const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = lb[i] > bound ? lb[i] : 0.0;
}

__global__ void element_wise_filterub_kernal(cupdlp_float *x,
                                             const cupdlp_float *ub,
                                             const cupdlp_float bound,
                                             const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = ub[i] < bound ? ub[i] : 0.0;
}

__global__ void init_cuda_vec_kernal(cupdlp_float *x, const cupdlp_float val,
                                     const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) x[i] = val;
}

//xUpdate = x - dPrimalStep * (cost - ATy)
__global__ void primal_grad_step_kernal(cupdlp_float *xUpdate,
                                        const cupdlp_float *x,
                                        const cupdlp_float *cost,
                                        const cupdlp_float *ATy,
                                        const cupdlp_float dPrimalStep,
                                        const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) xUpdate[i] = x[i] - dPrimalStep * (cost[i] - ATy[i]);
}

//yUpdate = y + dDualStep * (b -2AxUpdate + Ax)
__global__ void dual_grad_step_kernal(cupdlp_float *yUpdate,
                                      const cupdlp_float *y,
                                      const cupdlp_float *b,
                                      const cupdlp_float *Ax,
                                      const cupdlp_float *AxUpdate,
                                      const cupdlp_float dDualStep,
                                      const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x; 
  if (i < len) yUpdate[i] = y[i] + dDualStep * (b[i] - 2 * AxUpdate[i] + Ax[i]);
}

// z = x - y
__global__ void naive_sub_kernal(cupdlp_float *z, const cupdlp_float *x,
                                  const cupdlp_float *y, const cupdlp_int len) {
  cupdlp_int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) z[i] = x[i] - y[i];
}