// Example 2. Application Using C and cuBLAS: 0-based indexing
//-----------------------------------------------------------
#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include "cublas_v2.h"
#define N 5

int main(void) {
  cudaError_t cudaStat;
  cublasStatus_t stat;
  cublasHandle_t handle;
  int i;
  float* devPtrA;
  float* a = 0;
  a = (float*)malloc(N * sizeof(*a));
  if (!a) {
    printf("host memory allocation failed");
    return EXIT_FAILURE;
  }

  for (i = 0; i < N; i++) a[i] = (float)(i * N);
  for (i = 0; i < N; i++) printf("%7.0f", a[(i)]);
  std::cout << std::endl;

  cudaStat = cudaMalloc((void**)&devPtrA, N * sizeof(*a));
  if (cudaStat != cudaSuccess) {
    printf("device memory allocation failed");
    free(a);
    return EXIT_FAILURE;
  }
  stat = cublasCreate(&handle);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf("CUBLAS initialization failed\n");
    free(a);
    cudaFree(devPtrA);
    return EXIT_FAILURE;
  }

  stat = cublasSetVector(N, sizeof(*a), a, 1, devPtrA, 1);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf("data download failed");
    free(a);
    cudaFree(devPtrA);
    cublasDestroy(handle);
    return EXIT_FAILURE;
  }

  float r;

  cublasStatus_t status = cublasSnrm2(handle, N, devPtrA, 1, &r);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf("error at dnorm");
    free(a);
    cudaFree(devPtrA);
    cublasDestroy(handle);
    return EXIT_FAILURE;
  }

  std::cout << "Snorm: " << r << std::endl << std::endl;

  float* b;
  b = (float*)malloc(N * sizeof(*b));

  stat = cublasGetVector(N, sizeof(*a), devPtrA, 1, b, 1);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf("data upload failed");
    free(a);
    cudaFree(devPtrA);
    cublasDestroy(handle);
    return EXIT_FAILURE;
  }

  cudaFree(devPtrA);
  cublasDestroy(handle);

  for (i = 0; i < N; i++) {
    printf("%7.0f", b[(i)]);
  }
  printf("\n");

  free(a);
  return EXIT_SUCCESS;
}