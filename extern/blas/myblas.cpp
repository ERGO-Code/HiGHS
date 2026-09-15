#include "mycblas.h"

void highs_openblas_set_num_threads(int num_threads) {
#if defined(HIPO_USES_OPENBLAS)
  openblas_set_num_threads(num_threads);
#endif
}

int highs_openblas_get_num_threads(void) {
#if defined(HIPO_USES_OPENBLAS)
  return openblas_get_num_threads();
#else
  return -1;
#endif
}
