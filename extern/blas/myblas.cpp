#include "mycblas.h"

int highs_openblas_set_num_threads(int num_threads) {
#if defined(HIPO_USES_OPENBLAS)
  openblas_set_num_threads(num_threads);
  return openblas_get_num_threads();
#else
  return -1;
#endif
}
