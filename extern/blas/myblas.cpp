#include "mycblas.h"

void highs_openblas_set_num_threads(int num_threads) {
#if defined(HIPO_USES_OPENBLAS)
  openblas_set_num_threads(num_threads);
#endif
}
