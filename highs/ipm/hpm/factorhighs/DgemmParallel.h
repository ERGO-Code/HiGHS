#ifndef FACTORHIGHS_DGEMM_PARALLEL_H
#define FACTORHIGHS_DGEMM_PARALLEL_H

#include "ipm/hpm/auxiliary/IntConfig.h"

namespace highspm {

// parallelise dgemm for use within factorisation
// Performs Q <- Q - R P^T in hybrid format.
// Parallelised over the rows of R and Q.
class dgemmParalleliser {
  const double* P_;
  const double* R_;
  double* Q_;
  const Int col_;
  const Int jb_;

 public:
  dgemmParalleliser(const double* P, const double* R, double* Q, Int col,
                    Int jb);

  void run(Int start, Int end, double beta) const;
};

void dgemmParallel(const double* P, const double* R, double* Q, Int col, Int jb,
                   Int row, Int nb, double beta = 1.0);

}  // namespace highspm

#endif