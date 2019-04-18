#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "HMatrix.h"
#include "HVector.h"
#include "HighsQp.h"

#include "SparseMatrix.h"

class ConjugateGradient {
 private:

  SparseMatrix A;
  HVector b;

  HVector x;
  HVector p;
  HVector r;
  HVector rnew;

  double alpha;
  double beta;

  int k;

 public:
  ConjugateGradient(SparseMatrix& A, HVector& b, const HVector& x0);
  ~ConjugateGradient();

  void solve();

  void iterate();
};

#endif