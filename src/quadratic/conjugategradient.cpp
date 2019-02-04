#include "conjugategradient.h"
#include "HConst.h"

#include "HighsIO.h"
#include <math.h>

ConjugateGradient::~ConjugateGradient() {

}

ConjugateGradient::ConjugateGradient(SparseMatrix& A, HVector& b, const HVector& x0) {
  this->A = A;
  this->x.setup(3);
  this->x.copy(&x0);
  this->r.setup(A.numCol);
  this->b.setup(A.numCol);
  this->b.copy(&b);
  this->p.setup(A.numCol);
  this->rnew.setup(A.numCol);
}

void printVector(HVector& vec, const char* name) {
  HighsPrintMessage(ModelLogLevel::ML_DETAILED, "%s count: %d\n", name, vec.count);
  for(int i=0; i<vec.count; i++) {
    HighsPrintMessage(ModelLogLevel::ML_DETAILED, "%d %lf \n", vec.index[i],vec.array[vec.index[i]]);
  }
  HighsPrintMessage(ModelLogLevel::ML_DETAILED, "norm: %lf\n", sqrt(vec.norm2()));
}

void ConjugateGradient::iterate() {
  //alpha = ||rk||/(pkApk)
  HVector pkTA;
  pkTA.setup(this->A.numCol);
  this->A.vec_mat_prod(this->p, &pkTA);
  double normr = this->r.norm2();

  this->alpha = normr / pkTA.scalarProduct(&this->p);

  //xk+1 = xk + alpha * pk
  this->x.saxpy(this->alpha, &this->p);

  //tk+1 = rk + alpha * Apk
  this->r.saxpy(this->alpha, &pkTA);

  //beta = ||rk+1||/||rk||
  this->beta = this->r.norm2() / normr;

  //pk+1 = -rk+1 + beta*pk
  this->p.scale(beta);
  this->p.saxpy(-1.0, &this->r);

  //k = k + 1
  this->k++;
}

void ConjugateGradient::solve() {
  // initialize
  //r0 = Ax0 - b
  this->A.vec_mat_prod(this->x, &this->r);
  this->r.saxpy(-1.0, &this->b);

  //p0 = -r0
  this->p.saxpy(-1.0, &this->r);
  
  this->k = 0;

  printVector(this->x, "x");
  while(this->r.norm2() > HIGHS_CONST_TINY) {
    this->iterate();
    printVector(this->x, "x");
  }

  HighsPrintMessage(7, "Solved using Conjugate Gradient Method in %d iterations. Final residual: %lf.\n", this->k, this->r.norm2());
}