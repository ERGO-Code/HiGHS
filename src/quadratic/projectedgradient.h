#ifndef PROJECTEDGRADIENT_H
#define PROJECTEDGRADIENT_H

#include "HighsQp.h"
#include "HVector.h"
#include "SparseMatrix.h"

class ProjectedGradient {
 public:
  void solveQpPenalty(HighsQp& qp, double rho, double mu, double nu, HVector x0, HVector xbar);

  void solveLpPenalty(HighsLp& lp, double mu, HVector& x0);

 private:
  void computeGradientConstantPart(HVector& c, SparseMatrix& A, double mu, HVector& b, HVector& gradient);
  void computeGradient(HVector& gradientConstant, double mu, SparseMatrix& A, HVector& x, HVector& gradient);
  double computeObjectiveValue(HVector& c, SparseMatrix& A, double mu, HVector& b);
  void computeSecondDerivativeVector(HVector& x, SparseMatrix& A, double mu, HVector& result);
  
  void computeBreakpoints(HVector& gradient, HVector& x, HVector& u, HVector& l, std::vector<double>& t, std::vector<double>& tbar);
  void computeProjectedGradient(HVector& gradient, std::vector<double>& t, std::vector<double>& tbar, int bp, int numCol, HVector& projectedGradient);
  
  void projectGradient(HVector& gradient, HVector& x, HVector& l, HVector& u, int numCol, HVector& result);
  void projectIterate(HVector& x, HVector& l, HVector& u, int numCol);
};



#endif