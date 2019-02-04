#include "projectedgradient.h"

#include "HConst.h"
#include "HighsIO.h"

#include <math.h>
#include <algorithm>

void printVector(HVector& vec, const char* name);

void ProjectedGradient::projectIterate(HVector& x, HVector& l, HVector& u, int numCol) {
  int nz = 0;
  for (int i=0; i < numCol; i++) {
     if (x.array[i] <= l.array[i] + HIGHS_CONST_TINY) {
       // variable at or under lower bound
       if (fabs(l.array[i]) < HIGHS_CONST_TINY) {
          x.array[i] = 0.0;
       } else {
          x.array[i] = l.array[i];
          x.index[nz] = i;
          nz++;
       }
     //} else if (x.array[i] >= u.array[i] - HIGHS_CONST_TINY) {
     //   // variable is at or above upper bound
     //   if (fabs(u.array[i]) < HIGHS_CONST_TINY) {
     //      x.array[i] = 0.0;
     //   } else {
     //     x.array[i] = u.array[i];
     //     x.index[nz] = i;
     //     nz++;
     //   }
     } else {
        // variable is between bounds
        if (fabs(x.array[i]) < HIGHS_CONST_TINY) {
           x.array[i] = 0.0;
        } else {
          x.index[nz] = i;
          nz++;
        }
     }
  }

  x.count = nz;
}

void ProjectedGradient::projectGradient(HVector& gradient, HVector& x,
                                        HVector& l, HVector& u, int numCol,
                                        HVector& result) {
  int nz = 0;
  for (int i = 0; i < numCol; i++) {
     if (x.array[i] <= l.array[i] + HIGHS_CONST_TINY && gradient.array[i] > 0.0) {
        // variable is at lower bounds
        result.array[i] = 0.0;
     } else if (x.array[i] >= u.array[i] - HIGHS_CONST_TINY && gradient.array[i] < 0.0) {
        // variable is at upper bound
        result.array[i] = 0.0;
     } else {
        // variable is within bounds
        if (fabs(gradient.array[i]) >= HIGHS_CONST_TINY) {
           result.index[nz] = i;
           result.array[i] = gradient.array[i];
           nz++;
        }
     }
  }
  result.count = nz;
}

// computes c - mu * A'b
void ProjectedGradient::computeGradientConstantPart(HVector& c, SparseMatrix& A,
                                                    double mu, HVector& b,
                                                    HVector& gradient) {
  if (mu > 0.0) {
    A.vec_mat_prod(b, &gradient);
    gradient.scale(-mu);
  }
  gradient.saxpy(1.0, &c);
}

// computes [c - mu * A'b] + mu * A'Ax
void ProjectedGradient::computeGradient(HVector& gradientConstant, double mu,
                                        SparseMatrix& A, HVector& x,
                                        HVector& gradient) {
  HVector temp(A.numRow);

  if (mu > 0.0) {
    A.mat_vec_prod(x, &temp);
    A.vec_mat_prod(temp, &gradient);
    gradient.scale(mu);
  }

  gradient.saxpy(1.0, &gradientConstant);
}


void ProjectedGradient::computeSecondDerivativeVector(HVector& x, SparseMatrix& A,
                                                double mu, HVector& result) {
  HVector temp(A.numRow);

  if (mu > 0.0) {
    A.mat_vec_prod(x, &temp);
    A.vec_mat_prod(temp, &result);
    result.scale(mu);
  }
}

void ProjectedGradient::computeBreakpoints(HVector& gradient, HVector& x,
                                           HVector& u, HVector& l,
                                           std::vector<double>& t,
                                           std::vector<double>& tbar) {
  // find breakpoints
  for (int i = 0; i < gradient.count; i++) {
    int index = gradient.index[i];
    double breakpoint = HIGHS_CONST_INF;
    if (gradient.array[index] < -HIGHS_CONST_TINY &&
        u.array[index] < HIGHS_CONST_INF) {
      breakpoint = (x.array[index] - u.array[index]) / gradient.array[index];
    } else if (gradient.array[index] > HIGHS_CONST_TINY &&
               l.array[index] > -HIGHS_CONST_INF) {
      breakpoint = (x.array[index] - l.array[index]) / gradient.array[index];
    }
    tbar.push_back(breakpoint);
    t.push_back(breakpoint);
  }
  t.push_back(0.0);
  t.push_back(HIGHS_CONST_INF);

  // sort breakpoints
  std::sort(t.begin(), t.end());

  // remove duplicates
  t.erase(unique(t.begin(), t.end()), t.end());
}

void ProjectedGradient::computeProjectedGradient(HVector& gradient,
                                                 std::vector<double>& t,
                                                 std::vector<double>& tbar,
                                                 int bp, int numCol,
                                                 HVector& projectedGradient) {
  // compute pj
  int nz = 0;
  for (int i = 0; i < numCol; i++) {
    if (t[bp] < tbar[i]) {
      if (fabs(gradient.array[i]) > HIGHS_CONST_TINY) {
        projectedGradient.array[i] = -gradient.array[i];
        projectedGradient.index[nz] = i;
        nz++;
      } else {
        // TODO: remove?
        projectedGradient.array[i] = 0.0;
      }
    }
  }
  projectedGradient.count = nz;
}

// solves c'x + mu/2 (Ax-b)'(Ax-b) s.t. l <= x <= u
void ProjectedGradient::solveLpPenalty(HighsLp& lp, double mu, HVector& x) {
  HighsPrintMessage(
      ML_VERBOSE,
      "solving qp arrising from lp-penalty method. Cols: %d Rows: %d\n",
      lp.numCol_, lp.numRow_);
  // data
  SparseMatrix A(lp.numCol_, lp.numRow_, lp.Astart_, lp.Aindex_, lp.Avalue_);
  HVector b(lp.rowUpper_, lp.numRow_);
  HVector l(lp.colLower_, lp.numCol_);
  HVector u(lp.colUpper_, lp.numCol_);
  HVector c(lp.colCost_, lp.numCol_);

  // temp variables
  HVector gradientConstant(lp.numCol_);
  HVector gradient(lp.numCol_);
  HVector breakpoints(lp.numCol_ + 1);

  this->computeGradientConstantPart(c, A, mu, b, gradientConstant);

  //printVector(gradientConstant, "constant gradient");

  int iteration = 0;
  while (iteration < 100000) {
     //this->projectIterate(x, l, u, lp.numCol_);
    HighsPrintMessage(ML_VERBOSE, "Iteration %d\n", iteration);
    // compute search direction
    gradient.setup(lp.numCol_);
    this->computeGradient(gradientConstant, mu, A, x, gradient);
    //printVector(gradient, "gradient");

    // project gradient?
    HVector projectedGradient(lp.numCol_);
    this->projectGradient(gradient, x, l, u, lp.numCol_, projectedGradient);
    //printVector(projectedGradient, "Projected Gradient");

    double norm = projectedGradient.norm2();
    HighsPrintMessage(ModelLogLevel::ML_DETAILED, "%d: norm %lf\n", iteration, norm);
    //if (norm < 10E-7) {
    //  break;
    //}
    if (projectedGradient.count == 0) {
       break;
    }

    std::vector<double> tbar;
    std::vector<double> t;
    this->computeBreakpoints(gradient, x, u, l, t, tbar);

    // TODO: maybe to t.size()-1?
    for (int bp = 0; bp < t.size(); bp++) {
      HVector p;
      p.setup(lp.numCol_);

      this->computeProjectedGradient(gradient, t, tbar, bp, lp.numCol_, p);

      // compute hessian
      HVector Hp_vec(lp.numCol_);
      this->computeSecondDerivativeVector(p, A, mu, Hp_vec);

      // compute fd
      double fd = gradientConstant.scalarProduct(&p) + x.scalarProduct(&Hp_vec);
      HighsPrintMessage(ML_VERBOSE, "fd: %lf\n", fd);

      if (fd > HIGHS_CONST_TINY) {
        // first local minimum is at the current iterate.
        break;
      } else {
        // compute fdd
        double fdd = p.scalarProduct(&Hp_vec);
        HighsPrintMessage(ML_VERBOSE, "fdd: %lf\n", fdd);

        if (fabs(fdd) < HIGHS_CONST_TINY) {
           break;
        }

        double t_star = -fd / fdd;
        HighsPrintMessage(ML_VERBOSE, "tstar: %lf\n", t_star);
        if (t_star < t[bp + 1] - t[bp]) {
          x.saxpy(t_star, &p);
          break;
        } else {
          x.saxpy(t[bp + 1] - t[bp], &p);
        }
      }
    }

    // search in subspace. Variables at bounds may not be changed.
    // solve subproblem (approximately):
    // min  0.5 mu *x'A'Ax + c'x b mu * A'b
    // s.t. x_i = x_i^c, i \in Active set(x^c)
    // l <= x <= u

    printVector(x, "x");
    iteration++;
  }
  
}