#include "projectedgradient.h"

#include "HConst.h"
#include "HighsIO.h"

#include <math.h>
#include <algorithm>

#define TOLERANCE_GRADIENT_ZERO 10E-2
#define TOLERANCE_CG_OBJECTIVE_CHANGE 10E-4
#define TOLERANCE_CG_GRADIENT_SMALL 10E-2
#define LIMIT_PGM_ITERATIONS 200

double ProjectedGradient::computeObjectiveValue(HVector& c, SparseMatrix& A, double mu, HVector& b, HVector& x) {
  double obj = 0.0;

  obj += c.scalarProduct(&x);

  if (mu > 0.0) {
    HVector temp(A.numRow);
    A.mat_vec_prod(x, &temp);
    temp.saxpy(-1.0, &b);
    obj += mu/2 * temp.norm2();
  }

  return obj;
}

void ProjectedGradient::projectGradient(double sign, HVector& gradient,
                                        HVector& x, HVector& l, HVector& u,
                                        int numCol, HVector& result) {
  int nz = 0;
  for (int i = 0; i < gradient.count; i++) {
    int index = gradient.index[i];
    if (x.array[index] <= l.array[index] + HIGHS_CONST_TINY) {
      result.array[index] = sign * fmin(0.0, sign * gradient.array[index]);
    } else if (x.array[index] >= u.array[index] - HIGHS_CONST_TINY) {
      result.array[index] = sign * fmax(0.0, sign * gradient.array[index]);
    } else {
      result.array[index] = gradient.array[index];
    }
    if (fabs(result.array[index]) > HIGHS_CONST_TINY) {
      result.index[nz] = index;
      nz++;
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

void ProjectedGradient::computeSecondDerivativeVector(HVector& x,
                                                      SparseMatrix& A,
                                                      double mu,
                                                      HVector& result) {
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
        projectedGradient.array[i] = HIGHS_CONST_ZERO;
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

  this->computeGradientConstantPart(c, A, mu, b, gradientConstant);

  int iteration = 0;
  int matlab_iterations = 0;
  while (iteration < LIMIT_PGM_ITERATIONS) {
    // compute search direction
    gradient.setup(lp.numCol_);
    this->computeGradient(gradientConstant, mu, A, x, gradient);

    // project gradient?
    HVector projectedGradient(lp.numCol_);
    this->projectGradient(1.0, gradient, x, l, u, lp.numCol_,
                          projectedGradient);
    double norm = projectedGradient.norm2();
    HighsPrintMessage(ML_MINIMAL, "%d: gradient norm %lf\n",
                      iteration, norm);

    // TODO: I feel like this should not be a fixed value for all problems if used in the context of idiot crash
    // maybe divide norm by the number of columns?
    if (norm < TOLERANCE_GRADIENT_ZERO) {
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

      if (fd > HIGHS_CONST_TINY) {
        // first local minimum is at the current iterate.
        break;
      } else {
        // compute fdd
        double fdd = p.scalarProduct(&Hp_vec);

        if (fabs(fdd) < HIGHS_CONST_TINY) {
          break;
        }

        double t_star = -fd / fdd;
        if (t_star < t[bp + 1] - t[bp]) {
          // x.saxpy(t_star, &p);
          break;
        } else {
          x.saxpy(t[bp + 1] - t[bp], &p);
        }
      }
    }
    matlab_iterations++;
    // search in subspace. Variables at bounds may not be changed.
    // solve subproblem (approximately):
    // min  0.5 mu *x'A'Ax + c'x b mu * A'b
    // s.t. x_i = x_i^c, i \in Active set(x^c)
    // l <= x <= u

    // SPECULATIVE CODE START
    gradient.setup(lp.numCol_);
    this->computeGradient(gradientConstant, mu, A, x, gradient);
    this->projectGradient(1.0, gradient, x, l, u, lp.numCol_, gradient);

    HVector sk(lp.numCol_);
    sk.saxpy(-1.0, &gradient);

    int cgIteration = 0;
    double prevObjValue = this->computeObjectiveValue(c, A, mu, b, x);
    while (cgIteration < lp.numCol_) {
      matlab_iterations++;
      double norm_sk = sk.norm2();
      if (norm_sk < HIGHS_CONST_TINY) {
        HighsPrintMessage(ML_MINIMAL, "zero search direction in CG\n");
        break;
      }

      HVector w(lp.numRow_);
      A.mat_vec_prod(sk, &w);
      w.tight();
      double q1 = mu * w.norm2();
      double q2 = sk.scalarProduct(&gradient);

      if (fabs(q1) < HIGHS_CONST_TINY) {
        break;
      }
      double alpha = -q2 / q1;

      if (alpha <= HIGHS_CONST_TINY) {
        HighsPrintMessage(ML_MINIMAL, "alpha in CG too small %lf\n", alpha);
        break;
      }

      x.tight();
      x.saxpy(alpha, &sk);
      // project x on bounds
      for (int i = 0; i < lp.numCol_; i++) {
        x.array[i] = fmin(fmax(l.array[i], x.array[i]), u.array[i]);
      }

      double norm_gk = gradient.norm2();
      this->computeGradient(gradientConstant, mu, A, x, gradient);
      this->projectGradient(1.0, gradient, x, l, u, lp.numCol_, gradient);

      double norm_gkp1 = gradient.norm2();
      if (norm_gkp1 < TOLERANCE_CG_GRADIENT_SMALL) {
        HighsPrintMessage(ML_MINIMAL, "CG: norm of gradient too small: %lf\n", norm_gkp1);
        break;
      }

      // if objective change is too small, stop CG and continue with PGM 
      double newObjValue = this->computeObjectiveValue(c, A, mu, b, x);
      double objChange = fabs(prevObjValue - newObjValue);
      if (objChange < TOLERANCE_CG_OBJECTIVE_CHANGE) {
        HighsPrintMessage(ML_MINIMAL, "CG iter %d: objective change too small: %lf\n", cgIteration, objChange);
        break;
      }
      prevObjValue = newObjValue;

      double beta = (norm_gkp1) / (norm_gk);
      sk.scale(beta);
      sk.tight();
      sk.saxpy(-1.0, &gradient);
      projectGradient(-1.0, sk, x, l, u, lp.numCol_, sk); // replaced gradient by sk

      if (sk.scalarProduct(&gradient) > HIGHS_CONST_TINY) {
        sk.copy(&gradient);
        sk.scale(-1.0);
      }

      cgIteration++;
    }

    // SPECULATIVE CODE END
    iteration++;
  }
  // x.peak();
  HighsPrintMessage(ML_MINIMAL, "Matlab Iterations %d\n", matlab_iterations);
}
