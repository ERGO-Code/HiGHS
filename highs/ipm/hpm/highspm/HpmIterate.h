#ifndef HIPO_ITERATE_H
#define HIPO_ITERATE_H

#include <vector>

#include "HpmModel.h"
#include "ipm/hpm/auxiliary/IntConfig.h"
#include "HpmInfo.h"

namespace hipo {

// Holds the Newton direction Delta(x,y,xl,xu,zl,zu)
struct NewtonDir {
  std::vector<double> x{};
  std::vector<double> y{};
  std::vector<double> xl{};
  std::vector<double> xu{};
  std::vector<double> zl{};
  std::vector<double> zu{};

  NewtonDir(Int m, Int n);
};

struct Iterate {
  // lp model
  const Model* model;

  // ipm point
  std::vector<double> x, xl, xu, y, zl, zu;

  // residuals
  std::vector<double> res1, res2, res3, res4, res5, res6;

  // Newton direction
  NewtonDir delta;

  // indicators
  double pobj, dobj, pinf, dinf, pdgap;

  double mu;
  std::vector<double> scaling;

  // smallest value of mu seen so far
  double best_mu;

  // number of small/large complementarity products
  Int num_small, num_large;

  // ===================================================================================
  // Functions to construct, clear and check for nan or inf
  // ===================================================================================
  Iterate(const Model& model_input);

  // clear existing data
  void clearIter();
  void clearRes();
  void clearDir();

  // check if any component is nan or infinite
  bool isNan() const;
  bool isInf() const;
  bool isResNan() const;
  bool isResInf() const;
  bool isDirNan() const;
  bool isDirInf() const;

  // ===================================================================================
  // Compute:
  //  mu = \sum xl(i) * zl(i) + \sum xu(j) * zu(j)
  // for variables: i with a finite lower bound
  //                j with a finite upper bound
  // ===================================================================================
  void computeMu();

  // ===================================================================================
  // Compute diagonal scaling Theta^{-1}
  //  Theta^{-1}_{ii} = zl(i) / xl(i) + zu(i) / xu(i)
  // Theta^{-1} only considers the terms above if the corresponding upper/lower
  // bound is finite.
  // ===================================================================================
  void computeScaling();

  // ===================================================================================
  // Compute convergence indicators
  // - primal and dual objectives, and primal-dual relative gap
  // - primal and dual infeasibilities (either scaled or unscaled)
  // - complementairy products
  // ===================================================================================
  void indicators();

  // objectives
  void primalObj();
  void dualObj();
  void pdGap();

  // infeasibilities
  void primalInfeas();
  void dualInfeas();
  void primalInfeasUnscaled();
  void dualInfeasUnscaled();
  double infeasAfterDropping() const;

  // complementarity products
  void products();

  // ===================================================================================
  // Compute:
  //  res1 = rhs - A * x
  //  res2 = lower - x + xl
  //  res3 = upper - x - xu
  //  res4 = c - A^T * y - zl + zu
  // Components of residuals 2,3 are set to zero if the corresponding
  // upper/lower bound is not finite.
  // ===================================================================================
  void residual1234();

  // ===================================================================================
  // Compute:
  //  res5 = sigma * mu * e - Xl * Zl * e
  //  res6 = sigma * mu * e - Xu * Zu * e
  // Components of residuals 5,6 are set to zero if the corresponding
  // upper/lower bound is not finite.
  // ===================================================================================
  void residual56(double sigma);

  // ===================================================================================
  // Compute:
  //  res7 = res4 - Xl^{-1} * (res5 + Zl * res2) + Xu^{-1} * (res6 - Zu * res3)
  // (the computation of res7 takes into account only the components for which
  // the correspoding upper/lower bounds are finite)
  // ===================================================================================
  std::vector<double> residual7() const;

  // ===================================================================================
  // Compute:
  //  res8 = res1 + A * Theta * res7
  // ===================================================================================
  std::vector<double> residual8(const std::vector<double>& res7) const;

  // ===================================================================================
  // Extract solution to be returned to user:
  // - remove extra slacks from x, xl, xu, zl, zu
  // - adjust sign of y for inequality constraints
  // - compute and adjust sign of slacks
  // ===================================================================================
  void extract(std::vector<double>& x_user, std::vector<double>& xl_user,
               std::vector<double>& xu_user, std::vector<double>& slack_user,
               std::vector<double>& y_user, std::vector<double>& zl_user,
               std::vector<double>& zu_user) const;

  // ===================================================================================
  // Extract complementary solution to be used for crossover with IPX:
  // - drop variables to obtain complementary (x,y,z)
  // - adjust y based on z-slacks
  // - compute slacks
  // - remove extra slacks from x, z
  // ===================================================================================
  void extract(std::vector<double>& x_user, std::vector<double>& slack_user,
               std::vector<double>& y_user, std::vector<double>& z_user) const;

  // ===================================================================================
  // Construct a complementary point (x,y,z), such that for each j, either xj is
  // at one of the bounds (lower or upper), or zj is zero.
  // ===================================================================================
  void dropToComplementarity(std::vector<double>& x_cmp,
                             std::vector<double>& y_cmp,
                             std::vector<double>& z_cmp) const;

  // ===================================================================================
  // Compute residuals after solution has been found, postprocessed and
  // unscaled.
  // ===================================================================================
  void finalResiduals(Info& info) const;
};

}  // namespace hipo

#endif