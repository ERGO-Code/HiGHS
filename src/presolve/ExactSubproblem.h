#ifndef PRESOVLE_EXACTSUBPROBLEM_H_
#define PRESOVLE_EXACTSUBPROBLEM_H_

#include "io/HighsIO.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HConst.h"
#include "quadratic/projectedgradient.h"
#include "simplex/HVector.h"

void solve_exact(const HighsLp& lp, const double mu, std::vector<double>& col_value) {
  HVector vector(col_value, col_value.size());

  HighsLp& lp_non_const = const_cast<HighsLp&>(lp);

  HighsOptions options;
  options.messageLevel = ML_NONE;
  HighsSetIO(options);

  ProjectedGradient projected_gradient;
  projected_gradient.solveLpPenalty(lp_non_const, mu, vector);

  options.messageLevel = ML_MINIMAL;
  HighsSetIO(options);

  col_value = vector.array;
}

#endif