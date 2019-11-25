/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/ICrashUtil.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "ICrashUtil.h"

#include <algorithm>

#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsUtils.h"

void muptiplyByTranspose(const HighsLp& lp, const std::vector<double>& v,
                         std::vector<double>& result) {
  assert(result.size() == lp.numCol_);
  assert(v.size() == lp.numRow_);

  result.assign(lp.numCol_, 0);
  for (int col = 0; col < lp.numCol_; col++) {
    for (int k = lp.Astart_[col]; k < lp.Astart_[col + 1]; k++) {
      const int row = lp.Aindex_[k];
      result.at(col) += lp.Avalue_[k] * lp.rowUpper_[row];
    }
  }
}

void printMinorIterationDetails(const double iteration, const double col,
                                const double old_value, const double update,
                                const double ctx, const std::vector<double>& r,
                                const double quadratic_objective) {
  double rnorm = getNorm2(r);
  std::cout << "iter " << iteration;
  std::cout << ", col " << col;
  std::cout << ", update " << update;
  std::cout << ", old_value " << old_value;
  std::cout << ", new_value " << old_value + update;
  std::cout << ", ctx " << ctx;
  std::cout << ", r " << rnorm;
  std::cout << ", quadratic_objective " << quadratic_objective;
  std::cout << std::endl;
}

bool initialize(const HighsLp& lp, HighsSolution& solution,
                std::vector<double>& lambda) {
  if (!isSolutionConsistent(lp, solution)) {
    // clear and resize solution.
    solution.col_value.clear();
    solution.col_dual.clear();
    solution.row_value.clear();
    solution.row_dual.clear();

    solution.col_value.resize(lp.numCol_);
  }

  for (int col = 0; col < lp.numCol_; col++) {
    if (lp.colLower_[col] <= 0 && lp.colUpper_[col] >= 0)
      solution.col_value[col] = 0;
    else if (lp.colLower_[col] > 0)
      solution.col_value[col] = lp.colLower_[col];
    else if (lp.colUpper_[col] < 0)
      solution.col_value[col] = lp.colUpper_[col];
    else {
      HighsPrintMessage(
          ML_ALWAYS, "ICrash error: setting initial value for column %d", col);
      return false;
    }
  }

  lambda.resize(lp.numRow_);
  lambda.assign(lp.numRow_, 0);

  return true;
}

// Should only work with kLinearized type.
void minimize_exact_ica_admm() {
  // !!! do not modify LP cost ***
  // update();
}

void minimize_exact_penalty() {
  // update();
}

double minimize_component_ica(const int col, const double mu,
                              const std::vector<double>& lambda,
                              const HighsLp& lp, double& objective,
                              std::vector<double>& residual,
                              HighsSolution& sol) {
  // Minimize quadratic for column col.

  // Formulas for a and b when minimizing for x_j
  // a = (1/(2*mu)) * sum_i a_ij^2
  // b = -(1/(2*mu)) sum_i (2 * a_ij * (sum_{k!=j} a_ik * x_k - b_i)) + c_j (\)
  //     + sum_i a_ij * lambda_i
  // b / 2 = -(1/(2*mu)) sum_i (2 * a_ij
  double a = 0.0;
  double b = 0.0;

  for (int k = lp.Astart_[col]; k < lp.Astart_[col + 1]; k++) {
    int row = lp.Aindex_[k];
    a += lp.Avalue_[k] * lp.Avalue_[k];
    // matlab but with b = b / 2
    double bracket = -residual[row] - lp.Avalue_[k] * sol.col_value[col];
    bracket += lambda[row];
    // clp minimizing for delta_x
    // double bracket_clp = - residual_[row];
    b += lp.Avalue_[k] * bracket;
  }

  a = (0.5 / mu) * a;
  b = (0.5 / mu) * b + 0.5 * lp.colCost_[col];

  double theta = -b / a;
  double delta_x = 0;

  // matlab
  double new_x;
  if (theta > 0)
    new_x = std::min(theta, lp.colUpper_[col]);
  else
    new_x = std::max(theta, lp.colLower_[col]);
  delta_x = new_x - sol.col_value[col];

  // clp minimizing for delta_x
  // if (theta > 0)
  //   delta_x = std::min(theta, lp_.colUpper_[col] - col_value_[col]);
  // else
  //   delta_x = std::max(theta, lp_.colLower_[col] - col_value_[col]);

  sol.col_value[col] += delta_x;

  // std::cout << "col " << col << ": " << delta_x << std::endl;

  // Update objective, row_value, residual after each component update.
  objective += lp.colCost_[col] * delta_x;
  for (int k = lp.Astart_[col]; k < lp.Astart_[col + 1]; k++) {
    int row = lp.Aindex_[k];
    residual[row] -= lp.Avalue_[k] * delta_x;
    sol.row_value[row] += lp.Avalue_[k] * delta_x;
  }

  return delta_x;
}

void updateResidual(bool piecewise, const HighsLp& lp, const HighsSolution& sol,
                    std::vector<double>& residual) {
  residual.clear();
  residual.assign(lp.numRow_, 0);

  if (!piecewise) {
    assert(isEqualityProblem(lp));
    for (int row = 0; row < lp.numRow_; row++)
      residual[row] = std::fabs(lp.rowUpper_[row] - sol.row_value[row]);
  } else {
    // piecewise
    for (int row = 0; row < lp.numRow_; row++) {
      double value = 0;
      if (sol.row_value[row] <= lp.rowLower_[row])
        value = lp.rowLower_[row] - sol.row_value[row];
      else if (sol.row_value[row] >= lp.rowUpper_[row])
        value = sol.row_value[row] - lp.rowUpper_[row];

      residual[row] = value;
    }
  }
}

// Allows negative residuals
void updateResidualICA(const HighsLp& lp, const HighsSolution& sol,
                       std::vector<double>& residual) {
  assert(isEqualityProblem(lp));
  for (int row = 0; row < lp.numRow_; row++)
    residual[row] = lp.rowUpper_[row] - sol.row_value[row];
}
