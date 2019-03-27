#include "FindFeasibility.h"

#include <algorithm>
#include <sstream>
#include <cmath>
#include <iomanip>

#include "io/HighsIO.h"
#include "presolve/ExactSubproblem.h"

constexpr double kExitTolerance = 0.00000001;

bool isEqualityProblem(const HighsLp& lp) {
  for (int row = 0; row < lp.numRow_; row++)
    if (lp.rowLower_[row] != lp.rowUpper_[row])
      return false;

  return true;
}

class Quadratic
{
 public:
  Quadratic(const HighsLp& lp,
            std::vector<double>& primal_values) :
            lp_(lp), col_value_(primal_values) { update(); }

  const std::vector<double>& getResidual() const { return residual_; }
  double getResidualNorm2() const { return residual_norm_2_; }
  double getObjective() const { return objective_; }

  void getSolution(HighsSolution& solution) const {
    solution.col_value = col_value_;
    solution.row_value = row_value_;

    // check what solution looks like
    double max = *std::max_element(col_value_.begin(), col_value_.end());
    double min = *std::min_element(col_value_.begin(), col_value_.end());

    HighsPrintMessage(ML_ALWAYS, "\n");
    HighsPrintMessage(ML_ALWAYS, "Solution max element: %4.3f\n", max);
    HighsPrintMessage(ML_ALWAYS, "Solution min element: %4.3f\n", min);
  }

  void minimize_by_component(const double mu,
                             const std::vector<double>& lambda);

  void minimize_exact_penalty(const double mu);

 private:
  const HighsLp& lp_;
  std::vector<double> col_value_;
  std::vector<double> row_value_;

  double objective_;
  double residual_norm_1_;
  double residual_norm_2_;
  std::vector<double> residual_;

  void updateObjective();
  void updateRowValue();
  void updateResidual();

  void update();
};

void Quadratic::update() {
    updateObjective();
    updateRowValue();
    updateResidual();
}

void Quadratic::updateRowValue() {
  row_value_.clear();
  row_value_.assign(lp_.numRow_, 0);

  for (int col = 0; col < lp_.numCol_; col++) {
    for (int k = lp_.Astart_[col]; k < lp_.Astart_[col+1]; k++) {
      int row = lp_.Aindex_[k];
      row_value_[row] += lp_.Avalue_[k] * col_value_[col];
    }
  }
}

void Quadratic::updateResidual() {
  residual_.clear();
  residual_.assign(lp_.numRow_, 0);
  residual_norm_1_ = 0;
  residual_norm_2_ = 0;

  for (int row = 0; row  < lp_.numRow_; row++) {
    // for the moment assuming rowLower == rowUpper
    residual_[row] = lp_.rowUpper_[row] - row_value_[row];

    residual_norm_1_ += std::fabs(residual_[row]);
    residual_norm_2_ += residual_[row] * residual_[row];
  }

  residual_norm_2_ = std::sqrt(residual_norm_2_);
}

void Quadratic::updateObjective() {
  objective_ = 0;
  for (int col = 0; col < lp_.numCol_; col++)
    objective_ += lp_.colCost_[col] * col_value_[col];
}

void Quadratic::minimize_exact_penalty(const double mu) {
  double mu_penalty = 1.0 / mu;
  const HighsLp& lp_ref = lp_;

  solve_exact(lp_ref, mu_penalty, col_value_);

  update();
}

void Quadratic::minimize_by_component(const double mu,
                                      const std::vector<double>& lambda) {
  int iterations = 100;

  for (int iteration = 0; iteration < iterations; iteration++) {
    for (int col = 0; col < lp_.numCol_; col++) {
      // todo: determine whether to minimize for col.
      // Minimize quadratic for column col.

      // Formulas for a and b when minimizing for x_j
      // a = (1/(2*mu)) * sum_i a_ij^2
      // b = -(1/(2*mu)) sum_i (2 * a_ij * (sum_{k!=j} a_ik * x_k - b_i)) + c_j \
      //     + sum_i a_ij * lambda_i
      // b / 2 = -(1/(2*mu)) sum_i (2 * a_ij
      double a = 0.0;
      double b = 0.0;

      for (int k = lp_.Astart_[col]; k < lp_.Astart_[col+1]; k++) {
        int row = lp_.Aindex_[k];
        a += lp_.Avalue_[k] * lp_.Avalue_[k];
        // matlab but with b = b / 2
        double bracket = - residual_[row] - lp_.Avalue_[k] * col_value_[col];
        bracket += lambda[row];
        // clp minimizing for delta_x
        // double bracket_clp = - residual_[row];
        b += lp_.Avalue_[k] * bracket;
      }

      a = (0.5 / mu) * a;
      b = (0.5 / mu) * b + 0.5 * lp_.colCost_[col];

      double theta = -b / a;
      double delta_x = 0;

      // matlab
      double new_x;
      if (theta > 0)
        new_x = std::min(theta, lp_.colUpper_[col]);
      else
        new_x = std::max(theta, lp_.colLower_[col]);
      delta_x = new_x - col_value_[col];

      // clp minimizing for delta_x
      // if (theta > 0)
      //   delta_x = std::min(theta, lp_.colUpper_[col] - col_value_[col]);
      // else
      //   delta_x = std::max(theta, lp_.colLower_[col] - col_value_[col]);

      col_value_[col] += delta_x;

      // Update objective, row_value, residual after each component update.
      objective_ += lp_.colCost_[col] * delta_x;
      for (int k = lp_.Astart_[col]; k < lp_.Astart_[col+1]; k++) {
        int row = lp_.Aindex_[k];
        residual_[row] -= lp_.Avalue_[k] * delta_x;
        row_value_[row] += lp_.Avalue_[k] * delta_x;
      }
    }

    // Code below commented out because updating after each component
    // minimization.
    // update();

    // updateResidual();
    // todo: check for early exit
  }
  update();
}

double chooseStartingMu(const HighsLp& lp) {
  return 0.001;
}


HighsStatus initialize(const HighsLp& lp,
                       HighsSolution& solution,
                       double& mu,
                       std::vector<double>& lambda)
{
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
      HighsLogMessage(HighsMessageType::ERROR,
                      "Error setting initial value for column %d", col);
      return HighsStatus::Error;
    }
  }

  mu = chooseStartingMu(lp);

  lambda.resize(lp.numRow_);
  lambda.assign(lp.numRow_, 0);

  return HighsStatus::OK;
}

HighsStatus runFeasibility(const HighsLp& lp,
                           HighsSolution& solution,
                           const MinimizationType type) {
  if (!isEqualityProblem(lp))
    return HighsStatus::NotImplemented;

  // Initialize x_0 ≥ 0, μ_1, λ_1 = 0.
  double mu;
  std::vector<double> lambda;

  HighsStatus status = initialize(lp, solution, mu, lambda);
  Quadratic quadratic(lp, solution.col_value);

  if (type == MinimizationType::kComponentWise)
    HighsPrintMessage(ML_ALWAYS, "Minimizing quadratic subproblem component-wise...\n");
  else if (type == MinimizationType::kExact)
    HighsPrintMessage(ML_ALWAYS, "Minimizing quadratic subproblem exactly...\n");

  // Report values at start.
  std::stringstream ss;
  double residual_norm_2 = quadratic.getResidualNorm2();
  ss << "Iteration " << std::setw(3) << 0 << ": objective " << std::setw(3)
      << std::fixed << std::setprecision(2)
      << quadratic.getObjective() << " residual " << std::setw(5)
      << std::scientific << quadratic.getResidualNorm2() << std::endl;
  HighsPrintMessage(ML_ALWAYS, ss.str().c_str());

  // Minimize approximately for K iterations.
  int K = 30;
  for (int iteration = 1; iteration < K + 1; iteration++) {
    // Minimize quadratic function.
    if (type == MinimizationType::kComponentWise)
      quadratic.minimize_by_component(mu, lambda);
    else if (type == MinimizationType::kExact)
      quadratic.minimize_exact_penalty(mu);

    // Report outcome.
    residual_norm_2 = quadratic.getResidualNorm2();
    ss.str(std::string());
    ss << "Iteration " << std::setw(3) << iteration << ": objective " << std::setw(3)
       << std::fixed << std::setprecision(2)
       << quadratic.getObjective() << " residual " << std::setw(5)
       << std::scientific << residual_norm_2 << std::endl;
    HighsPrintMessage(ML_ALWAYS, ss.str().c_str());

    // Exit if feasible.
    if (residual_norm_2 < kExitTolerance)
      break;

    // Update mu every third iteration, otherwise update lambda.
    if (iteration % 3 == 2) {
      mu = 0.1 * mu;
    } else {
      lambda = quadratic.getResidual();
      for (int row = 0; row < lp.numRow_; row++)
        lambda[row] = mu * lambda[row];
    }
  }

  quadratic.getSolution(solution);
  HighsPrintMessage(ML_ALWAYS,
                    "\nSolution set at the end of feasibility search.\n");

  return HighsStatus::OK;
}
