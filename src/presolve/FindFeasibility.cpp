#include "FindFeasibility.h"

#include <sstream>

#include "io/HighsIO.h"

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
    solution.colValue_ = col_value_;
    solution.rowValue_ = row_value_;
  }

  void minimize_by_component(const double mu,
                             const std::vector<double>& lambda);

 private:
  const HighsLp& lp_;
  vector<double> col_value_;
  vector<double> row_value_;

  double objective_;
  double residual_norm_1_;
  double residual_norm_2_;
  vector<double> residual_;

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


void Quadratic::minimize_by_component(const double mu,
                                      const std::vector<double>& lambda) {
  int iterations = 15;

  // Formulas for a and b when minimizing for x_j
  // a = (1/(2*mu)) * sum_i a_ij^2
  // b = -(1/(2*mu) sum_i (2 * a_ij * (sum_{k!=j} a_ik * x_k - b_i)) + c_j \
  //     + sum_i a_ij * lambda_i

  double a = 0.0;
  double b = 0.0;

  for (int iteration = 0; iteration < iterations; iteration++) {
    for (int col = 0; col < lp_.numCol_; col++) {
      // todo: determine whether to minimize for col.
      // Minimize quadratic for column col.

      for (int k = lp_.Astart_[col]; k < lp_.Astart_[col+1]; k++) {
        int row = lp_.Aindex_[k];
        a += lp_.Avalue_[k] * lp_.Avalue_[k];
      // matlab
        double bracket = 2 * residual_[row] - lp_.Avalue_[k] * col_value_[col];
      // clp
      // double bracket = - residual[row];
        b += lp_.Avalue_[k] * bracket;
      }

      a = (0.5 / mu) * a;
      b = (0.5 / mu) * b + 0.5 * lp_.colCost_[col];

      double theta = -b / a;
      double delta_x = 0;

      if (theta > 0)
        delta_x = std::min(theta, lp_.colUpper_[col] - col_value_[col]);
      else
        delta_x = std::max(theta, lp_.colLower_[col] - col_value_[col]);

      col_value_[col] += delta_x;

      // Update objective, row_value, residual after each component update.
      objective_[col] += lp_.colCost_[col] * delta_x;
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
}

double chooseStartingMu(const HighsLp& lp) {
  return 1000;
}


HighsStatus initialize(const HighsLp& lp,
                       HighsSolution& solution,
                       double& mu,
                       std::vector<double>& lambda)
{
  if (!isSolutionConsistent(lp, solution)) {
    // clear and resize solution.
    solution.colValue_.clear();
    solution.colDual_.clear();
    solution.rowValue_.clear();
    solution.rowDual_.clear();

    solution.colValue_.resize(lp.numCol_);
  }

  for (int col = 0; col < lp.numCol_; col++) {
    if (lp.colLower_[col] <= 0 && lp.colUpper_[col] >= 0)
      solution.colValue_[col] = 0;
    else if (lp.colLower_[col] > 0)
      solution.colValue_[col] = lp.colLower_[col];
    else if (lp.colUpper_[col] < 0)
      solution.colValue_[col] = lp.colUpper_[col];
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

HighsStatus runFeasibility(const HighsLp& lp, HighsSolution& solution) {
  if (!isEqualityProblem(lp))
    return HighsStatus::NotImplemented;

  // Initialize x_0 ≥ 0, μ_1, λ_1 = 0.
  double mu;
  std::vector<double> lambda;

  HighsStatus status = initialize(lp, solution, mu, lambda);
  Quadratic quadratic(lp, solution.colValue_);

  int K = 10;
  for (int iteration = 0; iteration < K; iteration++) {
    // Minimize quadratic function.
    quadratic.minimize_by_component(mu, lambda);
    std::stringstream ss;

    // Report outcome.
    double residual_norm_2 = quadratic.getResidualNorm2();
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
