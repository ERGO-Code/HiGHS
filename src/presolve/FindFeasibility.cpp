#include "FindFeasibility.h"
#include "io/HighsIO.h"

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
            lp_(lp), col_value(primal_values) {
              updateObjective();
              updateRowValue();
              updateResidual();
            }

  void setSolution(const std::vector<double>& values) {
    col_value = (values);
    updateObjective();
    updateRowValue();
    updateResidual();
  }

 private:
  const HighsLp& lp_;
  vector<double> col_value;
  vector<double> row_value;

  double objective;
  double residual_norm_1;
  double residual_norm_2;
  vector<double> residual;

  void updateObjective();
  void updateRowValue();
  void updateResidual();
};

void Quadratic::updateRowValue() {
  row_value.clear();
  row_value.assign(lp_.numRow_, 0);

  for (int col = 0; col < lp_.numCol_; col++) {
    for (int k = lp_.Astart_[col]; k < lp_.Astart_[col+1]; k++) {
      int row = lp_.Aindex_[k];
      row_value[row] += lp_.Avalue_[k] * col_value[col];
    }
  }
}

void Quadratic::updateResidual() {
  residual.clear();
  residual.assign(lp_.numRow_, 0);
  residual_norm_1 = 0;
  residual_norm_2 = 0;

  for (int row = 0; row  < lp_.numRow_; row++) {
    // for the moment assuming rowLower == rowUpper
    residual[row] = lp_.rowUpper_[row] - row_value[row];

    residual_norm_1 += std::fabs(residual[row]);
    residual_norm_2 += residual[row] * residual[row];
  }

  residual_norm_2 = std::sqrt(residual_norm_2);
}

void Quadratic::updateObjective() {
  objective = 0;
  for (int col = 0; col < lp_.numCol_; col++)
    objective += lp_.colCost_[col] * col_value[col];
}

double chooseStartingMu(const HighsLp& lp) {
  return 1;
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

  int K = 15;
  for (int k = 0; k < K; k++)
    // Minimize quadratic function.

    // Possibly update mu.

    // Else update lambda.

    return HighsStatus::OK;
}

