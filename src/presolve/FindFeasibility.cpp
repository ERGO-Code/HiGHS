#include "FindFeasibility.h"
#include "io/HighsIO.h"

struct Quadratic
{
  double objective;
  double residual_norm_1;
  double residual_norm_2;

  vector<double> residual;
};

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

  mu = 1;

  lambda.resize(lp.numRow_);
  lambda.assign(lp.numRow_, 0);

  return HighsStatus::OK;
}

HighsStatus runIdiot(const HighsLp& lp, HighsSolution& solution) {
  // Initialize x_0 ≥ 0, μ_1, λ_1 = 0.
  double mu;
  std::vector<double> lambda;

  HighsStatus status = initialize(lp, solution, mu, lambda);

  int K = 15;
  for (int k = 0; k < K; k++)
    // Minimize quadratic function.

    // Possibly update mu.

    // Else update lambda.

    return HighsStatus::OK;
}