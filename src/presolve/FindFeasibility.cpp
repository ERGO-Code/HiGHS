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
            lp_(lp), col_value(primal_values) { update(); }

  const std::vector<double>& getResidual() const { return residual; }

  void getSolution(HighsSolution& solution) const {
    solution.colValue_ = col_value;
    solution.rowValue_ = row_value;
  }

  void minimize_by_component(const double mu,
                             const std::vector<double>& lambda);

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

  void update();
};

void Quadratic::update() {
    updateObjective();
    updateRowValue();
    updateResidual();
}

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


void Quadratic::minimize_by_component(const double mu,
                                      const std::vector<double>& lambda) {
  int iterations = 5;

  // todo: check notes for formulas for a and b and x.
  double a = 0.0;
  double b = 0.0;

  for (int iteration = 0; iteration < iterations; iteration++) {
    for (int col = 0; col < lp_.numCol_; col++) {
      // todo: determine whether to minimize for col.
      // Minimize quadratic for column col.
      for (int k = lp_.Astart_[col]; k < lp_.Astart_[col+1]; k++) {
        int row = lp_.Aindex_[k];
        a += lp_.Avalue_[k] * lp_.Avalue_[k];
        b += lp_.Avalue_[k] * col_value[col];
      }

      a = (0.5 / mu) * a;
      b = (0.5 / mu) * b + 0.5 * lp_.colCost_[col];

      double theta = -b / a; // todo: -b / 2a? see notes.
      double delta_x = 0;

      if (theta > 0)
        delta_x = std::min(theta, lp_.colUpper_[col] - col_value[col]);
      else
        delta_x = std::max(theta, lp_.colLower_[col] - col_value[col]);

      col_value[col] += delta_x;

    }

    // For the moment update only after each pass of columns. Later maybe
    // update objective, row_value, residual after each component update.
    update();
  }
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
  for (int iteration = 0; iteration < K; iteration++) {
    // Minimize quadratic function.
    quadratic.minimize_by_component(mu, lambda);

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
                    "Solution set at the end of feasibility search.\n");

  return HighsStatus::OK;
}
