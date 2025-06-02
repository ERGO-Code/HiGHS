#include "Numeric.h"

#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "HybridSolveHandler.h"
#include "Timing.h"
#include "ipm/hpm/auxiliary/Auxiliary.h"
#include "ipm/hpm/auxiliary/HpmLog.h"
#include "ipm/hpm/auxiliary/VectorOperations.h"
#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"

namespace highspm {

Numeric::Numeric(const Symbolic& S) : S_{S} {
  // initialise solve handler
  SH_.reset(new HybridSolveHandler(S_, sn_columns_, swaps_, pivot_2x2_));
}

Int Numeric::solve(std::vector<double>& x) const {
  // Return the number of solves performed

#if HPM_TIMING_LEVEL >= 1
  Clock clock{};
#endif

#if HPM_TIMING_LEVEL >= 2
  Clock clock_fine{};
#endif

  // permute rhs
  permuteVectorInverse(x, S_.iperm());

  // make a copy of permuted rhs, for refinement
  const std::vector<double> rhs(x);

#if HPM_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeSolvePrepare, clock_fine.stop());
  clock_fine.start();
#endif

  // solve
  SH_->forwardSolve(x);
  SH_->diagSolve(x);
  SH_->backwardSolve(x);
  DataCollector::get()->countSolves();

#if HPM_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeSolveSolve, clock_fine.stop());
#endif

  // iterative refinement
  Int refine_solves = refine(rhs, x);

#if HPM_TIMING_LEVEL >= 2
  clock_fine.start();
#endif

  // unpermute solution
  permuteVector(x, S_.iperm());

#if HPM_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeSolvePrepare, clock_fine.stop());
#endif

#if HPM_TIMING_LEVEL >= 1
  DataCollector::get()->sumTime(kTimeSolve, clock.stop());
#endif

  return refine_solves + 1;
}

std::vector<double> Numeric::residual(const std::vector<double>& rhs,
                                      const std::vector<double>& x) const {
  // Compute the residual rhs - A * x - Reg * x
  std::vector<double> res(rhs);
  symProduct(ptrA_, rowsA_, valA_, x, res, -1.0);
  for (Int i = 0; i < x.size(); ++i) res[i] -= total_reg_[i] * x[i];

  return res;
}

std::vector<double> Numeric::residualQuad(const std::vector<double>& rhs,
                                          const std::vector<double>& x) const {
  std::vector<HighsCDouble> res(rhs.size());
  for (Int i = 0; i < rhs.size(); ++i) res[i] = rhs[i];

  symProductQuad(ptrA_, rowsA_, valA_, x, res, -1.0);

  for (Int i = 0; i < x.size(); ++i)
    res[i] -= (HighsCDouble)total_reg_[i] * (HighsCDouble)x[i];

  std::vector<double> result(res.size());
  for (Int i = 0; i < res.size(); ++i) {
    result[i] = (double)res[i];
  }

  return result;
}

Int Numeric::refine(const std::vector<double>& rhs,
                    std::vector<double>& x) const {
  // Return the number of solver performed

  double old_omega{};
  Int solves_counter{};

#if HPM_TIMING_LEVEL >= 2
  Clock clock{};
#endif

  // compute residual
  std::vector<double> res = residualQuad(rhs, x);

#if HPM_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeSolveResidual, clock.stop());
  clock.start();
#endif

  double omega = computeOmega(rhs, x, res);

#if HPM_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeSolveOmega, clock.stop());
#endif

  Log::printDevVerbose("   # start  %.2e\n", omega);

  Int iter = 0;
  for (; iter < kMaxRefinementIter; ++iter) {
    // termination criterion
    if (omega < kRefinementTolerance) break;

#if HPM_TIMING_LEVEL >= 2
    clock.start();
#endif

    // compute correction
    SH_->forwardSolve(res);
    SH_->diagSolve(res);
    SH_->backwardSolve(res);
    DataCollector::get()->countSolves();
    ++solves_counter;

#if HPM_TIMING_LEVEL >= 2
    DataCollector::get()->sumTime(kTimeSolveSolve, clock.stop());
    clock.start();
#endif

    // add correction
    std::vector<double> temp(x);
    vectorAdd(temp, res);

#if HPM_TIMING_LEVEL >= 2
    DataCollector::get()->sumTime(kTimeSolvePrepare, clock.stop());
    clock.start();
#endif

    // compute new residual
    res = residualQuad(rhs, temp);

#if HPM_TIMING_LEVEL >= 2
    DataCollector::get()->sumTime(kTimeSolveResidual, clock.stop());
    clock.start();
#endif

    old_omega = omega;
    omega = computeOmega(rhs, temp, res);

    Log::printDevVerbose("   # refine %.2e\n", omega);

#if HPM_TIMING_LEVEL >= 2
    DataCollector::get()->sumTime(kTimeSolveOmega, clock.stop());
#endif

    if (omega < old_omega) {
      x = temp;
    } else {
      omega = old_omega;
      Log::printDevVerbose("   ## reject\n");
      break;
    }
  }

  DataCollector::get()->setOmega(omega);

  return solves_counter;
}

double Numeric::computeOmega(const std::vector<double>& b,
                             const std::vector<double>& x,
                             const std::vector<double>& res) const {
  // Termination of iterative refinement based on "Solving sparse linear systems
  // with sparse backward error", Arioli, Demmel, Duff.

  const Int n = x.size();

  // infinity norm of x
  const double inf_norm_x = infNorm(x);

  // |A|*|x|
  std::vector<double> abs_prod(n);
  for (Int col = 0; col < n; ++col) {
    for (Int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
      Int row = rowsA_[el];
      double val = valA_[el];
      abs_prod[row] += std::abs(val) * std::abs(x[col]);
      if (row != col) abs_prod[col] += std::abs(val) * std::abs(x[row]);
    }
  }

  double omega_1{};
  double omega_2{};

  for (Int i = 0; i < n; ++i) {
    // threshold 1000 * n * eps * (||Ai|| * ||x|| + |bi|)
    double tau =
        1000.0 * n * 1e-16 * (inf_norm_cols_[i] * inf_norm_x + std::abs(b[i]));

    if (abs_prod[i] + std::abs(b[i]) > tau) {
      // case 1, denominator is large enough
      double omega = std::abs(res[i]) / (abs_prod[i] + std::abs(b[i]));
      omega_1 = std::max(omega_1, omega);
    } else {
      // case 2, denominator would be small, change it
      double omega =
          std::abs(res[i]) / (abs_prod[i] + one_norm_cols_[i] * inf_norm_x);
      omega_2 = std::max(omega_2, omega);
    }
  }

  return omega_1 + omega_2;
}

void Numeric::conditionNumber() const {
  HighsRandom random;

  const Int n = S_.size();

  // estimate largest eigenvalue with power iteration:
  // x <- x / ||x||
  // y <- M * x
  // lambda = ||y||

  double lambda_large{};
  std::vector<double> x(n);
  for (Int i = 0; i < n; ++i) x[i] = random.fraction();

  for (Int iter = 0; iter < 10; ++iter) {
    // normalize x
    double x_norm = norm2(x);
    vectorScale(x, 1.0 / x_norm);

    // multiply by matrix
    std::vector<double> y(n);
    symProduct(ptrA_, rowsA_, valA_, x, y);

    double norm_y = norm2(y);

    if (std::abs(lambda_large - norm_y) / std::max(1.0, lambda_large) < 1e-2) {
      // converged
      break;
    }

    lambda_large = norm_y;
    x = std::move(y);
  }

  // estimate inverse of smallest eigenvalue with inverse power iteration:
  // x <- x / ||x||
  // y <- M^-1 * x
  // lambda_inv = ||y||

  double lambda_small_inv{};
  for (Int i = 0; i < n; ++i) x[i] = random.fraction();

  for (Int iter = 0; iter < 10; ++iter) {
    // normalize x
    double x_norm = norm2(x);
    vectorScale(x, 1.0 / x_norm);

    // solve with matrix
    std::vector<double> y(x);
    SH_->forwardSolve(y);
    SH_->diagSolve(y);
    SH_->backwardSolve(y);

    double norm_y = norm2(y);

    if (std::abs(lambda_small_inv - norm_y) / std::max(1.0, lambda_small_inv) <
        1e-2) {
      // converged
      break;
    }

    lambda_small_inv = norm_y;
    x = std::move(y);
  }

  // condition number: lambda_large / lambda_small
  double cond = lambda_large * lambda_small_inv;
}

}  // namespace highspm