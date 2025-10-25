#include "Numeric.h"

#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "HybridSolveHandler.h"
#include "ReturnValues.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"
#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"

namespace hipo {

Int Numeric::solve(std::vector<double>& x, Int* solve_count,
                   double* omega) const {
  // Return the number of solves performed

  if (!sn_columns_ || !S_) return kRetInvalidPointer;

  // initialise solve handler
  SH_.reset(new HybridSolveHandler(*S_, *sn_columns_, swaps_, pivot_2x2_));

  SH_->setData(data_);

#if HIPO_TIMING_LEVEL >= 1
  Clock clock{};
#endif

#if HIPO_TIMING_LEVEL >= 2
  Clock clock_fine{};
#endif

  // permute rhs
  permuteVectorInverse(x, S_->iperm());

  // make a copy of permuted rhs, for refinement
  const std::vector<double> rhs(x);

#if HIPO_TIMING_LEVEL >= 2
  if (data_) data_->sumTime(kTimeSolvePrepare, clock_fine.stop());
  clock_fine.start();
#endif

  // solve
  SH_->forwardSolve(x);
  SH_->diagSolve(x);
  SH_->backwardSolve(x);

#if HIPO_TIMING_LEVEL >= 2
  if (data_) data_->sumTime(kTimeSolveSolve, clock_fine.stop());
#endif

  // iterative refinement
  auto refine_data = refine(rhs, x);

#if HIPO_TIMING_LEVEL >= 2
  clock_fine.start();
#endif

  // unpermute solution
  permuteVector(x, S_->iperm());

#if HIPO_TIMING_LEVEL >= 2
  if (data_) data_->sumTime(kTimeSolvePrepare, clock_fine.stop());
#endif

#if HIPO_TIMING_LEVEL >= 1
  if (data_) data_->sumTime(kTimeSolve, clock.stop());
#endif

  if (solve_count) *solve_count = refine_data.first + 1;
  if (omega) *omega = refine_data.second;

  return kRetOk;
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

std::pair<Int, double> Numeric::refine(const std::vector<double>& rhs,
                                       std::vector<double>& x) const {
  // Return the number of solver performed

  double old_omega{};
  Int solves_counter{};

#if HIPO_TIMING_LEVEL >= 2
  Clock clock{};
#endif

  // compute residual
  std::vector<double> res = residualQuad(rhs, x);

#if HIPO_TIMING_LEVEL >= 2
  if (data_) data_->sumTime(kTimeSolveResidual, clock.stop());
  clock.start();
#endif

  double omega = computeOmega(rhs, x, res);

#if HIPO_TIMING_LEVEL >= 2
  if (data_) data_->sumTime(kTimeSolveOmega, clock.stop());
#endif

  // if(log_) log_->printDevVerbose("   # start  %.2e\n", omega);

  Int iter = 0;
  for (; iter < kMaxRefinementIter; ++iter) {
    // termination criterion
    if (omega < kRefinementTolerance) break;

#if HIPO_TIMING_LEVEL >= 2
    clock.start();
#endif

    // compute correction
    SH_->forwardSolve(res);
    SH_->diagSolve(res);
    SH_->backwardSolve(res);
    ++solves_counter;

#if HIPO_TIMING_LEVEL >= 2
    if (data_) data_->sumTime(kTimeSolveSolve, clock.stop());
    clock.start();
#endif

    // add correction
    std::vector<double> temp(x);
    vectorAdd(temp, res);

#if HIPO_TIMING_LEVEL >= 2
    if (data_) data_->sumTime(kTimeSolvePrepare, clock.stop());
    clock.start();
#endif

    // compute new residual
    res = residualQuad(rhs, temp);

#if HIPO_TIMING_LEVEL >= 2
    if (data_) data_->sumTime(kTimeSolveResidual, clock.stop());
    clock.start();
#endif

    old_omega = omega;
    omega = computeOmega(rhs, temp, res);

    // if(log_) log_->printDevVerbose("   # refine %.2e\n", omega);

#if HIPO_TIMING_LEVEL >= 2
    if (data_) data_->sumTime(kTimeSolveOmega, clock.stop());
#endif

    if (omega < old_omega) {
      x = temp;
    } else {
      omega = old_omega;
      // if(log_) log_->printDevVerbose("   ## reject\n");
      break;
    }
  }

  return {solves_counter, omega};
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

}  // namespace hipo