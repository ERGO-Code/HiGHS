#include "Numeric.h"

#include "DataCollector.h"
#include "FactorHighsSettings.h"
#include "HybridSolveHandler.h"
#include "ReturnValues.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"
#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"

namespace hipo {

Int Numeric::solve(double* x) const {
  // Return the number of solves performed

  if (!sn_columns_ || !S_ || !data_ || !options_) return kRetInvalidPointer;

  HIPO_CLOCK_CREATE;

  // initialise solve handler
  HybridSolveHandler SH(*S_, *sn_columns_, swaps_, pivot_2x2_, *data_,
                        *options_);

  // permute rhs
  HIPO_CLOCK_START(2);
  permuteVectorInverse(x, S_->iperm());
  HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare);

  // solve
  HIPO_CLOCK_START(2);
  SH.forwardSolve(x);
  SH.diagSolve(x);
  SH.backwardSolve(x);
  HIPO_CLOCK_STOP(2, *data_, kTimeSolveSolve);

  // unpermute solution
  HIPO_CLOCK_START(2);
  permuteVector(x, S_->iperm());
  HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare);

  HIPO_CLOCK_STOP(1, *data_, kTimeSolve);
  return kRetOk;
}

Int Numeric::forwardSolve(double* x) const {
  if (!sn_columns_ || !S_ || !data_ || !options_) return kRetInvalidPointer;
  HybridSolveHandler SH(*S_, *sn_columns_, swaps_, pivot_2x2_, *data_,
                        *options_);
  permuteVectorInverse(x, S_->iperm());
  SH.forwardSolve(x);
  return kRetOk;
}
Int Numeric::diagSolve(double* x) const {
  if (!sn_columns_ || !S_ || !data_ || !options_) return kRetInvalidPointer;
  HybridSolveHandler SH(*S_, *sn_columns_, swaps_, pivot_2x2_, *data_,
                        *options_);
  SH.diagSolve(x);
  return kRetOk;
}
Int Numeric::backwardSolve(double* x) const {
  if (!sn_columns_ || !S_ || !data_ || !options_) return kRetInvalidPointer;
  HybridSolveHandler SH(*S_, *sn_columns_, swaps_, pivot_2x2_, *data_,
                        *options_);
  SH.backwardSolve(x);
  permuteVector(x, S_->iperm());
  return kRetOk;
}

#define SOLVE_MULTIPLE(f)                                        \
  if (k == 1)                                                    \
    return f(x);                                                 \
  else {                                                         \
    highs::parallel::TaskGroup tg;                               \
    const Int n = S_->size();                                    \
    std::atomic<bool> fail{false};                               \
    for (Int i = 0; i < k; ++i) {                                \
      tg.spawn([=, &fail]() {                                    \
        Int status = f(&x[i * n]);                               \
        if (status) fail.store(true, std::memory_order_relaxed); \
      });                                                        \
    }                                                            \
    tg.taskWait();                                               \
    return fail;                                                 \
  }

Int Numeric::solve(double* x, Int k) const { SOLVE_MULTIPLE(solve); }
Int Numeric::forwardSolve(double* x, Int k) const {
  SOLVE_MULTIPLE(forwardSolve);
}
Int Numeric::diagSolve(double* x, Int k) const { SOLVE_MULTIPLE(diagSolve); }
Int Numeric::backwardSolve(double* x, Int k) const {
  SOLVE_MULTIPLE(backwardSolve);
}

void Numeric::getReg(double* reg) {
  std::memcpy(reg, total_reg_.data(), total_reg_.size() * sizeof(double));
}

void Numeric::inertia(Int& pos, Int& neg, Int& zero, double tol) const {
  if (!sn_columns_ || !S_ || !data_ || !options_) return;
  HybridSolveHandler SH(*S_, *sn_columns_, swaps_, pivot_2x2_, *data_,
                        *options_);
  SH.inertia(pos, neg, zero, tol);
}

}  // namespace hipo