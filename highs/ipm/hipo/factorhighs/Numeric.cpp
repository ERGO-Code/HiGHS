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

  if (!sn_columns_ || !S_) return kRetInvalidPointer;

  HIPO_CLOCK_CREATE;

  // initialise solve handler
  HybridSolveHandler SH(*S_, *sn_columns_, swaps_, pivot_2x2_, *data_,
                        pivoting_);

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

void Numeric::getReg(double* reg) {
  std::memcpy(reg, total_reg_.data(), total_reg_.size() * sizeof(double));
}

void Numeric::inertia(Int& pos, Int& neg, Int& zero, double tol) const {
  HybridSolveHandler SH(*S_, *sn_columns_, swaps_, pivot_2x2_, *data_,
                        pivoting_);
  SH.inertia(pos, neg, zero, tol);
}

}  // namespace hipo