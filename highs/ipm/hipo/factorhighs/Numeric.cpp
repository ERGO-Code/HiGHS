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

Int Numeric::solve(std::vector<double>& x) const {
  // Return the number of solves performed

  if (!sn_columns_ || !S_) return kRetInvalidPointer;

  HIPO_CLOCK_CREATE;

  // initialise solve handler
  HybridSolveHandler SH(*S_, *sn_columns_, swaps_, pivot_2x2_, *data_);

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

void Numeric::getReg(std::vector<double>& reg) {
  // unpermute regularisation
  permuteVector(total_reg_, S_->iperm());

  reg = std::move(total_reg_);
}

}  // namespace hipo