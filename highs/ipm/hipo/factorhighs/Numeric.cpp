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

  return kRetOk;
}

void Numeric::getReg(std::vector<double>& reg) {
  // unpermute regularisation
  permuteVector(total_reg_, S_->iperm());

  reg = std::move(total_reg_);
}

}  // namespace hipo