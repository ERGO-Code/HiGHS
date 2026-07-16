#include "Numeric.h"

#include <algorithm>

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

void Numeric::finaliseFactor() {
  // Compute the per-block swap flags and size the solve workspace.

  const Int nb = S_->blockSize();
  Int max_ld = 0;

  swap_flags_.assign(S_->sn(), {});
  for (Int sn = 0; sn < S_->sn(); ++sn) {
    const Int ldSn = S_->ptr(sn + 1) - S_->ptr(sn);
    max_ld = std::max(max_ld, ldSn);

    const Int sn_size = S_->snStart(sn + 1) - S_->snStart(sn);
    const Int n_blocks = (sn_size - 1) / nb + 1;
    std::vector<uint8_t>& flags = swap_flags_[sn];
    flags.assign(n_blocks, 0);

    // defensive: without pivoting there may be no swaps stored
    if ((Int)swaps_[sn].size() < sn_size) continue;

    for (Int j = 0; j < n_blocks; ++j) {
      const Int jb = std::min(nb, sn_size - nb * j);
      const Int* sw = &swaps_[sn][nb * j];
      for (Int i = 0; i < jb; ++i) {
        if (sw[i] != i) {
          flags[j] = 1;
          break;
        }
      }
    }
  }

  solve_work_.resize(max_ld);
}

Int Numeric::solve(std::vector<double>& x) const {
  // Return the number of solves performed

  if (!sn_columns_ || !S_) return kRetInvalidPointer;

  HIPO_CLOCK_CREATE;

  // initialise solve handler
  HybridSolveHandler SH(*S_, *sn_columns_, swaps_, swap_flags_, pivot_2x2_,
                        solve_work_, *data_);

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