#include "FormatHandler.h"

#include "CallAndTimeBlas.h"
#include "DenseFact.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

FormatHandler::FormatHandler(const Symbolic& S, Int sn, const Regul& regul,
                             std::vector<double>& frontal,
                             std::vector<double>& clique,
                             std::vector<Int>& swaps,
                             std::vector<double>& pivot_2x2)
    : S_{&S},
      regul_{regul},
      sn_{sn},
      nb_{S_->blockSize()},
      sn_size_{S_->snStart(sn_ + 1) - S_->snStart(sn_)},
      ldf_{S_->ptr(sn_ + 1) - S_->ptr(sn_)},
      ldc_{ldf_ - sn_size_},
      frontal_{frontal},
      clique_{clique},
      swaps_{swaps},
      pivot_2x2_{pivot_2x2} {
  local_reg_.assign(sn_size_, 0.0);
  swaps_.assign(sn_size_, 0);
  pivot_2x2_.assign(sn_size_, 0.0);
}

void FormatHandler::terminate(std::vector<double>& total_reg) {
  // Move local regularisation into total regularisation.
  for (Int i = 0; i < sn_size_; ++i)
    total_reg[S_->snStart(sn_) + i] = local_reg_[i];

  // This function should not require a lock, since all threads access different
  // locations within total_reg, sn_columns and schur_contribution.
}

}  // namespace hipo
