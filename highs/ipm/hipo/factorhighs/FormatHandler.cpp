#include "FormatHandler.h"

#include "CallAndTimeBlas.h"
#include "DenseFact.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

FormatHandler::FormatHandler(const Symbolic& S, Int sn,
                             std::vector<double>& frontal, double* clique_ptr,
                             const FHoptions& FH_opt, double* total_reg,
                             double* total_unreg_pivots, Int* total_swaps,
                             double* total_pivot_2x2)
    : S_{&S},
      sn_{sn},
      sn_size_{S_->snStart(sn_ + 1) - S_->snStart(sn_)},
      ldf_{(Int)(S_->ptr(sn_ + 1) - S_->ptr(sn_))},
      ldc_{ldf_ - sn_size_},
      FH_opt_{FH_opt},
      frontal_{frontal},
      clique_ptr_{clique_ptr},
      reg_{total_reg + S_->snStart(sn_)},
      unreg_pivots_{total_unreg_pivots + S_->snStart(sn_)},
      swaps_{total_swaps + S_->snStart(sn_)},
      pivot_2x2_{total_pivot_2x2 + S_->snStart(sn_)} {}

void FormatHandler::terminate(std::vector<double>& clique) {
  // Move the local clique into its final position.
  // In this way, the shared object schur_contribution_ is accessed only here,
  // while a local copy is used for the assembly and dense factorisation. This
  // should avoid the problem of false sharing.
  // This function should not require a lock, since all threads access different
  // locations within schur_contribution.

  clique = std::move(clique_);
}

}  // namespace hipo
