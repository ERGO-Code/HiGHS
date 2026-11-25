#ifndef FACTORHIGHS_HYBRID_HYBRID_FORMAT_H
#define FACTORHIGHS_HYBRID_HYBRID_FORMAT_H

#include "DataCollector.h"
#include "FormatHandler.h"

namespace hipo {

class HybridHybridFormatHandler : public FormatHandler {
  DataCollector& data_;
  const Int nb_;

  void initFrontal() override;
  void initClique() override;
  void assembleFrontal(Int i, Int j, double val) override;
  void assembleChild(Int child_sn, const double* child) override;
  Int denseFactorise(double reg_thresh, const double* static_reg,
                     ScratchArena& scratch) override;

  void assembleFrontalMultiple(Int& num, const double* child, Int nc,
                               Int child_sn, Int row, Int col, Int i, Int j);
  void assembleClique(const double* child, Int nc, Int child_sn);

 public:
  // number of entries of the frontal matrix of a supernode with sn_size
  // columns and ldf rows, stored in blocks of nb columns
  static Int64 frontalSize(Int ldf, Int sn_size, Int nb);

  HybridHybridFormatHandler(const Symbolic& S, Int sn, DataCollector& data,
                            double* frontal, double* clique_ptr,
                            const FHoptions& FH_opt, double* total_reg,
                            double* total_unreg_pivots, Int* total_swaps,
                            double* total_pivot_2x2);
};

}  // namespace hipo

#endif