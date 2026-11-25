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
  HybridHybridFormatHandler(const Symbolic& S, Int sn, DataCollector& data,
                            std::vector<double>& frontal, double* clique_ptr,
                            const FHoptions& FH_opt, double* total_reg,
                            double* total_unreg_pivots, std::vector<Int>& swaps,
                            std::vector<double>& pivot_2x2);
};

}  // namespace hipo

#endif