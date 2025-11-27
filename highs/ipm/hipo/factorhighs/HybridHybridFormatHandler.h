#ifndef FACTORHIGHS_HYBRID_HYBRID_FORMAT_H
#define FACTORHIGHS_HYBRID_HYBRID_FORMAT_H

#include "DataCollector.h"
#include "FormatHandler.h"

namespace hipo {

class HybridHybridFormatHandler : public FormatHandler {
  std::vector<Int64> diag_start_;
  DataCollector& data_;

  void initFrontal() override;
  void initClique() override;
  void assembleFrontal(Int i, Int j, double val) override;
  void assembleFrontalMultiple(Int num, const double* child, Int nc,
                               Int child_sn, Int row, Int col, Int i,
                               Int j) override;
  Int denseFactorise(double reg_thresh) override;
  void assembleClique(const double* child, Int nc, Int child_sn) override;
  void extremeEntries() override;

 public:
  HybridHybridFormatHandler(const Symbolic& S, Int sn, const Regul& regul,
                            DataCollector& data, std::vector<double>& frontal,
                            double* clique_ptr);
};

}  // namespace hipo

#endif