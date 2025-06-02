#ifndef HIGHSPM_FACTORHIGHS_SOLVER_H
#define HIGHSPM_FACTORHIGHS_SOLVER_H

#include <algorithm>

#include "HpmInfo.h"
#include "HpmModel.h"
#include "LinearSolver.h"
#include "ipm/hpm/auxiliary/IntConfig.h"
#include "ipm/hpm/factorhighs/FactorHiGHS.h"

namespace highspm {

class FactorHiGHSSolver : public LinearSolver {
  // symbolic factorisation
  Symbolic S_;

  // numeric factorisation
  Numeric N_;

  // keep track of whether as or ne is being factorised
  bool use_as_ = true;

  HpmInfo* info_ = nullptr;

  Int choose(const HpmModel& model, HpmOptions& options);
  Int setNla(const HpmModel& model, HpmOptions& options);
  void setParallel(HpmOptions& options);

 public:
  FactorHiGHSSolver(const HpmOptions& options, HpmInfo* info);

  // Override functions
  Int factorAS(const HighsSparseMatrix& A,
               const std::vector<double>& scaling) override;
  Int factorNE(const HighsSparseMatrix& A,
               const std::vector<double>& scaling) override;
  Int solveNE(const std::vector<double>& rhs,
              std::vector<double>& lhs) override;
  Int solveAS(const std::vector<double>& rhs_x,
              const std::vector<double>& rhs_y, std::vector<double>& lhs_x,
              std::vector<double>& lhs_y) override;
  Int setup(const HpmModel& model, HpmOptions& options) override;
  void clear() override;
  double flops() const override;
  double spops() const override;
  double nz() const override;
};

}  // namespace highspm

#endif
