#ifndef HIPO_FACTORHIGHS_SOLVER_H
#define HIPO_FACTORHIGHS_SOLVER_H

#include <algorithm>

#include "HpmInfo.h"
#include "HpmModel.h"
#include "LinearSolver.h"
#include "ipm/hpm/auxiliary/IntConfig.h"
#include "ipm/hpm/factorhighs/FactorHiGHS.h"

namespace hipo {

class FactorHiGHSSolver : public LinearSolver {
  // symbolic factorisation
  Symbolic S_;

  // numeric factorisation
  Numeric N_;

  // keep track of whether as or ne is being factorised
  bool use_as_ = true;

  Info* info_ = nullptr;

  Int choose(const Model& model, Options& options);
  Int setNla(const Model& model, Options& options);
  void setParallel(Options& options);

 public:
  FactorHiGHSSolver(const Options& options, Info* info);

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
  Int setup(const Model& model, Options& options) override;
  void clear() override;
  double flops() const override;
  double spops() const override;
  double nz() const override;
};

}  // namespace hipo

#endif
