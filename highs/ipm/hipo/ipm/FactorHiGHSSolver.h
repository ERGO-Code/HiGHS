#ifndef HIPO_FACTORHIGHS_SOLVER_H
#define HIPO_FACTORHIGHS_SOLVER_H

#include <algorithm>

#include "Info.h"
#include "IpmData.h"
#include "LinearSolver.h"
#include "Model.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/factorhighs/FactorHiGHS.h"

namespace hipo {

class FactorHiGHSSolver : public LinearSolver {
  // symbolic factorisation
  Symbolic S_;

  // numeric factorisation
  Numeric N_;

  // object to perform factorisation
  FHsolver FH_;

  // normal equations data
  std::vector<Int> ptrNE_, rowsNE_;
  std::vector<double> valNE_;
  HighsSparseMatrix AT_;

  const Regularisation& regul_;

  Info* info_ = nullptr;
  IpmData* data_ = nullptr;
  const LogHighs& log_;

  Int chooseNla(const Model& model, Options& options);
  Int setNla(const Model& model, Options& options);
  void setParallel(Options& options);
  Int buildNEstructureDense(
      const HighsSparseMatrix& A,
      int64_t max_num_nz = std::numeric_limits<Int>::max());
  Int buildNEstructureSparse(
      const HighsSparseMatrix& A,
      int64_t max_num_nz = std::numeric_limits<Int>::max());
  Int buildNEvalues(const HighsSparseMatrix& A,
                    const std::vector<double>& scaling);

 public:
  FactorHiGHSSolver(const Options& options, const Regularisation& regul, Info* info,
                    IpmData* record, const LogHighs& log);

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
