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
  // object to perform factorisation
  FHsolver FH_;

  // symbolic factorisation
  Symbolic S_;

  // normal equations data
  std::vector<Int> ptrNE_;
  std::vector<Int> rowsNE_;
  std::vector<double> valNE_;
  std::vector<Int> ptrA_rw_, idxA_rw_;
  std::vector<Int> corr_A_;

  // augmented system data
  std::vector<Int> ptrAS_;
  std::vector<Int> rowsAS_;
  std::vector<double> valAS_;

  const Regularisation& regul_;

  Info& info_;
  IpmData& data_;
  const LogHighs& log_;

  const Model& model_;
  const HighsSparseMatrix& A_;
  const Int mA_, nA_, nzA_;

  Options& options_;

  Int chooseNla();
  Int setNla();
  void setParallel();
  Int chooseOrdering(const std::vector<Int>& rows, const std::vector<Int>& ptr,
                      const std::vector<Int>& signs, Symbolic& S);

  Int buildNEstructure(Int64 nz_limit = kHighsIInf);
  Int buildNEvalues(const std::vector<double>& scaling);
  void freeNEmemory();

  Int buildASstructure(Int64 nz_limit = kHighsIInf);
  Int buildASvalues(const std::vector<double>& scaling);
  void freeASmemory();

  Int analyseAS(Symbolic& S);
  Int analyseNE(Symbolic& S, Int64 nz_limit = kHighsIInf);

 public:
  FactorHiGHSSolver(Options& options, const Model& model,
                    const Regularisation& regul, Info& info, IpmData& record,
                    const LogHighs& log);

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
  Int setup() override;
  void clear() override;
  double flops() const override;
  double spops() const override;
  double nz() const override;
  void getReg(std::vector<double>& reg) override;
};

}  // namespace hipo

#endif
