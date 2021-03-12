/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef MIP_HIGHS_MIP_SOLVER_H_
#define MIP_HIGHS_MIP_SOLVER_H_

#include "Highs.h"
#include "lp_data/HighsOptions.h"

struct HighsMipSolverData;
class HighsCutPool;

class HighsMipSolver {
 public:
  const HighsOptions* options_mip_;
  const HighsLp* model_;
  HighsModelStatus modelstatus_;
  bool submip;
  const HighsBasis* rootbasis;

  std::unique_ptr<HighsMipSolverData> mipdata_;

  void run();

  int numCol() const { return model_->numCol_; }

  int numRow() const { return model_->numRow_; }

  int numNonzero() const { return model_->Aindex_.size(); }

  const double* colCost() const { return model_->colCost_.data(); }

  double colCost(int col) const { return model_->colCost_[col]; }

  const double* rowLower() const { return model_->rowLower_.data(); }

  double rowLower(int col) const { return model_->rowLower_[col]; }

  const double* rowUpper() const { return model_->rowUpper_.data(); }

  double rowUpper(int col) const { return model_->rowUpper_[col]; }

  const HighsVarType* variableType() const {
    return model_->integrality_.data();
  }

  HighsVarType variableType(int col) const { return model_->integrality_[col]; }

  HighsMipSolver(const HighsOptions& options, const HighsLp& lp,
                 bool submip = false);

  ~HighsMipSolver();

  void setModel(const HighsLp& model) { model_ = &model; }

  HighsTimer timer_;
  PresolveComponent presolve_;
  HighsPresolveStatus runPresolve();
  HighsPostsolveStatus runPostsolve();
};

#endif
