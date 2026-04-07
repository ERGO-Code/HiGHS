/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsReformulation.cpp
 * @brief
 */
// #include <sstream>

#include "Highs.h"
// #include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
// #include "mip/HighsMipSolver.h"  // For getGapString
// #include "model/HighsHessianUtils.h"
// #include "simplex/HSimplex.h"
// #include "util/HighsMatrixUtils.h"
// #include "util/HighsSort.h"

HighsStatus Highs::doReformulation() {
  HighsStatus status = HighsStatus::kOk;
  assert(model_.lp_.has_infinite_cost_ ==
         model_.lp_.hasInfiniteCost(options_.infinite_cost));
  if (model_.lp_.has_infinite_cost_) {
    // If the model has infinite costs, then try to remove them. The
    // status indicates the success of this operation and, if
    // it's unsuccessful, the model will not have been modified.
    status = handleInfiniteCost();
    if (status != HighsStatus::kOk) {
      assert(status == HighsStatus::kError);
      setHighsModelStatusAndClearSolutionAndBasis(HighsModelStatus::kUnknown);
      return status;
    }
    assert(!model_.lp_.has_infinite_cost_);
  } else {
    assert(!model_.lp_.hasInfiniteCost(options_.infinite_cost));
  }
  return status;
}

void Highs::undoReformulation(HighsStatus& optimize_status) {
  this->restoreInfiniteCost(optimize_status);
  //  this->model_.lp_.unapplyMods();
}

HighsStatus Highs::handleInfiniteCost() {
  HighsLp& lp = this->model_.lp_;
  if (!lp.has_infinite_cost_) return HighsStatus::kOk;
  HighsLpMods& mods = lp.mods_;
  double inf_cost = this->options_.infinite_cost;
  for (HighsInt k = 0; k < 2; k++) {
    // Pass twice: first checking that infinite costs can be handled,
    // then handling them, so that model is unmodified if infinite
    // costs cannot be handled
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      double cost = lp.col_cost_[iCol];
      if (cost > -inf_cost && cost < inf_cost) continue;
      double lower = lp.col_lower_[iCol];
      double upper = lp.col_upper_[iCol];
      if (lp.isMip()) {
        if (lp.integrality_[iCol] == HighsVarType::kInteger) {
          lower = std::ceil(lower);
          upper = std::floor(upper);
        }
      }
      if (cost <= -inf_cost) {
        if (lp.sense_ == ObjSense::kMinimize) {
          // Minimizing with -inf cost so try to fix at upper bound
          if (upper < kHighsInf) {
            if (k) lp.col_lower_[iCol] = upper;
          } else {
            highsLogUser(options_.log_options, HighsLogType::kError,
                         "Cannot minimize with a cost on variable %d of %g and "
                         "upper bound of %g\n",
                         int(iCol), cost, upper);
            return HighsStatus::kError;
          }
        } else {
          // Maximizing with -inf cost so try to fix at lower bound
          if (lower > -kHighsInf) {
            if (k) lp.col_upper_[iCol] = lower;
          } else {
            highsLogUser(options_.log_options, HighsLogType::kError,
                         "Cannot maximize with a cost on variable %d of %g and "
                         "lower bound of %g\n",
                         int(iCol), cost, lower);
            return HighsStatus::kError;
          }
        }
      } else {
        if (lp.sense_ == ObjSense::kMinimize) {
          // Minimizing with inf cost so try to fix at lower bound
          if (lower > -kHighsInf) {
            if (k) lp.col_upper_[iCol] = lower;
          } else {
            highsLogUser(options_.log_options, HighsLogType::kError,
                         "Cannot minimize with a cost on variable %d of %g and "
                         "lower bound of %g\n",
                         int(iCol), cost, lower);
            return HighsStatus::kError;
          }
        } else {
          // Maximizing with inf cost so try to fix at upper bound
          if (upper < kHighsInf) {
            if (k) lp.col_lower_[iCol] = upper;
          } else {
            highsLogUser(options_.log_options, HighsLogType::kError,
                         "Cannot maximize with a cost on variable %d of %g and "
                         "upper bound of %g\n",
                         int(iCol), cost, upper);
            return HighsStatus::kError;
          }
        }
      }
      if (k) {
        mods.save_inf_cost_variable_index.push_back(iCol);
        mods.save_inf_cost_variable_cost.push_back(cost);
        mods.save_inf_cost_variable_lower.push_back(lower);
        mods.save_inf_cost_variable_upper.push_back(upper);
        lp.col_cost_[iCol] = 0;
      }
    }
  }
  // Infinite costs have been removed, but their presence in the
  // original model is known from mods.save_inf_cost_variable_*, so
  // set lp.has_infinite_cost_ to be false to avoid assert when run()
  // is called using copy of model in MIP solver (See #1446)
  lp.has_infinite_cost_ = false;

  return HighsStatus::kOk;
}

void Highs::restoreInfiniteCost(HighsStatus& optimize_status) {
  HighsLp& lp = this->model_.lp_;
  HighsBasis& basis = this->basis_;
  HighsLpMods& mods = lp.mods_;
  HighsInt num_inf_cost = mods.save_inf_cost_variable_index.size();
  if (num_inf_cost <= 0) return;
  assert(num_inf_cost);
  for (HighsInt ix = 0; ix < num_inf_cost; ix++) {
    HighsInt iCol = mods.save_inf_cost_variable_index[ix];
    double cost = mods.save_inf_cost_variable_cost[ix];
    double lower = mods.save_inf_cost_variable_lower[ix];
    double upper = mods.save_inf_cost_variable_upper[ix];
    double value = solution_.value_valid ? solution_.col_value[iCol] : 0;
    if (basis.valid) {
      assert(basis.col_status[iCol] != HighsBasisStatus::kBasic);
      if (lp.col_lower_[iCol] == lower) {
        basis.col_status[iCol] = HighsBasisStatus::kLower;
      } else {
        basis.col_status[iCol] = HighsBasisStatus::kUpper;
      }
    }
    assert(lp.col_cost_[iCol] == 0);
    if (value) this->info_.objective_function_value += value * cost;
    lp.col_cost_[iCol] = cost;
    lp.col_lower_[iCol] = lower;
    lp.col_upper_[iCol] = upper;
  }
  // Infinite costs have been reintroduced, so reset to true the flag
  // that was set false in Highs::handleInfiniteCost() (See #1446)
  lp.has_infinite_cost_ = true;
  mods.clearInfiniteCostRecord();

  if (this->model_status_ == HighsModelStatus::kInfeasible) {
    // Model is infeasible with the infinite cost variables fixed at
    // appropriate values, so model status cannot be determined
    this->model_status_ = HighsModelStatus::kUnknown;
    setHighsModelStatusAndClearSolutionAndBasis(this->model_status_);
    optimize_status = highsStatusFromHighsModelStatus(model_status_);
  }
}
