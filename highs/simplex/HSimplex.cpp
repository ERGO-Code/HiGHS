/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplex.cpp
 * @brief
 */

#include "simplex/HSimplex.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "lp_data/HighsLpUtils.h"
// #include "util/HighsSort.h"

using std::max;
using std::min;
// using std::runtime_error;

void accommodateAlienBasis(HighsLpSolverObject& solver_object) {
  HighsLp& lp = solver_object.lp_;
  HighsBasis& basis = solver_object.basis_;
  HighsOptions& options = solver_object.options_;
  assert(basis.alien);
  HighsInt num_row = lp.num_row_;
  HighsInt num_col = lp.num_col_;
  assert((int)basis.col_status.size() >= num_col);
  assert((int)basis.row_status.size() >= num_row);
  std::vector<HighsInt> basic_index;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (basis.col_status[iCol] == HighsBasisStatus::kBasic)
      basic_index.push_back(iCol);
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (basis.row_status[iRow] == HighsBasisStatus::kBasic)
      basic_index.push_back(num_col + iRow);
  }
  HighsInt num_basic_variables = basic_index.size();
  HFactor factor;
  factor.setupGeneral(&lp.a_matrix_, num_basic_variables, basic_index.data(),
                      kDefaultPivotThreshold, kDefaultPivotTolerance,
                      kHighsDebugLevelMin, &options.log_options);
  HighsInt rank_deficiency = factor.build();
  // Must not have timed out
  assert(rank_deficiency >= 0);
  // Deduce the basis from basic_index
  //
  // Set all basic variables to nonbasic
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (basis.col_status[iCol] == HighsBasisStatus::kBasic)
      basis.col_status[iCol] = HighsBasisStatus::kNonbasic;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (basis.row_status[iRow] == HighsBasisStatus::kBasic)
      basis.row_status[iRow] = HighsBasisStatus::kNonbasic;
  }
  // Set at most the first num_row variables in basic_index to basic
  const HighsInt use_basic_variables = std::min(num_row, num_basic_variables);
  // num_basic_variables is no longer needed, so can be used as a check
  num_basic_variables = 0;
  for (HighsInt iRow = 0; iRow < use_basic_variables; iRow++) {
    HighsInt iVar = basic_index[iRow];
    if (iVar < num_col) {
      basis.col_status[iVar] = HighsBasisStatus::kBasic;
    } else {
      basis.row_status[iVar - num_col] = HighsBasisStatus::kBasic;
    }
    num_basic_variables++;
  }
  // Complete the assignment of basic variables using the logicals of
  // non-pivotal rows
  const HighsInt num_missing = num_row - num_basic_variables;
  for (HighsInt k = 0; k < num_missing; k++) {
    HighsInt iRow = factor.row_with_no_pivot[rank_deficiency + k];
    basis.row_status[iRow] = HighsBasisStatus::kBasic;
    num_basic_variables++;
  }
  assert(num_basic_variables == num_row);
}

HighsStatus formSimplexLpBasisAndFactorReturn(
    const HighsStatus return_status, HighsLpSolverObject& solver_object) {
  HighsLp& lp = solver_object.lp_;
  HighsLp& ekk_lp = solver_object.ekk_instance_.lp_;
  if (lp.is_moved_) lp.moveBackLpAndUnapplyScaling(ekk_lp);
  return return_status;
}

HighsStatus formSimplexLpBasisAndFactor(HighsLpSolverObject& solver_object,
                                        const bool only_from_known_basis) {
  // Ideally, forms a SimplexBasis from the HighsBasis in the
  // HighsLpSolverObject
  //
  // If only_from_known_basis is true and
  // initialiseSimplexLpBasisAndFactor finds that there is no simplex
  // basis, then its error return is passed down
  //
  // If only_from_known_basis is false, then the basis is completed
  // with logicals if it is rank deficient (from singularity or being
  // incomplete)
  //
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsLp& lp = solver_object.lp_;
  HighsBasis& basis = solver_object.basis_;
  HighsOptions& options = solver_object.options_;
  HEkk& ekk_instance = solver_object.ekk_instance_;
  HighsSimplexStatus& ekk_status = ekk_instance.status_;
  lp.ensureColwise();
  const bool passed_scaled = lp.is_scaled_;
  // Consider scaling the LP
  if (!passed_scaled) considerSimplexScaling(options, lp);
  const bool check_basis = basis.alien || (!basis.valid && basis.useful);
  if (check_basis) {
    // The basis needs to be checked for rank deficiency, and possibly
    // completed if it is rectangular
    //
    // If it's not valid but useful, but not alien,
    // accommodateAlienBasis will assert, so make the basis alien
    basis.alien = true;
    assert(!only_from_known_basis);
    accommodateAlienBasis(solver_object);
    basis.alien = false;
    // Unapply any scaling used only for factorization to check and
    // complete the basis
    if (!passed_scaled) lp.unapplyScale();
    // Check that any scaling the LP arrived with has not been removed
    assert(lp.is_scaled_ == passed_scaled);
    return HighsStatus::kOk;
  }
  // Move the HighsLpSolverObject's LP to EKK
  ekk_instance.moveLp(solver_object);
  if (!ekk_status.has_basis) {
    // The Ekk instance has no simplex basis, so pass the HiGHS basis
    HighsStatus call_status = ekk_instance.setBasis(basis);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "setBasis");
    if (return_status == HighsStatus::kError)
      return formSimplexLpBasisAndFactorReturn(return_status, solver_object);
  }
  // Now form the invert
  assert(ekk_status.has_basis);
  call_status =
      ekk_instance.initialiseSimplexLpBasisAndFactor(only_from_known_basis);
  // If the current basis cannot be inverted, return an error
  if (call_status != HighsStatus::kOk)
    return formSimplexLpBasisAndFactorReturn(HighsStatus::kError,
                                             solver_object);
  // Once the invert is formed, move back the LP and remove any scaling.
  return formSimplexLpBasisAndFactorReturn(HighsStatus::kOk, solver_object);
}

void SimplexBasis::clear() {
  this->hash = 0;
  this->basicIndex_.clear();
  this->nonbasicFlag_.clear();
  this->nonbasicMove_.clear();
  this->debug_id = -1;
  this->debug_update_count = -1;
  this->debug_origin_name = "None";
}

void SimplexBasis::setup(const HighsInt num_col, const HighsInt num_row) {
  this->hash = 0;
  this->basicIndex_.resize(num_row);
  this->nonbasicFlag_.resize(num_col + num_row);
  this->nonbasicMove_.resize(num_col + num_row);
  this->debug_id = -1;
  this->debug_update_count = -1;
  this->debug_origin_name = "None";
}

void appendNonbasicColsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                               HighsInt XnumNewCol) {
  assert(highs_basis.valid);
  if (!highs_basis.valid) {
    printf("\n!!Appending columns to invalid basis!!\n\n");
  }
  // Add nonbasic structurals
  if (XnumNewCol == 0) return;
  HighsInt newNumCol = lp.num_col_ + XnumNewCol;
  highs_basis.col_status.resize(newNumCol);
  // Make any new columns nonbasic
  for (HighsInt iCol = lp.num_col_; iCol < newNumCol; iCol++) {
    if (!highs_isInfinity(-lp.col_lower_[iCol])) {
      highs_basis.col_status[iCol] = HighsBasisStatus::kLower;
    } else if (!highs_isInfinity(lp.col_upper_[iCol])) {
      highs_basis.col_status[iCol] = HighsBasisStatus::kUpper;
    } else {
      highs_basis.col_status[iCol] = HighsBasisStatus::kZero;
    }
  }
}

void appendNonbasicColsToBasis(HighsLp& lp, SimplexBasis& basis,
                               HighsInt XnumNewCol) {
  // Add nonbasic structurals
  if (XnumNewCol == 0) return;
  HighsInt newNumCol = lp.num_col_ + XnumNewCol;
  HighsInt newNumTot = newNumCol + lp.num_row_;
  basis.nonbasicFlag_.resize(newNumTot);
  basis.nonbasicMove_.resize(newNumTot);
  // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
  // necessary
  for (HighsInt iRow = lp.num_row_ - 1; iRow >= 0; iRow--) {
    HighsInt iCol = basis.basicIndex_[iRow];
    if (iCol >= lp.num_col_) {
      // This basic variable is a row, so shift its index
      basis.basicIndex_[iRow] += XnumNewCol;
    }
    basis.nonbasicFlag_[newNumCol + iRow] =
        basis.nonbasicFlag_[lp.num_col_ + iRow];
    basis.nonbasicMove_[newNumCol + iRow] =
        basis.nonbasicMove_[lp.num_col_ + iRow];
  }
  // Make any new columns nonbasic
  for (HighsInt iCol = lp.num_col_; iCol < newNumCol; iCol++) {
    basis.nonbasicFlag_[iCol] = kNonbasicFlagTrue;
    double lower = lp.col_lower_[iCol];
    double upper = lp.col_upper_[iCol];
    int8_t move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (fabs(lower) < fabs(upper)) {
          move = kNonbasicMoveUp;
        } else {
          move = kNonbasicMoveDn;
        }
      } else {
        // Lower (since upper bound is infinite)
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = kNonbasicMoveDn;
    } else {
      // FREE
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis.nonbasicMove_[iCol] = move;
  }
}

void appendBasicRowsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                            HighsInt XnumNewRow) {
  assert(highs_basis.valid);
  if (!highs_basis.valid) {
    printf("\n!!Appending columns to invalid basis!!\n\n");
  }
  // Add basic logicals
  if (XnumNewRow == 0) return;
  HighsInt newNumRow = lp.num_row_ + XnumNewRow;
  highs_basis.row_status.resize(newNumRow);
  // Make the new rows basic
  for (HighsInt iRow = lp.num_row_; iRow < newNumRow; iRow++) {
    highs_basis.row_status[iRow] = HighsBasisStatus::kBasic;
  }
}

void appendBasicRowsToBasis(HighsLp& lp, SimplexBasis& basis,
                            HighsInt XnumNewRow) {
  // Add basic logicals
  if (XnumNewRow == 0) return;

  HighsInt newNumRow = lp.num_row_ + XnumNewRow;
  HighsInt newNumTot = lp.num_col_ + newNumRow;
  basis.nonbasicFlag_.resize(newNumTot);
  basis.nonbasicMove_.resize(newNumTot);
  basis.basicIndex_.resize(newNumRow);
  // Make the new rows basic
  for (HighsInt iRow = lp.num_row_; iRow < newNumRow; iRow++) {
    basis.nonbasicFlag_[lp.num_col_ + iRow] = kNonbasicFlagFalse;
    basis.nonbasicMove_[lp.num_col_ + iRow] = 0;
    basis.basicIndex_[iRow] = lp.num_col_ + iRow;
  }
}

void getUnscaledInfeasibilities(const HighsOptions& options,
                                const HighsScale& scale,
                                const SimplexBasis& basis,
                                const HighsSimplexInfo& info,
                                HighsInfo& highs_info) {
  const double primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  const double dual_feasibility_tolerance = options.dual_feasibility_tolerance;

  HighsInt& num_primal_infeasibilities = highs_info.num_primal_infeasibilities;
  double& max_primal_infeasibility = highs_info.max_primal_infeasibility;
  double& sum_primal_infeasibilities = highs_info.sum_primal_infeasibilities;
  HighsInt& num_dual_infeasibilities = highs_info.num_dual_infeasibilities;
  double& max_dual_infeasibility = highs_info.max_dual_infeasibility;
  double& sum_dual_infeasibilities = highs_info.sum_dual_infeasibilities;

  // Zero the counts of unscaled primal and dual infeasibilities
  num_primal_infeasibilities = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibilities = 0;
  num_dual_infeasibilities = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibilities = 0;

  double scale_mu = 1.0;
  assert(int(scale.col.size()) == scale.num_col);
  assert(int(scale.row.size()) == scale.num_row);
  for (HighsInt iVar = 0; iVar < scale.num_col + scale.num_row; iVar++) {
    // Look at the dual infeasibilities of nonbasic variables
    if (basis.nonbasicFlag_[iVar] == kNonbasicFlagFalse) continue;
    // No dual infeasibility for fixed rows and columns
    if (info.workLower_[iVar] == info.workUpper_[iVar]) continue;
    bool col = iVar < scale.num_col;
    HighsInt iCol = 0;
    HighsInt iRow = 0;
    if (col) {
      iCol = iVar;
      assert(int(scale.col.size()) > iCol);
      scale_mu = 1 / (scale.col[iCol] / scale.cost);
    } else {
      iRow = iVar - scale.num_col;
      assert(int(scale.row.size()) > iRow);
      scale_mu = scale.row[iRow] * scale.cost;
    }
    const double dual = info.workDual_[iVar];
    const double lower = info.workLower_[iVar];
    const double upper = info.workUpper_[iVar];
    const double unscaled_dual = dual * scale_mu;

    double dual_infeasibility;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(unscaled_dual);
    } else {
      // Not fixed: any dual infeasibility is given by value signed by
      // nonbasicMove. This assumes that nonbasicMove=0 for fixed
      // variables
      dual_infeasibility = -basis.nonbasicMove_[iVar] * unscaled_dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
  // Look at the primal infeasibilities of basic variables
  for (HighsInt ix = 0; ix < scale.num_row; ix++) {
    HighsInt iVar = basis.basicIndex_[ix];
    bool col = iVar < scale.num_col;
    HighsInt iCol = 0;
    HighsInt iRow = 0;
    if (col) {
      iCol = iVar;
      scale_mu = scale.col[iCol];
    } else {
      iRow = iVar - scale.num_col;
      scale_mu = 1 / scale.row[iRow];
    }
    double unscaled_lower = info.baseLower_[ix] * scale_mu;
    double unscaled_value = info.baseValue_[ix] * scale_mu;
    double unscaled_upper = info.baseUpper_[ix] * scale_mu;
    // @primal_infeasibility calculation
    double primal_infeasibility = 0;
    if (unscaled_value < unscaled_lower - primal_feasibility_tolerance) {
      primal_infeasibility = unscaled_lower - unscaled_value;
    } else if (unscaled_value > unscaled_upper + primal_feasibility_tolerance) {
      primal_infeasibility = unscaled_value - unscaled_upper;
    }
    if (primal_infeasibility > 0) {
      num_primal_infeasibilities++;
      max_primal_infeasibility =
          max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibilities += primal_infeasibility;
    }
  }
  setSolutionStatus(highs_info);
}

void setSolutionStatus(HighsInfo& highs_info) {
  if (highs_info.num_primal_infeasibilities < 0) {
    highs_info.primal_solution_status = kSolutionStatusNone;
  } else if (highs_info.num_primal_infeasibilities == 0) {
    highs_info.primal_solution_status = kSolutionStatusFeasible;
  } else {
    highs_info.primal_solution_status = kSolutionStatusInfeasible;
  }
  if (highs_info.num_dual_infeasibilities < 0) {
    highs_info.dual_solution_status = kSolutionStatusNone;
  } else if (highs_info.num_dual_infeasibilities == 0) {
    highs_info.dual_solution_status = kSolutionStatusFeasible;
  } else {
    highs_info.dual_solution_status = kSolutionStatusInfeasible;
  }
}
// SCALING

bool considerSimplexScaling(const HighsOptions& options, HighsLp& lp) {
  // Indicate whether new scaling has been determined in the return value.
  bool new_scaling = false;
  // Consider scaling the LP - either by finding new factors or by
  // applying any existing factors
  const bool allow_scaling =
      lp.num_col_ > 0 &&
      options.simplex_scale_strategy != kSimplexScaleStrategyOff;
  if (lp.scale_.has_scaling && !allow_scaling) {
    // LP had scaling before, but now it is not permitted, clear any
    // scaling. Return true as scaling position has changed
    lp.clearScale();
    return true;
  }
  const bool scaling_not_tried = lp.scale_.strategy == kSimplexScaleStrategyOff;
  const bool new_scaling_strategy =
      options.simplex_scale_strategy != lp.scale_.strategy &&
      options.simplex_scale_strategy != kSimplexScaleStrategyChoose;
  const bool try_scaling =
      allow_scaling && (scaling_not_tried || new_scaling_strategy);
  if (try_scaling) {
    // Scaling will be tried, so ensure that any previous scaling is not applied
    lp.unapplyScale();
    const bool analyse_lp_data =
        kHighsAnalysisLevelModelData & options.highs_analysis_level;
    if (analyse_lp_data) analyseLp(options.log_options, lp);
    simplexScaleLp(options, lp);
    // If the LP is now scaled, then the scaling is new
    new_scaling = lp.is_scaled_;
    if (analyse_lp_data && lp.is_scaled_) analyseLp(options.log_options, lp);
  } else if (lp.scale_.has_scaling) {
    // Scaling factors are known, so ensure that they are applied
    lp.applyScale();
  }
  // Ensure that either the LP has scale factors and is scaled, or
  // it doesn't have scale factors and isn't scaled
  assert(lp.scale_.has_scaling == lp.is_scaled_);
  return new_scaling;
}

void simplexScaleLp(const HighsOptions& options, HighsLp& lp,
                    const bool force_scaling) {
  lp.clearScaling();
  HighsInt numCol = lp.num_col_;
  HighsInt numRow = lp.num_row_;
  // Scaling not well defined for models with no columns
  assert(numCol > 0);
  vector<double>& colCost = lp.col_cost_;
  vector<double>& colLower = lp.col_lower_;
  vector<double>& colUpper = lp.col_upper_;
  vector<double>& rowLower = lp.row_lower_;
  vector<double>& rowUpper = lp.row_upper_;

  // Save the simplex_scale_strategy so that the option can be
  // modified for the course of this method
  HighsInt simplex_scale_strategy = options.simplex_scale_strategy;
  // Determine the actual strategy to use
  HighsInt use_scale_strategy = simplex_scale_strategy;
  if (use_scale_strategy == kSimplexScaleStrategyChoose) {
    // HiGHS is left to choose: currently use forced equilibration, but maybe do
    // something more intelligent
    use_scale_strategy = kSimplexScaleStrategyForcedEquilibration;
  }
  // Find out range of matrix values and skip matrix scaling if all
  // |values| are in [0.2, 5]
  const double no_scaling_original_matrix_min_value = 0.2;
  const double no_scaling_original_matrix_max_value = 5.0;
  double original_matrix_min_value = kHighsInf;
  double original_matrix_max_value = 0;
  lp.a_matrix_.range(original_matrix_min_value, original_matrix_max_value);
  if (kSimplexScaleDev) {
    double min_cost = kHighsInf;
    double max_cost = -kHighsInf;
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      double abs_cost = std::fabs(lp.col_cost_[iCol]);
      if (abs_cost > 0) {
	min_cost = std::min(abs_cost, min_cost);
	max_cost = std::max(abs_cost, max_cost);
      }
    }
    if (kSimplexScaleDevReport) {
      printf("grepSimplexRangeTxt HighsScale: costs in [%g, %g]; matrix in [%g, %g]; %s: %s\n",
	     min_cost, max_cost,
	     original_matrix_min_value, original_matrix_max_value,
	     lp.model_name_.c_str(),
	     lp.origin_name_.c_str());
      printf("grepSimplexRangeCsv,%g,%g,%g,%g,%s,%s\n",
	     min_cost, max_cost,
	     original_matrix_min_value, original_matrix_max_value,
	     lp.model_name_.c_str(),
	     lp.origin_name_.c_str());
    }
  }
  // Possibly force scaling, otherwise base the decision on the range
  // of values in the matrix, values that will be used later for
  // reporting
  const bool no_scaling = force_scaling
                              ? false
                              : (original_matrix_min_value >=
                                 no_scaling_original_matrix_min_value) &&
                                    (original_matrix_max_value <=
                                     no_scaling_original_matrix_max_value);
  HighsScale& scale = lp.scale_;
  bool scaled_matrix = false;
  if (no_scaling) {
    // No matrix scaling
    highsLogDev(options.log_options, HighsLogType::kInfo,
		"Scaling: Matrix has [min, max] values of [%g, %g] within "
		"[%g, %g] so no scaling performed\n",
		original_matrix_min_value, original_matrix_max_value,
		no_scaling_original_matrix_min_value,
		no_scaling_original_matrix_max_value);
  } else {
    // Try scaling, so assign unit factors - partly because initial
    // factors may be assumed by the scaling method, but also because
    // scaling factors may not be computed for empty rows/columns
    scale.col.assign(numCol, 1);
    scale.row.assign(numRow, 1);
    const bool equilibration_scaling =
        use_scale_strategy == kSimplexScaleStrategyEquilibration ||
        use_scale_strategy == kSimplexScaleStrategyForcedEquilibration;
    // Try scaling. Value of scaled_matrix indicates whether scaling
    // was considered valuable (and performed). If it's not valuable
    // then the matrix remains unscaled
    if (equilibration_scaling) {
      scaled_matrix = equilibrationScaleMatrix(options, lp, use_scale_strategy);
    } else {
      scaled_matrix = maxValueScaleMatrix(options, lp, use_scale_strategy);
    }
    if (scaled_matrix) {
      // Matrix is scaled, so scale the bounds and costs
      for (HighsInt iCol = 0; iCol < numCol; iCol++) {
        colLower[iCol] /= scale.col[iCol];
        colUpper[iCol] /= scale.col[iCol];
        colCost[iCol] *= scale.col[iCol];
      }
      for (HighsInt iRow = 0; iRow < numRow; iRow++) {
        rowLower[iRow] *= scale.row[iRow];
        rowUpper[iRow] *= scale.row[iRow];
      }
      scale.has_scaling = true;
      scale.num_col = numCol;
      scale.num_row = numRow;
      scale.cost = 1.0;
      lp.is_scaled_ = true;
    }
  }
  bool scaled_costs = false;
  if (kSimplexScaleDev &&
      use_scale_strategy == kSimplexScaleStrategyMaxValueMatrixAndCost) {
    // Consider cost scaling
    simplexScaleCost(options, lp);
    scaled_costs = scale.cost != 1.0;
    scale.has_scaling = true;
    lp.is_scaled_ = true;
  }
  // If neither costs nor matrix are scaled, then clear the scaling
  if (!scaled_costs && !scaled_matrix) lp.clearScaling();
  // Record the scaling strategy used
  lp.scale_.strategy = use_scale_strategy;
}

bool equilibrationScaleMatrix(const HighsOptions& options, HighsLp& lp,
                              const HighsInt use_scale_strategy) {
  HighsInt numCol = lp.num_col_;
  HighsInt numRow = lp.num_row_;
  HighsScale& scale = lp.scale_;
  vector<double>& colScale = scale.col;
  vector<double>& rowScale = scale.row;
  vector<HighsInt>& Astart = lp.a_matrix_.start_;
  vector<HighsInt>& Aindex = lp.a_matrix_.index_;
  vector<double>& Avalue = lp.a_matrix_.value_;
  vector<double>& colCost = lp.col_cost_;

  HighsInt simplex_scale_strategy = use_scale_strategy;

  double original_matrix_min_value = kHighsInf;
  double original_matrix_max_value = 0;
  for (HighsInt k = 0, AnX = Astart[numCol]; k < AnX; k++) {
    double value = fabs(Avalue[k]);
    original_matrix_min_value = min(original_matrix_min_value, value);
    original_matrix_max_value = max(original_matrix_max_value, value);
  }

  // Include cost in scaling if minimum nonzero cost is less than 0.1
  double min_nonzero_cost = kHighsInf;
  for (HighsInt i = 0; i < numCol; i++) {
    if (colCost[i]) min_nonzero_cost = min(fabs(colCost[i]), min_nonzero_cost);
  }
  bool include_cost_in_scaling = false;
  include_cost_in_scaling = min_nonzero_cost < 0.1;

  // Limits on scaling factors
  double max_allow_scale;
  double min_allow_scale;
  // Now that kHighsInf =
  // std::numeric_limits<double>::infinity(), this Qi-trick doesn't
  // work so, in recognition, use the old value of kHighsInf
  const double finite_infinity = 1e200;
  max_allow_scale = pow(2.0, options.allowed_matrix_scale_factor);
  min_allow_scale = 1 / max_allow_scale;

  double min_allow_col_scale = min_allow_scale;
  double max_allow_col_scale = max_allow_scale;
  double min_allow_row_scale = min_allow_scale;
  double max_allow_row_scale = max_allow_scale;

  // Search up to 6 times
  vector<double> row_min_value(numRow, finite_infinity);
  vector<double> row_max_value(numRow, 1 / finite_infinity);
  for (HighsInt search_count = 0; search_count < 6; search_count++) {
    // Find column scale, prepare row data
    for (HighsInt iCol = 0; iCol < numCol; iCol++) {
      // For column scale (find)
      double col_min_value = finite_infinity;
      double col_max_value = 1 / finite_infinity;
      double abs_col_cost = fabs(colCost[iCol]);
      if (include_cost_in_scaling && abs_col_cost != 0) {
        col_min_value = min(col_min_value, abs_col_cost);
        col_max_value = max(col_max_value, abs_col_cost);
      }
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        double value = fabs(Avalue[k]) * rowScale[Aindex[k]];
        col_min_value = min(col_min_value, value);
        col_max_value = max(col_max_value, value);
      }
      double col_equilibration = 1 / sqrt(col_min_value * col_max_value);
      // Ensure that column scale factor is not excessively large or small
      colScale[iCol] =
          min(max(min_allow_col_scale, col_equilibration), max_allow_col_scale);
      // For row scale (only collect)
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        double value = fabs(Avalue[k]) * colScale[iCol];
        row_min_value[iRow] = min(row_min_value[iRow], value);
        row_max_value[iRow] = max(row_max_value[iRow], value);
      }
    }
    // For row scale (find)
    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
      double row_equilibration =
          1 / sqrt(row_min_value[iRow] * row_max_value[iRow]);
      // Ensure that row scale factor is not excessively large or small
      rowScale[iRow] =
          min(max(min_allow_row_scale, row_equilibration), max_allow_row_scale);
    }
    row_min_value.assign(numRow, finite_infinity);
    row_max_value.assign(numRow, 1 / finite_infinity);
  }
  // Make it numerically better
  // Also determine the max and min row and column scaling factors
  double min_col_scale = finite_infinity;
  double max_col_scale = 1 / finite_infinity;
  double min_row_scale = finite_infinity;
  double max_row_scale = 1 / finite_infinity;
  const double log2 = log(2.0);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / log2 + 0.5));
    min_col_scale = min(colScale[iCol], min_col_scale);
    max_col_scale = max(colScale[iCol], max_col_scale);
  }
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    rowScale[iRow] = pow(2.0, floor(log(rowScale[iRow]) / log2 + 0.5));
    min_row_scale = min(rowScale[iRow], min_row_scale);
    max_row_scale = max(rowScale[iRow], max_row_scale);
  }
  // Apply scaling to matrix and bounds
  double matrix_min_value = finite_infinity;
  double matrix_max_value = 0;
  double min_original_col_equilibration = finite_infinity;
  double sum_original_log_col_equilibration = 0;
  double max_original_col_equilibration = 0;
  double min_original_row_equilibration = finite_infinity;
  double sum_original_log_row_equilibration = 0;
  double max_original_row_equilibration = 0;
  double min_col_equilibration = finite_infinity;
  double sum_log_col_equilibration = 0;
  double max_col_equilibration = 0;
  double min_row_equilibration = finite_infinity;
  double sum_log_row_equilibration = 0;
  double max_row_equilibration = 0;
  vector<double> original_row_min_value(numRow, finite_infinity);
  vector<double> original_row_max_value(numRow, 1 / finite_infinity);
  row_min_value.assign(numRow, finite_infinity);
  row_max_value.assign(numRow, 1 / finite_infinity);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    double original_col_min_value = finite_infinity;
    double original_col_max_value = 1 / finite_infinity;
    double col_min_value = finite_infinity;
    double col_max_value = 1 / finite_infinity;
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      HighsInt iRow = Aindex[k];
      const double original_value = fabs(Avalue[k]);
      original_col_min_value = min(original_value, original_col_min_value);
      original_col_max_value = max(original_value, original_col_max_value);
      original_row_min_value[iRow] =
          min(original_row_min_value[iRow], original_value);
      original_row_max_value[iRow] =
          max(original_row_max_value[iRow], original_value);
      Avalue[k] *= (colScale[iCol] * rowScale[iRow]);
      const double value = fabs(Avalue[k]);
      col_min_value = min(value, col_min_value);
      col_max_value = max(value, col_max_value);
      row_min_value[iRow] = min(row_min_value[iRow], value);
      row_max_value[iRow] = max(row_max_value[iRow], value);
    }
    matrix_min_value = min(matrix_min_value, col_min_value);
    matrix_max_value = max(matrix_max_value, col_max_value);

    const double original_col_equilibration =
        1 / sqrt(original_col_min_value * original_col_max_value);
    min_original_col_equilibration =
        min(original_col_equilibration, min_original_col_equilibration);
    sum_original_log_col_equilibration += log(original_col_equilibration);
    max_original_col_equilibration =
        max(original_col_equilibration, max_original_col_equilibration);
    const double col_equilibration = 1 / sqrt(col_min_value * col_max_value);
    min_col_equilibration = min(col_equilibration, min_col_equilibration);
    sum_log_col_equilibration += log(col_equilibration);
    max_col_equilibration = max(col_equilibration, max_col_equilibration);
  }

  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    const double original_row_equilibration =
        1 / sqrt(original_row_min_value[iRow] * original_row_max_value[iRow]);
    min_original_row_equilibration =
        min(original_row_equilibration, min_original_row_equilibration);
    sum_original_log_row_equilibration += log(original_row_equilibration);
    max_original_row_equilibration =
        max(original_row_equilibration, max_original_row_equilibration);
    const double row_equilibration =
        1 / sqrt(row_min_value[iRow] * row_max_value[iRow]);
    min_row_equilibration = min(row_equilibration, min_row_equilibration);
    sum_log_row_equilibration += log(row_equilibration);
    max_row_equilibration = max(row_equilibration, max_row_equilibration);
  }
  const double geomean_original_col_equilibration =
      exp(sum_original_log_col_equilibration / numCol);
  const double geomean_original_row_equilibration =
      exp(sum_original_log_row_equilibration / numRow);
  const double geomean_col_equilibration =
      exp(sum_log_col_equilibration / numCol);
  const double geomean_row_equilibration =
      exp(sum_log_row_equilibration / numRow);
  if (options.log_dev_level) {
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Original equilibration: min/mean/max %11.4g/%11.4g/%11.4g "
        "(cols); min/mean/max %11.4g/%11.4g/%11.4g (rows)\n",
        min_original_col_equilibration, geomean_original_col_equilibration,
        max_original_col_equilibration, min_original_row_equilibration,
        geomean_original_row_equilibration, max_original_row_equilibration);
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Final    equilibration: min/mean/max %11.4g/%11.4g/%11.4g "
        "(cols); min/mean/max %11.4g/%11.4g/%11.4g (rows)\n",
        min_col_equilibration, geomean_col_equilibration, max_col_equilibration,
        min_row_equilibration, geomean_row_equilibration,
        max_row_equilibration);
  }

  // Compute the mean equilibration improvement
  const double geomean_original_col =
      max(geomean_original_col_equilibration,
          1 / geomean_original_col_equilibration);
  const double geomean_original_row =
      max(geomean_original_row_equilibration,
          1 / geomean_original_row_equilibration);
  const double geomean_col =
      max(geomean_col_equilibration, 1 / geomean_col_equilibration);
  const double geomean_row =
      max(geomean_row_equilibration, 1 / geomean_row_equilibration);
  const double mean_equilibration_improvement =
      sqrt((geomean_original_col * geomean_original_row) /
           (geomean_col * geomean_row));
  // Compute the extreme equilibration improvement
  const double original_col_ratio =
      max_original_col_equilibration / min_original_col_equilibration;
  const double original_row_ratio =
      max_original_row_equilibration / min_original_row_equilibration;
  const double col_ratio = max_col_equilibration / min_col_equilibration;
  const double row_ratio = max_row_equilibration / min_row_equilibration;
  const double extreme_equilibration_improvement =
      (original_col_ratio + original_row_ratio) / (col_ratio + row_ratio);
  // Compute the max/min matrix value improvement
  const double matrix_value_ratio = matrix_max_value / matrix_min_value;
  const double original_matrix_value_ratio =
      original_matrix_max_value / original_matrix_min_value;
  const double matrix_value_ratio_improvement =
      original_matrix_value_ratio / matrix_value_ratio;
  if (options.log_dev_level) {
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Extreme equilibration improvement =      ( %11.4g + "
        "%11.4g) / ( %11.4g + %11.4g)  =      %11.4g / %11.4g  = %11.4g\n",
        original_col_ratio, original_row_ratio, col_ratio, row_ratio,
        (original_col_ratio + original_row_ratio), (col_ratio + row_ratio),
        extreme_equilibration_improvement);
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Mean    equilibration improvement = sqrt(( %11.4g * "
        "%11.4g) / ( %11.4g * %11.4g)) = sqrt(%11.4g / %11.4g) = %11.4g\n",
        geomean_original_col, geomean_original_row, geomean_col, geomean_row,
        (geomean_original_col * geomean_original_row),
        (geomean_col * geomean_row), mean_equilibration_improvement);
    highsLogDev(
        options.log_options, HighsLogType::kInfo,
        "Scaling: Yields [min, max, ratio] matrix values of [%0.4g, %0.4g, "
        "%0.4g]; Originally [%0.4g, %0.4g, %0.4g]: Improvement of %0.4g\n",
        matrix_min_value, matrix_max_value, matrix_value_ratio,
        original_matrix_min_value, original_matrix_max_value,
        original_matrix_value_ratio, matrix_value_ratio_improvement);
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "Scaling: Improves    mean equilibration by a factor %0.4g\n",
                mean_equilibration_improvement);
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "Scaling: Improves extreme equilibration by a factor %0.4g\n",
                extreme_equilibration_improvement);
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "Scaling: Improves max/min matrix values by a factor %0.4g\n",
                matrix_value_ratio_improvement);
  }
  const bool possibly_abandon_scaling =
      simplex_scale_strategy != kSimplexScaleStrategyForcedEquilibration;
  const double improvement_factor = extreme_equilibration_improvement *
                                    mean_equilibration_improvement *
                                    matrix_value_ratio_improvement;

  const double improvement_factor_required = 1.0;
  const bool poor_improvement =
      improvement_factor < improvement_factor_required;

  // Possibly abandon scaling if it's not improved equilibration significantly
  if (possibly_abandon_scaling && poor_improvement) {
    // Unscale the matrix
    for (HighsInt iCol = 0; iCol < numCol; iCol++) {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        Avalue[k] /= (colScale[iCol] * rowScale[iRow]);
      }
    }
    if (options.log_dev_level)
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Improvement factor %0.4g < %0.4g required, so no "
                  "scaling applied\n",
                  improvement_factor, improvement_factor_required);
    return false;
  } else {
    if (options.log_dev_level) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Factors are in [%0.4g, %0.4g] for columns and in "
                  "[%0.4g, %0.4g] for rows\n",
                  min_col_scale, max_col_scale, min_row_scale, max_row_scale);
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Improvement factor is %0.4g >= %0.4g so scale LP\n",
                  improvement_factor, improvement_factor_required);
      if (extreme_equilibration_improvement < 1.0) {
        highsLogDev(
            options.log_options, HighsLogType::kWarning,
            "Scaling: Applying scaling with extreme improvement of %0.4g\n",
            extreme_equilibration_improvement);
      }
      if (mean_equilibration_improvement < 1.0) {
        highsLogDev(
            options.log_options, HighsLogType::kWarning,
            "Scaling: Applying scaling with mean improvement of %0.4g\n",
            mean_equilibration_improvement);
      }
      if (matrix_value_ratio_improvement < 1.0) {
        highsLogDev(options.log_options, HighsLogType::kWarning,
                    "Scaling: Applying scaling with matrix value ratio "
                    "improvement of %0.4g\n",
                    matrix_value_ratio_improvement);
      }
      if (improvement_factor < 10 * improvement_factor_required) {
        highsLogDev(options.log_options, HighsLogType::kWarning,
                    "Scaling: Applying scaling with improvement factor %0.4g "
                    "< 10*(%0.4g) improvement\n",
                    improvement_factor, improvement_factor_required);
      }
    }
  }
  return true;
}

bool maxValueScaleMatrix(const HighsOptions& options, HighsLp& lp,
                         const HighsInt use_scale_strategy) {
  HighsInt numCol = lp.num_col_;
  HighsInt numRow = lp.num_row_;
  HighsScale& scale = lp.scale_;
  vector<double>& colScale = scale.col;
  vector<double>& rowScale = scale.row;
  vector<HighsInt>& Astart = lp.a_matrix_.start_;
  vector<HighsInt>& Aindex = lp.a_matrix_.index_;
  vector<double>& Avalue = lp.a_matrix_.value_;

  assert(options.simplex_scale_strategy >= kSimplexScaleStrategyMaxValue);

  const double log2 = log(2.0);
  const double max_allow_scale = pow(2.0, options.allowed_matrix_scale_factor);
  const double min_allow_scale = 1 / max_allow_scale;

  const double min_allow_col_scale = min_allow_scale;
  const double max_allow_col_scale = max_allow_scale;
  const double min_allow_row_scale = min_allow_scale;
  const double max_allow_row_scale = max_allow_scale;

  double min_row_scale = kHighsInf;
  double max_row_scale = 0;
  double original_matrix_min_value = kHighsInf;
  double original_matrix_max_value = 0;
  // Determine the row scaling. Also determine the max/min row scaling
  // factors, and max/min original matrix values
  vector<double> row_max_value(numRow, 0);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      const HighsInt iRow = Aindex[k];
      const double value = fabs(Avalue[k]);
      row_max_value[iRow] = max(row_max_value[iRow], value);
      original_matrix_min_value = min(original_matrix_min_value, value);
      original_matrix_max_value = max(original_matrix_max_value, value);
    }
  }
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    if (row_max_value[iRow]) {
      double row_scale_value = 1 / row_max_value[iRow];
      // Convert the row scale factor to the nearest power of two, and
      // ensure that it is not excessively large or small
      row_scale_value = pow(2.0, floor(log(row_scale_value) / log2 + 0.5));
      row_scale_value =
          min(max(min_allow_row_scale, row_scale_value), max_allow_row_scale);
      min_row_scale = min(row_scale_value, min_row_scale);
      max_row_scale = max(row_scale_value, max_row_scale);
      rowScale[iRow] = row_scale_value;
    }
  }
  // Determine the column scaling, whilst applying the row scaling
  // Also determine the max/min column scaling factors, and max/min
  // matrix values
  double min_col_scale = kHighsInf;
  double max_col_scale = 0;
  double matrix_min_value = kHighsInf;
  double matrix_max_value = 0;
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    double col_max_value = 0;
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      const HighsInt iRow = Aindex[k];
      Avalue[k] *= rowScale[iRow];
      const double value = fabs(Avalue[k]);
      col_max_value = max(col_max_value, value);
    }
    if (col_max_value) {
      double col_scale_value = 1 / col_max_value;
      // Convert the col scale factor to the nearest power of two, and
      // ensure that it is not excessively large or small
      col_scale_value = pow(2.0, floor(log(col_scale_value) / log2 + 0.5));
      col_scale_value =
          min(max(min_allow_col_scale, col_scale_value), max_allow_col_scale);
      min_col_scale = min(col_scale_value, min_col_scale);
      max_col_scale = max(col_scale_value, max_col_scale);
      colScale[iCol] = col_scale_value;
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        Avalue[k] *= colScale[iCol];
        const double value = fabs(Avalue[k]);
        matrix_min_value = min(matrix_min_value, value);
        matrix_max_value = max(matrix_max_value, value);
      }
    }
  }
  const double matrix_value_ratio = matrix_max_value / matrix_min_value;
  const double original_matrix_value_ratio =
      original_matrix_max_value / original_matrix_min_value;
  const double matrix_value_ratio_improvement =
      original_matrix_value_ratio / matrix_value_ratio;

  const double improvement_factor = matrix_value_ratio_improvement;

  const double improvement_factor_required = 1.0;
  const bool poor_improvement =
      improvement_factor <= improvement_factor_required;

  if (poor_improvement) {
    // Unscale the matrix
    for (HighsInt iCol = 0; iCol < numCol; iCol++) {
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
        HighsInt iRow = Aindex[k];
        Avalue[k] /= (colScale[iCol] * rowScale[iRow]);
      }
    }
    if (options.log_dev_level)
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Improvement factor %0.4g < %0.4g required, so no "
                  "scaling applied\n",
                  improvement_factor, improvement_factor_required);
    return false;
  } else {
    if (options.log_dev_level) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Scaling: Factors are in [%0.4g, %0.4g] for columns and in "
                  "[%0.4g, %0.4g] for rows\n",
                  min_col_scale, max_col_scale, min_row_scale, max_row_scale);
      highsLogDev(
          options.log_options, HighsLogType::kInfo,
          "Scaling: Yields [min, max, ratio] matrix values of [%0.4g, %0.4g, "
          "%0.4g]; Originally [%0.4g, %0.4g, %0.4g]: Improvement of %0.4g\n",
          matrix_min_value, matrix_max_value, matrix_value_ratio,
          original_matrix_min_value, original_matrix_max_value,
          original_matrix_value_ratio, matrix_value_ratio_improvement);
    }
    return true;
  }
}

HighsStatus applyScalingToLpCol(HighsLp& lp, const HighsInt col,
                                const double colScale) {
  if (col < 0) return HighsStatus::kError;
  if (col >= lp.num_col_) return HighsStatus::kError;
  if (!colScale) return HighsStatus::kError;

  lp.a_matrix_.scaleCol(col, colScale);
  lp.col_cost_[col] *= colScale;
  if (colScale > 0) {
    lp.col_lower_[col] /= colScale;
    lp.col_upper_[col] /= colScale;
  } else {
    const double new_upper = lp.col_lower_[col] / colScale;
    lp.col_lower_[col] = lp.col_upper_[col] / colScale;
    lp.col_upper_[col] = new_upper;
  }
  return HighsStatus::kOk;
}

HighsStatus applyScalingToLpRow(HighsLp& lp, const HighsInt row,
                                const double rowScale) {
  if (row < 0) return HighsStatus::kError;
  if (row >= lp.num_row_) return HighsStatus::kError;
  if (!rowScale) return HighsStatus::kError;

  lp.a_matrix_.scaleRow(row, rowScale);
  if (rowScale > 0) {
    lp.row_lower_[row] *= rowScale;
    lp.row_upper_[row] *= rowScale;
  } else {
    const double new_upper = lp.row_lower_[row] * rowScale;
    lp.row_lower_[row] = lp.row_upper_[row] * rowScale;
    lp.row_upper_[row] = new_upper;
  }
  return HighsStatus::kOk;
}

void simplexUnscaleSolution(HighsSolution& solution, const HighsScale& scale) {
  assert(scale.has_scaling);
  assert(solution.col_value.size() == static_cast<size_t>(scale.num_col));
  assert(solution.col_dual.size() == static_cast<size_t>(scale.num_col));
  assert(solution.row_value.size() == static_cast<size_t>(scale.num_row));
  assert(solution.row_dual.size() == static_cast<size_t>(scale.num_row));

  for (HighsInt iCol = 0; iCol < scale.num_col; iCol++) {
    solution.col_value[iCol] *= scale.col[iCol];
    solution.col_dual[iCol] /= (scale.col[iCol] / scale.cost);
  }
  for (HighsInt iRow = 0; iRow < scale.num_row; iRow++) {
    solution.row_value[iRow] /= scale.row[iRow];
    solution.row_dual[iRow] *= (scale.row[iRow] * scale.cost);
  }
}

void simplexScaleCost(const HighsOptions& options, HighsLp& lp) {
  double& cost_scale = lp.scale_.cost;
  // Scale the costs by no less than minAlwCostScale
  double max_allowed_cost_scale = pow(2.0, options.allowed_cost_scale_factor);
  double max_nonzero_cost = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    if (lp.col_cost_[iCol]) 
      max_nonzero_cost = max(fabs(lp.col_cost_[iCol]), max_nonzero_cost);
  // Scaling the costs up effectively means that the problem is solved
  // to a smaller dual tolerance, which seems pointless. HiGHS will
  // warn about excessively small costs, and users can set
  // user_objective_scale to a positive value.
  //
  // Scaling the costs down effectively means that the problem is
  // solved to a larger dual tolerance, which may require clean-up
  // iterations when the scaling is removed after solving the scaled
  // problem
  cost_scale = 1.0;
  // Scale if the max cost is greater than 16
  if (max_nonzero_cost > 16) {
    cost_scale = max_nonzero_cost;
    cost_scale = pow(2.0, floor(log(cost_scale) / log(2.0) + 0.5));
    cost_scale = min(cost_scale, max_allowed_cost_scale);
  }
  if (cost_scale == 1.0) {
    highsLogDev(options.log_options, HighsLogType::kInfo,
		"LP cost vector not scaled down: max cost is %g\n",
		max_nonzero_cost);
  } else {
    // Scale the costs (and record of max_nonzero_cost) by cost_scale, being at
    // most max_allowed_cost_scale
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
      lp.col_cost_[iCol] /= cost_scale;
    highsLogDev(options.log_options, HighsLogType::kInfo,
		"LP cost vector scaled down by %g: max scaled cost is %g\n", cost_scale,
		max_nonzero_cost/cost_scale);
  }
}

void simplexUnscaleCost(HighsLp& lp) {
  double& cost_scale = lp.scale_.cost;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    lp.col_cost_[iCol] *= cost_scale;
}

bool isBasisRightSize(const HighsLp& lp, const SimplexBasis& basis) {
  bool right_size = true;
  right_size =
      (HighsInt)basis.nonbasicFlag_.size() == lp.num_col_ + lp.num_row_ &&
      right_size;
  right_size =
      (HighsInt)basis.nonbasicMove_.size() == lp.num_col_ + lp.num_row_ &&
      right_size;
  right_size = (HighsInt)basis.basicIndex_.size() == lp.num_row_ && right_size;
  return right_size;
}
