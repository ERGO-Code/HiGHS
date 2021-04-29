/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplexDebug.cpp
 * @brief
 */

#include "simplex/HSimplexDebug.h"

#include <string>

#include "lp_data/HighsDebug.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolutionDebug.h"
#include "simplex/HEkkDebug.h"
#include "simplex/HEkkDualRow.h"
#include "simplex/HFactorDebug.h"
#include "simplex/HSimplex.h"
#include "simplex/SimplexTimer.h"

// Methods for Ekk

HighsDebugStatus ekkDebugSimplexLp(const HighsModelObject& highs_model_object) {
  // Non-trivially expensive check that the .lp, if valid is .lp scaled
  // according to .scale
  const HEkk& ekk_instance = highs_model_object.ekk_instance_;
  const HighsSimplexStatus& status = ekk_instance.status_;
  if (!status.valid ||
      highs_model_object.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsOptions& options = ekk_instance.options_;
  const HighsLp& highs_lp = highs_model_object.lp_;
  const HighsLp& simplex_lp = ekk_instance.lp_;
  const HighsScale& scale = highs_model_object.scale_;
  const HFactor& factor = ekk_instance.factor_;

  bool right_size = true;
  right_size = (HighsInt)scale.col.size() == simplex_lp.numCol_ && right_size;
  right_size = (HighsInt)scale.row.size() == simplex_lp.numRow_ && right_size;
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "scale size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  // Take a copy of the original LP
  HighsLp check_lp = highs_lp;
  if (applyScalingToLp(options.log_options, check_lp, scale) !=
      HighsStatus::kOk) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "debugSimplexLp: Error scaling check LP\n");
    return HighsDebugStatus::kLogicalError;
  }
  const bool lp_data_ok = check_lp == simplex_lp;
  if (!lp_data_ok) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "debugSimplexLp: Check LP and simplex LP not equal\n");
    assert(lp_data_ok);
    return_status = HighsDebugStatus::kLogicalError;
  }

  if (status.has_basis) {
    const bool basis_correct =
        debugDebugToHighsStatus(ekkDebugBasisCorrect(ekk_instance)) !=
        HighsStatus::kError;
    if (!basis_correct) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Supposed to be a Simplex basis, but incorrect\n");
      assert(basis_correct);
      return_status = HighsDebugStatus::kLogicalError;
    }
  }

  if (status.has_invert) {
    const bool invert_ok = debugDebugToHighsStatus(debugCheckInvert(
                               options, factor)) != HighsStatus::kError;
    if (!invert_ok) {
      highsLogUser(
          options.log_options, HighsLogType::kError,
          "Supposed to be a Simplex basis inverse, but too inaccurate\n");
      assert(invert_ok);
      return_status = HighsDebugStatus::kLogicalError;
    }
  }
  return return_status;
}

HighsDebugStatus debugBasisConsistent(const HighsOptions& options,
                                      const HighsLp& lp,
                                      const SimplexBasis& basis) {
  // Cheap analysis of a Simplex basis, checking vector sizes, numbers
  // of basic/nonbasic variables and non-repetition of basic variables
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  // Check consistency of nonbasicFlag
  if (debugNonbasicFlagConsistent(options, lp, basis) ==
      HighsDebugStatus::kLogicalError) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicFlag inconsistent\n");
    return_status = HighsDebugStatus::kLogicalError;
  }
  const bool right_size = (HighsInt)basis.basicIndex_.size() == lp.numRow_;
  // Check consistency of basicIndex
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "basicIndex size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  // Use localNonbasicFlag so that duplicate entries in basicIndex can
  // be spotted
  vector<int8_t> localNonbasicFlag = basis.nonbasicFlag_;
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    HighsInt iCol = basis.basicIndex_[iRow];
    HighsInt flag = localNonbasicFlag[iCol];
    // Indicate that this column has been found in basicIndex
    localNonbasicFlag[iCol] = -1;
    if (flag) {
      // Nonzero value for localNonbasicFlag entry means that column is either
      if (flag == kNonbasicFlagTrue) {
        // Nonbasic...
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Entry basicIndex_[%" HIGHSINT_FORMAT
                     "] = %" HIGHSINT_FORMAT " is not basic\n",
                     iRow, iCol);
      } else {
        // .. or is -1 since it has already been found in basicIndex
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Entry basicIndex_[%" HIGHSINT_FORMAT
                     "] = %" HIGHSINT_FORMAT " is already basic\n",
                     iRow, iCol);
        assert(flag == -1);
      }
      assert(!flag);
      return_status = HighsDebugStatus::kLogicalError;
    }
  }
  return return_status;
}

void debugDualChuzcFailNorms(
    const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    double& workDataNorm, const HighsInt numVar, const double* workDual,
    double& workDualNorm) {
  workDataNorm = 0;
  for (HighsInt i = 0; i < workCount; i++) {
    double value = workData[i].second;
    workDataNorm += value * value;
  }
  workDataNorm = sqrt(workDataNorm);
  workDualNorm = 0;
  for (HighsInt iVar = 0; iVar < numVar; iVar++) {
    double value = workDual[iVar];
    workDualNorm += value * value;
  }
  workDualNorm = sqrt(workDualNorm);
}

HighsDebugStatus debugDualChuzcFailQuad0(
    const HighsOptions& options, const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    const HighsInt numVar, const double* workDual, const double selectTheta,
    const double remainTheta, const bool force) {
  // Non-trivially expensive assessment of CHUZC failure
  if (options.highs_debug_level < kHighsDebugLevelCostly && !force)
    return HighsDebugStatus::kNotChecked;

  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     No change in loop 2 so return error\n");
  double workDataNorm;
  double workDualNorm;
  debugDualChuzcFailNorms(workCount, workData, workDataNorm, numVar, workDual,
                          workDualNorm);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workCount = %" HIGHSINT_FORMAT
              "; selectTheta=%g; remainTheta=%g\n",
              workCount, selectTheta, remainTheta);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workDataNorm = %g; workDualNorm = %g\n",
              workDataNorm, workDualNorm);
  return HighsDebugStatus::kOk;
}

HighsDebugStatus debugDualChuzcFailQuad1(
    const HighsOptions& options, const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    const HighsInt numVar, const double* workDual, const double selectTheta,
    const bool force) {
  // Non-trivially expensive assessment of CHUZC failure
  if (options.highs_debug_level < kHighsDebugLevelCostly && !force)
    return HighsDebugStatus::kNotChecked;

  highsLogDev(
      options.log_options, HighsLogType::kInfo,
      "DualChuzC:     No group identified in quad search so return error\n");
  double workDataNorm;
  double workDualNorm;
  debugDualChuzcFailNorms(workCount, workData, workDataNorm, numVar, workDual,
                          workDualNorm);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workCount = %" HIGHSINT_FORMAT
              "; selectTheta=%g\n",
              workCount, selectTheta);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workDataNorm = %g; workDualNorm = %g\n",
              workDataNorm, workDualNorm);
  return HighsDebugStatus::kOk;
}

HighsDebugStatus debugDualChuzcFailHeap(
    const HighsOptions& options, const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    const HighsInt numVar, const double* workDual, const double selectTheta,
    const bool force) {
  // Non-trivially expensive assessment of CHUZC failure
  if (options.highs_debug_level < kHighsDebugLevelCostly && !force)
    return HighsDebugStatus::kNotChecked;

  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     No entries in heap so return error\n");
  double workDataNorm;
  double workDualNorm;
  debugDualChuzcFailNorms(workCount, workData, workDataNorm, numVar, workDual,
                          workDualNorm);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workCount = %" HIGHSINT_FORMAT
              "; selectTheta=%g\n",
              workCount, selectTheta);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workDataNorm = %g; workDualNorm = %g\n",
              workDataNorm, workDualNorm);
  return HighsDebugStatus::kOk;
}

HighsDebugStatus debugNonbasicFlagConsistent(const HighsOptions& options,
                                             const HighsLp& lp,
                                             const SimplexBasis& basis) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  HighsInt numTot = lp.numCol_ + lp.numRow_;
  const bool right_size = (HighsInt)basis.nonbasicFlag_.size() == numTot;
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicFlag size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  HighsInt num_basic_variables = 0;
  for (HighsInt var = 0; var < numTot; var++) {
    if (basis.nonbasicFlag_[var] == kNonbasicFlagFalse) {
      num_basic_variables++;
    } else {
      assert(basis.nonbasicFlag_[var] == kNonbasicFlagTrue);
    }
  }
  bool right_num_basic_variables = num_basic_variables == lp.numRow_;
  if (!right_num_basic_variables) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicFlag has %" HIGHSINT_FORMAT ", not %" HIGHSINT_FORMAT
                 " basic variables\n",
                 num_basic_variables, lp.numRow_);
    assert(right_num_basic_variables);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

// Methods for HMO

const double excessive_absolute_primal_norm = 1e12;
const double excessive_relative_primal_norm = 1e6;
const double large_absolute_primal_norm = sqrt(excessive_absolute_primal_norm);
const double large_relative_primal_norm = sqrt(excessive_relative_primal_norm);

const double excessive_absolute_nonbasic_dual_norm = 1e12;
const double excessive_relative_nonbasic_dual_norm = 1e6;
const double large_absolute_nonbasic_dual_norm =
    sqrt(excessive_absolute_nonbasic_dual_norm);
const double large_relative_nonbasic_dual_norm =
    sqrt(excessive_relative_nonbasic_dual_norm);

const double large_absolute_basic_dual_norm = 1e-12;
const double large_relative_basic_dual_norm = 1e-14;
const double excessive_absolute_basic_dual_norm =
    sqrt(large_absolute_basic_dual_norm);
const double excessive_relative_basic_dual_norm =
    sqrt(large_relative_basic_dual_norm);

const double computed_primal_excessive_absolute_norm =
    excessive_absolute_primal_norm;
const double computed_primal_excessive_relative_norm =
    excessive_relative_primal_norm;
const double computed_primal_large_absolute_norm = large_absolute_primal_norm;
const double computed_primal_large_relative_norm = large_relative_primal_norm;

const double computed_dual_excessive_absolute_nonbasic_dual_norm =
    excessive_absolute_nonbasic_dual_norm;
const double computed_dual_excessive_relative_nonbasic_dual_norm =
    excessive_relative_nonbasic_dual_norm;
const double computed_dual_large_absolute_nonbasic_dual_norm =
    large_absolute_nonbasic_dual_norm;
const double computed_dual_large_relative_nonbasic_dual_norm =
    large_relative_nonbasic_dual_norm;

const double computed_dual_excessive_absolute_basic_dual_norm =
    excessive_absolute_basic_dual_norm;
const double computed_dual_excessive_relative_basic_dual_norm =
    excessive_relative_basic_dual_norm;
const double computed_dual_large_absolute_basic_dual_norm =
    large_absolute_basic_dual_norm;
const double computed_dual_large_relative_basic_dual_norm =
    large_relative_basic_dual_norm;

const double computed_dual_small_relative_nonbasic_dual_change_norm = 1e-12;
const double computed_dual_large_relative_nonbasic_dual_change_norm =
    sqrt(computed_dual_small_relative_nonbasic_dual_change_norm);
const double computed_dual_small_absolute_nonbasic_dual_change_norm = 1e-6;
const double computed_dual_large_absolute_nonbasic_dual_change_norm =
    sqrt(computed_dual_small_absolute_nonbasic_dual_change_norm);

const double updated_objective_small_relative_error = 1e-12;
const double updated_objective_large_relative_error =
    sqrt(updated_objective_small_relative_error);
const double updated_objective_small_absolute_error = 1e-6;
const double updated_objective_large_absolute_error =
    sqrt(updated_objective_small_absolute_error);

const double excessive_basis_condition = 1e16;
const double large_basis_condition = sqrt(excessive_basis_condition);
const double fair_basis_condition = sqrt(large_basis_condition);

const double cleanup_large_absolute_nonbasic_dual_change_norm = 1e-12;
const double cleanup_large_relative_nonbasic_dual_change_norm = 1e-6;
const double cleanup_excessive_absolute_nonbasic_dual_change_norm =
    sqrt(cleanup_large_absolute_nonbasic_dual_change_norm);
const double cleanup_excessive_relative_nonbasic_dual_change_norm =
    sqrt(cleanup_large_relative_nonbasic_dual_change_norm);

const double freelist_excessive_pct_num_entries = 25.0;
const double freelist_large_pct_num_entries = 10.0;
const double freelist_fair_pct_num_entries = 1.0;

/*
HighsDebugStatus debugSimplexLp(const HighsModelObject& highs_model_object) {
  // Non-trivially expensive check that the .lp, if valid is .lp scaled
  // according to .scale
  const HighsSimplexStatus& status =
      highs_model_object.status_;
  if (!status.valid ||
      highs_model_object.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsOptions& options = highs_model_object.options_;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsScale& scale = highs_model_object.scale_;
  const HFactor& factor = highs_model_object.factor_;

  bool right_size = true;
  right_size = (HighsInt)scale.col.size() == lp.numCol_ && right_size;
  right_size = (HighsInt)scale.row.size() == lp.numRow_ && right_size;
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                    "scale size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  // Take a copy of the original LP
  HighsLp check_lp = lp;
  if (applyScalingToLp(options, check_lp, scale) != HighsStatus::kOk) {
    highsLogUser(options.log_options, HighsLogType::kError,
                    "debugSimplexLp: Error scaling check LP\n");
    return HighsDebugStatus::kLogicalError;
  }
  const bool lp_data_ok = check_lp == lp;
  if (!lp_data_ok) {
    highsLogUser(options.log_options, HighsLogType::kError,
                    "debugSimplexLp: Check LP and simplex LP not equal\n");
    assert(lp_data_ok);
    return_status = HighsDebugStatus::kLogicalError;
  }

  if (status.has_basis) {
    const bool basis_correct =
        debugDebugToHighsStatus(debugSimplexBasisCorrect(highs_model_object)) !=
        HighsStatus::kError;
    if (!basis_correct) {
      highsLogUser(options.log_options, HighsLogType::kError,
                      "Supposed to be a Simplex basis, but incorrect\n");
      assert(basis_correct);
      return_status = HighsDebugStatus::kLogicalError;
    }
  }

  if (status.has_invert) {
    const bool invert_ok = debugDebugToHighsStatus(debugCheckInvert(
                               options, factor)) != HighsStatus::kError;
    if (!invert_ok) {
      highsLogUser(options.log_options, HighsLogType::kError,
          "Supposed to be a Simplex basis inverse, but too inaccurate\n");
      assert(invert_ok);
      return_status = HighsDebugStatus::kLogicalError;
    }
  }
  return return_status;
}

HighsDebugStatus debugBasisRightSize(const HighsOptions& options,
                                     const HighsLp& lp,
                                     const SimplexBasis& basis) {
  // Cheap analysis of a Simplex basis, checking vector sizes
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  bool right_size = isBasisRightSize(lp, basis);
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                    "Simplex basis size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugSimplexInfoBasisRightSize(
    const HighsModelObject& highs_model_object) {
  // Trivially cheap check of dimensions and sizes
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;

  const HighsOptions& options = highs_model_object.options_;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsSimplexInfo& info = highs_model_object.info_;
  const SimplexBasis& basis = highs_model_object.basis_;

  HighsInt numCol = lp.numCol_;
  HighsInt numRow = lp.numRow_;
  HighsInt numTot = numCol + numRow;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  bool dimension_ok =
      numCol == lp.numCol_ && numRow == lp.numRow_;
  assert(dimension_ok);
  if (!dimension_ok) {
    highsLogUser(options.log_options, HighsLogType::kError,
        "LP-SimplexLP dimension incompatibility (%" HIGHSINT_FORMAT ", %"
HIGHSINT_FORMAT ") != (%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT ")\n", numCol,
        lp.numCol_, numRow, lp.numRow_);
    return_status = HighsDebugStatus::kLogicalError;
  }
  bool right_size = true;
  right_size = (HighsInt)info.workCost_.size() == numTot && right_size;
  right_size = (HighsInt)info.workDual_.size() == numTot && right_size;
  right_size = (HighsInt)info.workShift_.size() == numTot && right_size;
  right_size = (HighsInt)info.workLower_.size() == numTot && right_size;
  right_size = (HighsInt)info.workUpper_.size() == numTot && right_size;
  right_size = (HighsInt)info.workRange_.size() == numTot && right_size;
  right_size = (HighsInt)info.workValue_.size() == numTot && right_size;
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                    "info work vector size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  if (debugBasisRightSize(options, lp, basis) !=
      HighsDebugStatus::kOk)
    return_status = HighsDebugStatus::kLogicalError;

  return return_status;
}

HighsDebugStatus debugComputePrimal(const HighsModelObject& highs_model_object,
                                    const std::vector<double>& primal_rhs) {
  // Non-trivially expensive analysis of computed primal values.
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kNotChecked;
  const std::vector<double>& primal_value =
      highs_model_object.info_.baseValue_;

  HighsInt num_row = highs_model_object.lp_.numRow_;

  // Use the size of the RHS to determine whether to use it
  const bool have_primal_rhs = (HighsInt)primal_rhs.size() == num_row;

  double primal_rhs_norm = 0;
  if (have_primal_rhs) {
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      primal_rhs_norm += fabs(primal_rhs[iRow]);
  }
  double computed_absolute_primal_norm = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    computed_absolute_primal_norm += fabs(primal_value[iRow]);

  std::string value_adjective;
   HighsLogType report_level;
  return_status = HighsDebugStatus::kOk;
  double computed_relative_primal_norm;
  if (primal_rhs_norm) {
    computed_relative_primal_norm =
        computed_absolute_primal_norm / primal_rhs_norm;
  } else {
    computed_relative_primal_norm = -1;
  }
  if (computed_relative_primal_norm > computed_primal_excessive_relative_norm ||
      computed_absolute_primal_norm > computed_primal_excessive_absolute_norm) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (computed_relative_primal_norm >
                 computed_primal_large_relative_norm ||
             computed_absolute_primal_norm >
                 computed_primal_large_absolute_norm) {
    value_adjective = "Large";
    report_level = HighsLogType::kWarning;
    return_status = HighsDebugStatus::kWarning;
  } else {
    value_adjective = "SMALL";
    report_level = HighsLogType::kVerbose;
  }
  highsLogDev(highs_model_object.options_.log_options, report_level,
      "ComputePrimal: %-9s absolute (%9.4g) or relative (%9.4g) norm of "
      "primal values\n",
      value_adjective.c_str(), computed_absolute_primal_norm,
      computed_relative_primal_norm);
  return return_status;
}
HighsDebugStatus debugComputeDual(const HighsModelObject& highs_model_object,
                                  const std::vector<double>& previous_dual,
                                  const std::vector<double>& basic_costs,
                                  const std::vector<double>& row_dual) {
  // Non-trivially expensive analysis of computed dual values.
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kNotChecked;
  const std::vector<double>& new_dual =
      highs_model_object.info_.workDual_;

  HighsInt num_row = highs_model_object.lp_.numRow_;
  HighsInt num_col = highs_model_object.lp_.numCol_;

  const bool have_basic_costs = (HighsInt)basic_costs.size() == num_row;
  const bool have_row_dual = (HighsInt)row_dual.size() == num_row;
  const bool have_previous_dual =
      (HighsInt)previous_dual.size() == num_col + num_row;

  double basic_costs_norm = 0;
  if (have_basic_costs) {
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      basic_costs_norm += fabs(basic_costs[iRow]);
  }
  double row_dual_norm = 0;
  if (have_row_dual) {
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      row_dual_norm += fabs(row_dual[iRow]);
  }
  double computed_dual_absolute_basic_dual_norm = 0;
  double computed_dual_absolute_nonbasic_dual_norm = 0;
  for (HighsInt iVar = 0; iVar < num_row + num_col; iVar++) {
    if (!highs_model_object.basis_.nonbasicFlag_[iVar]) {
      computed_dual_absolute_basic_dual_norm += fabs(new_dual[iVar]);
      continue;
    }
    computed_dual_absolute_nonbasic_dual_norm += fabs(new_dual[iVar]);
  }
  std::string value_adjective;
   HighsLogType report_level;
  return_status = HighsDebugStatus::kOk;
  // Comment on the norm of the basic costs being zero
  if (have_basic_costs && !basic_costs_norm) {
    highsLogUser(highs_model_object.options_.log_options,
HighsLogType::kWarning, "ComputeDual:   basic cost norm is = %9.4g\n",
basic_costs_norm); return_status = HighsDebugStatus::kWarning;
  }
  // Comment on the norm of the nonbasic duals being zero
  if (!computed_dual_absolute_nonbasic_dual_norm) {
    highsLogUser(highs_model_object.options_.log_options,
                    HighsLogType::kWarning,
                    "ComputeDual:   nonbasic dual norm is = %9.4g\n",
                    computed_dual_absolute_nonbasic_dual_norm);
    return_status = HighsDebugStatus::kWarning;
  }

  // Comment on the norm of basic duals (relative to the norm of the
  // basic costs) which, as c_B-BB^{-1}c_B, should be zero
  double computed_dual_relative_basic_dual_norm;
  if (basic_costs_norm) {
    computed_dual_relative_basic_dual_norm =
        computed_dual_absolute_basic_dual_norm / basic_costs_norm;
  } else {
    computed_dual_relative_basic_dual_norm = -1;
  }
  if (computed_dual_relative_basic_dual_norm >
          computed_dual_excessive_relative_basic_dual_norm ||
      computed_dual_absolute_basic_dual_norm >
          computed_dual_excessive_absolute_basic_dual_norm) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (computed_dual_relative_basic_dual_norm >
                 computed_dual_large_relative_basic_dual_norm ||
             computed_dual_absolute_basic_dual_norm >
                 computed_dual_large_absolute_basic_dual_norm) {
    value_adjective = "Large";
    report_level = HighsLogType::kWarning;
    return_status = HighsDebugStatus::kWarning;
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
  }
  highsLogDev(highs_model_object.options_.log_options, report_level,
      "ComputeDual:   %-9s absolute (%9.4g) or relative (%9.4g) norm of "
      "   basic dual values\n",
      value_adjective.c_str(), computed_dual_absolute_basic_dual_norm,
      computed_dual_relative_basic_dual_norm);
  // Comment on the norm of nonbasic duals relative to the norm of the
  // basic costs
  double computed_dual_relative_nonbasic_dual_norm;
  if (basic_costs_norm) {
    computed_dual_relative_nonbasic_dual_norm =
        computed_dual_absolute_nonbasic_dual_norm / basic_costs_norm;
  } else {
    computed_dual_relative_nonbasic_dual_norm = -1;
  }
  if (computed_dual_relative_nonbasic_dual_norm >
          computed_dual_excessive_relative_nonbasic_dual_norm ||
      computed_dual_absolute_nonbasic_dual_norm >
          computed_dual_excessive_absolute_nonbasic_dual_norm) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (computed_dual_relative_nonbasic_dual_norm >
                 computed_dual_large_relative_nonbasic_dual_norm ||
             computed_dual_absolute_nonbasic_dual_norm >
                 computed_dual_large_absolute_nonbasic_dual_norm) {
    value_adjective = "Large";
    report_level = HighsLogType::kWarning;
    return_status = HighsDebugStatus::kWarning;
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
  }
  highsLogDev(highs_model_object.options_.log_options, report_level,
      "ComputeDual:   %-9s absolute (%9.4g) or relative (%9.4g) norm of "
      "nonbasic dual values\n",
      value_adjective.c_str(), computed_dual_absolute_nonbasic_dual_norm,
      computed_dual_relative_nonbasic_dual_norm);
  double report_basic_costs_norm = -1;
  if (basic_costs_norm) report_basic_costs_norm = basic_costs_norm;
  double report_row_dual_norm = -1;
  if (row_dual_norm) report_row_dual_norm = row_dual_norm;
  highsLogDev(highs_model_object.options_.log_options, report_level,
                    "ComputeDual:   B.pi=c_B has |c_B| = %9.4g; |pi| = %9.4g; "
                    "|pi^TA-c| = [basic %9.4g; nonbasic %9.4g]\n",
                    report_basic_costs_norm, report_row_dual_norm,
                    computed_dual_absolute_basic_dual_norm,
                    computed_dual_absolute_nonbasic_dual_norm);
  if (have_previous_dual) {
    // Comment on the change in the dual values
    std::string change_adjective;
    double computed_dual_absolute_nonbasic_dual_change_norm = 0;
    for (HighsInt iVar = 0; iVar < num_row + num_col; iVar++) {
      if (!highs_model_object.basis_.nonbasicFlag_[iVar]) continue;
      computed_dual_absolute_nonbasic_dual_change_norm +=
          fabs(new_dual[iVar] - previous_dual[iVar]);
    }
    double computed_dual_relative_nonbasic_dual_change_norm;
    if (computed_dual_absolute_nonbasic_dual_norm) {
      computed_dual_relative_nonbasic_dual_change_norm =
          computed_dual_absolute_nonbasic_dual_change_norm /
          computed_dual_absolute_nonbasic_dual_norm;
    } else {
      computed_dual_relative_nonbasic_dual_change_norm = -1;
    }
    if (computed_dual_relative_nonbasic_dual_change_norm >
            computed_dual_large_relative_nonbasic_dual_change_norm ||
        computed_dual_absolute_nonbasic_dual_change_norm >
            computed_dual_large_absolute_nonbasic_dual_change_norm) {
      change_adjective = "Large";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kWarning;
    } else if (computed_dual_relative_nonbasic_dual_change_norm >
                   computed_dual_small_relative_nonbasic_dual_change_norm ||
               computed_dual_absolute_nonbasic_dual_change_norm >
                   computed_dual_small_absolute_nonbasic_dual_change_norm) {
      change_adjective = "Small";
      report_level = HighsLogType::kWarning;
      return_status = HighsDebugStatus::kWarning;
    } else {
      change_adjective = "OK";
      report_level = HighsLogType::kVerbose;
    }
    highsLogDev(highs_model_object.options_.log_options, report_level,
                      "ComputeDual:   %-9s absolute (%9.4g) or relative "
                      "(%9.4g) nonbasic dual change\n",
                      change_adjective.c_str(),
                      computed_dual_absolute_nonbasic_dual_change_norm,
                      computed_dual_relative_nonbasic_dual_change_norm);
  }
  return return_status;
}

HighsDebugStatus debugSimplexDualFeasibility(
    const HighsModelObject& highs_model_object, const std::string message,
    const bool force) {
  // Non-trivially expensive check of dual feasibility.
  if (highs_model_object.options_.highs_debug_level <
          kHighsDebugLevelCostly &&
      !force)
    return HighsDebugStatus::kNotChecked;
  if (force)
    highsLogDev(highs_model_object.options_.log_options, HighsLogType::kInfo,
                      "SmplxDuFeas:   Forcing debug\n");

  const HighsLp& lp = highs_model_object.lp_;
  const HighsSimplexInfo& info = highs_model_object.info_;
  const SimplexBasis& basis = highs_model_object.basis_;
  double scaled_dual_feasibility_tolerance =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;

  HighsInt num_dual_infeasibilities = 0;
  double sum_dual_infeasibilities = 0;
  double max_dual_infeasibility = 0;
  for (HighsInt iVar = 0; iVar < lp.numCol_ + lp.numRow_;
iVar++) { if (!basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double dual = info.workDual_[iVar];
    const double lower = info.workLower_[iVar];
    const double upper = info.workUpper_[iVar];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else {
      // Not free: any dual infeasibility is given by the dual value
      // signed by nonbasicMove
      dual_infeasibility = -basis.nonbasicMove_[iVar] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
  if (num_dual_infeasibilities) {
    highsLogDev(highs_model_object.options_.log_options, HighsLogType::kError,
                      "SmplxDuFeas:   num/max/sum simplex dual infeasibilities "
                      "= %" HIGHSINT_FORMAT " / %g / %g - %s\n",
                      num_dual_infeasibilities, max_dual_infeasibility,
                      sum_dual_infeasibilities, message.c_str());
    return HighsDebugStatus::kLogicalError;
  }
  return HighsDebugStatus::kOk;
}

HighsDebugStatus debugUpdatedObjectiveValue(
    HighsModelObject& highs_model_object, const SimplexAlgorithm algorithm,
    const HighsInt phase, const std::string message, const bool force) {
  // Non-trivially expensive check of updated objective value. Computes the
  // exact objective value
  if (highs_model_object.options_.highs_debug_level <
          kHighsDebugLevelCostly &&
      !force)
    return HighsDebugStatus::kNotChecked;
  HighsSimplexInfo& info = highs_model_object.info_;

  static bool have_previous_exact_primal_objective_value;
  static double previous_exact_primal_objective_value;
  static double previous_updated_primal_objective_value;
  static double updated_primal_objective_correction;

  static bool have_previous_exact_dual_objective_value;
  static double previous_exact_dual_objective_value;
  static double previous_updated_dual_objective_value;
  static double updated_dual_objective_correction;
  if (phase < 0) {
    if (algorithm == SimplexAlgorithm::kPrimal) {
      have_previous_exact_primal_objective_value = false;
    } else {
      have_previous_exact_dual_objective_value = false;
    }
    return HighsDebugStatus::kOk;
  }
  double exact_objective_value;
  double updated_objective_value;
  bool have_previous_exact_objective_value;
  // Assign values to prevent compiler warning
  double previous_exact_objective_value = 0;
  double previous_updated_objective_value = 0;
  double updated_objective_correction = 0;
  std::string algorithm_name;
  if (algorithm == SimplexAlgorithm::kPrimal) {
    algorithm_name = "primal";
    have_previous_exact_objective_value =
        have_previous_exact_primal_objective_value;
    if (have_previous_exact_objective_value) {
      previous_exact_objective_value = previous_exact_primal_objective_value;
      previous_updated_objective_value =
          previous_updated_primal_objective_value;
      updated_objective_correction = updated_primal_objective_correction;
    }
    updated_objective_value = info.updated_primal_objective_value;
    // Save the current objective value so that it can be recovered
    // after calling computePrimalObjectiveValue
    double save_objective_value = info.primal_objective_value;
    computePrimalObjectiveValue(highs_model_object);
    exact_objective_value = info.primal_objective_value;
    info.primal_objective_value = save_objective_value;
  } else {
    algorithm_name = "dual";
    have_previous_exact_objective_value =
        have_previous_exact_dual_objective_value;
    if (have_previous_exact_objective_value) {
      previous_exact_objective_value = previous_exact_dual_objective_value;
      previous_updated_objective_value = previous_updated_dual_objective_value;
      updated_objective_correction = updated_dual_objective_correction;
    }
    updated_objective_value = info.updated_dual_objective_value;
    // Save the current objective value so that it can be recovered
    // after calling computeDualObjectiveValue
    double save_objective_value = info.dual_objective_value;
    computeDualObjectiveValue(highs_model_object, phase);
    exact_objective_value = info.dual_objective_value;
    info.dual_objective_value = save_objective_value;
  }
  double change_in_objective_value = 0;
  double change_in_updated_objective_value = 0;
  if (have_previous_exact_objective_value) {
    change_in_objective_value =
        exact_objective_value - previous_exact_objective_value;
    change_in_updated_objective_value =
        updated_objective_value - previous_updated_objective_value;
    updated_objective_value += updated_objective_correction;
  } else {
    updated_objective_correction = 0;
  }
  const double updated_objective_error =
      exact_objective_value - updated_objective_value;
  const double updated_objective_absolute_error = fabs(updated_objective_error);
  const double updated_objective_relative_error =
      updated_objective_absolute_error / max(1.0, fabs(exact_objective_value));
  updated_objective_correction += updated_objective_error;

  // Now update the records of previous objective value
  if (algorithm == SimplexAlgorithm::kPrimal) {
    have_previous_exact_primal_objective_value = true;
    previous_exact_primal_objective_value = exact_objective_value;
    previous_updated_primal_objective_value = updated_objective_value;
    updated_primal_objective_correction = updated_objective_correction;
  } else {
    have_previous_exact_dual_objective_value = true;
    previous_exact_dual_objective_value = exact_objective_value;
    previous_updated_dual_objective_value = updated_objective_value;
    updated_dual_objective_correction = updated_objective_correction;
  }

  // Now analyse the error
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  std::string error_adjective;
   HighsLogType report_level;
  bool at_least_small_error =
      updated_objective_relative_error >
          updated_objective_small_relative_error ||
      updated_objective_absolute_error > updated_objective_small_absolute_error;
  if (!at_least_small_error) return return_status;
  if (updated_objective_relative_error >
          updated_objective_large_relative_error ||
      updated_objective_absolute_error >
          updated_objective_large_absolute_error) {
    error_adjective = "Large";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kLargeError;
  } else if (updated_objective_relative_error >
                 updated_objective_small_relative_error ||
             updated_objective_absolute_error >
                 updated_objective_small_absolute_error) {
    error_adjective = "Small";
    report_level = HighsLogType::kWarning;
    return_status = HighsDebugStatus::kSmallError;
  } else {
    error_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  highsLogDev(highs_model_object.options_.log_options, report_level,
      "UpdateObjVal:  %-9s absolute (%9.4g) or relative (%9.4g) error in "
      "updated %s objective value"
      " - objective change - exact (%9.4g) updated (%9.4g) | %s\n",
      error_adjective.c_str(), updated_objective_error,
      updated_objective_relative_error, algorithm_name.c_str(),
      change_in_objective_value, change_in_updated_objective_value,
      message.c_str());
  return return_status;
}

HighsDebugStatus debugUpdatedObjectiveValue(
    const HighsModelObject& highs_model_object,
    const SimplexAlgorithm algorithm) {
  // Cheap check of updated objective value - assumes that the
  // objective value computed directly is correct, so only call after
  // this has been done
  if (highs_model_object.options_.highs_debug_level == kHighsDebugLevelNone)
    return HighsDebugStatus::kNotChecked;
  const HighsSimplexInfo& info = highs_model_object.info_;
  std::string algorithm_name = "dual";
  if (algorithm == SimplexAlgorithm::kPrimal) algorithm_name = "primal";
  double exact_objective_value;
  double updated_objective_value;
  if (algorithm == SimplexAlgorithm::kPrimal) {
    assert(highs_model_object.status_.has_primal_objective_value);
    exact_objective_value = info.primal_objective_value;
    updated_objective_value = info.updated_primal_objective_value;
  } else {
    assert(highs_model_object.status_.has_dual_objective_value);
    exact_objective_value = info.dual_objective_value;
    updated_objective_value = info.updated_dual_objective_value;
  }
  const double updated_objective_error =
      exact_objective_value - updated_objective_value;
  const double updated_objective_absolute_error = fabs(updated_objective_error);
  const double updated_objective_relative_error =
      updated_objective_absolute_error / max(1.0, fabs(exact_objective_value));

  // Now analyse the error
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  std::string error_adjective;
   HighsLogType report_level;
  if (updated_objective_relative_error >
          updated_objective_large_relative_error ||
      updated_objective_absolute_error >
          updated_objective_large_absolute_error) {
    error_adjective = "Large";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kLargeError;
  } else if (updated_objective_relative_error >
                 updated_objective_small_relative_error ||
             updated_objective_absolute_error >
                 updated_objective_small_absolute_error) {
    error_adjective = "Small";
    report_level = HighsLogType::kWarning;
    return_status = HighsDebugStatus::kSmallError;
  } else {
    error_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  highsLogDev(highs_model_object.options_.log_options, report_level,
                    "UpdateObjVal:  %-9s large absolute (%9.4g) or relative "
                    "(%9.4g) error in updated %s objective value\n",
                    error_adjective.c_str(), updated_objective_error,
                    updated_objective_relative_error, algorithm_name.c_str());
  return return_status;
}

HighsDebugStatus debugFixedNonbasicMove(
    const HighsModelObject& highs_model_object) {
  // Non-trivially expensive check of nonbasicMove for fixed variables
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsSimplexInfo& info = highs_model_object.info_;
  const SimplexBasis& basis = highs_model_object.basis_;
  HighsInt num_fixed_variable_move_errors = 0;
  for (HighsInt iVar = 0; iVar < lp.numCol_ + lp.numRow_;
iVar++) { if (!basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    if (info.workLower_[iVar] == info.workUpper_[iVar] &&
        basis.nonbasicMove_[iVar])
      num_fixed_variable_move_errors++;
  }
  assert(num_fixed_variable_move_errors == 0);
  if (num_fixed_variable_move_errors) {
    highsLogDev(highs_model_object.options_.log_options, HighsLogType::kError,
                      "There are %" HIGHSINT_FORMAT " fixed nonbasicMove
errors", num_fixed_variable_move_errors); return
HighsDebugStatus::kLogicalError;
  }
  return HighsDebugStatus::kOk;
}

HighsDebugStatus debugNonbasicMove(const HighsModelObject& highs_model_object) {
  // Non-trivially expensive check of NonbasicMove
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsOptions& options = highs_model_object.options_;
  const HighsLp& lp = highs_model_object.lp_;
  const SimplexBasis& basis = highs_model_object.basis_;
  HighsInt num_free_variable_move_errors = 0;
  HighsInt num_lower_bounded_variable_move_errors = 0;
  HighsInt num_upper_bounded_variable_move_errors = 0;
  HighsInt num_boxed_variable_move_errors = 0;
  HighsInt num_fixed_variable_move_errors = 0;
  const HighsInt numTot = lp.numCol_ + lp.numRow_;
  bool right_size = (HighsInt)basis.nonbasicMove_.size() == numTot;
  // Check consistency of nonbasicMove
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                    "nonbasicMove size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  double lower;
  double upper;

  for (HighsInt iVar = 0; iVar < numTot; iVar++) {
    if (!basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic variable
    if (iVar < lp.numCol_) {
      lower = lp.colLower_[iVar];
      upper = lp.colUpper_[iVar];
    } else {
      HighsInt iRow = iVar - lp.numCol_;
      lower = -lp.rowUpper_[iRow];
      upper = -lp.rowLower_[iRow];
    }

    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free
        if (basis.nonbasicMove_[iVar]) {
          num_free_variable_move_errors++;
        }
      } else {
        // Only lower bounded
        if (basis.nonbasicMove_[iVar] != kNonbasicMoveUp) {
          num_lower_bounded_variable_move_errors++;
        }
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded
        if (basis.nonbasicMove_[iVar] != kNonbasicMoveDn) {
          num_upper_bounded_variable_move_errors++;
        }
      } else {
        // Boxed or fixed
        if (lower != upper) {
          // Boxed
          if (!basis.nonbasicMove_[iVar]) {
            num_boxed_variable_move_errors++;
          }
        } else {
          // Fixed
          if (basis.nonbasicMove_[iVar]) {
            num_fixed_variable_move_errors++;
          }
        }
      }
    }
  }
  HighsInt num_errors =
      num_free_variable_move_errors + num_lower_bounded_variable_move_errors +
      num_upper_bounded_variable_move_errors + num_boxed_variable_move_errors +
      num_fixed_variable_move_errors;

  if (num_errors) {
    highsLogUser(options.log_options, HighsLogType::kError,
        "There are %" HIGHSINT_FORMAT " nonbasicMove errors: %" HIGHSINT_FORMAT
" free; %" HIGHSINT_FORMAT " lower; %" HIGHSINT_FORMAT " upper; %"
HIGHSINT_FORMAT " \n" "boxed; %" HIGHSINT_FORMAT " fixed", num_errors,
num_free_variable_move_errors, num_lower_bounded_variable_move_errors,
        num_upper_bounded_variable_move_errors, num_boxed_variable_move_errors,
        num_fixed_variable_move_errors);
    assert(num_errors == 0);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugBasisCondition(const HighsModelObject& highs_model_object,
                                     const std::string message) {
  // Non-trivially expensive assessment of basis condition
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  double basis_condition = computeBasisCondition(highs_model_object);
  std::string value_adjective;
   HighsLogType report_level;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  if (basis_condition > excessive_basis_condition) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (basis_condition > large_basis_condition) {
    value_adjective = "Large";
    report_level = HighsLogType::kWarning;
    return_status = HighsDebugStatus::kWarning;
  } else if (basis_condition > fair_basis_condition) {
    value_adjective = "Fair";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  highsLogDev(highs_model_object.options_.log_options, report_level,
      "BasisCond:     %-9s basis condition estimate (%9.4g) - %s\n",
      value_adjective.c_str(), basis_condition, message.c_str());
  return return_status;
}

HighsDebugStatus debugCleanup(HighsModelObject& highs_model_object,
                              const std::vector<double>& original_dual) {
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsSimplexInfo& info = highs_model_object.info_;
  const SimplexBasis& basis = highs_model_object.basis_;
  // Make sure that the original_dual has been set up
  assert((HighsInt)original_dual.size() == lp.numCol_ +
lp.numRow_); const std::vector<double>& new_dual =
info.workDual_;

  const double dual_feasibility_tolerance =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;
  HighsInt num_dual_sign_change = 0;
  double cleanup_absolute_nonbasic_dual_change_norm = 0;
  double cleanup_absolute_nonbasic_dual_norm = 0;
  for (HighsInt iVar = 0; iVar < lp.numCol_ + lp.numRow_;
iVar++) { if (!basis.nonbasicFlag_[iVar]) continue;
    cleanup_absolute_nonbasic_dual_norm += std::fabs(new_dual[iVar]);
#ifdef HiGHSDEV
    const double nonbasic_dual_change =
        std::fabs(new_dual[iVar] - original_dual[iVar]);
    cleanup_absolute_nonbasic_dual_change_norm += nonbasic_dual_change;
#endif
    const double max_dual =
        std::max(std::fabs(new_dual[iVar]), std::fabs(original_dual[iVar]));
    if (max_dual > dual_feasibility_tolerance &&
        new_dual[iVar] * original_dual[iVar] < 0)
      num_dual_sign_change++;
  }
  // Comment on the norm of the nonbasic duals being zero
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  if (!cleanup_absolute_nonbasic_dual_norm) {
    highsLogUser(highs_model_object.options_.log_options,
                    HighsLogType::kWarning,
                    "DualCleanup:   dual norm is = %9.4g\n",
                    cleanup_absolute_nonbasic_dual_norm);
    return_status = HighsDebugStatus::kWarning;
  }
  // Comment on the norm of the change being zero
  if (!cleanup_absolute_nonbasic_dual_change_norm) {
    highsLogUser(highs_model_object.options_.log_options,
                    HighsLogType::kWarning,
                    "DualCleanup:   dual norm is = %9.4g\n",
                    cleanup_absolute_nonbasic_dual_change_norm);
    return_status = HighsDebugStatus::kWarning;
  }
  double cleanup_relative_nonbasic_dual_change_norm;
  if (cleanup_absolute_nonbasic_dual_norm) {
    cleanup_relative_nonbasic_dual_change_norm =
        cleanup_absolute_nonbasic_dual_change_norm /
        cleanup_absolute_nonbasic_dual_norm;
  } else {
    cleanup_relative_nonbasic_dual_change_norm = -1;
  }
  std::string value_adjective;
   HighsLogType report_level;
  if (cleanup_absolute_nonbasic_dual_change_norm >
          cleanup_excessive_absolute_nonbasic_dual_change_norm ||
      cleanup_relative_nonbasic_dual_change_norm >
          cleanup_excessive_relative_nonbasic_dual_change_norm) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (cleanup_absolute_nonbasic_dual_change_norm >
                 cleanup_large_absolute_nonbasic_dual_change_norm ||
             cleanup_relative_nonbasic_dual_change_norm >
                 cleanup_large_relative_nonbasic_dual_change_norm) {
    value_adjective = "Large";
    report_level = HighsLogType::kWarning;
    return_status = HighsDebugStatus::kWarning;
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  highsLogDev(highs_model_object.options_.log_options, report_level,
      "DualCleanup:   %-9s absolute (%9.4g) or relative (%9.4g) dual change, "
      "with %" HIGHSINT_FORMAT " meaningful sign change(s)\n",
      value_adjective.c_str(), cleanup_absolute_nonbasic_dual_change_norm,
      cleanup_relative_nonbasic_dual_change_norm, num_dual_sign_change);
  return return_status;
}

HighsDebugStatus debugFreeListNumEntries(
    const HighsModelObject& highs_model_object, const std::set<HighsInt>&
freeList) { if (highs_model_object.options_.highs_debug_level <
kHighsDebugLevelCheap) return HighsDebugStatus::kNotChecked;

  HighsInt freelist_num_entries = 0;
  if (freeList.size() > 0) {
    std::set<HighsInt>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++)
      freelist_num_entries++;
  }

  const HighsInt numTot = highs_model_object.lp_.numCol_ +
                     highs_model_object.lp_.numRow_;
  double pct_freelist_num_entries = (100.0 * freelist_num_entries) / numTot;

  std::string value_adjective;
   HighsLogType report_level;
  HighsDebugStatus return_status = HighsDebugStatus::kNotChecked;

  if (pct_freelist_num_entries > freelist_excessive_pct_num_entries) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
  } else if (pct_freelist_num_entries > freelist_large_pct_num_entries) {
    value_adjective = "Large";
    report_level = HighsLogType::kWarning;
  } else if (pct_freelist_num_entries > freelist_fair_pct_num_entries) {
    value_adjective = "Fair";
    report_level = HighsLogType::kVerbose;
  } else {
    value_adjective = "OK";
    if (freelist_num_entries) {
      report_level = HighsLogType::kError;
    } else {
      report_level = HighsLogType::kVerbose;
    }
    return_status = HighsDebugStatus::kOk;
  }

  highsLogDev(highs_model_object.options_.log_options, report_level,
      "FreeList   :   %-9s percentage (%6.2g) of %" HIGHSINT_FORMAT " variables
on free list\n", value_adjective.c_str(), pct_freelist_num_entries, numTot);

  return return_status;
}

void debugDualChuzcWorkDataAndGroupReport(
    const HighsModelObject& highs_model_object, const double workDelta,
    const double workTheta, const std::string message,
    const HighsInt report_workCount,
    const std::vector<std::pair<HighsInt, double>>& report_workData,
    const std::vector<HighsInt>& report_workGroup) {
  const HighsOptions& options = highs_model_object.options_;
  const std::vector<HighsInt>& workMove =
      highs_model_object.basis_.nonbasicMove_;
  const std::vector<double>& workDual =
      highs_model_object.info_.workDual_;
  const std::vector<double>& workRange =
      highs_model_object.info_.workRange_;
  const double Td =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;
  double totalChange = initial_total_change;
  const double totalDelta = fabs(workDelta);
  highsLogDev(options.log_options, HighsLogType::kInfo,
  "\n%s: totalDelta = %10.4g\nworkData\n  En iCol       Dual      Value    "
      "  Ratio     Change\n",
      message.c_str(), totalDelta);
  for (HighsInt i = 0; i < report_workCount; i++) {
    HighsInt iCol = report_workData[i].first;
    double value = report_workData[i].second;
    double dual = workMove[iCol] * workDual[iCol];
    totalChange += value * (workRange[iCol]);
    highsLogDev(options.log_options, HighsLogType::kInfo,
                      "%4" HIGHSINT_FORMAT " %4" HIGHSINT_FORMAT " %10.4g %10.4g
%10.4g %10.4g\n", i, iCol, dual, value, dual / value, totalChange);
  }
  double selectTheta = workTheta;
  highsLogDev(options.log_options, HighsLogType::kInfo,
                    "workGroup\n  Ix:   selectTheta Entries\n");
  for (HighsInt group = 0; group < (HighsInt)report_workGroup.size() - 1;
group++) { highsLogDev(options.log_options, HighsLogType::kInfo,
                      "%4" HIGHSINT_FORMAT ": selectTheta = %10.4g ", group,
selectTheta); for (HighsInt en = report_workGroup[group]; en <
report_workGroup[group + 1]; en++) { highsLogDev(options.log_options,
HighsLogType::kInfo,
                        "%4" HIGHSINT_FORMAT " ", en);
    }
    highsLogDev(options.log_options, HighsLogType::kInfo, "\n");
    HighsInt en = report_workGroup[group + 1];
    HighsInt iCol = report_workData[en].first;
    double value = report_workData[en].second;
    double dual = workMove[iCol] * workDual[iCol];
    selectTheta = (dual + Td) / value;
  }
}

HighsDebugStatus debugDualChuzcWorkDataAndGroup(
    const HighsModelObject& highs_model_object, const double workDelta,
    const double workTheta, const HighsInt workCount, const HighsInt
alt_workCount, const HighsInt breakIndex, const HighsInt alt_breakIndex, const
std::vector<std::pair<HighsInt, double>>& workData, const
std::vector<std::pair<HighsInt, double>>& sorted_workData, const
std::vector<HighsInt>& workGroup, const std::vector<HighsInt>& alt_workGroup) {
  // Cheap comparison and possible non-trivially expensive reporting
  // of the two sorting methods for BFRT nodes in dual CHUZC
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  const HighsOptions& options = highs_model_object.options_;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  HighsInt workPivot = workData[breakIndex].first;
  HighsInt alt_workPivot = sorted_workData[alt_breakIndex].first;
  if (alt_workPivot != workPivot) {
    highsLogDev(options.log_options, HighsLogType::kInfo,
                      "Quad workPivot = %" HIGHSINT_FORMAT "; Heap workPivot =
%" HIGHSINT_FORMAT "\n", workPivot, alt_workPivot); return_status =
HighsDebugStatus::kWarning; if (highs_model_object.options_.highs_debug_level <
        kHighsDebugLevelCostly)
      return return_status;
    debugDualChuzcWorkDataAndGroupReport(highs_model_object, workDelta,
                                         workTheta, "Original", workCount,
                                         workData, workGroup);
    debugDualChuzcWorkDataAndGroupReport(
        highs_model_object, workDelta, workTheta, "Heap-derived", alt_workCount,
        sorted_workData, alt_workGroup);
  }
  return return_status;
}

HighsDebugStatus debugSimplexBasicSolution(
    const string message, const HighsModelObject& highs_model_object) {
  // Non-trivially expensive analysis of a simplex basic solution, starting from
  // solution_params
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;

  if (highsStatusFromHighsModelStatus(
          highs_model_object.scaled_model_status_) != HighsStatus::kOk)
    return HighsDebugStatus::kOk;
  HighsDebugStatus return_status = HighsDebugStatus::kNotChecked;

  const HighsLp& lp = highs_model_object.lp_;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsScale& scale = highs_model_object.scale_;
  const HighsSimplexInfo& info = highs_model_object.info_;
  const SimplexBasis& basis = highs_model_object.basis_;

  return_status = debugSimplexInfoBasisRightSize(highs_model_object);
  if (return_status == HighsDebugStatus::kLogicalError) return return_status;

  // Determine a HiGHS basis from the simplex basis. Only basic/nonbasic is
  // needed
  HighsBasis highs_basis;
  highs_basis.col_status.resize(lp.numCol_);
  highs_basis.row_status.resize(lp.numRow_);
  // Now scatter the indices of basic variables
  for (HighsInt iVar = 0; iVar < lp.numCol_ + lp.numRow_; iVar++) {
    if (iVar < lp.numCol_) {
      HighsInt iCol = iVar;
      if (basis.nonbasicFlag_[iVar] == kNonbasicFlagTrue) {
        highs_basis.col_status[iCol] = HighsBasisStatus::kNonbasic;
      } else {
        highs_basis.col_status[iCol] = HighsBasisStatus::kBasic;
      }
    } else {
      HighsInt iRow = iVar - lp.numCol_;
      if (basis.nonbasicFlag_[iVar] == kNonbasicFlagTrue) {
        highs_basis.row_status[iRow] = HighsBasisStatus::kNonbasic;
      } else {
        highs_basis.row_status[iRow] = HighsBasisStatus::kBasic;
      }
    }
  }
  highs_basis.valid_ = true;
  // Possibly scaled model
  // Determine a HiGHS solution simplex solution
  HighsSolution solution;
  solution.col_value.resize(lp.numCol_);
  solution.col_dual.resize(lp.numCol_);
  solution.row_value.resize(lp.numRow_);
  solution.row_dual.resize(lp.numRow_);

  for (HighsInt iVar = 0; iVar < lp.numCol_ + lp.numRow_; iVar++) {
    if (iVar < lp.numCol_) {
      HighsInt iCol = iVar;
      solution.col_value[iCol] = info.workValue_[iVar];
      solution.col_dual[iCol] =
          (HighsInt)lp.sense_ * info.workDual_[iVar];
    } else {
      HighsInt iRow = iVar - lp.numCol_;
      solution.row_value[iRow] = -info.workValue_[iVar];
      solution.row_dual[iRow] =
          (HighsInt)lp.sense_ * info.workDual_[iVar];
    }
  }
  // Now insert the basic values
  for (HighsInt ix = 0; ix < lp.numRow_; ix++) {
    HighsInt iVar = basis.basicIndex_[ix];
    if (iVar < lp.numCol_) {
      solution.col_value[iVar] = info.baseValue_[ix];
      solution.col_dual[iVar] = 0;
    } else {
      HighsInt iRow = iVar - lp.numCol_;
      solution.row_value[iRow] = -info.baseValue_[ix];
      solution.row_dual[iRow] = 0;
    }
  }

  const std::string message_scaled = message + " - scaled";
  return_status = debugWorseStatus(
      debugHighsBasicSolution(message_scaled, highs_model_object.options_,
                              lp, basis, solution,
                              highs_model_object.scaled_solution_params_,
                              highs_model_object.scaled_model_status_),
      return_status);

  if (!highs_model_object.scale_.is_scaled) return return_status;

  // Doesn't work if simplex LP has permuted columns
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    solution.col_value[iCol] *= scale.col[iCol];
    solution.col_dual[iCol] /= (scale.col[iCol] / scale.cost);
  }
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    solution.row_value[iRow] /= scale.row[iRow];
    solution.row_dual[iRow] *= (scale.row[iRow] * scale.cost);
  }
  // Cannot assume unscaled solution params or unscaled model status are known
  const std::string message_unscaled = message + " - unscaled";
  return_status = debugWorseStatus(
      debugHighsBasicSolution(message_unscaled, highs_model_object.options_, lp,
                              basis, solution),
      return_status);

  // Scaled model
  return return_status;
}

HighsDebugStatus debugSimplexHighsSolutionDifferences(
    const HighsModelObject& highs_model_object) {
  // Nontrivially expensive check of dimensions and sizes
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;

  const HighsOptions& options = highs_model_object.options_;
  const HighsSolution& solution = highs_model_object.solution_;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsSimplexInfo& info = highs_model_object.info_;
  const SimplexBasis& basis = highs_model_object.basis_;
  const HighsScale& scale = highs_model_object.scale_;

  HighsDebugStatus return_status = HighsDebugStatus::kNotChecked;

  // Go through the columns, finding the differences in nonbasic column values
  // and duals
  double max_nonbasic_col_value_difference = 0;
  double max_nonbasic_col_dual_difference = 0;
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    HighsInt iVar = iCol;
    if (basis.nonbasicFlag_[iVar] == kNonbasicFlagTrue) {
      // Consider this nonbasic column
      double local_col_value = info.workValue_[iVar] * scale.col[iCol];
      double local_col_dual = (HighsInt)lp.sense_ *
                              info.workDual_[iVar] /
                              (scale.col[iCol] / scale.cost);
      double value_difference =
          fabs(local_col_value - solution.col_value[iCol]);
      double dual_difference = fabs(local_col_dual - solution.col_dual[iCol]);
      max_nonbasic_col_value_difference =
          std::max(value_difference, max_nonbasic_col_value_difference);
      max_nonbasic_col_dual_difference =
          std::max(dual_difference, max_nonbasic_col_dual_difference);
    }
  }
  // Go through the rows, finding the differences in nonbasic and
  // basic row values and duals, as well as differences in basic
  // column values and duals
  double max_nonbasic_row_value_difference = 0;
  double max_nonbasic_row_dual_difference = 0;
  double max_basic_col_value_difference = 0;
  double max_basic_col_dual_difference = 0;
  double max_basic_row_value_difference = 0;
  double max_basic_row_dual_difference = 0;

  for (HighsInt ix = 0; ix < lp.numRow_; ix++) {
    HighsInt iRow = ix;
    HighsInt iVar = lp.numCol_ + iRow;
    if (basis.nonbasicFlag_[iVar] == kNonbasicFlagTrue) {
      // Consider this nonbasic row
      double local_row_value =
          -info.workValue_[iVar] / scale.row[iRow];
      double local_row_dual = (HighsInt)lp.sense_ *
                              info.workDual_[iVar] *
                              (scale.row[iRow] * scale.cost);
      double value_difference =
          fabs(local_row_value - solution.row_value[iRow]);
      double dual_difference = fabs(local_row_dual - solution.row_dual[iRow]);
      max_nonbasic_row_value_difference =
          std::max(value_difference, max_nonbasic_row_value_difference);
      max_nonbasic_row_dual_difference =
          std::max(dual_difference, max_nonbasic_row_dual_difference);
    }
    // Consider the basic variable associated with this row index
    iVar = basis.basicIndex_[ix];
    if (iVar < lp.numCol_) {
      // Consider this basic column
      HighsInt iCol = iVar;
      double local_col_value = info.baseValue_[ix] * scale.col[iCol];
      double local_col_dual = 0;
      double value_difference =
          fabs(local_col_value - solution.col_value[iCol]);
      double dual_difference = fabs(local_col_dual - solution.col_dual[iCol]);
      max_basic_col_value_difference =
          std::max(value_difference, max_basic_col_value_difference);
      max_basic_col_dual_difference =
          std::max(dual_difference, max_basic_col_dual_difference);
    } else {
      // Consider this basic row
      iRow = iVar - lp.numCol_;
      double local_row_value = -info.baseValue_[ix] / scale.row[iRow];
      double local_row_dual = 0;
      double value_difference =
          fabs(local_row_value - solution.row_value[iRow]);
      double dual_difference = fabs(local_row_dual - solution.row_dual[iRow]);
      max_basic_row_value_difference =
          std::max(value_difference, max_basic_row_value_difference);
      max_basic_row_dual_difference =
          std::max(dual_difference, max_basic_row_dual_difference);
    }
  }

  highsLogDev(options.log_options, HighsLogType::kInfo,
                    "\nHiGHS-simplex solution differences\n");
  std::string value_adjective;
   HighsLogType report_level;
  return_status = HighsDebugStatus::kOk;
  if (max_nonbasic_col_value_difference > 0) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
    highsLogDev(options.log_options, report_level,
        "HighsSimplexD: %-9s Nonbasic column value difference: %9.4g\n",
        value_adjective.c_str(), max_nonbasic_col_value_difference);
  }
  if (max_nonbasic_row_value_difference > 0) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
    highsLogDev(options.log_options, report_level,
        "HighsSimplexD: %-9s Nonbasic row    value difference: %9.4g\n",
        value_adjective.c_str(), max_nonbasic_row_value_difference);
  }

  return_status = debugWorseStatus(
      debugAssessSolutionNormDifference(options, "Basic   column value",
                                        max_basic_col_value_difference),
      return_status);
  return_status = debugWorseStatus(
      debugAssessSolutionNormDifference(options, "Basic      row value",
                                        max_basic_row_value_difference),
      return_status);
  return_status = debugWorseStatus(
      debugAssessSolutionNormDifference(options, "Nonbasic column dual",
                                        max_nonbasic_col_dual_difference),
      return_status);
  return_status = debugWorseStatus(
      debugAssessSolutionNormDifference(options, "Nonbasic    row dual",
                                        max_nonbasic_row_dual_difference),
      return_status);

  if (max_basic_col_dual_difference > 0) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
    highsLogDev(options.log_options, report_level,
        "HighsSimplexD: %-9s Basic    column dual difference: %9.4g\n",
        value_adjective.c_str(), max_basic_col_dual_difference);
  }
  if (max_basic_row_dual_difference > 0) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
    highsLogDev(options.log_options, report_level,
        "HighsSimplexD: %-9s Basic    row     dual difference: %9.4g\n",
        value_adjective.c_str(), max_basic_row_dual_difference);
  }

  return return_status;
}

HighsDebugStatus debugAssessSolutionNormDifference(const HighsOptions& options,
                                                   const std::string type,
                                                   const double difference) {
  const double small_difference = 1e-12;
  const double large_difference = 1e-8;
  const double excessive_difference = 1e-4;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  if (difference <= small_difference) return return_status;
  std::string value_adjective;
   HighsLogType report_level;

  if (difference > excessive_difference) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (difference > large_difference) {
    value_adjective = "Large";
    report_level = HighsLogType::kWarning;
    return_status = HighsDebugStatus::kWarning;
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
  }
  highsLogDev(options.log_options,report_level,
                    "HighsSimplexD: %-9s %s difference: %9.4g\n",
                    value_adjective.c_str(), type.c_str(), difference);
  return return_status;
}

HighsDebugStatus debugSimplexBasisCorrect(
    const HighsModelObject& highs_model_object) {
  // Nontrivially expensive analysis of a Simplex basis, checking
  // consistency, and then correctness of nonbasicMove
  const HighsOptions& options = highs_model_object.options_;
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  const HighsLp& lp = highs_model_object.lp_;
  const SimplexBasis& basis = highs_model_object.basis_;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const bool consistent =
      debugBasisConsistent(options, lp, basis) !=
      HighsDebugStatus::kLogicalError;
  if (!consistent) {
    highsLogUser(options.log_options, HighsLogType::kError,
                    "Supposed to be a Simplex basis, but not consistent\n");
    assert(consistent);
    return_status = HighsDebugStatus::kLogicalError;
  }
  if (options.highs_debug_level < kHighsDebugLevelCostly)
    return return_status;
  const bool correct_nonbasicMove =
      debugNonbasicMove(highs_model_object) != HighsDebugStatus::kLogicalError;
  if (!correct_nonbasicMove) {
    highsLogUser(options.log_options, HighsLogType::kError,
        "Supposed to be a Simplex basis, but nonbasicMove is incorrect\n");
    assert(correct_nonbasicMove);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugOkForSolve(const HighsModelObject& highs_model_object,
                                 const HighsInt phase) {
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  const HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsSimplexStatus& status =
      highs_model_object.status_;
  const SimplexBasis& basis = highs_model_object.basis_;
  const HighsOptions& options = highs_model_object.options_;
  bool ok;
  // Minimal check - just look at flags. This means we trust them!
  ok = status.has_basis && status.has_matrix &&
       status.has_factor_arrays &&
       status.has_dual_steepest_edge_weights &&
       status.has_invert;
  if (!ok) {
    if (!status.has_basis)
      highsLogUser(options.log_options, HighsLogType::kError,
                      "Not OK to solve since status.has_basis =
%" HIGHSINT_FORMAT "\n", status.has_basis); if
(!status.has_matrix) highsLogUser(options.log_options,
HighsLogType::kError, "Not OK to solve since status.has_matrix =
%" HIGHSINT_FORMAT "\n", status.has_matrix);
    //    if (!status.has_factor_arrays)
    //      highsLogUser(options.log_options, HighsLogType::kError,
    //                  "Not OK to solve since
    //      status.has_factor_arrays = %" HIGHSINT_FORMAT "\n",
    //             status.has_factor_arrays);
    if (!status.has_dual_steepest_edge_weights)
      highsLogUser(options.log_options, HighsLogType::kError,
                      "Not OK to solve since \n"
                      "status.has_dual_steepest_edge_weights = %"
HIGHSINT_FORMAT "", status.has_dual_steepest_edge_weights); if
(!status.has_invert) highsLogUser(options.log_options,
HighsLogType::kError, "Not OK to solve since status.has_invert =
%" HIGHSINT_FORMAT "\n", status.has_invert);
  }
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCostly)
    return return_status;
  // Basis and data check
  if (debugBasisConsistent(highs_model_object.options_, lp,
                           highs_model_object.basis_) ==
      HighsDebugStatus::kLogicalError)
    return HighsDebugStatus::kLogicalError;
  if (!debugWorkArraysOk(highs_model_object, phase))
    return HighsDebugStatus::kLogicalError;
  const HighsInt numTot = lp.numCol_ + lp.numRow_;
  for (HighsInt var = 0; var < numTot; ++var) {
    if (basis.nonbasicFlag_[var]) {
      // Nonbasic variable
      if (!debugOneNonbasicMoveVsWorkArraysOk(highs_model_object, var))
        return HighsDebugStatus::kLogicalError;
    }
  }
  return return_status;
}

// Methods below are not called externally

bool debugWorkArraysOk(const HighsModelObject& highs_model_object,
                       const HighsInt phase) {
  const HighsLp& lp = highs_model_object.lp_;
  const HighsSimplexInfo& info = highs_model_object.info_;
  const HighsOptions& options = highs_model_object.options_;
  bool ok = true;
  // Only check phase 2 bounds: others will have been set by solve() so can be
  // trusted
  if (phase == 2) {
    for (HighsInt col = 0; col < lp.numCol_; ++col) {
      HighsInt var = col;
      if (!highs_isInfinity(-info.workLower_[var])) {
        ok = info.workLower_[var] == lp.colLower_[col];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
              "For col %" HIGHSINT_FORMAT ", info.workLower_ should be
%g but is %g\n", col, lp.colLower_[col], info.workLower_[var]);
return ok;
        }
      }
      if (!highs_isInfinity(info.workUpper_[var])) {
        ok = info.workUpper_[var] == lp.colUpper_[col];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
              "For col %" HIGHSINT_FORMAT ", info.workUpper_ should be
%g but is %g\n", col, lp.colUpper_[col], info.workUpper_[var]);
return ok;
        }
      }
    }
    for (HighsInt row = 0; row < lp.numRow_; ++row) {
      HighsInt var = lp.numCol_ + row;
      if (!highs_isInfinity(-info.workLower_[var])) {
        ok = info.workLower_[var] == -lp.rowUpper_[row];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
              "For row %" HIGHSINT_FORMAT ", info.workLower_ should be
%g but is %g\n", row, -lp.rowUpper_[row], info.workLower_[var]);
return ok;
        }
      }
      if (!highs_isInfinity(info.workUpper_[var])) {
        ok = info.workUpper_[var] == -lp.rowLower_[row];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
              "For row %" HIGHSINT_FORMAT ", info.workUpper_ should be
%g but is %g\n", row, -lp.rowLower_[row], info.workUpper_[var]);
return ok;
        }
      }
    }
  }
  const HighsInt numTot = lp.numCol_ + lp.numRow_;
  for (HighsInt var = 0; var < numTot; ++var) {
    ok = info.workRange_[var] ==
         (info.workUpper_[var] - info.workLower_[var]);
    if (!ok) {
      highsLogUser(options.log_options, HighsLogType::kError,
          "For variable %" HIGHSINT_FORMAT ", info.workRange_ should be
%g = %g - %g \n" "but is %g", var, info.workUpper_[var] -
info.workLower_[var], info.workUpper_[var],
info.workLower_[var], info.workRange_[var]); return ok;
    }
  }
  // Don't check perturbed costs: these will have been set by solve() so can be
  // trusted
  if (!info.costs_perturbed) {
    for (HighsInt col = 0; col < lp.numCol_; ++col) {
      HighsInt var = col;
      ok = info.workCost_[var] ==
           (HighsInt)lp.sense_ * lp.colCost_[col];
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
            "For col %" HIGHSINT_FORMAT ", info.workLower_ should be %g
but is %g\n", col, lp.colLower_[col], info.workCost_[var]);
        return ok;
      }
    }
    for (HighsInt row = 0; row < lp.numRow_; ++row) {
      HighsInt var = lp.numCol_ + row;
      ok = info.workCost_[var] == 0.;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
            "For row %" HIGHSINT_FORMAT ", info.workCost_ should be zero
but is %g\n", row, info.workCost_[var]); return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool debugOneNonbasicMoveVsWorkArraysOk(
    const HighsModelObject& highs_model_object, const HighsInt var) {
  const HighsLp& lp = highs_model_object.lp_;
  const HighsSimplexInfo& info = highs_model_object.info_;
  const SimplexBasis& basis = highs_model_object.basis_;
  const HighsOptions& options = highs_model_object.options_;
  assert(var >= 0);
  assert(var < lp.numCol_ + lp.numRow_);
  // Make sure we're not checking a basic variable
  if (!basis.nonbasicFlag_[var]) return true;
  bool ok;
  if (!highs_isInfinity(-info.workLower_[var])) {
    if (!highs_isInfinity(info.workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (info.workLower_[var] == info.workUpper_[var]) {
        // Fixed variable
        ok = basis.nonbasicMove_[var] == kNonbasicMoveZe;
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
              "Fixed variable %" HIGHSINT_FORMAT " (lp.numCol_ = %"
HIGHSINT_FORMAT ") [%11g, %11g, \n"
              "%11g] so nonbasic "
              "move should be zero but is %" HIGHSINT_FORMAT "",
              var, lp.numCol_, info.workLower_[var],
              info.workValue_[var], info.workUpper_[var],
              basis.nonbasicMove_[var]);
          return ok;
        }
        ok = info.workValue_[var] == info.workLower_[var];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                          "Fixed variable %" HIGHSINT_FORMAT "
(lp.numCol_ = %" HIGHSINT_FORMAT ") so \n" "info.work value
should be %g but " "is %g", var, lp.numCol_,
info.workLower_[var], info.workValue_[var]); return ok;
        }
      } else {
        // Boxed variable
        ok = (basis.nonbasicMove_[var] == kNonbasicMoveUp) ||
             (basis.nonbasicMove_[var] == kNonbasicMoveDn);
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
              "Boxed variable %" HIGHSINT_FORMAT " (lp.numCol_ = %"
HIGHSINT_FORMAT ") [%11g, %11g, \n"
              "%11g] range %g so "
              "nonbasic move should be up/down but is  %" HIGHSINT_FORMAT "",
              var, lp.numCol_, info.workLower_[var],
              info.workValue_[var], info.workUpper_[var],
              info.workUpper_[var] - info.workLower_[var],
              basis.nonbasicMove_[var]);
          return ok;
        }
        if (basis.nonbasicMove_[var] == kNonbasicMoveUp) {
          ok = info.workValue_[var] == info.workLower_[var];
          if (!ok) {
            highsLogUser(options.log_options, HighsLogType::kError,
                            "Boxed variable %" HIGHSINT_FORMAT "
(lp.numCol_ = %" HIGHSINT_FORMAT ") with \n" "kNonbasicMoveUp so work "
"value should be %g but is %g", var, lp.numCol_,
info.workLower_[var], info.workValue_[var]); return ok;
          }
        } else {
          ok = info.workValue_[var] == info.workUpper_[var];
          if (!ok) {
            highsLogUser(options.log_options, HighsLogType::kError,
                            "Boxed variable %" HIGHSINT_FORMAT "
(lp.numCol_ = %" HIGHSINT_FORMAT ") with \n" "kNonbasicMoveDn so work "
"value should be %g but is %g", var, lp.numCol_,
info.workUpper_[var], info.workValue_[var]); return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = basis.nonbasicMove_[var] == kNonbasicMoveUp;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
            "Finite lower bound and infinite upper bound variable %"
HIGHSINT_FORMAT " \n"
            "(lp.numCol_ = "
            "%" HIGHSINT_FORMAT ") [%11g, %11g, %11g] so nonbasic move should be
up=%2" HIGHSINT_FORMAT " but is  "
            "%" HIGHSINT_FORMAT "",
            var, lp.numCol_, info.workLower_[var],
            info.workValue_[var], info.workUpper_[var],
            kNonbasicMoveUp, basis.nonbasicMove_[var]);
        return ok;
      }
      ok = info.workValue_[var] == info.workLower_[var];
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
            "Finite lower bound and infinite upper bound variable %"
HIGHSINT_FORMAT " \n"
            "(lp.numCol_ = "
            "%" HIGHSINT_FORMAT ") so work value should be %g but is %g",
            var, lp.numCol_, info.workLower_[var],
            info.workValue_[var]);
        return ok;
      }
    }
  } else {
    // Infinite lower bound
    if (!highs_isInfinity(info.workUpper_[var])) {
      ok = basis.nonbasicMove_[var] == kNonbasicMoveDn;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
            "Finite upper bound and infinite lower bound variable %"
HIGHSINT_FORMAT " \n"
            "(lp.numCol_ = "
            "%" HIGHSINT_FORMAT ") [%11g, %11g, %11g] so nonbasic move should be
down but is  "
            "%" HIGHSINT_FORMAT "",
            var, lp.numCol_, info.workLower_[var],
            info.workValue_[var], info.workUpper_[var],
            basis.nonbasicMove_[var]);
        return ok;
      }
      ok = info.workValue_[var] == info.workUpper_[var];
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
            "Finite upper bound and infinite lower bound variable %"
HIGHSINT_FORMAT " \n"
            "(lp.numCol_ = "
            "%" HIGHSINT_FORMAT ") so work value should be %g but is %g",
            var, lp.numCol_, info.workUpper_[var],
            info.workValue_[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = basis.nonbasicMove_[var] == kNonbasicMoveZe;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
            "Free variable %" HIGHSINT_FORMAT " (lp.numCol_ = %"
HIGHSINT_FORMAT ") [%11g, %11g, %11g] \n" "so nonbasic " "move should be zero
but is  %" HIGHSINT_FORMAT "", var, lp.numCol_,
info.workLower_[var], info.workValue_[var],
info.workUpper_[var], basis.nonbasicMove_[var]); return ok;
      }
      ok = info.workValue_[var] == 0.0;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
            "Free variable %" HIGHSINT_FORMAT " (lp.numCol_ = %"
HIGHSINT_FORMAT ") so work value should \n" "be zero but " "is %g", var,
lp.numCol_, info.workValue_[var]); return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool debugAllNonbasicMoveVsWorkArraysOk(
    const HighsModelObject& highs_model_object) {
  const HighsLp& lp = highs_model_object.lp_;
  //    HighsSimplexInfo &info = highs_model_object.info_;
  const SimplexBasis& basis = highs_model_object.basis_;
  const HighsOptions& options = highs_model_object.options_;
  bool ok;
  const HighsInt numTot = lp.numCol_ + lp.numRow_;
  for (HighsInt var = 0; var < numTot; ++var) {
    highsLogUser(options.log_options, HighsLogType::kError,
        "NonbasicMoveVsWorkArrays: var = %2" HIGHSINT_FORMAT ";
basis.nonbasicFlag_[var] \n"
        "= %2" HIGHSINT_FORMAT "",
        var, basis.nonbasicFlag_[var]);
    if (!basis.nonbasicFlag_[var]) continue;
    ok = debugOneNonbasicMoveVsWorkArraysOk(highs_model_object, var);
    if (!ok) {
      highsLogUser(options.log_options, HighsLogType::kError,
          "Error in NonbasicMoveVsWorkArrays for nonbasic variable %"
HIGHSINT_FORMAT "\n", var); assert(ok); return ok;
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

void debugReportReinvertOnNumericalTrouble(
    const std::string method_name, const HighsModelObject& highs_model_object,
    const double numerical_trouble_measure, const double alpha_from_col,
    const double alpha_from_row, const double numerical_trouble_tolerance,
    const bool reinvert) {
  if (highs_model_object.options_.highs_debug_level < kHighsDebugLevelCheap)
    return;
  const double abs_alpha_from_col = fabs(alpha_from_col);
  const double abs_alpha_from_row = fabs(alpha_from_row);
  const double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  const HighsInt iteration_count = highs_model_object.iteration_counts_.simplex;
  const HighsInt update_count = highs_model_object.info_.update_count;
  const std::string model_name = highs_model_object.lp_.model_name_;

  const bool numerical_trouble =
      numerical_trouble_measure > numerical_trouble_tolerance;
  const bool near_numerical_trouble =
      10 * numerical_trouble_measure > numerical_trouble_tolerance;

  const bool wrong_sign = alpha_from_col * alpha_from_row <= 0;
  if (!near_numerical_trouble && !wrong_sign) return;
  std::string adjective;
  if (numerical_trouble) {
    adjective = "       exceeds";
  } else if (near_numerical_trouble) {
    adjective = "almost exceeds";
  } else {
    adjective = "clearly satisfies";
  }
  highsLogUser(highs_model_object.options_.log_options,
                  HighsLogType::kWarning,
                  "%s (%s) [Iter %" HIGHSINT_FORMAT "; Update %" HIGHSINT_FORMAT
"] Col: %11.4g; Row: %11.4g; Diff \n"
                  "= %11.4g: Measure %11.4g %s %11.4g",
                  method_name.c_str(), model_name.c_str(), iteration_count,
                  update_count, abs_alpha_from_col, abs_alpha_from_row,
                  abs_alpha_diff, numerical_trouble_measure, adjective.c_str(),
                  numerical_trouble_tolerance);
  if (wrong_sign) {
    highsLogUser(highs_model_object.options_.log_options,
                    HighsLogType::kWarning,
                    "   Incompatible signs for Col: %11.4g and Row: %11.4g\n",
                    alpha_from_col, alpha_from_row);
  }
  if ((numerical_trouble || wrong_sign) && !reinvert) {
    highsLogUser(highs_model_object.options_.log_options,
                    HighsLogType::kWarning,
                    "   Numerical trouble or wrong sign and not reinverting\n");
  }
}
*/
