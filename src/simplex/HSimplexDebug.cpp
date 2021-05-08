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

