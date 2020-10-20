/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HEkkDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include <string>

#include "simplex/HEkkDebug.h"

HighsDebugStatus ekkDebugBasisConsistent(const HighsOptions& options,
                                      const HighsLp& simplex_lp,
                                      const SimplexBasis& simplex_basis) {
  // Cheap analysis of a Simplex basis, checking vector sizes, numbers
  // of basic/nonbasic variables and non-repetition of basic variables
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  // Check consistency of nonbasicFlag
  if (ekkDebugNonbasicFlagConsistent(options, simplex_lp, simplex_basis) ==
      HighsDebugStatus::LOGICAL_ERROR) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "nonbasicFlag inconsistent");
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  const bool right_size =
      (int)simplex_basis.basicIndex_.size() == simplex_lp.numRow_;
  // Check consistency of basicIndex
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "basicIndex size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  // Use localNonbasicFlag so that duplicate entries in basicIndex can
  // be spotted
  vector<int> localNonbasicFlag = simplex_basis.nonbasicFlag_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    int iCol = simplex_basis.basicIndex_[iRow];
    int flag = localNonbasicFlag[iCol];
    // Indicate that this column has been found in basicIndex
    localNonbasicFlag[iCol] = -1;
    if (flag) {
      // Nonzero value for localNonbasicFlag entry means that column is either
      if (flag == NONBASIC_FLAG_TRUE) {
        // Nonbasic...
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Entry basicIndex_[%d] = %d is not basic", iRow, iCol);
      } else {
        // .. or is -1 since it has already been found in basicIndex
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Entry basicIndex_[%d] = %d is already basic", iRow,
                        iCol);
        assert(flag == -1);
      }
      assert(!flag);
      return_status = HighsDebugStatus::LOGICAL_ERROR;
    }
  }
  return return_status;
}

HighsDebugStatus ekkDebugNonbasicFlagConsistent(
    const HighsOptions& options, const HighsLp& simplex_lp,
    const SimplexBasis& simplex_basis) {
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  const bool right_size = (int)simplex_basis.nonbasicFlag_.size() == numTot;
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "nonbasicFlag size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  int num_basic_variables = 0;
  for (int var = 0; var < numTot; var++) {
    if (simplex_basis.nonbasicFlag_[var] == NONBASIC_FLAG_FALSE) {
      num_basic_variables++;
    } else {
      assert(simplex_basis.nonbasicFlag_[var] == NONBASIC_FLAG_TRUE);
    }
  }
  bool right_num_basic_variables = num_basic_variables == simplex_lp.numRow_;
  if (!right_num_basic_variables) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "nonbasicFlag has %d, not %d basic variables",
                    num_basic_variables, simplex_lp.numRow_);
    assert(right_num_basic_variables);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  return return_status;
}

HighsDebugStatus ekkDebugOkForSolve(const HEkk& ekk_instance,
				 const SimplexAlgorithm algorithm,
                                 const int phase,
				 const bool perturbed) {
  if (ekk_instance.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsDebugStatus return_status = HighsDebugStatus::OK;
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  const HighsSimplexLpStatus& simplex_lp_status =
      ekk_instance.simplex_lp_status_;
  const SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;
  const HighsOptions& options = ekk_instance.options_;
  bool ok;
  // Minimal check - just look at flags. This means we trust them!
  ok = simplex_lp_status.has_basis && simplex_lp_status.has_matrix_col_wise &&
       simplex_lp_status.has_matrix_row_wise &&
       simplex_lp_status.has_factor_arrays &&
       simplex_lp_status.has_dual_steepest_edge_weights &&
       simplex_lp_status.has_invert;
  if (!ok) {
    if (!simplex_lp_status.has_basis)
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Not OK to solve since simplex_lp_status.has_basis = %d",
                      simplex_lp_status.has_basis);
    if (!simplex_lp_status.has_matrix_col_wise)
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "Not OK to solve since simplex_lp_status.has_matrix_col_wise "
          "= %d",
          simplex_lp_status.has_matrix_col_wise);
    if (!simplex_lp_status.has_matrix_row_wise)
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "Not OK to solve since simplex_lp_status.has_matrix_row_wise "
          "= %d",
          simplex_lp_status.has_matrix_row_wise);
    //    if (!simplex_lp_status.has_factor_arrays)
    //      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
    //                  "Not OK to solve since
    //      simplex_lp_status.has_factor_arrays = %d",
    //             simplex_lp_status.has_factor_arrays);
    if (!simplex_lp_status.has_dual_steepest_edge_weights)
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Not OK to solve since "
                      "simplex_lp_status.has_dual_steepest_edge_weights = %d",
                      simplex_lp_status.has_dual_steepest_edge_weights);
    if (!simplex_lp_status.has_invert)
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Not OK to solve since simplex_lp_status.has_invert = %d",
                      simplex_lp_status.has_invert);
  }
  if (ekk_instance.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return return_status;
  // Basis and data check
  if (ekkDebugBasisConsistent(ekk_instance.options_, simplex_lp,
                           ekk_instance.simplex_basis_) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return HighsDebugStatus::LOGICAL_ERROR;
  if (!ekkDebugWorkArraysOk(ekk_instance, phase))
    return HighsDebugStatus::LOGICAL_ERROR;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; ++var) {
    if (simplex_basis.nonbasicFlag_[var]) {
      // Nonbasic variable
      if (!ekkDebugOneNonbasicMoveVsWorkArraysOk(ekk_instance, var))
        return HighsDebugStatus::LOGICAL_ERROR;
    }
  }
  return return_status;
}

// Methods below are not called externally

bool ekkDebugWorkArraysOk(const HEkk& ekk_instance,
                       const int phase) {
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  const HighsSimplexInfo& simplex_info = ekk_instance.simplex_info_;
  const HighsOptions& options = ekk_instance.options_;
  bool ok = true;
  // Only check phase 2 bounds: others will have been set by solve() so can be
  // trusted
  if (phase == 2) {
    for (int col = 0; col < simplex_lp.numCol_; ++col) {
      int var = col;
      if (!highs_isInfinity(-simplex_info.workLower_[var])) {
        ok = simplex_info.workLower_[var] == simplex_lp.colLower_[col];
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "For col %d, simplex_info.workLower_ should be %g but is %g", col,
              simplex_lp.colLower_[col], simplex_info.workLower_[var]);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info.workUpper_[var])) {
        ok = simplex_info.workUpper_[var] == simplex_lp.colUpper_[col];
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "For col %d, simplex_info.workUpper_ should be %g but is %g", col,
              simplex_lp.colUpper_[col], simplex_info.workUpper_[var]);
          return ok;
        }
      }
    }
    for (int row = 0; row < simplex_lp.numRow_; ++row) {
      int var = simplex_lp.numCol_ + row;
      if (!highs_isInfinity(-simplex_info.workLower_[var])) {
        ok = simplex_info.workLower_[var] == -simplex_lp.rowUpper_[row];
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "For row %d, simplex_info.workLower_ should be %g but is %g", row,
              -simplex_lp.rowUpper_[row], simplex_info.workLower_[var]);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info.workUpper_[var])) {
        ok = simplex_info.workUpper_[var] == -simplex_lp.rowLower_[row];
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "For row %d, simplex_info.workUpper_ should be %g but is %g", row,
              -simplex_lp.rowLower_[row], simplex_info.workUpper_[var]);
          return ok;
        }
      }
    }
  }
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; ++var) {
    ok = simplex_info.workRange_[var] ==
         (simplex_info.workUpper_[var] - simplex_info.workLower_[var]);
    if (!ok) {
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "For variable %d, simplex_info.workRange_ should be %g = %g - %g "
          "but is %g",
          var, simplex_info.workUpper_[var] - simplex_info.workLower_[var],
          simplex_info.workUpper_[var], simplex_info.workLower_[var],
          simplex_info.workRange_[var]);
      return ok;
    }
  }
  // Don't check perturbed costs: these will have been set by solve() so can be
  // trusted
  if (!simplex_info.costs_perturbed) {
    for (int col = 0; col < simplex_lp.numCol_; ++col) {
      int var = col;
      ok = simplex_info.workCost_[var] ==
           (int)simplex_lp.sense_ * simplex_lp.colCost_[col];
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "For col %d, simplex_info.workLower_ should be %g but is %g", col,
            simplex_lp.colLower_[col], simplex_info.workCost_[var]);
        return ok;
      }
    }
    for (int row = 0; row < simplex_lp.numRow_; ++row) {
      int var = simplex_lp.numCol_ + row;
      ok = simplex_info.workCost_[var] == 0.;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "For row %d, simplex_info.workCost_ should be zero but is %g", row,
            simplex_info.workCost_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool ekkDebugOneNonbasicMoveVsWorkArraysOk(
    const HEkk& ekk_instance, const int var) {
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  const HighsSimplexInfo& simplex_info = ekk_instance.simplex_info_;
  const SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;
  const HighsOptions& options = ekk_instance.options_;
  assert(var >= 0);
  assert(var < simplex_lp.numCol_ + simplex_lp.numRow_);
  // Make sure we're not checking a basic variable
  if (!simplex_basis.nonbasicFlag_[var]) return true;
  bool ok;
  if (!highs_isInfinity(-simplex_info.workLower_[var])) {
    if (!highs_isInfinity(simplex_info.workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (simplex_info.workLower_[var] == simplex_info.workUpper_[var]) {
        // Fixed variable
        ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "Fixed variable %d (simplex_lp.numCol_ = %d) [%11g, %11g, "
              "%11g] so nonbasic "
              "move should be zero but is %d",
              var, simplex_lp.numCol_, simplex_info.workLower_[var],
              simplex_info.workValue_[var], simplex_info.workUpper_[var],
              simplex_basis.nonbasicMove_[var]);
          return ok;
        }
        ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
        if (!ok) {
          HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                          "Fixed variable %d (simplex_lp.numCol_ = %d) so "
                          "simplex_info.work value should be %g but "
                          "is %g",
                          var, simplex_lp.numCol_, simplex_info.workLower_[var],
                          simplex_info.workValue_[var]);
          return ok;
        }
      } else {
        // Boxed variable
        ok = (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) ||
             (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_DN);
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "Boxed variable %d (simplex_lp.numCol_ = %d) [%11g, %11g, "
              "%11g] range %g so "
              "nonbasic move should be up/down but is  %d",
              var, simplex_lp.numCol_, simplex_info.workLower_[var],
              simplex_info.workValue_[var], simplex_info.workUpper_[var],
              simplex_info.workUpper_[var] - simplex_info.workLower_[var],
              simplex_basis.nonbasicMove_[var]);
          return ok;
        }
        if (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
          ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
          if (!ok) {
            HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                            "Boxed variable %d (simplex_lp.numCol_ = %d) with "
                            "NONBASIC_MOVE_UP so work "
                            "value should be %g but is %g",
                            var, simplex_lp.numCol_,
                            simplex_info.workLower_[var],
                            simplex_info.workValue_[var]);
            return ok;
          }
        } else {
          ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
          if (!ok) {
            HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                            "Boxed variable %d (simplex_lp.numCol_ = %d) with "
                            "NONBASIC_MOVE_DN so work "
                            "value should be %g but is %g",
                            var, simplex_lp.numCol_,
                            simplex_info.workUpper_[var],
                            simplex_info.workValue_[var]);
            return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Finite lower bound and infinite upper bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be up=%2d but is  "
            "%d",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            NONBASIC_MOVE_UP, simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Finite lower bound and infinite upper bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) so work value should be %g but is %g",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var]);
        return ok;
      }
    }
  } else {
    // Infinite lower bound
    if (!highs_isInfinity(simplex_info.workUpper_[var])) {
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_DN;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Finite upper bound and infinite lower bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be down but is  "
            "%d",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Finite upper bound and infinite lower bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) so work value should be %g but is %g",
            var, simplex_lp.numCol_, simplex_info.workUpper_[var],
            simplex_info.workValue_[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Free variable %d (simplex_lp.numCol_ = %d) [%11g, %11g, %11g] "
            "so nonbasic "
            "move should be zero but is  %d",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == 0.0;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Free variable %d (simplex_lp.numCol_ = %d) so work value should "
            "be zero but "
            "is %g",
            var, simplex_lp.numCol_, simplex_info.workValue_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool ekkDebugAllNonbasicMoveVsWorkArraysOk(
    const HEkk& ekk_instance) {
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  //    HighsSimplexInfo &simplex_info = ekk_instance.simplex_info_;
  const SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;
  const HighsOptions& options = ekk_instance.options_;
  bool ok;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; ++var) {
    HighsLogMessage(
        options.logfile, HighsMessageType::ERROR,
        "NonbasicMoveVsWorkArrays: var = %2d; simplex_basis.nonbasicFlag_[var] "
        "= %2d",
        var, simplex_basis.nonbasicFlag_[var]);
    if (!simplex_basis.nonbasicFlag_[var]) continue;
    ok = ekkDebugOneNonbasicMoveVsWorkArraysOk(ekk_instance, var);
    if (!ok) {
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "Error in NonbasicMoveVsWorkArrays for nonbasic variable %d", var);
      assert(ok);
      return ok;
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

void ekkDebugReportReinvertOnNumericalTrouble(
    const std::string method_name, const HEkk& ekk_instance,
    const double numerical_trouble_measure, const double alpha_from_col,
    const double alpha_from_row, const double numerical_trouble_tolerance,
    const bool reinvert) {
  if (ekk_instance.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return;
  const double abs_alpha_from_col = fabs(alpha_from_col);
  const double abs_alpha_from_row = fabs(alpha_from_row);
  const double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  const int iteration_count = ekk_instance.iteration_count_;
  const int update_count = ekk_instance.simplex_info_.update_count;
  const std::string model_name = ekk_instance.simplex_lp_.model_name_;

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
  HighsLogMessage(ekk_instance.options_.logfile,
                  HighsMessageType::WARNING,
                  "%s (%s) [Iter %d; Update %d] Col: %11.4g; Row: %11.4g; Diff "
                  "= %11.4g: Measure %11.4g %s %11.4g",
                  method_name.c_str(), model_name.c_str(), iteration_count,
                  update_count, abs_alpha_from_col, abs_alpha_from_row,
                  abs_alpha_diff, numerical_trouble_measure, adjective.c_str(),
                  numerical_trouble_tolerance);
  if (wrong_sign) {
    HighsLogMessage(ekk_instance.options_.logfile,
                    HighsMessageType::WARNING,
                    "   Incompatible signs for Col: %11.4g and Row: %11.4g",
                    alpha_from_col, alpha_from_row);
  }
  if ((numerical_trouble || wrong_sign) && !reinvert) {
    HighsLogMessage(ekk_instance.options_.logfile,
                    HighsMessageType::WARNING,
                    "   Numerical trouble or wrong sign and not reinverting");
  }
}
