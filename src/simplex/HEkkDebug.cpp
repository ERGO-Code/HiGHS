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
/**@file lp_data/HEkkDebug.cpp
 * @brief
 */

#include "simplex/HEkkDebug.h"

#include <math.h>

#include <cassert>
#include <string>

#include "lp_data/HighsDebug.h"
#include "lp_data/HighsModelUtils.h"

using std::fabs;
using std::max;

const double ok_feasibility_difference = 1e-3;

const double large_basic_dual = 1e-12;
const double excessive_basic_dual = sqrt(large_basic_dual);

const double large_residual_error = 1e-12;
const double excessive_residual_error = sqrt(large_residual_error);

const double updated_dual_small_relative_error = 1e-12;
const double updated_dual_large_relative_error =
    sqrt(updated_dual_small_relative_error);
const double updated_dual_small_absolute_error = 1e-6;
const double updated_dual_large_absolute_error =
    sqrt(updated_dual_small_absolute_error);

HighsDebugStatus ekkDebugSimplex(const std::string message,
                                 const HEkk& ekk_instance,
                                 const SimplexAlgorithm algorithm,
                                 const HighsInt phase, const bool initialise) {
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  static double max_max_basic_dual;
  static double max_max_primal_residual;
  static double max_max_dual_residual;
  if (initialise) {
    max_max_basic_dual = 0;
    max_max_primal_residual = 0;
    max_max_dual_residual = 0;
    return ::HighsDebugStatus::kOk;
  }
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  const HighsSimplexInfo& simplex_info = ekk_instance.simplex_info_;
  const SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;
  const HighsOptions& options = ekk_instance.options_;

  const HighsInt num_col = simplex_lp.numCol_;
  const HighsInt num_row = simplex_lp.numRow_;
  const HighsInt num_tot = num_col + num_row;
  const HighsInt iteration_count = ekk_instance.iteration_count_;
  std::string value_adjective;
  HighsLogType report_level;

  // Check the nonbasic flags are all NONBASIC_FLAG_TRUE or NONBASIC_FLAG_FALSE
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    HighsInt flag = simplex_basis.nonbasicFlag_[iVar];
    bool flag_error = flag != NONBASIC_FLAG_TRUE && flag != NONBASIC_FLAG_FALSE;
    if (flag_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Variable %" HIGHSINT_FORMAT
                   " has "
                   "nonbasic flag = %" HIGHSINT_FORMAT "\n",
                   message.c_str(), iteration_count, iVar, flag);
      assert(!flag_error);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  const double dual_feasibility_tolerance = options.dual_feasibility_tolerance;
  HighsInt num_dual_infeasibility = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibility = 0;
  HighsInt num_primal_infeasibility = 0;
  double max_primal_infeasibility = 0;
  double sum_primal_infeasibility = 0;
  // Check the nonbasic variables
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE) continue;
    // For nonbasic variables, check that they are on a bound (or free
    // at 0 with correct nonbasic move. Determine dual infeasibilities
    double dual = simplex_info.workDual_[iVar];
    double lower = simplex_info.workLower_[iVar];
    double upper = simplex_info.workUpper_[iVar];
    double value = simplex_info.workValue_[iVar];
    double primal_error = 0;
    double dual_infeasibility = 0;
    HighsInt move;
    if (lower == upper) {
      primal_error = fabs(lower - value);
      move = NONBASIC_MOVE_ZE;
    } else if (value == lower) {
      move = NONBASIC_MOVE_UP;
      dual_infeasibility = max(-dual, 0.);
    } else if (value == upper) {
      move = NONBASIC_MOVE_DN;
      dual_infeasibility = max(dual, 0.);
    } else {
      // If not fixed or at a bound, only valid status can be zero and free
      primal_error = fabs(value);
      move = NONBASIC_MOVE_ZE;
      dual_infeasibility = fabs(dual);
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility > dual_feasibility_tolerance) {
        num_dual_infeasibility++;
      }
      max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
    if (primal_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Nonbasic variable %" HIGHSINT_FORMAT
                   " "
                   "has primal error "
                   "= %g for [%g, %g, %g]\n",
                   message.c_str(), iteration_count, iVar, primal_error, lower,
                   value, upper);
      assert(!primal_error);
      return HighsDebugStatus::kLogicalError;
    }
    bool move_error = move != simplex_basis.nonbasicMove_[iVar];
    if (move_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Nonbasic variable %" HIGHSINT_FORMAT
                   " "
                   "has move error "
                   "[%" HIGHSINT_FORMAT " <> %" HIGHSINT_FORMAT
                   "] for [%g, %g, %g]\n",
                   message.c_str(), iteration_count, iVar, move,
                   simplex_basis.nonbasicMove_[iVar], lower, value, upper);
      assert(!move_error);
      return HighsDebugStatus::kLogicalError;
    }
  }
  // Check the basic variables
  double max_basic_dual = 0;
  const double base =
      simplex_info.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = simplex_basis.basicIndex_[iRow];
    // For basic variables, check that the nonbasic flag isn't set,
    // that nonbasicMove is zero, that baseLower/Upper are correct,
    // that the dual is zero and, in primal phase 1, that the cost is
    // correct. Determine primal infeasibilities
    bool nonbasicFlag_error =
        simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE;
    if (nonbasicFlag_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Basic variable %" HIGHSINT_FORMAT
                   " "
                   "has nonbasicFlag = %" HIGHSINT_FORMAT "\n",
                   message.c_str(), iteration_count, iVar,
                   simplex_basis.nonbasicFlag_[iVar]);
      assert(!nonbasicFlag_error);
      return HighsDebugStatus::kLogicalError;
    }
    bool nonbasicMove_error = simplex_basis.nonbasicMove_[iVar];
    if (nonbasicMove_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Basic variable %" HIGHSINT_FORMAT
                   " "
                   "has nonbasicMove = %" HIGHSINT_FORMAT "\n",
                   message.c_str(), iteration_count, iVar,
                   simplex_basis.nonbasicMove_[iVar]);
      assert(!nonbasicMove_error);
      return HighsDebugStatus::kLogicalError;
    }
    double workLower = simplex_info.workLower_[iVar];
    double workUpper = simplex_info.workUpper_[iVar];
    double cost = simplex_info.workCost_[iVar];
    double dual = simplex_info.workDual_[iVar];
    double lower = simplex_info.baseLower_[iRow];
    double upper = simplex_info.baseUpper_[iRow];
    double value = simplex_info.baseValue_[iRow];
    bool baseBound_error = workLower != lower || workUpper != upper;
    if (baseBound_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Basic variable %" HIGHSINT_FORMAT
                   " "
                   "(in row %" HIGHSINT_FORMAT
                   ") has "
                   "baseBound [%g, %g] and workBound [%g, %g]\n",
                   message.c_str(), iteration_count, iVar, iRow, lower, upper,
                   workLower, workUpper);
      assert(!baseBound_error);
      return HighsDebugStatus::kLogicalError;
    }
    max_basic_dual = max(fabs(dual), max_basic_dual);
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1;
    }
    if (algorithm == SimplexAlgorithm::kPrimal && phase == 1) {
      double primal_phase1_cost = bound_violated;
      if (base)
        primal_phase1_cost *= 1 + base * simplex_info.numTotRandomValue_[iRow];
      bool primal_phase1_cost_error = fabs(cost - primal_phase1_cost);
      if (primal_phase1_cost_error) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                     " Basic variable %" HIGHSINT_FORMAT
                     " "
                     "(in row %" HIGHSINT_FORMAT
                     ") has "
                     "primal phase 1 cost %g for [%g, %g, %g]\n",
                     message.c_str(), iteration_count, iVar, iRow, cost, lower,
                     value, upper);
        assert(!primal_phase1_cost_error);
        return HighsDebugStatus::kLogicalError;
      }
    }
    if (!bound_violated) continue;
    // @primal_infeasibility calculation
    double primal_infeasibility = 0;
    if (bound_violated < 0) {
      primal_infeasibility = lower - value;
    } else {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibility++;
    max_primal_infeasibility =
        max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibility += primal_infeasibility;
  }
  // Report on basic dual values
  if (max_basic_dual > excessive_basic_dual) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kInfo;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
  } else if (max_basic_dual > large_basic_dual) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = debugWorseStatus(HighsDebugStatus::kWarning, return_status);
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = debugWorseStatus(HighsDebugStatus::kOk, return_status);
  }
  if (max_basic_dual > 2 * max_max_basic_dual) {
    highsLogDev(options.log_options, report_level,
                "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                " %-9s max   basic dual = %9.4g\n",
                message.c_str(), iteration_count, value_adjective.c_str(),
                max_basic_dual);
    max_max_basic_dual = max_basic_dual;
  }
  // Check that the number, max and sums of primal and dual infeasibilities (if
  // known) are correct
  const HighsInt info_num_primal_infeasibility =
      ekk_instance.simplex_info_.num_primal_infeasibility;
  if (info_num_primal_infeasibility >= 0) {
    const bool illegal_num_primal_infeasibility =
        num_primal_infeasibility != info_num_primal_infeasibility;
    if (illegal_num_primal_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %" HIGHSINT_FORMAT
                   " not "
                   "%" HIGHSINT_FORMAT " primal infeasibilities\n",
                   message.c_str(), iteration_count, num_primal_infeasibility,
                   info_num_primal_infeasibility);
      assert(!illegal_num_primal_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double info_max_primal_infeasibility =
      ekk_instance.simplex_info_.max_primal_infeasibility;
  if (info_max_primal_infeasibility >= 0) {
    const bool illegal_max_primal_infeasibility =
        fabs(max_primal_infeasibility - info_max_primal_infeasibility) >
        ok_feasibility_difference;
    if (illegal_max_primal_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %g not "
                   "%g max primal infeasibility\n",
                   message.c_str(), iteration_count, max_primal_infeasibility,
                   info_max_primal_infeasibility);
      assert(!illegal_max_primal_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double info_sum_primal_infeasibility =
      ekk_instance.simplex_info_.sum_primal_infeasibility;
  if (info_sum_primal_infeasibility >= 0) {
    const bool illegal_sum_primal_infeasibility =
        fabs(sum_primal_infeasibility - info_sum_primal_infeasibility) >
        ok_feasibility_difference;
    if (illegal_sum_primal_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %g not "
                   "%g sum primal infeasibilities\n",
                   message.c_str(), iteration_count, sum_primal_infeasibility,
                   info_sum_primal_infeasibility);
      assert(!illegal_sum_primal_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const HighsInt info_num_dual_infeasibility =
      ekk_instance.simplex_info_.num_dual_infeasibility;
  if (info_num_dual_infeasibility >= 0) {
    const bool illegal_num_dual_infeasibility =
        num_dual_infeasibility != info_num_dual_infeasibility;
    if (illegal_num_dual_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %" HIGHSINT_FORMAT
                   " not "
                   "%" HIGHSINT_FORMAT " dual infeasibilities\n",
                   message.c_str(), iteration_count, num_dual_infeasibility,
                   info_num_dual_infeasibility);
      assert(!illegal_num_dual_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double info_max_dual_infeasibility =
      ekk_instance.simplex_info_.max_dual_infeasibility;
  if (info_max_dual_infeasibility >= 0) {
    const bool illegal_max_dual_infeasibility =
        fabs(max_dual_infeasibility - info_max_dual_infeasibility) >
        ok_feasibility_difference;
    if (illegal_max_dual_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %g not "
                   "%g max dual infeasibility\n",
                   message.c_str(), iteration_count, max_dual_infeasibility,
                   info_max_dual_infeasibility);
      assert(!illegal_max_dual_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double info_sum_dual_infeasibility =
      ekk_instance.simplex_info_.sum_dual_infeasibility;
  if (info_sum_dual_infeasibility >= 0) {
    const bool illegal_sum_dual_infeasibility =
        fabs(sum_dual_infeasibility - info_sum_dual_infeasibility) >
        ok_feasibility_difference;
    if (illegal_sum_dual_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %g not "
                   "%g sum dual infeasibilities\n",
                   message.c_str(), iteration_count, sum_dual_infeasibility,
                   info_sum_dual_infeasibility);
      assert(!illegal_sum_dual_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  // Check any assumed feasibility
  bool require_primal_feasible_in_primal_simplex =
      algorithm == SimplexAlgorithm::kPrimal && (phase == 0 || phase == 2);
  bool require_primal_feasible_in_dual_simplex =
      algorithm == SimplexAlgorithm::kDual && phase == 0;
  bool require_primal_feasible = require_primal_feasible_in_primal_simplex ||
                                 require_primal_feasible_in_dual_simplex;
  bool illegal_primal_infeasibility =
      require_primal_feasible && num_primal_infeasibility > 0;
  if (illegal_primal_infeasibility) {
    // Should be primal feasible, but isn't
    highsLogUser(options.log_options, HighsLogType::kError,
                 "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                 " Should be primal "
                 "feasible, but num / max "
                 "/ sum primal infeasibility is %" HIGHSINT_FORMAT
                 " / %g / %g\n",
                 message.c_str(), iteration_count, num_primal_infeasibility,
                 max_primal_infeasibility, sum_primal_infeasibility);
    assert(!illegal_primal_infeasibility);
    return HighsDebugStatus::kLogicalError;
  }
  bool require_dual_feasible_in_dual_simplex =
      algorithm == SimplexAlgorithm::kDual &&
      ekk_instance.simplex_lp_status_.has_fresh_rebuild &&
      ekk_instance.simplex_info_.allow_cost_perturbation &&
      ekk_instance.scaled_model_status_ != HighsModelStatus::kDualInfeasible &&
      ekk_instance.scaled_model_status_ !=
          HighsModelStatus::kPrimalDualInfeasible;

  bool illegal_dual_infeasibility =
      (require_dual_feasible_in_dual_simplex || phase == 0) &&
      num_dual_infeasibility > 0;
  if (illegal_dual_infeasibility) {
    // Dual simplex or optimal but has dual infeasibilities
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
        " Should be dual "
        "feasible, but num / max / "
        "sum dual infeasibility is %" HIGHSINT_FORMAT
        " / %g / %g; Phase = %" HIGHSINT_FORMAT "; status = %s\n",
        message.c_str(), iteration_count, num_dual_infeasibility,
        max_dual_infeasibility, sum_dual_infeasibility, phase,
        utilModelStatusToString(ekk_instance.scaled_model_status_).c_str());
    assert(!illegal_dual_infeasibility);
    return HighsDebugStatus::kLogicalError;
  }

  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCostly)
    return return_status;
  // Now determine the primal and dual residuals.
  //
  // This uses the primal values for the columns to determine row
  // activities that are checked against the primal values for the
  // rows. It uses the pi vector to determine column duals. The
  // entries of the pi vector are the negated duals for nonbasic rows,
  // and costs for basic rows. The latter are normally zero, but will
  // be nonzero if the constraint is violated in primal phase 1, or if
  // the row cost is a perturbed zero in dual simplex.
  vector<double> primal_value(num_tot);
  vector<double> dual_value(num_tot);
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    primal_value[iVar] = simplex_info.workValue_[iVar];
    dual_value[iVar] = simplex_info.workDual_[iVar];
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = simplex_basis.basicIndex_[iRow];
    primal_value[iVar] = simplex_info.baseValue_[iRow];
    dual_value[iVar] = -simplex_info.workCost_[iVar];
  }
  // Accumulate primal_activities
  double max_dual_residual = 0;
  vector<double> primal_activity(num_row, 0);
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    double dual = simplex_info.workCost_[iCol];
    double value = primal_value[iCol];
    for (HighsInt iEl = simplex_lp.Astart_[iCol];
         iEl < simplex_lp.Astart_[iCol + 1]; iEl++) {
      HighsInt iRow = simplex_lp.Aindex_[iEl];
      HighsInt iVar = num_col + iRow;
      double Avalue = simplex_lp.Avalue_[iEl];
      primal_activity[iRow] += value * Avalue;
      dual += dual_value[iVar] * Avalue;
    }
    double dual_residual = fabs(dual - simplex_info.workDual_[iCol]);
    max_dual_residual = max(dual_residual, max_dual_residual);
  }
  // Remember that simplex row values are the negated row activities
  double max_primal_residual = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    double primal_residual = fabs(primal_activity[iRow] + primal_value[iVar]);
    max_primal_residual = max(primal_residual, max_primal_residual);
  }
  if (max_primal_residual > excessive_residual_error) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kInfo;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
  } else if (max_primal_residual > large_residual_error) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = debugWorseStatus(HighsDebugStatus::kWarning, return_status);
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = debugWorseStatus(HighsDebugStatus::kOk, return_status);
  }
  if (max_primal_residual > 2 * max_max_primal_residual) {
    highsLogDev(options.log_options, report_level,
                "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                " %-9s max primal residual = %9.4g\n",
                message.c_str(), iteration_count, value_adjective.c_str(),
                max_primal_residual);
    max_max_primal_residual = max_primal_residual;
  }
  if (max_dual_residual > excessive_residual_error) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kInfo;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
  } else if (max_dual_residual > large_residual_error) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = debugWorseStatus(HighsDebugStatus::kWarning, return_status);
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = debugWorseStatus(HighsDebugStatus::kOk, return_status);
  }
  if (max_dual_residual > 2 * max_max_dual_residual) {
    highsLogDev(options.log_options, report_level,
                "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                " %-9s max   dual residual = %9.4g\n",
                message.c_str(), iteration_count, value_adjective.c_str(),
                max_dual_residual);
    max_max_dual_residual = max_dual_residual;
  }
  return return_status;
}

HighsDebugStatus ekkDebugBasisCorrect(const HEkk& ekk_instance) {
  // Nontrivially expensive analysis of a Simplex basis, checking
  // consistency, and then correctness of nonbasicMove
  const HighsOptions& options = ekk_instance.options_;
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const bool consistent =
      ekkDebugBasisConsistent(ekk_instance) != HighsDebugStatus::kLogicalError;
  if (!consistent) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Supposed to be a Simplex basis, but not consistent\n");
    assert(consistent);
    return_status = HighsDebugStatus::kLogicalError;
  }
  if (options.highs_debug_level < kHighsDebugLevelCostly) return return_status;
  const bool correct_nonbasicMove =
      ekkDebugNonbasicMove(ekk_instance) != HighsDebugStatus::kLogicalError;
  if (!correct_nonbasicMove) {
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "Supposed to be a Simplex basis, but nonbasicMove is incorrect\n");
    assert(correct_nonbasicMove);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus ekkDebugNonbasicMove(const HEkk& ekk_instance) {
  // Non-trivially expensive check of NonbasicMove
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsOptions& options = ekk_instance.options_;
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  const SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;
  HighsInt num_free_variable_move_errors = 0;
  HighsInt num_lower_bounded_variable_move_errors = 0;
  HighsInt num_upper_bounded_variable_move_errors = 0;
  HighsInt num_boxed_variable_move_errors = 0;
  HighsInt num_fixed_variable_move_errors = 0;
  const HighsInt numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  bool right_size = (HighsInt)simplex_basis.nonbasicMove_.size() == numTot;
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
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic variable
    if (iVar < simplex_lp.numCol_) {
      lower = simplex_lp.colLower_[iVar];
      upper = simplex_lp.colUpper_[iVar];
    } else {
      HighsInt iRow = iVar - simplex_lp.numCol_;
      lower = -simplex_lp.rowUpper_[iRow];
      upper = -simplex_lp.rowLower_[iRow];
    }

    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free
        if (simplex_basis.nonbasicMove_[iVar]) {
          num_free_variable_move_errors++;
        }
      } else {
        // Only lower bounded
        if (simplex_basis.nonbasicMove_[iVar] != NONBASIC_MOVE_UP) {
          num_lower_bounded_variable_move_errors++;
        }
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded
        if (simplex_basis.nonbasicMove_[iVar] != NONBASIC_MOVE_DN) {
          num_upper_bounded_variable_move_errors++;
        }
      } else {
        // Boxed or fixed
        if (lower != upper) {
          // Boxed
          if (!simplex_basis.nonbasicMove_[iVar]) {
            num_boxed_variable_move_errors++;
          }
        } else {
          // Fixed
          if (simplex_basis.nonbasicMove_[iVar]) {
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
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "There are %" HIGHSINT_FORMAT " nonbasicMove errors: %" HIGHSINT_FORMAT
        " free; %" HIGHSINT_FORMAT " lower; %" HIGHSINT_FORMAT
        " upper; %" HIGHSINT_FORMAT
        " "
        "boxed; %" HIGHSINT_FORMAT " fixed\n",
        num_errors, num_free_variable_move_errors,
        num_lower_bounded_variable_move_errors,
        num_upper_bounded_variable_move_errors, num_boxed_variable_move_errors,
        num_fixed_variable_move_errors);
    assert(num_errors == 0);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus ekkDebugBasisConsistent(const HEkk& ekk_instance) {
  // Cheap analysis of a Simplex basis, checking vector sizes, numbers
  // of basic/nonbasic variables and non-repetition of basic variables
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsOptions& options = ekk_instance.options_;
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  const SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;
  // Check consistency of nonbasicFlag
  if (ekkDebugNonbasicFlagConsistent(ekk_instance) ==
      HighsDebugStatus::kLogicalError) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicFlag inconsistent\n");
    return_status = HighsDebugStatus::kLogicalError;
  }
  const bool right_size =
      (HighsInt)simplex_basis.basicIndex_.size() == simplex_lp.numRow_;
  // Check consistency of basicIndex
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "basicIndex size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  // Use localNonbasicFlag so that duplicate entries in basicIndex can
  // be spotted
  vector<int8_t> localNonbasicFlag = simplex_basis.nonbasicFlag_;
  for (HighsInt iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    HighsInt iCol = simplex_basis.basicIndex_[iRow];
    HighsInt flag = localNonbasicFlag[iCol];
    // Indicate that this column has been found in basicIndex
    localNonbasicFlag[iCol] = -1;
    if (flag) {
      // Nonzero value for localNonbasicFlag entry means that column is either
      if (flag == NONBASIC_FLAG_TRUE) {
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

HighsDebugStatus ekkDebugNonbasicFlagConsistent(const HEkk& ekk_instance) {
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsOptions& options = ekk_instance.options_;
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  const SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;
  HighsInt numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  const bool right_size =
      (HighsInt)simplex_basis.nonbasicFlag_.size() == numTot;
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicFlag size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  HighsInt num_basic_variables = 0;
  for (HighsInt var = 0; var < numTot; var++) {
    if (simplex_basis.nonbasicFlag_[var] == NONBASIC_FLAG_FALSE) {
      num_basic_variables++;
    } else {
      assert(simplex_basis.nonbasicFlag_[var] == NONBASIC_FLAG_TRUE);
    }
  }
  bool right_num_basic_variables = num_basic_variables == simplex_lp.numRow_;
  if (!right_num_basic_variables) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicFlag has %" HIGHSINT_FORMAT ", not %" HIGHSINT_FORMAT
                 " basic variables\n",
                 num_basic_variables, simplex_lp.numRow_);
    assert(right_num_basic_variables);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus ekkDebugOkForSolve(
    const HEkk& ekk_instance, const SimplexAlgorithm algorithm,
    const HighsInt phase, const HighsModelStatus scaled_model_status) {
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  const HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  const HighsSimplexLpStatus& simplex_lp_status =
      ekk_instance.simplex_lp_status_;
  const SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;
  const HighsOptions& options = ekk_instance.options_;
  bool ok;
  // Minimal check - just look at flags. This means we trust them!
  ok = simplex_lp_status.has_basis && simplex_lp_status.has_matrix &&
       simplex_lp_status.has_factor_arrays &&
       //       simplex_lp_status.has_dual_steepest_edge_weights &&
       simplex_lp_status.has_invert;
  if (!ok) {
    if (!simplex_lp_status.has_basis)
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Not OK to solve since simplex_lp_status.has_basis = "
                   "%" HIGHSINT_FORMAT "\n",
                   simplex_lp_status.has_basis);
    if (!simplex_lp_status.has_matrix)
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Not OK to solve since simplex_lp_status.has_matrix = "
                   "%" HIGHSINT_FORMAT "\n",
                   simplex_lp_status.has_matrix);
    if (!simplex_lp_status.has_factor_arrays)
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Not OK to solve since simplex_lp_status.has_factor_arrays "
                   "= %" HIGHSINT_FORMAT "\n",
                   simplex_lp_status.has_factor_arrays);
    if (!simplex_lp_status.has_dual_steepest_edge_weights)
      highsLogUser(
          options.log_options, HighsLogType::kError,
          "Not OK to solve since "
          "simplex_lp_status.has_dual_steepest_edge_weights = %" HIGHSINT_FORMAT
          "\n",
          simplex_lp_status.has_dual_steepest_edge_weights);
    if (!simplex_lp_status.has_invert)
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Not OK to solve since simplex_lp_status.has_invert = "
                   "%" HIGHSINT_FORMAT "\n",
                   simplex_lp_status.has_invert);
  }
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCostly)
    return return_status;
  // Basis and data check
  if (ekkDebugBasisConsistent(ekk_instance) == HighsDebugStatus::kLogicalError)
    return HighsDebugStatus::kLogicalError;
  // Check work cost, lower, upper and range
  if (!ekkDebugWorkArraysOk(ekk_instance, algorithm, phase,
                            scaled_model_status))
    return HighsDebugStatus::kLogicalError;
  const HighsInt numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  // Check nonbasic move against work cost, lower, upper and range
  for (HighsInt var = 0; var < numTot; ++var) {
    if (simplex_basis.nonbasicFlag_[var]) {
      // Nonbasic variable
      if (!ekkDebugOneNonbasicMoveVsWorkArraysOk(ekk_instance, var))
        return HighsDebugStatus::kLogicalError;
    }
  }
  return return_status;
}

// Methods below are not called externally

bool ekkDebugWorkArraysOk(const HEkk& ekk_instance,
                          const SimplexAlgorithm algorithm,
                          const HighsInt phase,
                          const HighsModelStatus scaled_model_status) {
  const HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  const HighsSimplexInfo& simplex_info = ekk_instance.simplex_info_;
  const HighsOptions& options = ekk_instance.options_;
  bool ok = true;
  // Don't check dual simplex phase 1 bounds or perturbed bounds
  const bool dual_phase1 = algorithm == SimplexAlgorithm::kDual && phase == 1;
  const bool primal_phase1 =
      algorithm == SimplexAlgorithm::kPrimal && phase == 1;
  if (!(dual_phase1 || simplex_info.bounds_perturbed)) {
    for (HighsInt col = 0; col < simplex_lp.numCol_; ++col) {
      HighsInt var = col;
      if (!highs_isInfinity(-simplex_info.workLower_[var])) {
        double lp_lower = simplex_info.workLower_[var];
        ok = lp_lower == simplex_lp.colLower_[col];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "For col %" HIGHSINT_FORMAT
                       ", simplex_info.workLower_ should be %g but is %g\n",
                       col, simplex_lp.colLower_[col], lp_lower);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info.workUpper_[var])) {
        double lp_upper = simplex_info.workUpper_[var];
        ok = lp_upper == simplex_lp.colUpper_[col];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "For col %" HIGHSINT_FORMAT
                       ", simplex_info.workUpper_ should be %g but is %g\n",
                       col, simplex_lp.colUpper_[col], lp_upper);
          return ok;
        }
      }
    }
    for (HighsInt row = 0; row < simplex_lp.numRow_; ++row) {
      HighsInt var = simplex_lp.numCol_ + row;
      if (!highs_isInfinity(-simplex_info.workLower_[var])) {
        double lp_lower = simplex_info.workLower_[var];
        ok = lp_lower == -simplex_lp.rowUpper_[row];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "For row %" HIGHSINT_FORMAT
                       ", simplex_info.workLower_ should be %g but is %g\n",
                       row, -simplex_lp.rowUpper_[row], lp_lower);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info.workUpper_[var])) {
        double lp_upper = simplex_info.workUpper_[var];
        ok = lp_upper == -simplex_lp.rowLower_[row];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "For row %" HIGHSINT_FORMAT
                       ", simplex_info.workUpper_ should be %g but is %g\n",
                       row, -simplex_lp.rowLower_[row], lp_upper);
          return ok;
        }
      }
    }
    const HighsInt numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
    for (HighsInt var = 0; var < numTot; ++var) {
      ok = simplex_info.workRange_[var] ==
           (simplex_info.workUpper_[var] - simplex_info.workLower_[var]);
      if (!ok) {
        highsLogUser(
            options.log_options, HighsLogType::kError,
            "For variable %" HIGHSINT_FORMAT
            ", simplex_info.workRange_ should be %g = %g - %g "
            "but is %g\n",
            var, simplex_info.workUpper_[var] - simplex_info.workLower_[var],
            simplex_info.workUpper_[var], simplex_info.workLower_[var],
            simplex_info.workRange_[var]);
        return ok;
      }
    }
  }
  // Don't check costs against the LP, when using primal simplex in
  // primal phase 1, if the LP is primal infeasible, or if the costs
  // have been perturbed
  if (!(primal_phase1 ||
        scaled_model_status == HighsModelStatus::kPrimalInfeasible ||
        simplex_info.costs_perturbed)) {
    for (HighsInt col = 0; col < simplex_lp.numCol_; ++col) {
      HighsInt var = col;
      double work_cost = simplex_info.workCost_[var];
      double ok_cost = (HighsInt)simplex_lp.sense_ * simplex_lp.colCost_[col];
      ok = work_cost == ok_cost;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "For col %" HIGHSINT_FORMAT
                     ", simplex_info.workCost_ should be %g but is %g\n",
                     col, ok_cost, simplex_info.workCost_[var]);
        return ok;
      }
    }
    for (HighsInt row = 0; row < simplex_lp.numRow_; ++row) {
      HighsInt var = simplex_lp.numCol_ + row;
      ok = simplex_info.workCost_[var] == 0.;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "For row %" HIGHSINT_FORMAT
                     ", simplex_info.workCost_ should be zero but is %g\n",
                     row, simplex_info.workCost_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool ekkDebugOneNonbasicMoveVsWorkArraysOk(const HEkk& ekk_instance,
                                           const HighsInt var) {
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
          highsLogUser(options.log_options, HighsLogType::kError,
                       "Fixed variable %" HIGHSINT_FORMAT
                       " (simplex_lp.numCol_ = %" HIGHSINT_FORMAT
                       ") [%11g, %11g, "
                       "%11g] so nonbasic "
                       "move should be zero but is %" HIGHSINT_FORMAT "\n",
                       var, simplex_lp.numCol_, simplex_info.workLower_[var],
                       simplex_info.workValue_[var],
                       simplex_info.workUpper_[var],
                       simplex_basis.nonbasicMove_[var]);
          return ok;
        }
        ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "Fixed variable %" HIGHSINT_FORMAT
                       " (simplex_lp.numCol_ = %" HIGHSINT_FORMAT
                       ") so "
                       "simplex_info.work value should be %g but "
                       "is %g\n",
                       var, simplex_lp.numCol_, simplex_info.workLower_[var],
                       simplex_info.workValue_[var]);
          return ok;
        }
      } else {
        // Boxed variable
        ok = (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) ||
             (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_DN);
        if (!ok) {
          highsLogUser(
              options.log_options, HighsLogType::kError,
              "Boxed variable %" HIGHSINT_FORMAT
              " (simplex_lp.numCol_ = %" HIGHSINT_FORMAT
              ") [%11g, %11g, "
              "%11g] range %g so "
              "nonbasic move should be up/down but is  %" HIGHSINT_FORMAT "\n",
              var, simplex_lp.numCol_, simplex_info.workLower_[var],
              simplex_info.workValue_[var], simplex_info.workUpper_[var],
              simplex_info.workUpper_[var] - simplex_info.workLower_[var],
              simplex_basis.nonbasicMove_[var]);
          return ok;
        }
        if (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
          ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
          if (!ok) {
            highsLogUser(options.log_options, HighsLogType::kError,
                         "Boxed variable %" HIGHSINT_FORMAT
                         " (simplex_lp.numCol_ = %" HIGHSINT_FORMAT
                         ") with "
                         "NONBASIC_MOVE_UP so work "
                         "value should be %g but is %g\n",
                         var, simplex_lp.numCol_, simplex_info.workLower_[var],
                         simplex_info.workValue_[var]);
            return ok;
          }
        } else {
          ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
          if (!ok) {
            highsLogUser(options.log_options, HighsLogType::kError,
                         "Boxed variable %" HIGHSINT_FORMAT
                         " (simplex_lp.numCol_ = %" HIGHSINT_FORMAT
                         ") with "
                         "NONBASIC_MOVE_DN so work "
                         "value should be %g but is %g\n",
                         var, simplex_lp.numCol_, simplex_info.workUpper_[var],
                         simplex_info.workValue_[var]);
            return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Finite lower bound and infinite upper bound variable "
                     "%" HIGHSINT_FORMAT
                     " "
                     "(simplex_lp.numCol_ = "
                     "%" HIGHSINT_FORMAT
                     ") [%11g, %11g, %11g] so nonbasic move should be "
                     "up=%2" HIGHSINT_FORMAT
                     " but is  "
                     "%" HIGHSINT_FORMAT "\n",
                     var, simplex_lp.numCol_, simplex_info.workLower_[var],
                     simplex_info.workValue_[var], simplex_info.workUpper_[var],
                     NONBASIC_MOVE_UP, simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Finite lower bound and infinite upper bound variable "
                     "%" HIGHSINT_FORMAT
                     " "
                     "(simplex_lp.numCol_ = "
                     "%" HIGHSINT_FORMAT
                     ") so work value should be %g but is %g\n",
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
        highsLogUser(
            options.log_options, HighsLogType::kError,
            "Finite upper bound and infinite lower bound variable "
            "%" HIGHSINT_FORMAT
            " "
            "(simplex_lp.numCol_ = "
            "%" HIGHSINT_FORMAT
            ") [%11g, %11g, %11g] so nonbasic move should be down but is  "
            "%" HIGHSINT_FORMAT "\n",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Finite upper bound and infinite lower bound variable "
                     "%" HIGHSINT_FORMAT
                     " "
                     "(simplex_lp.numCol_ = "
                     "%" HIGHSINT_FORMAT
                     ") so work value should be %g but is %g\n",
                     var, simplex_lp.numCol_, simplex_info.workUpper_[var],
                     simplex_info.workValue_[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Free variable %" HIGHSINT_FORMAT
                     " (simplex_lp.numCol_ = %" HIGHSINT_FORMAT
                     ") [%11g, %11g, %11g] "
                     "so nonbasic "
                     "move should be zero but is  %" HIGHSINT_FORMAT "\n",
                     var, simplex_lp.numCol_, simplex_info.workLower_[var],
                     simplex_info.workValue_[var], simplex_info.workUpper_[var],
                     simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == 0.0;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Free variable %" HIGHSINT_FORMAT
                     " (simplex_lp.numCol_ = %" HIGHSINT_FORMAT
                     ") so work value should "
                     "be zero but "
                     "is %g\n",
                     var, simplex_lp.numCol_, simplex_info.workValue_[var]);
        return ok;
      }
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
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap) return;
  const double abs_alpha_from_col = fabs(alpha_from_col);
  const double abs_alpha_from_row = fabs(alpha_from_row);
  const double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  const HighsInt iteration_count = ekk_instance.iteration_count_;
  const HighsInt update_count = ekk_instance.simplex_info_.update_count;
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
  highsLogUser(ekk_instance.options_.log_options, HighsLogType::kWarning,
               "%s (%s) [Iter %" HIGHSINT_FORMAT "; Update %" HIGHSINT_FORMAT
               "] Col: %11.4g; Row: %11.4g; Diff "
               "= %11.4g: Measure %11.4g %s %11.4g\n",
               method_name.c_str(), model_name.c_str(), iteration_count,
               update_count, abs_alpha_from_col, abs_alpha_from_row,
               abs_alpha_diff, numerical_trouble_measure, adjective.c_str(),
               numerical_trouble_tolerance);
  if (wrong_sign) {
    highsLogUser(ekk_instance.options_.log_options, HighsLogType::kWarning,
                 "   Incompatible signs for Col: %11.4g and Row: %11.4g\n",
                 alpha_from_col, alpha_from_row);
  }
  if ((numerical_trouble || wrong_sign) && !reinvert) {
    highsLogUser(ekk_instance.options_.log_options, HighsLogType::kWarning,
                 "   Numerical trouble or wrong sign and not reinverting\n");
  }
}

HighsDebugStatus ekkDebugUpdatedDual(const HighsOptions& options,
                                     const double updated_dual,
                                     const double computed_dual) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  std::string error_adjective;
  HighsLogType report_level;
  double updated_dual_absolute_error = fabs(updated_dual - computed_dual);
  double updated_dual_relative_error =
      updated_dual_absolute_error / max(abs(computed_dual), 1.0);
  bool sign_error = updated_dual * computed_dual <= 0;
  bool at_least_small_error =
      sign_error ||
      updated_dual_absolute_error > updated_dual_small_absolute_error ||
      updated_dual_relative_error > updated_dual_small_relative_error;
  if (!at_least_small_error) return return_status;

  if (updated_dual_relative_error > updated_dual_large_relative_error ||
      updated_dual_absolute_error > updated_dual_large_absolute_error) {
    error_adjective = "Large";
    report_level = HighsLogType::kInfo;
    return_status = HighsDebugStatus::kLargeError;
  } else if (updated_dual_relative_error > updated_dual_small_relative_error ||
             updated_dual_absolute_error > updated_dual_small_absolute_error) {
    error_adjective = "Small";
    report_level = HighsLogType::kDetailed;
    return_status = HighsDebugStatus::kSmallError;
  } else {
    error_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  if (sign_error) {
    report_level = HighsLogType::kInfo;
    return_status = HighsDebugStatus::kLargeError;
  }
  highsLogDev(
      options.log_options, report_level,
      "UpdatedDual:  %-9s absolute (%9.4g) or relative (%9.4g) error in "
      "updated dual value",
      error_adjective.c_str(), updated_dual_absolute_error,
      updated_dual_relative_error);
  if (sign_error) {
    highsLogDev(options.log_options, report_level,
                ": Also sign error with (%9.4g, %9.4g)\n", updated_dual,
                computed_dual);
  } else {
    highsLogDev(options.log_options, report_level, "\n");
  }
  return return_status;
}

HighsDebugStatus ekkDebugNonbasicFreeColumnSet(
    const HEkk& ekk_instance, const HighsInt num_free_col,
    const HSet nonbasic_free_col_set) {
  const HighsOptions& options = ekk_instance.options_;
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  const HighsLp& lp = ekk_instance.simplex_lp_;
  const HighsSimplexInfo& simplex_info = ekk_instance.simplex_info_;
  const SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;
  HighsInt num_tot = lp.numCol_ + lp.numRow_;

  // Check the number of free columns
  HighsInt check_num_free_col = 0;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (simplex_info.workLower_[iVar] <= -kHighsInf &&
        simplex_info.workUpper_[iVar] >= kHighsInf)
      check_num_free_col++;
  }
  if (check_num_free_col != num_free_col) {
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "NonbasicFreeColumnData: Number of free columns should be "
                "%" HIGHSINT_FORMAT ", not %" HIGHSINT_FORMAT "\n",
                check_num_free_col, num_free_col);
    return HighsDebugStatus::kLogicalError;
  }
  if (!num_free_col) return HighsDebugStatus::kOk;
  // Debug HSet nonbasic_free_col
  bool nonbasic_free_col_ok = nonbasic_free_col_set.debug();
  if (!nonbasic_free_col_ok) {
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "NonbasicFreeColumnData: HSet error\n");
    return HighsDebugStatus::kLogicalError;
  }

  // Check that we have the right number of nonbasic free columns
  const HighsInt& num_nonbasic_free_col = nonbasic_free_col_set.count();
  HighsInt check_num_nonbasic_free_col = 0;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    bool nonbasic_free =
        simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE &&
        simplex_info.workLower_[iVar] <= -kHighsInf &&
        simplex_info.workUpper_[iVar] >= kHighsInf;
    if (nonbasic_free) check_num_nonbasic_free_col++;
  }
  if (check_num_nonbasic_free_col != num_nonbasic_free_col) {
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "NonbasicFreeColumnData: Set should have %" HIGHSINT_FORMAT
                " entries, not %" HIGHSINT_FORMAT "\n",
                check_num_nonbasic_free_col, num_nonbasic_free_col);
    return HighsDebugStatus::kLogicalError;
  }
  // Check that all in the set are nonbasic free columns
  const vector<HighsInt>& nonbasic_free_col_set_entry =
      nonbasic_free_col_set.entry();
  for (HighsInt ix = 0; ix < num_nonbasic_free_col; ix++) {
    HighsInt iVar = nonbasic_free_col_set_entry[ix];
    bool nonbasic_free =
        simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE &&
        simplex_info.workLower_[iVar] <= -kHighsInf &&
        simplex_info.workUpper_[iVar] >= kHighsInf;
    if (!nonbasic_free) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "NonbasicFreeColumnData: Variable %" HIGHSINT_FORMAT
                  " in nonbasic free "
                  "set has nonbasicFlag = %" HIGHSINT_FORMAT
                  " and bounds [%g, %g]\n",
                  iVar, simplex_basis.nonbasicFlag_[iVar],
                  simplex_info.workLower_[iVar], simplex_info.workUpper_[iVar]);
      return HighsDebugStatus::kLogicalError;
    }
  }
  return HighsDebugStatus::kOk;
}

HighsDebugStatus ekkDebugRowMatrix(const HEkk& ekk_instance) {
  /*
  printf("Checking row-wise matrix\n");
  for (HighsInt row = 0; row < numRow; row++) {
    for (HighsInt el = ARstart[row]; el < AR_Nend[row]; el++) {
      HighsInt col = ARindex[el];
      if (!nonbasicFlag_[col]) {
        printf("Row-wise matrix error: col %" HIGHSINT_FORMAT ", (el = %"
  HIGHSINT_FORMAT " for row %" HIGHSINT_FORMAT ") is basic\n", col, el, row);
        return false;
      }
    }
    for (HighsInt el = AR_Nend[row]; el < ARstart[row + 1]; el++) {
      HighsInt col = ARindex[el];
      if (nonbasicFlag_[col]) {
        printf(
            "Row-wise matrix error: col %" HIGHSINT_FORMAT ", (el = %"
  HIGHSINT_FORMAT " for row %" HIGHSINT_FORMAT ") is nonbasic\n", col, el, row);
        return false;
      }
    }
  }
  */
  return HighsDebugStatus::kOk;
}
