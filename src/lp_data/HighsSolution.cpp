/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolution.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsSolution.h"

#include <string>
#include <vector>

#include "ipm/IpxSolution.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsSolutionDebug.h"
#include "util/HighsUtils.h"

#ifdef IPX_ON
#include "ipm/IpxStatus.h"
#include "ipm/ipx/include/ipx_status.h"
#include "ipm/ipx/src/lp_solver.h"
#endif

void getPrimalDualInfeasibilitiesFromHighsBasicSolution(
    const HighsLp& lp, const HighsBasis& basis, const HighsSolution& solution,
    HighsSolutionParams& solution_params) {
  HighsPrimalDualErrors primal_dual_errors;
  double primal_objective_value;
  double dual_objective_value;
  const int report_level = -1;
  getPrimalDualInfeasibilitiesAndErrorsFromHighsBasicSolution(
      lp, basis, solution, solution_params, primal_dual_errors,
      primal_objective_value, dual_objective_value, report_level);
}

void getPrimalDualInfeasibilitiesAndErrorsFromHighsBasicSolution(
    const HighsLp& lp, const HighsBasis& basis, const HighsSolution& solution,
    HighsSolutionParams& solution_params,
    HighsPrimalDualErrors& primal_dual_errors, double& primal_objective_value,
    double& dual_objective_value, const int report_level) {
  double primal_feasibility_tolerance =
      solution_params.primal_feasibility_tolerance;
  double dual_feasibility_tolerance =
      solution_params.dual_feasibility_tolerance;

  // solution_params are the values computed in this method.
  int& num_primal_infeasibilities = solution_params.num_primal_infeasibilities;
  double& max_primal_infeasibility = solution_params.max_primal_infeasibility;
  double& sum_primal_infeasibilities =
      solution_params.sum_primal_infeasibilities;
  int& num_dual_infeasibilities = solution_params.num_dual_infeasibilities;
  double& max_dual_infeasibility = solution_params.max_dual_infeasibility;
  double& sum_dual_infeasibilities = solution_params.sum_dual_infeasibilities;

  num_primal_infeasibilities = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibilities = 0;
  num_dual_infeasibilities = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibilities = 0;

  vector<double> primal_activities;
  vector<double> dual_activities;
  primal_activities.assign(lp.numRow_, 0);
  dual_activities.resize(lp.numCol_);
  int num_non_basic_var = 0;
  int num_basic_var = 0;

  int& num_nonzero_basic_duals = primal_dual_errors.num_nonzero_basic_duals;
  int& num_large_nonzero_basic_duals =
      primal_dual_errors.num_large_nonzero_basic_duals;
  double& max_nonzero_basic_dual = primal_dual_errors.max_nonzero_basic_dual;
  double& sum_nonzero_basic_duals = primal_dual_errors.sum_nonzero_basic_duals;

  int& num_off_bound_nonbasic = primal_dual_errors.num_off_bound_nonbasic;
  double& max_off_bound_nonbasic = primal_dual_errors.max_off_bound_nonbasic;
  double& sum_off_bound_nonbasic = primal_dual_errors.sum_off_bound_nonbasic;

  int& num_primal_residual = primal_dual_errors.num_primal_residual;
  double& max_primal_residual = primal_dual_errors.max_primal_residual;
  double& sum_primal_residual = primal_dual_errors.sum_primal_residual;

  int& num_dual_residual = primal_dual_errors.num_dual_residual;
  double& max_dual_residual = primal_dual_errors.max_dual_residual;
  double& sum_dual_residual = primal_dual_errors.sum_dual_residual;

  num_nonzero_basic_duals = 0;
  num_large_nonzero_basic_duals = 0;
  max_nonzero_basic_dual = 0;
  sum_nonzero_basic_duals = 0;

  num_off_bound_nonbasic = 0;
  max_off_bound_nonbasic = 0;
  sum_off_bound_nonbasic = 0;
  num_primal_residual = 0;
  max_primal_residual = 0;
  sum_primal_residual = 0;
  num_dual_residual = 0;
  max_dual_residual = 0;
  sum_dual_residual = 0;

  // Initialise the objective value calculations. Done using
  // HighsSolution so offset is vanilla
  primal_objective_value = lp.offset_;
  dual_objective_value = lp.offset_;

  bool header_written = false;
  double off_bound_nonbasic;
  double primal_infeasibility;
  double dual_infeasibility;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double lower = lp.colLower_[iCol];
    double upper = lp.colUpper_[iCol];
    double value = solution.col_value[iCol];
    double dual = solution.col_dual[iCol];
    HighsBasisStatus status = basis.col_status[iCol];
    primal_objective_value += lp.colCost_[iCol] * value;
    if (status != HighsBasisStatus::BASIC) dual_objective_value += value * dual;
    // Flip dual according to lp.sense_
    dual *= (int)lp.sense_;
    bool report = false;
    bool query = analyseVarBasicSolution(
        report, primal_feasibility_tolerance, dual_feasibility_tolerance,
        status, lower, upper, value, dual, num_non_basic_var, num_basic_var,
        off_bound_nonbasic, primal_infeasibility, dual_infeasibility);
    if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
    max_off_bound_nonbasic =
        std::max(off_bound_nonbasic, max_off_bound_nonbasic);
    sum_off_bound_nonbasic += off_bound_nonbasic;
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibilities++;
    max_primal_infeasibility =
        std::max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibilities += primal_infeasibility;
    if (status == HighsBasisStatus::BASIC) {
      double abs_basic_dual = dual_infeasibility;
      if (abs_basic_dual > 0) {
        num_nonzero_basic_duals++;
        if (abs_basic_dual > dual_feasibility_tolerance)
          num_large_nonzero_basic_duals++;
        max_nonzero_basic_dual =
            std::max(abs_basic_dual, max_nonzero_basic_dual);
        sum_nonzero_basic_duals += abs_basic_dual;
      }
    } else {
      if (dual_infeasibility > dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
    report = report_level == 3 || (report_level == 2 && query);
    if (report) {
      if (!header_written) {
        printf(
            "\nColumns\nIndex NonBs Mv [          LB,           UB]       "
            "Primal         Dual    PrimalIfs      DualIfs\n");
        header_written = true;
      }
      printf("%5d %5d [%12g, %12g] %12g %12g", iCol, (int)status, lower, upper,
             value, dual);
      printf(" %12g %12g", primal_infeasibility, dual_infeasibility);
      analyseVarBasicSolution(
          report, primal_feasibility_tolerance, dual_feasibility_tolerance,
          status, lower, upper, value, dual, num_non_basic_var, num_basic_var,
          off_bound_nonbasic, primal_infeasibility, dual_infeasibility);
      printf("\n");
    }
    dual_activities[iCol] = lp.colCost_[iCol];
    for (int el = lp.Astart_[iCol]; el < lp.Astart_[iCol + 1]; el++) {
      int iRow = lp.Aindex_[el];
      double Avalue = lp.Avalue_[el];
      primal_activities[iRow] += value * Avalue;
      dual_activities[iCol] += solution.row_dual[iRow] * Avalue;
    }
  }
  bool report = report_level >= 2;
  header_written = false;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    double primal_residual =
        fabs(primal_activities[iRow] - solution.row_value[iRow]);
    if (primal_residual > primal_feasibility_tolerance) {
      if (report) {
        if (!header_written) {
          printf(
              "\nRow primal residuals\nIndex     Activity     Solution     "
              "Residual\n");
          header_written = true;
        }
        printf("%5d %12g %12g %12g\n", iRow, primal_activities[iRow],
               solution.row_value[iRow], primal_residual);
      }
      num_primal_residual++;
    }
    max_primal_residual = std::max(primal_residual, max_primal_residual);
    sum_primal_residual += primal_residual;
  }
  header_written = false;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double dual_residual =
        fabs(dual_activities[iCol] - solution.col_dual[iCol]);
    if (dual_residual > dual_feasibility_tolerance) {
      if (report) {
        if (!header_written) {
          printf(
              "\nRow dual residuals\nIndex     Activity     Solution     "
              "Residual\n");
          header_written = true;
        }
        printf("%5d %12g %12g %12g\n", iCol, dual_activities[iCol],
               solution.col_dual[iCol], dual_residual);
      }
      num_dual_residual++;
    }
    max_dual_residual = std::max(dual_residual, max_dual_residual);
    sum_dual_residual += dual_residual;
  }
  header_written = false;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    double lower = lp.rowLower_[iRow];
    double upper = lp.rowUpper_[iRow];
    double value = solution.row_value[iRow];
    double dual = -solution.row_dual[iRow];
    HighsBasisStatus status = basis.row_status[iRow];
    if (status != HighsBasisStatus::BASIC) dual_objective_value += value * dual;
    // Flip dual according to lp.sense_
    dual *= (int)lp.sense_;
    bool report = false;
    bool query = analyseVarBasicSolution(
        report, primal_feasibility_tolerance, dual_feasibility_tolerance,
        status, lower, upper, value, dual, num_non_basic_var, num_basic_var,
        off_bound_nonbasic, primal_infeasibility, dual_infeasibility);
    if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
    max_off_bound_nonbasic =
        std::max(off_bound_nonbasic, max_off_bound_nonbasic);
    sum_off_bound_nonbasic += off_bound_nonbasic;
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibilities++;
    max_primal_infeasibility =
        std::max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibilities += primal_infeasibility;
    if (status == HighsBasisStatus::BASIC) {
      double abs_basic_dual = dual_infeasibility;
      if (abs_basic_dual > 0) {
        num_nonzero_basic_duals++;
        if (abs_basic_dual > dual_feasibility_tolerance)
          num_large_nonzero_basic_duals++;
        max_nonzero_basic_dual =
            std::max(abs_basic_dual, max_nonzero_basic_dual);
        sum_nonzero_basic_duals += abs_basic_dual;
      }
    } else {
      if (dual_infeasibility > dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
    report = report_level == 3 || (report_level == 2 && query);
    if (report) {
      if (!header_written) {
        printf(
            "Rows\nIndex NonBs Mv [          LB,           UB]       Primal    "
            "     Dual    PrimalIfs      DualIfs\n");
        header_written = true;
      }
      printf("%5d %5d [%12g, %12g] %12g %12g", iRow, (int)status, lower, upper,
             value, dual);
      printf(" %12g %12g", primal_infeasibility, dual_infeasibility);
      analyseVarBasicSolution(
          report, primal_feasibility_tolerance, dual_feasibility_tolerance,
          status, lower, upper, value, dual, num_non_basic_var, num_basic_var,
          off_bound_nonbasic, primal_infeasibility, dual_infeasibility);
      printf("\n");
    }
  }
}

bool analyseVarBasicSolution(bool report,
                             const double primal_feasibility_tolerance,
                             const double dual_feasibility_tolerance,
                             const HighsBasisStatus status, const double lower,
                             const double upper, const double value,
                             const double dual, int& num_non_basic_var,
                             int& num_basic_var, double& off_bound_nonbasic,
                             double& primal_infeasibility,
                             double& dual_infeasibility) {
  double middle = (lower + upper) * 0.5;

  bool query = false;
  bool count = !report;
  off_bound_nonbasic = 0;
  double primal_residual = std::max(lower - value, value - upper);
  primal_infeasibility = std::max(primal_residual, 0.);
  // ToDo Strange: nonbasic_flag seems to be inverted???
  if (status == HighsBasisStatus::BASIC) {
    // Basic variable: look for primal infeasibility
    if (count) num_basic_var++;
    if (primal_infeasibility > primal_feasibility_tolerance) {
      // Outside a bound
      if (value < lower) {
        query = true;
        if (report)
          printf(": Basic below lower bound by %12g", primal_residual);
      } else {
        query = true;
        if (report)
          printf(": Basic above upper bound by %12g", primal_residual);
      }
    }
    dual_infeasibility = fabs(dual);
    if (dual_infeasibility > dual_feasibility_tolerance) {
      query = true;
      if (report) printf(": Dual infeasibility of %12g", dual_infeasibility);
    }
  } else {
    // Nonbasic variable: look for primal and dual infeasibility
    if (count) num_non_basic_var++;

    if (primal_infeasibility > primal_feasibility_tolerance) {
      // Outside a bound
      dual_infeasibility = 0;
      if (value < lower) {
        query = true;
        if (report)
          printf(": Nonbasic below lower bound by %12g", primal_residual);
      } else {
        query = true;
        if (report)
          printf(": Nonbasic above upper bound by %12g", primal_residual);
      }
    } else if (primal_residual >= -primal_feasibility_tolerance) {
      // At a bound: check for dual feasibility
      if (lower < upper) {
        // Non-fixed variable
        if (value < middle) {
          // At lower
          dual_infeasibility = std::max(-dual, 0.);
          if (dual_infeasibility > dual_feasibility_tolerance) {
            // Dual infeasiblility
            query = true;
            if (report)
              printf(": Dual infeasibility of %12g", dual_infeasibility);
          }
        } else {
          // At Upper
          dual_infeasibility = std::max(dual, 0.);
          if (dual_infeasibility > dual_feasibility_tolerance) {
            // Dual infeasiblility
            query = true;
            if (report)
              printf(": Dual infeasibility of %12g", dual_infeasibility);
          }
        }
      } else {
        // Fixed variable
        dual_infeasibility = 0;
      }
    } else {
      // Between bounds (or free)
      if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
        // Free
        if (report) printf(": Nonbasic free");
      } else {
        query = true;
        if (report) printf(": Nonbasic off bound by %12g", -primal_residual);
        off_bound_nonbasic = -primal_residual;
      }
      dual_infeasibility = fabs(dual);
      if (dual_infeasibility > dual_feasibility_tolerance) {
        query = true;
        if (report) printf(": Dual infeasibility of %12g", dual_infeasibility);
      }
    }
  }
  query = false;
  return query;
}

#ifdef HiGHSDEV
void analyseSimplexAndHighsSolutionDifferences(
    const HighsModelObject& highs_model_object) {
  const HighsSolution& solution = highs_model_object.solution_;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const HighsSolutionParams& scaled_solution_params =
      highs_model_object.scaled_solution_params_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  const HighsScale& scale = highs_model_object.scale_;

  const double scaled_primal_feasibility_tolerance =
      scaled_solution_params.primal_feasibility_tolerance;
  const double scaled_dual_feasibility_tolerance =
      scaled_solution_params.dual_feasibility_tolerance;

  // Go through the columns, finding the differences in nonbasic column values
  // and duals
  int num_nonbasic_col_value_differences = 0;
  double sum_nonbasic_col_value_differences = 0;
  int num_nonbasic_col_dual_differences = 0;
  double sum_nonbasic_col_dual_differences = 0;
  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    int iVar = iCol;
    if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE) {
      // Consider this nonbasic column
      double local_col_value = simplex_info.workValue_[iVar] * scale.col_[iCol];
      double local_col_dual = (int)simplex_lp.sense_ *
                              simplex_info.workDual_[iVar] /
                              (scale.col_[iCol] / scale.cost_);
      double value_difference =
          fabs(local_col_value - solution.col_value[iCol]);
      double dual_difference = fabs(local_col_dual - solution.col_dual[iCol]);
      if (value_difference > scaled_primal_feasibility_tolerance)
        num_nonbasic_col_value_differences++;
      sum_nonbasic_col_value_differences += value_difference;
      if (value_difference > scaled_dual_feasibility_tolerance)
        num_nonbasic_col_dual_differences++;
      sum_nonbasic_col_dual_differences += dual_difference;
    }
  }
  // Go through the rows, finding the differences in nonbasic and
  // basic row values and duals, as well as differences in basic
  // column values and duals
  int num_nonbasic_row_value_differences = 0;
  double sum_nonbasic_row_value_differences = 0;
  int num_nonbasic_row_dual_differences = 0;
  double sum_nonbasic_row_dual_differences = 0;
  int num_basic_col_value_differences = 0;
  double sum_basic_col_value_differences = 0;
  int num_basic_col_dual_differences = 0;
  double sum_basic_col_dual_differences = 0;
  int num_basic_row_value_differences = 0;
  double sum_basic_row_value_differences = 0;
  int num_basic_row_dual_differences = 0;
  double sum_basic_row_dual_differences = 0;

  for (int ix = 0; ix < simplex_lp.numRow_; ix++) {
    int iRow = ix;
    int iVar = simplex_lp.numCol_ + iRow;
    if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE) {
      // Consider this nonbasic row
      double local_row_value =
          -simplex_info.workValue_[iVar] / scale.row_[iRow];
      double local_row_dual = (int)simplex_lp.sense_ *
                              simplex_info.workDual_[iVar] *
                              (scale.row_[iRow] * scale.cost_);
      double value_difference =
          fabs(local_row_value - solution.row_value[iRow]);
      double dual_difference = fabs(local_row_dual - solution.row_dual[iRow]);
      if (value_difference > scaled_primal_feasibility_tolerance)
        num_nonbasic_row_value_differences++;
      sum_nonbasic_row_value_differences += value_difference;
      if (value_difference > scaled_dual_feasibility_tolerance)
        num_nonbasic_row_dual_differences++;
      sum_nonbasic_row_dual_differences += dual_difference;
    }
    // Consider the basic variable associated with this row index
    iVar = simplex_basis.basicIndex_[ix];
    if (iVar < simplex_lp.numCol_) {
      // Consider this basic column
      int iCol = iVar;
      double local_col_value = simplex_info.baseValue_[ix] * scale.col_[iCol];
      double local_col_dual = 0;
      double value_difference =
          fabs(local_col_value - solution.col_value[iCol]);
      double dual_difference = fabs(local_col_dual - solution.col_dual[iCol]);
      if (value_difference > scaled_primal_feasibility_tolerance)
        num_basic_col_value_differences++;
      sum_basic_col_value_differences += value_difference;
      if (value_difference > scaled_dual_feasibility_tolerance)
        num_basic_col_dual_differences++;
      sum_basic_col_dual_differences += dual_difference;
    } else {
      // Consider this basic row
      iRow = iVar - simplex_lp.numCol_;
      double local_row_value = -simplex_info.baseValue_[ix] / scale.row_[iRow];
      double local_row_dual = 0;
      double value_difference =
          fabs(local_row_value - solution.row_value[iRow]);
      double dual_difference = fabs(local_row_dual - solution.row_dual[iRow]);
      if (value_difference > scaled_primal_feasibility_tolerance)
        num_basic_row_value_differences++;
      sum_basic_row_value_differences += value_difference;
      if (value_difference > scaled_dual_feasibility_tolerance)
        num_basic_row_dual_differences++;
      sum_basic_row_dual_differences += dual_difference;
    }
  }
  double acceptable_difference_sum =
      scaled_primal_feasibility_tolerance + scaled_dual_feasibility_tolerance;
  bool significant_nonbasic_value_differences =
      sum_nonbasic_col_value_differences + sum_nonbasic_row_value_differences >
      0;
  bool significant_basic_value_differences =
      sum_basic_col_value_differences + sum_basic_row_value_differences >
      2 * acceptable_difference_sum;
  bool significant_nonbasic_col_dual_differences =
      sum_nonbasic_col_dual_differences > acceptable_difference_sum;
  bool significant_nonbasic_row_dual_differences =
      sum_nonbasic_row_dual_differences > acceptable_difference_sum;
  bool significant_basic_dual_differences =
      sum_basic_col_dual_differences + sum_basic_row_dual_differences > 0;
  if (significant_nonbasic_value_differences ||
      significant_basic_value_differences ||
      significant_nonbasic_col_dual_differences ||
      significant_nonbasic_row_dual_differences ||
      significant_basic_dual_differences) {
    printf(
        "In transition(): There are significant value and dual differences\n");
    /*
      printf("   nonbasic_value_differences = %d\n",
      significant_nonbasic_value_differences); printf(" basic_value_differences
      = %d\n", significant_basic_value_differences); printf("
      nonbasic_col_dual_differences = %d\n",
      significant_nonbasic_col_dual_differences); printf("
      nonbasic_row_dual_differences = %d\n",
      significant_nonbasic_row_dual_differences); printf("
      basic_dual_differences = %d\n", significant_basic_dual_differences);
      */
  } else {
    printf(
        "In transition(): There are no significant value and dual "
        "differences\n");
  }
  if (significant_nonbasic_value_differences) {
    if (sum_nonbasic_col_value_differences > 0)
      printf("Nonbasic column value differences: %6d (%11.4g)\n",
             num_nonbasic_col_value_differences,
             sum_nonbasic_col_value_differences);
    if (sum_nonbasic_row_value_differences > 0)
      printf("Nonbasic row    value differences: %6d (%11.4g)\n",
             num_nonbasic_row_value_differences,
             sum_nonbasic_row_value_differences);
  }
  if (significant_basic_value_differences) {
    if (sum_basic_col_value_differences > acceptable_difference_sum)
      printf("Basic    column value differences: %6d (%11.4g)\n",
             num_basic_col_value_differences, sum_basic_col_value_differences);
    if (sum_basic_row_value_differences > acceptable_difference_sum)
      printf("Basic    row    value differences: %6d (%11.4g)\n",
             num_basic_row_value_differences, sum_basic_row_value_differences);
  }
  if (significant_nonbasic_col_dual_differences)
    printf("Nonbasic column  dual differences: %6d (%11.4g)\n",
           num_nonbasic_col_dual_differences,
           sum_nonbasic_col_dual_differences);
  if (significant_nonbasic_row_dual_differences)
    printf("Nonbasic row     dual differences: %6d (%11.4g)\n",
           num_nonbasic_row_dual_differences,
           sum_nonbasic_row_dual_differences);
  if (significant_basic_dual_differences) {
    if (sum_basic_col_dual_differences > 0)
      printf("Basic    column  dual differences: %6d (%11.4g)\n",
             num_basic_col_dual_differences, sum_basic_col_dual_differences);
    if (sum_basic_row_dual_differences > 0)
      printf("Basic    row     dual differences: %6d (%11.4g)\n",
             num_basic_row_dual_differences, sum_basic_row_dual_differences);
  }
  printf(
      "grep_transition,%s,%.15g,%d,%g,%d,%g,%s,%d,%g,%d,%g,%d,%g,%d,%g,Primal,%"
      "d,%g,%d,%g,Dual,%d,%g,%d,%g\n",
      simplex_lp.model_name_.c_str(), simplex_info.primal_objective_value,
      scaled_solution_params.num_primal_infeasibilities,
      scaled_solution_params.sum_primal_infeasibilities,
      scaled_solution_params.num_dual_infeasibilities,
      scaled_solution_params.sum_dual_infeasibilities,
      utilHighsModelStatusToString(highs_model_object.scaled_model_status_)
          .c_str(),
      num_nonbasic_col_value_differences, sum_nonbasic_col_value_differences,
      num_nonbasic_row_value_differences, sum_nonbasic_row_value_differences,
      num_basic_col_value_differences, sum_basic_col_value_differences,
      num_basic_row_value_differences, sum_basic_row_value_differences,
      num_nonbasic_col_dual_differences, sum_nonbasic_col_dual_differences,
      num_nonbasic_row_dual_differences, sum_nonbasic_row_dual_differences,
      num_basic_col_dual_differences, sum_basic_col_dual_differences,
      num_basic_row_dual_differences, sum_basic_row_dual_differences);
}
#endif

#ifdef IPX_ON
HighsStatus ipxSolutionToHighsSolution(
    FILE* logfile, const HighsLp& lp, const std::vector<double>& rhs,
    const std::vector<char>& constraint_type, const int ipx_num_col,
    const int ipx_num_row, const std::vector<double>& ipx_x,
    const std::vector<double>& ipx_slack_vars,
    // const std::vector<double>& ipx_y,
    HighsSolution& highs_solution) {
  // Resize the HighsSolution
  highs_solution.col_value.resize(lp.numCol_);
  highs_solution.row_value.resize(lp.numRow_);
  //  highs_solution.col_dual.resize(lp.numCol_);
  //  highs_solution.row_dual.resize(lp.numRow_);

  const std::vector<double>& ipx_col_value = ipx_x;
  const std::vector<double>& ipx_row_value = ipx_slack_vars;
  //  const std::vector<double>& ipx_col_dual = ipx_x;
  //  const std::vector<double>& ipx_row_dual = ipx_y;

  // Row activities are needed to set activity values of free rows -
  // which are ignored by IPX
  vector<double> row_activity;
  bool get_row_activities = ipx_num_row < lp.numRow_;
#ifdef HiGHSDEV
  // For debugging, get the row activities if there are any boxed
  // constraints
  get_row_activities = get_row_activities || ipx_num_col > lp.numCol_;
#endif
  if (get_row_activities) row_activity.assign(lp.numRow_, 0);
  for (int col = 0; col < lp.numCol_; col++) {
    highs_solution.col_value[col] = ipx_col_value[col];
    //    highs_solution.col_dual[col] = ipx_col_dual[col];
    if (get_row_activities) {
      // Accumulate row activities to assign value to free rows
      for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
        int row = lp.Aindex_[el];
        row_activity[row] += highs_solution.col_value[col] * lp.Avalue_[el];
      }
    }
  }
  int ipx_row = 0;
  int ipx_slack = lp.numCol_;
  int num_boxed_rows = 0;
  for (int row = 0; row < lp.numRow_; row++) {
    double lower = lp.rowLower_[row];
    double upper = lp.rowUpper_[row];
    if (lower <= -HIGHS_CONST_INF && upper >= HIGHS_CONST_INF) {
      // Free row - removed by IPX so set it to its row activity
      highs_solution.row_value[row] = row_activity[row];
      //      highs_solution.row_dual[row] = 0;
    } else {
      // Non-free row, so IPX will have it
      if ((lower > -HIGHS_CONST_INF && upper < HIGHS_CONST_INF) &&
          (lower < upper)) {
        // Boxed row - look at its slack
        num_boxed_rows++;
        highs_solution.row_value[row] = ipx_col_value[ipx_slack];
        //	highs_solution.row_dual[row] = -ipx_col_dual[ipx_slack];
        // Update the slack to be used for boxed rows
        ipx_slack++;
      } else {
        highs_solution.row_value[row] = rhs[ipx_row] - ipx_row_value[ipx_row];
        //        highs_solution.row_dual[row] = -ipx_row_dual[ipx_row];
      }
      // Update the IPX row index
      ipx_row++;
    }
  }
  assert(ipx_row == ipx_num_row);
  assert(ipx_slack == ipx_num_col);

  // Flip dual according to lp.sense_
  /*
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    highs_solution.col_dual[iCol] *= (int)lp.sense_;
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    highs_solution.row_dual[iRow] *= (int)lp.sense_;
  }
  */
  return HighsStatus::OK;
}

HighsStatus ipxBasicSolutionToHighsBasicSolution(
    FILE* logfile, const HighsLp& lp, const std::vector<double>& rhs,
    const std::vector<char>& constraint_type, const IpxSolution& ipx_solution,
    HighsBasis& highs_basis, HighsSolution& highs_solution) {
  // Resize the HighsSolution and HighsBasis
  highs_solution.col_value.resize(lp.numCol_);
  highs_solution.row_value.resize(lp.numRow_);
  highs_solution.col_dual.resize(lp.numCol_);
  highs_solution.row_dual.resize(lp.numRow_);
  highs_basis.col_status.resize(lp.numCol_);
  highs_basis.row_status.resize(lp.numRow_);

  const std::vector<double>& ipx_col_value = ipx_solution.ipx_col_value;
  const std::vector<double>& ipx_row_value = ipx_solution.ipx_row_value;
  const std::vector<double>& ipx_col_dual = ipx_solution.ipx_col_dual;
  const std::vector<double>& ipx_row_dual = ipx_solution.ipx_row_dual;
  const std::vector<ipx::Int>& ipx_col_status = ipx_solution.ipx_col_status;
  const std::vector<ipx::Int>& ipx_row_status = ipx_solution.ipx_row_status;

  // Set up meaningful names for values of ipx_col_status and ipx_row_status to
  // be used later in comparisons
  const ipx::Int ipx_basic = 0;
  const ipx::Int ipx_nonbasic_at_lb = -1;
  const ipx::Int ipx_nonbasic_at_ub = -2;
  const ipx::Int ipx_superbasic = -3;
  // Row activities are needed to set activity values of free rows -
  // which are ignored by IPX
  vector<double> row_activity;
  bool get_row_activities = ipx_solution.num_row < lp.numRow_;
#ifdef HiGHSDEV
  // For debugging, get the row activities if there are any boxed
  // constraints
  get_row_activities = get_row_activities || ipx_solution.num_col > lp.numCol_;
#endif
  if (get_row_activities) row_activity.assign(lp.numRow_, 0);
  int num_basic_variables = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    bool unrecognised = false;
    if (ipx_col_status[col] == ipx_basic) {
      // Column is basic
      highs_basis.col_status[col] = HighsBasisStatus::BASIC;
      highs_solution.col_value[col] = ipx_col_value[col];
      highs_solution.col_dual[col] = 0;
    } else if (ipx_col_status[col] == ipx_nonbasic_at_lb) {
      // Column is nonbasic at lower bound
      highs_basis.col_status[col] = HighsBasisStatus::LOWER;
      highs_solution.col_value[col] = ipx_col_value[col];
      highs_solution.col_dual[col] = ipx_col_dual[col];
    } else if (ipx_col_status[col] == ipx_nonbasic_at_ub) {
      // Column is nonbasic at upper bound
      highs_basis.col_status[col] = HighsBasisStatus::UPPER;
      highs_solution.col_value[col] = ipx_col_value[col];
      highs_solution.col_dual[col] = ipx_col_dual[col];
    } else if (ipx_col_status[col] == ipx_superbasic) {
      // Column is superbasic
      highs_basis.col_status[col] = HighsBasisStatus::ZERO;
      highs_solution.col_value[col] = ipx_col_value[col];
      highs_solution.col_dual[col] = ipx_col_dual[col];
    } else {
      unrecognised = true;
#ifdef HiGHSDEV
      printf(
          "\nError in IPX conversion: Unrecognised value ipx_col_status[%2d] = "
          "%d\n",
          col, (int)ipx_col_status[col]);
#endif
    }
#ifdef HiGHSDEV
    if (unrecognised)
      printf("Bounds [%11.4g, %11.4g]\n", lp.colLower_[col], lp.colUpper_[col]);
    if (unrecognised)
      printf(
          "Col %2d ipx_col_status[%2d] = %2d; x[%2d] = %11.4g; z[%2d] = "
          "%11.4g\n",
          col, col, (int)ipx_col_status[col], col, ipx_col_value[col], col,
          ipx_col_dual[col]);
#endif
    assert(!unrecognised);
    if (unrecognised) {
      HighsLogMessage(logfile, HighsMessageType::ERROR,
                      "Unrecognised ipx_col_status value from IPX");
      return HighsStatus::Error;
    }
    if (get_row_activities) {
      // Accumulate row activities to assign value to free rows
      for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
        int row = lp.Aindex_[el];
        row_activity[row] += highs_solution.col_value[col] * lp.Avalue_[el];
      }
    }
    if (highs_basis.col_status[col] == HighsBasisStatus::BASIC)
      num_basic_variables++;
  }
  int ipx_row = 0;
  int ipx_slack = lp.numCol_;
  int num_boxed_rows = 0;
  int num_boxed_rows_basic = 0;
  int num_boxed_row_slacks_basic = 0;
  for (int row = 0; row < lp.numRow_; row++) {
    bool unrecognised = false;
    double lower = lp.rowLower_[row];
    double upper = lp.rowUpper_[row];
#ifdef HiGHSDEV
    int this_ipx_row = ipx_row;
#endif
    if (lower <= -HIGHS_CONST_INF && upper >= HIGHS_CONST_INF) {
      // Free row - removed by IPX so make it basic at its row activity
      highs_basis.row_status[row] = HighsBasisStatus::BASIC;
      highs_solution.row_value[row] = row_activity[row];
      highs_solution.row_dual[row] = 0;
    } else {
      // Non-free row, so IPX will have it
      if ((lower > -HIGHS_CONST_INF && upper < HIGHS_CONST_INF) &&
          (lower < upper)) {
        // Boxed row - look at its slack
        num_boxed_rows++;
        double slack_value = ipx_col_value[ipx_slack];
        double slack_dual = ipx_col_dual[ipx_slack];
        double value = slack_value;
        double dual = -slack_dual;
        if (ipx_row_status[ipx_row] == ipx_basic) {
          // Row is basic
          num_boxed_rows_basic++;
          highs_basis.row_status[row] = HighsBasisStatus::BASIC;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = 0;
        } else if (ipx_col_status[ipx_slack] == ipx_basic) {
          // Slack is basic
          num_boxed_row_slacks_basic++;
          highs_basis.row_status[row] = HighsBasisStatus::BASIC;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = 0;
        } else if (ipx_col_status[ipx_slack] == ipx_nonbasic_at_lb) {
          // Slack at lower bound
          highs_basis.row_status[row] = HighsBasisStatus::LOWER;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else if (ipx_col_status[ipx_slack] == ipx_nonbasic_at_ub) {
          // Slack is at its upper bound
          assert(ipx_col_status[ipx_slack] == ipx_nonbasic_at_ub);
          highs_basis.row_status[row] = HighsBasisStatus::UPPER;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else {
          unrecognised = true;
#ifdef HiGHSDEV
          printf(
              "\nError in IPX conversion: Row %2d (IPX row %2d) has "
              "unrecognised value ipx_col_status[%2d] = %d\n",
              row, ipx_row, ipx_slack, (int)ipx_col_status[ipx_slack]);
#endif
        }
        // Update the slack to be used for boxed rows
        ipx_slack++;
      } else if (ipx_row_status[ipx_row] == ipx_basic) {
        // Row is basic
        highs_basis.row_status[row] = HighsBasisStatus::BASIC;
        highs_solution.row_value[row] = rhs[ipx_row] - ipx_row_value[ipx_row];
        highs_solution.row_dual[row] = 0;
      } else {
        // Nonbasic row at fixed value, lower bound or upper bound
        assert(ipx_row_status[ipx_row] ==
               -1);  // const ipx::Int ipx_nonbasic_row = -1;
        double value = rhs[ipx_row] - ipx_row_value[ipx_row];
        double dual = -ipx_row_dual[ipx_row];
        if (constraint_type[ipx_row] == '>') {
          // Row is at its lower bound
          highs_basis.row_status[row] = HighsBasisStatus::LOWER;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else if (constraint_type[ipx_row] == '<') {
          // Row is at its upper bound
          highs_basis.row_status[row] = HighsBasisStatus::UPPER;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else if (constraint_type[ipx_row] == '=') {
          // Row is at its fixed value
          highs_basis.row_status[row] = HighsBasisStatus::LOWER;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else {
          unrecognised = true;
#ifdef HiGHSDEV
          printf(
              "\nError in IPX conversion: Row %2d: cannot handle "
              "constraint_type[%2d] = %d\n",
              row, ipx_row, constraint_type[ipx_row]);
#endif
        }
      }
      // Update the IPX row index
      ipx_row++;
    }
#ifdef HiGHSDEV
    if (unrecognised)
      printf("Bounds [%11.4g, %11.4g]\n", lp.rowLower_[row], lp.rowUpper_[row]);
    if (unrecognised)
      printf(
          "Row %2d ipx_row_status[%2d] = %2d; s[%2d] = %11.4g; y[%2d] = "
          "%11.4g\n",
          row, this_ipx_row, (int)ipx_row_status[this_ipx_row], this_ipx_row,
          ipx_row_value[this_ipx_row], this_ipx_row,
          ipx_row_dual[this_ipx_row]);
#endif
    assert(!unrecognised);
    if (unrecognised) {
      HighsLogMessage(logfile, HighsMessageType::ERROR,
                      "Unrecognised ipx_row_status value from IPX");
      return HighsStatus::Error;
    }
    if (highs_basis.row_status[row] == HighsBasisStatus::BASIC)
      num_basic_variables++;
  }
  assert(num_basic_variables == lp.numRow_);
  highs_basis.valid_ = true;
  assert(ipx_row == ipx_solution.num_row);
  assert(ipx_slack == ipx_solution.num_col);

  // Flip dual according to lp.sense_
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    highs_solution.col_dual[iCol] *= (int)lp.sense_;
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    highs_solution.row_dual[iRow] *= (int)lp.sense_;
  }

#ifdef HiGHSDEV
  if (num_boxed_rows)
    printf("Of %d boxed rows: %d are basic and %d have basic slacks\n",
           num_boxed_rows, num_boxed_rows_basic, num_boxed_row_slacks_basic);
#endif
  return HighsStatus::OK;
}
#endif

std::string iterationsToString(const HighsIterationCounts& iterations_counts) {
  std::string iteration_statement = "";
  bool not_first = false;
  int num_positive_count = 0;
  if (iterations_counts.simplex) num_positive_count++;
  if (iterations_counts.ipm) num_positive_count++;
  if (iterations_counts.crossover) num_positive_count++;
  if (num_positive_count == 0) {
    iteration_statement += "0 iterations; ";
    return iteration_statement;
  }
  if (num_positive_count > 1) iteration_statement += "(";
  int count;
  std::string count_str;
  count = iterations_counts.simplex;
  if (count) {
    count_str = std::to_string(count);
    if (not_first) iteration_statement += "; ";
    iteration_statement += count_str + " " + "Simplex";
    not_first = true;
  }
  count = iterations_counts.ipm;
  if (count) {
    count_str = std::to_string(count);
    if (not_first) iteration_statement += "; ";
    iteration_statement += count_str + " " + "IPM";
    not_first = true;
  }
  count = iterations_counts.crossover;
  if (count) {
    count_str = std::to_string(count);
    if (not_first) iteration_statement += "; ";
    iteration_statement += count_str + " " + "Crossover";
    not_first = true;
  }
  if (num_positive_count > 1) {
    iteration_statement += ") Iterations; ";
  } else {
    iteration_statement += " iterations; ";
  }
  return iteration_statement;
}

void resetModelStatusAndSolutionParams(HighsModelObject& highs_model_object) {
  resetModelStatusAndSolutionParams(
      highs_model_object.unscaled_model_status_,
      highs_model_object.unscaled_solution_params_,
      highs_model_object.options_);
  resetModelStatusAndSolutionParams(highs_model_object.scaled_model_status_,
                                    highs_model_object.scaled_solution_params_,
                                    highs_model_object.options_);
}

void resetModelStatusAndSolutionParams(HighsModelStatus& model_status,
                                       HighsSolutionParams& solution_params,
                                       const HighsOptions& options) {
  model_status = HighsModelStatus::NOTSET;
  resetSolutionParams(solution_params, options);
}

void resetSolutionParams(HighsSolutionParams& solution_params,
                         const HighsOptions& options) {
  // Set the feasibility tolerances - not affected by invalidateSolutionParams
  solution_params.primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  solution_params.dual_feasibility_tolerance =
      options.dual_feasibility_tolerance;

  // Save a copy of the unscaled solution params to recover the iteration counts
  // and objective
  HighsSolutionParams save_solution_params;
  copySolutionObjectiveParams(solution_params, save_solution_params);
  // Invalidate the solution params then reset the feasibility
  // tolerances and recover the objective
  invalidateSolutionParams(solution_params);
  copySolutionObjectiveParams(save_solution_params, solution_params);
}

// Invalidate a HighsSolutionParams instance
void invalidateSolutionParams(HighsSolutionParams& solution_params) {
  solution_params.objective_function_value = 0;
  invalidateSolutionStatusParams(solution_params);
  invalidateSolutionInfeasibilityParams(solution_params);
}

// Invalidate the solution status values in a HighsSolutionParams
// instance.
void invalidateSolutionStatusParams(HighsSolutionParams& solution_params) {
  solution_params.primal_status = PrimalDualStatus::STATUS_NOTSET;
  solution_params.dual_status = PrimalDualStatus::STATUS_NOTSET;
}

// Invalidate the infeasibility values in a HighsSolutionParams
// instance. Setting the number of infeasibilities to negative values
// indicates that they aren't known
void invalidateSolutionInfeasibilityParams(
    HighsSolutionParams& solution_params) {
  solution_params.num_primal_infeasibilities = -1;
  solution_params.sum_primal_infeasibilities = 0;
  solution_params.max_primal_infeasibility = 0;
  solution_params.num_dual_infeasibilities = -1;
  solution_params.sum_dual_infeasibilities = 0;
  solution_params.max_dual_infeasibility = 0;
}

void copySolutionObjectiveParams(
    const HighsSolutionParams& from_solution_params,
    HighsSolutionParams& to_solution_params) {
  to_solution_params.objective_function_value =
      from_solution_params.objective_function_value;
}

void copyFromSolutionParams(HighsInfo& highs_info,
                            const HighsSolutionParams& solution_params) {
  highs_info.primal_status = solution_params.primal_status;
  highs_info.dual_status = solution_params.dual_status;
  highs_info.objective_function_value =
      solution_params.objective_function_value;
  highs_info.num_primal_infeasibilities =
      solution_params.num_primal_infeasibilities;
  highs_info.max_primal_infeasibility =
      solution_params.max_primal_infeasibility;
  highs_info.sum_primal_infeasibilities =
      solution_params.sum_primal_infeasibilities;
  highs_info.num_dual_infeasibilities =
      solution_params.num_dual_infeasibilities;
  highs_info.max_dual_infeasibility = solution_params.max_dual_infeasibility;
  highs_info.sum_dual_infeasibilities =
      solution_params.sum_dual_infeasibilities;
}
