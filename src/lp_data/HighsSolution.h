/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolution.h
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHSSOLUTION_H_
#define LP_DATA_HIGHSSOLUTION_H_

#include <vector>

#include "HConfig.h"
#include "lp_data/HighsLp.h"
//#include "lp_data/HighsOptions.h"

struct HighsSolutionParams {
  double primal_objective_value;
  double dual_objective_value;
  int iteration_count;
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  int num_primal_infeasibilities;
  double sum_primal_infeasibilities;
  double max_primal_infeasibility;
  int num_dual_infeasibilities;
  double sum_dual_infeasibilities;
  double max_dual_infeasibility;
};

bool analyseVarSolution(
			bool report,
			const double primal_feasibility_tolerance,
			const double dual_feasibility_tolerance,
			const HighsBasisStatus status,
			const double lower,
			const double upper,
			const double value,
			const double dual,
			int& num_non_basic_var,
			int& num_basic_var,
			double& off_bound_nonbasic,
			double& primal_infeasibility,
			double& dual_infeasibility) {
  double middle = (lower + upper) * 0.5;

  bool query = false;
  bool count = !report;
  off_bound_nonbasic = 0;
  double primal_residual = max(lower - value, value - upper);
  primal_infeasibility = max(primal_residual, 0.);
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
          dual_infeasibility = max(-dual, 0.);
          if (dual_infeasibility > dual_feasibility_tolerance) {
            // Dual infeasiblility
            query = true;
            if (report)
              printf(": Dual infeasibility of %12g", dual_infeasibility);
          }
        } else {
          // At Upper
          dual_infeasibility = max(dual, 0.);
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

HighsModelStatus analyseHighsSolution(const HighsLp& lp,
				      const HighsBasis& basis,
				      const HighsSolution& solution,
				      HighsSolutionParams& solution_params,
				      const int report_level,
				      const string message) {

  HighsLogMessage(HighsMessageType::INFO, "HiGHS basic solution: Analysis %s", message.c_str());
  double primal_feasibility_tolerance =
      solution_params.primal_feasibility_tolerance;
  double dual_feasibility_tolerance =
      solution_params.dual_feasibility_tolerance;
  vector<double> primal_activities;
  vector<double> dual_activities;
  primal_activities.assign(lp.numRow_, 0);
  dual_activities.resize(lp.numCol_);
  int num_non_basic_var = 0;
  int num_basic_var = 0;
  int num_off_bound_nonbasic = 0;
  double max_off_bound_nonbasic = 0;
  double sum_off_bound_nonbasic = 0;
  bool header_written = false;
  double off_bound_nonbasic;
  double primal_infeasibility;
  double dual_infeasibility;
  int num_primal_infeasibilities = 0;
  double max_primal_infeasibility = 0;
  double sum_primal_infeasibilities = 0;
  int num_dual_infeasibilities = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibilities = 0;
  int num_nonzero_basic_duals = 0;
  int num_large_nonzero_basic_duals = 0;
  double max_nonzero_basic_dual = 0;
  double sum_nonzero_basic_duals = 0;
  // Initialise the objective value calculations
  double primal_objective_value = lp.offset_;
  double dual_objective_value = lp.offset_;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double lower = lp.colLower_[iCol];
    double upper = lp.colUpper_[iCol];
    double value = solution.col_value[iCol];
    double dual = solution.col_dual[iCol];
    HighsBasisStatus status = basis.col_status[iCol];
    primal_objective_value += lp.colCost_[iCol] * value;
    if (status != HighsBasisStatus::BASIC)
      dual_objective_value += value * dual;
    bool report = false;
    bool query = analyseVarSolution(
        report,
	primal_feasibility_tolerance,
	dual_feasibility_tolerance,
	status, lower, upper, value, dual, num_non_basic_var,
        num_basic_var, off_bound_nonbasic, primal_infeasibility,
        dual_infeasibility);
    if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
    max_off_bound_nonbasic = max(off_bound_nonbasic, max_off_bound_nonbasic);
    sum_off_bound_nonbasic += off_bound_nonbasic;
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibilities++;
    max_primal_infeasibility =
        max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibilities += primal_infeasibility;
    if (status == HighsBasisStatus::BASIC) {
      double abs_basic_dual = dual_infeasibility;
      if (abs_basic_dual > 0) {
	num_nonzero_basic_duals++;
	if (abs_basic_dual > dual_feasibility_tolerance) num_large_nonzero_basic_duals++;
	max_nonzero_basic_dual = max(abs_basic_dual, max_nonzero_basic_dual);
	sum_nonzero_basic_duals += abs_basic_dual;
      }
    } else {
      if (dual_infeasibility > dual_feasibility_tolerance)
	num_dual_infeasibilities++;
      max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
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
      analyseVarSolution(
          report,
	  primal_feasibility_tolerance,
	  dual_feasibility_tolerance,
	  status, lower, upper, value, dual, num_non_basic_var,
          num_basic_var, off_bound_nonbasic, primal_infeasibility,
          dual_infeasibility);
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
  int num_primal_residual = 0;
  double max_primal_residual = 0;
  double sum_primal_residual = 0;
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
    max_primal_residual = max(primal_residual, max_primal_residual);
    sum_primal_residual += primal_residual;
  }
  header_written = false;
  int num_dual_residual = 0;
  double max_dual_residual = 0;
  double sum_dual_residual = 0;
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
    max_dual_residual = max(dual_residual, max_dual_residual);
    sum_dual_residual += dual_residual;
  }
  header_written = false;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    double lower = lp.rowLower_[iRow];
    double upper = lp.rowUpper_[iRow];
    double value = solution.row_value[iRow];
    double dual = -solution.row_dual[iRow];
    HighsBasisStatus status = basis.row_status[iRow];
    if (status != HighsBasisStatus::BASIC)
      dual_objective_value += value * dual;
    bool report = false;
    bool query = analyseVarSolution(
        report,
	primal_feasibility_tolerance,
	dual_feasibility_tolerance,
	status, lower, upper, value, dual, num_non_basic_var,
        num_basic_var, off_bound_nonbasic, primal_infeasibility,
        dual_infeasibility);
    if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
    max_off_bound_nonbasic = max(off_bound_nonbasic, max_off_bound_nonbasic);
    sum_off_bound_nonbasic += off_bound_nonbasic;
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibilities++;
    max_primal_infeasibility =
        max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibilities += primal_infeasibility;
    if (status == HighsBasisStatus::BASIC) {
      double abs_basic_dual = dual_infeasibility;
      if (abs_basic_dual > 0) {
	num_nonzero_basic_duals++;
	if (abs_basic_dual > dual_feasibility_tolerance) num_large_nonzero_basic_duals++;
	max_nonzero_basic_dual = max(abs_basic_dual, max_nonzero_basic_dual);
	sum_nonzero_basic_duals += abs_basic_dual;
      }
    } else {
      if (dual_infeasibility > dual_feasibility_tolerance)
	num_dual_infeasibilities++;
      max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
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
      analyseVarSolution(
          report,
	  primal_feasibility_tolerance,
	  dual_feasibility_tolerance,
	  status, lower, upper, value, dual, num_non_basic_var,
          num_basic_var, off_bound_nonbasic, primal_infeasibility,
          dual_infeasibility);
      printf("\n");
    }
  }
  // Save the solution data
  solution_params.primal_objective_value = primal_objective_value;
  solution_params.dual_objective_value = dual_objective_value;
  solution_params.num_primal_infeasibilities = num_primal_infeasibilities;
  solution_params.max_primal_infeasibility = max_primal_infeasibility;
  solution_params.sum_primal_infeasibilities = sum_primal_infeasibilities;
  solution_params.num_dual_infeasibilities = num_dual_infeasibilities;
  solution_params.max_dual_infeasibility = max_dual_infeasibility;
  solution_params.sum_dual_infeasibilities = sum_dual_infeasibilities;

  HighsModelStatus model_status;
  bool primal_feasible =
      num_primal_infeasibilities ==
      0;  // && max_primal_residual < primal_feasibility_tolerance;
  bool dual_feasible = num_dual_infeasibilities ==
                       0;  // && max_dual_residual < dual_feasibility_tolerance;
  if (primal_feasible) {
    if (dual_feasible) {
      model_status = HighsModelStatus::OPTIMAL;
    } else {
      model_status = HighsModelStatus::PRIMAL_FEASIBLE;
    }
  } else {
    if (dual_feasible) {
      model_status = HighsModelStatus::DUAL_FEASIBLE;
    } else {
      model_status = HighsModelStatus::NOTSET;
    }
  }
  if (num_nonzero_basic_duals) {
    HighsLogMessage(
		    HighsMessageType::WARNING,
		    "HiGHS basic solution: %d (%d large) nonzero basic duals; max = %g; sum = %g",
		    num_nonzero_basic_duals, num_large_nonzero_basic_duals, max_nonzero_basic_dual, sum_nonzero_basic_duals);
  }
  if (num_off_bound_nonbasic)  {
    HighsLogMessage(HighsMessageType::WARNING,
                    "Off-bound num/max/sum           %6d/%11.4g/%11.4g",
                    num_off_bound_nonbasic, max_off_bound_nonbasic, sum_off_bound_nonbasic);
  }
  if (report_level>0) {
    HighsLogMessage(HighsMessageType::INFO,
                    "Primal    num/max/sum residuals %6d/%11.4g/%11.4g: num/max/sum "
                    "infeasibilities %6d/%11.4g/%11.4g",
                    num_primal_residual, max_primal_residual,
                    sum_primal_residual, num_primal_infeasibilities,
                    max_primal_infeasibility, sum_primal_infeasibilities);
    HighsLogMessage(HighsMessageType::INFO,
                    "Dual      num/max/sum residuals %6d/%11.4g/%11.4g: num/max/sum "
                    "infeasibilities %6d/%11.4g/%11.4g",
                    num_dual_residual, max_dual_residual, sum_dual_residual,
                    num_dual_infeasibilities, max_dual_infeasibility,
                    sum_dual_infeasibilities);
    double relative_objective_difference =
      fabs(primal_objective_value-dual_objective_value)/
      max(max(1.0, fabs(primal_objective_value)), fabs(dual_objective_value));
    HighsLogMessage(HighsMessageType::INFO,
                    "Relative objective difference = %.4g",
                    relative_objective_difference);
  } 
  std::string model_status_str = utilHighsModelStatusToString(model_status);
  HighsLogMessage(HighsMessageType::INFO,
      "HiGHS basic solution: Iterations = %d; Objective = %.15g; ",
      solution_params.iteration_count, primal_objective_value);
  HighsLogMessage(HighsMessageType::INFO,
		  "Infeasibilities Pr %d(%g); Du %d(%g); Status: %s",
		  solution_params.num_primal_infeasibilities,
		  solution_params.sum_primal_infeasibilities,
		  solution_params.num_dual_infeasibilities,
		  solution_params.sum_dual_infeasibilities,
		  model_status_str.c_str());


  /*
  HighsLogMessage(HighsMessageType::INFO,
      "HiGHS basic solution: Iterations = %d; Objective = %.15g; Infeasibilities Pr %d(%g); Du %d(%g); Status: %s",
      solution_params.iteration_count, primal_objective_value,
      solution_params.num_primal_infeasibilities,
      solution_params.sum_primal_infeasibilities,
      solution_params.num_dual_infeasibilities,
      solution_params.sum_dual_infeasibilities,
      model_status_str.c_str());
  */

#ifdef HiGHSDEV
  printf("grep_AnBsSol,%s,%s,%.15g,%s,%d,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g\n",
	 lp.model_name_.c_str(), message.c_str(), primal_objective_value,
	 utilHighsModelStatusToString(model_status_).c_str(),
	 num_nonzero_basic_duals, num_large_nonzero_basic_duals, max_nonzero_basic_dual, sum_nonzero_basic_duals,
	 num_off_bound_nonbasic, max_off_bound_nonbasic, sum_off_bound_nonbasic,
	 num_primal_residual, max_primal_residual, sum_primal_residual,
	 num_primal_infeasibilities, max_primal_infeasibility, sum_primal_infeasibilities,
	 num_dual_residual, max_dual_residual, sum_dual_residual,
	 num_dual_infeasibilities, max_dual_infeasibility, sum_dual_infeasibilities);
#endif
  return model_status;
}
#endif  // LP_DATA_HIGHSSOLUTION_H_
