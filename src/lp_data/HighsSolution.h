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
//#include <stdlib.h>
#include <string>

#include "HConfig.h"
#include "lp_data/HighsLp.h"
#include "ipm/IpxSolution.h"
#include "lp_data/HighsOptions.h"
#ifdef IPX_ON
#include "ipm/IpxStatus.h"
#include "ipm/ipx/include/ipx_status.h"
#include "ipm/ipx/src/lp_solver.h"
#endif

void copyToSolutionParams(HighsSolutionParams& solution_params, const HighsOptions& options, const HighsSimplexInfo& simplex_info) {
  solution_params.primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  solution_params.dual_feasibility_tolerance = options.dual_feasibility_tolerance;
  solution_params.simplex_iteration_count = simplex_info.iteration_count;
  solution_params.ipm_iteration_count = 0;
  solution_params.crossover_iteration_count = 0;
}

#ifdef IPX_ON
void copyToSolutionParams(HighsSolutionParams& solution_params, const HighsOptions& options, const ipx::Info& ipx_info) {
  solution_params.primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  solution_params.dual_feasibility_tolerance = options.dual_feasibility_tolerance;
  solution_params.simplex_iteration_count = 0;
  solution_params.ipm_iteration_count = (int)ipx_info.iter;
  solution_params.crossover_iteration_count = 0;
}
#endif

void copyFromSolutionParams(HighsSimplexInfo& simplex_info, const HighsSolutionParams& solution_params) {
  simplex_info.iteration_count = solution_params.simplex_iteration_count;
  simplex_info.primal_status = solution_params.primal_status;
  simplex_info.dual_status = solution_params.dual_status;
  simplex_info.primal_objective_value = solution_params.objective_function_value;
  simplex_info.dual_objective_value = solution_params.objective_function_value;
  simplex_info.num_primal_infeasibilities = solution_params.num_primal_infeasibilities;
  simplex_info.max_primal_infeasibility = solution_params.max_primal_infeasibility;
  simplex_info.sum_primal_infeasibilities = solution_params.sum_primal_infeasibilities;
  simplex_info.num_dual_infeasibilities = solution_params.num_dual_infeasibilities;
  simplex_info.max_dual_infeasibility = solution_params.max_dual_infeasibility;
  simplex_info.sum_dual_infeasibilities = solution_params.sum_dual_infeasibilities;
}

void copyFromSolutionParams(HighsInfo& highs_info, const HighsSolutionParams& solution_params) {
  highs_info.simplex_iteration_count = solution_params.simplex_iteration_count;
  highs_info.ipm_iteration_count = solution_params.ipm_iteration_count;
  highs_info.crossover_iteration_count = solution_params.crossover_iteration_count;
  highs_info.primal_status = solution_params.primal_status;
  highs_info.dual_status = solution_params.dual_status;
  highs_info.objective_function_value = solution_params.objective_function_value;
  highs_info.num_primal_infeasibilities = solution_params.num_primal_infeasibilities;
  highs_info.max_primal_infeasibility = solution_params.max_primal_infeasibility;
  highs_info.sum_primal_infeasibilities = solution_params.sum_primal_infeasibilities;
  highs_info.num_dual_infeasibilities = solution_params.num_dual_infeasibilities;
  highs_info.max_dual_infeasibility = solution_params.max_dual_infeasibility;
  highs_info.sum_dual_infeasibilities = solution_params.sum_dual_infeasibilities;
}

#ifdef IPX_ON
HighsStatus ipxToHighsBasicSolution(const HighsLp& lp,
				    const std::vector<double>& rhs,
				    const std::vector<char>& constraint_type,
				    const IpxSolution& ipx_solution,
				    HighsBasis& highs_basis,
				    HighsSolution& highs_solution) {
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

  // Set up meaningful names for values of ipx_col_status and ipx_row_status to be
  // used later in comparisons
  const ipx::Int ipx_basic = 0;
  const ipx::Int ipx_nonbasic_at_lb = -1;
  const ipx::Int ipx_nonbasic_at_ub = -2;
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
    } else {
      unrecognised = true;
#ifdef HiGHSDEV
      printf("\nError in IPX conversion: Unrecognised value ipx_col_status[%2d] = %d\n", col, (int)ipx_col_status[col]);
#endif	    
    }
#ifdef HiGHSDEV
      if (unrecognised) printf("Bounds [%11.4g, %11.4g]\n", lp.colLower_[col], lp.colUpper_[col]);
      if (unrecognised)
	printf("Col %2d ipx_col_status[%2d] = %2d; x[%2d] = %11.4g; z[%2d] = %11.4g\n",
	       col, col, (int)ipx_col_status[col], col, ipx_col_value[col], col, ipx_col_dual[col]);
#endif
      assert(!unrecognised);
      if (unrecognised) {
	HighsLogMessage(HighsMessageType::ERROR, "Unrecognised ipx_col_status value from IPX");
	return HighsStatus::Error;
      }
      if (get_row_activities) {
	// Accumulate row activities to assign value to free rows
	for (int el = lp.Astart_[col]; el < lp.Astart_[col+1]; el++) {
	  int row = lp.Aindex_[el];
	  row_activity[row] += highs_solution.col_value[col]*lp.Avalue_[el];
	}
      }
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
	if ((lower > -HIGHS_CONST_INF && upper < HIGHS_CONST_INF) && (lower < upper)) {
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
	    printf("\nError in IPX conversion: Row %2d (IPX row %2d) has unrecognised value ipx_col_status[%2d] = %d\n",
		   row, ipx_row, ipx_slack, (int)ipx_col_status[ipx_slack]);
#endif	    
	  }
	  // Update the slack to be used for boxed rows
	  ipx_slack++;
	} else if (ipx_row_status[ipx_row] == ipx_basic) {
	  // Row is basic
	  highs_basis.row_status[row] = HighsBasisStatus::BASIC;
	  highs_solution.row_value[row] = rhs[ipx_row]-ipx_row_value[ipx_row];
	  highs_solution.row_dual[row] = 0;
	} else {
	  // Nonbasic row at fixed value, lower bound or upper bound
	  assert(ipx_row_status[ipx_row] == -1);// const ipx::Int ipx_nonbasic_row = -1;
	  double value = rhs[ipx_row]-ipx_row_value[ipx_row];
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
	    printf("\nError in IPX conversion: Row %2d: cannot handle constraint_type[%2d] = %d\n", row, ipx_row, constraint_type[ipx_row]);
#endif	    
	  }
	}
	  // Update the IPX row index
	ipx_row++;
      }
#ifdef HiGHSDEV
      if (unrecognised) printf("Bounds [%11.4g, %11.4g]\n", lp.rowLower_[row], lp.rowUpper_[row]);
      if (unrecognised)
	printf("Row %2d ipx_row_status[%2d] = %2d; s[%2d] = %11.4g; y[%2d] = %11.4g\n",
	       row, this_ipx_row, (int)ipx_row_status[this_ipx_row], this_ipx_row, ipx_row_value[this_ipx_row], this_ipx_row, ipx_row_dual[this_ipx_row]);
#endif      
      assert(!unrecognised);
      if (unrecognised) {
	HighsLogMessage(HighsMessageType::ERROR, "Unrecognised ipx_row_status value from IPX");
	return HighsStatus::Error;
      }
    }
    assert(ipx_row == ipx_solution.num_row);
    assert(ipx_slack == ipx_solution.num_col);

#ifdef HiGHSDEV
    if (num_boxed_rows)
      printf("Of %d boxed rows: %d are basic and %d have basic slacks\n", 
	     num_boxed_rows, num_boxed_rows_basic, num_boxed_row_slacks_basic);
#endif
    return HighsStatus::OK; 
}
#endif    

bool analyseVarBasicSolution(
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

std::string iterationsToString(const HighsSolutionParams& solution_params) {
  std::string iteration_statement = "";
  bool not_first = false;
  int num_positive_count = 0;
  if (solution_params.simplex_iteration_count) num_positive_count++;
  if (solution_params.ipm_iteration_count) num_positive_count++;
  if (solution_params.crossover_iteration_count) num_positive_count++;
  if (num_positive_count == 0) {
    iteration_statement += "0 iterations; ";
    return iteration_statement;
  }
  if (num_positive_count > 1) iteration_statement += "(";
  int count;
  std::string count_str;
  count = solution_params.simplex_iteration_count;
  if (count) {
    count_str = std::to_string(count);
    if (not_first) iteration_statement += "; ";
    iteration_statement += count_str + " " + "Simplex";
    not_first = true;
  }
  count = solution_params.ipm_iteration_count;
  if (count) {
    count_str = std::to_string(count);
    if (not_first) iteration_statement += "; ";
    iteration_statement += count_str + " " + "IPM";
    not_first = true;
  }
  count = solution_params.crossover_iteration_count;
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

// Analyse the HiGHS basic solution of the given LP. Currently only
// used with the unscaled LP, but would work just as well with a
// scaled LP. The primal and dual feasibility tolerances are passed in
// via solution_params, which returns the int and double data obtained
// about the solution. The overall model status is returned in the
// argument.
HighsModelStatus analyseHighsBasicSolution(const HighsLp& lp,
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
    bool query = analyseVarBasicSolution(
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
      analyseVarBasicSolution(
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
    bool query = analyseVarBasicSolution(
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
      analyseVarBasicSolution(
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
  solution_params.objective_function_value = primal_objective_value;
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
  // Determine the model status
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
  // Determine the primal status
  if (primal_feasible) {
    solution_params.primal_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
  } else {
    solution_params.primal_status = PrimalDualStatus::STATUS_NO_SOLUTION;
  }
  // Determine the dual status
  if (dual_feasible) {
    solution_params.dual_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
  } else {
    solution_params.dual_status = PrimalDualStatus::STATUS_NO_SOLUTION;
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
  HighsLogMessage(HighsMessageType::INFO,
		  "HiGHS basic solution: %sObjective = %.15g; Status: %s",
		  iterationsToString(solution_params).c_str(),
		  primal_objective_value,
		  utilHighsModelStatusToString(model_status).c_str());

  HighsLogMessage(HighsMessageType::INFO,
		  "HiGHS basic solution: Infeasibilities Pr %d(%g); Du %d(%g)",
		  solution_params.num_primal_infeasibilities,
		  solution_params.sum_primal_infeasibilities,
		  solution_params.num_dual_infeasibilities,
		  solution_params.sum_dual_infeasibilities);

#ifdef HiGHSDEV
  printf("grep_AnBsSol,%s,%s,%d,%d,%d,%.15g,%s,%d,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g\n",
	 lp.model_name_.c_str(), message.c_str(),
	 solution_params.simplex_iteration_count,
	 solution_params.ipm_iteration_count,
	 solution_params.crossover_iteration_count,
	 primal_objective_value,
	 utilHighsModelStatusToString(model_status).c_str(),
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
