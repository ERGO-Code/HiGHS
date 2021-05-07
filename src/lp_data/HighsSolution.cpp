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
/**@file lp_data/HighsSolution.cpp
 * @brief Class-independent utilities for HiGHS
 */
#include "lp_data/HighsSolution.h"

#include <string>
#include <vector>

#include "io/HighsIO.h"
#include "ipm/IpxSolution.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsSolutionDebug.h"
#include "util/HighsUtils.h"

#ifdef IPX_ON
#include "ipm/ipx/include/ipx_status.h"
#include "ipm/ipx/src/lp_solver.h"
#endif

void getKktFailures(const HighsLp& lp, const HighsSolution& solution,
                    const HighsBasis& basis,
                    HighsSolutionParams& solution_params) {
  HighsPrimalDualErrors primal_dual_errors;
  getKktFailures(lp, solution, basis, solution_params, primal_dual_errors);
}

void getKktFailures(const HighsLp& lp, const HighsSolution& solution,
                    const HighsBasis& basis,
                    HighsSolutionParams& solution_params,
                    HighsPrimalDualErrors& primal_dual_errors,
		    const bool get_residuals) {
  double primal_feasibility_tolerance =
      solution_params.primal_feasibility_tolerance;
  double dual_feasibility_tolerance =
      solution_params.dual_feasibility_tolerance;
  // solution_params are the values computed in this method.
  HighsInt& num_primal_infeasibility = solution_params.num_primal_infeasibility;
  double& max_primal_infeasibility = solution_params.max_primal_infeasibility;
  double& sum_primal_infeasibility = solution_params.sum_primal_infeasibility;
  HighsInt& num_dual_infeasibility = solution_params.num_dual_infeasibility;
  double& max_dual_infeasibility = solution_params.max_dual_infeasibility;
  double& sum_dual_infeasibility = solution_params.sum_dual_infeasibility;

  num_primal_infeasibility = kHighsIllegalInfeasibilityCount;
  max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  sum_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  solution_params.primal_status = kHighsPrimalDualStatusNoSolution;

  num_dual_infeasibility = kHighsIllegalInfeasibilityCount;
  max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  sum_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  solution_params.dual_status = kHighsPrimalDualStatusNoSolution;

  const bool& have_primal_solution = solution.value_valid;
  const bool& have_dual_solution = solution.dual_valid;
  const bool& have_basis = basis.valid;
  // Check that there is no dual solution if there's no primal solution
  assert(have_primal_solution || !have_dual_solution);
  // Check that there is no basis if there's no dual solution
  assert(have_dual_solution || !have_basis);

  if (have_primal_solution) {
    // There's a primal solution, so check its size and initialise the
    // infeasiblilty counts
    assert((int)solution.col_value.size() >= lp.numCol_);
    assert((int)solution.row_value.size() >= lp.numRow_);
    num_primal_infeasibility = 0;
    max_primal_infeasibility = 0;
    sum_primal_infeasibility = 0;
    if (have_dual_solution) {
      // There's a dual solution, so check its size and initialise the
      // infeasiblilty counts
      assert((int)solution.col_dual.size() >= lp.numCol_);
      assert((int)solution.row_dual.size() >= lp.numRow_);
      num_dual_infeasibility = 0;
      max_dual_infeasibility = 0;
      sum_dual_infeasibility = 0;
    }
  }

  HighsInt& num_primal_residual = primal_dual_errors.num_primal_residual;
  double& max_primal_residual = primal_dual_errors.max_primal_residual;
  double& sum_primal_residual = primal_dual_errors.sum_primal_residual;

  HighsInt& num_dual_residual = primal_dual_errors.num_dual_residual;
  double& max_dual_residual = primal_dual_errors.max_dual_residual;
  double& sum_dual_residual = primal_dual_errors.sum_dual_residual;

  HighsInt& num_nonzero_basic_duals =
      primal_dual_errors.num_nonzero_basic_duals;
  HighsInt& num_large_nonzero_basic_duals =
      primal_dual_errors.num_large_nonzero_basic_duals;
  double& max_nonzero_basic_dual = primal_dual_errors.max_nonzero_basic_dual;
  double& sum_nonzero_basic_duals = primal_dual_errors.sum_nonzero_basic_duals;

  HighsInt& num_off_bound_nonbasic = primal_dual_errors.num_off_bound_nonbasic;
  double& max_off_bound_nonbasic = primal_dual_errors.max_off_bound_nonbasic;
  double& sum_off_bound_nonbasic = primal_dual_errors.sum_off_bound_nonbasic;

  // Initialise HighsPrimalDualErrors

  if (have_primal_solution && get_residuals) {
    num_primal_residual = 0;
    max_primal_residual = 0;
    sum_primal_residual = 0;
  } else {
    num_primal_residual = kHighsIllegalInfeasibilityCount;
    max_primal_residual = kHighsIllegalInfeasibilityMeasure;
    sum_primal_residual = kHighsIllegalInfeasibilityMeasure;
  }
  if (have_dual_solution && get_residuals) {
    num_dual_residual = 0;
    max_dual_residual = 0;
    sum_dual_residual = 0;
  } else {
    num_dual_residual = kHighsIllegalInfeasibilityCount;
    max_dual_residual = kHighsIllegalInfeasibilityMeasure;
    sum_dual_residual = kHighsIllegalInfeasibilityMeasure;
  }
  if (have_basis) {
    num_nonzero_basic_duals = 0;
    num_large_nonzero_basic_duals = 0;
    max_nonzero_basic_dual = 0;
    sum_nonzero_basic_duals = 0;

    num_off_bound_nonbasic = 0;
    max_off_bound_nonbasic = 0;
    sum_off_bound_nonbasic = 0;
  } else {
    num_nonzero_basic_duals = kHighsIllegalInfeasibilityCount;
    num_large_nonzero_basic_duals = kHighsIllegalInfeasibilityCount;
    max_nonzero_basic_dual = kHighsIllegalInfeasibilityMeasure;
    sum_nonzero_basic_duals = kHighsIllegalInfeasibilityMeasure;

    num_off_bound_nonbasic = kHighsIllegalInfeasibilityCount;
    max_off_bound_nonbasic = kHighsIllegalInfeasibilityMeasure;
    sum_off_bound_nonbasic = kHighsIllegalInfeasibilityMeasure;
  }
  // Without a primal solution, nothing can be done!
  if (!have_primal_solution) return;
  std::vector<double> primal_activities;
  std::vector<double> dual_activities;
  if (get_residuals) {
    primal_activities.assign(lp.numRow_, 0);
    if (have_dual_solution) 
      dual_activities.resize(lp.numCol_);
  }
  HighsInt num_basic_var = 0;
  HighsInt num_non_basic_var = 0;

  // Set status to a value so the compiler doesn't think it might be unassigned.
  HighsBasisStatus status = HighsBasisStatus::kNonbasic;
  // Set status_pointer to be NULL unless we have a basis, in which
  // case it's the pointer to status. Can't make it a pointer to the
  // entry of basis since this is const.
  HighsBasisStatus* status_pointer = NULL;
  if (have_basis) status_pointer = &status;

  double primal_infeasibility;
  double dual_infeasibility;
  double value_residual;
  double lower;
  double upper;
  double value;
  double dual = 0;
  for (HighsInt iVar = 0; iVar < lp.numCol_ + lp.numRow_; iVar++) {
    if (iVar < lp.numCol_) {
      HighsInt iCol = iVar;
      lower = lp.colLower_[iCol];
      upper = lp.colUpper_[iCol];
      value = solution.col_value[iCol];
      if (have_dual_solution) dual = solution.col_dual[iCol];
      if (have_basis) status = basis.col_status[iCol];
    } else {
      HighsInt iRow = iVar - lp.numCol_;
      lower = lp.rowLower_[iRow];
      upper = lp.rowUpper_[iRow];
      value = solution.row_value[iRow];
      if (have_dual_solution) dual = -solution.row_dual[iRow];
      if (have_basis) status = basis.row_status[iRow];
    }
    // Flip dual according to lp.sense_
    dual *= (HighsInt)lp.sense_;
    getVariableKktFailures(  // const HighsLogOptions& log_options,
        primal_feasibility_tolerance, dual_feasibility_tolerance, lower, upper,
        value, dual, status_pointer, primal_infeasibility, dual_infeasibility,
        value_residual);
    // Accumulate primal infeasiblilties
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibility++;
    max_primal_infeasibility =
        std::max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibility += primal_infeasibility;

    if (have_dual_solution) {
      // Accumulate dual infeasiblilties
      if (dual_infeasibility > dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
    if (have_basis) {
      if (status == HighsBasisStatus::kBasic) {
        num_basic_var++;
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
        num_non_basic_var++;
        double off_bound_nonbasic = value_residual;
        if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
        max_off_bound_nonbasic =
            std::max(off_bound_nonbasic, max_off_bound_nonbasic);
        sum_off_bound_nonbasic += off_bound_nonbasic;
      }
    }
    if (iVar < lp.numCol_ && get_residuals) {
      HighsInt iCol = iVar;
      if (have_dual_solution) dual_activities[iCol] = lp.colCost_[iCol];
      for (HighsInt el = lp.Astart_[iCol]; el < lp.Astart_[iCol + 1]; el++) {
        HighsInt iRow = lp.Aindex_[el];
        double Avalue = lp.Avalue_[el];
        primal_activities[iRow] += value * Avalue;
        if (have_dual_solution)
          dual_activities[iCol] += solution.row_dual[iRow] * Avalue;
      }
    }
  }
  if (get_residuals) {
    const double large_residual_error = 1e-12;
    for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
      double primal_residual_error =
        std::fabs(primal_activities[iRow] - solution.row_value[iRow]);
      if (primal_residual_error > large_residual_error) num_primal_residual++;
      max_primal_residual = std::max(primal_residual_error, max_primal_residual);
      sum_primal_residual += primal_residual_error;
    }
    if (have_dual_solution) {
      for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
	double dual_residual_error =
          std::fabs(dual_activities[iCol] - solution.col_dual[iCol]);
	if (dual_residual_error > large_residual_error) num_dual_residual++;
	max_dual_residual = std::max(dual_residual_error, max_dual_residual);
	sum_dual_residual += dual_residual_error;
      }
    }
  }
  // Assign primal solution status
  if (num_primal_infeasibility) {
    solution_params.primal_status = kHighsPrimalDualStatusInfeasiblePoint;
  } else {
    solution_params.primal_status = kHighsPrimalDualStatusFeasiblePoint;
  }
  if (have_dual_solution) {
    // Assign dual solution status
    if (num_dual_infeasibility) {
      solution_params.dual_status = kHighsPrimalDualStatusInfeasiblePoint;
    } else {
      solution_params.dual_status = kHighsPrimalDualStatusFeasiblePoint;
    }
  }
}
// Gets the KKT failures for a variable. The lack of a basis status is
// indicated by status_pointer being null.
//
// Value and dual are used compute the primal and dual infeasibility -
// according to the basis status (if valid) or primal value.It's up to
// the calling method to ignore these if the value or dual are not
// valid.
//
// If the basis status is valid, then the numbers of basic and
// nonbasic variables are updated, and the extent to which a nonbasic
// variable is off its bound is returned.
void getVariableKktFailures(
    const double primal_feasibility_tolerance,
    const double dual_feasibility_tolerance, const double lower,
    const double upper, const double value, const double dual,
    HighsBasisStatus* status_pointer, double& primal_infeasibility,
    double& dual_infeasibility, double& value_residual) {
  const double middle = (lower + upper) * 0.5;
  // @primal_infeasibility calculation
  primal_infeasibility = 0;
  if (value < lower - primal_feasibility_tolerance) {
    // Below lower
    primal_infeasibility = lower - value;
  } else if (value > upper + primal_feasibility_tolerance) {
    // Above upper
    primal_infeasibility = value - upper;
  }
  value_residual = std::min(std::fabs(lower - value), std::fabs(value - upper));
  // Determine whether the variable is at a bound using the basis
  // status (if valid) or value residual.
  bool at_a_bound = value_residual <= primal_feasibility_tolerance;
  if (status_pointer != NULL) {
    // If the variable is basic, then consider it not to be at a bound
    // so that any dual value yields an infeasibility value
    if (*status_pointer == HighsBasisStatus::kBasic) at_a_bound = false;
  }
  if (at_a_bound) {
    // At a bound
    double middle = (lower + upper) * 0.5;
    if (lower < upper) {
      // Non-fixed variable
      if (value < middle) {
        // At lower
        dual_infeasibility = std::max(-dual, 0.);
      } else {
        // At upper
        dual_infeasibility = std::max(dual, 0.);
      }
    } else {
      // Fixed variable
      dual_infeasibility = 0;
    }
  } else {
    // Off bounds (or free)
    dual_infeasibility = fabs(dual);
  }
}

void getReportKktFailures(const HighsOptions& options, const HighsLp& lp,
                          const HighsSolution& solution,
                          const HighsBasis& basis) {
  const HighsLogOptions& log_options = options.log_options;
  HighsSolutionParams solution_params;
  solution_params.primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  solution_params.dual_feasibility_tolerance =
      options.dual_feasibility_tolerance;
  HighsPrimalDualErrors primal_dual_errors;
  getKktFailures(lp, solution, basis, solution_params, primal_dual_errors, true);

  bool have_kkt_failure = false;
  highsLogDev(log_options, HighsLogType::kInfo, "KKT Failures\n");
  have_kkt_failure =
      reportKktFailureInt("Number of primal infeasibilities         ",
                          log_options,
                          solution_params.num_primal_infeasibility) ||
      have_kkt_failure;
  have_kkt_failure =
      reportKktFailureDouble("Maximum   primal infeasibility           ",
                             log_options,
                             solution_params.max_primal_infeasibility) ||
      have_kkt_failure;
  have_kkt_failure =
      reportKktFailureDouble("Sum of    primal infeasibilities         ",
                             log_options,
                             solution_params.sum_primal_infeasibility) ||
      have_kkt_failure;

  have_kkt_failure = reportKktFailureInt(
                         "Number of dual infeasibilities           ",
                         log_options, solution_params.num_dual_infeasibility) ||
                     have_kkt_failure;
  have_kkt_failure = reportKktFailureDouble(
                         "Maximum   dual infeasibility             ",
                         log_options, solution_params.max_dual_infeasibility) ||
                     have_kkt_failure;
  have_kkt_failure = reportKktFailureDouble(
                         "Sum of    dual infeasibilities           ",
                         log_options, solution_params.sum_dual_infeasibility) ||
                     have_kkt_failure;

  have_kkt_failure =
      reportKktFailureInt("Number of       nonzero basic dual values",
                          log_options,
                          primal_dual_errors.num_nonzero_basic_duals) ||
      have_kkt_failure;
  have_kkt_failure =
      reportKktFailureInt("Number of large nonzero basic dual values",
                          log_options,
                          primal_dual_errors.num_large_nonzero_basic_duals) ||
      have_kkt_failure;
  have_kkt_failure =
      reportKktFailureDouble("Maximum         nonzero basic dual value ",
                             log_options,
                             primal_dual_errors.max_nonzero_basic_dual) ||
      have_kkt_failure;
  have_kkt_failure =
      reportKktFailureDouble("Sum of          nonzero basic dual values",
                             log_options,
                             primal_dual_errors.sum_nonzero_basic_duals) ||
      have_kkt_failure;

  have_kkt_failure =
      reportKktFailureInt("Number of off-bound nonbasic values      ",
                          log_options,
                          primal_dual_errors.num_off_bound_nonbasic) ||
      have_kkt_failure;
  have_kkt_failure =
      reportKktFailureDouble("Maximum   off-bound nonbasic value       ",
                             log_options,
                             primal_dual_errors.max_off_bound_nonbasic) ||
      have_kkt_failure;
  have_kkt_failure =
      reportKktFailureDouble("Sum of    off-bound nonbasic values      ",
                             log_options,
                             primal_dual_errors.sum_off_bound_nonbasic) ||
      have_kkt_failure;
  have_kkt_failure = reportKktFailureInt(
                         "Number of primal residuals               ",
                         log_options, primal_dual_errors.num_primal_residual) ||
                     have_kkt_failure;
  have_kkt_failure = reportKktFailureDouble(
                         "Maximum   primal residual                ",
                         log_options, primal_dual_errors.max_primal_residual) ||
                     have_kkt_failure;
  have_kkt_failure = reportKktFailureDouble(
                         "Sum of    primal residuals               ",
                         log_options, primal_dual_errors.sum_primal_residual) ||
                     have_kkt_failure;
  have_kkt_failure =
      reportKktFailureInt("Number of dual residuals                 ",
                          log_options, primal_dual_errors.num_dual_residual) ||
      have_kkt_failure;
  have_kkt_failure = reportKktFailureDouble(
                         "Maximum   dual residual                  ",
                         log_options, primal_dual_errors.max_dual_residual) ||
                     have_kkt_failure;
  have_kkt_failure = reportKktFailureDouble(
                         "Sum of    dual residuals                 ",
                         log_options, primal_dual_errors.sum_dual_residual) ||
                     have_kkt_failure;
  if (!have_kkt_failure)
    highsLogDev(log_options, HighsLogType::kInfo, "   None!\n");
}

bool reportKktFailureInt(const std::string message,
                         const HighsLogOptions& log_options,
                         const HighsInt value) {
  if (value > 0) {
    highsLogDev(log_options, HighsLogType::kInfo, "%s: %" HIGHSINT_FORMAT "\n",
                message.c_str(), value);
    return true;
  }
  return false;
}

bool reportKktFailureDouble(const std::string message,
                            const HighsLogOptions& log_options,
                            const double value) {
  if (value > 0) {
    highsLogDev(log_options, HighsLogType::kInfo, "%s: %g\n", message.c_str(),
                value);
    return true;
  }
  return false;
}

void getPrimalDualInfeasibilities(const HighsLp& lp,
                                  const HighsSolution& solution,
                                  HighsSolutionParams& solution_params) {
  double primal_feasibility_tolerance =
      solution_params.primal_feasibility_tolerance;
  double dual_feasibility_tolerance =
      solution_params.dual_feasibility_tolerance;

  // solution_params are the values computed in this method.
  HighsInt& num_primal_infeasibility = solution_params.num_primal_infeasibility;
  double& max_primal_infeasibility = solution_params.max_primal_infeasibility;
  double& sum_primal_infeasibility = solution_params.sum_primal_infeasibility;
  HighsInt& num_dual_infeasibility = solution_params.num_dual_infeasibility;
  double& max_dual_infeasibility = solution_params.max_dual_infeasibility;
  double& sum_dual_infeasibility = solution_params.sum_dual_infeasibility;

  num_primal_infeasibility = kHighsIllegalInfeasibilityCount;
  max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  sum_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  solution_params.primal_status = kHighsPrimalDualStatusNoSolution;

  num_dual_infeasibility = kHighsIllegalInfeasibilityCount;
  max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  sum_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  solution_params.dual_status = kHighsPrimalDualStatusNoSolution;

  const bool have_primal_solution = solution.value_valid;
  if (have_primal_solution) {
    // There's a primal solution, so check its size and initialise the
    // infeasiblilty counts
    assert((int)solution.col_value.size() >= lp.numCol_);
    assert((int)solution.row_value.size() >= lp.numRow_);
    num_primal_infeasibility = 0;
    max_primal_infeasibility = 0;
    sum_primal_infeasibility = 0;
  } else {
    // Nothing can be done, so return
    return;
  }
  const bool have_dual_solution = solution.dual_valid;
  if (have_dual_solution) {
    // There's a dual solution, so check its size and initialise the
    // infeasiblilty counts
    assert((int)solution.col_dual.size() >= lp.numCol_);
    assert((int)solution.row_dual.size() >= lp.numRow_);
    num_dual_infeasibility = 0;
    max_dual_infeasibility = 0;
    sum_dual_infeasibility = 0;
  }

  double primal_infeasibility;
  double dual_infeasibility;
  double lower;
  double upper;
  double value;
  double dual;
  for (HighsInt iVar = 0; iVar < lp.numCol_ + lp.numRow_; iVar++) {
    if (iVar < lp.numCol_) {
      HighsInt iCol = iVar;
      lower = lp.colLower_[iCol];
      upper = lp.colUpper_[iCol];
      value = solution.col_value[iCol];
      if (have_dual_solution) dual = solution.col_dual[iCol];
    } else {
      HighsInt iRow = iVar - lp.numCol_;
      lower = lp.rowLower_[iRow];
      upper = lp.rowUpper_[iRow];
      value = solution.row_value[iRow];
      if (have_dual_solution) dual = -solution.row_dual[iRow];
    }
    // Flip dual according to lp.sense_
    if (have_dual_solution) dual *= (HighsInt)lp.sense_;

    // @primal_infeasibility calculation
    primal_infeasibility = 0;
    if (value < lower - primal_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + primal_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > 0) {
      if (primal_infeasibility > primal_feasibility_tolerance)
        num_primal_infeasibility++;
      max_primal_infeasibility =
          std::max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibility += primal_infeasibility;
    }
    if (have_dual_solution) {
      double primal_residual =
          std::min(std::fabs(lower - value), std::fabs(value - upper));
      if (primal_residual <= primal_feasibility_tolerance) {
        // At a bound
        double middle = (lower + upper) * 0.5;
        if (lower < upper) {
          // Non-fixed variable
          if (value < middle) {
            // At lower
            dual_infeasibility = std::max(-dual, 0.);
          } else {
            // At Upper
            dual_infeasibility = std::max(dual, 0.);
          }
        } else {
          // Fixed variable
          dual_infeasibility = 0;
        }
      } else {
        // Off bounds (or free)
        dual_infeasibility = fabs(dual);
      }
      if (dual_infeasibility > dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  if (num_primal_infeasibility == 0) {
    solution_params.primal_status = kHighsPrimalDualStatusFeasiblePoint;
  } else {
    solution_params.primal_status = kHighsPrimalDualStatusInfeasiblePoint;
  }
  if (have_dual_solution) {
    if (num_dual_infeasibility == 0) {
      solution_params.dual_status = kHighsPrimalDualStatusFeasiblePoint;
    } else {
      solution_params.dual_status = kHighsPrimalDualStatusInfeasiblePoint;
    }
  }
}

double computeObjectiveValue(const HighsLp& lp, const HighsSolution& solution) {
  double objective_value = 0;
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++)
    objective_value += lp.colCost_[iCol] * solution.col_value[iCol];
  objective_value += lp.offset_;
  return objective_value;
}

// Refine any HighsBasisStatus::kNonbasic settings according to the LP
// and any solution values
void refineBasis(const HighsLp& lp, const HighsSolution& solution,
                 HighsBasis& basis) {
  assert(basis.valid);
  assert(isBasisRightSize(lp, basis));
  const bool have_highs_solution = isSolutionRightSize(lp, solution);

  const HighsInt num_col = lp.numCol_;
  const HighsInt num_row = lp.numRow_;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (basis.col_status[iCol] != HighsBasisStatus::kNonbasic) continue;
    const double lower = lp.colLower_[iCol];
    const double upper = lp.colUpper_[iCol];
    HighsBasisStatus status = HighsBasisStatus::kNonbasic;
    if (lower == upper) {
      status = HighsBasisStatus::kLower;
    } else if (!highs_isInfinity(-lower)) {
      if (!highs_isInfinity(upper)) {
        if (have_highs_solution) {
          if (solution.col_value[iCol] < 0.5 * (lower + upper)) {
            status = HighsBasisStatus::kLower;
          } else {
            status = HighsBasisStatus::kUpper;
          }
        } else {
          if (fabs(lower) < fabs(upper)) {
            status = HighsBasisStatus::kLower;
          } else {
            status = HighsBasisStatus::kUpper;
          }
        }
      } else {
        status = HighsBasisStatus::kLower;
      }
    } else if (!highs_isInfinity(upper)) {
      status = HighsBasisStatus::kUpper;
    } else {
      status = HighsBasisStatus::kZero;
    }
    assert(status != HighsBasisStatus::kNonbasic);
    basis.col_status[iCol] = status;
  }

  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (basis.row_status[iRow] != HighsBasisStatus::kNonbasic) continue;
    const double lower = lp.rowLower_[iRow];
    const double upper = lp.rowUpper_[iRow];
    HighsBasisStatus status = HighsBasisStatus::kNonbasic;
    if (lower == upper) {
      status = HighsBasisStatus::kLower;
    } else if (!highs_isInfinity(-lower)) {
      if (!highs_isInfinity(upper)) {
        if (have_highs_solution) {
          if (solution.row_value[iRow] < 0.5 * (lower + upper)) {
            status = HighsBasisStatus::kLower;
          } else {
            status = HighsBasisStatus::kUpper;
          }
        } else {
          if (fabs(lower) < fabs(upper)) {
            status = HighsBasisStatus::kLower;
          } else {
            status = HighsBasisStatus::kUpper;
          }
        }
      } else {
        status = HighsBasisStatus::kLower;
      }
    } else if (!highs_isInfinity(upper)) {
      status = HighsBasisStatus::kUpper;
    } else {
      status = HighsBasisStatus::kZero;
    }
    assert(status != HighsBasisStatus::kNonbasic);
    basis.row_status[iRow] = status;
  }
}

#ifdef IPX_ON
HighsStatus ipxSolutionToHighsSolution(
    const HighsLogOptions& log_options, const HighsLp& lp,
    const std::vector<double>& rhs, const std::vector<char>& constraint_type,
    const HighsInt ipx_num_col, const HighsInt ipx_num_row,
    const std::vector<double>& ipx_x, const std::vector<double>& ipx_slack_vars,
    // const std::vector<double>& ipx_y,
    HighsSolution& highs_solution) {
  // Resize the HighsSolution
  highs_solution.col_value.resize(lp.numCol_);
  highs_solution.row_value.resize(lp.numRow_);

  const std::vector<double>& ipx_col_value = ipx_x;
  const std::vector<double>& ipx_row_value = ipx_slack_vars;
  //  const std::vector<double>& ipx_col_dual = ipx_x;
  //  const std::vector<double>& ipx_row_dual = ipx_y;

  // Row activities are needed to set activity values of free rows -
  // which are ignored by IPX
  vector<double> row_activity;
  bool get_row_activities = ipx_num_row < lp.numRow_;
  if (get_row_activities) row_activity.assign(lp.numRow_, 0);
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    highs_solution.col_value[col] = ipx_col_value[col];
    if (get_row_activities) {
      // Accumulate row activities to assign value to free rows
      for (HighsInt el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
        HighsInt row = lp.Aindex_[el];
        row_activity[row] += highs_solution.col_value[col] * lp.Avalue_[el];
      }
    }
  }
  HighsInt ipx_row = 0;
  HighsInt ipx_slack = lp.numCol_;
  HighsInt num_boxed_rows = 0;
  for (HighsInt row = 0; row < lp.numRow_; row++) {
    double lower = lp.rowLower_[row];
    double upper = lp.rowUpper_[row];
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free row - removed by IPX so set it to its row activity
      highs_solution.row_value[row] = row_activity[row];
    } else {
      // Non-free row, so IPX will have it
      if ((lower > -kHighsInf && upper < kHighsInf) && (lower < upper)) {
        // Boxed row - look at its slack
        num_boxed_rows++;
        highs_solution.row_value[row] = ipx_col_value[ipx_slack];
        // Update the slack to be used for boxed rows
        ipx_slack++;
      } else {
        highs_solution.row_value[row] = rhs[ipx_row] - ipx_row_value[ipx_row];
      }
      // Update the IPX row index
      ipx_row++;
    }
  }
  assert(ipx_row == ipx_num_row);
  assert(ipx_slack == ipx_num_col);
  // Indicate that the primal, but not dual solution is known
  highs_solution.value_valid = true;
  highs_solution.dual_valid = false;
  return HighsStatus::kOk;
}

HighsStatus ipxBasicSolutionToHighsBasicSolution(
    const HighsLogOptions& log_options, const HighsLp& lp,
    const std::vector<double>& rhs, const std::vector<char>& constraint_type,
    const IpxSolution& ipx_solution, HighsBasis& highs_basis,
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
  if (get_row_activities) row_activity.assign(lp.numRow_, 0);
  HighsInt num_basic_variables = 0;
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    bool unrecognised = false;
    if (ipx_col_status[col] == ipx_basic) {
      // Column is basic
      highs_basis.col_status[col] = HighsBasisStatus::kBasic;
      highs_solution.col_value[col] = ipx_col_value[col];
      highs_solution.col_dual[col] = 0;
    } else if (ipx_col_status[col] == ipx_nonbasic_at_lb) {
      // Column is nonbasic at lower bound
      highs_basis.col_status[col] = HighsBasisStatus::kLower;
      highs_solution.col_value[col] = ipx_col_value[col];
      highs_solution.col_dual[col] = ipx_col_dual[col];
    } else if (ipx_col_status[col] == ipx_nonbasic_at_ub) {
      // Column is nonbasic at upper bound
      highs_basis.col_status[col] = HighsBasisStatus::kUpper;
      highs_solution.col_value[col] = ipx_col_value[col];
      highs_solution.col_dual[col] = ipx_col_dual[col];
    } else if (ipx_col_status[col] == ipx_superbasic) {
      // Column is superbasic
      highs_basis.col_status[col] = HighsBasisStatus::kZero;
      highs_solution.col_value[col] = ipx_col_value[col];
      highs_solution.col_dual[col] = ipx_col_dual[col];
    } else {
      unrecognised = true;
      highsLogDev(log_options, HighsLogType::kError,
                  "\nError in IPX conversion: Unrecognised value "
                  "ipx_col_status[%2" HIGHSINT_FORMAT
                  "] = "
                  "%" HIGHSINT_FORMAT "\n",
                  col, (HighsInt)ipx_col_status[col]);
    }
    if (unrecognised) {
      highsLogDev(log_options, HighsLogType::kError,
                  "Bounds [%11.4g, %11.4g]\n", lp.colLower_[col],
                  lp.colUpper_[col]);
      highsLogDev(log_options, HighsLogType::kError,
                  "Col %2" HIGHSINT_FORMAT " ipx_col_status[%2" HIGHSINT_FORMAT
                  "] = %2" HIGHSINT_FORMAT "; x[%2" HIGHSINT_FORMAT
                  "] = %11.4g; z[%2" HIGHSINT_FORMAT
                  "] = "
                  "%11.4g\n",
                  col, col, (HighsInt)ipx_col_status[col], col,
                  ipx_col_value[col], col, ipx_col_dual[col]);
      assert(!unrecognised);
      highsLogUser(log_options, HighsLogType::kError,
                   "Unrecognised ipx_col_status value from IPX\n");
      return HighsStatus::kError;
    }
    if (get_row_activities) {
      // Accumulate row activities to assign value to free rows
      for (HighsInt el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
        HighsInt row = lp.Aindex_[el];
        row_activity[row] += highs_solution.col_value[col] * lp.Avalue_[el];
      }
    }
    if (highs_basis.col_status[col] == HighsBasisStatus::kBasic)
      num_basic_variables++;
  }
  HighsInt ipx_row = 0;
  HighsInt ipx_slack = lp.numCol_;
  HighsInt num_boxed_rows = 0;
  HighsInt num_boxed_rows_basic = 0;
  HighsInt num_boxed_row_slacks_basic = 0;
  for (HighsInt row = 0; row < lp.numRow_; row++) {
    bool unrecognised = false;
    double lower = lp.rowLower_[row];
    double upper = lp.rowUpper_[row];
    HighsInt this_ipx_row = ipx_row;
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free row - removed by IPX so make it basic at its row activity
      highs_basis.row_status[row] = HighsBasisStatus::kBasic;
      highs_solution.row_value[row] = row_activity[row];
      highs_solution.row_dual[row] = 0;
    } else {
      // Non-free row, so IPX will have it
      if ((lower > -kHighsInf && upper < kHighsInf) && (lower < upper)) {
        // Boxed row - look at its slack
        num_boxed_rows++;
        double slack_value = ipx_col_value[ipx_slack];
        double slack_dual = ipx_col_dual[ipx_slack];
        double value = slack_value;
        double dual = -slack_dual;
        if (ipx_row_status[ipx_row] == ipx_basic) {
          // Row is basic
          num_boxed_rows_basic++;
          highs_basis.row_status[row] = HighsBasisStatus::kBasic;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = 0;
        } else if (ipx_col_status[ipx_slack] == ipx_basic) {
          // Slack is basic
          num_boxed_row_slacks_basic++;
          highs_basis.row_status[row] = HighsBasisStatus::kBasic;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = 0;
        } else if (ipx_col_status[ipx_slack] == ipx_nonbasic_at_lb) {
          // Slack at lower bound
          highs_basis.row_status[row] = HighsBasisStatus::kLower;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else if (ipx_col_status[ipx_slack] == ipx_nonbasic_at_ub) {
          // Slack is at its upper bound
          assert(ipx_col_status[ipx_slack] == ipx_nonbasic_at_ub);
          highs_basis.row_status[row] = HighsBasisStatus::kUpper;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else {
          unrecognised = true;
          highsLogDev(log_options, HighsLogType::kError,
                      "Error in IPX conversion: Row %2" HIGHSINT_FORMAT
                      " (IPX row %2" HIGHSINT_FORMAT
                      ") has "
                      "unrecognised value ipx_col_status[%2" HIGHSINT_FORMAT
                      "] = %" HIGHSINT_FORMAT "\n",
                      row, ipx_row, ipx_slack,
                      (HighsInt)ipx_col_status[ipx_slack]);
        }
        // Update the slack to be used for boxed rows
        ipx_slack++;
      } else if (ipx_row_status[ipx_row] == ipx_basic) {
        // Row is basic
        highs_basis.row_status[row] = HighsBasisStatus::kBasic;
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
          highs_basis.row_status[row] = HighsBasisStatus::kLower;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else if (constraint_type[ipx_row] == '<') {
          // Row is at its upper bound
          highs_basis.row_status[row] = HighsBasisStatus::kUpper;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else if (constraint_type[ipx_row] == '=') {
          // Row is at its fixed value
          highs_basis.row_status[row] = HighsBasisStatus::kLower;
          highs_solution.row_value[row] = value;
          highs_solution.row_dual[row] = dual;
        } else {
          unrecognised = true;
          highsLogDev(log_options, HighsLogType::kError,
                      "Error in IPX conversion: Row %2" HIGHSINT_FORMAT
                      ": cannot handle "
                      "constraint_type[%2" HIGHSINT_FORMAT
                      "] = %" HIGHSINT_FORMAT "\n",
                      row, ipx_row, constraint_type[ipx_row]);
        }
      }
      // Update the IPX row index
      ipx_row++;
    }
    if (unrecognised) {
      highsLogDev(log_options, HighsLogType::kError,
                  "Bounds [%11.4g, %11.4g]\n", lp.rowLower_[row],
                  lp.rowUpper_[row]);
      highsLogDev(log_options, HighsLogType::kError,
                  "Row %2" HIGHSINT_FORMAT " ipx_row_status[%2" HIGHSINT_FORMAT
                  "] = %2" HIGHSINT_FORMAT "; s[%2" HIGHSINT_FORMAT
                  "] = %11.4g; y[%2" HIGHSINT_FORMAT
                  "] = "
                  "%11.4g\n",
                  row, this_ipx_row, (HighsInt)ipx_row_status[this_ipx_row],
                  this_ipx_row, ipx_row_value[this_ipx_row], this_ipx_row,
                  ipx_row_dual[this_ipx_row]);
      assert(!unrecognised);
      highsLogUser(log_options, HighsLogType::kError,
                   "Unrecognised ipx_row_status value from IPX\n");
      return HighsStatus::kError;
    }
    if (highs_basis.row_status[row] == HighsBasisStatus::kBasic)
      num_basic_variables++;
  }
  assert(num_basic_variables == lp.numRow_);
  assert(ipx_row == ipx_solution.num_row);
  assert(ipx_slack == ipx_solution.num_col);

  // Flip dual according to lp.sense_
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    highs_solution.col_dual[iCol] *= (HighsInt)lp.sense_;
  }
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    highs_solution.row_dual[iRow] *= (HighsInt)lp.sense_;
  }

  if (num_boxed_rows)
    highsLogDev(log_options, HighsLogType::kInfo,
                "Of %" HIGHSINT_FORMAT " boxed rows: %" HIGHSINT_FORMAT
                " are basic and %" HIGHSINT_FORMAT " have basic slacks\n",
                num_boxed_rows, num_boxed_rows_basic,
                num_boxed_row_slacks_basic);
  // Indicate that the primal solution, dual solution and basis are valid
  highs_solution.value_valid = true;
  highs_solution.dual_valid = true;
  highs_basis.valid = true;
  return HighsStatus::kOk;
}
#endif

std::string iterationsToString(const HighsIterationCounts& iterations_counts) {
  std::string iteration_statement = "";
  bool not_first = false;
  HighsInt num_positive_count = 0;
  if (iterations_counts.simplex) num_positive_count++;
  if (iterations_counts.ipm) num_positive_count++;
  if (iterations_counts.crossover) num_positive_count++;
  if (num_positive_count == 0) {
    iteration_statement += "0 iterations; ";
    return iteration_statement;
  }
  if (num_positive_count > 1) iteration_statement += "(";
  HighsInt count;
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
  resetModelStatusAndSolutionParams(highs_model_object.unscaled_model_status_,
                                    highs_model_object.solution_params_,
                                    highs_model_object.options_);
}

void resetModelStatusAndSolutionParams(HighsModelStatus& model_status,
                                       HighsSolutionParams& solution_params,
                                       const HighsOptions& options) {
  model_status = HighsModelStatus::kNotset;
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
  solution_params.primal_status = kHighsPrimalDualStatusNotset;
  solution_params.dual_status = kHighsPrimalDualStatusNotset;
}

// Invalidate the infeasibility values in a HighsSolutionParams
// instance.
void invalidateSolutionInfeasibilityParams(
    HighsSolutionParams& solution_params) {
  solution_params.num_primal_infeasibility = kHighsIllegalInfeasibilityCount;
  solution_params.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  solution_params.sum_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  solution_params.num_dual_infeasibility = kHighsIllegalInfeasibilityCount;
  solution_params.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  solution_params.sum_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
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
      solution_params.num_primal_infeasibility;
  highs_info.max_primal_infeasibility =
      solution_params.max_primal_infeasibility;
  highs_info.sum_primal_infeasibilities =
      solution_params.sum_primal_infeasibility;
  highs_info.num_dual_infeasibilities = solution_params.num_dual_infeasibility;
  highs_info.max_dual_infeasibility = solution_params.max_dual_infeasibility;
  highs_info.sum_dual_infeasibilities = solution_params.sum_dual_infeasibility;
}

bool isBasisConsistent(const HighsLp& lp, const HighsBasis& basis) {
  bool consistent = true;
  consistent = isBasisRightSize(lp, basis) && consistent;
  HighsInt num_basic_variables = 0;
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    if (basis.col_status[iCol] == HighsBasisStatus::kBasic)
      num_basic_variables++;
  }
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    if (basis.row_status[iRow] == HighsBasisStatus::kBasic)
      num_basic_variables++;
  }
  bool right_num_basic_variables = num_basic_variables == lp.numRow_;
  consistent = right_num_basic_variables && consistent;
  return consistent;
}

bool isPrimalSolutionRightSize(const HighsLp& lp,
                               const HighsSolution& solution) {
  return (HighsInt)solution.col_value.size() == lp.numCol_ &&
         (HighsInt)solution.row_value.size() == lp.numRow_;
}

bool isDualSolutionRightSize(const HighsLp& lp, const HighsSolution& solution) {
  return (HighsInt)solution.col_dual.size() == lp.numCol_ &&
         (HighsInt)solution.row_dual.size() == lp.numRow_;
}

bool isSolutionRightSize(const HighsLp& lp, const HighsSolution& solution) {
  return isPrimalSolutionRightSize(lp, solution) &&
         isDualSolutionRightSize(lp, solution);
}

bool isBasisRightSize(const HighsLp& lp, const HighsBasis& basis) {
  return (HighsInt)basis.col_status.size() == lp.numCol_ &&
         (HighsInt)basis.row_status.size() == lp.numRow_;
}

void clearPrimalSolutionUtil(HighsSolution& solution) {
  solution.col_value.clear();
  solution.row_value.clear();
  solution.value_valid = false;
}

void clearDualSolutionUtil(HighsSolution& solution) {
  solution.col_dual.clear();
  solution.row_dual.clear();
  solution.dual_valid = false;
}

void clearSolutionUtil(HighsSolution& solution) {
  clearPrimalSolutionUtil(solution);
  clearDualSolutionUtil(solution);
}

void clearBasisUtil(HighsBasis& basis) {
  basis.row_status.clear();
  basis.col_status.clear();
  basis.valid = false;
}
