/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
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
#include "ipm/ipx/ipx_status.h"
#include "ipm/ipx/lp_solver.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolutionDebug.h"

const uint8_t kHighsSolutionLo = -1;
const uint8_t kHighsSolutionNo = 0;
const uint8_t kHighsSolutionUp = 1;

const bool printf_kkt = false;  // true;  //

void getKktFailures(const HighsOptions& options, const HighsModel& model,
                    const HighsSolution& solution, const HighsBasis& basis,
                    HighsInfo& highs_info) {
  HighsPrimalDualErrors primal_dual_errors;
  getKktFailures(options, model, solution, basis, highs_info,
                 primal_dual_errors);
}

void getKktFailures(const HighsOptions& options, const HighsModel& model,
                    const HighsSolution& solution, const HighsBasis& basis,
                    HighsInfo& highs_info,
                    HighsPrimalDualErrors& primal_dual_errors,
                    const bool get_residuals) {
  vector<double> gradient;
  model.objectiveGradient(solution.col_value, gradient);
  const HighsLp& lp = model.lp_;
  getKktFailures(options, model.isQp(), lp, gradient, solution, highs_info,
                 get_residuals);
  getPrimalDualBasisErrors(options, lp, solution, basis, primal_dual_errors);
  getPrimalDualGlpsolErrors(options, lp, gradient, solution,
                            primal_dual_errors);
}

void getLpKktFailures(const HighsOptions& options, const HighsLp& lp,
                      const HighsSolution& solution, const HighsBasis& basis,
                      HighsInfo& highs_info) {
  HighsPrimalDualErrors primal_dual_errors;
  getLpKktFailures(options, lp, solution, basis, highs_info,
                   primal_dual_errors);
}

void getLpKktFailures(const HighsOptions& options, const HighsLp& lp,
                      const HighsSolution& solution, const HighsBasis& basis,
                      HighsInfo& highs_info,
                      HighsPrimalDualErrors& primal_dual_errors,
                      const bool get_residuals) {
  getKktFailures(options, false, lp, lp.col_cost_, solution, highs_info,
                 get_residuals);
  getPrimalDualBasisErrors(options, lp, solution, basis, primal_dual_errors);
  getPrimalDualGlpsolErrors(options, lp, lp.col_cost_, solution,
                            primal_dual_errors);
}

void getKktFailures(const HighsOptions& options, const bool is_qp,
                    const HighsLp& lp, const std::vector<double>& gradient,
                    const HighsSolution& solution, HighsInfo& highs_info,
                    const bool get_residuals) {
  double primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  double dual_feasibility_tolerance = options.dual_feasibility_tolerance;
  double primal_residual_tolerance = options.primal_residual_tolerance;
  double dual_residual_tolerance = options.dual_residual_tolerance;
  double optimality_tolerance = options.optimality_tolerance;
  if (options.kkt_tolerance != kDefaultKktTolerance) {
    primal_feasibility_tolerance = options.kkt_tolerance;
    dual_feasibility_tolerance = options.kkt_tolerance;
    primal_residual_tolerance = options.kkt_tolerance;
    dual_residual_tolerance = options.kkt_tolerance;
    optimality_tolerance = options.kkt_tolerance;
  }
  // highs_info are the values computed in this method

  HighsInt& num_primal_infeasibility = highs_info.num_primal_infeasibilities;
  double& max_primal_infeasibility = highs_info.max_primal_infeasibility;
  double& sum_primal_infeasibility = highs_info.sum_primal_infeasibilities;

  HighsInt& num_dual_infeasibility = highs_info.num_dual_infeasibilities;
  double& max_dual_infeasibility = highs_info.max_dual_infeasibility;
  double& sum_dual_infeasibility = highs_info.sum_dual_infeasibilities;

  HighsInt& num_relative_primal_infeasibility =
      highs_info.num_relative_primal_infeasibilities;
  double& max_relative_primal_infeasibility =
      highs_info.max_relative_primal_infeasibility;

  HighsInt& num_relative_dual_infeasibility =
      highs_info.num_relative_dual_infeasibilities;
  double& max_relative_dual_infeasibility =
      highs_info.max_relative_dual_infeasibility;

  HighsInt& num_primal_residual_error = highs_info.num_primal_residual_errors;
  double& max_primal_residual_error = highs_info.max_primal_residual_error;

  HighsInt& num_dual_residual_error = highs_info.num_dual_residual_errors;
  double& max_dual_residual_error = highs_info.max_dual_residual_error;

  HighsInt& num_relative_primal_residual_error =
      highs_info.num_relative_primal_residual_errors;
  double& max_relative_primal_residual_error =
      highs_info.max_relative_primal_residual_error;

  HighsInt& num_relative_dual_residual_error =
      highs_info.num_relative_dual_residual_errors;
  double& max_relative_dual_residual_error =
      highs_info.max_relative_dual_residual_error;

  HighsInt& num_complementarity_violation =
      highs_info.num_complementarity_violations;
  double& max_complementarity_violation =
      highs_info.max_complementarity_violation;

  double& primal_dual_objective_error = highs_info.primal_dual_objective_error;

  const bool& have_primal_solution = solution.value_valid;
  const bool& have_dual_solution = solution.dual_valid;
  const bool have_integrality = (lp.integrality_.size() != 0);

  // Primal/dual solution status on entry can be feasible or
  // infeasible, because it refers to the previous problem solved, so
  // can't check whether something meaningful is being changed
  highs_info.primal_solution_status = kSolutionStatusNone;
  highs_info.dual_solution_status = kSolutionStatusNone;

  // Check that there is no dual solution if there's no primal solution
  assert(have_primal_solution || !have_dual_solution);

  // Invalidate all the KKT measures
  highs_info.invalidateKkt();

  if (have_primal_solution) {
    // There's a primal solution, so check its size and initialise the
    // infeasibility counts
    assert((int)solution.col_value.size() >= lp.num_col_);
    assert((int)solution.row_value.size() >= lp.num_row_);
    num_primal_infeasibility = 0;
    max_primal_infeasibility = 0;
    sum_primal_infeasibility = 0;
    num_relative_primal_infeasibility = 0;
    max_relative_primal_infeasibility = 0;
    if (get_residuals) {
      num_primal_residual_error = 0;
      max_primal_residual_error = 0;
      num_relative_primal_residual_error = 0;
      max_relative_primal_residual_error = 0;
    }
  }
  if (have_dual_solution) {
    // There's a dual solution, so check its size and initialise the
    // infeasibility counts
    assert((int)solution.col_dual.size() >= lp.num_col_);
    assert((int)solution.row_dual.size() >= lp.num_row_);
    num_dual_infeasibility = 0;
    max_dual_infeasibility = 0;
    sum_dual_infeasibility = 0;
    num_relative_dual_infeasibility = 0;
    max_relative_dual_infeasibility = 0;
    if (get_residuals) {
      num_dual_residual_error = 0;
      max_dual_residual_error = 0;
      num_relative_dual_residual_error = 0;
      max_relative_dual_residual_error = 0;
    }
  }
  // Without a primal solution, nothing can be done!
  if (!have_primal_solution) return;

  // Possibly compute Ax and A^Ty (to form A^Ty-c) as residual check
  // for solution.row_value and -solution.col_dual
  std::vector<double> primal_activity;
  std::vector<double> dual_activity;
  if (get_residuals)
    lp.a_matrix_.productQuad(primal_activity, solution.col_value);
  if (get_residuals && have_dual_solution)
    lp.a_matrix_.productTransposeQuad(dual_activity, solution.row_dual);

  double max_col_primal_infeasibility = 0;
  double max_col_dual_infeasibility = 0;
  double max_relative_col_primal_infeasibility = 0;
  double max_relative_col_dual_infeasibility = 0;

  double primal_infeasibility;
  double relative_primal_infeasibility;
  double dual_infeasibility;
  double cost = 0.0;
  double lower;
  double upper;
  double value;
  double dual = 0;
  uint8_t at_status;
  uint8_t mid_status;
  HighsVarType integrality = HighsVarType::kContinuous;
  HighsBasisStatus status;
  // Compute the infinity norm of all near-active column and row
  // bounds, since they contribute to the magnitude of the row values,
  // and dividing their residual error by the norm gives a relative
  // residual error: don't consider large inactive bounds, since they
  // don't affect the model
  //
  // In IPM without crossover, the columns and rows with values close
  // to their bounds are marginal in driving the algorithm so, by
  // computing the norm as we do, we capture the active bounds of all
  // these variables.
  //
  // In PDLP, only the constraint bounds contribute to the norm, but
  // do the variable bounds not (also) drive the algorithm
  //
  // In simplex the values of the nonbasic variables would be used,
  // since they determine the RHS of the system solved for the values
  // of the basic variables values. That said, by computing the norm
  // as we do, we capture the active bounds of all the nonbasic variables
  double highs_norm_bounds = 0.0;
  // Compute the infinity norm of all near-active column duals, since
  // they contribute to the magnitude of the row values, and dividing
  // their residual error by the norm gives a relative residual error:
  // don't consider large inactive bounds, since they don't affect the
  // model
  double highs_norm_costs = 0.0;

  // Pass twice through this loop, once to determine the bound and
  // cost norms, and once to use them to assess relative
  // infeasibilities and residual errors
  for (HighsInt pass = 0; pass < 2; pass++) {
    for (HighsInt iVar = 0; iVar < lp.num_col_ + lp.num_row_; iVar++) {
      const bool is_col = iVar < lp.num_col_;
      if (is_col) {
        HighsInt iCol = iVar;
        cost = gradient[iCol];
        lower = lp.col_lower_[iCol];
        upper = lp.col_upper_[iCol];
        value = solution.col_value[iCol];
        if (have_dual_solution) dual = solution.col_dual[iCol];
        if (have_integrality) integrality = lp.integrality_[iCol];
        if (pass == 0) {
          if (dual * dual < dual_feasibility_tolerance) {
            // Dual close to zero
            highs_norm_costs = std::max(std::fabs(cost), highs_norm_costs);
          }
          if (get_residuals && have_dual_solution) {
            // Subtract off the gradient value
            HighsCDouble q_dual_activity = dual_activity[iCol];
            q_dual_activity -= gradient[iCol];
            dual_activity[iCol] = double(q_dual_activity);
          }
        }
        //
      } else {
        HighsInt iRow = iVar - lp.num_col_;
        lower = lp.row_lower_[iRow];
        upper = lp.row_upper_[iRow];
        value = solution.row_value[iRow];
        if (have_dual_solution) dual = solution.row_dual[iRow];
        integrality = HighsVarType::kContinuous;
      }
      // Flip dual according to lp.sense_
      dual *= (HighsInt)lp.sense_;
      // Get the primal and dual infeasibility for this variable, and
      //
      // at_status: Indicates whether the variable is close to its
      // lower bound, upper bound or not at all. Use this as a proxy
      // for being non-basic, so the active bound contributes to
      // highs_norm_bounds for calculating relative primal measures.
      //
      // mid_status: Indicates whether a variable with meaningful
      // bound interval length is below the midpoint of its bound
      // interval, or above the midpoint. The midpoint is infinite for
      // one-sided variables. For free variables, fixed variables, or
      // variables with small positive bound interval lengths,
      // mid_status is returned as kHighsSolutionNo.
      getVariableKktFailures(primal_feasibility_tolerance,
                             dual_feasibility_tolerance, lower, upper, value,
                             dual, integrality, primal_infeasibility,
                             dual_infeasibility, at_status, mid_status);
      if (pass == 0) {
        // If the primal value is close to a bound then include the bound
        // in the active bound norm
        if (at_status == kHighsSolutionLo) {
          highs_norm_bounds = std::max(std::fabs(lower), highs_norm_bounds);
        } else if (at_status == kHighsSolutionUp) {
          highs_norm_bounds = std::max(std::fabs(upper), highs_norm_bounds);
        }
      } else {
        if (primal_infeasibility > 0) {
          // Accumulate primal infeasibilities
          if (primal_infeasibility > primal_feasibility_tolerance)
            num_primal_infeasibility++;
          if (max_primal_infeasibility < primal_infeasibility)
            max_primal_infeasibility = primal_infeasibility;
          sum_primal_infeasibility += primal_infeasibility;
          // Determine the denominator for the relative primal
          // infeasiblilty
          double relative_bound_measure = highs_norm_bounds;
          if (at_status == kHighsSolutionNo) {
            // Primal value is infeasible, but not close to a bound:
            // unusual, but possible if absolute primal infeasiblities
            // are not small. Bound has not been included in
            // highs_norm_bounds, but should be used for local
            // relative infeasiblility
            if (mid_status == kHighsSolutionNo ||
                mid_status == kHighsSolutionLo) {
              // Bound interval is short, or variable is below the
              // midpoint (which is only possible if lower is finite)
              relative_bound_measure =
                  std::max(std::fabs(lower), relative_bound_measure);
            } else {
              assert(mid_status == kHighsSolutionUp);
              // Variable is above the midpoint of the bound interval
              // (which is only possible if upper is finite)
              relative_bound_measure =
                  std::max(std::fabs(upper), relative_bound_measure);
            }
          }
          double relative_primal_infeasibility =
              primal_infeasibility / (1.0 + relative_bound_measure);

          if (relative_primal_infeasibility > primal_feasibility_tolerance)
            num_relative_primal_infeasibility++;
          if (max_relative_primal_infeasibility < relative_primal_infeasibility)
            max_relative_primal_infeasibility = relative_primal_infeasibility;
        }
        if (have_dual_solution) {
          if (dual_infeasibility > 0) {
            // Accumulate dual infeasibilities
            if (dual_infeasibility > dual_feasibility_tolerance)
              num_dual_infeasibility++;
            if (max_dual_infeasibility < dual_infeasibility)
              max_dual_infeasibility = dual_infeasibility;
            sum_dual_infeasibility += dual_infeasibility;

            // Determine the denominator for the relative dual
            // infeasiblilty
            double relative_cost_measure = highs_norm_costs;
            if (is_col && cost && dual * dual >= dual_feasibility_tolerance) {
              // Dual value is infeasible, but not close to zero:
              // unusual, but possible if absolute dual infeasiblities
              // are not small. Hence the cost has not been included in
              // highs_norm_costs, but should be used for local relative
              // infeasiblility.
              //
              // updateRelativeMeasure(cost, relative_cost_measure);
              relative_cost_measure =
                  std::max(std::fabs(cost), relative_cost_measure);
            }
            double relative_dual_infeasibility =
                dual_infeasibility / (1.0 + relative_cost_measure);

            if (relative_dual_infeasibility > dual_feasibility_tolerance)
              num_relative_dual_infeasibility++;
            if (max_relative_dual_infeasibility < relative_dual_infeasibility)
              max_relative_dual_infeasibility = relative_dual_infeasibility;
          }
        }
        if (!is_col && get_residuals) {
          HighsInt iRow = iVar - lp.num_col_;
          assert(iRow >= 0);
          double primal_residual_error =
              std::fabs(primal_activity[iRow] - solution.row_value[iRow]);
          double relative_primal_residual_error =
              primal_residual_error / (1.0 + highs_norm_bounds);

          if (primal_residual_error > primal_residual_tolerance)
            num_primal_residual_error++;
          if (max_primal_residual_error < primal_residual_error)
            max_primal_residual_error = primal_residual_error;

          if (relative_primal_residual_error > primal_residual_tolerance)
            num_relative_primal_residual_error++;
          if (max_relative_primal_residual_error <
              relative_primal_residual_error)
            max_relative_primal_residual_error = relative_primal_residual_error;
        }
        if (is_col && get_residuals && have_dual_solution) {
          HighsInt iCol = iVar;
          assert(iCol < lp.num_col_);
          double dual_residual_error =
              std::fabs(dual_activity[iCol] + solution.col_dual[iCol]);
          double relative_dual_residual_error =
              dual_residual_error / (1.0 + highs_norm_costs);

          if (dual_residual_error > dual_residual_tolerance)
            num_dual_residual_error++;
          if (max_dual_residual_error < dual_residual_error)
            max_dual_residual_error = dual_residual_error;

          if (relative_dual_residual_error > dual_residual_tolerance)
            num_relative_dual_residual_error++;
          if (max_relative_dual_residual_error < relative_dual_residual_error)
            max_relative_dual_residual_error = relative_dual_residual_error;
        }
      }
      if (pass == 1 && iVar == lp.num_col_ - 1) {
        max_col_primal_infeasibility = max_primal_infeasibility;
        max_col_dual_infeasibility = max_dual_infeasibility;

        max_relative_col_primal_infeasibility =
            max_relative_primal_infeasibility;
        max_relative_col_dual_infeasibility = max_relative_dual_infeasibility;

        max_primal_infeasibility = 0;
        max_dual_infeasibility = 0;

        max_relative_primal_infeasibility = 0;
        max_relative_dual_infeasibility = 0;
      }
    }
  }

  double max_row_primal_infeasibility = max_primal_infeasibility;
  double max_row_dual_infeasibility = max_dual_infeasibility;

  double max_relative_row_primal_infeasibility =
      max_relative_primal_infeasibility;
  double max_relative_row_dual_infeasibility = max_relative_dual_infeasibility;

  max_primal_infeasibility =
      std::max(max_col_primal_infeasibility, max_row_primal_infeasibility);
  max_dual_infeasibility =
      std::max(max_col_dual_infeasibility, max_row_dual_infeasibility);

  max_relative_primal_infeasibility =
      std::max(max_relative_col_primal_infeasibility,
               max_relative_row_primal_infeasibility);
  max_relative_dual_infeasibility = std::max(
      max_relative_col_dual_infeasibility, max_relative_row_dual_infeasibility);

  if (have_dual_solution) {
    // Determine the sum of complementarity violations
    const bool have_values = getComplementarityViolations(
        lp, solution, optimality_tolerance, num_complementarity_violation,
        max_complementarity_violation);
    assert(have_values);
  }

  if (have_dual_solution) {
    // Determine the primal-dual objective error
    //
    // IPX computes objective_gap = (pobjective-dobjective) / (1.0 +
    // 0.5*std::fabs(pobjective+dobjective));
    //
    // PDLP computes dRelObjGap = fabs(dPrimalObj - dDualObj) / (1.0 +
    // fabs(dPrimalObj) + fabs(dDualObj));
    //
    // Use the PDLP relative primal-dual objective error
    //
    // The dual objective for a QP has a -(1/2)x^TQx term, and this
    // can be computed from the gradient (g = Qx + c) as
    // -(1/2)(g-c)^Tx = (1/2)(c-g)^Tx, so pass a pointer to the
    // gradient data if this is necessary
    const double* gradient_p = is_qp ? gradient.data() : nullptr;
    double dual_objective_value;
    bool dual_objective_status = computeDualObjectiveValue(
        gradient_p, lp, solution, dual_objective_value);
    assert(dual_objective_status);
    const double abs_objective_difference =
        std::fabs(highs_info.objective_function_value - dual_objective_value);
    primal_dual_objective_error = abs_objective_difference;
    const double denominator = 1.0 +
                               std::fabs(highs_info.objective_function_value) +
                               std::fabs(dual_objective_value);
    primal_dual_objective_error = primal_dual_objective_error / denominator;
  }

  if (printf_kkt || options.log_dev_level > 0) {
    printf("getKktFailures:: cost norm = %8.3g; bound norm = %8.3g\n",
           highs_norm_costs, highs_norm_bounds);
    printf(
        "getKktFailures:                      LP  (abs / rel)    "
        "     Col (abs / rel)         Row (abs / rel)\n");

    printf(
        "getKktFailures: primal infeasibility %8.3g / %8.3g     %8.3g / "
        "%8.3g     %8.3g / %8.3g\n",
        highs_info.max_primal_infeasibility,
        highs_info.max_relative_primal_infeasibility,
        max_col_primal_infeasibility, max_relative_col_primal_infeasibility,
        max_row_primal_infeasibility, max_relative_row_primal_infeasibility);
    if (have_dual_solution)
      printf(
          "getKktFailures:   dual infeasibility %8.3g / %8.3g     %8.3g / "
          "%8.3g     %8.3g / %8.3g\n",
          highs_info.max_dual_infeasibility,
          highs_info.max_relative_dual_infeasibility,
          max_col_dual_infeasibility, max_relative_col_dual_infeasibility,
          max_row_dual_infeasibility, max_relative_row_dual_infeasibility);
    if (get_residuals) {
      printf("getKktFailures: primal residual      %8.3g / %8.3g\n",
             highs_info.max_primal_residual_error,
             highs_info.max_relative_primal_residual_error);
      if (have_dual_solution)
        printf("getKktFailures:   dual residual      %8.3g / %8.3g\n",
               highs_info.max_dual_residual_error,
               highs_info.max_relative_dual_residual_error);
    }
    if (!is_qp && have_dual_solution)
      printf("getKktFailures: objective gap        %8.3g\n",
             highs_info.primal_dual_objective_error);
  }

  // Assign primal (and possibly dual) solution status according to
  // existence of primal and dual feasibilities
  //
  // For LP this may mean that primal/dual feasibility of the solution
  // is claimed when there are residual errors, or denied when PDLP or
  // IPX without crossover are used, since they are deemed feasible
  // according to relative tolerances. This will be cleaned up in
  // Highs::lpKktCheck, when the existence of a basis determines
  // whether absolute or relative measures are used.

  if (num_primal_infeasibility) {
    highs_info.primal_solution_status = kSolutionStatusInfeasible;
  } else {
    highs_info.primal_solution_status = kSolutionStatusFeasible;
  }
  if (have_dual_solution) {
    if (num_dual_infeasibility) {
      highs_info.dual_solution_status = kSolutionStatusInfeasible;
    } else {
      highs_info.dual_solution_status = kSolutionStatusFeasible;
    }
  }
}

// Gets the KKT failures for a variable.
//
// Value and dual are used compute the primal and dual infeasibility
// It's up to the calling method to ignore these if the value or dual
// are not valid.
void getVariableKktFailures(const double primal_feasibility_tolerance,
                            const double dual_feasibility_tolerance,
                            const double lower, const double upper,
                            const double value, const double dual,
                            const HighsVarType integrality,
                            double& primal_infeasibility,
                            double& dual_infeasibility, uint8_t& at_status,
                            uint8_t& mid_status) {
  const HighsInt feasibility_tolerance_mu = 0.0;
  // @primal_infeasibility calculation
  primal_infeasibility = 0;
  if (value < lower - primal_feasibility_tolerance * feasibility_tolerance_mu) {
    // Below lower
    primal_infeasibility = lower - value;
  } else if (value >
             upper + primal_feasibility_tolerance * feasibility_tolerance_mu) {
    // Above upper
    primal_infeasibility = value - upper;
  }
  // Determine whether this value is close to a bound
  at_status = kHighsSolutionNo;
  double residual = std::fabs(lower - value);
  if (residual * residual <= primal_feasibility_tolerance) {
    // Close to lower bound
    at_status = kHighsSolutionLo;
  } else {
    // Not close to lower bound: maybe close to upper bound
    residual = std::fabs(value - upper);
    if (residual * residual <= primal_feasibility_tolerance)
      at_status = kHighsSolutionUp;
  }
  // Look for dual sign errors that exceed the tolerance. For boxed
  // variables the test is discontinuous at the midpoint, but any
  // meaningful dual value on a meaningful interval will show up as a
  // large complementarity error
  mid_status = kHighsSolutionNo;
  if (lower < upper) {
    double length = upper - lower;
    // Non-fixed variable
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free variable
      dual_infeasibility = fabs(dual);
    } else if (length * length > primal_feasibility_tolerance) {
      // Interval length is meaningful
      //
      // Compute the mid-point of the bound interval. This will be
      // +infinity for LB variables; -infinity for UB variables;
      // finite for boxed and fixed variables
      const double middle = (lower + upper) * 0.5;
      if (value < middle) {
        // Below the mid-point, so use lower bound optimality condition:
        // feasiblity is dual >= -tolerance
        mid_status = kHighsSolutionLo;
        dual_infeasibility = std::max(-dual, 0.);
      } else {
        // Below the mid-point, so use upper bound optimality condition:
        // feasiblity is dual <= tolerance
        mid_status = kHighsSolutionUp;
        dual_infeasibility = std::max(dual, 0.);
      }
    } else {
      // Interval length is less than
      // sqrt(primal_feasibility_tolerance) so dual infeasibility is
      // hard to define
      dual_infeasibility = 0;
    }
  } else {
    // Fixed variable
    dual_infeasibility = 0;
  }
  // Account for semi-variables
  const bool semi_variable = integrality == HighsVarType::kSemiContinuous ||
                             integrality == HighsVarType::kSemiInteger;
  if (semi_variable && std::fabs(value) < primal_feasibility_tolerance)
    primal_infeasibility = 0;
}

void getPrimalDualGlpsolErrors(const HighsOptions& options, const HighsLp& lp,
                               const std::vector<double>& gradient,
                               const HighsSolution& solution,
                               HighsPrimalDualErrors& primal_dual_errors) {
  double primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  double dual_feasibility_tolerance = options.dual_feasibility_tolerance;
  double primal_residual_tolerance = options.primal_residual_tolerance;
  double dual_residual_tolerance = options.dual_residual_tolerance;
  double optimality_tolerance = options.optimality_tolerance;

  // clang-format off
  HighsInt& num_primal_residual_error = primal_dual_errors.glpsol_num_primal_residual_errors;
  double&   sum_primal_residual_error = primal_dual_errors.glpsol_sum_primal_residual_errors;
  
  HighsInt& num_dual_residual_error = primal_dual_errors.glpsol_num_dual_residual_errors;
  double&   sum_dual_residual_error = primal_dual_errors.glpsol_sum_dual_residual_errors;
  
  double&   max_primal_residual_error = primal_dual_errors.glpsol_max_primal_residual.absolute_value;
  HighsInt& max_primal_residual_index = primal_dual_errors.glpsol_max_primal_residual.absolute_index;

  double&   max_relative_primal_residual_error = primal_dual_errors.glpsol_max_primal_residual.relative_value;
  HighsInt& max_relative_primal_residual_index = primal_dual_errors.glpsol_max_primal_residual.relative_index;

  double&   max_primal_infeasibility =       primal_dual_errors.glpsol_max_primal_infeasibility.absolute_value;
  HighsInt& max_primal_infeasibility_index = primal_dual_errors.glpsol_max_primal_infeasibility.absolute_index;

  double&   max_relative_primal_infeasibility =       primal_dual_errors.glpsol_max_primal_infeasibility.relative_value;
  HighsInt& max_relative_primal_infeasibility_index = primal_dual_errors.glpsol_max_primal_infeasibility.relative_index;

  double&   max_dual_residual_error = primal_dual_errors.glpsol_max_dual_residual.absolute_value;
  HighsInt& max_dual_residual_index = primal_dual_errors.glpsol_max_dual_residual.absolute_index;

  double&   max_relative_dual_residual_error = primal_dual_errors.glpsol_max_dual_residual.relative_value;
  HighsInt& max_relative_dual_residual_index = primal_dual_errors.glpsol_max_dual_residual.relative_index;

  double&   max_dual_infeasibility =       primal_dual_errors.glpsol_max_dual_infeasibility.absolute_value;
  HighsInt& max_dual_infeasibility_index = primal_dual_errors.glpsol_max_dual_infeasibility.absolute_index;
  // clang-format on

  // Relative dual infeasiblities are same as absolute

  primal_dual_errors.glpsol_max_primal_infeasibility.invalidate();
  primal_dual_errors.glpsol_max_dual_infeasibility.invalidate();

  const bool have_primal_solution = solution.value_valid;
  const bool have_dual_solution = solution.dual_valid;
  const bool have_integrality = lp.integrality_.size() != 0;

  // Check that there is no dual solution if there's no primal solution
  assert(have_primal_solution || !have_dual_solution);

  if (have_primal_solution) {
    // There's a primal solution, so check its size and initialise the
    // infeasibility counts
    assert((int)solution.col_value.size() >= lp.num_col_);
    assert((int)solution.row_value.size() >= lp.num_row_);
    max_primal_infeasibility = 0;
    primal_dual_errors.glpsol_max_primal_infeasibility.reset();
    if (have_dual_solution) {
      // There's a dual solution, so check its size and initialise the
      // infeasibility counts
      assert((int)solution.col_dual.size() >= lp.num_col_);
      assert((int)solution.row_dual.size() >= lp.num_row_);
      max_dual_infeasibility = 0;
      primal_dual_errors.glpsol_max_dual_infeasibility.reset();
    }
  }

  if (have_primal_solution) {
    num_primal_residual_error = 0;
    max_primal_residual_error = 0;
    max_relative_primal_residual_error = 0;
    primal_dual_errors.glpsol_max_primal_residual.reset();

  } else {
    num_primal_residual_error = kHighsIllegalResidualCount;
    max_primal_residual_error = kHighsIllegalResidualMeasure;
    max_relative_primal_residual_error = kHighsIllegalResidualMeasure;
    primal_dual_errors.glpsol_max_primal_residual.invalidate();
  }
  if (have_dual_solution) {
    num_dual_residual_error = 0;
    max_dual_residual_error = 0;
    max_relative_dual_residual_error = 0;
    primal_dual_errors.glpsol_max_dual_residual.reset();
  } else {
    num_dual_residual_error = kHighsIllegalResidualCount;
    max_dual_residual_error = kHighsIllegalResidualMeasure;
    max_relative_dual_residual_error = kHighsIllegalResidualMeasure;
    primal_dual_errors.glpsol_max_dual_residual.invalidate();
  }
  // Without a primal solution, nothing can be done!
  if (!have_primal_solution) return;

  // Residuals are formed for Glpsol via their positive and negative
  // terms so that meaningful relative values can be computed
  std::vector<double> primal_positive_sum;
  std::vector<double> primal_negative_sum;
  std::vector<double> dual_positive_sum;
  std::vector<double> dual_negative_sum;
  primal_positive_sum.assign(lp.num_row_, 0);
  primal_negative_sum.assign(lp.num_row_, 0);
  if (have_dual_solution) {
    dual_positive_sum.resize(lp.num_col_);
    dual_negative_sum.resize(lp.num_col_);
  }

  double primal_infeasibility;
  double relative_primal_infeasibility;
  double dual_infeasibility;
  double lower;
  double upper;
  double value;
  double dual = 0;
  uint8_t at_status;
  uint8_t mid_status;
  HighsVarType integrality = HighsVarType::kContinuous;
  for (HighsInt iVar = 0; iVar < lp.num_col_ + lp.num_row_; iVar++) {
    if (iVar < lp.num_col_) {
      HighsInt iCol = iVar;
      lower = lp.col_lower_[iCol];
      upper = lp.col_upper_[iCol];
      value = solution.col_value[iCol];
      if (have_dual_solution) dual = solution.col_dual[iCol];
      if (have_integrality) integrality = lp.integrality_[iCol];
    } else {
      HighsInt iRow = iVar - lp.num_col_;
      lower = lp.row_lower_[iRow];
      upper = lp.row_upper_[iRow];
      value = solution.row_value[iRow];
      if (have_dual_solution) dual = solution.row_dual[iRow];
      integrality = HighsVarType::kContinuous;
    }
    // Flip dual according to lp.sense_
    dual *= (HighsInt)lp.sense_;

    getVariableKktFailures(primal_feasibility_tolerance,
                           dual_feasibility_tolerance, lower, upper, value,
                           dual, integrality, primal_infeasibility,
                           dual_infeasibility, at_status, mid_status);

    relative_primal_infeasibility = 0;
    if (mid_status == kHighsSolutionLo) {
      relative_primal_infeasibility =
          primal_infeasibility / (1 + std::fabs(lower));
    } else if (mid_status == kHighsSolutionUp) {
      relative_primal_infeasibility =
          primal_infeasibility / (1 + std::fabs(upper));
    } else if (lower > -kHighsInf) {
      // Variable has small bound length
      relative_primal_infeasibility =
          primal_infeasibility / (1 + std::fabs(lower));
    }

    if (max_primal_infeasibility < primal_infeasibility) {
      max_primal_infeasibility = primal_infeasibility;
      max_primal_infeasibility_index = iVar;
    }
    if (max_relative_primal_infeasibility < relative_primal_infeasibility) {
      max_relative_primal_infeasibility = relative_primal_infeasibility;
      max_relative_primal_infeasibility_index = iVar;
    }

    if (have_dual_solution) {
      // Accumulate dual infeasibilities
      if (max_dual_infeasibility < dual_infeasibility) {
        max_dual_infeasibility = dual_infeasibility;
        max_dual_infeasibility_index = iVar;
      }
    }
    if (iVar < lp.num_col_) {
      HighsInt iCol = iVar;
      if (have_dual_solution) {
        if (gradient[iCol] > 0) {
          dual_positive_sum[iCol] = gradient[iCol];
        } else {
          dual_negative_sum[iCol] = -gradient[iCol];
        }
      }
      for (HighsInt el = lp.a_matrix_.start_[iCol];
           el < lp.a_matrix_.start_[iCol + 1]; el++) {
        HighsInt iRow = lp.a_matrix_.index_[el];
        double Avalue = lp.a_matrix_.value_[el];
        double term = value * Avalue;
        if (term > 0) {
          primal_positive_sum[iRow] += term;
        } else {
          primal_negative_sum[iRow] -= term;
        }
        // @FlipRowDual += became -=
        if (have_dual_solution) {
          double term = -solution.row_dual[iRow] * Avalue;
          if (term > 0) {
            dual_positive_sum[iCol] += term;
          } else {
            dual_negative_sum[iCol] -= term;
          }
        }
      }
    }
  }

  // Relative dual infeasiblities are same as absolute
  primal_dual_errors.glpsol_max_dual_infeasibility.relative_value =
      primal_dual_errors.glpsol_max_dual_infeasibility.absolute_value;

  primal_dual_errors.glpsol_max_dual_infeasibility.relative_index =
      primal_dual_errors.glpsol_max_dual_infeasibility.absolute_index;
}

void getPrimalDualBasisErrors(const HighsOptions& options, const HighsLp& lp,
                              const HighsSolution& solution,
                              const HighsBasis& basis,
                              HighsPrimalDualErrors& primal_dual_errors) {
  double primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  double dual_feasibility_tolerance = options.dual_feasibility_tolerance;
  const bool& have_primal_solution = solution.value_valid;
  const bool& have_dual_solution = solution.dual_valid;
  const bool& have_basis = basis.valid;

  // Check that there is no dual solution if there's no primal solution
  assert(have_primal_solution || !have_dual_solution);

  // Check that there is no basis if there's no dual solution
  assert(have_dual_solution || !have_basis);

  HighsInt& num_nonzero_basic_dual = primal_dual_errors.num_nonzero_basic_duals;
  double& max_nonzero_basic_dual = primal_dual_errors.max_nonzero_basic_dual;
  double& sum_nonzero_basic_dual = primal_dual_errors.sum_nonzero_basic_duals;

  HighsInt& num_off_bound_nonbasic = primal_dual_errors.num_off_bound_nonbasic;
  double& max_off_bound_nonbasic = primal_dual_errors.max_off_bound_nonbasic;
  double& sum_off_bound_nonbasic = primal_dual_errors.sum_off_bound_nonbasic;

  if (have_basis) {
    num_nonzero_basic_dual = 0;
    max_nonzero_basic_dual = 0;
    sum_nonzero_basic_dual = 0;

    num_off_bound_nonbasic = 0;
    max_off_bound_nonbasic = 0;
    sum_off_bound_nonbasic = 0;
  } else {
    num_nonzero_basic_dual = kHighsIllegalInfeasibilityCount;
    max_nonzero_basic_dual = kHighsIllegalInfeasibilityMeasure;
    sum_nonzero_basic_dual = kHighsIllegalInfeasibilityMeasure;

    num_off_bound_nonbasic = kHighsIllegalInfeasibilityCount;
    max_off_bound_nonbasic = kHighsIllegalInfeasibilityMeasure;
    sum_off_bound_nonbasic = kHighsIllegalInfeasibilityMeasure;
  }
  // Without a primal solution or a basis, nothing can be done!
  if (!have_primal_solution || !have_basis) return;
  // Makes no sense to get primal dual basis failures without a primal
  // solution, dual solution and basis
  assert(have_primal_solution && have_dual_solution && have_basis);
  double primal_infeasibility;
  double relative_primal_infeasibility;
  double dual_infeasibility;
  double value_residual;
  double lower;
  double upper;
  double value;
  double dual = 0;
  HighsBasisStatus status;
  bool status_value_ok;
  for (HighsInt iVar = 0; iVar < lp.num_col_ + lp.num_row_; iVar++) {
    if (iVar < lp.num_col_) {
      HighsInt iCol = iVar;
      lower = lp.col_lower_[iCol];
      upper = lp.col_upper_[iCol];
      value = solution.col_value[iCol];
      dual = solution.col_dual[iCol];
      status = basis.col_status[iCol];
    } else {
      HighsInt iRow = iVar - lp.num_col_;
      lower = lp.row_lower_[iRow];
      upper = lp.row_upper_[iRow];
      value = solution.row_value[iRow];
      dual = solution.row_dual[iRow];
      status = basis.row_status[iRow];
    }
    value_residual =
        std::min(std::fabs(lower - value), std::fabs(value - upper));
    // Flip dual according to lp.sense_
    dual *= (HighsInt)lp.sense_;
    status_value_ok = true;
    // Check that kLower and kUpper are consistent with value and
    // bounds - for debugging QP basis errors
    //
    // With very large values, accuracy is lost in adding/subtracting
    // the feasibility tolerance from the bounds, so skip if this may
    // occur
    if (status == HighsBasisStatus::kLower) {
      if (std::fabs(lower) / primal_feasibility_tolerance < 1e-16)
        status_value_ok = value >= lower - primal_feasibility_tolerance &&
                          value <= lower + primal_feasibility_tolerance;
    } else if (status == HighsBasisStatus::kUpper) {
      if (std::fabs(upper) / primal_feasibility_tolerance < 1e-16)
        status_value_ok = value >= upper - primal_feasibility_tolerance &&
                          value <= upper + primal_feasibility_tolerance;
    }

    if (!status_value_ok)
      highsLogUser(
          options.log_options, HighsLogType::kError,
          "getPrimalDualBasisErrors: %s %d status-value error: [%23.18g; "
          "%23.18g; %23.18g] has "
          "residual %g\n",
          iVar < lp.num_col_ ? "Column" : "Row   ",
          iVar < lp.num_col_ ? int(iVar) : int(iVar - lp.num_col_), lower,
          value, upper, value_residual);
    assert(status_value_ok);

    if (status == HighsBasisStatus::kBasic) {
      double abs_basic_dual = std::fabs(dual);
      if (abs_basic_dual > 0) {
        num_nonzero_basic_dual++;
        max_nonzero_basic_dual =
            std::max(abs_basic_dual, max_nonzero_basic_dual);
        sum_nonzero_basic_dual += abs_basic_dual;
      }
    } else {
      double off_bound_nonbasic = value_residual;
      if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
      max_off_bound_nonbasic =
          std::max(off_bound_nonbasic, max_off_bound_nonbasic);
      sum_off_bound_nonbasic += off_bound_nonbasic;
    }
  }
}

bool getComplementarityViolations(const HighsLp& lp,
                                  const HighsSolution& solution,
                                  const double optimality_tolerance,
                                  HighsInt& num_complementarity_violation,
                                  double& max_complementarity_violation) {
  num_complementarity_violation = kHighsIllegalComplementarityCount;
  max_complementarity_violation = kHighsIllegalComplementarityViolation;
  if (!solution.dual_valid) return false;

  num_complementarity_violation = 0;
  max_complementarity_violation = 0;
  double primal_residual = 0;
  for (HighsInt iVar = 0; iVar < lp.num_col_ + lp.num_row_; iVar++) {
    const bool is_col = iVar < lp.num_col_;
    const HighsInt iRow = iVar - lp.num_col_;
    const double primal =
        is_col ? solution.col_value[iVar] : solution.row_value[iRow];
    const double dual =
        is_col ? solution.col_dual[iVar] : solution.row_dual[iRow];
    const double lower = is_col ? lp.col_lower_[iVar] : lp.row_lower_[iRow];
    const double upper = is_col ? lp.col_upper_[iVar] : lp.row_upper_[iRow];
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free
      primal_residual = 1;
    } else {
      const double mid = (lower + upper) * 0.5;
      primal_residual =
          primal < mid ? std::fabs(lower - primal) : std::fabs(upper - primal);
    }
    const double dual_residual = std::fabs(dual);
    const double complementarity_violation = primal_residual * dual_residual;
    if (complementarity_violation > optimality_tolerance)
      num_complementarity_violation++;
    max_complementarity_violation =
        std::max(complementarity_violation, max_complementarity_violation);
  }
  return true;
}

bool computeDualObjectiveValue(const HighsModel& model,
                               const HighsSolution& solution,
                               double& dual_objective_value) {
  const HighsLp& lp = model.lp_;
  if (!model.isQp())
    return computeDualObjectiveValue(nullptr, lp, solution,
                                     dual_objective_value);
  assert(solution.col_value.size() == static_cast<size_t>(lp.num_col_));
  // Model is QP, so compute gradient Qx + c so generic
  // computeDualObjectiveValue can be used
  std::vector<double> gradient;
  model.objectiveGradient(solution.col_value, gradient);
  return computeDualObjectiveValue(gradient.data(), lp, solution,
                                   dual_objective_value);
}

bool computeDualObjectiveValue(const double* gradient, const HighsLp& lp,
                               const HighsSolution& solution,
                               double& dual_objective_value) {
  dual_objective_value = 0;
  if (!solution.dual_valid) return false;
  // #2184 Make sure that the solution corresponds to this LP
  assert(solution.col_value.size() == static_cast<size_t>(lp.num_col_));
  assert(solution.col_dual.size() == static_cast<size_t>(lp.num_col_));
  assert(solution.row_value.size() == static_cast<size_t>(lp.num_row_));
  assert(solution.row_dual.size() == static_cast<size_t>(lp.num_row_));

  dual_objective_value = lp.offset_;
  if (gradient) {
    // The dual objective for a QP has a -(1/2)x^TQx term, and this
    // can be computed from the gradient (g = Qx + c) as
    // -(1/2)(g-c)^Tx = (1/2)(c-g)^Tx, a pointer to the gradient data
    // is passed if this is necessary
    double quad_value = 0;
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      quad_value +=
          (lp.col_cost_[iCol] - gradient[iCol]) * solution.col_value[iCol];
    }
    dual_objective_value += 0.5 * quad_value;
  }
  double bound = 0;
  for (HighsInt iVar = 0; iVar < lp.num_col_ + lp.num_row_; iVar++) {
    const bool is_col = iVar < lp.num_col_;
    const HighsInt iRow = iVar - lp.num_col_;
    const double primal =
        is_col ? solution.col_value[iVar] : solution.row_value[iRow];
    const double dual =
        is_col ? solution.col_dual[iVar] : solution.row_dual[iRow];
    const double lower = is_col ? lp.col_lower_[iVar] : lp.row_lower_[iRow];
    const double upper = is_col ? lp.col_upper_[iVar] : lp.row_upper_[iRow];
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free
      bound = 1;
    } else {
      const double mid = (lower + upper) * 0.5;
      bound = primal < mid ? lower : upper;
    }
    dual_objective_value += bound * dual;
  }
  return true;
}

void HighsError::print(std::string message) {
  printf(
      "\n%s\nAbsolute value = %11.4g; index = %9d\nRelative value = %11.4g; "
      "index = %9d\n",
      message.c_str(), this->absolute_value, (int)this->absolute_index,
      this->relative_value, (int)this->relative_index);
}

void HighsError::reset() {
  this->absolute_value = 0;
  this->absolute_index = 0;
  this->relative_value = 0;
  this->relative_index = 0;
}

void HighsError::invalidate() {
  this->absolute_value = kHighsIllegalErrorValue;
  this->absolute_index = kHighsIllegalErrorIndex;
  this->relative_value = kHighsIllegalErrorValue;
  this->relative_index = kHighsIllegalErrorIndex;
}

double computeObjectiveValue(const HighsLp& lp, const HighsSolution& solution) {
  double objective_value = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    objective_value += lp.col_cost_[iCol] * solution.col_value[iCol];
  objective_value += lp.offset_;
  return objective_value;
}

// Refine any HighsBasisStatus::kNonbasic settings according to the LP
// and any solution values
void refineBasis(const HighsLp& lp, const HighsSolution& solution,
                 HighsBasis& basis) {
  assert(basis.useful);
  assert(isBasisRightSize(lp, basis));
  const bool have_highs_solution = solution.value_valid;

  const HighsInt num_col = lp.num_col_;
  const HighsInt num_row = lp.num_row_;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (basis.col_status[iCol] != HighsBasisStatus::kNonbasic) continue;
    const double lower = lp.col_lower_[iCol];
    const double upper = lp.col_upper_[iCol];
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
    const double lower = lp.row_lower_[iRow];
    const double upper = lp.row_upper_[iRow];
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

HighsStatus ipxSolutionToHighsSolution(
    const HighsOptions& options, const HighsLp& lp,
    const std::vector<double>& rhs, const std::vector<char>& constraint_type,
    const HighsInt ipx_num_col, const HighsInt ipx_num_row,
    const std::vector<double>& ipx_x, const std::vector<double>& ipx_slack_vars,
    const std::vector<double>& ipx_y, const std::vector<double>& ipx_zl,
    const std::vector<double>& ipx_zu, HighsSolution& highs_solution) {
  // Resize the HighsSolution
  highs_solution.col_value.resize(lp.num_col_);
  highs_solution.row_value.resize(lp.num_row_);
  highs_solution.col_dual.resize(lp.num_col_);
  highs_solution.row_dual.resize(lp.num_row_);

  const std::vector<double>& ipx_col_value = ipx_x;
  const std::vector<double>& ipx_row_value = ipx_slack_vars;

  // Row activities are needed to set activity values of free rows -
  // which are ignored by IPX
  vector<double> row_activity;
  const bool get_row_activities = ipx_num_row < lp.num_row_;
  if (get_row_activities) row_activity.assign(lp.num_row_, 0);
  HighsInt ipx_slack = lp.num_col_;
  assert(ipx_num_row == lp.num_row_);
  HighsInt dual_infeasibility_count = 0;
  double primal_infeasibility;
  double relative_primal_infeasibility;
  double dual_infeasibility;
  double value_residual;
  uint8_t at_status;   // Not used
  uint8_t mid_status;  // Not used
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    double value = ipx_col_value[iCol];
    if (get_row_activities) {
      // Accumulate row activities to assign value to free rows
      for (HighsInt el = lp.a_matrix_.start_[iCol];
           el < lp.a_matrix_.start_[iCol + 1]; el++) {
        HighsInt row = lp.a_matrix_.index_[el];
        row_activity[row] += value * lp.a_matrix_.value_[el];
      }
    }
    double dual = ipx_zl[iCol] - ipx_zu[iCol];
    highs_solution.col_value[iCol] = value;
    highs_solution.col_dual[iCol] = dual;
  }
  HighsInt ipx_row = 0;
  ipx_slack = lp.num_col_;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    double lower = lp.row_lower_[iRow];
    double upper = lp.row_upper_[iRow];
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free row - removed by IPX so set it to its row activity
      highs_solution.row_value[iRow] = row_activity[iRow];
      highs_solution.row_dual[iRow] = 0;
      continue;
    }
    // Non-free row, so IPX will have it
    double value = 0;
    double dual = 0;
    if ((lower > -kHighsInf && upper < kHighsInf) && (lower < upper)) {
      assert(constraint_type[ipx_row] == '=');
      // Boxed row - look at its slack
      value = ipx_col_value[ipx_slack];
      dual = ipx_zl[ipx_slack] - ipx_zu[ipx_slack];
      // Update the slack to be used for boxed rows
      ipx_slack++;
    } else {
      value = rhs[ipx_row] - ipx_row_value[ipx_row];
      dual = ipx_y[ipx_row];
    }
    highs_solution.row_value[iRow] = value;
    highs_solution.row_dual[iRow] = dual;
    // Update the IPX row index
    ipx_row++;
  }
  assert(ipx_row == ipx_num_row);
  assert(ipx_slack == ipx_num_col);
  if (lp.sense_ == ObjSense::kMaximize) {
    // Flip dual values since original LP is maximization
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
      highs_solution.col_dual[iCol] *= -1;
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
      highs_solution.row_dual[iRow] *= -1;
  }

  // Indicate that the primal and dual solution are known
  highs_solution.value_valid = true;
  highs_solution.dual_valid = true;
  return HighsStatus::kOk;
}

HighsStatus ipxBasicSolutionToHighsBasicSolution(
    const HighsLogOptions& log_options, const HighsLp& lp,
    const std::vector<double>& rhs, const std::vector<char>& constraint_type,
    const IpxSolution& ipx_solution, HighsBasis& highs_basis,
    HighsSolution& highs_solution) {
  // Resize the HighsSolution and HighsBasis
  highs_solution.col_value.resize(lp.num_col_);
  highs_solution.row_value.resize(lp.num_row_);
  highs_solution.col_dual.resize(lp.num_col_);
  highs_solution.row_dual.resize(lp.num_row_);
  highs_basis.col_status.resize(lp.num_col_);
  highs_basis.row_status.resize(lp.num_row_);

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
  bool get_row_activities = ipx_solution.num_row < lp.num_row_;
  if (get_row_activities) row_activity.assign(lp.num_row_, 0);
  HighsInt num_basic_variables = 0;
  for (HighsInt col = 0; col < lp.num_col_; col++) {
    bool unrecognised = false;
    if (ipx_col_status[col] == ipx_basic) {
      // Column is basic
      highs_basis.col_status[col] = HighsBasisStatus::kBasic;
      highs_solution.col_value[col] = ipx_col_value[col];
      highs_solution.col_dual[col] = 0;
    } else {
      // Column is nonbasic. Setting of ipx_col_status is consistent
      // with dual value for fixed columns
      if (ipx_col_status[col] == ipx_nonbasic_at_lb) {
        // Column is at lower bound
        highs_basis.col_status[col] = HighsBasisStatus::kLower;
        highs_solution.col_value[col] = ipx_col_value[col];
        highs_solution.col_dual[col] = ipx_col_dual[col];
      } else if (ipx_col_status[col] == ipx_nonbasic_at_ub) {
        // Column is at upper bound
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
    }
    if (unrecognised) {
      highsLogDev(log_options, HighsLogType::kError,
                  "Bounds [%11.4g, %11.4g]\n", lp.col_lower_[col],
                  lp.col_upper_[col]);
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
      for (HighsInt el = lp.a_matrix_.start_[col];
           el < lp.a_matrix_.start_[col + 1]; el++) {
        HighsInt row = lp.a_matrix_.index_[el];
        row_activity[row] +=
            highs_solution.col_value[col] * lp.a_matrix_.value_[el];
      }
    }
    if (highs_basis.col_status[col] == HighsBasisStatus::kBasic)
      num_basic_variables++;
  }
  HighsInt ipx_row = 0;
  HighsInt ipx_slack = lp.num_col_;
  HighsInt num_boxed_rows = 0;
  HighsInt num_boxed_rows_basic = 0;
  HighsInt num_boxed_row_slacks_basic = 0;
  for (HighsInt row = 0; row < lp.num_row_; row++) {
    bool unrecognised = false;
    double lower = lp.row_lower_[row];
    double upper = lp.row_upper_[row];
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
        // @FlipRowDual -slack_dual became slack_dual
        double dual = slack_dual;
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
        // @FlipRowDual -ipx_row_dual[ipx_row]; became ipx_row_dual[ipx_row];
        double dual = ipx_row_dual[ipx_row];
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
          // Row is at its fixed value: set HighsBasisStatus according
          // to sign of dual.
          //
          // Don't worry about maximization problems. IPX solves them
          // as minimizations with negated costs, so a negative dual
          // yields HighsBasisStatus::kUpper here, and dual signs are
          // then flipped below, so HighsBasisStatus::kUpper will have
          // corresponding positive dual.
          highs_basis.row_status[row] =
              dual >= 0 ? HighsBasisStatus::kLower : HighsBasisStatus::kUpper;
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
                  "Bounds [%11.4g, %11.4g]\n", lp.row_lower_[row],
                  lp.row_upper_[row]);
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
  assert(num_basic_variables == lp.num_row_);
  assert(ipx_row == ipx_solution.num_row);
  assert(ipx_slack == ipx_solution.num_col);

  if (lp.sense_ == ObjSense::kMaximize) {
    // Flip dual values since original LP is maximization
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
      highs_solution.col_dual[iCol] *= -1;
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
      highs_solution.row_dual[iRow] *= -1;
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
  highs_basis.useful = true;
  return HighsStatus::kOk;
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
  if (!passed_scaled) considerScaling(options, lp);
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

void resetModelStatusAndHighsInfo(HighsLpSolverObject& solver_object) {
  resetModelStatusAndHighsInfo(solver_object.model_status_,
                               solver_object.highs_info_);
}

void resetModelStatusAndHighsInfo(HighsModelStatus& model_status,
                                  HighsInfo& highs_info) {
  model_status = HighsModelStatus::kNotset;
  highs_info.objective_function_value = 0;
  highs_info.primal_solution_status = kSolutionStatusNone;
  highs_info.dual_solution_status = kSolutionStatusNone;
  highs_info.invalidateKkt();
}

bool isBasisConsistent(const HighsLp& lp, const HighsBasis& basis) {
  if (!isBasisRightSize(lp, basis)) return false;

  HighsInt num_basic_variables = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (basis.col_status[iCol] == HighsBasisStatus::kBasic)
      num_basic_variables++;
  }
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (basis.row_status[iRow] == HighsBasisStatus::kBasic)
      num_basic_variables++;
  }
  return num_basic_variables == lp.num_row_;
}

bool isColPrimalSolutionRightSize(const HighsLp& lp,
                                  const HighsSolution& solution) {
  return solution.col_value.size() == static_cast<size_t>(lp.num_col_);
}

bool isRowPrimalSolutionRightSize(const HighsLp& lp,
                                  const HighsSolution& solution) {
  return solution.row_value.size() == static_cast<size_t>(lp.num_row_);
}

bool isPrimalSolutionRightSize(const HighsLp& lp,
                               const HighsSolution& solution) {
  return isColPrimalSolutionRightSize(lp, solution) &&
         isRowPrimalSolutionRightSize(lp, solution);
}

bool isColDualSolutionRightSize(const HighsLp& lp,
                                const HighsSolution& solution) {
  return solution.col_dual.size() == static_cast<size_t>(lp.num_col_);
}

bool isRowDualSolutionRightSize(const HighsLp& lp,
                                const HighsSolution& solution) {
  return solution.row_dual.size() == static_cast<size_t>(lp.num_row_);
}

bool isDualSolutionRightSize(const HighsLp& lp, const HighsSolution& solution) {
  return isColDualSolutionRightSize(lp, solution) &&
         isRowDualSolutionRightSize(lp, solution);
}

bool isSolutionRightSize(const HighsLp& lp, const HighsSolution& solution) {
  return isPrimalSolutionRightSize(lp, solution) &&
         isDualSolutionRightSize(lp, solution);
}

bool isBasisRightSize(const HighsLp& lp, const HighsBasis& basis) {
  return basis.col_status.size() == static_cast<size_t>(lp.num_col_) &&
         basis.row_status.size() == static_cast<size_t>(lp.num_row_);
}

void reportLpKktFailures(const HighsLp& lp, const HighsOptions& options,
                         const HighsInfo& info, const std::string& message) {
  const HighsLogOptions& log_options = options.log_options;
  double primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  double dual_feasibility_tolerance = options.dual_feasibility_tolerance;
  double primal_residual_tolerance = options.primal_residual_tolerance;
  double dual_residual_tolerance = options.dual_residual_tolerance;
  double optimality_tolerance = options.optimality_tolerance;
  if (options.kkt_tolerance != kDefaultKktTolerance) {
    primal_feasibility_tolerance = options.kkt_tolerance;
    dual_feasibility_tolerance = options.kkt_tolerance;
    primal_residual_tolerance = options.kkt_tolerance;
    dual_residual_tolerance = options.kkt_tolerance;
    optimality_tolerance = options.kkt_tolerance;
  }

  const bool force_report = false;
  const bool has_kkt_failures =
      info.num_primal_infeasibilities > 0 ||
      info.num_dual_infeasibilities > 0 ||
      info.num_primal_residual_errors > 0 ||
      info.num_dual_residual_errors > 0 ||
      info.primal_dual_objective_error > optimality_tolerance;
  if (!has_kkt_failures && !force_report) return;

  HighsLogType log_type =
      has_kkt_failures ? HighsLogType::kWarning : HighsLogType::kInfo;

  highsLogUser(log_options, log_type, "LP solution KKT conditions%s%s\n",
               message == "" ? "" : ": ", message == "" ? "" : message.c_str());

  highsLogUser(
      log_options, HighsLogType::kInfo,
      "num/max %6d / %8.3g (relative %6d / %8.3g) primal "
      "infeasibilities     (tolerance = %4.0e)\n",
      int(info.num_primal_infeasibilities), info.max_primal_infeasibility,
      int(info.num_relative_primal_infeasibilities),
      info.max_relative_primal_infeasibility, primal_feasibility_tolerance);
  highsLogUser(log_options, HighsLogType::kInfo,
               "num/max %6d / %8.3g (relative %6d / %8.3g)   dual "
               "infeasibilities     (tolerance = %4.0e)\n",
               int(info.num_dual_infeasibilities), info.max_dual_infeasibility,
               int(info.num_relative_dual_infeasibilities),
               info.max_relative_dual_infeasibility,
               dual_feasibility_tolerance);
  if (info.num_primal_residual_errors >= 0)
    highsLogUser(
        log_options, HighsLogType::kInfo,
        "num/max %6d / %8.3g (relative %6d / %8.3g) primal residual "
        "errors     (tolerance = %4.0e)\n",
        int(info.num_primal_residual_errors), info.max_primal_residual_error,
        int(info.num_relative_primal_residual_errors),
        info.max_relative_primal_residual_error, primal_residual_tolerance);
  if (info.num_dual_residual_errors >= 0)
    highsLogUser(
        log_options, HighsLogType::kInfo,
        "num/max %6d / %8.3g (relative %6d / %8.3g)   dual residual "
        "errors     (tolerance = %4.0e)\n",
        int(info.num_dual_residual_errors), info.max_dual_residual_error,
        int(info.num_relative_dual_residual_errors),
        info.max_relative_dual_residual_error, dual_residual_tolerance);
  if (info.primal_dual_objective_error !=
      kHighsIllegalComplementarityViolation) {
    highsLogUser(
        log_options, HighsLogType::kInfo,
        "                          "
        "               %1d / %8.3g  P-D objective error        "
        "(tolerance = %4.0e)\n",
        info.primal_dual_objective_error > optimality_tolerance ? 1 : 0,
        info.primal_dual_objective_error, optimality_tolerance);
  }
  if (printf_kkt) {
    printf("grepLpKktFailures,%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%g\n",
           options.solver.c_str(), lp.model_name_.c_str(),
           lp.origin_name_.c_str(), int(info.num_primal_infeasibilities),
           int(info.num_dual_infeasibilities),
           int(info.num_primal_residual_errors),
           int(info.num_dual_residual_errors),
           int(info.num_relative_primal_infeasibilities),
           int(info.num_relative_dual_infeasibilities),
           int(info.num_relative_primal_residual_errors),
           int(info.num_relative_dual_residual_errors),
           info.primal_dual_objective_error);
  }
}

bool HighsSolution::hasUndefined() const {
  for (double value : this->col_value)
    if (value == kHighsUndefined) return true;
  return false;
}

void HighsSolution::invalidate() {
  this->value_valid = false;
  this->dual_valid = false;
}

void HighsSolution::clear() {
  this->invalidate();
  this->col_value.clear();
  this->row_value.clear();
  this->col_dual.clear();
  this->row_dual.clear();
}

void HighsObjectiveSolution::clear() { this->col_value.clear(); }

void HighsBasis::print(std::string message) const {
  if (!this->useful) return;
  this->printScalars(message);
  for (HighsInt iCol = 0; iCol < HighsInt(this->col_status.size()); iCol++)
    printf("Basis: col_status[%2d] = %d\n", int(iCol),
           int(this->col_status[iCol]));
  for (HighsInt iRow = 0; iRow < HighsInt(this->row_status.size()); iRow++)
    printf("Basis: row_status[%2d] = %d\n", int(iRow),
           int(this->row_status[iRow]));
}

void HighsBasis::printScalars(std::string message) const {
  printf("\nBasis: %s\n", message.c_str());
  printf(" valid = %d\n", this->valid);
  printf(" alien = %d\n", this->alien);
  printf(" useful = %d\n", this->useful);
  printf(" was_alien = %d\n", this->was_alien);
  printf(" debug_id = %d\n", int(this->debug_id));
  printf(" debug_update_count = %d\n", int(this->debug_update_count));
  printf(" debug_origin_name = %s\n", this->debug_origin_name.c_str());
}

void HighsBasis::invalidate() {
  this->valid = false;
  this->alien = true;
  this->useful = false;
  this->was_alien = true;
  this->debug_id = -1;
  this->debug_update_count = -1;
  this->debug_origin_name = "None";
}

void HighsBasis::clear() {
  this->invalidate();
  this->row_status.clear();
  this->col_status.clear();
}
