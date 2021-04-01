/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
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
#include "io/HighsIO.h"
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

void getPrimalDualInfeasibilities(const HighsLp& lp, const HighsBasis& basis,
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

  num_primal_infeasibility = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibility = 0;
  num_dual_infeasibility = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibility = 0;

  double primal_infeasibility;
  double dual_infeasibility;
  double lower;
  double upper;
  double value;
  double dual;
  HighsBasisStatus status;
  for (HighsInt iVar = 0; iVar < lp.numCol_ + lp.numRow_; iVar++) {
    if (iVar < lp.numCol_) {
      HighsInt iCol = iVar;
      lower = lp.colLower_[iCol];
      upper = lp.colUpper_[iCol];
      value = solution.col_value[iCol];
      dual = solution.col_dual[iCol];
      status = basis.col_status[iCol];
    } else {
      HighsInt iRow = iVar - lp.numCol_;
      lower = lp.rowLower_[iRow];
      upper = lp.rowUpper_[iRow];
      value = solution.row_value[iRow];
      dual = -solution.row_dual[iRow];
      status = basis.row_status[iRow];
    }
    // Flip dual according to lp.sense_
    dual *= (HighsInt)lp.sense_;

    double primal_residual = std::max(lower - value, value - upper);
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
    //    primal_infeasibility = std::max(primal_residual, 0.);
    //    if (primal_infeasibility > primal_feasibility_tolerance)
    //      num_primal_infeasibility++;
    //    max_primal_infeasibility =
    //        std::max(primal_infeasibility, max_primal_infeasibility);
    //    sum_primal_infeasibility += primal_infeasibility;

    if (status != HighsBasisStatus::BASIC) {
      // Nonbasic variable: look for dual infeasibility
      if (primal_residual >= -primal_feasibility_tolerance) {
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
        // Between bounds (or free)
        dual_infeasibility = fabs(dual);
      }
      if (dual_infeasibility > dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
}

// Refine any HighsBasisStatus::NONBASIC settings according to the LP
// and any solution values
void refineBasis(const HighsLp& lp, const HighsSolution& solution,
                 HighsBasis& basis) {
  assert(basis.valid_);
  assert(isBasisRightSize(lp, basis));
  const bool have_highs_solution = isSolutionRightSize(lp, solution);

  const HighsInt num_col = lp.numCol_;
  const HighsInt num_row = lp.numRow_;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (basis.col_status[iCol] != HighsBasisStatus::NONBASIC) continue;
    const double lower = lp.colLower_[iCol];
    const double upper = lp.colUpper_[iCol];
    HighsBasisStatus status = HighsBasisStatus::NONBASIC;
    if (lower == upper) {
      status = HighsBasisStatus::LOWER;
    } else if (!highs_isInfinity(-lower)) {
      if (!highs_isInfinity(upper)) {
        if (have_highs_solution) {
          if (solution.col_value[iCol] < 0.5 * (lower + upper)) {
            status = HighsBasisStatus::LOWER;
          } else {
            status = HighsBasisStatus::UPPER;
          }
        } else {
          if (fabs(lower) < fabs(upper)) {
            status = HighsBasisStatus::LOWER;
          } else {
            status = HighsBasisStatus::UPPER;
          }
        }
      } else {
        status = HighsBasisStatus::LOWER;
      }
    } else if (!highs_isInfinity(upper)) {
      status = HighsBasisStatus::UPPER;
    } else {
      status = HighsBasisStatus::ZERO;
    }
    assert(status != HighsBasisStatus::NONBASIC);
    basis.col_status[iCol] = status;
  }

  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (basis.row_status[iRow] != HighsBasisStatus::NONBASIC) continue;
    const double lower = lp.rowLower_[iRow];
    const double upper = lp.rowUpper_[iRow];
    HighsBasisStatus status = HighsBasisStatus::NONBASIC;
    if (lower == upper) {
      status = HighsBasisStatus::LOWER;
    } else if (!highs_isInfinity(-lower)) {
      if (!highs_isInfinity(upper)) {
        if (have_highs_solution) {
          if (solution.row_value[iRow] < 0.5 * (lower + upper)) {
            status = HighsBasisStatus::LOWER;
          } else {
            status = HighsBasisStatus::UPPER;
          }
        } else {
          if (fabs(lower) < fabs(upper)) {
            status = HighsBasisStatus::LOWER;
          } else {
            status = HighsBasisStatus::UPPER;
          }
        }
      } else {
        status = HighsBasisStatus::LOWER;
      }
    } else if (!highs_isInfinity(upper)) {
      status = HighsBasisStatus::UPPER;
    } else {
      status = HighsBasisStatus::ZERO;
    }
    assert(status != HighsBasisStatus::NONBASIC);
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
  // No dual values are known, but ensure that the vectors are sized
  // and assigned
  highs_solution.col_dual.assign(lp.numCol_, 0.0);
  highs_solution.row_dual.assign(lp.numRow_, 0.0);

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
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    highs_solution.col_value[col] = ipx_col_value[col];
    //    highs_solution.col_dual[col] = ipx_col_dual[col];
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
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    highs_solution.col_dual[iCol] *= (HighsInt)lp.sense_;
  }
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    highs_solution.row_dual[iRow] *= (HighsInt)lp.sense_;
  }
  */
  return HighsStatus::OK;
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
#ifdef HiGHSDEV
  // For debugging, get the row activities if there are any boxed
  // constraints
  get_row_activities = get_row_activities || ipx_solution.num_col > lp.numCol_;
#endif
  if (get_row_activities) row_activity.assign(lp.numRow_, 0);
  HighsInt num_basic_variables = 0;
  for (HighsInt col = 0; col < lp.numCol_; col++) {
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
          col, (HighsInt)ipx_col_status[col]);
#endif
    }
#ifdef HiGHSDEV
    if (unrecognised)
      printf("Bounds [%11.4g, %11.4g]\n", lp.colLower_[col], lp.colUpper_[col]);
    if (unrecognised)
      printf(
          "Col %2d ipx_col_status[%2d] = %2d; x[%2d] = %11.4g; z[%2d] = "
          "%11.4g\n",
          col, col, (HighsInt)ipx_col_status[col], col, ipx_col_value[col], col,
          ipx_col_dual[col]);
#endif
    assert(!unrecognised);
    if (unrecognised) {
      highsLogUser(log_options, HighsLogType::ERROR,
                   "Unrecognised ipx_col_status value from IPX\n");
      return HighsStatus::Error;
    }
    if (get_row_activities) {
      // Accumulate row activities to assign value to free rows
      for (HighsInt el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
        HighsInt row = lp.Aindex_[el];
        row_activity[row] += highs_solution.col_value[col] * lp.Avalue_[el];
      }
    }
    if (highs_basis.col_status[col] == HighsBasisStatus::BASIC)
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
#ifdef HiGHSDEV
    HighsInt this_ipx_row = ipx_row;
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
              row, ipx_row, ipx_slack, (HighsInt)ipx_col_status[ipx_slack]);
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
          row, this_ipx_row, (HighsInt)ipx_row_status[this_ipx_row],
          this_ipx_row, ipx_row_value[this_ipx_row], this_ipx_row,
          ipx_row_dual[this_ipx_row]);
#endif
    assert(!unrecognised);
    if (unrecognised) {
      highsLogUser(log_options, HighsLogType::ERROR,
                   "Unrecognised ipx_row_status value from IPX\n");
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
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    highs_solution.col_dual[iCol] *= (HighsInt)lp.sense_;
  }
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    highs_solution.row_dual[iRow] *= (HighsInt)lp.sense_;
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
// instance.
void invalidateSolutionInfeasibilityParams(
    HighsSolutionParams& solution_params) {
  solution_params.num_primal_infeasibility = illegal_infeasibility_count;
  solution_params.max_primal_infeasibility = illegal_infeasibility_measure;
  solution_params.sum_primal_infeasibility = illegal_infeasibility_measure;
  solution_params.num_dual_infeasibility = illegal_infeasibility_count;
  solution_params.max_dual_infeasibility = illegal_infeasibility_measure;
  solution_params.sum_dual_infeasibility = illegal_infeasibility_measure;
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
    if (basis.col_status[iCol] == HighsBasisStatus::BASIC)
      num_basic_variables++;
  }
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    if (basis.row_status[iRow] == HighsBasisStatus::BASIC)
      num_basic_variables++;
  }
  bool right_num_basic_variables = num_basic_variables == lp.numRow_;
  consistent = right_num_basic_variables && consistent;
  return consistent;
}

bool isSolutionRightSize(const HighsLp& lp, const HighsSolution& solution) {
  bool right_size = true;
  right_size = (HighsInt)solution.col_value.size() == lp.numCol_ && right_size;
  right_size = (HighsInt)solution.col_dual.size() == lp.numCol_ && right_size;
  right_size = (HighsInt)solution.row_value.size() == lp.numRow_ && right_size;
  right_size = (HighsInt)solution.row_dual.size() == lp.numRow_ && right_size;
  return right_size;
}

bool isBasisRightSize(const HighsLp& lp, const HighsBasis& basis) {
  bool right_size = true;
  right_size = (HighsInt)basis.col_status.size() == lp.numCol_ && right_size;
  right_size = (HighsInt)basis.row_status.size() == lp.numRow_ && right_size;
  return right_size;
}

void clearSolutionUtil(HighsSolution& solution) {
  solution.col_dual.clear();
  solution.col_value.clear();
  solution.row_dual.clear();
  solution.row_value.clear();
}

void clearBasisUtil(HighsBasis& basis) {
  basis.row_status.clear();
  basis.col_status.clear();
  basis.valid_ = false;
}
