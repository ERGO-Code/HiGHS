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
/**@file ipm/IpxWrapper.h
 * @brief
 */
#ifndef IPM_IPX_WRAPPER_H_
#define IPM_IPX_WRAPPER_H_

#include <algorithm>

#include "ipm/IpxSolution.h"
#include "ipm/IpxStatus.h"
#include "ipm/ipx/include/ipx_status.h"
#include "ipm/ipx/src/lp_solver.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"

IpxStatus fillInIpxData(const HighsLp& lp, ipx::Int& num_col,
                        std::vector<double>& obj, std::vector<double>& col_lb,
                        std::vector<double>& col_ub, ipx::Int& num_row,
                        std::vector<ipx::Int>& Ap, std::vector<ipx::Int>& Ai,
                        std::vector<double>& Ax, std::vector<double>& rhs,
                        std::vector<char>& constraint_type) {
  num_col = lp.numCol_;
  num_row = lp.numRow_;

  // For each row with both a lower and an upper bound introduce one new column
  // so num_col may increase. Ignore each free row so num_row may decrease.
  // lba <= a'x <= uba becomes
  // a'x-s = 0 and lba <= s <= uba.

  // For each row with bounds on both sides introduce explicit slack and
  // transfer bounds.
  assert((HighsInt)lp.rowLower_.size() == num_row);
  assert((HighsInt)lp.rowUpper_.size() == num_row);

  std::vector<HighsInt> general_bounded_rows;
  std::vector<HighsInt> free_rows;

  for (HighsInt row = 0; row < num_row; row++)
    if (lp.rowLower_[row] < lp.rowUpper_[row] &&
        lp.rowLower_[row] > -kHighsInf && lp.rowUpper_[row] < kHighsInf)
      general_bounded_rows.push_back(row);
    else if (lp.rowLower_[row] <= -kHighsInf && lp.rowUpper_[row] >= kHighsInf)
      free_rows.push_back(row);

  const HighsInt num_slack = general_bounded_rows.size();

  // For each row except free rows add entry to char array and set up rhs
  // vector
  rhs.reserve(num_row);
  constraint_type.reserve(num_row);

  for (int row = 0; row < num_row; row++) {
    if (lp.rowLower_[row] > -kHighsInf && lp.rowUpper_[row] >= kHighsInf) {
      rhs.push_back(lp.rowLower_[row]);
      constraint_type.push_back('>');
    } else if (lp.rowLower_[row] <= -kHighsInf &&
               lp.rowUpper_[row] < kHighsInf) {
      rhs.push_back(lp.rowUpper_[row]);
      constraint_type.push_back('<');
    } else if (lp.rowLower_[row] == lp.rowUpper_[row]) {
      rhs.push_back(lp.rowUpper_[row]);
      constraint_type.push_back('=');
    } else if (lp.rowLower_[row] > -kHighsInf &&
               lp.rowUpper_[row] < kHighsInf) {
      // general bounded
      rhs.push_back(0);
      constraint_type.push_back('=');
    }
  }

  std::vector<HighsInt> reduced_rowmap(lp.numRow_, -1);
  if (free_rows.size() > 0) {
    HighsInt counter = 0;
    HighsInt findex = 0;
    for (HighsInt row = 0; row < lp.numRow_; row++) {
      if (free_rows[findex] == row) {
        findex++;
        continue;
      } else {
        reduced_rowmap[row] = counter;
        counter++;
      }
    }
  } else {
    for (HighsInt k = 0; k < lp.numRow_; k++) reduced_rowmap[k] = k;
  }
  num_row -= free_rows.size();
  num_col += num_slack;

  std::vector<HighsInt> sizes(num_col, 0);

  for (HighsInt col = 0; col < lp.numCol_; col++)
    for (HighsInt k = lp.Astart_[col]; k < lp.Astart_[col + 1]; k++) {
      HighsInt row = lp.Aindex_[k];
      if (lp.rowLower_[row] > -kHighsInf || lp.rowUpper_[row] < kHighsInf)
        sizes[col]++;
    }
  // Copy Astart and Aindex to ipx::Int array.
  HighsInt nnz = lp.Aindex_.size();
  Ap.resize(num_col + 1);
  Ai.reserve(nnz + num_slack);
  Ax.reserve(nnz + num_slack);

  // Set starting points of original and newly introduced columns.
  Ap[0] = 0;
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    Ap[col + 1] = Ap[col] + sizes[col];
    //    printf("Struc Ap[%2" HIGHSINT_FORMAT "] = %2" HIGHSINT_FORMAT ";
    //    Al[%2" HIGHSINT_FORMAT "] = %2" HIGHSINT_FORMAT "\n", col,
    //    (int)Ap[col], col, (int)sizes[col]);
  }
  for (HighsInt col = lp.numCol_; col < (HighsInt)num_col; col++) {
    Ap[col + 1] = Ap[col] + 1;
    //    printf("Slack Ap[%2" HIGHSINT_FORMAT "] = %2" HIGHSINT_FORMAT "\n",
    //    col, (int)Ap[col]);
  }
  //  printf("Fictn Ap[%2" HIGHSINT_FORMAT "] = %2" HIGHSINT_FORMAT "\n",
  //  (int)num_col, (int)Ap[num_col]);
  for (HighsInt k = 0; k < nnz; k++) {
    HighsInt row = lp.Aindex_[k];
    if (lp.rowLower_[row] > -kHighsInf || lp.rowUpper_[row] < kHighsInf) {
      Ai.push_back(reduced_rowmap[lp.Aindex_[k]]);
      Ax.push_back(lp.Avalue_[k]);
    }
  }

  for (HighsInt k = 0; k < num_slack; k++) {
    Ai.push_back((ipx::Int)general_bounded_rows[k]);
    Ax.push_back(-1);
  }

  // Column bound vectors.
  col_lb.resize(num_col);
  col_ub.resize(num_col);
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    if (lp.colLower_[col] <= -kHighsInf)
      col_lb[col] = -INFINITY;
    else
      col_lb[col] = lp.colLower_[col];

    if (lp.colUpper_[col] >= kHighsInf)
      col_ub[col] = INFINITY;
    else
      col_ub[col] = lp.colUpper_[col];
  }
  for (HighsInt slack = 0; slack < num_slack; slack++) {
    const int row = general_bounded_rows[slack];
    col_lb[lp.numCol_ + slack] = lp.rowLower_[row];
    col_ub[lp.numCol_ + slack] = lp.rowUpper_[row];
  }

  obj.resize(num_col);
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    obj[col] = (HighsInt)lp.sense_ * lp.colCost_[col];
  }
  obj.insert(obj.end(), num_slack, 0);
  /*
  for (int col = 0; col < num_col; col++)
    printf("Col %2" HIGHSINT_FORMAT ": [%11.4g, %11.4g] Cost = %11.4g; Start =
  %" HIGHSINT_FORMAT "\n", col, col_lb[col], col_ub[col], obj[col],
  (int)Ap[col]); for (int row = 0; row < num_row; row++) printf("Row %2"
  HIGHSINT_FORMAT ": RHS = %11.4g; Type = %" HIGHSINT_FORMAT "\n", row,
  rhs[row], constraint_type[row]); for (int col = 0; col < num_col; col++) { for
  (int el = Ap[col]; el < Ap[col+1]; el++) { printf("El %2" HIGHSINT_FORMAT ":
  [%2" HIGHSINT_FORMAT ", %11.4g]\n", el, (int)Ai[el], Ax[el]);
    }
  }
  */

  return IpxStatus::OK;
}

HighsStatus reportIpxSolveStatus(const HighsOptions& options,
                                 const ipx::Int solve_status,
                                 const ipx::Int error_flag) {
  if (solve_status == IPX_STATUS_solved) {
    highsLogUser(options.log_options, HighsLogType::kInfo, "Ipx: Solved\n");
    return HighsStatus::kOk;
  } else if (solve_status == IPX_STATUS_stopped) {
    highsLogUser(options.log_options, HighsLogType::kWarning, "Ipx: Stopped\n");
    return HighsStatus::kWarning;
  } else if (solve_status == IPX_STATUS_no_model) {
    if (error_flag == IPX_ERROR_argument_null) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - argument_null\n");
      return HighsStatus::kError;
    } else if (error_flag == IPX_ERROR_invalid_dimension) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - invalid dimension\n");
      return HighsStatus::kError;
    } else if (error_flag == IPX_ERROR_invalid_matrix) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - invalid matrix\n");
      return HighsStatus::kError;
    } else if (error_flag == IPX_ERROR_invalid_vector) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - invalid vector\n");
      return HighsStatus::kError;
    } else if (error_flag == IPX_ERROR_invalid_basis) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - invalid basis\n");
      return HighsStatus::kError;
    } else {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - unrecognised error\n");
      return HighsStatus::kError;
    }
  } else if (solve_status == IPX_STATUS_out_of_memory) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Ipx: Out of memory\n");
    return HighsStatus::kError;
  } else if (solve_status == IPX_STATUS_internal_error) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Ipx: Internal error %" HIGHSINT_FORMAT "\n", (int)error_flag);
    return HighsStatus::kError;
  } else {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Ipx: unrecognised solve status = %" HIGHSINT_FORMAT "\n",
                 (int)solve_status);
    return HighsStatus::kError;
  }
  return HighsStatus::kError;
}

HighsStatus reportIpxIpmCrossoverStatus(const HighsOptions& options,
                                        const ipx::Int status,
                                        const bool ipm_status) {
  std::string method_name;
  if (ipm_status)
    method_name = "IPM      ";
  else
    method_name = "Crossover";
  if (status == IPX_STATUS_not_run) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s not run\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_optimal) {
    highsLogUser(options.log_options, HighsLogType::kInfo, "Ipx: %s optimal\n",
                 method_name.c_str());
    return HighsStatus::kOk;
  } else if (status == IPX_STATUS_imprecise) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s imprecise\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_primal_infeas) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s primal infeasible\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_dual_infeas) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s dual infeasible\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_time_limit) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s reached time limit\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_iter_limit) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s reached iteration limit\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_no_progress) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s no progress\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_failed) {
    highsLogUser(options.log_options, HighsLogType::kError, "Ipx: %s failed\n",
                 method_name.c_str());
    return HighsStatus::kError;
  } else if (status == IPX_STATUS_debug) {
    highsLogUser(options.log_options, HighsLogType::kError, "Ipx: %s debug\n",
                 method_name.c_str());
    return HighsStatus::kError;
  } else {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Ipx: %s unrecognised status\n", method_name.c_str());
    return HighsStatus::kError;
  }
  return HighsStatus::kError;
}

bool ipxStatusError(const bool status_error, const HighsOptions& options,
                    std::string message, const int value = -1) {
  if (status_error) {
    if (value < 0) {
      highsLogUser(options.log_options, HighsLogType::kError, "Ipx: %s\n",
                   message.c_str());
    } else {
      highsLogUser(options.log_options, HighsLogType::kError, "Ipx: %s %d\n",
                   message.c_str(), value);
    }
    fflush(NULL);
  }
  assert(!status_error);
  return status_error;
}

bool illegalIpxSolvedStatus(ipx::Info& ipx_info, const HighsOptions& options) {
  bool found_illegal_status = false;
  //========
  // For IPX
  //========
  // Can solve and be optimal
  // Can solve and be imprecise
  // Can solve and be primal_infeas
  // Can solve and be dual_infeas
  // Cannot solve and reach time limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_time_limit, options,
                     "solved  status_ipm should not be IPX_STATUS_time_limit");
  // Cannot solve and reach iteration limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_iter_limit, options,
                     "solved  status_ipm should not be IPX_STATUS_iter_limit");
  // Cannot solve and make no progress
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_no_progress, options,
                     "solved  status_ipm should not be IPX_STATUS_no_progress");
  // Cannot solve and failed
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_failed, options,
                     "solved  status_ipm should not be IPX_STATUS_failed");
  // Cannot solve and debug
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_debug, options,
                     "solved  status_ipm should not be IPX_STATUS_debug");
  //==============
  // For crossover
  //==============
  // Can solve and be optimal
  // Can solve and be imprecise
  // Cannot solve with primal infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_primal_infeas, options,
          "solved  status_crossover should not be IPX_STATUS_primal_infeas");
  // Cannot solve with dual infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_dual_infeas, options,
          "solved  status_crossover should not be IPX_STATUS_dual_infeas");
  // Cannot solve and reach time limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_time_limit, options,
          "solved  status_crossover should not be IPX_STATUS_time_limit");
  // Cannot solve and reach time limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_iter_limit, options,
          "solved  status_crossover should not be IPX_STATUS_iter_limit");
  // Cannot solve and make no progress
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_no_progress, options,
          "solved  status_crossover should not be IPX_STATUS_no_progress");
  // Cannot solve and failed
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_failed, options,
          "solved  status_crossover should not be IPX_STATUS_failed");
  // Cannot solve and debug
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_crossover == IPX_STATUS_debug, options,
                     "solved  status_crossover should not be IPX_STATUS_debug");
  return found_illegal_status;
}

bool illegalIpxStoppedIpmStatus(ipx::Info& ipx_info,
                                const HighsOptions& options) {
  bool found_illegal_status = false;
  // Cannot stop and be optimal
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_optimal, options,
                     "stopped status_ipm should not be IPX_STATUS_optimal");
  // Cannot stop and be imprecise
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_imprecise, options,
                     "stopped status_ipm should not be IPX_STATUS_imprecise");
  // Cannot stop with primal infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_ipm == IPX_STATUS_primal_infeas, options,
          "stopped status_ipm should not be IPX_STATUS_primal_infeas");
  // Cannot stop with dual infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_dual_infeas, options,
                     "stopped status_ipm should not be IPX_STATUS_dual_infeas");
  // Can stop with time limit
  // Can stop with iter limit
  // Can stop with no progress
  // Cannot stop and failed - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_failed, options,
                     "stopped status_ipm should not be IPX_STATUS_failed");
  // Cannot stop and debug - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_debug, options,
                     "stopped status_ipm should not be IPX_STATUS_debug");
  return found_illegal_status;
}

bool illegalIpxStoppedCrossoverStatus(ipx::Info& ipx_info,
                                      const HighsOptions& options) {
  bool found_illegal_status = false;
  // Cannot stop and be optimal
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_optimal, options,
          "stopped status_crossover should not be IPX_STATUS_optimal");
  // Cannot stop and be imprecise
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_imprecise, options,
          "stopped status_crossover should not be IPX_STATUS_imprecise");
  // Cannot stop with primal infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_primal_infeas, options,
          "stopped status_crossover should not be IPX_STATUS_primal_infeas");
  // Cannot stop with dual infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_dual_infeas, options,
          "stopped status_crossover should not be IPX_STATUS_dual_infeas");
  // Cannot stop and reach iteration limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_iter_limit, options,
          "stopped status_crossover should not be IPX_STATUS_iter_limit");
  // Can stop and reach time limit
  // Cannot stop with no_progress
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_no_progress, options,
          "stopped status_crossover should not be IPX_STATUS_no_progress");
  // Cannot stop and failed - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_failed, options,
          "stopped status_crossover should not be IPX_STATUS_failed");
  // Cannot stop and debug - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_crossover == IPX_STATUS_debug, options,
                     "stopped status_crossover should not be IPX_STATUS_debug");
  return found_illegal_status;
}

void reportIpmNoProgress(const HighsOptions& options,
                         const ipx::Info& ipx_info) {
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "No progress: primal objective value       = %11.4g\n",
               ipx_info.pobjval);
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "No progress: max absolute primal residual = %11.4g\n",
               ipx_info.abs_presidual);
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "No progress: max absolute   dual residual = %11.4g\n",
               ipx_info.abs_dresidual);
}

HighsStatus analyseIpmNoProgress(const ipx::Info& ipx_info,
                                 const ipx::Parameters& parameters,
                                 HighsModelStatus& model_status) {
  if (ipx_info.abs_presidual > parameters.ipm_feasibility_tol) {
    // Looks like the LP is infeasible
    model_status = HighsModelStatus::kPrimalInfeasible;
    return HighsStatus::kOk;
  } else if (ipx_info.abs_dresidual > parameters.ipm_optimality_tol) {
    // Looks like the LP is unbounded
    model_status = HighsModelStatus::kPrimalUnbounded;
    return HighsStatus::kOk;
  } else if (ipx_info.pobjval < -kHighsInf) {
    // Looks like the LP is unbounded
    model_status = HighsModelStatus::kPrimalUnbounded;
    return HighsStatus::kOk;
  } else {
    // Don't know
    model_status = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }
  return HighsStatus::kWarning;
}

HighsStatus solveLpIpx(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, bool& imprecise_solution,
                       HighsBasis& highs_basis, HighsSolution& highs_solution,
                       HighsIterationCounts& iteration_counts,
                       HighsModelStatus& model_status,
                       HighsSolutionParams& solution_params) {
  imprecise_solution = false;
  resetModelStatusAndSolutionParams(model_status, solution_params, options);
  // Create the LpSolver instance
  ipx::LpSolver lps;
  // Set IPX parameters
  //
  // Cannot set internal IPX parameters directly since they are
  // private, so create instance of parameters
  ipx::Parameters parameters;
  // Set IPX parameters from options
  //
  // Set display according to output
  parameters.display = 1;
  if (!options.output_flag) parameters.display = 0;
  // Modify parameters.debug according to log_dev_level
  parameters.debug = 0;
  if (options.log_dev_level == kHighsLogDevLevelDetailed) {
    // Default options.log_dev_level setting is kHighsLogDevLevelNone, yielding
    // default setting debug = 0
    parameters.debug = 0;
  } else if (options.log_dev_level == kHighsLogDevLevelInfo) {
    parameters.debug = 3;
  } else if (options.log_dev_level == kHighsLogDevLevelVerbose) {
    parameters.debug = 4;
  }
  // Just test feasibility and optimality tolerances for now
  // ToDo Set more parameters
  parameters.ipm_feasibility_tol =
      min(solution_params.primal_feasibility_tolerance,
          solution_params.dual_feasibility_tolerance);

  parameters.ipm_optimality_tol = options.ipm_optimality_tolerance;
  parameters.crossover_start = options.start_crossover_tolerance;
  // Determine the run time allowed for IPX
  parameters.time_limit = options.time_limit - timer.readRunHighsClock();
  parameters.ipm_maxiter = options.ipm_iteration_limit - iteration_counts.ipm;
  // Determine if crossover is to be run or not
  parameters.crossover = options.run_crossover;
  if (!parameters.crossover) {
    // If crossover is not run, then set crossover_start to -1 so that
    // IPX can terminate according to its feasibility and optimality
    // tolerances
    parameters.crossover_start = -1;
  }

  // Set the internal IPX parameters
  lps.SetParameters(parameters);

  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;
  IpxStatus result = fillInIpxData(lp, num_col, objective, col_lb, col_ub,
                                   num_row, Ap, Ai, Av, rhs, constraint_type);
  if (result != IpxStatus::OK) return HighsStatus::kError;
  highsLogUser(options.log_options, HighsLogType::kInfo,
               "IPX model has %" HIGHSINT_FORMAT " rows, %" HIGHSINT_FORMAT
               " columns and %" HIGHSINT_FORMAT " nonzeros\n",
               num_row, num_col, Ap[num_col]);

  ipx::Int load_status =
      lps.LoadModel(num_col, &objective[0], &col_lb[0], &col_ub[0], num_row,
                    &Ap[0], &Ai[0], &Av[0], &rhs[0], &constraint_type[0]);

  // todo: handle load error
  ipx::Int solve_status = lps.Solve();

  // Get solver and solution information.
  // Struct ipx_info defined in ipx/include/ipx_info.h
  ipx::Info ipx_info = lps.GetInfo();
  iteration_counts.ipm += (HighsInt)ipx_info.iter;

  iteration_counts.crossover += (HighsInt)ipx_info.updates_crossover;
  // iteration_counts.crossover += (int)ipx_info.pushes_crossover;

  // If not solved...
  if (solve_status != IPX_STATUS_solved) {
    const HighsStatus solve_return_status =
        reportIpxSolveStatus(options, solve_status, ipx_info.errflag);
    // Return error if IPX solve error has occurred
    if (solve_return_status == HighsStatus::kError) {
      model_status = HighsModelStatus::kSolveError;
      return HighsStatus::kError;
    }
  }
  bool ipm_status = true;
  const HighsStatus ipm_return_status =
      reportIpxIpmCrossoverStatus(options, ipx_info.status_ipm, ipm_status);
  ipm_status = false;
  const HighsStatus crossover_return_status = reportIpxIpmCrossoverStatus(
      options, ipx_info.status_crossover, ipm_status);
  // Return error if IPX IPM or crossover error has occurred
  if (ipm_return_status == HighsStatus::kError ||
      crossover_return_status == HighsStatus::kError) {
    model_status = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }
  // Should only reach here if Solve() returned IPX_STATUS_solved or
  // IPX_STATUS_stopped
  if (ipxStatusError(
          solve_status != IPX_STATUS_solved &&
              solve_status != IPX_STATUS_stopped,
          options, "solve_status should be solved or stopped here but value is",
          (int)solve_status))
    return HighsStatus::kError;

  if (solve_status == IPX_STATUS_stopped) {
    //
    // Look at the reason why IPX stopped
    //
    // Return error if stopped status settings occur that JAJH doesn't
    // think should happen
    //
    //==============
    // For crossover
    //==============
    if (illegalIpxStoppedCrossoverStatus(ipx_info, options))
      return HighsStatus::kError;
    // Can stop and reach time limit
    if (ipx_info.status_crossover == IPX_STATUS_time_limit) {
      model_status = HighsModelStatus::kReachedTimeLimit;
      return HighsStatus::kWarning;
    }
    //========
    // For IPM
    //========
    //
    // Note that IPX can stop with IPM optimal, imprecise,
    // primal_infeas or dual_infeas, due to crossover stopping with
    // time limit, and this is why crossover returns are tested first
    if (illegalIpxStoppedIpmStatus(ipx_info, options))
      return HighsStatus::kError;
    // Can stop with time limit
    // Can stop with iter limit
    // Can stop with no progress
    if (ipx_info.status_ipm == IPX_STATUS_time_limit) {
      model_status = HighsModelStatus::kReachedTimeLimit;
      return HighsStatus::kWarning;
    } else if (ipx_info.status_ipm == IPX_STATUS_iter_limit) {
      model_status = HighsModelStatus::kReachedIterationLimit;
      return HighsStatus::kWarning;
    } else if (ipx_info.status_ipm == IPX_STATUS_no_progress) {
      reportIpmNoProgress(options, ipx_info);
      return analyseIpmNoProgress(ipx_info, lps.GetParameters(), model_status);
    }
  }
  // Should only reach here if Solve() returned IPX_STATUS_solved
  if (ipxStatusError(solve_status != IPX_STATUS_solved, options,
                     "solve_status should be solved here but value is",
                     (int)solve_status))
    return HighsStatus::kError;
  // Return error if solved status settings occur that JAJH doesn't
  // think should happen
  if (illegalIpxSolvedStatus(ipx_info, options)) return HighsStatus::kError;
  //==============
  // For crossover
  //==============
  // Can be not run
  // Can solve and be optimal
  // Can solve and be imprecise
  //========
  // For IPX
  //========
  // Can solve and be optimal
  // Can solve and be imprecise
  // Can solve and be primal_infeas
  // Can solve and be dual_infeas
  if (ipx_info.status_ipm == IPX_STATUS_primal_infeas) {
    model_status = HighsModelStatus::kPrimalInfeasible;
    return HighsStatus::kOk;
  } else if (ipx_info.status_ipm == IPX_STATUS_dual_infeas) {
    model_status = HighsModelStatus::kPrimalUnbounded;
    return HighsStatus::kOk;
  }

  // Should only reach here if crossover is not run, optimal or imprecise
  if (ipxStatusError(ipx_info.status_crossover != IPX_STATUS_not_run &&
                         ipx_info.status_crossover != IPX_STATUS_optimal &&
                         ipx_info.status_crossover != IPX_STATUS_imprecise,
                     options,
                     "crossover status should be not run, optimal or imprecise "
                     "but value is",
                     (int)ipx_info.status_crossover))
    return HighsStatus::kError;

  // Get the interior solution (available if IPM was started).
  // GetInteriorSolution() returns the final IPM iterate, regardless if the
  // IPM terminated successfully or not. (Only in case of out-of-memory no
  // solution exists.)
  std::vector<double> x(num_col);
  std::vector<double> xl(num_col);
  std::vector<double> xu(num_col);
  std::vector<double> zl(num_col);
  std::vector<double> zu(num_col);
  std::vector<double> slack(num_row);
  std::vector<double> y(num_row);

  lps.GetInteriorSolution(&x[0], &xl[0], &xu[0], &slack[0], &y[0], &zl[0],
                          &zu[0]);

  // Basic solution depends on crossover being run
  const bool have_basic_solution =
      ipx_info.status_crossover != IPX_STATUS_not_run;
  imprecise_solution = ipx_info.status_crossover == IPX_STATUS_imprecise;
  if (have_basic_solution) {
    IpxSolution ipx_solution;
    ipx_solution.num_col = num_col;
    ipx_solution.num_row = num_row;
    ipx_solution.ipx_col_value.resize(num_col);
    ipx_solution.ipx_row_value.resize(num_row);
    ipx_solution.ipx_col_dual.resize(num_col);
    ipx_solution.ipx_row_dual.resize(num_row);
    ipx_solution.ipx_row_status.resize(num_row);
    ipx_solution.ipx_col_status.resize(num_col);
    lps.GetBasicSolution(
        &ipx_solution.ipx_col_value[0], &ipx_solution.ipx_row_value[0],
        &ipx_solution.ipx_row_dual[0], &ipx_solution.ipx_col_dual[0],
        &ipx_solution.ipx_row_status[0], &ipx_solution.ipx_col_status[0]);

    // Convert the IPX basic solution to a HiGHS basic solution
    ipxBasicSolutionToHighsBasicSolution(options.log_options, lp, rhs,
                                         constraint_type, ipx_solution,
                                         highs_basis, highs_solution);
  } else {
    ipxSolutionToHighsSolution(options.log_options, lp, rhs, constraint_type,
                               num_col, num_row, x, slack, highs_solution);
    highs_basis.valid_ = false;
  }
  HighsStatus return_status;
  if (imprecise_solution) {
    model_status = HighsModelStatus::kNotset;
    return_status = HighsStatus::kWarning;
  } else {
    model_status = HighsModelStatus::kOptimal;
    solution_params.primal_status = kHighsPrimalDualStatusFeasiblePoint;
    // Currently only have a dual solution if there is a basic solution
    if (have_basic_solution)
      solution_params.dual_status = kHighsPrimalDualStatusFeasiblePoint;
    return_status = HighsStatus::kOk;
  }
  double objective_function_value = lp.offset_;
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++)
    objective_function_value +=
        highs_solution.col_value[iCol] * lp.colCost_[iCol];
  solution_params.objective_function_value = objective_function_value;
  if (highs_basis.valid_)
    getPrimalDualInfeasibilities(lp, highs_basis, highs_solution,
                                 solution_params);
  return return_status;
}
#endif
