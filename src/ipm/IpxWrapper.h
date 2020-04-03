/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/IpxWrapper.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
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
  assert(lp.rowLower_.size() == (unsigned int)num_row);
  assert(lp.rowUpper_.size() == (unsigned int)num_row);

  std::vector<int> general_bounded_rows;
  std::vector<int> free_rows;

  for (int row = 0; row < num_row; row++)
    if (lp.rowLower_[row] < lp.rowUpper_[row] &&
        lp.rowLower_[row] > -HIGHS_CONST_INF &&
        lp.rowUpper_[row] < HIGHS_CONST_INF)
      general_bounded_rows.push_back(row);
    else if (lp.rowLower_[row] == -HIGHS_CONST_INF &&
             lp.rowUpper_[row] == HIGHS_CONST_INF)
      free_rows.push_back(row);

  const int num_slack = general_bounded_rows.size();

  // For each row except free rows add entry to char array and set up rhs
  // vector
  rhs.reserve(num_row);
  constraint_type.reserve(num_row);

  for (int row = 0; row < num_row; row++) {
    if (lp.rowLower_[row] > -HIGHS_CONST_INF &&
        lp.rowUpper_[row] == HIGHS_CONST_INF) {
      rhs.push_back(lp.rowLower_[row]);
      constraint_type.push_back('>');
    } else if (lp.rowLower_[row] == -HIGHS_CONST_INF &&
               lp.rowUpper_[row] < HIGHS_CONST_INF) {
      rhs.push_back(lp.rowUpper_[row]);
      constraint_type.push_back('<');
    } else if (lp.rowLower_[row] == lp.rowUpper_[row]) {
      rhs.push_back(lp.rowUpper_[row]);
      constraint_type.push_back('=');
    } else if (lp.rowLower_[row] > -HIGHS_CONST_INF &&
               lp.rowUpper_[row] < HIGHS_CONST_INF) {
      // general bounded
      rhs.push_back(0);
      constraint_type.push_back('=');
    }
  }

  std::vector<int> reduced_rowmap(lp.numRow_, -1);
  if (free_rows.size() > 0) {
    int counter = 0;
    int findex = 0;
    for (int row = 0; row < lp.numRow_; row++) {
      if (free_rows[findex] == row) {
        findex++;
        continue;
      } else {
        reduced_rowmap[row] = counter;
        counter++;
      }
    }
  } else {
    for (int k = 0; k < lp.numRow_; k++) reduced_rowmap[k] = k;
  }
  num_row -= free_rows.size();
  num_col += num_slack;

  std::vector<int> sizes(num_col, 0);

  for (int col = 0; col < lp.numCol_; col++)
    for (int k = lp.Astart_[col]; k < lp.Astart_[col + 1]; k++) {
      int row = lp.Aindex_[k];
      if (lp.rowLower_[row] > -HIGHS_CONST_INF ||
          lp.rowUpper_[row] < HIGHS_CONST_INF)
        sizes[col]++;
    }
  // Copy Astart and Aindex to ipx::Int array.
  int nnz = lp.Aindex_.size();
  Ap.resize(num_col + 1);
  Ai.reserve(nnz + num_slack);
  Ax.reserve(nnz + num_slack);

  // Set starting points of original and newly introduced columns.
  Ap[0] = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    Ap[col + 1] = Ap[col] + sizes[col];
    //    printf("Struc Ap[%2d] = %2d; Al[%2d] = %2d\n", col, (int)Ap[col], col,
    //    (int)sizes[col]);
  }
  for (int col = lp.numCol_; col < (int)num_col; col++) {
    Ap[col + 1] = Ap[col] + 1;
    //    printf("Slack Ap[%2d] = %2d\n", col, (int)Ap[col]);
  }
  //  printf("Fictn Ap[%2d] = %2d\n", (int)num_col, (int)Ap[num_col]);
  for (int k = 0; k < nnz; k++) {
    int row = lp.Aindex_[k];
    if (lp.rowLower_[row] > -HIGHS_CONST_INF ||
        lp.rowUpper_[row] < HIGHS_CONST_INF) {
      Ai.push_back(reduced_rowmap[lp.Aindex_[k]]);
      Ax.push_back(lp.Avalue_[k]);
    }
  }

  for (int k = 0; k < num_slack; k++) {
    Ai.push_back((ipx::Int)general_bounded_rows[k]);
    Ax.push_back(-1);
  }

  // Column bound vectors.
  col_lb.resize(num_col);
  col_ub.resize(num_col);
  for (int col = 0; col < lp.numCol_; col++) {
    if (lp.colLower_[col] == -HIGHS_CONST_INF)
      col_lb[col] = -INFINITY;
    else
      col_lb[col] = lp.colLower_[col];

    if (lp.colUpper_[col] == HIGHS_CONST_INF)
      col_ub[col] = INFINITY;
    else
      col_ub[col] = lp.colUpper_[col];
  }
  for (int slack = 0; slack < num_slack; slack++) {
    const int row = general_bounded_rows[slack];
    col_lb[lp.numCol_ + slack] = lp.rowLower_[row];
    col_ub[lp.numCol_ + slack] = lp.rowUpper_[row];
  }

  obj.resize(num_col);
  for (int col = 0; col < lp.numCol_; col++) {
    obj[col] = (int)lp.sense_ * lp.colCost_[col];
  }
  obj.insert(obj.end(), num_slack, 0);
#ifdef HiGHSDEV
  printf("IPX model has %d columns, %d rows and %d nonzeros\n", (int)num_col,
         (int)num_row, (int)Ap[num_col]);
#endif
  /*
  for (int col = 0; col < num_col; col++)
    printf("Col %2d: [%11.4g, %11.4g] Cost = %11.4g; Start = %d\n", col,
  col_lb[col], col_ub[col], obj[col], (int)Ap[col]); for (int row = 0; row <
  num_row; row++) printf("Row %2d: RHS = %11.4g; Type = %d\n", row, rhs[row],
  constraint_type[row]); for (int col = 0; col < num_col; col++) { for (int el =
  Ap[col]; el < Ap[col+1]; el++) { printf("El %2d: [%2d, %11.4g]\n", el,
  (int)Ai[el], Ax[el]);
    }
  }
  */

  return IpxStatus::OK;
}

HighsStatus solveLpIpx(const HighsLp& lp, const HighsOptions& options,
                       HighsBasis& highs_basis, HighsSolution& highs_solution,
                       HighsModelStatus& unscaled_model_status,
                       HighsSolutionParams& unscaled_solution_params) {
  int debug = 0;
  resetModelStatusAndSolutionParams(unscaled_model_status,
                                    unscaled_solution_params, options);
#ifdef CMAKE_BUILD_TYPE
  debug = 1;
#endif

  // Create the LpSolver instance
  ipx::LpSolver lps;
  // Set IPX parameters
  //
  // Cannot set internal IPX parameters directly since they are
  // private, so create instance of parameters
  ipx::Parameters parameters;
  // parameters.crossover = 1; by default
  if (debug) parameters.debug = 1;
  // Set IPX parameters from options
  // Just test feasibility and optimality tolerances for now
  // ToDo Set more parameters
  parameters.ipm_feasibility_tol =
      unscaled_solution_params.primal_feasibility_tolerance;
  parameters.ipm_optimality_tol =
      unscaled_solution_params.dual_feasibility_tolerance;
  // Set the internal IPX parameters
  lps.SetParameters(parameters);

  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;
  IpxStatus result = fillInIpxData(lp, num_col, objective, col_lb, col_ub,
                                   num_row, Ap, Ai, Av, rhs, constraint_type);
  if (result != IpxStatus::OK) return HighsStatus::Error;

  ipx::Int status =
      lps.Solve(num_col, &objective[0], &col_lb[0], &col_ub[0], num_row, &Ap[0],
                &Ai[0], &Av[0], &rhs[0], &constraint_type[0]);

#ifdef HiGHSDEV
  int int_status = status;
  if (status != 1000) printf("IPX Solve: status = %d\n", int_status);
#endif
  if (status != IPX_STATUS_solved) {
    unscaled_model_status = HighsModelStatus::SOLVE_ERROR;
    // fatal error (invalid input, out of memory, etc.)
    std::cout << " status: " << status << ','
              << " errflag: " << lps.GetInfo().errflag << '\n';
    // return IpxStatus::ErrorOrNotOptimal;
    return HighsStatus::Error;
  }

  // Get solver and solution information.
  ipx::Info ipx_info = lps.GetInfo();
  // Struct ipx_info defined in ipx/include/ipx_info.h
#ifdef HiGHSDEV
  int int_status_ipm = ipx_info.status_ipm;
  if (ipx_info.status_ipm != 1)
    printf("IPX Solve: status_ipm = %d\n", int_status_ipm);
#endif

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

#ifdef HiGHSDEV
  int int_status_crossover = ipx_info.status_crossover;
  if (int_status_crossover != 1)
    printf("IPX GetInteriorSolution: status_crossover = %d\n",
           int_status_crossover);
#endif

  if (ipx_info.status_crossover == IPX_STATUS_optimal ||
      ipx_info.status_crossover == IPX_STATUS_imprecise) {
    if (ipx_info.status_crossover == IPX_STATUS_imprecise) {
      HighsLogMessage(
          options.logfile, HighsMessageType::WARNING,
          "Ipx Crossover status imprecise: at least one of primal and dual "
          "infeasibilities of basic solution is not within parameters pfeastol "
          "and dfeastol. Simplex clean up will be required");
      // const double abs_presidual = ipx_info.abs_presidual;
      // const double abs_dresidual = ipx_info.abs_dresidual;
      // const double rel_presidual = ipx_info.rel_presidual;
      // const double rel_dresidual = ipx_info.rel_dresidual;
      // const double rel_objgap = ipx_info.rel_objgap;
    }

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

#ifdef HiGHSDEV
    int int_status_crossover = ipx_info.status_crossover;
    if (int_status_crossover != 1)
      printf("IPX GetBasicSolution: status_crossover = %d\n",
             int_status_crossover);
#endif

    // Convert the IPX basic solution to a HiGHS basic solution
    ipxToHighsBasicSolution(options.logfile, lp, rhs, constraint_type,
                            ipx_solution, highs_basis, highs_solution);

    // Set optimal
#ifdef HiGHSDEV
    if (ipx_info.status_crossover != IPX_STATUS_optimal)
      printf("IPX: Setting unscaled model status erroneously to OPTIMAL\n");
#endif
    unscaled_model_status = HighsModelStatus::OPTIMAL;
    unscaled_solution_params.ipm_iteration_count = (int)ipx_info.iter;
    unscaled_solution_params.objective_function_value = ipx_info.objval;
    getPrimalDualInfeasibilitiesFromHighsBasicSolution(
        lp, highs_basis, highs_solution, unscaled_solution_params);
  }
  //  return IpxStatus::OK;
  return HighsStatus::OK;
}
#endif
