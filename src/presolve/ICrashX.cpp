/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/ICrashX.h"

#include <algorithm>
#include <iostream>

#include "HConfig.h"
#include "ipm/IpxWrapper.h"

#ifndef IPX_ON
bool callCrossover(const HighsLp& lp, const HighsOptions& options,
                   HighsSolution& solution, HighsBasis& highs_basis) {
  return false;
}
#else

#include "ipm/IpxWrapper.h"

bool callCrossover(const HighsLp& lp, const HighsOptions& options,
                   HighsSolution& solution, HighsBasis& basis) {
  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;

  fillInIpxData(lp, num_col, num_row, objective, col_lb, col_ub, Ap, Ai, Av,
                rhs, constraint_type);
  // if (res != IpxStatus::OK) return false;

  const HighsLogOptions& log_options = options.log_options;

  ipx::Parameters parameters;
  parameters.crossover = true;
  parameters.crash_basis = 1;  // 0 = slack basis; 1 = crash basis
  parameters.display = 1;
  if (!options.output_flag) parameters.display = 0;
  // Modify parameters.debug according to log_dev_level
  parameters.debug = 0;
  if (options.log_dev_level == kHighsLogDevLevelDetailed) {
    parameters.debug = 0;
  } else if (options.log_dev_level == kHighsLogDevLevelInfo) {
    parameters.debug = 2;
  } else if (options.log_dev_level == kHighsLogDevLevelVerbose) {
    parameters.debug = 4;
  }

  ipx::LpSolver lps;
  lps.SetParameters(parameters);

  ipx::Int errflag =
      lps.LoadModel(num_col, &objective[0], &col_lb[0], &col_ub[0], num_row,
                    &Ap[0], &Ai[0], &Av[0], &rhs[0], &constraint_type[0]);
  if (errflag != 0) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Error loading ipx model\n");
    return false;
  }

  // Set x values within bounds.
  std::vector<double> x(solution.col_value);
  for (int i = 0; i < num_col; i++) {
    x[i] = std::max(x[i], col_lb[i]);
    x[i] = std::min(x[i], col_ub[i]);
  }

  // Build slack variables from rhs-A*x but subject to sign conditions.
  std::vector<double> slack(rhs);
  for (int i = 0; i < num_col; i++) {
    for (int p = Ap[i]; p < Ap[i + 1]; ++p) slack[Ai[p]] -= Av[p] * x[i];
  }
  for (int i = 0; i < num_row; i++) {
    switch (constraint_type[i]) {
      case '=':
        slack[i] = 0.0;
        break;
      case '<':
        slack[i] = std::max(slack[i], 0.0);
        break;
      case '>':
        slack[i] = std::min(slack[i], 0.0);
        break;
    }
  }
  if (solution.dual_valid && solution.col_dual.size() == num_col &&
      solution.row_dual.size() == num_row) {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Calling IPX crossover with primal and dual values\n");
    errflag = lps.CrossoverFromStartingPoint(
        &x[0], &slack[0], &solution.row_dual[0], &solution.col_dual[0]);
  } else {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Calling IPX crossover with only primal values\n");
    errflag =
        lps.CrossoverFromStartingPoint(&x[0], &slack[0], nullptr, nullptr);
  }

  if (errflag != 0) {
    highsLogUser(log_options, HighsLogType::kError,
                 "IPX crossover error: flag = %d\n", (int)errflag);
    return false;
  }
  ipx::Info info = lps.GetInfo();
  if (info.status_crossover != IPX_STATUS_optimal &&
      info.status_crossover != IPX_STATUS_imprecise) {
    highsLogUser(log_options, HighsLogType::kError,
                 "IPX crossover failed: status = %d\n",
                 (int)info.status_crossover);
    return false;
  }

  // Get basis
  IpxSolution ipx_solution;
  ipx_solution.num_col = num_col;
  ipx_solution.num_row = num_row;
  ipx_solution.ipx_col_value.resize(num_col);
  ipx_solution.ipx_row_value.resize(num_row);
  ipx_solution.ipx_col_dual.resize(num_col);
  ipx_solution.ipx_row_dual.resize(num_row);
  ipx_solution.ipx_row_status.resize(num_row);
  ipx_solution.ipx_col_status.resize(num_col);
  errflag = lps.GetBasicSolution(
      &ipx_solution.ipx_col_value[0], &ipx_solution.ipx_row_value[0],
      &ipx_solution.ipx_row_dual[0], &ipx_solution.ipx_col_dual[0],
      &ipx_solution.ipx_row_status[0], &ipx_solution.ipx_col_status[0]);
  if (errflag != 0) {
    highsLogUser(log_options, HighsLogType::kError,
                 "IPX crossover getting basic solution: flag = %d\n",
                 (int)errflag);
    return false;
  }

  // Convert the IPX basic solution to a HiGHS basic solution
  HighsStatus status = ipxBasicSolutionToHighsBasicSolution(
      options.log_options, lp, rhs, constraint_type, ipx_solution, basis,
      solution);

  if (status != HighsStatus::kOk) {
    highsLogUser(
        log_options, HighsLogType::kError,
        "Failed to convert IPX basic solution to Highs basic solution\n");
    return false;
  }

  return true;
}

#endif
