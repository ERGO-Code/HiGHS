#include <iostream>

#include "presolve/ICrash.h"
#include "HConfig.h"

#ifndef IPX_ON
bool callCrossover(const HighsLp& lp, const HighsOptions& options, const std::vector<double>& x_values,
                   HighsSolution& solution, HighsBasis& highs_basis) {
  return false;
}
#else

#include "ipm/IpxWrapper.h"

bool callCrossover(const HighsLp& lp, const HighsOptions& options, const std::vector<double>& x_values,
                   HighsSolution& solution, HighsBasis& basis) {
  std::cout << "Calling ipx crossover after icrash...";

  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;

  IpxStatus res = fillInIpxData(lp, num_col, objective, col_lb, col_ub, num_row,
                                Ap, Ai, Av, rhs, constraint_type);

  if (res != IpxStatus::OK) return false;

  ipx::Parameters parameters;
  parameters.crossover = true;

  ipx::LpSolver lps;
  lps.SetParameters(parameters);

  ipx::Int load_status =
      lps.LoadModel(num_col, &objective[0], &col_lb[0], &col_ub[0], num_row,
                    &Ap[0], &Ai[0], &Av[0], &rhs[0], &constraint_type[0]);

  // run crossover
  // lps.RunCrossover_();

  // specify primal values coming from icrash
  const double* x = &x_values[0];

  std::vector<double> collb(num_col, 0);
  std::vector<double> colub(num_col, 0);
  std::vector<double> zsy(num_row, 0);

  for (int i = 0; i < num_col; i++) {
    if (col_lb[i] == -INFINITY) collb[i] = INFINITY;
    if (col_ub[i] == INFINITY) colub[i] = INFINITY;
  }

  const double* zl = &collb[0];
  const double* zu = &colub[0];
  const double* y = &zsy[0];

  const int flag = lps.LoadIPMStartingPoint(x, zl, zu, y, y, y, y);
  if (flag)
    return false;
  lps.RunCrossover_X();

  // get basis
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
  ipxBasicSolutionToHighsBasicSolution(options.logfile, lp, rhs,
                                       constraint_type, ipx_solution,
                                       basis, solution);

  std::cout << "Crossover basic solution >>>" << std::endl;

  return true;
}

#endif