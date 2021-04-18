#include "presolve/ICrash.h"

#include <iostream>


#ifndef IPX_ON
bool callCrossover(const HighsLp& lp, ICrashInfo& result) {return false;}
#else

#include "ipm/IpxWrapper.h"

bool callCrossover(const HighsLp& lp, ICrashInfo& result) {
  std::cout << "Calling ipx crossover after icrash...";

  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;

  IpxStatus res = fillInIpxData(lp, num_col, objective, col_lb, col_ub,
                                   num_row, Ap, Ai, Av, rhs, constraint_type);

  if (res != IpxStatus::OK) return false;

  ipx::Parameters parameters;
  parameters.crossover = true;

  ipx::LpSolver lps;
  lps.SetParameters(parameters);

  ipx::Int load_status =
      lps.LoadModel(num_col, &objective[0], &col_lb[0], &col_ub[0], num_row,
                    &Ap[0], &Ai[0], &Av[0], &rhs[0], &constraint_type[0]);
  // todo: specify primal values coming from icrash

  // run crossover 
  // lps.RunCrossover_();

  const double* x = &result.x_values[0];

  std::vector<double> collb(num_col,0);
  std::vector<double> colub(num_col,0);
  std::vector<double> zsy(num_row,0);

  for (int i=0; i<num_col; i++) {
    if (col_lb[i] == -INFINITY)
      collb[i] = INFINITY;
    if (col_ub[i] == INFINITY)
      colub[i] = INFINITY;
  }

  const double* zl = &collb[0];
  const double* zu = &colub[0];
  const double* y = &zsy[0];

  lps.LoadIPMStartingPoint(x,zl,zu,y,y,y,y);
  lps.RunCrossover_X();

  return true;
}

#endif