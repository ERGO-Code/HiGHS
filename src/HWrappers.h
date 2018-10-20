#ifndef HWRAPPERS_H_
#define HWRAPPERS_H_

#include "HModel.h"
#include "HConst.h"

#ifdef IPX
#include "lp_solver.h"
#endif

#ifdef IPX
int solveModelWithIpx(HModel &model)
{
  int debug = 0;

#ifdef CMAKE_BUILD_TYPE
  debug = 1;
#endif

  ipx::LpSolver lps;
  ipx::Parameters parameters;

  // parameters.crossover = 1; by default
  if (debug)
    parameters.debug = 1;

  ipx::Int num_col = model.numCol;
  ipx::Int num_row = model.numRow;

  // For each row with both a lower and an upper bound introduce one new column
  // so num_col may increase. Ignore each free row so num_row may decrease.

  // For each row with bounds on both sides introduce explicit slack and
  // transfer bounds.
  assert(model.rowLower.size() == (unsigned int) num_row);
  assert(model.rowUpper.size() == (unsigned int) num_row);
  std::vector<int> general_bounded_rows;
  std::vector<int> free_rows;

  for (int i = 0; i < num_row; i++)
    if (model.rowLower[i] > -HSOL_CONST_INF && model.rowUpper[i] < HSOL_CONST_INF)
      general_bounded_rows.push_back(i);

    // For each row (excluding free rows) add entry to char array and set up rhs
    // vector
    std::vector<double> rhs;
    std::vector<char> constraint_type;
    rhs.reserve(num_row);
    constraint_type.reserve(num_row);
    
    for (int i = 0; i < num_row; i++)
    {
      if (model.rowLower[i] > -HSOL_CONST_INF && model.rowUpper[i] == HSOL_CONST_INF)
      {
        rhs.push_back(model.rowLower[i]);
        constraint_type.push_back( '>');
      }
      else if (model.rowLower[i] == -HSOL_CONST_INF && model.rowUpper[i] < HSOL_CONST_INF)
      {
        rhs.push_back(model.rowUpper[i]);
        constraint_type.push_back('<');
      }
      else if (model.rowLower[i] == model.rowUpper[i])
      {
        rhs.push_back(model.rowUpper[i]);
        constraint_type.push_back('=');
      }
      else
      {
        free_rows.push_back(i);
      }
    }

    num_col += general_bounded_rows.size();
    num_row -= free_rows.size();

    // Copy Astart and Aindex to ipx::Int array. 
    int nnz = model.Aindex.size();
    std::vector<ipx::Int> Ap(num_col + 1);
    std::vector<ipx::Int> Ai(nnz);
    std::vector<double> Av(nnz);

    std::vector<double> col_lb();

    // Set starting points of original and newly introduced columns.
    for (int j = 0; j <= model.numCol; j++) Ap[j] = model.Astart[j];
    for (int j = model.numCol; j <= (int) num_col; j++) Ap[j] = Ap[j-1] + 1;

    for (int k = 0; k < nnz; k++)
    { 
      Ai[k] = model.Aindex[k];
      Av[k] = model.Avalue[k];
    }
    for (int k = 0; k < (int) general_bounded_rows.size(); k++)
    {
      Ai[nnz + k] = (ipx::Int) general_bounded_rows[k];
      Av[nnz + k] = 1;
    }

    // Copy row bounds and truncate if needed.
    // TODO

    // Transfer bounds from general bounded rows to slacks in rhs. 
    // TODO

    ipx::Int status = lps.Solve(
        num_col,
        &model.colCost[0],
        &model.colLower[0],
        &model.colUpper[0],
        model.numRow,
        &Ap[0],
        &Ai[0],
        &Av[0],
        &rhs[0],
        &constraint_type[0]);

   if (status != IPX_STATUS_ok)
    {
      // fatal error (invalid input, out of memory, etc.)
      std::cout << " status: " << status << ','
                << " errflag: " << lps.GetInfo().errflag << '\n';
      return 1;
    }

    // Get solver and solution information.
    ipx::Info info = lps.GetInfo();
    /*
    // Get the interior solution (available if IPM was started).
    double x[num_var], xl[num_var], xu[num_var], slack[num_constr];
    double y[num_constr], zl[num_var], zu[num_var];
    lps.GetInteriorSolution(x, xl, xu, slack, y, zl, zu);
    */

    // Get the basic solution (available if crossover terminated without error).
    double xbasic[num_col], sbasic[num_row];
    double ybasic[num_row], zbasic[num_col];
    ipx::Int cbasis[num_row], vbasis[num_col];

    lps.GetBasicSolution(xbasic, sbasic, ybasic, zbasic, cbasis, vbasis);

    // Set solution in model.
    // TODO

  return 0;
}
#endif

#endif