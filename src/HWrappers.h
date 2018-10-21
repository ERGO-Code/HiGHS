#ifndef HWRAPPERS_H_
#define HWRAPPERS_H_

#include "HModel.h"
#include "HConst.h"

#ifdef IPX
#include "lp_solver.h"
#endif

#ifdef IPX

int fillInIpxData(const HModel &model, ipx::Int &num_col,
                  std::vector<double> &obj, std::vector<double> &col_lb,
                  std::vector<double> &col_ub, ipx::Int &num_row,
                  std::vector<ipx::Int> &Ap,
                  std::vector<ipx::Int> &Ai, std::vector<double> &Ax,
                  std::vector<double> &rhs,
                  std::vector<char> &constraint_type)
{
  num_col = model.numCol;
  num_row = model.numRow;

  // For each row with both a lower and an upper bound introduce one new column
  // so num_col may increase. Ignore each free row so num_row may decrease.
  // lba <= a'x <= uba becomes
  // a'x-s = 0 and lba <= s <= uba.

  // For each row with bounds on both sides introduce explicit slack and
  // transfer bounds.
  assert(model.rowLower.size() == (unsigned int)num_row);
  assert(model.rowUpper.size() == (unsigned int)num_row);
  std::vector<int> general_bounded_rows;
  std::vector<int> free_rows;

  for (int row = 0; row < num_row; row++)
    if (model.rowLower[row] > -HSOL_CONST_INF &&
        model.rowUpper[row] < HSOL_CONST_INF)
      general_bounded_rows.push_back(row);

  // For each row (excluding free rows) add entry to char array and set up rhs
  // vector
  rhs.reserve(num_row);
  constraint_type.reserve(num_row);

  for (int row = 0; row < num_row; row++)
  {
    if (model.rowLower[row] > -HSOL_CONST_INF &&
        model.rowUpper[row] == HSOL_CONST_INF)
    {
      rhs.push_back(model.rowLower[row]);
      constraint_type.push_back('>');
    }
    else if (model.rowLower[row] == -HSOL_CONST_INF &&
             model.rowUpper[row] < HSOL_CONST_INF)
    {
      rhs.push_back(model.rowUpper[row]);
      constraint_type.push_back('<');
    }
    else if (model.rowLower[row] == model.rowUpper[row])
    {
      rhs.push_back(model.rowUpper[row]);
      constraint_type.push_back('=');
    }
    else
    {
      free_rows.push_back(row);
    }
  }

  // num_row is still not decreased so it can be used as below.
  std::vector<int> truncated_row_index(num_row);
  if (free_rows.size() == 0)
  {
    for (int row = 0; row < num_row; row++)
      truncated_row_index[row] = row;
  }
  else
  {
    truncated_row_index.assign(num_row, -1);
    int current = 0;
    int index = 0;
    for (int row = 0; row < num_row; row++)
    {
      if (row != free_rows[current])
      {
        truncated_row_index[row] = index;
        index++;
      }
    }
  }

  num_col += general_bounded_rows.size();
  num_row -= free_rows.size();

  // Copy Astart and Aindex to ipx::Int array.
  int nnz = model.Aindex.size();
  Ap.resize(num_col + 1);
  Ai.resize(nnz);
  Ax.resize(nnz);

  // Set starting points of original and newly introduced columns.
  for (int col = 0; col <= model.numCol; col++)
    Ap[col] = model.Astart[col];
  for (int col = model.numCol; col <= (int)num_col; col++)
  {
    Ap[col] = Ap[col - 1] + 1;
  }
  for (int k = 0; k < nnz; k++)
  {
    Ai[k] = truncated_row_index[model.Aindex[k]];
    Ax[k] = model.Avalue[k];
  }
  for (int k = 0; k < (int)general_bounded_rows.size(); k++)
  {
    Ai[nnz + k] = (ipx::Int)general_bounded_rows[k];
    Ax[nnz + k] = 1;
  }

  // Column bound vectors.
  col_lb.reserve(num_col);
  col_ub.resize(num_col);
  for (int col = 0; col < model.numCol; col++)
  {
    if (model.colLower[col] == -HSOL_CONST_INF)
      col_lb[col] = -INFINITY;
    else
      col_lb[col] = model.colLower[col];

    if (model.colUpper[col] == HSOL_CONST_INF)
      col_ub[col] = INFINITY;
    else
      col_ub[col] = model.colUpper[col];
  }
  for (int slack = 0; slack < (int)general_bounded_rows.size(); slack++)
  {
    const int row = general_bounded_rows[slack];
    col_lb[model.numCol + slack] = model.rowLower[row];
    col_ub[model.numCol + slack] = model.rowUpper[row];
  }

  obj = model.colCost;

  return 0;
}

int solveModelWithIpx(HModel model)
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

  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;
  fillInIpxData(model, num_col, objective, col_lb, col_ub, num_row, Ap, Ai, Av,
                rhs, constraint_type);

  ipx::Int status = lps.Solve(
      num_col,
      &objective[0],
      &col_lb[0],
      &col_ub[0],
      num_row,
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