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
  std::vector<double> rhs;
  std::vector<char> constraint_type;
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

  std::vector<int> truncated_row_index(num_row, -1);
  int current = 0;
  int index = 0;
  // num_row is still not decreased so the criterion below is OK.
  for (int row = 0; row < num_row; row++)
  {
    if (row != free_rows[current])
    {
      truncated_row_index[row] = index;
      index++;
    }
  }

  num_col += general_bounded_rows.size();
  num_row -= free_rows.size();

  // Copy Astart and Aindex to ipx::Int array.
  int nnz = model.Aindex.size();
  std::vector<ipx::Int> Ap(num_col + 1);
  std::vector<ipx::Int> Ai(nnz);
  std::vector<double> Av(nnz);

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
    Av[k] = model.Avalue[k];
  }
  for (int k = 0; k < (int)general_bounded_rows.size(); k++)
  {
    Ai[nnz + k] = (ipx::Int)general_bounded_rows[k];
    Av[nnz + k] = 1;
  }

  // Column bound vectors.
  std::vector<double> col_lb(num_col);
  std::vector<double> col_ub(num_col);
  for (int col = 0; col < model.numCol; col++)
  {
    if (model.colLower == -HSOL_CONST_INF)
      col_lb[col] = -INFINITY;
    else
      col_lb[col] = model.colLower[col];

    if (model.colUpper == HSOL_CONST_INF)
      col_ub[col] = INFINITY;
    else
      col_ub[col] = model.colUpper[col];
  }
  for (int slack = 0; slack < general_bounded_rows.size(); slack++)
  {
    const int row = general_bounded_rows[slack];
    col_lb[model.numCol + slack] = model.rowLower[row];
    col_ub[model.numCol + slack] = model.rowUpper[row];
  }

  // Copy row bounds and truncate if needed.
  std::vector<double> row_lb;
  std::vector<double> row_ub;

  if (free_rows.size() == 0)
  {
    row_lb = model.rowLower;
    row_ub = model.rowUpper;
  }
  else
  {
    // There is free rows to remove.
    row_lb.reserve(num_row);
    row_ub.reserve(num_row);
    int current = 0;
    int num_free = free_rows.size();
    for (int i = 0; i < model.numRow; i++)
    {
        if (current == num_free || 
            (current < num_free && i != free_rows[current])
        {
        // Row is not free.
        if (row_lb == -HSOL_CONST_INF)
        {
          row_lb.push_back(-INFINITY);
          row_ub.push_back(model.rowUpper[i]);
        }
        else if (row_ub == HSOL_CONST_INF)
        {
          row_lb.push_back(model.rowLower[i]);
          row_ub.push_back(INFINITY);
        }
        else
        {
          // Row is general bounded. Bounds are transfered onto slack.
          row_lb.push_back(0);;
          row_ub.push_back(0);;
        }
        } if (i == free_rows[current])
        {
        current++;
        }
    }
  }

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