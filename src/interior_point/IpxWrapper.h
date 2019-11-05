#ifndef INTERIOR_POINT_IPX_WRAPPER_H_
#define INTERIOR_POINT_IPX_WRAPPER_H_

#include <algorithm>

#include "interior_point/ipx/include/ipx_status.h"
#include "interior_point/ipx/src/lp_solver.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"

enum class IpxStatus { OK, Error, ErrorFreeRow, ErrorOrNotOptimal };

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
  for (int col = 1; col <= lp.numCol_; col++)
    Ap[col] = lp.Astart_[col - 1] + sizes[col - 1];
  for (int col = lp.numCol_; col <= (int)num_col; col++) {
    Ap[col] = Ap[col - 1] + 1;
  }

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

  obj = lp.colCost_;
  obj.insert(obj.end(), num_slack, 0);
  return IpxStatus::OK;
}

IpxStatus solveModelWithIpx(const HighsLp& lp, HighsSolution& solution) {
  int debug = 0;

#ifdef CMAKE_BUILD_TYPE
  debug = 1;
#endif

  ipx::LpSolver lps;
  ipx::Parameters parameters;

  // parameters.crossover = 1; by default
  if (debug) parameters.debug = 1;

  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;
  IpxStatus result = fillInIpxData(lp, num_col, objective, col_lb, col_ub,
                                   num_row, Ap, Ai, Av, rhs, constraint_type);
  if (result != IpxStatus::OK) return result;

  ipx::Int status =
      lps.Solve(num_col, &objective[0], &col_lb[0], &col_ub[0], num_row, &Ap[0],
                &Ai[0], &Av[0], &rhs[0], &constraint_type[0]);

  if (status != IPX_STATUS_solved) {
    // fatal error (invalid input, out of memory, etc.)
    std::cout << " status: " << status << ','
              << " errflag: " << lps.GetInfo().errflag << '\n';
    return IpxStatus::ErrorOrNotOptimal;
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

  // Set solution in lp.
  // TODO

  return IpxStatus::OK;
}
#endif