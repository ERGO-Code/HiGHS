#ifndef INTERIOR_POINT_IPX_WRAPPER_H_
#define INTERIOR_POINT_IPX_WRAPPER_H_

#include <algorithm>

#include "interior_point/IpxStatus.h"
#include "interior_point/ipx/include/ipx_status.h"
#include "interior_point/ipx/src/lp_solver.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"

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

IpxStatus solveModelWithIpx(const HighsLp& lp, HighsSolution& solution,
                            HighsBasis& basis    // todo: fill up HighsBasis
) {
  int debug = 0;

#ifdef CMAKE_BUILD_TYPE
  debug = 1;
#endif

  ipx::LpSolver lps;
  ipx::Parameters parameters;

  // todo: set tolerances

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
  // todo: fill up HighsSolution

  if (info.status_crossover == IPX_STATUS_optimal ||
      info.status_crossover == IPX_STATUS_imprecise) {
    if (info.status_crossover == IPX_STATUS_imprecise) {
      HighsPrintMessage(
          ML_ALWAYS,
          "Ipx Crossover status imprecise: at least one of primal and dual "
          "infeasibilities of basic solution is not within parameters pfeastol "
          "and dfeastol. Simplex clean up will be required.\n");
      // const double abs_presidual = info.abs_presidual;
      // const double abs_dresidual = info.abs_dresidual;
      // const double rel_presidual = info.rel_presidual;
      // const double rel_dresidual = info.rel_dresidual;
      // const double rel_objgap = info.rel_objgap;
    }

    std::vector<double> xbasic(num_col);
    std::vector<double> sbasic(num_row);
    std::vector<double> ybasic(num_row);
    std::vector<double> zbasic(num_col);
    std::vector<ipx::Int> cbasis(num_row);
    std::vector<ipx::Int> vbasis(num_col);

    lps.GetBasicSolution(&xbasic[0], &sbasic[0], &ybasic[0], &zbasic[0], &cbasis[0], &vbasis[0]);
    
    // Extract the basis and primal/dual solution from IPX into the
    // HighsSolution and HighsBasis
    // 
    // Check that IPX has the same number of rows and columns
    // Sometimes it won't, but can't handle that yet!
    if (lp.numCol_ != num_col || lp.numRow_ != num_row) {
      printf("Cannot convert IPX basis and solution with different numbers of columns[%d, %d] or rows[%d, %d]\n",
	     lp.numCol_, (int)num_col, lp.numRow_, (int)num_row);
      return IpxStatus::OK;
    }
    basis.col_status.resize(num_col);
    basis.row_status.resize(num_row);
    solution.col_value.resize(num_col);
    solution.col_dual.resize(num_col);
    solution.row_value.resize(num_row);
    solution.row_dual.resize(num_row);
    const ipx::Int ipx_basic = 0;
    const ipx::Int ipx_nonbasic_at_lb = -1;
    const ipx::Int ipx_nonbasic_at_ub = -2;
    double ipx_objective_value = lp.offset_;
    std::vector<double> row_activity;
    row_activity.assign(num_row, 0);
    
    int num_primal_infeasibilities = 0;
    int num_dual_infeasibilities = 0;
    double max_primal_infeasibility = 0;
    double max_dual_infeasibility = 0;
    double sum_primal_infeasibilities = 0;
    double sum_dual_infeasibilities = 0;
    for (int col=0; col<num_col; col++) {
      bool unrecognised = false;
      double primal_infeasibility = 0;
      double dual_infeasibility = 0;
      if (vbasis[col] == ipx_basic) {
	// Column is basic
	basis.col_status[col] = HighsBasisStatus::BASIC;
	solution.col_value[col] = xbasic[col];
	solution.col_dual[col] = 0;
      } else if (vbasis[col] == ipx_nonbasic_at_lb) {
	// Column is nonbasic at lower bound
	basis.col_status[col] = HighsBasisStatus::LOWER;
	solution.col_value[col] = xbasic[col];
	solution.col_dual[col] = zbasic[col];
	dual_infeasibility = max(0.0, -solution.col_dual[col]);
      } else if (vbasis[col] == ipx_nonbasic_at_ub) {
	// Column is nonbasic at upper bound
	basis.col_status[col] = HighsBasisStatus::UPPER;
	solution.col_value[col] = xbasic[col];
	solution.col_dual[col] = zbasic[col];
	dual_infeasibility = max(0.0, solution.col_dual[col]);
      } else {
	unrecognised = true;
	solution.col_value[col] = 0;
	solution.col_dual[col] = 0;
	dual_infeasibility = 1e200;
      }
      primal_infeasibility = max(lp.colLower_[col] - solution.col_value[col], solution.col_value[col] - lp.colUpper_[col]);
      if (unrecognised) printf("Unrecognised vbasis value from IPX: ");
      //      if (unrecognised)
	printf("Col %2d vbasis[%2d] = %2d; x[%2d] = %11.4g; z[%2d] = %11.4g\n",
	       col, col, (int)vbasis[col], col, xbasic[col], col, zbasic[col]);
      for (int el=lp.Astart_[col]; el<lp.Astart_[col+1]; el++) {
	int row = lp.Aindex_[el];
	row_activity[row] += solution.col_value[col]*lp.Avalue_[el];
      }
      ipx_objective_value += solution.col_value[col]*lp.colCost_[col];
      if (primal_infeasibility > 1e-7) {
	num_primal_infeasibilities++;
	max_primal_infeasibility = max(primal_infeasibility, max_primal_infeasibility);
	sum_primal_infeasibilities += primal_infeasibility;
      }
      if (dual_infeasibility > 1e-7) {
	num_dual_infeasibilities++;
	max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
	sum_dual_infeasibilities += dual_infeasibility;
      }
    }
    double norm_row_activity_error = 0;
    for (int row=0; row<num_row; row++) {
      bool unrecognised = false;
      double primal_infeasibility = 0;
      double dual_infeasibility = 0;
      if (cbasis[row] == ipx_basic) {
	// Row is basic
	basis.row_status[row] = HighsBasisStatus::BASIC;
	solution.row_value[row] = -sbasic[row];
	solution.row_dual[row] = 0;
      } else if (cbasis[row] == ipx_nonbasic_at_lb) {
	// Row is nonbasic at lower bound
	basis.row_status[row] = HighsBasisStatus::LOWER;
	solution.row_value[row] = lp.rowLower_[row]-sbasic[row];
	solution.row_dual[row] = -ybasic[row];
	dual_infeasibility = max(0.0, solution.row_dual[row]);
      } else if (cbasis[row] == ipx_nonbasic_at_ub) {
	// Row is nonbasic at upper bound
	basis.row_status[row] = HighsBasisStatus::UPPER;
	solution.row_value[row] = lp.rowUpper_[row]-sbasic[row];
	solution.row_dual[row] = -ybasic[row];
	dual_infeasibility = max(0.0, -solution.row_dual[row]);
      } else {
	unrecognised = true;
	solution.row_value[row] = 0;
	solution.row_dual[row] = 0;
	dual_infeasibility = 1e200;
      }
      if (unrecognised) printf("Unrecognised cbasis value from IPX: ");
      //      if (unrecognised)
      printf("Row %2d cbasis[%2d] = %2d; s[%2d] = %11.4g; y[%2d] = %11.4g\n",
	     row, row, (int)cbasis[row], row, sbasic[row], row, ybasic[row]);
      double row_activity_error = fabs(row_activity[row]-solution.row_value[row]);
      //      if (row_activity_error>1e-7) printf("Row activity error: %2d, %11.4g [%11.4g, %11.4g]\n", row,
      //					  row_activity_error, row_activity[row], solution.row_value[row]);  
      norm_row_activity_error += row_activity_error;
      if (primal_infeasibility > 1e-7) {
	num_primal_infeasibilities++;
	max_primal_infeasibility = max(primal_infeasibility, max_primal_infeasibility);
	sum_primal_infeasibilities += primal_infeasibility;
      }
      if (dual_infeasibility > 1e-7) {
	num_dual_infeasibilities++;
	max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
	sum_dual_infeasibilities += dual_infeasibility;
      }
    }
    printf("IPX objective =            %g\n", ipx_objective_value);
    printf("  ||row activity error|| = %g\n", norm_row_activity_error);
  }
  return IpxStatus::OK;
}
#endif
