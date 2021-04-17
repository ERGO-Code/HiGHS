#include "presolve/ICrash.h"

#include <iostream>

#include "ipm/IpxSolution.h"
#include "ipm/IpxStatus.h"
#include "ipm/ipx/include/ipx_status.h"
#include "ipm/ipx/src/lp_solver.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsOptions.h"

#ifndef IPX_ON
bool callCrossover(const HighsLp& lp, ICrashInfo& result) {return false;}

#else

IpxStatus fillInIpxData_(const HighsLp& lp, ipx::Int& num_col,
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
    else if (lp.rowLower_[row] <= -HIGHS_CONST_INF &&
             lp.rowUpper_[row] >= HIGHS_CONST_INF)
      free_rows.push_back(row);

  const int num_slack = general_bounded_rows.size();

  // For each row except free rows add entry to char array and set up rhs
  // vector
  rhs.reserve(num_row);
  constraint_type.reserve(num_row);

  for (int row = 0; row < num_row; row++) {
    if (lp.rowLower_[row] > -HIGHS_CONST_INF &&
        lp.rowUpper_[row] >= HIGHS_CONST_INF) {
      rhs.push_back(lp.rowLower_[row]);
      constraint_type.push_back('>');
    } else if (lp.rowLower_[row] <= -HIGHS_CONST_INF &&
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
    if (lp.colLower_[col] <= -HIGHS_CONST_INF)
      col_lb[col] = -INFINITY;
    else
      col_lb[col] = lp.colLower_[col];

    if (lp.colUpper_[col] >= HIGHS_CONST_INF)
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

bool callCrossover(const HighsLp& lp, ICrashInfo& result) {
  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;

  IpxStatus res = fillInIpxData_(lp, num_col, objective, col_lb, col_ub,
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
  lps.RunCrossover_();
  return true;
}

bool callCrossover_(const HighsLp& lp, ICrashInfo& result) {
  std::cout << "Calling ipx crossover after icrash...";
  
  // First equality problems only: qaps

  // assert(basis_);
  //  const Int m = model_.rows();
  //  const Int n = model_.cols();
  // const Vector& lb = model_.lb();
  // const Vector& ub = model_.ub();
  // basic_statuses_.clear();

  // // Construct a complementary primal-dual point from the final IPM iterate.
  // // This usually increases the residuals to Ax=b and A'y+z=c.
  // x_crossover_.resize(n+m);
  // y_crossover_.resize(m);
  // z_crossover_.resize(n+m);
  // iterate_->DropToComplementarity(x_crossover_, y_crossover_, z_crossover_);

  // // Run crossover. Perform dual pushes in increasing order and primal pushes
  // // in decreasing order of the scaling factors from the final IPM iterate.
  // {
  //     Vector weights(n+m);
  //     for (Int j = 0; j < n+m; j++)
  //         weights[j] = iterate_->ScalingFactor(j);
  //     Crossover crossover(control_);
  //     crossover.PushAll(basis_.get(), x_crossover_, y_crossover_,
  //                       z_crossover_, &weights[0], &info_);
  //     info_.time_crossover =
  //         crossover.time_primal() + crossover.time_dual();
  //     info_.updates_crossover =
  //         crossover.primal_pivots() + crossover.dual_pivots();
  //     if (info_.status_crossover != IPX_STATUS_optimal) {
  //         // Crossover failed. Discard solution.
  //         x_crossover_.resize(0);
  //         y_crossover_.resize(0);
  //         z_crossover_.resize(0);
  //         return;
  //     }
  // }

  // // Recompute vertex solution and set basic statuses.
  // basis_->ComputeBasicSolution(x_crossover_, y_crossover_, z_crossover_);
  // basic_statuses_.resize(n+m);
  // for (Int j = 0; j < basic_statuses_.size(); j++) {
  //     if (basis_->IsBasic(j)) {
  //         basic_statuses_[j] = IPX_basic;
  //     } else {
  //         if (lb[j] == ub[j])
  //             basic_statuses_[j] = z_crossover_[j] >= 0.0 ?
  //                 IPX_nonbasic_lb : IPX_nonbasic_ub;
  //         else if (x_crossover_[j] == lb[j])
  //             basic_statuses_[j] = IPX_nonbasic_lb;
  //         else if (x_crossover_[j] == ub[j])
  //             basic_statuses_[j] = IPX_nonbasic_ub;
  //         else
  //             basic_statuses_[j] = IPX_superbasic;
  //     }
  // }
  // control_.Debug()
  //     << Textline("Bound violation of basic solution:")
  //     << sci2(PrimalInfeasibility(model_, x_crossover_)) << '\n'
  //     << Textline("Dual sign violation of basic solution:")
  //     << sci2(DualInfeasibility(model_, x_crossover_, z_crossover_)) << '\n';
  // control_.Debug()
  //     << Textline("Minimum singular value of basis matrix:")
  //     << sci2(basis_->MinSingularValue()) << '\n';

  // // Declare crossover status "imprecise" if the vertex solution defined by
  // // the final basis does not satisfy tolerances.
  // model_.EvaluateBasicSolution(x_crossover_, y_crossover_, z_crossover_,
  //                              basic_statuses_, &info_);
  // if (info_.primal_infeas > control_.pfeasibility_tol() ||
  //     info_.dual_infeas > control_.dfeasibility_tol())
  //     info_.status_crossover = IPX_STATUS_imprecise;
  return false;

  return true;
}

#endif