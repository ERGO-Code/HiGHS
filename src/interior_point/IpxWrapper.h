/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file interior_point/IpxWrapper.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef INTERIOR_POINT_IPX_WRAPPER_H_
#define INTERIOR_POINT_IPX_WRAPPER_H_

#include <algorithm>

#include "interior_point/IpxStatus.h"
#include "interior_point/ipx/include/ipx_status.h"
#include "interior_point/ipx/src/lp_solver.h"
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
    Ap[col+1] = Ap[col] + sizes[col];
    //    printf("Struc Ap[%2d] = %2d; Al[%2d] = %2d\n", col, (int)Ap[col], col, (int)sizes[col]);
  }
  for (int col = lp.numCol_; col < (int)num_col; col++) {
    Ap[col+1] = Ap[col] + 1;
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

  obj = lp.colCost_;
  obj.insert(obj.end(), num_slack, 0);
#ifdef HiGHSDEV
  printf("IPX model has %d columns, %d rows and %d nonzeros\n", (int)num_col, (int)num_row, (int)Ap[num_col]);
#endif  
  /*
  for (int col = 0; col < num_col; col++)
    printf("Col %2d: [%11.4g, %11.4g] Cost = %11.4g; Start = %d\n", col, col_lb[col], col_ub[col], obj[col], (int)Ap[col]);
  for (int row = 0; row < num_row; row++) 
    printf("Row %2d: RHS = %11.4g; Type = %d\n", row, rhs[row], constraint_type[row]);
  for (int col = 0; col < num_col; col++) {
    for (int el = Ap[col]; el < Ap[col+1]; el++) {
      printf("El %2d: [%2d, %11.4g]\n", el, (int)Ai[el], Ax[el]);
    }
  }
  */

  return IpxStatus::OK;
}

IpxStatus solveModelWithIpx(const HighsLp& lp,
			    const HighsOptions& options,
			    HighsModelStatus& highs_model_status,
			    HighsInfo& highs_info,
			    HighsSolution& highs_solution,
                            HighsBasis& highs_basis) {
  int debug = 0;

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
  parameters.ipm_feasibility_tol = options.primal_feasibility_tolerance;
  parameters.ipm_optimality_tol = options.dual_feasibility_tolerance;
  // Set the internal IPX parameters
  lps.SetParameters(parameters);

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
  ipx::Info ipx_info = lps.GetInfo();
  // Struct ipx_info defined in ipx/include/ipx_info.h

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

  if (ipx_info.status_crossover == IPX_STATUS_optimal ||
      ipx_info.status_crossover == IPX_STATUS_imprecise) {
    if (ipx_info.status_crossover == IPX_STATUS_imprecise) {
      HighsPrintMessage(
          ML_ALWAYS,
          "Ipx Crossover status imprecise: at least one of primal and dual "
          "infeasibilities of basic solution is not within parameters pfeastol "
          "and dfeastol. Simplex clean up will be required.\n");
      // const double abs_presidual = ipx_info.abs_presidual;
      // const double abs_dresidual = ipx_info.abs_dresidual;
      // const double rel_presidual = ipx_info.rel_presidual;
      // const double rel_dresidual = ipx_info.rel_dresidual;
      // const double rel_objgap = ipx_info.rel_objgap;
    }

    IpxSolution ipx_solution;
    ipx_solution.num_col = num_col;
    ipx_solution.num_row = num_row;
    ipx_solution.xbasic.resize(num_col);
    ipx_solution.sbasic.resize(num_row);
    ipx_solution.ybasic.resize(num_row);
    ipx_solution.zbasic.resize(num_col);
    ipx_solution.cbasis.resize(num_row);
    ipx_solution.vbasis.resize(num_col);

    std::vector<double>& xbasic = ipx_solution.xbasic;
    std::vector<double>& sbasic = ipx_solution.sbasic;
    std::vector<double>& ybasic = ipx_solution.ybasic;
    std::vector<double>& zbasic = ipx_solution.zbasic;
    std::vector<ipx::Int>& cbasis = ipx_solution.cbasis;
    std::vector<ipx::Int>& vbasis = ipx_solution.vbasis;

    lps.GetBasicSolution(&xbasic[0], &sbasic[0], &ybasic[0], &zbasic[0], &cbasis[0], &vbasis[0]);
    
    // Extract the basis and primal/dual solution from IPX into the
    // HighsSolution and HighsBasis
    //
    // Resize the HighsSolution and HighsBasis
    highs_basis.col_status.resize(lp.numCol_);
    highs_basis.row_status.resize(lp.numRow_);
    highs_solution.col_value.resize(lp.numCol_);
    highs_solution.col_dual.resize(lp.numCol_);
    highs_solution.row_value.resize(lp.numRow_);
    highs_solution.row_dual.resize(lp.numRow_);

    ipxToHighsBasicSolution(lp, rhs, constraint_type, ipx_solution, highs_basis, highs_solution);

    /*
    // Set up meaningful names for values of vbasis and cbasis to be
    // used later in comparisons
    const ipx::Int ipx_basic = 0;
    const ipx::Int ipx_nonbasic_at_lb = -1;
    const ipx::Int ipx_nonbasic_at_ub = -2;

    // Row activities are needed to set activity values of free rows -
    // which are ignored by IPX
    vector<double> row_activity;
    bool get_row_activities = num_row < lp.numRow_;
#ifdef HiGHSDEV
    // For debugging, get the row activities if there are any boxed
    // constraints
    get_row_activities = get_row_activities || num_col > lp.numCol_;
#endif    
    if (get_row_activities) row_activity.assign(lp.numRow_, 0);
    for (int col = 0; col < lp.numCol_; col++) {
      bool unrecognised = false;
      if (vbasis[col] == ipx_basic) {
	// Column is basic
	highs_basis.col_status[col] = HighsBasisStatus::BASIC;
	highs_solution.col_value[col] = xbasic[col];
	highs_solution.col_dual[col] = 0;
      } else if (vbasis[col] == ipx_nonbasic_at_lb) {
	// Column is nonbasic at lower bound
	highs_basis.col_status[col] = HighsBasisStatus::LOWER;
	highs_solution.col_value[col] = xbasic[col];
	highs_solution.col_dual[col] = zbasic[col];
      } else if (vbasis[col] == ipx_nonbasic_at_ub) {
	// Column is nonbasic at upper bound
	highs_basis.col_status[col] = HighsBasisStatus::UPPER;
	highs_solution.col_value[col] = xbasic[col];
	highs_solution.col_dual[col] = zbasic[col];
      } else {
	unrecognised = true;
#ifdef HiGHSDEV
	printf("\nError in IPX conversion: Unrecognised value vbasis[%2d] = %d\n", col, (int)vbasis[col]);
#endif	    
      }
#ifdef HiGHSDEV
      if (unrecognised) printf("Bounds [%11.4g, %11.4g]\n", lp.colLower_[col], lp.colUpper_[col]);
      if (unrecognised)
	printf("Col %2d vbasis[%2d] = %2d; x[%2d] = %11.4g; z[%2d] = %11.4g\n",
	       col, col, (int)vbasis[col], col, xbasic[col], col, zbasic[col]);
#endif
      assert(!unrecognised);
      if (unrecognised) {
	HighsLogMessage(HighsMessageType::ERROR, "Unrecognised vbasis value from IPX");
	return IpxStatus::Error;
      }
      if (get_row_activities) {
	// Accumulate row activities to assign value to free rows
	for (int el = lp.Astart_[col]; el < lp.Astart_[col+1]; el++) {
	  int row = lp.Aindex_[el];
	  row_activity[row] += highs_solution.col_value[col]*lp.Avalue_[el];
	}
      }
    }
    int ipx_row = 0;
    int ipx_slack = lp.numCol_;
    int num_boxed_rows = 0;
    int num_boxed_rows_basic = 0;
    int num_boxed_row_slacks_basic = 0;
    for (int row = 0; row < lp.numRow_; row++) {
      bool unrecognised = false;
      double lower = lp.rowLower_[row];
      double upper = lp.rowUpper_[row];
#ifdef HiGHSDEV
      int this_ipx_row = ipx_row;
#endif      
      if (lower <= -HIGHS_CONST_INF && upper >= HIGHS_CONST_INF) {
	// Free row - removed by IPX so make it basic at its row activity
	highs_basis.row_status[row] = HighsBasisStatus::BASIC;
	highs_solution.row_value[row] = row_activity[row];
	highs_solution.row_dual[row] = 0;
      } else {
	// Non-free row, so IPX will have it
	if ((lower > -HIGHS_CONST_INF && upper < HIGHS_CONST_INF) && (lower < upper)) {
	  // Boxed row - look at its slack
	  num_boxed_rows++;
	  double slack_value = xbasic[ipx_slack];
	  double slack_dual = zbasic[ipx_slack];
	  /*
	  double row_value = rhs[ipx_row]-sbasic[ipx_row];
	  double row_dual = -ybasic[ipx_row];
	  printf("Boxed row: %2d [%11.4g, %11.4g]\n", row, lower, upper);
	  printf("   Value = %11.4g; RHS = %11.4g; Activity = %11.4g; Dual = %11.4g; cbasis = %d\n",
		 row_value, row_activity[row], rhs[ipx_row], row_dual, (int)cbasis[ipx_row]);
	  printf("   Slack = %11.4g;                Dual = %11.4g [%11.4g, %11.4g] vbasis = %d\n",
		 slack_value, slack_dual, xl[ipx_slack], xu[ipx_slack], (int)vbasis[ipx_slack]);
	  */
	  /*
	  double value = slack_value;
	  double dual = -slack_dual;
	  if (cbasis[ipx_row] == ipx_basic) {
	    // Row is basic
	    num_boxed_rows_basic++;
	    highs_basis.row_status[row] = HighsBasisStatus::BASIC;
	    highs_solution.row_value[row] = value;
	    highs_solution.row_dual[row] = 0;
	  } else if (vbasis[ipx_slack] == ipx_basic) {
	    // Slack is basic
	    num_boxed_row_slacks_basic++;
	    highs_basis.row_status[row] = HighsBasisStatus::BASIC;
	    highs_solution.row_value[row] = value;
	    highs_solution.row_dual[row] = 0;
	  } else if (vbasis[ipx_slack] == ipx_nonbasic_at_lb) {
	    // Slack at lower bound
	    highs_basis.row_status[row] = HighsBasisStatus::LOWER;
	    highs_solution.row_value[row] = value;
	    highs_solution.row_dual[row] = dual;
	  } else if (vbasis[ipx_slack] == ipx_nonbasic_at_ub) {
	    // Slack is at its upper bound
	    assert(vbasis[ipx_slack] == ipx_nonbasic_at_ub);
	    highs_basis.row_status[row] = HighsBasisStatus::UPPER;
	    highs_solution.row_value[row] = value;
	    highs_solution.row_dual[row] = dual;
	  } else {
	    unrecognised = true;
#ifdef HiGHSDEV
	    printf("\nError in IPX conversion: Row %2d (IPX row %2d) has unrecognised value vbasis[%2d] = %d\n",
		   row, ipx_row, ipx_slack, (int)vbasis[ipx_slack]);
	    double row_value = rhs[ipx_row]-sbasic[ipx_row];
	    double row_dual = -ybasic[ipx_row];
	    printf("Boxed row: %2d [%11.4g, %11.4g]\n", row, lower, upper);
	    printf("   Value = %11.4g; RHS = %11.4g; Activity = %11.4g; Dual = %11.4g; cbasis = %d\n",
		   row_value, row_activity[row], rhs[ipx_row], row_dual, (int)cbasis[ipx_row]);
	    printf("   Slack = %11.4g;                Dual = %11.4g [%11.4g, %11.4g] vbasis = %d\n",
		   slack_value, slack_dual, xl[ipx_slack], xu[ipx_slack], (int)vbasis[ipx_slack]);
#endif	    
	  }
	  // Update the slack to be used for boxed rows
	  ipx_slack++;
	} else if (cbasis[ipx_row] == ipx_basic) {
	  // Row is basic
	  highs_basis.row_status[row] = HighsBasisStatus::BASIC;
	  highs_solution.row_value[row] = rhs[ipx_row]-sbasic[ipx_row];
	  highs_solution.row_dual[row] = 0;
	} else {
	  // Nonbasic row at fixed value, lower bound or upper bound
	  assert(cbasis[ipx_row] == -1);// const ipx::Int ipx_nonbasic_row = -1;
	  double value = rhs[ipx_row]-sbasic[ipx_row];
	  double dual = -ybasic[ipx_row];
	  if (constraint_type[ipx_row] == '>') {
	    // Row is at its lower bound
	    highs_basis.row_status[row] = HighsBasisStatus::LOWER;
	    highs_solution.row_value[row] = value;
	    highs_solution.row_dual[row] = dual;
	  } else if (constraint_type[ipx_row] == '<') {
	    // Row is at its upper bound
	    highs_basis.row_status[row] = HighsBasisStatus::UPPER;
	    highs_solution.row_value[row] = value;
	    highs_solution.row_dual[row] = dual;
	  } else if (constraint_type[ipx_row] == '=') {
	    // Row is at its fixed value
	    highs_basis.row_status[row] = HighsBasisStatus::LOWER;
	    highs_solution.row_value[row] = value;
	    highs_solution.row_dual[row] = dual;
	  } else {
	    unrecognised = true;
#ifdef HiGHSDEV
	    printf("\nError in IPX conversion: Row %2d: cannot handle constraint_type[%2d] = %d\n", row, ipx_row, constraint_type[ipx_row]);
#endif	    
	  }
	}
	  // Update the IPX row index
	ipx_row++;
      }
#ifdef HiGHSDEV
      if (unrecognised) printf("Bounds [%11.4g, %11.4g]\n", lp.rowLower_[row], lp.rowUpper_[row]);
      if (unrecognised)
	printf("Row %2d cbasis[%2d] = %2d; s[%2d] = %11.4g; y[%2d] = %11.4g\n",
	       row, this_ipx_row, (int)cbasis[this_ipx_row], this_ipx_row, sbasic[this_ipx_row], this_ipx_row, ybasic[this_ipx_row]);
#endif      
      assert(!unrecognised);
      if (unrecognised) {
	HighsLogMessage(HighsMessageType::ERROR, "Unrecognised cbasis value from IPX");
	return IpxStatus::Error;
      }
    }
    assert(ipx_row == num_row);
    assert(ipx_slack == num_col);

#ifdef HiGHSDEV
    if (num_boxed_rows)
      printf("Of %d boxed rows: %d are basic and %d have basic slacks\n", 
	     num_boxed_rows, num_boxed_rows_basic, num_boxed_row_slacks_basic);
#endif
	  */
    int report_level = -1;
#ifdef HiGHSDEV
    report_level = 1;
#endif
    HighsSolutionParams solution_params;
    copyToSolutionParams(solution_params, options, ipx_info);
    highs_model_status = analyseHighsBasicSolution(lp, highs_basis, highs_solution,
						   solution_params, report_level, "after IPX");
    copyFromSolutionParams(highs_info, solution_params);
  }
  return IpxStatus::OK;
}
#endif
