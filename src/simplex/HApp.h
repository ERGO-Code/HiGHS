/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef SIMPLEX_HAPP_H_
#define SIMPLX_HAPP_H_

// todo: clear includes.
#include <getopt.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "HConfig.h"
#include "HCrash.h"
#include "HDual.h"
#include "HighsLp.h"
#include "HighsLpUtils.h"
#include "HighsModelObject.h"
#include "HighsStatus.h"
#include "HighsUtils.h"
//#include "HRanging.h"
#include "HSimplex.h"
#include "HighsSimplexInterface.h"
#include "SimplexConst.h"

using std::cout;
using std::endl;
using std::flush;

HighsStatus LpStatusToHighsStatus(SimplexSolutionStatus simplex_solution_status) {
  switch (simplex_solution_status) {
  case SimplexSolutionStatus::OUT_OF_TIME:
      return HighsStatus::Timeout;
  case SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      return HighsStatus::ReachedDualObjectiveUpperBound;
  case SimplexSolutionStatus::FAILED:
      return HighsStatus::SolutionError;
  case SimplexSolutionStatus::SINGULAR:
      return HighsStatus::SolutionError;
  case SimplexSolutionStatus::UNBOUNDED:
      return HighsStatus::Unbounded;
  case SimplexSolutionStatus::INFEASIBLE:
      return HighsStatus::Infeasible;
  case SimplexSolutionStatus::OPTIMAL:
      return HighsStatus::Optimal;
  default:
      return HighsStatus::NotImplemented;
  }
}

HighsStatus solveSimplex(
			 const HighsOptions& opt,
                         HighsModelObject& highs_model_object
			 ) {
  // Just solves the LP in highs_model_object.scaled_lp_
  HighsTimer &timer = highs_model_object.timer_;
  HighsSimplexInfo &simplex_info_ = highs_model_object.simplex_info_;

  bool ranging = true;
  // Initialize solver and set dual solver options from simplex options
  HDual dual_solver(highs_model_object);
  dual_solver.options();
  
  // If after postsolve. todo: advanced basis start here.
  if (opt.clean_up) {
    initialise_from_nonbasic(highs_model_object); // initFromNonbasic();
    timer.start(timer.solve_clock);
    dual_solver.solve();
    timer.stop(timer.solve_clock);
    return LpStatusToHighsStatus(simplex_info_.solution_status);
  }

  // Crash, if HighsModelObject has basis information.
  if (simplex_info_.crash_strategy != SimplexCrashStrategy::OFF) {
    HCrash crash;
    timer.start(timer.crash_clock);
    crash.crash(highs_model_object, 0);
    timer.stop(timer.crash_clock);
  }

  timer.start(timer.solve_clock);
  // Solve, depending on the options.
  // Parallel.
  if (simplex_info_.simplex_strategy == SimplexStrategy::DUAL_TASKS) {
    dual_solver.solve(8);
  } else if (simplex_info_.simplex_strategy == SimplexStrategy::DUAL_MULTI) {
    //    if (opt.partitionFile.size() > 0) {model.strOption[STROPT_PARTITION_FILE] = opt.partitionFile;}
    dual_solver.solve(8);
#ifdef HiGHSDEV
    // Reinstate this once simplex::writePivots is written
    //    if (simplex_info_.simplex_strategy == SimplexStrategy::DUAL_MULTI) writePivots("multi");
    //    if (simplex_info_.simplex_strategy == SimplexStrategy::DUAL_TASKS) writePivots("tasks");
#endif
  } else {
    // Serial. Based on previous solvePlainJAJH.

    //    double lcSolveTime;

    //  bool FourThreads = true;
    bool FourThreads = false;
    //  bool EightThreads = true;
    bool EightThreads = false;

    if (FourThreads)
      dual_solver.solve(4);
    else if (EightThreads)
      dual_solver.solve(8);
    else
      dual_solver.solve();

#ifdef HiGHSDEV
    double currentRunHighsTime = highs_model_object.timer_.readRunHighsClock();
    printf(
        "\nBnchmkHsol01 After presolve        ,hsol,%3d,%16s, %d,%d,"
        "%10.3f,%20.10e,%10d,%10d,%10d\n",
        (int) simplex_info_.solution_status,
	highs_model_object.lp_.model_name_.c_str(),
	highs_model_object.lp_.numRow_,
        highs_model_object.lp_.numCol_,
	currentRunHighsTime,
	simplex_info_.dualObjectiveValue,
	simplex_info_.dual_phase1_iteration_count,
        simplex_info_.dual_phase2_iteration_count,
	simplex_info_.primal_phase1_iteration_count);
#endif
    //    reportLp(highs_model_object.lp_);
    //    reportLpSolution(highs_model_object);
    HighsStatus result = LpStatusToHighsStatus(simplex_info_.solution_status);

    timer.stop(timer.solve_clock);


    if (result != HighsStatus::Optimal) return result;


    // TODO Reinstate this once solve after postsolve is performed
    vector<double> colPrAct;
    vector<double> colDuAct;
    vector<double> rowPrAct;
    vector<double> rowDuAct;

    HighsSimplexInterface simplex_interface(highs_model_object);
    simplex_interface.get_primal_dual_values(colPrAct, colDuAct, rowPrAct, rowDuAct);
    double lp_objective_value = simplex_interface.get_lp_objective_value(colPrAct);
#ifdef HiGHSDEV
    printf("Computed LP objective value = %g\n", lp_objective_value);
#endif

  }
  return HighsStatus::Optimal;
}

HighsStatus solveScip(const HighsOptions& opt, HighsModelObject& highs_model_object) {
  printf("Called solveScip.\n");

  const HighsLp &lp = highs_model_object.lp_;
  HighsBasis &basis = highs_model_object.basis_;

  HighsSimplexInterface simplex_interface(highs_model_object);

  // Extract columns numCol-3..numCol-1
  int FmCol = lp.numCol_ - 3;
  int ToCol = lp.numCol_ - 1;
  int numExtractCols = ToCol - FmCol + 1;
  vector<double> XcolCost;
  vector<double> XcolLower;
  vector<double> XcolUpper;
  vector<int> XAstart;
  vector<int> XAindex;
  vector<double> XAvalue;
  //  simplex_interface.util_extractCols(FmCol, ToCol, XcolCost, XcolLower, XcolUpper,
  //			 XAstart, XAindex, XAvalue);

  //  printf("Returned from simplex_interface.util_extractCols with\n");
  //  simplex_interface.util_reportColVec(numExtractCols, XcolCost, XcolLower, XcolUpper);
  //  simplex_interface.util_reportColMtx(numExtractCols, XAstart, XAindex, XAvalue);

  // Delete the columns just extracted
  simplex_interface.util_delete_cols(FmCol, ToCol);
  //  simplex_interface.util_reportModel();

  // Extract rows numRow-3..numRow-1
  int FmRow = lp.numRow_ - 3;
  int ToRow = lp.numRow_ - 1;
  int numExtractRows = ToRow - FmRow + 1;
  vector<double> XrowLower;
  vector<double> XrowUpper;
  vector<int> XARstart;
  vector<int> XARindex;
  vector<double> XARvalue;
  //  simplex_interface.util_extract_rows(FmRow, ToRow, &(*XrowLower.begin()),
  //  &(*XrowUpper.begin()),
  // &(*XARstart.begin()), &(*XARindex.begin()), &(*XARvalue.begin()));

  //  printf("Returned from simplex_interface.util_extractRows with\n");
  //  simplex_interface.util_reportRowVec(numExtractRows, XrowLower, XrowUpper);
  //  report_row_matrix(numExtractRows, XARstart, XARindex, XARvalue);

  // Delete the rows just extracted
  simplex_interface.util_delete_rows(FmRow, ToRow);
  //  simplex_interface.util_reportModel();

  // Extract all remaining rows
  FmRow = 0;
  ToRow = lp.numRow_ - 1;
  int num0ExtractRows = ToRow - FmRow + 1;
  vector<double> X0rowLower;
  vector<double> X0rowUpper;
  vector<int> X0ARstart;
  vector<int> X0ARindex;
  vector<double> X0ARvalue;

  // simplex_interface.util_extract_rows(FmRow, ToRow, &(*X0rowLower.begin()),
  // &(*X0rowUpper.begin()),
  //			 &(*X0ARstart.begin()), &(*X0ARindex.begin()),
  //&(*X0ARvalue.begin()));

  // Delete the rows just extracted
  simplex_interface.util_delete_rows(FmRow, ToRow);
  //  simplex_interface.util_reportModel();

  // Extract all remaining columns
  FmCol = 0;
  ToCol = lp.numCol_ - 1;
  int num0ExtractCols = ToCol - FmCol + 1;
  vector<double> X0colCost;
  vector<double> X0colLower;
  vector<double> X0colUpper;
  vector<int> X0Astart;
  vector<int> X0Aindex;
  vector<double> X0Avalue;
  //  simplex_interface.util_extract_cols(FmCol, ToCol, X0colCost, X0colLower, X0colUpper,
  //			 X0Astart, X0Aindex, X0Avalue);

  // Delete the columns just extracted
  simplex_interface.util_delete_cols(FmCol, ToCol);
  //  simplex_interface.util_reportModel();

  int nnonz = 0;
  simplex_interface.util_add_cols(num0ExtractCols, &X0colCost[0], &X0colLower[0],
                     &X0colUpper[0], nnonz, &X0Astart[0], &X0Aindex[0],
                     &X0Avalue[0]);
  //  simplex_interface.util_reportModel();

  nnonz = X0ARstart[num0ExtractRows];
  simplex_interface.util_add_rows(num0ExtractRows, &X0rowLower[0], &X0rowUpper[0], nnonz,
                     &X0ARstart[0], &X0ARindex[0], &X0ARvalue[0]);
  //  simplex_interface.util_reportModel();

  nnonz = XARstart[numExtractRows];
  simplex_interface.util_add_rows(numExtractRows, &XrowLower[0], &XrowUpper[0], nnonz,
                     &XARstart[0], &XARindex[0], &XARvalue[0]);
  //  simplex_interface.util_reportModel();

  nnonz = XAstart[numExtractCols];
  simplex_interface.util_add_cols(numExtractCols, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                     nnonz, &XAstart[0], &XAindex[0], &XAvalue[0]);
  //  simplex_interface.util_reportModel();

  scale_solver_lp(highs_model_object);
  HDual dual_solver(highs_model_object);
  dual_solver.solve();
  //  reportLpSolution(highs_model_object);
  simplex_interface.report_simplex_outcome("SCIP 1");

  vector<double> colPrimal(lp.numCol_);
  vector<double> colDual(lp.numCol_);
  vector<double> colLower(lp.numCol_);
  vector<double> colUpper(lp.numCol_);
  vector<double> rowPrimal(lp.numRow_);
  vector<double> rowDual(lp.numRow_);
  vector<double> rowLower(lp.numRow_);
  vector<double> rowUpper(lp.numRow_);
  simplex_interface.get_primal_dual_values(colPrimal, colDual, rowPrimal, rowDual);
  getLpColBounds(highs_model_object.solver_lp_, 0, lp.numCol_ - 1, &colLower[0], &colUpper[0]);
  getLpRowBounds(highs_model_object.solver_lp_, 0, lp.numRow_ - 1, &rowLower[0], &rowUpper[0]);

  double og_colLower;
  double og_colUpper;
  int colBoundIndex;
  double nw_colLower;
  double nw_colUpper;

  int num_resolve = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    getLpColBounds(lp, col, col, &og_colLower, &og_colUpper);
    printf("\nColumn %2d has primal value %11g and bounds [%11g, %11g]", col,
           colPrimal[col], og_colLower, og_colUpper);
    if (basis.nonbasicFlag_[col]) {
      printf(": nonbasic so don't branch\n");
      continue;
    } else {
      double rsdu =
          min(colPrimal[col] - og_colLower, og_colUpper - colPrimal[col]);
      if (rsdu < 0.1) {
        printf(": basic but rsdu = %11g so don't branch\n", rsdu);
        continue;
      }
      printf(": basic with rsdu = %11g so branch\n\n", rsdu);
      num_resolve++;
      colBoundIndex = col;
      if (highs_isInfinity(og_colUpper))
        nw_colLower = colPrimal[col] + 1;
      else
        nw_colLower = og_colUpper;
      nw_colUpper = og_colUpper;
      printf("Calling simplex_interface.change_col_bounds_set(1, %d, %g, %g)\n", colBoundIndex, nw_colLower, nw_colUpper);
      simplex_interface.change_col_bounds_set(1, &colBoundIndex, &nw_colLower, &nw_colUpper);
      scale_solver_lp(highs_model_object);
      dual_solver.solve();
      simplex_interface.report_simplex_outcome("SCIP 2");
      // Was &nw_colLower, &nw_colUpper); and might be more interesting for
      // avgas
      simplex_interface.change_col_bounds_set(1, &colBoundIndex, &og_colLower, &og_colUpper);
      if (num_resolve >= 10) break;
    }
  }
  printf("Returning from solveSCIP\n");
  cout << flush;
  return HighsStatus::OK;
}

// Single function to solve an lp according to options and fill
// solution in solution.
HighsStatus runSimplexSolver(const HighsOptions& opt,
                             HighsModelObject& highs_model_object) {
  // For the moment handle scip case separately.
  if (opt.scip) return solveScip(opt, highs_model_object);

  HighsTimer &timer = highs_model_object.timer_;

  // Set up aliases
  const HighsLp &lp_ = highs_model_object.lp_;
  HighsBasis &basis_ = highs_model_object.basis_;
  HighsScale &scale_ = highs_model_object.scale_;
  HighsLp &solver_lp_ = highs_model_object.solver_lp_;
  HighsSimplexInfo &simplex_info_ = highs_model_object.simplex_info_;
  HMatrix &matrix_ = highs_model_object.matrix_;
  HFactor &factor_ = highs_model_object.factor_;

  // Copy the LP to the structure to be used by the solver
  solver_lp_ = lp_;

  // Set simplex options from HiGHS options
  options(highs_model_object, opt);

  // Possibly transpose the LP to be solved. This will change the
  // numbers of rows and columns in the LP to be solved
  if (simplex_info_.transpose_solver_lp) transpose_solver_lp(highs_model_object);

  // Now that the numbers of rows and columns in the LP to be solved
  // are fixed, initialise the real and integer random vectors
  initialise_solver_lp_random_vectors(highs_model_object);
  //
  // Allocate memory for the basis
  // assignBasis();
  const int numTot = highs_model_object.lp_.numCol_ + highs_model_object.lp_.numRow_;
  basis_.basicIndex_.resize(highs_model_object.lp_.numRow_);
  basis_.nonbasicFlag_.assign(numTot, 0);
  basis_.nonbasicMove_.resize(numTot);
  //
  // Possibly scale the LP to be used by the solver
  //
  // Initialise unit scaling factors, to simplify things is no scaling
  // is performed
  scaleHighsModelInit(highs_model_object);
  if (simplex_info_.scale_solver_lp)
    scale_solver_lp(highs_model_object);
  //
  // Possibly permute the columns of the LP to be used by the solver. 
  if (simplex_info_.permute_solver_lp)
    permute_solver_lp(highs_model_object);
  //
  // Possibly tighten the bounds of LP to be used by the solver. 
  if (simplex_info_.tighten_solver_lp)
    tighten_solver_lp(highs_model_object);
  //
#ifdef HiGHSDEV
  // Analyse the scaled LP
  if (simplex_info_.analyseLp) {
    util_analyseLp(lp_, "Unscaled");
    if (simplex_info_.solver_lp_is_scaled) {
      util_analyseVectorValues("Column scaling factors", lp_.numCol_, scale_.col_, false);
      util_analyseVectorValues("Row    scaling factors", lp_.numRow_, scale_.row_, false);
      util_analyseLp(solver_lp_, "Scaled");
    }
  }
  //  report_solver_lp_status_flags(highs_model_object);
#endif

  initialise_with_logical_basis(highs_model_object); // initWithLogicalBasis();

  matrix_.setup_lgBs(solver_lp_.numCol_, solver_lp_.numRow_,
		     &solver_lp_.Astart_[0],
		     &solver_lp_.Aindex_[0],
		     &solver_lp_.Avalue_[0]);
  
  factor_.setup(solver_lp_.numCol_, solver_lp_.numRow_,
		&solver_lp_.Astart_[0],
		&solver_lp_.Aindex_[0],
		&solver_lp_.Avalue_[0],
		&basis_.basicIndex_[0]);

  // Crash, if needed.

  HighsStatus result = solveSimplex(opt, highs_model_object);

  // todo uncomment line below.
  if (result != HighsStatus::Optimal) return result;

  // HighsSolution set values in highs_model_object.
  HighsSolution& solution = highs_model_object.solution_;
  HighsSimplexInterface simplex_interface(highs_model_object);
  simplex_interface.get_primal_dual_values(solution.colValue_,
					   solution.colDual_,
					   solution.rowValue_,
					   solution.rowDual_);
  simplex_interface.get_basicIndex_nonbasicFlag(highs_model_object.basis_info_.basis_index,
						highs_model_object.basis_info_.nonbasic_flag);

  highs_model_object.basis_info_.nonbasic_move = basis_.nonbasicMove_;

  return result;
}

#endif
