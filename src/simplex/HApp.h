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

//#include "HAPI.h"
#include "HConfig.h"
#include "HConst.h"
#include "HDual.h"
#include "HTimer.h"
#include "HighsLp.h"
#include "HighsModelObject.h"
#include "HCrash.h"
#include "HRanging.h"
#include "Scaling.h"

HighsStatus LpStatusToHighsStatus(const int lp_status) {
  switch (lp_status) {
    case LP_Status_OutOfTime:
      return HighsStatus::Timeout;
    case LP_Status_Failed:
      return HighsStatus::SolutionError;
    case LP_Status_Infeasible:
      return HighsStatus::Infeasible;
    case LP_Status_Unbounded:
      return HighsStatus::Unbounded;
    case LP_Status_Optimal:
      return HighsStatus::Optimal;
    default:
      return HighsStatus::NotImplemented;
  }
}

HighsStatus solveSimplex(const HighsOptions& opt,
                         HighsModelObject& highs_model) {
  HModel& model = highs_model.hmodel_[0];

  bool ranging = true;
  // Initialize solver.
  HDual solver;

  // If after postsolve. todo: advanced basis start here.
  if (opt.clean_up) {
    model.initFromNonbasic();

    solver.solve(highs_model);
    return LpStatusToHighsStatus(model.problemStatus);
  }

  solver.setCrash(opt.crashMode);
  // Crash, if HighsModelObject has basis information.
  if (opt.crash) {
    HCrash crash;
    crash.crash(highs_model, solver.Crash_Mode);
  }

  // Solve, depending on the options.
  // Parallel.
  if (opt.sip) {
    model.intOption[INTOPT_PERMUTE_FLAG] = 1;
    solver.solve(highs_model, HDUAL_VARIANT_TASKS, 8);
  } else if (opt.pami) {
    model.intOption[INTOPT_PERMUTE_FLAG] = 1;
    if (opt.partitionFile) {
      model.strOption[STROPT_PARTITION_FILE] = opt.partitionFile;
    }
    if (opt.cut) {
      model.dblOption[DBLOPT_PAMI_CUTOFF] = opt.cut;
    } else {
      solver.solve(highs_model, HDUAL_VARIANT_MULTI, 8);
    }
#ifdef HiGHSDEV
    if (opt.pami) model.writePivots("multi");
    if (opt.sip) model.writePivots("tasks");
#endif
  } else {
    // Serial. Based on previous solvePlainJAJH.

// todo: do setup and presolve timings elsewhere.

    double crashTime = 0; 
#ifdef HiGHSDEV
    double crossoverTime = 0;
    double presolve2Time = 0;
#endif
    double solveTime = 0;
    int solveIt = 0;
#ifdef HiGHSDEV
    int solvePh1DuIt = 0;
    int solvePh2DuIt = 0;
    int solvePrIt = 0;
#endif
    double lcSolveTime;

    HDual solver;

    vector<double> colPrAct;
    vector<double> colDuAct;
    vector<double> rowPrAct;
    vector<double> rowDuAct;

    solver.setPrice(opt.priceMode);
    solver.setEdWt(opt.edWtMode);
    solver.setTimeLimit(opt.timeLimit);

    model.timer.reset();

    //  bool FourThreads = true;
    bool FourThreads = false;
    //  bool EightThreads = true;
    bool EightThreads = false;

    if (FourThreads)
      solver.solve(highs_model, HDUAL_VARIANT_MULTI, 4);
    else if (EightThreads)
      solver.solve(highs_model, HDUAL_VARIANT_MULTI, 8);
    else
      solver.solve(highs_model);

    lcSolveTime = model.timer.getTime();
    solveTime += lcSolveTime;
    solveIt += model.numberIteration;

#ifdef HiGHSDEV
    solvePh1DuIt += solver.n_ph1_du_it;
    solvePh2DuIt += solver.n_ph2_du_it;
    solvePrIt += solver.n_pr_it;
    printf(
        "\nBnchmkHsol01 After presolve        ,hsol,%3d,%16s, %d,%d,"
        "%10.3f,%20.10e,%10d,%10d,%10d\n",
        model.getPrStatus(), model.modelName.c_str(), highs_model.lp_.numRow_,
        highs_model.lp_.numCol_, lcSolveTime, model.dualObjectiveValue, solver.n_ph1_du_it,
        solver.n_ph2_du_it, solver.n_pr_it);
#endif

    // Possibly recover bounds after presolve (after using bounds tightened by
    // presolve)
    if (model.usingImpliedBoundsPresolve) {
      //		Recover the true bounds overwritten by the implied
      // bounds
#ifdef HiGHSDEV
      printf("\nRecovering bounds after using implied bounds and resolving\n");
#endif
      if (model.problemStatus != LP_Status_OutOfTime) {
        model.copy_savedBoundsToModelBounds();

        model.timer.reset();
        solver.solve(highs_model);
 /*       lcSolveTime = model.timer.getTime();
        solveTime += lcSolveTime;
        solveIt += model.numberIteration;  */
        model.util_reportSolverOutcome("After recover:   ");
#ifdef HiGHSDEV
        solvePh1DuIt += solver.n_ph1_du_it;
        solvePh2DuIt += solver.n_ph2_du_it;
        solvePrIt += solver.n_pr_it;
        printf(
            "\nBnchmkHsol02 After restoring bounds,hsol,%3d,%16s, %d,%d,"
            "%10.3f,%20.10e,%10d,%10d,%10d\n",
            model.getPrStatus(), model.modelName.c_str(), highs_model.lp_.numRow_,
            highs_model.lp_.numCol_, lcSolveTime, model.dualObjectiveValue, solver.n_ph1_du_it,
            solver.n_ph2_du_it, solver.n_pr_it);
#endif
      }
    }
    //    HighsUtils utils; utils.reportLp(highs_model.lp_);
    //    utils.reportLpSolution(highs_model);
    HighsStatus result = LpStatusToHighsStatus(model.problemStatus);
    if (result != HighsStatus::Optimal) return result;


/* todo: do elsewhere once timing is added.
#ifdef HiGHSDEV
    double sumTime =
        crashTime + solveTime;
        // setupTime + presolve1Time + crashTime + solveTime + postsolveTime;
    printf(
        "Time: setup = %10.3f; presolve = %10.3f; crash = %10.3f; solve = "
        "%10.3f; postsolve = %10.3f; sum = %10.3f; total = %10.3f\n",
        setupTime, presolve1Time, crashTime, solveTime, postsolveTime, sumTime,
        model.totalTime);
    cout << flush;
    double errTime = abs(sumTime - model.totalTime);
    if (errTime > 1e-3) printf("!! Sum-Total time error of %g\n", errTime);
#endif
*/

    // TODO Reinstate this once solve after postsolve is performed
    //  model.util_getPrimalDualValues(colPrAct, colDuAct, rowPrAct, rowDuAct);
    //  double Ph2Objective = model.computePh2Objective(colPrAct);
    //  printf("Computed Phase 2 objective = %g\n", Ph2Objective);

/* todo: do elsewhere once timing is added.
#ifdef HiGHSDEV
    bool rpBnchmk = false;
    if (rpBnchmk) {
      int numCol = highs_model.lp_.numCol_;
      int numRow = highs_model.lp_.numRow_;
      printf(
          "\nBnchmkHsol99,hsol,%3d,%16s,Presolve %s,"
          "Crash %s,EdWt %s,Price %s,%d,%d,%10.3f,%10.3f,"
          "%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,"
          "%20.10e,%10d,%10.3f,"
          "%d\n",
          model.getPrStatus(), model.modelName.c_str(), Presolve_ArgV,
          Crash_ArgV, EdWt_ArgV, Price_ArgV, numRow, numCol, setupTime,
          presolve1Time, crashTime, crossoverTime, presolve2Time, solveTime,
          postsolveTime, model.dualObjective, model.numberIteration,
          model.totalTime, solver.n_wg_DSE_wt);
      cout << flush;
    }
#endif
*/

  }
  return HighsStatus::Optimal;
}

HighsStatus solveScip(const HighsOptions& opt, HighsModelObject& highs_model) {
  printf("Called solveScip.\n");

  HModel& model = highs_model.hmodel_[0];
  cout << flush;
  //  model.util_reportModel();

  // Extract columns numCol-3..numCol-1
  int FmCol = highs_model.lp_.numCol_ - 3;
  int ToCol = highs_model.lp_.numCol_ - 1;
  int numExtractCols = ToCol - FmCol + 1;
  vector<double> XcolCost;
  vector<double> XcolLower;
  vector<double> XcolUpper;
  vector<int> XAstart;
  vector<int> XAindex;
  vector<double> XAvalue;
  //  model.util_extractCols(FmCol, ToCol, XcolCost, XcolLower, XcolUpper,
  //			 XAstart, XAindex, XAvalue);

  //  printf("Returned from model.util_extractCols with\n");
  //  model.util_reportColVec(numExtractCols, XcolCost, XcolLower, XcolUpper);
  //  model.util_reportColMtx(numExtractCols, XAstart, XAindex, XAvalue);

  // Delete the columns just extracted
  model.util_deleteCols(FmCol, ToCol);
  //  model.util_reportModel();

  // Extract rows numRow-3..numRow-1
  int FmRow = highs_model.lp_.numRow_ - 3;
  int ToRow = highs_model.lp_.numRow_ - 1;
  int numExtractRows = ToRow - FmRow + 1;
  vector<double> XrowLower;
  vector<double> XrowUpper;
  vector<int> XARstart;
  vector<int> XARindex;
  vector<double> XARvalue;
  //  model.util_extractRows(FmRow, ToRow, &(*XrowLower.begin()),
  //  &(*XrowUpper.begin()),
  // &(*XARstart.begin()), &(*XARindex.begin()), &(*XARvalue.begin()));

  //  printf("Returned from model.util_extractRows with\n");
  //  model.util_reportRowVec(numExtractRows, XrowLower, XrowUpper);
  //  model.util_reportRowMtx(numExtractRows, XARstart, XARindex, XARvalue);

  // Delete the rows just extracted
  model.util_deleteRows(FmRow, ToRow);
  //  model.util_reportModel();

  // Extract all remaining rows
  FmRow = 0;
  ToRow = highs_model.lp_.numRow_ - 1;
  int num0ExtractRows = ToRow - FmRow + 1;
  vector<double> X0rowLower;
  vector<double> X0rowUpper;
  vector<int> X0ARstart;
  vector<int> X0ARindex;
  vector<double> X0ARvalue;

  // model.util_extractRows(FmRow, ToRow, &(*X0rowLower.begin()),
  // &(*X0rowUpper.begin()),
  //			 &(*X0ARstart.begin()), &(*X0ARindex.begin()),
  //&(*X0ARvalue.begin()));

  // Delete the rows just extracted
  model.util_deleteRows(FmRow, ToRow);
  //  model.util_reportModel();

  // Extract all remaining columns
  FmCol = 0;
  ToCol = highs_model.lp_.numCol_ - 1;
  int num0ExtractCols = ToCol - FmCol + 1;
  vector<double> X0colCost;
  vector<double> X0colLower;
  vector<double> X0colUpper;
  vector<int> X0Astart;
  vector<int> X0Aindex;
  vector<double> X0Avalue;
  //  model.util_extractCols(FmCol, ToCol, X0colCost, X0colLower, X0colUpper,
  //			 X0Astart, X0Aindex, X0Avalue);

  // Delete the columns just extracted
  model.util_deleteCols(FmCol, ToCol);
  //  model.util_reportModel();

  int nnonz = 0;
  model.util_addCols(num0ExtractCols, &X0colCost[0], &X0colLower[0],
                     &X0colUpper[0], nnonz, &X0Astart[0], &X0Aindex[0],
                     &X0Avalue[0]);
  //  model.util_reportModel();

  nnonz = X0ARstart[num0ExtractRows];
  model.util_addRows(num0ExtractRows, &X0rowLower[0], &X0rowUpper[0], nnonz,
                     &X0ARstart[0], &X0ARindex[0], &X0ARvalue[0]);
  //  model.util_reportModel();

  nnonz = XARstart[numExtractRows];
  model.util_addRows(numExtractRows, &XrowLower[0], &XrowUpper[0], nnonz,
                     &XARstart[0], &XARindex[0], &XARvalue[0]);
  //  model.util_reportModel();

  nnonz = XAstart[numExtractCols];
  model.util_addCols(numExtractCols, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                     nnonz, &XAstart[0], &XAindex[0], &XAvalue[0]);
  //  model.util_reportModel();

  model.scaleModel();
  HDual solver;
  solver.solve(highs_model);
  model.util_reportModelSolution(model.lp_scaled_);
  model.util_reportSolverOutcome("SCIP 1");

  vector<double> colPrimal(highs_model.lp_.numCol_);
  vector<double> colDual(highs_model.lp_.numCol_);
  vector<double> colLower(highs_model.lp_.numCol_);
  vector<double> colUpper(highs_model.lp_.numCol_);
  vector<double> rowPrimal(highs_model.lp_.numRow_);
  vector<double> rowDual(highs_model.lp_.numRow_);
  vector<double> rowLower(highs_model.lp_.numRow_);
  vector<double> rowUpper(highs_model.lp_.numRow_);
  model.util_getPrimalDualValues(colPrimal, colDual, rowPrimal, rowDual);
  model.util_getColBounds(model.lp_scaled_, 0, highs_model.lp_.numCol_ - 1, &colLower[0], &colUpper[0]);
  model.util_getRowBounds(model.lp_scaled_, 0, highs_model.lp_.numRow_ - 1, &rowLower[0], &rowUpper[0]);

  double og_colLower;
  double og_colUpper;
  int colBoundIndex;
  double nw_colLower;
  double nw_colUpper;

  int num_resolve = 0;
  for (int col = 0; col < highs_model.lp_.numCol_; col++) {
    model.util_getColBounds(model.lp_scaled_, col, col, &og_colLower, &og_colUpper);
    printf("\nColumn %2d has primal value %11g and bounds [%11g, %11g]", col,
           colPrimal[col], og_colLower, og_colUpper);
    if (model.basis_->nonbasicFlag_[col]) {
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
      printf("Calling model.util_chgColBounds(1, %d, %g, %g)\n", colBoundIndex,
             nw_colLower, nw_colUpper);
      model.util_chgColBoundsSet(1, &colBoundIndex, &nw_colLower, &nw_colUpper);
      printf("Calling model.scaleModel()\n");
      model.scaleModel();
      //      printf("Calling solver.solve(highs_model)\n");
      solver.solve(highs_model);
      //      printf("Called solver.solve(highs_model)\n");
      model.util_reportSolverOutcome("SCIP 2");
      // Was &nw_colLower, &nw_colUpper); and might be more interesting for
      // avgas
      model.util_chgColBoundsSet(1, &colBoundIndex, &og_colLower, &og_colUpper);
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
                             HighsModelObject& highs_model) {
  cout << "=================================================================="
       << endl;

  // For the moment handle scip case separately.
  if (opt.scip) return solveScip(opt, highs_model);

  // When runSimplexSolver is called initialize an instance of HModel inside the
  // HighsModelObject. This will then be passed to HDual.
  highs_model.hmodel_.push_back(HModel());

  HModel& model = highs_model.hmodel_[0];
  const HighsLp& lp_ = highs_model.lp_;

  // Allocate memory for the basis
  // assignBasis();
  const int numTot = highs_model.lp_.numCol_ + highs_model.lp_.numRow_;
  highs_model.basis_.basicIndex_.resize(highs_model.lp_.numRow_);
  highs_model.basis_.nonbasicFlag_.assign(numTot, 0);
  highs_model.basis_.nonbasicMove_.resize(numTot);

  // Set pointers within HModel for the basis, scaling data structure and simplex information data structure
  model.basis_ = &highs_model.basis_;
  model.scale_ = &highs_model.scale_;
  model.simplex_ = &highs_model.simplex_;

  bool load_fromArrays = false;
  if (load_fromArrays) {
    model.load_fromArrays(lp_.numCol_, lp_.sense_, &lp_.colCost_[0],
			  &lp_.colLower_[0], &lp_.colUpper_[0], lp_.numRow_,
			  &lp_.rowLower_[0], &lp_.rowUpper_[0], lp_.nnz_,
			  &lp_.Astart_[0], &lp_.Aindex_[0], &lp_.Avalue_[0]);
    // Scaling: Separate from simplex.
    model.scaleModel();

  } else {
    // Copy the LP to the structure to be scaled and then scale it
    scaleHighsModel(highs_model);
    model.lp_scaled_ = highs_model.lp_scaled_;
    model.initWithLogicalBasis();

  }

  const HighsLp& lp_scaled_ = model.lp_scaled_;
  highs_model.matrix_.setup_lgBs(lp_scaled_.numCol_, lp_scaled_.numRow_,
				  &lp_scaled_.Astart_[0],
				  &lp_scaled_.Aindex_[0],
				  &lp_scaled_.Avalue_[0]);

  highs_model.factor_.setup(lp_scaled_.numCol_, lp_scaled_.numRow_,
			   &lp_scaled_.Astart_[0],
			   &lp_scaled_.Aindex_[0],
			   &lp_scaled_.Avalue_[0],
			   &highs_model.basis_.basicIndex_[0]);

  // Set pointers within HModel for the matrix and factor data structure
  model.matrix_ = &highs_model.matrix_;
  model.factor_ = &highs_model.factor_;


  // Crash, if needed.

  HighsStatus result = solveSimplex(opt, highs_model);

  // todo uncomment line below.
  if (result != HighsStatus::Optimal) return result;

  // HighsSolution set values in highs_model.
  HighsSolution& solution = highs_model.solution_;
  highs_model.hmodel_[0].util_getPrimalDualValues(
      solution.colValue_, solution.colDual_, solution.rowValue_,
      solution.rowDual_);
  model.util_getBasicIndexNonbasicFlag(highs_model.basis_info_.basis_index,
                                       highs_model.basis_info_.nonbasic_flag);

  highs_model.basis_info_.nonbasic_move = model.basis_->nonbasicMove_;

  cout << "==================================================================\n";

  return result;
}

#endif
