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

#include "HAPI.h"
#include "HConfig.h"
#include "HConst.h"
#include "HDual.h"
#include "HTimer.h"
#include "HighsLp.h"
#include "HighsModelObject.h"

HighsStatus solveSimplex(const HighsOptions& opt,
                         HighsModelObject& highs_model) {
  HModel& model = highs_model.hmodel_[0];

  // todo: cange with Julian's HighsModelObject parameter.
  HDual solver;

  // Solve, depending on the options.
  // Parallel.
  if (opt.sip) {
    model.intOption[INTOPT_PERMUTE_FLAG] = 1;
    solver.solve(&model, HDUAL_VARIANT_TASKS, 8);
    model.util_reportSolverOutcome("Solve tasks");
  } else if (opt.pami) {
    model.intOption[INTOPT_PERMUTE_FLAG] = 1;
    if (opt.partitionFile) {
      model.strOption[STROPT_PARTITION_FILE] = opt.partitionFile;
    }
    if (opt.cut) {
      model.dblOption[DBLOPT_PAMI_CUTOFF] = opt.cut;
    } else {
      solver.solve(&model, HDUAL_VARIANT_MULTI, 8);
      model.util_reportSolverOutcome("multi");
    }
#ifdef HiGHSDEV
    if (opt.pami) model.writePivots("multi");
    if (opt.sip) model.writePivots("tasks");
#endif
  } else {
    // Serial. Based on previous solvePlainJAJH.

    double setupTime = 0;
    double presolve1Time = 0;
    double crashTime = 0;
#ifdef HiGHSDEV
    double crossoverTime = 0;
    double presolve2Time = 0;
#endif
    double solveTime = 0;
    double postsolveTime = 0;
    int solveIt = 0;
#ifdef HiGHSDEV
    int solvePh1DuIt = 0;
    int solvePh2DuIt = 0;
    int solvePrIt = 0;
#endif
    double lcSolveTime;
    HDual solver;

    const bool presolveNoScale = false;

    vector<double> colPrAct;
    vector<double> colDuAct;
    vector<double> rowPrAct;
    vector<double> rowDuAct;

    //	printf("model.intOption[INTOPT_PRINT_FLAG] = %d\n",
    //model.intOption[INTOPT_PRINT_FLAG]);
    solver.setPresolve(opt.presolveMode);
    solver.setPrice(opt.priceMode);
    solver.setEdWt(opt.edWtMode);
    solver.setCrash(opt.crashMode);
    solver.setTimeLimit(opt.timeLimit);

    model.timer.reset();
    bool with_presolve = solver.Presolve_Mode == Presolve_Mode_On;
    //  bool FourThreads = true;
    bool FourThreads = false;
    //  bool EightThreads = true;
    bool EightThreads = false;

    printf("solvePlainJAJH: with_presolve = %d\n", with_presolve);
    if (with_presolve) {
      if (solver.Crash_Mode > 0) {
        HCrash crash;
        crash.crash(&model, solver.Crash_Mode);
        crashTime += model.timer.getTime();
      }

      if (FourThreads)
        solver.solve(&model, HDUAL_VARIANT_MULTI, 4);
      else if (EightThreads)
        solver.solve(&model, HDUAL_VARIANT_MULTI, 8);
      else
        solver.solve(&model);

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
          model.getPrStatus(), model.modelName.c_str(), model.numRow,
          model.numCol, lcSolveTime, model.objective, solver.n_ph1_du_it,
          solver.n_ph2_du_it, solver.n_pr_it);
#endif

      // Possibly recover bounds after presolve (after using bounds tightened by
      // presolve)
      if (model.usingImpliedBoundsPresolve) {
        //		Recover the true bounds overwritten by the implied
        //bounds
#ifdef HiGHSDEV
        printf(
            "\nRecovering bounds after using implied bounds and resolving\n");
#endif
        if (model.problemStatus != LP_Status_OutOfTime) {
          model.copy_savedBoundsToModelBounds();

          model.timer.reset();
          solver.solve(&model);
          lcSolveTime = model.timer.getTime();
          solveTime += lcSolveTime;
          solveIt += model.numberIteration;
          model.util_reportSolverOutcome("After recover:   ");
#ifdef HiGHSDEV
          solvePh1DuIt += solver.n_ph1_du_it;
          solvePh2DuIt += solver.n_ph2_du_it;
          solvePrIt += solver.n_pr_it;
          printf(
              "\nBnchmkHsol02 After restoring bounds,hsol,%3d,%16s, %d,%d,"
              "%10.3f,%20.10e,%10d,%10d,%10d\n",
              model.getPrStatus(), model.modelName.c_str(), model.numRow,
              model.numCol, lcSolveTime, model.objective, solver.n_ph1_du_it,
              solver.n_ph2_du_it, solver.n_pr_it);
#endif
        }
      }

      if (model.problemStatus != LP_Status_OutOfTime) {
        // Perform postsolve
      }
    } else {
      setupTime += model.timer.getTime();
      if (solver.Crash_Mode > 0) {
        HCrash crash;
        //      printf("Calling crash.crash(&model,
        //      solver.Crash_Mode);\n");cout<<flush;
        crash.crash(&model, solver.Crash_Mode);
        // printf("Called  crash.crash(&model,
        // solver.Crash_Mode);\n");cout<<flush;
        crashTime += model.timer.getTime();
      }
      //		printf("model.intOption[INTOPT_PRINT_FLAG] = %d\n",
      //model.intOption[INTOPT_PRINT_FLAG]);
      model.scaleModel();
      if (FourThreads)
        solver.solve(&model, HDUAL_VARIANT_MULTI, 4);
      else if (EightThreads)
        solver.solve(&model, HDUAL_VARIANT_MULTI, 8);
      else
        solver.solve(&model);
      solveTime += model.timer.getTime();
      int problemStatus = model.getPrStatus();
      //    printf("After solve() model status is %d\n", problemStatus);
      if (problemStatus == LP_Status_Unset) {
        HCrash crash;
        crash.crash(&model, Crash_Mode_Bs);
        solver.solve(&model);
        solveTime += model.timer.getTime();
        //   int problemStatus = model.getPrStatus(); printf("After solve()
        //   model status is %d\n", problemStatus);
      }
    }
#ifdef HiGHSDEV
    double sumTime =
        setupTime + presolve1Time + crashTime + solveTime + postsolveTime;
    printf(
        "Time: setup = %10.3f; presolve = %10.3f; crash = %10.3f; solve = "
        "%10.3f; postsolve = %10.3f; sum = %10.3f; total = %10.3f\n",
        setupTime, presolve1Time, crashTime, solveTime, postsolveTime, sumTime,
        model.totalTime);
    cout << flush;
    double errTime = abs(sumTime - model.totalTime);
    if (errTime > 1e-3) printf("!! Sum-Total time error of %g\n", errTime);
#endif
    // TODO Reinstate this once solve after postsolve is performed
    //  model.util_getPrimalDualValues(colPrAct, colDuAct, rowPrAct, rowDuAct);
    //  double Ph2Objective = model.computePh2Objective(colPrAct);
    //  printf("Computed Phase 2 objective = %g\n", Ph2Objective);
    model.util_reportSolverOutcome("Final:           ");
#ifdef HiGHSDEV
    bool rpBnchmk = false;
    if (rpBnchmk) {
      int numCol = model.numCol;
      int numRow = model.numRow;
      printf(
          "\nBnchmkHsol99,hsol,%3d,%16s,Presolve %s,"
          "Crash %s,EdWt %s,Price %s,%d,%d,%10.3f,%10.3f,"
          "%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,"
          "%20.10e,%10d,%10.3f,"
          "%d\n",
          model.getPrStatus(), model.modelName.c_str(), Presolve_ArgV,
          Crash_ArgV, EdWt_ArgV, Price_ArgV, numRow, numCol, setupTime,
          presolve1Time, crashTime, crossoverTime, presolve2Time, solveTime,
          postsolveTime, model.objective, model.numberIteration,
          model.totalTime, solver.n_wg_DSE_wt);
      cout << flush;
    }
#endif
  }
  return HighsStatus::OK;
}

HighsStatus afterPostsolve(const HighsOptions& opt, HighsModelObject& model) {
  // Extract the model from what's recreated in postsolve
  // printf("\nload_fromPostsolve\n");
  model.load_fromPostsolve(pre);
  model.shiftObjectiveValue(pre->objShift);
  postsolveTime += model.timer.getTime();
  // Save the solved results
  model.totalTime += model.timer.getTime();
#ifdef HiGHSDEV
  model.util_reportModelSolution();
#endif

#ifdef HiGHSDEV
  printf("\nBefore solve after Postsolve\n");
  cout << flush;
#endif
  model.timer.reset();
  solver.solve(&model);
  lcSolveTime = model.timer.getTime();
  solveTime += lcSolveTime;
  solveIt += model.numberIteration;
  model.util_reportSolverOutcome("After postsolve: ");
#ifdef HiGHSDEV
  solvePh1DuIt += solver.n_ph1_du_it;
  solvePh2DuIt += solver.n_ph2_du_it;
  solvePrIt += solver.n_pr_it;
  printf(
      "\nBnchmkHsol03 After postsolve       ,hsol,%3d,%16s, %d,%d,"
      "%10.3f,%20.10e,%10d,%10d,%10d\n",
      model.getPrStatus(), model.modelName.c_str(), model.numRow, model.numCol,
      lcSolveTime, model.objective, solver.n_ph1_du_it, solver.n_ph2_du_it,
      solver.n_pr_it);
  cout << flush;
#endif
  return HighsStatus::OK;
}

HighsStatus solveScip(const HighsOptions& opt, HighsModelObject& highs_model) {}

// Single function to solve an lp according to options and fill
// solution in solution.
HighsStatus runSimplexSolver(const HighsOptions& opt,
                             HighsModelObject& highs_model) {
  cout << "=================================================================="
       << endl;

  // For the mo ent handle scip case separately.
  if (opt.scip) return solveScip(opt, highs_model);

  // When runSimplexSolver is called initialize an instance of HModel inside the
  // HighsModelObject. This will then be passed to HDual.
  highs_model.hmodel_.push_back(HModel());

  HModel& model = highs_model.hmodel_[0];
  const HighsLp& lp = highs_model.lp_;

  model.load_fromArrays(lp.numCol_, lp.sense_, &lp.colCost_[0],
                        &lp.colLower_[0], &lp.colUpper_[0], lp.numRow_,
                        &lp.rowLower_[0], &lp.rowUpper_[0], lp.nnz_,
                        &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);

  model.intOption[INTOPT_PRINT_FLAG] = 1;

  // Scaling: Separate from simplex.
  model.scaleModel();

  // Crash, if needed.

  HighsStatus result = solveSimplex(opt, highs_model);

  // todo uncomment line below.
  // if (result == optimal)
  // HighsSolution set values in highs_model.
  HighsSolution& solution = highs_model.solution_;
  highs_model.hmodel_[0].util_getPrimalDualValues(
      solution.colValue_, solution.colDual_, solution.rowValue_,
      solution.rowDual_);
  model.util_getBasicIndexNonbasicFlag(highs_model.basis_info_.basis_index,
                                       highs_model.basis_info_.nonbasic_flag);

  cout << "==================================================================";

  return result;
}

#endif
