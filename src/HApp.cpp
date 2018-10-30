#include "HApp.h"

using namespace std;

int solvePlain(HModel &model) {
  model.intOption[INTOPT_PRINT_FLAG] = 1;
  if (model.intOption[INTOPT_PRINT_FLAG]) model.util_reportModelBrief();
#ifdef HiGHSDEV
    //  cout << "\n Using solvePlain() - Calling model.scaleModel()\n" <<
    //  endl;
#endif
  model.scaleModel();
  HDual solver;
#ifdef HiGHSDEV
  //  cout << "\n Using solvePlain() - Calling solver.solve(&model)\n" <<
  //  endl;
#endif
  solver.solve(&model);
  model.util_reportSolverOutcome("Solve plain");
#ifdef HiGHSDEV
  // model.util_reportModelStatus();
  // model.util_reportModelBrief();
  // model.util_reportModelDense();
  // Possibly analyse the degeneracy of the primal and dual activities
  // model.util_anPrDuDgn();
  // model.util_reportModelSolution();
#endif
  return 0;
}

int solvePlainAPI(HModel &model) {
  model.intOption[INTOPT_PRINT_FLAG] = 1;
  int RpRtVec = 0;
  printf("\nUsing SolvePlainAPI\n\n");
  //  model.intOption[INTOPT_PRINT_FLAG] = 1;

  // Use the arrays read from an MPS file to test the routine to
  // solve a model passed by arrays. First copy the data.
  int XnumCol = model.numCol;
  int XnumRow = model.numRow;
  int XnumNz = model.Astart[model.numCol];
  int XobjSense = model.objSense;
  int XobjOffset = model.objOffset;
  double *XcolCost;
  double *XcolLower;
  double *XcolUpper;
  double *XrowLower;
  double *XrowUpper;
  int *XAstart;
  int *XAindex;
  double *XAvalue;

  XcolCost = (double *)malloc(sizeof(double) * XnumCol);
  XcolLower = (double *)malloc(sizeof(double) * XnumCol);
  XcolUpper = (double *)malloc(sizeof(double) * XnumCol);
  XrowLower = (double *)malloc(sizeof(double) * XnumRow);
  XrowUpper = (double *)malloc(sizeof(double) * XnumRow);
  XAstart = (int *)malloc(sizeof(int) * (XnumCol + 1));
  XAindex = (int *)malloc(sizeof(int) * XnumNz);
  XAvalue = (double *)malloc(sizeof(double) * XnumNz);

  memcpy(XcolCost, &(model.colCost[0]), sizeof(double) * model.numCol);
  memcpy(XcolLower, &(model.colLower[0]), sizeof(double) * model.numCol);
  memcpy(XcolUpper, &(model.colUpper[0]), sizeof(double) * model.numCol);
  memcpy(XrowLower, &(model.rowLower[0]), sizeof(double) * model.numRow);
  memcpy(XrowUpper, &(model.rowUpper[0]), sizeof(double) * model.numRow);
  memcpy(XAstart, &(model.Astart[0]), sizeof(int) * (XnumCol + 1));
  memcpy(XAindex, &(model.Aindex[0]), sizeof(int) * XnumNz);
  memcpy(XAvalue, &(model.Avalue[0]), sizeof(double) * XnumNz);

  model.clearModel();

  int probStatus, basisStatus;
  double *colPrimalValues;
  double *colDualValues;
  double *rowPrimalValues;
  double *rowDualValues;
  int *basicVariables;

  colPrimalValues = (double *)malloc(sizeof(double) * XnumCol);
  rowPrimalValues = (double *)malloc(sizeof(double) * XnumRow);
  colDualValues = (double *)malloc(sizeof(double) * XnumCol);
  rowDualValues = (double *)malloc(sizeof(double) * XnumRow);
  basicVariables = (int *)malloc(sizeof(int) * XnumRow);

  probStatus = 0;
  basisStatus = HiGHS_basisStatus_no;
  solve_fromArrays(&probStatus, &basisStatus, XnumCol, XnumRow, XnumNz,
                   XobjSense, XobjOffset, XcolCost, XcolLower, XcolUpper,
                   XrowLower, XrowUpper, XAstart, XAindex, XAvalue,
                   colPrimalValues, colDualValues, rowPrimalValues,
                   rowDualValues, basicVariables);
  printf("GOT back from solve_fromArrays: probStatus = %d\n", probStatus);
  cout << flush;
  if (RpRtVec) {
    for (int col = 0; col < XnumCol; col++) {
      printf("Col %3d: Pr, Du = %11.4g, %11.4g\n", col, colPrimalValues[col],
             colDualValues[col]);
    }
    for (int row = 0; row < XnumRow; row++) {
      printf("Row %3d: Bc = %3d; Pr, Du = %11.4g, %11.4g\n", row,
             basicVariables[row], rowPrimalValues[row], rowDualValues[row]);
    }
  }
  basisStatus = HiGHS_basisStatus_yes;
  solve_fromArrays(&probStatus, &basisStatus, XnumCol, XnumRow, XnumNz,
                   XobjSense, XobjOffset, XcolCost, XcolLower, XcolUpper,
                   XrowLower, XrowUpper, XAstart, XAindex, XAvalue,
                   colPrimalValues, colDualValues, rowPrimalValues,
                   rowDualValues, basicVariables);
  printf("GOT back from solve_fromArrays: probStatus = %d\n", probStatus);
  cout << flush;
  if (RpRtVec) {
    for (int col = 0; col < XnumCol; col++) {
      printf("Col %3d: Pr, Du = %11.4g, %11.4g\n", col, colPrimalValues[col],
             colDualValues[col]);
    }
    for (int row = 0; row < XnumRow; row++) {
      printf("Row %3d: Bc = %3d; Pr, Du = %11.4g, %11.4g\n", row,
             basicVariables[row], rowPrimalValues[row], rowDualValues[row]);
    }
  }

  /*
  //TEMP CODE TO GENERATE BOXED PROBLEMS
  for (int col = 0; col < XnumCol; col++) {
    //    printf("Col %3d has (%11.4g, %11.4g, %11.4g)", col, XcolLower[col],
  colPrimalValues[col], XcolUpper[col]); if (-XcolLower[col] >=
  HSOL_CONST_INF) { XcolLower[col] = colPrimalValues[col] - 0.1*fmax(1,
  fabs(colPrimalValues[col]));
    }
    if (XcolUpper[col] >= HSOL_CONST_INF) {
      XcolUpper[col] = colPrimalValues[col] + 0.1*fmax(1,
  fabs(colPrimalValues[col]));
    }
    //    printf(": now (%11.4g, %11.4g)\n", XcolLower[col], XcolUpper[col]);
  }
  const char *fileName = "WrMl.mps";
  //  int* integerColumn; integerColumn = (int *) malloc(sizeof(int)*XnumCol);
  for (int col=0; col< XnumCol; col++) {integerColumn[col]=0;} vector<int>
  integerColumn; integerColumn.assign(XnumCol,0); writeMPS(fileName, XnumRow,
  XnumCol, XobjSense, XobjOffset, XAstart, XAindex, XAvalue, XcolCost,
  XcolLower, XcolUpper, XrowLower, XrowUpper, integerColumn);
  //TEMP CODE TO GENERATE BOXED PROBLEMS
*/
  return 0;
}

// Ivet
int solvePlainWithPresolve(HModel &model) {
  double time1;

  double obj1 = presolve(model, time1);
  (void)obj1;

  // to test singularity of basis matrix after postsolve
  /*
        HModel model2;

        model2.load_fromMPS(filename);
        model2.scaleModel();

        HDual solver2;
        solver2.solve(&model2);
        model2.util_reportSolverOutcome("Presolve 1");
        double obj2 = model2.util_getObjectiveValue();

        //testing
        int * x = model2.getNonbasicFlag();
        int len = model2.getNumTot();
        vector<int> v(len);
        for (int i=0;i<len;i++)
                v[i] = (x[i]);

        int * xx = model2.getBaseIndex();
        len = model2.getNumRow();
        vector<int> vv(len);
        for (int i=0;i<len;i++)
                vv[i] = (xx[i]);

        HModel mod;
        mod.load_fromMPS(filename);
        double time1;
        HPresolve * pre = new HPresolve();
        mod.copy_fromHModelToHPresolve(pre);
        pre->presolve();
        mod.load_fromPresolve(pre);

        HDual solver;
        mod.scaleModel();
    solver.solve(&mod);
    mod.timer.reset();

    pre->setProblemStatus(mod.getPrStatus());
    mod.util_getPrimalDualValues(pre->colValue, pre->colDual, pre->rowValue,
    pre->rowDual); mod.util_getBasicIndexNonbasicFlag(pre->basicIndex,
    pre->nonbasicFlag); pre->setNBFfullproblem(v, vv);  //testing
    pre->postsolve();
    mod.shiftObjectiveValue(pre->objShift);

    mod.util_reportSolverOutcome("Postsolve 1");

    time1 = mod.totalTime ;
        double obj1 = mod.util_getObjectiveValue();
        delete pre;


        if (abs(obj1-obj2)>0.000001)
                cout<<"OBJECTIVE FAIL: diff = "<<obj1-obj2<<endl;
        */

  /*	// to compare data presolve receives and returns after postsolve
        HModel model;
        model.load_fromMPS(filename);

        HinOut test("","");
        test.getData(model);


        double time1;
        double obj1 = presolve(model, time1);

        test.readDataPostsolve(model);
        test.compareData(2);
*/
  return 0;
}

// Julian
int solveSCIP(HModel &model) {
  printf("Called solveSCIP\n");
  cout << flush;
  //  model.util_reportModel();

  // Extract columns numCol-3..numCol-1
  int FmCol = model.numCol - 3;
  int ToCol = model.numCol - 1;
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
  int FmRow = model.numRow - 3;
  int ToRow = model.numRow - 1;
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
  ToRow = model.numRow - 1;
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
  ToCol = model.numCol - 1;
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

  model.numTot = model.numCol + model.numRow;
  model.scaleModel();
  HDual solver;
  solver.solve(&model);
  model.util_reportModelSolution();
  model.util_reportSolverOutcome("SCIP 1");

  vector<double> colPrimal(model.numCol);
  vector<double> colDual(model.numCol);
  vector<double> colLower(model.numCol);
  vector<double> colUpper(model.numCol);
  vector<double> rowPrimal(model.numRow);
  vector<double> rowDual(model.numRow);
  vector<double> rowLower(model.numRow);
  vector<double> rowUpper(model.numRow);
  model.util_getPrimalDualValues(colPrimal, colDual, rowPrimal, rowDual);
  model.util_getColBounds(0, model.numCol - 1, &colLower[0], &colUpper[0]);
  model.util_getRowBounds(0, model.numRow - 1, &rowLower[0], &rowUpper[0]);

  double og_colLower;
  double og_colUpper;
  int colBoundIndex;
  double nw_colLower;
  double nw_colUpper;

  int num_resolve = 0;
  for (int col = 0; col < model.numCol; col++) {
    model.util_getColBounds(col, col, &og_colLower, &og_colUpper);
    printf("\nColumn %2d has primal value %11g and bounds [%11g, %11g]", col,
           colPrimal[col], og_colLower, og_colUpper);
    if (model.nonbasicFlag[col]) {
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
      if (model.hsol_isInfinity(og_colUpper))
        nw_colLower = colPrimal[col] + 1;
      else
        nw_colLower = og_colUpper;
      nw_colUpper = og_colUpper;
      printf("Calling model.util_chgColBounds(1, %d, %g, %g)\n", colBoundIndex,
             nw_colLower, nw_colUpper);
      model.util_chgColBoundsSet(1, &colBoundIndex, &nw_colLower, &nw_colUpper);
      printf("Calling model.scaleModel()\n");
      model.scaleModel();
      //      printf("Calling solver.solve(&model)\n");
      solver.solve(&model);
      //      printf("Called solver.solve(&model)\n");
      model.util_reportSolverOutcome("SCIP 2");
      // Was &nw_colLower, &nw_colUpper); and might be more interesting for
      // avgas
      model.util_chgColBoundsSet(1, &colBoundIndex, &og_colLower, &og_colUpper);
      if (num_resolve >= 10) break;
    }
  }
  printf("Returning from solveSCIP\n");
  cout << flush;
  return 0;
}

int solvePlainJAJH(HModel &model, const char *Price_ArgV, const char *EdWt_ArgV,
                   const char *Crash_ArgV, const char *Presolve_ArgV,
                   double TimeLimit_ArgV) {
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
  model.intOption[INTOPT_PRINT_FLAG] = 1;
  HDual solver;

  const bool presolveNoScale = false;

  vector<double> colPrAct;
  vector<double> colDuAct;
  vector<double> rowPrAct;
  vector<double> rowDuAct;

  //	printf("model.intOption[INTOPT_PRINT_FLAG] = %d\n",
  // model.intOption[INTOPT_PRINT_FLAG]);
  solver.setPresolve(Presolve_ArgV);
  solver.setPrice(Price_ArgV);
  solver.setEdWt(EdWt_ArgV);
  solver.setCrash(Crash_ArgV);
  solver.setTimeLimit(TimeLimit_ArgV);

  model.timer.reset();
  bool with_presolve = solver.Presolve_Mode == Presolve_Mode_On;
  //  bool FourThreads = true;
  bool FourThreads = false;
  //  bool EightThreads = true;
  bool EightThreads = false;

  printf("solvePlainJAJH: with_presolve = %d\n", with_presolve);
  if (with_presolve) {
    // Check size
    if (model.numRow == 0) return 1;
    HPresolve *pre = new HPresolve();
    model.copy_fromHModelToHPresolve(pre);
    setupTime += model.timer.getTime();
    model.timer.reset();
    pre->presolve();
    //		For consistency, the following should be done within
    // presolve
    model.totalTime += model.timer.getTime();
    pre->reportTimes();
    model.load_fromPresolve(pre);
    presolve1Time += model.timer.getTime();
    setupTime += model.timer.getTime();

    if (solver.Crash_Mode > 0) {
      HCrash crash;
      crash.crash(&model, solver.Crash_Mode);
      crashTime += model.timer.getTime();
    }

    //    printf("model.intOption[INTOPT_PRINT_FLAG] = %d\n",
    //    model.intOption[INTOPT_PRINT_FLAG]);
    if (presolveNoScale)
      printf(
          "*****************************\n* !!Not currently scaling!! "
          "*\n*****************************\n");
    else {
      model.scaleModel();
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
    model.util_reportSolverOutcome("After presolve:  ");
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
    // bounds
#ifdef HiGHSDEV
      printf("\nRecovering bounds after using implied bounds and resolving\n");
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
#ifdef HiGHSDEV
      printf("\nPostsolving\n");
#endif
      model.timer.reset();
      // Copy model's problem status into presolve's problem status
      pre->setProblemStatus(model.getPrStatus());
      // Extract solution into primal and dual, row and column
      // arrays. Undos scaling since presolve works on the unscaled
      // model
      model.util_getPrimalDualValues(pre->colValue, pre->colDual, pre->rowValue,
                                     pre->rowDual);
      model.util_getBasicIndexNonbasicFlag(pre->basicIndex, pre->nonbasicFlag);

      // Perform postsolve
      pre->postsolve();
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
          model.getPrStatus(), model.modelName.c_str(), model.numRow,
          model.numCol, lcSolveTime, model.objective, solver.n_ph1_du_it,
          solver.n_ph2_du_it, solver.n_pr_it);
      cout << flush;
#endif
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
    // model.intOption[INTOPT_PRINT_FLAG]);
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
        model.getPrStatus(), model.modelName.c_str(), Presolve_ArgV, Crash_ArgV,
        EdWt_ArgV, Price_ArgV, numRow, numCol, setupTime, presolve1Time,
        crashTime, crossoverTime, presolve2Time, solveTime, postsolveTime,
        model.objective, model.numberIteration, model.totalTime,
        solver.n_wg_DSE_wt);
    cout << flush;
  }
#endif
  return 0;
}

Status presolve(const HighsLp &lp, HighsLp &reduced_lp) {
  return Status::NotImplemented;
}

double presolve(HModel &mod, double &time) {
  cout << "------\n";

  HPresolve *pre = new HPresolve();
  mod.copy_fromHModelToHPresolve(pre);
  int status = pre->presolve();
  if (status == HPresolve::Unset) {
    // pre->reportTimes();
    mod.load_fromPresolve(pre);

    HDual solver;
    mod.scaleModel();
    solver.solve(&mod);
    pre->setProblemStatus(mod.getPrStatus());
    mod.util_getPrimalDualValues(pre->colValue, pre->colDual, pre->rowValue,
                                 pre->rowDual);
    mod.util_getBasicIndexNonbasicFlag(pre->basicIndex, pre->nonbasicFlag);
    pre->postsolve();
    mod.load_fromPostsolve(pre);
    solver.solve(&mod);
    mod.util_reportSolverOutcome("Postsolve");
    time = mod.totalTime;
  } else if (status == HPresolve::Empty) {
    pre->postsolve();
    mod.load_fromPostsolve(pre);
    HDual solver;

    solver.solve(&mod);
    mod.util_reportSolverOutcome("Postsolve");
    time = mod.totalTime;
  } else {
    if (status == HPresolve::Infeasible)
      mod.problemStatus = LP_Status_Infeasible;
    else if (status == HPresolve::Unbounded)
      mod.problemStatus = LP_Status_Unbounded;
    else {
      std::cout << "Unknown, status=" << status << std::endl;
      mod.problemStatus = LP_Status_Failed;
      delete pre;
      return 0;
    }

    mod.util_reportSolverOutcome("Presolve");
  }

  delete pre;
  return mod.util_getObjectiveValue();
}

int solvePlainExperiments(const char *filename) {
  bool exp = true;
  ofstream myfile;
  if (exp) {
    std::string crName(filename);

    std::string sub2 = crName;
    if (sub2[0] == '.' && sub2[1] == '.')
      sub2 = crName.substr(10, crName.size());
    myfile.open("../experiments/out", ios::app);
    myfile << sub2;
    myfile << " &  ";
    myfile.close();
    myfile.open("../experiments/t1", ios::app);
    myfile << sub2;
    myfile << " &  ";
    myfile.close();
    //		myfile.open("../experiments/t2", ios::app);
    //		myfile << sub2;
    //		myfile << " &  ";
    //		myfile.close();
    myfile.open("../experiments/t3", ios::app);
    myfile << sub2;
    myfile << " &  ";
    myfile.close();
  }

  HModel model;
  int RtCd;
  RtCd = model.load_fromMPS(filename);
  if (RtCd) return RtCd;
  // Check size
  if (model.numRow == 0) return 1;

  double time1;
  double obj1 = presolve(model, time1);

  cout << "----------\n";

  HModel model2;
  RtCd = model2.load_fromMPS(filename);
  if (RtCd) return RtCd;
  model2.scaleModel();

  HDual solver2;
  solver2.solve(&model2);
  model2.util_reportSolverOutcome("SolvePlainExperiments");
  double obj2 = model2.util_getObjectiveValue();

  // get primal column solution
  vector<double> colValue2, colDual2, rowValue2, rowDual2;
  model.util_getPrimalDualValues(colValue2, colDual2, rowValue2, rowDual2);

  if (exp) {
    ofstream myfile;
    myfile.open("../experiments/out", ios::app);
    if (abs(obj1 - obj2) <= 0.000001)
      myfile << " obj pass" << endl;
    else
      myfile << " obj fail" << endl;
    myfile.close();
  }

  return 0;
}

int solveTasks(HModel &model) {
  model.intOption[INTOPT_PRINT_FLAG] = 1;
  model.intOption[INTOPT_PERMUTE_FLAG] = 1;

  model.scaleModel();
  HDual solver;
  solver.solve(&model, HDUAL_VARIANT_TASKS, 8);

  model.util_reportSolverOutcome("Solve tasks");
#ifdef HiGHSDEV
  model.writePivots("tasks");
#endif
  return 0;
}

int solveMulti(HModel &model, const char *partitionfile) {
  model.intOption[INTOPT_PRINT_FLAG] = 1;
  model.intOption[INTOPT_PERMUTE_FLAG] = 1;
  if (partitionfile) {
    model.strOption[STROPT_PARTITION_FILE] = partitionfile;
  }

  model.scaleModel();
  HDual solver;
  //    solver.solve(&model, HDUAL_VARIANT_MULTI, 1);
  //    solver.solve(&model, HDUAL_VARIANT_MULTI, 2);
  //    solver.solve(&model, HDUAL_VARIANT_MULTI, 4);
  solver.solve(&model, HDUAL_VARIANT_MULTI, 8);

  model.util_reportSolverOutcome("Solve multi");
#ifdef HiGHSDEV
  model.writePivots("multi");
#endif
  return 0;
}

#ifdef EXT_PRESOLVE

void copyMatrix(const Problem<double> &problem, vector<int> &Astart,
                vector<int> &Aend, vector<int> &Aindex,
                vector<double> &Avalue) {
  int numCol = problem.getNCols();
  int nnz = problem.getConstraintMatrix().getNnz();

  assert((unsigned int)numCol + 1 == Aend.size());
  assert((unsigned int)numCol + 1 == Astart.size());

  vector<pair<int, size_t>> vp;
  vp.reserve(numCol);

  for (int i = 0; i != numCol; ++i) {
    vp.push_back(make_pair(Astart.at(i), i));
  }

  // Sorting will put lower values ahead of larger ones,
  // resolving ties using the original index
  sort(vp.begin(), vp.end());

  vector<int> Aendtmp;
  Aendtmp = Aend;
  const int *Aindex_ = problem.getConstraintMatrix().getTransposeRowIndices();
  const double *Avalue_ = problem.getConstraintMatrix().getTransposeValues();

  int iPut = 0;
  for (size_t i = 0; i != vp.size(); ++i) {
    int col = vp.at(i).second;
    int k = vp.at(i).first;
    Astart.at(col) = iPut;
    while (k < Aendtmp.at(col)) {
      Avalue[iPut] = Avalue_[k];
      Aindex[iPut] = Aindex_[k];
      iPut++;
      k++;
    }
    Aend.at(col) = iPut;
  }

  assert(iPut == nnz);

  return;
}

int solveExternalPresolve(const char *fileName) {
  HModel model;
  int RtCd = model.load_fromMPS(fileName);
  if (RtCd) return RtCd;

  // Now we got a loaded model that we will pass to external presolve
  // set up data
  SparseStorage<double> matrix_transpose(&model.Avalue[0], &model.Astart[0],
                                         &model.Aindex[0], model.numCol,
                                         model.numRow, model.Avalue.size());

  // code below makes a row wise copy and passes that
  // HPresolve *pre = new HPresolve();
  // model.copy_fromHModelToHPresolve(pre);
  // pre->makeARCopy();
  // SparseStorage<double> matrix(&pre->ARvalue[0], &pre->ARstart[0],
  // &pre->ARindex[0],
  //                            model.numRow, model.numCol,
  //                            model.Avalue.size());
  // problem.setConstraintMatrix(matrix, model.rowLower, model.rowUpper);

  Problem<double> problem;
  problem.setObjective(model.colCost);
  problem.setName(string(fileName));
  problem.setConstraintMatrix(matrix_transpose, model.rowLower, model.rowUpper,
                              true);
  problem.setVariableDomainsLP(model.colLower, model.colUpper);
  problem.fixInfiniteBounds(HSOL_CONST_INF);

  // presolve
  Presolve<double> presolve;
  // presolve.addPresolveMethod(...);
  presolve.apply(problem);

  // Load presolved problem in solver

  // Update old HModel and set up solver to solve
  vector<int> Astart = problem.getConstraintMatrix().getTransposeColStart();
  vector<int> Aend = problem.getConstraintMatrix().getTransposeColEnd();

  int *Aindex_p = NULL;
  double *Avalue_p = NULL;

  vector<int> Aindex;
  vector<double> Avalue;

  // check if matrix copy is necessary
  int numCol = problem.getNCols();
  bool isNeeded = false;
  for (int i = 0; i < numCol; ++i) {
    if (Astart[i + 1] != Aend[i]) {
      isNeeded = true;
      break;
    }
  }

  // if we need matrix copy
  if (isNeeded) {
    int nnz = problem.getConstraintMatrix().getNnz();
    Aindex.resize(nnz);
    Avalue.resize(nnz);
    copyMatrix(problem, Astart, Aend, Aindex, Avalue);
    Astart.at(numCol) = nnz;
  }
  // if not needed do not make matrix copy
  else {
    Aindex_p = (int *)problem.getConstraintMatrix().getTransposeRowIndices();
    Avalue_p = (double *)problem.getConstraintMatrix().getTransposeValues();
  }

  vector<double> colLower = problem.getLowerBounds();
  vector<double> colUpper = problem.getUpperBounds();
  vector<double> rowLower = problem.getConstraintMatrix().getLeftHandSides();
  vector<double> rowUpper = problem.getConstraintMatrix().getRightHandSides();

  int nCols = problem.getNCols();
  int nRows = problem.getNRows();

  for (int i = 0; i < nCols; i++) {
    if (colLower.at(i) <= -HSOL_CONST_INF) colLower.at(i) = -HSOL_CONST_INF;
    if (colUpper.at(i) >= HSOL_CONST_INF) colUpper.at(i) = HSOL_CONST_INF;
  }

  for (int i = 0; i < nRows; i++) {
    if (rowLower.at(i) <= -HSOL_CONST_INF) rowLower.at(i) = -HSOL_CONST_INF;
    if (rowUpper.at(i) >= HSOL_CONST_INF) rowUpper.at(i) = HSOL_CONST_INF;
  }

  if (!Aindex_p) Aindex_p = &Aindex[0];
  if (!Avalue_p) Avalue_p = &Avalue[0];

  model.load_fromArrays(nCols, 1, &(problem.getObjective().coefficients[0]),
                        &colLower[0], &colUpper[0], nRows, &rowLower[0],
                        &rowUpper[0], problem.getConstraintMatrix().getNnz(),
                        &Astart[0], Aindex_p, Avalue_p);

  model.scaleModel();

  // solve
  HDual solver;
  solver.solve(&model);
  double obj = model.util_getObjectiveValue();

  // pass reduced solutions back to external presolve class for postsolve
  vector<double> colValue, colDual, rowValue, rowDual;
  model.util_getPrimalDualValues(colValue, colDual, rowValue, rowDual);
  problem.setPrimalValues(colValue);
  problem.setDualValues(colDual);
  problem.setRowValues(rowValue);
  problem.setRowDuals(rowDual);

  // postsolve
  problem.postsolve();

  // get solution from postsolve
  colValue = problem.getPrimalValues();
  colDual = problem.getDualValues();
  rowValue = problem.getRowValues();
  rowDual = problem.getRowDuals();

  double objective = problem.getOptimalObjective();

  // report solution
  cout << "Optimal objecive value after postsolve: " << objective << endl;

  // check that obj value is the same by solving the original problemV

  HModel model2;
  RtCd = model2.load_fromMPS(fileName);
  if (RtCd) return RtCd;
  model2.scaleModel();

  // HPresolve *pre3 = new HPresolve();
  // model2.copy_fromHModelToHPresolve(pre3);
  // pre3->initializeVectors();
  // pre3->print(0);

  // check
  HDual solver2;
  solver2.solve(&model2);
  double obj2 = model2.util_getObjectiveValue();

  if (abs(obj2 - obj) > 0.00001)
    cout << "OBJECTIVES DIFFER" << endl;
  else
    cout << "Objectives match." << endl;
  return 0;
}
#endif
