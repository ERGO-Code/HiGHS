/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDual.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HDual.h"
#include "HConst.h"
#include "HCrash.h"
#include "HMatrix.h"
#include "HPrimal.h"
#include "HTimer.h"

//#include "HighsSetup.h"
#include "HighsLp.h"
#include "HighsIO.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>
using namespace std;

void HDual::solve(HModel *ptr_model, int variant, int num_threads) {
  assert(ptr_model != NULL);
  dual_variant = variant;
  model = ptr_model;
  // Setup two work buffers in model required for solve()
  model->buffer.setup(model->numRow);
  model->bufferLong.setup(model->numCol);
  // Setup aspects of the model data which are needed for solve() but better
  // left until now for efficiency reasons.
  model->setup_for_solve();
  model->problemStatus = LP_Status_Unset;
  model->numberIteration = 0;
#ifdef HiGHSDEV
  // Initialise numbers and times of rebuilds and inverts.
  totalRebuilds = 0;
  totalRebuildTime = 0;
  model->totalInverts = 0;
  model->totalInvertTime = 0;
#endif
  // Cannot solve box-constrained LPs TODO: Fix this
  if (model->numRow == 0) return;

  model->timer.reset();

  n_ph1_du_it = 0;
  n_ph2_du_it = 0;
  n_pr_it = 0;
  // Set SolveBailout to be true if control is to be returned immediately to
  // calling function
  SolveBailout = false;

  HighsPrintMessage(HighsMessageType::INFO, "Using HighsPrintMessage to report TimeLimitValue   on entry to HDual::solve() as %12g\n", TimeLimitValue);
  //  HighsPrintMessage(HighsMessageType::INFO, "Using HighsPrintMessage to report option.timeLimit on entry to HDual::solve() as %12g\n", option.timeLimit);

  if (TimeLimitValue == 0) {
    TimeLimitValue = 1000000.0;
#ifdef HiGHSDEV
    HighsPrintMessage(HighsMessageType::WARNING, "Changing TimeLimitValue from 0 to %g\n", TimeLimitValue);
#endif
  }

  // Initialise working environment
  // Does LOTS, including initialisation of edge weights. Should only
  // be called if model dimension changes
  init(num_threads);

  model->initCost(1);
  if (!model->mlFg_haveFreshInvert) {
    int rankDeficiency = model->computeFactor();
    if (rankDeficiency) {
      throw runtime_error("Dual initialise: singular-basis-matrix");
    }
#ifdef HiGHSDEV
    double bsCond = an_bs_cond(model);
    HighsPrintMessage(HighsMessageType::INFO, "Initial basis condition estimate of %11.4g is", bsCond);
    if (bsCond > 1e12) {
      HighsPrintMessage(HighsMessageType::INFO, " excessive\n");
      return;
    } else {
      HighsPrintMessage(HighsMessageType::INFO, " OK\n");
    }
#endif
  }
  // Consider initialising edge weights
  //
  // NB workEdWt is assigned and initialised to 1s in
  // dualRHS.setup(model) so that CHUZR is well defined, even for
  // Dantzig pricing
  //
#ifdef HiGHSDEV
  //  printf("model->mlFg_haveEdWt 2 = %d; EdWt_Mode = %d; EdWt_Mode_DSE =
  //  %d\n",
  //	 model->mlFg_haveEdWt, EdWt_Mode, EdWt_Mode_DSE);cout<<flush;
  //  printf("Edge weights known? %d\n", !model->mlFg_haveEdWt);cout<<flush;
#endif
  if (!model->mlFg_haveEdWt) {
    // Edge weights are not known
    // Set up edge weights according to EdWt_Mode and iz_DSE_wt
    if (EdWt_Mode == EdWt_Mode_Dvx) {
      // Using dual Devex edge weights
      // Zero the number of Devex frameworks used and set up the first one
      n_dvx_fwk = 0;
      dvx_ix.assign(numTot, 0);
      iz_dvx_fwk();
    } else if (EdWt_Mode == EdWt_Mode_DSE) {
      // Using dual steepest edge (DSE) weights
      int numBasicStructurals = numRow - model->numBasicLogicals;
#ifdef HiGHSDEV
      n_wg_DSE_wt = 0;
      printf(
          "If (0<numBasicStructurals = %d) && %d = iz_DSE_wt: Compute exact "
          "DSE weights\n",
          numBasicStructurals, iz_DSE_wt);
#endif
      if (numBasicStructurals > 0 && iz_DSE_wt) {
        // Basis is not logical and DSE weights are to be initialised
#ifdef HiGHSDEV
        printf("Compute exact DSE weights\n");  // int RpI = 1;
        double IzDseEdWtTT = model->timer.getTime();
#endif
        for (int i = 0; i < numRow; i++) {
#ifdef HiGHSDEV
          //	  if (i==RpI) {printf("Computing exact DSE weight %d\n", i); RpI
          //= RpI*2;}
#endif
          row_ep.clear();
          row_ep.count = 1;
          row_ep.index[0] = i;
          row_ep.array[i] = 1;
          row_ep.packFlag = false;
          factor->btran(row_ep, row_epDensity);
          dualRHS.workEdWt[i] = row_ep.norm2();
          double lc_OpRsDensity = (double)row_ep.count / numRow;
          uOpRsDensityRec(lc_OpRsDensity, row_epDensity);
        }
#ifdef HiGHSDEV
        IzDseEdWtTT = model->timer.getTime() - IzDseEdWtTT;
        printf("Computed %d initial DSE weights in %gs\n", numRow, IzDseEdWtTT);
        if (model->intOption[INTOPT_PRINT_FLAG])
          printf(
              "solve:: %d basic structurals: computed %d initial DSE weights "
              "in %gs, %d, %d, %g\n",
              numBasicStructurals, numRow, IzDseEdWtTT, numBasicStructurals,
              numRow, IzDseEdWtTT);
#endif
      }
#ifdef HiGHSDEV
      else {
        if (model->intOption[INTOPT_PRINT_FLAG])
          printf(
              "solve:: %d basic structurals: starting from B=I so unit initial "
              "DSE weights\n",
              numBasicStructurals);
      }
#endif
    }
    // Indicate that edge weights are known
    model->mlFg_haveEdWt = 1;
  }

#ifdef HiGHSDEV
  //  printf("model->mlFg_haveEdWt 3 = %d\n", model->mlFg_haveEdWt);cout<<flush;
  bool rp_bs_cond = false;
  if (rp_bs_cond) {
    double bs_cond = an_bs_cond(ptr_model);
    printf("Initial basis condition estimate is %g\n", bs_cond);
  }
#endif

  model->computeDual();

  model->computeDualInfeasInDual(&dualInfeasCount);
  solvePhase = dualInfeasCount > 0 ? 1 : 2;

  // Find largest dual. No longer adjust the dual tolerance accordingly
  double largeDual = 0;
  for (int i = 0; i < numTot; i++) {
    if (model->getNonbasicFlag()[i]) {
      double myDual = fabs(workDual[i] * jMove[i]);
      if (largeDual < myDual) largeDual = myDual;
    }
  }
#ifdef HiGHSDEV
  //  printf("Solve: Large dual = %g\n", largeDual);cout<<flush;
#endif

  // Check that the model is OK to solve:
  //
  // Level 0 just checks the flags
  //
  // Level 1 also checks that the basis is OK and that the necessary
  // data in work* is populated.
  //
  // Level 2 (will) checks things like the nonbasic duals and basic
  // primal values
  //
  // Level 3 (will) checks expensive things like the INVERT and
  // steepeest edge weights
  //
#ifdef HiGHSDEV
  //  if ((solvePhase != 1) && (solvePhase != 2)) {printf("In solve():
  //  solvePhase = %d\n", solvePhase);cout<<flush;}
#endif
  bool ok = model->OKtoSolve(1, solvePhase);
  if (!ok) {
    printf("NOT OK TO SOLVE???\n");
    cout << flush;
  }
  //  assert(ok);

#ifdef HiGHSDEV
  //  Analyse the initial values of primal and dual variables
  //  an_iz_vr_v();
#endif

  // The major solving loop

  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiGHSDEV
  iterateIzAn();

  while (solvePhase) {
#ifdef HiGHSDEV
    int it0 = model->numberIteration;
    // printf("HDual::solve Phase %d: Iteration %d; totalTime = %g; timer.getTime = %g\n",
    // solvePhase, model->numberIteration, model->totalTime, model->timer.getTime());cout<<flush;
#endif
    switch (solvePhase) {
      case 1:
        solve_phase1();
#ifdef HiGHSDEV
        n_ph1_du_it += (model->numberIteration - it0);
#endif
        break;
      case 2:
        solve_phase2();
#ifdef HiGHSDEV
        n_ph2_du_it += (model->numberIteration - it0);
#endif
        break;
      case 4:
        break;
      default:
        solvePhase = 0;
        break;
    }
    // Jump for primal
    if (solvePhase == 4) break;
    // Possibly bail out
    if (SolveBailout) break;
  }

#ifdef HiGHSDEV
  if (AnIterLg) iterateRpAn();
  // Report the ticks before primal
  if (dual_variant == HDUAL_VARIANT_PLAIN) {
    int reportList[] = {
        HTICK_INVERT,        HTICK_PERM_WT,        HTICK_COMPUTE_DUAL,
        HTICK_CORRECT_DUAL,  HTICK_COMPUTE_PRIMAL, HTICK_COLLECT_PR_IFS,
        HTICK_COMPUTE_DUOBJ, HTICK_REPORT_INVERT,  HTICK_CHUZR1,
        HTICK_BTRAN,         HTICK_PRICE,          HTICK_CHUZC0,
        HTICK_CHUZC1,        HTICK_CHUZC2,         HTICK_CHUZC3,
        HTICK_CHUZC4,        HTICK_DEVEX_WT,       HTICK_FTRAN,
        HTICK_FTRAN_BFRT,    HTICK_FTRAN_DSE,      HTICK_UPDATE_DUAL,
        HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,  HTICK_DEVEX_IZ,
        HTICK_UPDATE_PIVOTS, HTICK_UPDATE_FACTOR,  HTICK_UPDATE_MATRIX};
    int reportCount = sizeof(reportList) / sizeof(int);
    model->timer.report(reportCount, reportList, 1.0);
    bool rpIterate = false;
    if (rpIterate) {
      int reportList[] = {HTICK_ITERATE};
      int reportCount = sizeof(reportList) / sizeof(int);
      model->timer.report(reportCount, reportList, 1.0);
    }
    if (rpIterate) {
      int reportList[] = {
          HTICK_ITERATE_REBUILD, HTICK_ITERATE_CHUZR,    HTICK_ITERATE_CHUZC,
          HTICK_ITERATE_FTRAN,   HTICK_ITERATE_VERIFY,   HTICK_ITERATE_DUAL,
          HTICK_ITERATE_PRIMAL,  HTICK_ITERATE_DEVEX_IZ, HTICK_ITERATE_PIVOTS};
      int reportCount = sizeof(reportList) / sizeof(int);
      model->timer.report(reportCount, reportList, 1.0);
    }
  }

  if (dual_variant == HDUAL_VARIANT_TASKS) {
    int reportList[] = {
        HTICK_INVERT,        HTICK_CHUZR1,        HTICK_BTRAN,
        HTICK_PRICE,         HTICK_CHUZC1,        HTICK_CHUZC2,
        HTICK_CHUZC3,        HTICK_DEVEX_WT,      HTICK_FTRAN,
        HTICK_FTRAN_BFRT,    HTICK_FTRAN_DSE,     HTICK_UPDATE_DUAL,
        HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT, HTICK_UPDATE_FACTOR,
        HTICK_GROUP1};
    int reportCount = sizeof(reportList) / sizeof(int);
    model->timer.report(reportCount, reportList, 1.0);
  }

  if (dual_variant == HDUAL_VARIANT_MULTI) {
    int reportList[] = {
        HTICK_INVERT,        HTICK_CHUZR1,        HTICK_BTRAN,
        HTICK_PRICE,         HTICK_CHUZC1,        HTICK_CHUZC2,
        HTICK_CHUZC3,        HTICK_DEVEX_WT,      HTICK_FTRAN,
        HTICK_FTRAN_BFRT,    HTICK_FTRAN_DSE,     HTICK_UPDATE_DUAL,
        HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT, HTICK_UPDATE_FACTOR,
        HTICK_UPDATE_ROW_EP};
    int reportCount = sizeof(reportList) / sizeof(int);
    model->timer.report(reportCount, reportList, 1.0);
    printf("PAMI   %-20s    CUTOFF  %6g    PERSISTENSE  %6g\n",
           model->modelName.c_str(), model->dblOption[DBLOPT_PAMI_CUTOFF],
           model->numberIteration / (1.0 + multi_iteration));
  }
#endif

  if (model->problemStatus != LP_Status_OutOfTime) {
    // Use primal to clean up if not out of time
#ifdef HiGHSDEV
    int it0 = model->numberIteration;
#endif
    if (solvePhase == 4) {
      HPrimal hPrimal;
      hPrimal.TimeLimitValue = TimeLimitValue;
      hPrimal.solvePhase2(model);
      // Add in the count and time for any primal rebuilds
#ifdef HiGHSDEV
      totalRebuildTime += hPrimal.totalRebuildTime;
      totalRebuilds += hPrimal.totalRebuilds;
#endif
    }
#ifdef HiGHSDEV
    n_pr_it += (model->numberIteration - it0);
#endif
  }
  // Save the solved results
  model->totalTime += model->timer.getTime();

#ifdef HiGHSDEV
  if (n_ph1_du_it + n_ph2_du_it + n_pr_it != model->numberIteration) {
    printf("Iteration total error \n");
  }
  printf("Iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d\n", n_ph1_du_it,
         n_ph2_du_it, n_pr_it, model->numberIteration);
  if (EdWt_Mode == EdWt_Mode_Dvx) {
    printf("Devex: n_dvx_fwk = %d; Average n_dvx_it = %d\n", n_dvx_fwk,
           model->numberIteration / n_dvx_fwk);
  }
  if (rp_bs_cond) {
    double bs_cond = an_bs_cond(ptr_model);
    printf("Optimal basis condition estimate is %g\n", bs_cond);
  }
#endif
#ifdef HiGHSDEV
  //  if ((solvePhase != 1) && (solvePhase != 2)) {printf("In solve():
  //  solvePhase = %d\n", solvePhase);cout<<flush;}
#endif
  ok = model->OKtoSolve(1, solvePhase);
  if (!ok) {
    printf("NOT OK After Solve???\n");
    cout << flush;
  }
  //  assert(ok);
#ifdef HiGHSDEV
  //  printf("model->mlFg_Report() 9\n");cout<<flush;
  //  model->mlFg_Report();cout<<flush;
#endif

#ifdef HiGHSDEV
  if (model->anInvertTime) {
    printf(
        "Time: Total inverts =  %4d; Total invert  time = %11.4g of Total time "
        "= %11.4g",
        model->totalInverts, model->totalInvertTime, model->totalTime);
    if (model->totalTime > 0.001) {
      printf(" (%6.2f%%)\n", (100 * model->totalInvertTime) / model->totalTime);
    } else {
      printf("\n");
    }
    cout << flush;
    printf(
        "Time: Total rebuilds = %4d; Total rebuild time = %11.4g of Total time "
        "= %11.4g",
        totalRebuilds, totalRebuildTime, model->totalTime);
    if (model->totalTime > 0.001) {
      printf(" (%6.2f%%)\n", (100 * totalRebuildTime) / model->totalTime);
    } else {
      printf("\n");
    }
    cout << flush;
  }
  model->util_anMlSol();

#endif
}

void HDual::init(int num_threads) {
  // Copy size, matrix and factor
  numCol = model->getNumCol();
  numRow = model->getNumRow();
  numTot = model->getNumTot();
  matrix = model->getMatrix();
  factor = model->getFactor();

  // Copy pointers
  jMove = model->getNonbasicMove();
  workDual = model->getWorkDual();
  //    JAJH: Only because I can't get this from HModel.h
  workValue = model->getWorkValue();
  workRange = model->getWorkRange();
  baseLower = model->getBaseLower();
  baseUpper = model->getBaseUpper();
  baseValue = model->getBaseValue();

  // Copy tolerances
  Tp = model->dblOption[DBLOPT_PRIMAL_TOL];
  Td = model->dblOption[DBLOPT_DUAL_TOL];

  // Setup local vectors
  columnDSE.setup(numRow);
  columnBFRT.setup(numRow);
  column.setup(numRow);
  row_ep.setup(numRow);
  row_ap.setup(numCol);
  //  row_ap_ultra.setup(numCol);
  columnDensity = 0;
  row_epDensity = 0;
  row_apDensity = 0;
  rowdseDensity = 0;
  // Setup other buffers
  dualRow.setup(model);
  dualRHS.setup(model);

  // Initialize for tasks
  if (dual_variant == HDUAL_VARIANT_TASKS) {
    init_slice(num_threads - 2);
  }

  // Initialize for multi
  if (dual_variant == HDUAL_VARIANT_MULTI) {
    multi_num = num_threads;
    if (multi_num < 1) multi_num = 1;
    if (multi_num > HSOL_THREAD_LIMIT) multi_num = HSOL_THREAD_LIMIT;
    for (int i = 0; i < multi_num; i++) {
      multi_choice[i].row_ep.setup(numRow);
      multi_choice[i].column.setup(numRow);
      multi_choice[i].columnBFRT.setup(numRow);
    }
    init_slice(multi_num - 1);
  }
  multi_iteration = 0;
  //  string partitionFile = model->strOption[STROPT_PARTITION_FILE];
  //  if (partitionFile.size())
  //  {
  //    dualRHS.setup_partition(partitionFile.c_str());
  //  }
}

void HDual::init_slice(int init_sliced_num) {
  // Number of slices
  slice_num = init_sliced_num;
  if (slice_num < 1) slice_num = 1;
  if (slice_num > HSOL_SLICED_LIMIT) slice_num = HSOL_SLICED_LIMIT;

  // Alias to the matrix
  const int *Astart = matrix->getAstart();
  const int *Aindex = matrix->getAindex();
  const double *Avalue = matrix->getAvalue();
  const int AcountX = Astart[numCol];

  // Figure out partition weight
  double sliced_countX = AcountX / slice_num;
  slice_start[0] = 0;
  for (int i = 0; i < slice_num - 1; i++) {
    int endColumn = slice_start[i] + 1;  // At least one column
    int endX = Astart[endColumn];
    int stopX = (i + 1) * sliced_countX;
    while (endX < stopX) {
      endX = Astart[++endColumn];
    }
    slice_start[i + 1] = endColumn;
    if (endColumn >= numCol) {
      slice_num = i;  // SHRINK
      break;
    }
  }
  slice_start[slice_num] = numCol;

  // Partition the matrix, row_ap and related packet
  vector<int> sliced_Astart;
  for (int i = 0; i < slice_num; i++) {
    // The matrix
    int mystart = slice_start[i];
    int mycount = slice_start[i + 1] - mystart;
    int mystartX = Astart[mystart];
    sliced_Astart.resize(mycount + 1);
    for (int k = 0; k <= mycount; k++)
      sliced_Astart[k] = Astart[k + mystart] - mystartX;
    // TODO generalise this call so slice can be used with non-logical initial
    // basis
    slice_matrix[i].setup_lgBs(mycount, numRow, &sliced_Astart[0],
                               Aindex + mystartX, Avalue + mystartX);

    // The row_ap and its packages
    slice_row_ap[i].setup(mycount);
    slice_dualRow[i].setupSlice(model, mycount);
  }
}

void HDual::solve_phase1() {
  model->util_reportMessage("dual-phase-1-start");
  // Switch to dual phase 1 bounds
  model->initBound(1);
  model->initValue();
  double lc_totalTime = model->totalTime + model->timer.getTime();
#ifdef HiGHSDEV
  // int lc_totalTime_rp_n = 0; printf("DualPh1: lc_totalTime = %5.2f; Record
  // %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
  // Main solving structure
  model->timer.recordStart(HTICK_ITERATE);
  for (;;) {
    model->timer.recordStart(HTICK_ITERATE_REBUILD);
    rebuild();
    model->timer.recordFinish(HTICK_ITERATE_REBUILD);
    for (;;) {
      switch (dual_variant) {
        default:
        case HDUAL_VARIANT_PLAIN:
          iterate();
          break;
        case HDUAL_VARIANT_TASKS:
          iterate_tasks();
          break;
        case HDUAL_VARIANT_MULTI:
          iterate_multi();
          break;
      }
      if (invertHint) break;
      // printf("HDual::solve_phase1: Iter = %d; Objective = %g\n",
      // model->numberIteration, model->dualObjectiveValue);
      /*
      if (model->dualObjectiveValue > model->dblOption[DBLOPT_OBJ_UB]) {
#ifdef SCIP_DEV
        printf("HDual::solve_phase1: %g = Objective >
dblOption[DBLOPT_OBJ_UB]\n", model->dualObjectiveValue, model->dblOption[DBLOPT_OBJ_UB]);
#endif
        model->problemStatus = LP_Status_ObjUB;
        break;
      }
      */
    }
    lc_totalTime = model->totalTime + model->timer.getTime();
#ifdef HiGHSDEV
    //      lc_totalTime_rp_n += 1; printf("DualPh1: lc_totalTime = %5.2f;
    //      Record %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
    if (lc_totalTime > TimeLimitValue) {
      SolveBailout = true;
      model->problemStatus = LP_Status_OutOfTime;
      break;
    }
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (model->countUpdate == 0) break;
    if (model->mlFg_haveFreshRebuild) break;
  }

  model->timer.recordFinish(HTICK_ITERATE);
  if (SolveBailout) return;

  if (rowOut == -1) {
    model->util_reportMessage("dual-phase-1-optimal");
    // Go to phase 2
    if (model->dualObjectiveValue == 0) {
      solvePhase = 2;
    } else {
      // We still have dual infeasible
      if (model->problemPerturbed) {
        // Clean up perturbation and go on
        cleanup();
        if (dualInfeasCount == 0) solvePhase = 2;
      } else {
        // Report dual infeasible
        solvePhase = -1;
        model->util_reportMessage("dual-infeasible");
        model->setProblemStatus(LP_Status_Unbounded);
      }
    }
  } else if (invertHint == invertHint_chooseColumnFail) {
    // chooseColumn has failed
    // Behave as "Report strange issues" below
    solvePhase = -1;
    model->util_reportMessage("dual-phase-1-not-solved");
    model->setProblemStatus(LP_Status_Failed);
  } else if (columnIn == -1) {
    // We got dual phase 1 unbounded - strange
    model->util_reportMessage("dual-phase-1-unbounded");
    if (model->problemPerturbed) {
      // Clean up perturbation and go on
      cleanup();
      if (dualInfeasCount == 0) solvePhase = 2;
    } else {
      // Report strange issues
      solvePhase = -1;
      model->util_reportMessage("dual-phase-1-not-solved");
      model->setProblemStatus(LP_Status_Failed);
    }
  }

  if (solvePhase == 2) {
    model->initBound();
    model->initValue();
  }
}

void HDual::solve_phase2() {
  model->util_reportMessage("dual-phase-2-start");

  // Collect free variables
  dualRow.create_Freelist();
  double lc_totalTime = model->totalTime + model->timer.getTime();
#ifdef HiGHSDEV
  //  int lc_totalTime_rp_n = 0; printf("DualPh2: lc_totalTime = %5.2f; Record
  //  %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
  // Main solving structure
  model->timer.recordStart(HTICK_ITERATE);
  for (;;) {
    // Outer loop of solve_phase2()
    // Rebuild all values, reinverting B if updates have been performed
    model->timer.recordStart(HTICK_ITERATE_REBUILD);
    rebuild();
    model->timer.recordFinish(HTICK_ITERATE_REBUILD);
    if (dualInfeasCount > 0) break;
    for (;;) {
      // Inner loop of solve_phase2()
      // Performs one iteration in case HDUAL_VARIANT_PLAIN:
      model->util_reportSolverProgress();
      switch (dual_variant) {
        default:
        case HDUAL_VARIANT_PLAIN:
          iterate();
          break;
        case HDUAL_VARIANT_TASKS:
          iterate_tasks();
          break;
        case HDUAL_VARIANT_MULTI:
          iterate_multi();
          break;
      }
      // invertHint can be true for various reasons see HModel.h
      if (invertHint) break;
        // Need the dual objective value to check for exceeding the
        // upper bound set by SCIP, but can't afford to compute it
        // every iteration!
        // Killer line for speed of HiGHS on hyper-sparse LPs!
        // Comment out when not working with SCIP!!
        //	  model->computeDualObjectiveValue();
#ifdef HiGHSDEV
        //      model->computeDualObjectiveValue();
        //      double pr_obj_v = model->computePrObj();
        //      printf("HDual::solve_phase2: Iter = %4d; Pr Obj = %.11g; Du Obj
        //      = %.11g\n",
        //	     model->numberIteration, pr_obj_v, model->dualObjectiveValue);
#endif
      if (model->dualObjectiveValue > model->dblOption[DBLOPT_OBJ_UB]) {
#ifdef SCIP_DEV
        printf(
            "HDual::solve_phase2: Objective = %g > %g = "
            "dblOption[DBLOPT_OBJ_UB]\n",
            model->dualObjectiveValue, model->dblOption[DBLOPT_OBJ_UB]);
#endif
        model->problemStatus = LP_Status_ObjUB;
        SolveBailout = true;
        break;
      }
    }
    lc_totalTime = model->totalTime + model->timer.getTime();
    if (model->problemStatus == LP_Status_ObjUB) {
      SolveBailout = true;
      break;
    }
#ifdef HiGHSDEV
    //      lc_totalTime_rp_n += 1; printf("DualPh2: lc_totalTime = %5.2f;
    //      Record %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
    if (lc_totalTime > TimeLimitValue) {
      model->problemStatus = LP_Status_OutOfTime;
      SolveBailout = true;
      break;
    }
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (model->countUpdate == 0) break;
    if (model->mlFg_haveFreshRebuild) break;
  }
  model->timer.recordFinish(HTICK_ITERATE);

  if (SolveBailout) {
    return;
  }
  if (dualInfeasCount > 0) {
    // There are dual infeasiblities so switch to Phase 1 and return
    model->util_reportMessage("dual-phase-2-found-free");
    solvePhase = 1;
  } else if (rowOut == -1) {
    // There is no candidate in CHUZR, even after rebuild so probably optimal
    model->util_reportMessage("dual-phase-2-optimal");
    //		printf("Rebuild: cleanup()\n");
    cleanup();
    if (dualInfeasCount > 0) {
      // There are dual infeasiblities after cleanup() so switch to primal
      // simplex
      solvePhase = 4;  // Do primal
    } else {
      // There are no dual infeasiblities after cleanup() so optimal!
      solvePhase = 0;
      model->util_reportMessage("problem-optimal");
      model->setProblemStatus(LP_Status_Optimal);
    }
  } else if (invertHint == invertHint_chooseColumnFail) {
    // chooseColumn has failed
    // Behave as "Report strange issues" below
    solvePhase = -1;
    model->util_reportMessage("dual-phase-2-not-solved");
    model->setProblemStatus(LP_Status_Failed);
  } else if (columnIn == -1) {
    // There is no candidate in CHUZC, so probably dual unbounded
    model->util_reportMessage("dual-phase-2-unbounded");
    if (model->problemPerturbed) {
      // If the costs have been perturbed, clean up and return
      cleanup();
    } else {
      // If the costs have not been perturbed, so dual unbounded---and hence
      // primal infeasible
      solvePhase = -1;
      model->util_reportMessage("problem-infeasible");
      model->setProblemStatus(LP_Status_Infeasible);
    }
  }
}

void HDual::rebuild() {
  // Save history information
  model->recordPivots(-1, -1, 0);  // Indicate REINVERT
#ifdef HiGHSDEV
  double tt0 = 0;
  if (anRebuildTime) tt0 = model->timer.getTime();
#endif
  int sv_invertHint = invertHint;
  invertHint = invertHint_no;  // Was 0

  // Possibly Rebuild model->factor
  bool reInvert = model->countUpdate > 0;
  if (!model->InvertIfRowOutNeg) {
    // Don't reinvert if rowOut is negative [equivalently, if sv_invertHint ==
    // invertHint_possiblyOptimal]
    if (sv_invertHint == invertHint_possiblyOptimal) {
      assert(rowOut == -1);
      reInvert = false;
    }
  }
  if (reInvert) {
    const int *baseIndex = model->getBaseIndex();
    // Scatter the edge weights so that, after INVERT,
    // they can be gathered according to the new

    // permutation of baseIndex
    model->timer.recordStart(HTICK_PERM_WT);
    for (int i = 0; i < numRow; i++)
      dualRHS.workEdWtFull[baseIndex[i]] = dualRHS.workEdWt[i];
    model->timer.recordFinish(HTICK_PERM_WT);

    model->timer.recordStart(HTICK_INVERT);

    // Call computeFactor to perform INVERT
    int rankDeficiency = model->computeFactor();
    model->timer.recordFinish(HTICK_INVERT);

    if (rankDeficiency)
      throw runtime_error("Dual reInvert: singular-basis-matrix");
    // Gather the edge weights according to the
    // permutation of baseIndex after INVERT
    model->timer.recordStart(HTICK_PERM_WT);
    for (int i = 0; i < numRow; i++)
      dualRHS.workEdWt[i] = dualRHS.workEdWtFull[baseIndex[i]];
    model->timer.recordFinish(HTICK_PERM_WT);

    // Possibly look at the basis condition
    //		double bsCond = an_bs_cond(model);
  }

  // Recompute dual solution
  model->timer.recordStart(HTICK_COMPUTE_DUAL);
  model->computeDual();
  model->timer.recordFinish(HTICK_COMPUTE_DUAL);

  model->timer.recordStart(HTICK_CORRECT_DUAL);
  model->correctDual(&dualInfeasCount);
  model->timer.recordFinish(HTICK_CORRECT_DUAL);

  // Recompute primal solution
  model->timer.recordStart(HTICK_COMPUTE_PRIMAL);
  model->computePrimal();
  model->timer.recordFinish(HTICK_COMPUTE_PRIMAL);

  // Collect primal infeasible as a list
  model->timer.recordStart(HTICK_COLLECT_PR_IFS);
  dualRHS.create_infeasArray();
  dualRHS.create_infeasList(columnDensity);
  model->timer.recordFinish(HTICK_COLLECT_PR_IFS);

  // Compute the objective value
  model->timer.recordStart(HTICK_COMPUTE_DUOBJ);
  model->computeDualObjectiveValue(solvePhase);

#ifdef HiGHSDEV
  model->checkDualObjectiveValue("After model->computeDualObjectiveValue");
  printf("Checking INVERT in rebuild()\n");
  model->factor.checkInvert();
#endif

  model->timer.recordFinish(HTICK_COMPUTE_DUOBJ);
  //	model->util_reportNumberIterationObjectiveValue(sv_invertHint);

  model->timer.recordStart(HTICK_REPORT_INVERT);
  iterateRpInvert(sv_invertHint);
  model->timer.recordFinish(HTICK_REPORT_INVERT);

  total_INVERT_TICK = factor->build_syntheticTick;  // Was factor->pseudoTick
  total_FT_inc_TICK = 0;
#ifdef HiGHSDEV
  total_fake = 0;
#endif
  total_syntheticTick = 0;

#ifdef HiGHSDEV
  if (anRebuildTime) {
    double rebuildTime = model->timer.getTime() - tt0;
    totalRebuilds++;
    totalRebuildTime += rebuildTime;
    printf(
        "Dual  Ph%-2d rebuild %4d (%1d) on iteration %9d: Rebuild time = "
        "%11.4g; Total rebuild time = %11.4g\n",
        solvePhase, totalRebuilds, sv_invertHint, model->numberIteration,
        rebuildTime, totalRebuildTime);
  }
#endif
  // Data are fresh from rebuild
  model->mlFg_haveFreshRebuild = 1;
}

void HDual::cleanup() {
  // Remove perturbation and recompute the dual solution
  model->util_reportMessage("dual-cleanup-shift");
  model->initCost();
  model->initBound();
  model->computeDual();
  model->computeDualObjectiveValue(solvePhase);
  //	model->util_reportNumberIterationObjectiveValue(-1);
  iterateRpInvert(-1);

  model->computeDualInfeasInPrimal(&dualInfeasCount);
}

void HDual::iterate() {
  // This is the main teration loop for dual revised simplex. All the
  // methods have as their first line if (invertHint) return;, where
  // invertHint is, for example, set to 1 when CHUZR finds no
  // candidate. This causes a break from the inner loop of
  // solve_phase% and, hence, a call to rebuild()

  //	Reporting:
  //	hsol row-wise matrix after update in updateMatrix(columnIn, columnOut);
  model->timer.recordStart(HTICK_ITERATE_CHUZR);
  chooseRow();
  model->timer.recordFinish(HTICK_ITERATE_CHUZR);

  model->timer.recordStart(HTICK_ITERATE_CHUZC);
  chooseColumn(&row_ep);
  model->timer.recordFinish(HTICK_ITERATE_CHUZC);

  model->timer.recordStart(HTICK_ITERATE_FTRAN);
  updateFtranBFRT();
  // updateFtran(); computes the pivotal column in the data structure "column"
  updateFtran();

  // updateFtranDSE performs the DSE FTRAN on pi_p
  if (EdWt_Mode == EdWt_Mode_DSE) updateFtranDSE(&row_ep);
  model->timer.recordFinish(HTICK_ITERATE_FTRAN);

  // updateVerify() Checks row-wise pivot against column-wise pivot for
  // numerical trouble
  model->timer.recordStart(HTICK_ITERATE_VERIFY);
  updateVerify();
  model->timer.recordFinish(HTICK_ITERATE_VERIFY);

  // updateDual() Updates the dual values
  model->timer.recordStart(HTICK_ITERATE_DUAL);
  updateDual();
  model->timer.recordFinish(HTICK_ITERATE_DUAL);

  // updatePrimal(&row_ep); Updates the primal values and the edge weights
  model->timer.recordStart(HTICK_ITERATE_PRIMAL);
  updatePrimal(&row_ep);
  model->timer.recordFinish(HTICK_ITERATE_PRIMAL);

  if ((EdWt_Mode == EdWt_Mode_Dvx) && (nw_dvx_fwk)) {
    model->timer.recordStart(HTICK_ITERATE_DEVEX_IZ);
    iz_dvx_fwk();
    model->timer.recordFinish(HTICK_ITERATE_DEVEX_IZ);
  }

  // Update the basis representation
  model->timer.recordStart(HTICK_ITERATE_PIVOTS);
  updatePivots();
  model->timer.recordFinish(HTICK_ITERATE_PIVOTS);

  // Analyse the iteration: possibly report; possibly switch strategy
  iterateAn();
}

void HDual::iterate_tasks() {
  slice_PRICE = 1;

  // Group 1
  chooseRow();

  // Disable slice when too sparse
  if (1.0 * row_ep.count / numRow < 0.01) slice_PRICE = 0;

  model->timer.recordStart(HTICK_GROUP1);
#pragma omp parallel
#pragma omp single
  {
#pragma omp task
    {
      columnDSE.copy(&row_ep);
      updateFtranDSE(&columnDSE);
    }
#pragma omp task
    {
      if (slice_PRICE)
        chooseColumn_slice(&row_ep);
      else
        chooseColumn(&row_ep);
#pragma omp task
      updateFtranBFRT();
#pragma omp task
      updateFtran();
#pragma omp taskwait
    }
  }
  model->timer.recordFinish(HTICK_GROUP1);

  updateVerify();
  updateDual();
  updatePrimal(&columnDSE);
  updatePivots();
}

void HDual::iterateIzAn() {
  AnIterIt0 = model->numberIteration;
  AnIterCostlyDseFq = 0;
#ifdef HiGHSDEV
  AnIterPrevRpNumCostlyDseIt = 0;
  AnIterPrevIt = 0;
  AnIterOpRec *AnIter;
  AnIter = &AnIterOp[AnIterOpTy_Btran];
  AnIter->AnIterOpName = "Btran";
  AnIter = &AnIterOp[AnIterOpTy_Price];
  AnIter->AnIterOpName = "Price";
  AnIter = &AnIterOp[AnIterOpTy_Ftran];
  AnIter->AnIterOpName = "Ftran";
  AnIter = &AnIterOp[AnIterOpTy_FtranBFRT];
  AnIter->AnIterOpName = "FtranBFRT";
  AnIter = &AnIterOp[AnIterOpTy_FtranDSE];
  AnIter->AnIterOpName = "FtranDSE";
  for (int k = 0; k < NumAnIterOpTy; k++) {
    AnIter = &AnIterOp[k];
    AnIter->AnIterOpLog10RsDsty = 0;
    AnIter->AnIterOpSuLog10RsDsty = 0;
    if (k == AnIterOpTy_Price) {
      AnIter->AnIterOpHyperCANCEL = 1.0;
      AnIter->AnIterOpHyperTRAN = 1.0;
      AnIter->AnIterOpRsDim = numCol;
    } else {
      if (k == AnIterOpTy_Btran) {
        AnIter->AnIterOpHyperCANCEL = hyperCANCEL;
        AnIter->AnIterOpHyperTRAN = hyperBTRANU;
      } else {
        AnIter->AnIterOpHyperCANCEL = hyperCANCEL;
        AnIter->AnIterOpHyperTRAN = hyperFTRANL;
      }
      AnIter->AnIterOpRsDim = numRow;
    }
    AnIter->AnIterOpNumCa = 0;
    AnIter->AnIterOpNumHyperOp = 0;
    AnIter->AnIterOpNumHyperRs = 0;
    AnIter->AnIterOpRsMxNNZ = 0;
    AnIter->AnIterOpSuNumCa = 0;
    AnIter->AnIterOpSuNumHyperOp = 0;
    AnIter->AnIterOpSuNumHyperRs = 0;
  }
  for (int k = 1; k <= AnIterNumInvertHint; k++) AnIterNumInvert[k] = 0;
  AnIterNumPrDgnIt = 0;
  AnIterNumDuDgnIt = 0;
  AnIterNumColPrice = 0;
  AnIterNumRowPrice = 0;
  AnIterNumRowPriceWSw = 0;
  AnIterNumRowPriceUltra = 0;
  for (int k = 0; k <= EdWt_Mode_Dan; k++) AnIterNumEdWtIt[k] = 0;
  AnIterNumCostlyDseIt = 0;
  AnIterTraceNumRec = 0;
  AnIterTraceIterDl = 1;
  AnIterTraceRec *lcAnIter = &AnIterTrace[0];
  lcAnIter->AnIterTraceIter = AnIterIt0;
  lcAnIter->AnIterTraceTime = model->timer.getTime();
#endif
}

void HDual::iterateAn() {
  // Possibly report on the iteration
  iterateRp();

  // Possibly switch from DSE to Dvx
  if (EdWt_Mode == EdWt_Mode_DSE) {
    double AnIterCostlyDseMeasureDen;
    //    AnIterCostlyDseMeasureDen = row_epDensity*columnDensity;
    AnIterCostlyDseMeasureDen =
        max(max(row_epDensity, columnDensity), row_apDensity);
    if (AnIterCostlyDseMeasureDen > 0) {
      AnIterCostlyDseMeasure = rowdseDensity / AnIterCostlyDseMeasureDen;
      AnIterCostlyDseMeasure = AnIterCostlyDseMeasure * AnIterCostlyDseMeasure;
    } else {
      AnIterCostlyDseMeasure = 0;
    }
    bool CostlyDseIt = AnIterCostlyDseMeasure > AnIterCostlyDseMeasureLimit &&
                       rowdseDensity > AnIterCostlyDseMnDensity;
    AnIterCostlyDseFq = (1 - runningAverageMu) * AnIterCostlyDseFq;
    if (CostlyDseIt) {
      AnIterNumCostlyDseIt++;
      AnIterCostlyDseFq += runningAverageMu * 1.0;
      int lcNumIter = model->numberIteration - AnIterIt0;
      if (alw_DSE2Dvx_sw &&
          (AnIterNumCostlyDseIt > lcNumIter * AnIterFracNumCostlyDseItbfSw) &&
          (lcNumIter > AnIterFracNumTot_ItBfSw * numTot)) {
        // At least 5% of the (at least) 0.1NumTot iterations have been costly
        // DSE so switch to Devex
#ifdef HiGHSDEV
        printf(
            "Switch from DSE to Dvx after %d costly DSE iterations of %d: "
            "Col_Dsty = %11.4g; R_Ep_Dsty = %11.4g; DSE_Dsty = %11.4g\n",
            AnIterNumCostlyDseIt, lcNumIter, rowdseDensity, row_epDensity,
            columnDensity);
#endif
        EdWt_Mode = EdWt_Mode_Dvx;
        // Zero the number of Devex frameworks used and set up the first one
        n_dvx_fwk = 0;
        dvx_ix.assign(numTot, 0);
        iz_dvx_fwk();
      }
    }
  }

#ifdef HiGHSDEV
  int AnIterCuIt = model->numberIteration;
  bool iterLg = AnIterCuIt % 100 == 0;
  iterLg = false;
  if (iterLg) {
    int lc_NumCostlyDseIt = AnIterNumCostlyDseIt - AnIterPrevRpNumCostlyDseIt;
    AnIterPrevRpNumCostlyDseIt = AnIterNumCostlyDseIt;
    printf("Iter %10d: ", AnIterCuIt);
    iterateRpDsty(true);
    iterateRpDsty(false);
    if (EdWt_Mode == EdWt_Mode_DSE) {
      int lc_pct = (100 * AnIterNumCostlyDseIt) / (AnIterCuIt - AnIterIt0);
      printf("| Fq = %4.2f; Su =%5d (%3d%%)", AnIterCostlyDseFq,
             AnIterNumCostlyDseIt, lc_pct);

      if (lc_NumCostlyDseIt > 0) printf("; LcNum =%3d", lc_NumCostlyDseIt);
    }
    printf("\n");
  }

  for (int k = 0; k < NumAnIterOpTy; k++) {
    AnIterOpRec *lcAnIterOp = &AnIterOp[k];
    if (lcAnIterOp->AnIterOpNumCa) {
      lcAnIterOp->AnIterOpSuNumCa += lcAnIterOp->AnIterOpNumCa;
      lcAnIterOp->AnIterOpSuNumHyperOp += lcAnIterOp->AnIterOpNumHyperOp;
      lcAnIterOp->AnIterOpSuNumHyperRs += lcAnIterOp->AnIterOpNumHyperRs;
      lcAnIterOp->AnIterOpSuLog10RsDsty += lcAnIterOp->AnIterOpLog10RsDsty;
    }
    lcAnIterOp->AnIterOpNumCa = 0;
    lcAnIterOp->AnIterOpNumHyperOp = 0;
    lcAnIterOp->AnIterOpNumHyperRs = 0;
    lcAnIterOp->AnIterOpLog10RsDsty = 0;
  }
  if (invertHint > 0) AnIterNumInvert[invertHint]++;
  if (thetaDual <= 0) AnIterNumDuDgnIt++;
  if (thetaPrimal <= 0) AnIterNumPrDgnIt++;
  if (AnIterCuIt > AnIterPrevIt)
    AnIterNumEdWtIt[EdWt_Mode] += (AnIterCuIt - AnIterPrevIt);

  AnIterTraceRec *lcAnIter = &AnIterTrace[AnIterTraceNumRec];
  //  if (model->numberIteration ==
  //  AnIterTraceIterRec[AnIterTraceNumRec]+AnIterTraceIterDl) {
  if (model->numberIteration == lcAnIter->AnIterTraceIter + AnIterTraceIterDl) {
    if (AnIterTraceNumRec == AnIterTraceMxNumRec) {
      for (int rec = 1; rec <= AnIterTraceMxNumRec / 2; rec++)
        AnIterTrace[rec] = AnIterTrace[2 * rec];
      AnIterTraceNumRec = AnIterTraceNumRec / 2;
      AnIterTraceIterDl = AnIterTraceIterDl * 2;
    } else {
      AnIterTraceNumRec++;
      lcAnIter = &AnIterTrace[AnIterTraceNumRec];
      lcAnIter->AnIterTraceIter = model->numberIteration;
      lcAnIter->AnIterTraceTime = model->timer.getTime();
      lcAnIter->AnIterTraceDsty[AnIterOpTy_Btran] = row_epDensity;
      lcAnIter->AnIterTraceDsty[AnIterOpTy_Price] = row_apDensity;
      lcAnIter->AnIterTraceDsty[AnIterOpTy_Ftran] = columnDensity;
      lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranBFRT] = columnDensity;
      if (EdWt_Mode == EdWt_Mode_DSE) {
        lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE] = rowdseDensity;
        lcAnIter->AnIterTraceAux0 = AnIterCostlyDseMeasure;
      } else {
        lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE] = 0;
        lcAnIter->AnIterTraceAux0 = 0;
      }
      lcAnIter->AnIterTraceEdWt_Mode = EdWt_Mode;
    }
  }
  AnIterPrevIt = AnIterCuIt;
#endif
}

void HDual::iterateRp() {
  if (model->intOption[INTOPT_PRINT_FLAG] != 4) return;
  int numIter = model->numberIteration;
  bool header = numIter % 10 == 1;
  header = true;  // JAJH10/10
  if (header) iterateRpFull(header);
  iterateRpFull(false);
}

void HDual::iterateRpFull(bool header) {
  if (header) {
    iterateRpIterPh(true);
    iterateRpDuObj(true);
#ifdef HiGHSDEV
    iterateRpIterDa(true);
    iterateRpDsty(true);
    printf(" FreeLsZ");
#endif
    printf("\n");
  } else {
    iterateRpIterPh(false);
    iterateRpDuObj(false);
#ifdef HiGHSDEV
    iterateRpIterDa(false);
    iterateRpDsty(false);
    printf(" %7d", dualRow.freeListSize);
#endif
    printf("\n");
  }
}

void HDual::iterateRpIterPh(bool header) {
  if (header) {
    printf(" Iteration Ph");
  } else {
    int numIter = model->numberIteration;
    printf(" %9d %2d", numIter, solvePhase);
  }
}
void HDual::iterateRpDuObj(bool header) {
  if (header) {
    printf("    DualObjective    ");
  } else {
    model->computeDualObjectiveValue(solvePhase);
    printf(" %20.10e", model->dualObjectiveValue);
  }
}

void HDual::iterateRpInvert(int i_v) {
  if (model->intOption[INTOPT_PRINT_FLAG] != 1 &&
      model->intOption[INTOPT_PRINT_FLAG] != 4)
    return;
  printf("Iter %10d:", model->numberIteration);
#ifdef HiGHSDEV
  iterateRpDsty(true);
  iterateRpDsty(false);
#endif
  iterateRpDuObj(false);
  printf(" %2d\n", i_v);
}

void HDual::uOpRsDensityRec(double lc_OpRsDensity, double &opRsDensity) {
  // Update an average density record for BTRAN, an FTRAN or PRICE
  opRsDensity = (1 - runningAverageMu) * (opRsDensity) +
                runningAverageMu * lc_OpRsDensity;
}

void HDual::chooseRow() {
  // Choose the index of a row to leave the basis (CHUZR)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  // Choose candidates repeatedly until candidate is OK or optimality is
  // detected
  for (;;) {
    // Choose the index of a good row to leave the basis
    dualRHS.choose_normal(&rowOut);
    if (rowOut == -1) {
      // No index found so may be dual optimal. By setting
      // invertHint>0 all subsequent methods in the iteration will
      // be skipped until reinversion and rebuild have taken place
      invertHint = invertHint_possiblyOptimal;
      return;
    }
    // Compute pi_p = B^{-T}e_p in row_ep
    model->timer.recordStart(HTICK_BTRAN);
    // Set up RHS for BTRAN
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = rowOut;
    row_ep.array[rowOut] = 1;
    row_ep.packFlag = true;
#ifdef HiGHSDEV
    if (AnIterLg) iterateOpRecBf(AnIterOpTy_Btran, row_ep, row_epDensity);
#endif
    // Perform BTRAN
    //      printf("\nBTRAN\n");
    factor->btran(row_ep, row_epDensity);
#ifdef HiGHSDEV
    if (AnIterLg) iterateOpRecAf(AnIterOpTy_Btran, row_ep);
#endif
    model->timer.recordFinish(HTICK_BTRAN);
    // Verify DSE weight
    if (EdWt_Mode == EdWt_Mode_DSE) {
      // For DSE, see how accurate the updated weight is
      // Save the updated weight
      double u_weight = dualRHS.workEdWt[rowOut];
      // Compute the weight from row_ep and over-write the updated weight
      double c_weight = dualRHS.workEdWt[rowOut] = row_ep.norm2();
      // If the weight error is acceptable then break out of the
      // loop. All we worry about is accepting rows with weights
      // which are not too small, since this can make the row look
      // unreasonably attractive
      if (u_weight >= 0.25 * c_weight) break;
#ifdef HiGHSDEV
      // Count the number of wrong DSE weights for reporting
      n_wg_DSE_wt += 1;
#endif
      // Weight error is unacceptable so look for another
      // candidate. Of course, it's possible that the same
      // candidate is chosen, but the weight will be correct (so
      // no infinite loop).
    } else {
      // If not using DSE then accept the row by breaking out of
      // the loop
      break;
    }
  }
  // Index of row to leave the basis has been found
  //
  // Assign basic info:
  //
  // Record the column (variable) associated with the leaving row
  columnOut = model->getBaseIndex()[rowOut];
  // Record the change in primal variable associated with the move to the bound
  // being violated
  if (baseValue[rowOut] < baseLower[rowOut]) {
    // Below the lower bound so set deltaPrimal = value - LB < 0
    deltaPrimal = baseValue[rowOut] - baseLower[rowOut];
  } else {
    // Above the upper bound so set deltaPrimal = value - UB > 0
    deltaPrimal = baseValue[rowOut] - baseUpper[rowOut];
  }
  // Set sourceOut to be -1 if deltaPrimal<0, otherwise +1 (since deltaPrimal>0)
  sourceOut = deltaPrimal < 0 ? -1 : 1;
  // Update the record of average row_ep (pi_p) density. This ignores
  // any BTRANs done for skipped candidates
  double lc_OpRsDensity = (double)row_ep.count / numRow;
  uOpRsDensityRec(lc_OpRsDensity, row_epDensity);
}

void HDual::chooseColumn(HVector *row_ep) {
  // Compute pivot row (PRICE) and choose the index of a column to enter the
  // basis (CHUZC)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  //
  // PRICE
  //
  model->timer.recordStart(HTICK_PRICE);
  row_ap.clear();

#ifdef HiGHSDEV
  bool anPriceEr = false;
  bool useUltraPrice = alw_price_ultra &&
                       row_apDensity * numCol * 10 < row_ap.ilP2 &&
                       row_apDensity < 1e-3;
#endif
  if (Price_Mode == Price_Mode_Col) {
    // Column-wise PRICE
#ifdef HiGHSDEV
    if (AnIterLg) {
      iterateOpRecBf(AnIterOpTy_Price, *row_ep, 0.0);
      AnIterNumColPrice++;
    }
#endif
    // Perform column-wise PRICE
    matrix->price_by_col(row_ap, *row_ep);
  } else {
    // By default, use row-wise PRICE, but possibly use column-wise
    // PRICE if the density of row_ep is too high
    double lc_dsty = (double)(*row_ep).count / numRow;
    if (alw_price_by_col_sw && (lc_dsty > dstyColPriceSw)) {
      // Use column-wise PRICE due to density of row_ep
#ifdef HiGHSDEV
      if (AnIterLg) {
        iterateOpRecBf(AnIterOpTy_Price, *row_ep, 0.0);
        AnIterNumColPrice++;
      }
#endif
      // Perform column-wise PRICE
      matrix->price_by_col(row_ap, *row_ep);
      // Zero the components of row_ap corresponding to basic variables
      // (nonbasicFlag[*]=0)
      for (int col = 0; col < numCol; col++) {
        row_ap.array[col] = model->nonbasicFlag[col] * row_ap.array[col];
      }
#ifdef HiGHSDEV
      // Ultra-sparse PRICE is in development
    } else if (useUltraPrice) {
      if (AnIterLg) {
        iterateOpRecBf(AnIterOpTy_Price, *row_ep, row_apDensity);
        AnIterNumRowPriceUltra++;
      }
      // Perform ultra-sparse row-wise PRICE
      matrix->price_by_row_ultra(row_ap, *row_ep);
      if (anPriceEr) {
        bool price_er;
        price_er = matrix->price_er_ck(row_ap, *row_ep);
        if (!price_er) printf("No ultra PRICE error\n");
      }
#endif
    } else if (alw_price_by_row_sw) {
      // Avoid hyper-sparse PRICE on current density of result or
      // switch if the density of row_ap becomes extreme
#ifdef HiGHSDEV
      if (AnIterLg) {
        iterateOpRecBf(AnIterOpTy_Price, *row_ep, row_apDensity);
        AnIterNumRowPriceWSw++;
      }
#endif
      // Set the value of the density of row_ap at which the switch to
      // sparse row-wise PRICE should be made
      const double sw_dsty = matrix->price_by_row_sw_dsty;
      // Perform hyper-sparse row-wise PRICE with switching
      matrix->price_by_row_w_sw(row_ap, *row_ep, row_apDensity, 0, sw_dsty);
    } else {
      // No avoiding hyper-sparse PRICE on current density of result
      // or switch if the density of row_ap becomes extreme
#ifdef HiGHSDEV
      if (AnIterLg) {
        iterateOpRecBf(AnIterOpTy_Price, *row_ep, 0.0);
        AnIterNumRowPrice++;
      }
#endif
      // Perform hyper-sparse row-wise PRICE
      matrix->price_by_row(row_ap, *row_ep);
    }
  }
#ifdef HiGHSDEV
  // Possibly analyse the error in the result of PRICE
  if (anPriceEr) {
    matrix->price_er_ck(row_ap, *row_ep);
  }
#endif
  // Update the record of average row_ap density
  double lc_OpRsDensity = (double)row_ap.count / numCol;
  uOpRsDensityRec(lc_OpRsDensity, row_apDensity);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecAf(AnIterOpTy_Price, row_ap);
#endif
  model->timer.recordFinish(HTICK_PRICE);
  //
  // CHUZC
  //
  // Section 0: Clear data and call create_Freemove to set a value of
  // nonbasicMove for all free columns to prevent their dual values
  // from being changed.
  model->timer.recordStart(HTICK_CHUZC0);
  dualRow.clear();
  dualRow.workDelta = deltaPrimal;
  dualRow.create_Freemove(row_ep);
  model->timer.recordFinish(HTICK_CHUZC0);
  //
  // Section 1: Pack row_ap and row_ep, then determine the possible
  // variables - candidates for CHUZC
  model->timer.recordStart(HTICK_CHUZC1);
  dualRow.choose_makepack(
      &row_ap, 0);  // Pack row_ap into the packIndex/Value of HDualRow
  dualRow.choose_makepack(
      row_ep, numCol);  // Pack row_ep into the packIndex/Value of HDualRow
  dualRow.choose_possible();  // Determine the possible variables - candidates
                              // for CHUZC
  model->timer.recordFinish(HTICK_CHUZC1);
  //
  // Take action if the step to an expanded bound is not positive, or
  // there are no candidates for CHUZC
  columnIn = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
    invertHint = invertHint_possiblyDualUnbounded;  // Was 1
    return;
  }
  //
  // Sections 2 and 3: Perform (bound-flipping) ratio test. This can
  // fail if the dual values are excessively large
  bool chooseColumnFail = dualRow.choose_final();
  if (chooseColumnFail) {
    invertHint = invertHint_chooseColumnFail;
    return;
  }
  //
  // Section 4: Reset the nonbasicMove values for free columns
  model->timer.recordStart(HTICK_CHUZC4);
  dualRow.delete_Freemove();
  model->timer.recordFinish(HTICK_CHUZC4);
  // Record values for basis change, checking for numerical problems and update
  // of dual variables
  columnIn = dualRow.workPivot;   // Index of the column entering the basis
  alphaRow = dualRow.workAlpha;   // Pivot value computed row-wise - used for
                                  // numerical checking
  thetaDual = dualRow.workTheta;  // Dual step length

  if (EdWt_Mode == EdWt_Mode_Dvx) {
    model->timer.recordStart(HTICK_DEVEX_WT);
    // Determine the exact Devex weight
    double og_dvx_wt_o_rowOut = dualRHS.workEdWt[rowOut];
    double tru_dvx_wt_o_rowOut = 0;
    // Loop over [row_ap; row_ep] using the packed values
    for (int el_n = 0; el_n < dualRow.packCount; el_n++) {
      int vr_n = dualRow.packIndex[el_n];
      double pv = dvx_ix[vr_n] * dualRow.packValue[el_n];
      tru_dvx_wt_o_rowOut += pv * pv;
    }
    tru_dvx_wt_o_rowOut = max(1.0, tru_dvx_wt_o_rowOut);
    // Analyse the Devex weight to determine whether a new framework
    // should be set up
    double dvx_rao = max(og_dvx_wt_o_rowOut / tru_dvx_wt_o_rowOut,
                         tru_dvx_wt_o_rowOut / og_dvx_wt_o_rowOut);
    int i_te = numRow / minRlvNumberDevexIterations;
    i_te = max(minAbsNumberDevexIterations, i_te);
    // Square maxAllowedDevexWeightRatio due to keeping squared
    // weights
    nw_dvx_fwk =
        dvx_rao > maxAllowedDevexWeightRatio * maxAllowedDevexWeightRatio ||
        n_dvx_it > i_te;
    dualRHS.workEdWt[rowOut] = tru_dvx_wt_o_rowOut;
    model->timer.recordFinish(HTICK_DEVEX_WT);
  }
  return;
}

void HDual::chooseColumn_slice(HVector *row_ep) {
  // Choose the index of a column to enter the basis (CHUZC) by
  // exploiting slices of the pivotal row - for SIP and PAMI
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  model->timer.recordStart(HTICK_CHUZC1);
  dualRow.clear();
  dualRow.workDelta = deltaPrimal;
  dualRow.create_Freemove(row_ep);

  // Row_ep:         PACK + CC1
#pragma omp task
  {
    dualRow.choose_makepack(row_ep, numCol);
    dualRow.choose_possible();
  }

  // Row_ap: PRICE + PACK + CC1
  for (int i = 0; i < slice_num; i++) {
#pragma omp task
    {
      slice_row_ap[i].clear();
      slice_matrix[i].price_by_row(slice_row_ap[i], *row_ep);

      slice_dualRow[i].clear();
      slice_dualRow[i].workDelta = deltaPrimal;
      slice_dualRow[i].choose_makepack(&slice_row_ap[i], slice_start[i]);
      slice_dualRow[i].choose_possible();
    }
  }
#pragma omp taskwait

  // Join CC1 results here
  for (int i = 0; i < slice_num; i++)
    dualRow.choose_joinpack(&slice_dualRow[i]);

  // Infeasible we created before
  columnIn = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
    invertHint = invertHint_possiblyDualUnbounded;  // Was 1
    model->timer.recordFinish(HTICK_CHUZC1);
    return;
  }
  model->timer.recordFinish(HTICK_CHUZC1);

  // Choose column 2, This only happens if didn't go out
  dualRow.choose_final();
  dualRow.delete_Freemove();
  columnIn = dualRow.workPivot;
  alphaRow = dualRow.workAlpha;
  thetaDual = dualRow.workTheta;
}

void HDual::updateFtran() {
  // Compute the pivotal column (FTRAN)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  model->timer.recordStart(HTICK_FTRAN);
  // Clear the picotal column and indicate that its values should be packed
  column.clear();
  column.packFlag = true;
  // Get the constraint matrix column by combining just one column
  // with unit multiplier
  matrix->collect_aj(column, columnIn, 1);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecBf(AnIterOpTy_Ftran, column, columnDensity);
#endif
  // Perform FTRAN
  //  printf("\nFTRAN\n");
  factor->ftran(column, columnDensity);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecAf(AnIterOpTy_Ftran, column);
#endif
  // Save the pivot value computed column-wise - used for numerical checking
  alpha = column.array[rowOut];
  model->timer.recordFinish(HTICK_FTRAN);
}

void HDual::updateFtranBFRT() {
  // Compute the RHS changes corresponding to the BFRT (FTRAN-BFRT)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  // Only time updateFtranBFRT if dualRow.workCount > 0;
  // If dualRow.workCount = 0 then dualRow.update_flip(&columnBFRT)
  // merely clears columnBFRT so no FTRAN is performed
  bool time_updateFtranBFRT = dualRow.workCount > 0;

  if (time_updateFtranBFRT) model->timer.recordStart(HTICK_FTRAN_BFRT);

  dualRow.update_flip(&columnBFRT);

  if (columnBFRT.count) {
#ifdef HiGHSDEV
    if (AnIterLg)
      iterateOpRecBf(AnIterOpTy_FtranBFRT, columnBFRT, columnDensity);
#endif
    // Perform FTRAN BFRT
    //    printf("\nFTRAN BFRT\n");
    factor->ftran(columnBFRT, columnDensity);
#ifdef HiGHSDEV
    if (AnIterLg) iterateOpRecAf(AnIterOpTy_FtranBFRT, columnBFRT);
#endif
  }
  if (time_updateFtranBFRT) model->timer.recordFinish(HTICK_FTRAN_BFRT);
}

void HDual::updateFtranDSE(HVector *DSE_Vector) {
  // Compute the vector required to update DSE weights - being FTRAN
  // applied to the pivotal column (FTRAN-DSE)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  model->timer.recordStart(HTICK_FTRAN_DSE);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecBf(AnIterOpTy_FtranDSE, *DSE_Vector, rowdseDensity);
#endif
  // Perform FTRAN DSE
  //  printf("\nFTRAN DSE\n");
  factor->ftran(*DSE_Vector, rowdseDensity);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecAf(AnIterOpTy_FtranDSE, *DSE_Vector);
#endif
  model->timer.recordFinish(HTICK_FTRAN_DSE);
}

void HDual::updateVerify() {
  // Compare the pivot value computed row-wise and column-wise and
  // determine whether reinversion is advisable
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  // Look at the relative difference between the absolute values of the two
  // pivot values
  double aCol = fabs(alpha);
  double aRow = fabs(alphaRow);
  double aDiff = fabs(aCol - aRow);
  numericalTrouble = aDiff / min(aCol, aRow);
  // Reinvert if the relative difference is large enough, and updates hav ebeen
  // performed
  if (numericalTrouble > 1e-7 && model->countUpdate > 0) {
    invertHint = invertHint_possiblySingularBasis;
  }
}

void HDual::updateDual() {
  // Update the dual values
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  // Update - dual (shift and back)
  if (thetaDual == 0)
    // Little to do if thetaDual is zero
    model->shiftCost(columnIn, -workDual[columnIn]);
  else {
    // Update the dual values (if packCount>0)
    dualRow.update_dual(thetaDual, columnOut);
    if (dual_variant != HDUAL_VARIANT_PLAIN && slice_PRICE) {
      // Update the dual variables slice-by-slice [presumably
      // nothing is done in the previous call to
      // dualRow.update_dual. TODO: Check with Qi
      for (int i = 0; i < slice_num; i++)
        slice_dualRow[i].update_dual(thetaDual, columnOut);
    }
  }
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;
  model->shiftBack(columnOut);
}

void HDual::updatePrimal(HVector *DSE_Vector) {
  // Update the primal values and any edge weights
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  // NB DSE_Vector is only computed if EdWt_Mode == EdWt_Mode_DSE
  // Update - primal and weight
  dualRHS.update_primal(&columnBFRT, 1);
  dualRHS.update_infeasList(&columnBFRT);
  double x_out = baseValue[rowOut];
  double l_out = baseLower[rowOut];
  double u_out = baseUpper[rowOut];
  thetaPrimal = (x_out - (deltaPrimal < 0 ? l_out : u_out)) / alpha;
  dualRHS.update_primal(&column, thetaPrimal);
  if (EdWt_Mode == EdWt_Mode_DSE) {
    double thisEdWt = dualRHS.workEdWt[rowOut] / (alpha * alpha);
    dualRHS.update_weight_DSE(&column, thisEdWt, -2 / alpha,
                              &DSE_Vector->array[0]);
    dualRHS.workEdWt[rowOut] = thisEdWt;
  } else if (EdWt_Mode == EdWt_Mode_Dvx) {
    // Pivotal row is for the current basis: weights are required for
    // the next basis so have to divide the current (exact) weight by
    // the pivotal value
    double thisEdWt = dualRHS.workEdWt[rowOut] / (alpha * alpha);
    double dvx_wt_o_rowOut = max(1.0, thisEdWt);
    // nw_wt is max(workEdWt[iRow], NewExactWeight*columnArray[iRow]^2);
    //
    // But NewExactWeight is dvx_wt_o_rowOut = max(1.0, dualRHS.workEdWt[rowOut] / (alpha * alpha))
    //
    // so nw_wt = max(workEdWt[iRow], dvx_wt_o_rowOut*columnArray[iRow]^2);
    //
    // Update rest of weights
    dualRHS.update_weight_Dvx(&column, dvx_wt_o_rowOut);
    dualRHS.workEdWt[rowOut] = dvx_wt_o_rowOut;
    n_dvx_it += 1;
  }
  dualRHS.update_infeasList(&column);

  if (EdWt_Mode == EdWt_Mode_DSE) {
    double lc_OpRsDensity = (double)DSE_Vector->count / numRow;
    uOpRsDensityRec(lc_OpRsDensity, rowdseDensity);
  }
  double lc_OpRsDensity = (double)column.count / numRow;
  uOpRsDensityRec(lc_OpRsDensity, columnDensity);

  //  total_fake += column.fakeTick;
  total_syntheticTick += column.syntheticTick;
  if (EdWt_Mode == EdWt_Mode_DSE) {
    //      total_fake += DSE_Vector->fakeTick;
    total_syntheticTick += DSE_Vector->syntheticTick;
  }
  total_FT_inc_TICK += column.syntheticTick;  // Was .pseudoTick
  if (EdWt_Mode == EdWt_Mode_DSE) {
    total_FT_inc_TICK += DSE_Vector->syntheticTick;  // Was .pseudoTick
  }
}

void HDual::updatePivots() {
  // UPDATE
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  //
  // Update the sets of indices of basic and nonbasic variables
  model->updatePivots(columnIn, rowOut, sourceOut);
  //  model->checkDualObjectiveValue("After  model->updatePivots");
  //
  // Update the iteration count and store the basis change if HiGHSDEV
  // is defined
  model->recordPivots(columnIn, columnOut, alpha);
  //
  // Update the invertible representation of the basis matrix
  model->updateFactor(&column, &row_ep, &rowOut, &invertHint);
  //
  // Update the row-wise representation of the nonbasic columns
  model->updateMatrix(columnIn, columnOut);
  //
  // Delete Freelist entry for columnIn
  dualRow.delete_Freelist(columnIn);
  //
  // Update the primal value for the row where the basis change has
  // occurred, and set the corresponding squared primal infeasibility
  // value in dualRHS.workArray
  dualRHS.update_pivots(rowOut, model->getWorkValue()[columnIn] + thetaPrimal);
  // Determine whether to reinvert based on the synthetic clock
  bool reinvert_syntheticClock =
      total_syntheticTick >= factor->build_syntheticTick;
#ifdef HiGHSDEV
  //    bool reinvert_syntheticClock = total_fake >=
  //    factor->build_syntheticTick;
#endif
  if (reinvert_syntheticClock && model->countUpdate >= 50) {
    invertHint = invertHint_syntheticClockSaysInvert;
  }
}

void HDual::iz_dvx_fwk() {
  // Initialise the Devex framework: reference set is all basic
  // variables
  model->timer.recordStart(HTICK_DEVEX_IZ);
  const int *NonbasicFlag = model->getNonbasicFlag();
  for (int vr_n = 0; vr_n < numTot; vr_n++) {
    //      if (model->getNonbasicFlag()[vr_n])
    //      if (NonbasicFlag[vr_n])
    //			Nonbasic variables not in reference set
    //	dvx_ix[vr_n] = dvx_not_in_R;
    //      else
    //			Basic variables in reference set
    //	dvx_ix[vr_n] = dvx_in_R;
    //
    // Assume all nonbasic variables have |NonbasicFlag[vr_n]|=1
    //
    // NonbasicFlag[vr_n]*NonbasicFlag[vr_n] evaluates faster than
    // abs(NonbasicFlag[vr_n])
    //
    // int dvx_ix_o_vr = 1-abs(NonbasicFlag[vr_n]);
    //
    int dvx_ix_o_vr = 1 - NonbasicFlag[vr_n] * NonbasicFlag[vr_n];
    dvx_ix[vr_n] = dvx_ix_o_vr;
  }
  dualRHS.workEdWt.assign(numRow, 1.0);  // Set all initial weights to 1
  n_dvx_it = 0;    // Zero the count of iterations with this Devex framework
  n_dvx_fwk += 1;  // Increment the number of Devex frameworks
  nw_dvx_fwk =
      false;  // Indicate that there's no need for a new Devex framework
  model->timer.recordFinish(HTICK_DEVEX_IZ);
}

void HDual::setCrash(const char *Crash_ArgV) {
  //	cout << "HDual::setCrash Crash_ArgV = " << Crash_ArgV << endl;
  if (strcmp(Crash_ArgV, "Off") == 0)
    Crash_Mode = Crash_Mode_No;
  else if (strcmp(Crash_ArgV, "LTSSF") == 0)
    Crash_Mode = Crash_Mode_Df;
  else if (strcmp(Crash_ArgV, "LTSSF1") == 0)
    Crash_Mode = Crash_Mode_LTSSF_k;
  else if (strcmp(Crash_ArgV, "LTSSF2") == 0)
    Crash_Mode = Crash_Mode_LTSSF_pri;
  else if (strcmp(Crash_ArgV, "LTSSF3") == 0)
    Crash_Mode = Crash_Mode_LTSF_k;
  else if (strcmp(Crash_ArgV, "LTSSF4") == 0)
    Crash_Mode = Crash_Mode_LTSF_pri;
  else if (strcmp(Crash_ArgV, "LTSSF5") == 0)
    Crash_Mode = Crash_Mode_LTSF;
  else if (strcmp(Crash_ArgV, "LTSSF6") == 0)
    Crash_Mode = Crash_Mode_Bixby;
  else if (strcmp(Crash_ArgV, "LTSSF7") == 0)
    Crash_Mode = Crash_Mode_BixbyNoNzCCo;
  else if (strcmp(Crash_ArgV, "Bs") == 0)
    Crash_Mode = Crash_Mode_Bs;
#ifdef HiGHSDEV
  else if (strcmp(Crash_ArgV, "TsSing") == 0)
    Crash_Mode = Crash_Mode_TsSing;
#endif
  else {
    cout << "HDual::setCrash unrecognised CrashArgV = " << Crash_ArgV
         << " - using No crash" << endl;
    Crash_Mode = Crash_Mode_No;
  }
  //	if (Crash_Mode == Crash_Mode_LTSSF) {
  //		crash.ltssf_iz_mode(Crash_Mode);
  //	}
  //		cout<<"HDual::setCrash Crash_Mode = " << Crash_Mode << endl;
}

void HDual::setPrice(const char *Price_ArgV) {
  //  cout<<"HDual::setPrice Price_ArgV = "<<Price_ArgV<<endl;
  alw_price_by_col_sw = false;
  alw_price_by_row_sw = false;
  alw_price_ultra = false;
  if (strcmp(Price_ArgV, "Col") == 0) {
    Price_Mode = Price_Mode_Col;
  } else if (strcmp(Price_ArgV, "Row") == 0) {
    Price_Mode = Price_Mode_Row;
  } else if (strcmp(Price_ArgV, "RowSw") == 0) {
    Price_Mode = Price_Mode_Row;
    alw_price_by_row_sw = true;
  } else if (strcmp(Price_ArgV, "RowSwColSw") == 0) {
    Price_Mode = Price_Mode_Row;
    alw_price_by_col_sw = true;
    alw_price_by_row_sw = true;
  } else if (strcmp(Price_ArgV, "RowUltra") == 0) {
    Price_Mode = Price_Mode_Row;
    alw_price_by_col_sw = true;
    alw_price_by_row_sw = true;
    alw_price_ultra = true;
  } else {
    cout << "HDual::setPrice unrecognised PriceArgV = " << Price_ArgV
         << " - using row Price with switch or colump price switch" << endl;
    Price_Mode = Price_Mode_Row;
    alw_price_by_col_sw = true;
    alw_price_by_row_sw = true;
  }
}

void HDual::setEdWt(const char *EdWt_ArgV) {
  //	cout<<"HDual::setEdWt EdWt_ArgV = "<<EdWt_ArgV<<endl;
  if (strcmp(EdWt_ArgV, "Dan") == 0)
    EdWt_Mode = EdWt_Mode_Dan;
  else if (strcmp(EdWt_ArgV, "Dvx") == 0)
    EdWt_Mode = EdWt_Mode_Dvx;
  else if (strcmp(EdWt_ArgV, "DSE") == 0) {
    EdWt_Mode = EdWt_Mode_DSE;
    iz_DSE_wt = true;
    alw_DSE2Dvx_sw = false;
  } else if (strcmp(EdWt_ArgV, "DSE0") == 0) {
    EdWt_Mode = EdWt_Mode_DSE;
    iz_DSE_wt = false;
    alw_DSE2Dvx_sw = false;
  } else if (strcmp(EdWt_ArgV, "DSE2Dvx") == 0) {
    EdWt_Mode = EdWt_Mode_DSE;
    iz_DSE_wt = true;
    alw_DSE2Dvx_sw = true;
  } else {
    cout << "HDual::setEdWt unrecognised EdWtArgV = " << EdWt_ArgV
         << " - using DSE with possible switch to Devex" << endl;
    EdWt_Mode = EdWt_Mode_DSE;
    iz_DSE_wt = true;
    alw_DSE2Dvx_sw = true;
  }
  //	cout<<"HDual::setEdWt iz_DSE_wt = " << iz_DSE_wt << endl;
}

void HDual::setTimeLimit(double TimeLimit_ArgV) {
  //	cout<<"HDual::setTimeLimit TimeLimit_ArgV = "<<TimeLimit_ArgV<<endl;
  TimeLimitValue = TimeLimit_ArgV;
}

void HDual::setPresolve(const char *Presolve_ArgV) {
  //	cout<<"HDual::setPresolve Presolve_ArgV = "<<Presolve_ArgV<<endl;
  if (strcmp(Presolve_ArgV, "Off") == 0)
    Presolve_Mode = Presolve_Mode_Off;
  else if (strcmp(Presolve_ArgV, "On") == 0)
    Presolve_Mode = Presolve_Mode_On;
  else {
    cout << "HDual::setPresolve unrecognised PresolveArgV = " << Presolve_ArgV
         << " - setting presolve off" << endl;
    Presolve_Mode = Presolve_Mode_Off;
  }
}

// Utility to get a row of the inverse of B for SCIP
int HDual::util_getBasisInvRow(int r, double *coef, int *inds, int *ninds) {
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = r;
  row_ep.array[r] = 1;
  row_ep.packFlag = true;
  factor->btran(row_ep, row_epDensity);
  //  printf("util_getBasisInvRow: nnz = %4d/%4d\n", row_ep.count, numRow);
  for (int row = 0; row < numRow; row++) {
    //    printf("BasisInvRow(%4d) = %11g\n", row,  row_ep.array[row]);
    coef[row] = row_ep.array[row];
  }
  if (0 <= row_ep.count && row_ep.count <= numRow) {
    for (int ix = 0; ix < row_ep.count; ix++) inds[ix] = row_ep.index[ix];
    ninds[0] = row_ep.count;
  } else {
    printf(
        "util_getBasisInvRow: row_ep.count < 0 or row_ep.count > numRow: %4d; "
        "%4d\n",
        row_ep.count, numRow);
    ninds[0] = -1;
  }
  cout << flush;
  return 0;
}

double HDual::an_bs_cond(HModel *ptr_model) {
  model = ptr_model;
  // Alias to the matrix
  matrix = model->getMatrix();
  const int *Astart = matrix->getAstart();
  const double *Avalue = matrix->getAvalue();
  // Compute the Hager condition number estimate for the basis matrix
  double NoDensity = 1;
  bs_cond_x.resize(numRow);
  bs_cond_y.resize(numRow);
  bs_cond_z.resize(numRow);
  bs_cond_w.resize(numRow);
  // x = ones(n,1)/n;
  // y = A\x;
  double mu = 1.0 / numRow;
  double norm_Binv;
  for (int r_n = 0; r_n < numRow; r_n++) bs_cond_x[r_n] = mu;
  row_ep.clear();
  row_ep.count = numRow;
  for (int r_n = 0; r_n < numRow; r_n++) {
    row_ep.index[r_n] = r_n;
    row_ep.array[r_n] = bs_cond_x[r_n];
  }
  for (int ps_n = 1; ps_n <= 5; ps_n++) {
    row_ep.packFlag = false;
    factor->ftran(row_ep, NoDensity);
    // zeta = sign(y);
    for (int r_n = 0; r_n < numRow; r_n++) {
      bs_cond_y[r_n] = row_ep.array[r_n];
      if (bs_cond_y[r_n] > 0)
        bs_cond_w[r_n] = 1.0;
      else if (bs_cond_y[r_n] < 0)
        bs_cond_w[r_n] = -1.0;
      else
        bs_cond_w[r_n] = 0.0;
    }
    // z=A'\zeta;
    row_ep.clear();
    row_ep.count = numRow;
    for (int r_n = 0; r_n < numRow; r_n++) {
      row_ep.index[r_n] = r_n;
      row_ep.array[r_n] = bs_cond_w[r_n];
    }
    row_ep.packFlag = false;
    factor->btran(row_ep, NoDensity);
    // norm_z = norm(z,'inf');
    // ztx = z'*x ;
    // NormEst = norm(y,1);
    // fd_i = 0;
    // for i=1:n
    //    if abs(z(i)) == norm_z
    //        fd_i = i;
    //        break
    //    end
    // end
    double norm_z = 0.0;
    double ztx = 0.0;
    norm_Binv = 0.0;
    int argmax_z = -1;
    for (int r_n = 0; r_n < numRow; r_n++) {
      bs_cond_z[r_n] = row_ep.array[r_n];
      double abs_z_v = abs(bs_cond_z[r_n]);
      if (abs_z_v > norm_z) {
        norm_z = abs_z_v;
        argmax_z = r_n;
      }
      ztx += bs_cond_z[r_n] * bs_cond_x[r_n];
      norm_Binv += abs(bs_cond_y[r_n]);
    }
    // printf("%2d: ||z||_inf = %8.2g; z^T*x = %8.2g; ||y||_1 = %g\n", ps_n,
    // norm_z, ztx, norm_Binv);
    if (norm_z <= ztx) break;
    // x = zeros(n,1);
    // x(fd_i) = 1;
    for (int r_n = 0; r_n < numRow; r_n++) bs_cond_x[r_n] = 0.0;
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = argmax_z;
    row_ep.array[argmax_z] = 1.0;
    bs_cond_x[argmax_z] = 1.0;
  }
  double norm_B = 0.0;
  for (int r_n = 0; r_n < numRow; r_n++) {
    int vr_n = model->getBaseIndex()[r_n];
    double c_norm = 0.0;
    if (vr_n < numCol)
      for (int el_n = Astart[vr_n]; el_n < Astart[vr_n + 1]; el_n++)
        c_norm += abs(Avalue[el_n]);
    else
      c_norm += 1.0;
    norm_B = max(c_norm, norm_B);
  }
  double cond_B = norm_Binv * norm_B;
  //  printf("Hager estimate of ||B^{-1}||_1 = %g; ||B||_1 = %g so cond_1(B)
  //  estimate is %g\n", norm_Binv, norm_B, cond_B);
  return cond_B;
}

#ifdef HiGHSDEV
void HDual::iterateRpIterDa(bool header) {
  if (header) {
    printf(
        " Inv       NumCk     LvR     LvC     EnC        DlPr        ThDu      "
        "  ThPr          Aa");
  } else {
    printf(" %3d %11.4g %7d %7d %7d %11.4g %11.4g %11.4g %11.4g", invertHint,
           numericalTrouble, rowOut, columnOut, columnIn, deltaPrimal,
           thetaDual, thetaPrimal, alpha);
  }
}

void HDual::iterateRpDsty(bool header) {
  if (header) {
    printf("  Col R_Ep R_Ap  DSE");
  } else {
    int l10ColDse = intLog10(columnDensity);
    int l10REpDse = intLog10(row_epDensity);
    int l10RapDse = intLog10(row_apDensity);
    double lc_rowdseDensity;
    if (EdWt_Mode == EdWt_Mode_DSE) {
      lc_rowdseDensity = rowdseDensity;
    } else {
      lc_rowdseDensity = 0;
    }
    int l10DseDse = intLog10(lc_rowdseDensity);
    printf(" %4d %4d %4d %4d", l10ColDse, l10REpDse, l10RapDse, l10DseDse);
  }
}

int HDual::intLog10(double v) {
  int intLog10V = -99;
  if (v > 0) intLog10V = log(v) / log(10.0);
  return intLog10V;
}

void HDual::iterateOpRecBf(int opTy, HVector &vector, double hist_dsty) {
  AnIterOpRec *AnIter = &AnIterOp[opTy];
  AnIter->AnIterOpNumCa++;
  double curr_dsty = 1.0 * vector.count / numRow;
  //  printf("%10s: %g<= %g;  %g<= %g\n", AnIter->AnIterOpName.c_str(),
  //	 curr_dsty, AnIter->AnIterOpHyperCANCEL,
  //	 hist_dsty, AnIter->AnIterOpHyperTRAN);
  if (curr_dsty <= AnIter->AnIterOpHyperCANCEL &&
      hist_dsty <= AnIter->AnIterOpHyperTRAN)
    AnIter->AnIterOpNumHyperOp++;
}

void HDual::iterateOpRecAf(int opTy, HVector &vector) {
  AnIterOpRec *AnIter = &AnIterOp[opTy];
  double rsDsty = 1.0 * vector.count / AnIter->AnIterOpRsDim;
  if (rsDsty <= hyperRESULT) AnIter->AnIterOpNumHyperRs++;
  AnIter->AnIterOpRsMxNNZ = max(vector.count, AnIter->AnIterOpRsMxNNZ);
  if (rsDsty > 0) {
    AnIter->AnIterOpLog10RsDsty += log(rsDsty) / log(10.0);
  } else {
    double vectorNorm = 0;
    for (int index = 0; index < AnIter->AnIterOpRsDim; index++) {
      double vectorValue = vector.array[index];
      vectorNorm += vectorValue * vectorValue;
    }
    vectorNorm = sqrt(vectorNorm);
    printf("Strange: operation %s has result density = %g: ||vector|| = %g\n",
           AnIter->AnIterOpName.c_str(), rsDsty, vectorNorm);
  }
}

void HDual::iterateRpAn() {
  int AnIterNumIter = model->numberIteration - AnIterIt0;
  printf("\nAnalysis of %d iterations (%d to %d)\n", AnIterNumIter,
         AnIterIt0 + 1, model->numberIteration);
  if (AnIterNumIter <= 0) return;
  int lc_EdWtNumIter;
  lc_EdWtNumIter = AnIterNumEdWtIt[EdWt_Mode_DSE];
  if (lc_EdWtNumIter > 0)
    printf("DSE for %7d (%3d%%) iterations\n", lc_EdWtNumIter,
           (100 * lc_EdWtNumIter) / AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[EdWt_Mode_Dvx];
  if (lc_EdWtNumIter > 0)
    printf("Dvx for %7d (%3d%%) iterations\n", lc_EdWtNumIter,
           (100 * lc_EdWtNumIter) / AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[EdWt_Mode_Dan];
  if (lc_EdWtNumIter > 0)
    printf("Dan for %7d (%3d%%) iterations\n", lc_EdWtNumIter,
           (100 * lc_EdWtNumIter) / AnIterNumIter);
  printf("\n");
  for (int k = 0; k < NumAnIterOpTy; k++) {
    AnIterOpRec *AnIter = &AnIterOp[k];
    int lcNumCa = AnIter->AnIterOpSuNumCa;
    printf("\n%-9s performed %7d times\n", AnIter->AnIterOpName.c_str(),
           AnIter->AnIterOpSuNumCa);
    if (lcNumCa > 0) {
      int lcHyperOp = AnIter->AnIterOpSuNumHyperOp;
      int lcHyperRs = AnIter->AnIterOpSuNumHyperRs;
      int pctHyperOp = (100 * lcHyperOp) / lcNumCa;
      int pctHyperRs = (100 * lcHyperRs) / lcNumCa;
      double lcRsDsty = pow(10.0, AnIter->AnIterOpSuLog10RsDsty / lcNumCa);
      int lcNumNNz = lcRsDsty * AnIter->AnIterOpRsDim;
      int lcMxNNz = AnIter->AnIterOpRsMxNNZ;
      double lcMxNNzDsty = (1.0 * lcMxNNz) / AnIter->AnIterOpRsDim;
      printf("   %11d hyper-sparse operations (%3d%%)\n", lcHyperOp,
             pctHyperOp);
      printf("   %11d hyper-sparse results    (%3d%%)\n", lcHyperRs,
             pctHyperRs);
      printf("   %11.4g density of result (%d nonzeros)\n", lcRsDsty, lcNumNNz);
      printf("   %11.4g density of result with max (%d) nonzeros\n",
             lcMxNNzDsty, lcMxNNz);
    }
  }
  int NumInvert = 0;
  for (int k = 1; k <= AnIterNumInvertHint; k++)
    NumInvert += AnIterNumInvert[k];
  if (NumInvert > 0) {
    int lcNumInvert = 0;
    printf("\nInvert    performed %7d times: average frequency = %d\n",
           NumInvert, AnIterNumIter / NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_updateLimitReached];
    if (lcNumInvert > 0)
      printf("%7d (%3d%%) Invert operations due to update limit reached\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_syntheticClockSaysInvert];
    if (lcNumInvert > 0)
      printf("%7d (%3d%%) Invert operations due to pseudo-clock\n", lcNumInvert,
             (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_possiblyOptimal];
    if (lcNumInvert > 0)
      printf("%7d (%3d%%) Invert operations due to possibly optimal\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_possiblyPrimalUnbounded];
    if (lcNumInvert > 0)
      printf("%7d (%3d%%) Invert operations due to possibly primal unbounded\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_possiblyDualUnbounded];
    if (lcNumInvert > 0)
      printf("%7d (%3d%%) Invert operations due to possibly dual unbounded\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_possiblySingularBasis];
    if (lcNumInvert > 0)
      printf("%7d (%3d%%) Invert operations due to possibly singular basis\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_primalInfeasibleInPrimalSimplex];
    if (lcNumInvert > 0)
      printf(
          "%7d (%3d%%) Invert operations due to primal infeasible in primal "
          "simplex\n",
          lcNumInvert, (100 * lcNumInvert) / NumInvert);
  }
  printf("\n%7d (%3d%%) primal degenerate iterations\n", AnIterNumPrDgnIt,
         (100 * AnIterNumPrDgnIt) / AnIterNumIter);
  printf("%7d (%3d%%)   dual degenerate iterations\n", AnIterNumDuDgnIt,
         (100 * AnIterNumDuDgnIt) / AnIterNumIter);
  int suPrice = AnIterNumColPrice + AnIterNumRowPrice + AnIterNumRowPriceWSw +
                AnIterNumRowPriceUltra;
  if (suPrice > 0) {
    printf("\n%7d Price operations:\n", suPrice);
    printf("%7d Col Price      (%3d%%)\n", AnIterNumColPrice,
           (100 * AnIterNumColPrice) / suPrice);
    printf("%7d Row Price      (%3d%%)\n", AnIterNumRowPrice,
           (100 * AnIterNumRowPrice) / suPrice);
    printf("%7d Row PriceWSw   (%3d%%)\n", AnIterNumRowPriceWSw,
           (100 * AnIterNumRowPriceWSw / suPrice));
    printf("%7d Row PriceUltra (%3d%%)\n", AnIterNumRowPriceUltra,
           (100 * AnIterNumRowPriceUltra / suPrice));
  }
  printf("\n%7d (%3d%%) costly DSE        iterations\n", AnIterNumCostlyDseIt,
         (100 * AnIterNumCostlyDseIt) / AnIterNumIter);

  //
  // Add a record for the final iterations: may end up with one more
  // than AnIterTraceMxNumRec records, so ensure that there is enough
  // space in the arrays
  //
  AnIterTraceNumRec++;
  AnIterTraceRec *lcAnIter;
  lcAnIter = &AnIterTrace[AnIterTraceNumRec];
  lcAnIter->AnIterTraceIter = model->numberIteration;
  lcAnIter->AnIterTraceTime = model->timer.getTime();
  lcAnIter->AnIterTraceDsty[AnIterOpTy_Btran] = row_epDensity;
  lcAnIter->AnIterTraceDsty[AnIterOpTy_Price] = row_apDensity;
  lcAnIter->AnIterTraceDsty[AnIterOpTy_Ftran] = columnDensity;
  lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranBFRT] = columnDensity;
  if (EdWt_Mode == EdWt_Mode_DSE) {
    lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE] = rowdseDensity;
    lcAnIter->AnIterTraceAux0 = AnIterCostlyDseMeasure;
  } else {
    lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE] = 0;
    lcAnIter->AnIterTraceAux0 = 0;
  }
  lcAnIter->AnIterTraceEdWt_Mode = EdWt_Mode;

  if (AnIterTraceIterDl >= 100) {
    printf("\n Iteration speed analysis\n");
    lcAnIter = &AnIterTrace[0];
    int fmIter = lcAnIter->AnIterTraceIter;
    double fmTime = lcAnIter->AnIterTraceTime;
    printf(
        "   Iter ( FmIter: ToIter)     Time Iter/sec |  Col R_Ep R_Ap  DSE | "
        "EdWt | Aux0\n");
    for (int rec = 1; rec <= AnIterTraceNumRec; rec++) {
      lcAnIter = &AnIterTrace[rec];
      int toIter = lcAnIter->AnIterTraceIter;
      double toTime = lcAnIter->AnIterTraceTime;
      int dlIter = toIter - fmIter;
      if (rec < AnIterTraceNumRec && dlIter != AnIterTraceIterDl)
        printf("STRANGE: %d = dlIter != AnIterTraceIterDl = %d\n", dlIter,
               AnIterTraceIterDl);
      double dlTime = toTime - fmTime;
      int iterSpeed = 0;
      if (dlTime > 0) iterSpeed = dlIter / dlTime;
      int lcEdWt_Mode = lcAnIter->AnIterTraceEdWt_Mode;
      int l10ColDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_Ftran]);
      int l10REpDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_Btran]);
      int l10RapDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_Price]);
      int l10DseDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE]);
      int l10Aux0 = intLog10(lcAnIter->AnIterTraceAux0);
      string strEdWt_Mode;
      if (lcEdWt_Mode == EdWt_Mode_DSE)
        strEdWt_Mode = "DSE";
      else if (lcEdWt_Mode == EdWt_Mode_Dvx)
        strEdWt_Mode = "Dvx";
      else if (lcEdWt_Mode == EdWt_Mode_Dan)
        strEdWt_Mode = "Dan";
      else
        strEdWt_Mode = "XXX";
      printf("%7d (%7d:%7d) %8.4f  %7d | %4d %4d %4d %4d |  %3s | %4d\n",
             dlIter, fmIter, toIter, dlTime, iterSpeed, l10ColDse, l10REpDse,
             l10RapDse, l10DseDse, strEdWt_Mode.c_str(), l10Aux0);
      fmIter = toIter;
      fmTime = toTime;
    }
    printf("\n");
  }
}

// TODO Put this in the right place - try to identify when workValue is set
// [nonbasic primals are set to the right bound for dual feasibility?]
void HDual::an_iz_vr_v() {
  double norm_bc_pr_vr = 0;
  double norm_bc_du_vr = 0;
  for (int r_n = 0; r_n < numRow; r_n++) {
    int vr_n = model->getBaseIndex()[r_n];
    norm_bc_pr_vr += baseValue[r_n] * baseValue[r_n];
    norm_bc_du_vr += workDual[vr_n] * workDual[vr_n];
  }
  double norm_nonbc_pr_vr = 0;
  double norm_nonbc_du_vr = 0;
  for (int vr_n = 0; vr_n < numTot; vr_n++) {
    if (model->getNonbasicFlag()[vr_n]) {
      double pr_act_v = model->getWorkValue()[vr_n];
      norm_nonbc_pr_vr += pr_act_v * pr_act_v;
      norm_nonbc_du_vr += workDual[vr_n] * workDual[vr_n];
    }
  }
  // printf("Initial point has %d dual infeasibilities\n", dualInfeasCount);
  // printf("Norm of the basic    primal activites is %g\n",
  // sqrt(norm_bc_pr_vr)); printf("Norm of the basic    dual   activites is
  // %g\n", sqrt(norm_bc_du_vr)); printf("Norm of the nonbasic primal activites
  // is %g\n", sqrt(norm_nonbc_pr_vr)); printf("Norm of the nonbasic dual
  // activites is %g\n", sqrt(norm_nonbc_du_vr));
}
#endif
