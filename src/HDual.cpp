#include "HDual.h"
#include "HConst.h"
#include "HTimer.h"
#include "HPrimal.h"
#include "HCrash.h"
#include "HMatrix.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>
#include <set>
#include <stdexcept>
using namespace std;

void HDual::solve(HModel *ptr_model, int variant, int num_threads)
{
  //  printf("\nEntering solve(HModel *ptr_model, int variant, int num_threads)\n");cout<<flush;
  assert(ptr_model != NULL);
  dual_variant = variant;
  model = ptr_model;
#ifdef HiGHSDEV
  //  printf("model->mlFg_Report() 1\n"); cout << flush; model->mlFg_Report(); cout << flush;
#endif
  // Setup two work buffers in model required for solve()
  model->buffer.setup(model->numRow);
  model->bufferLong.setup(model->numCol);
  
  //  printf("model->mlFg_haveEdWt 0 = %d\n", model->mlFg_haveEdWt);cout<<flush;
  
  // Setup aspects of the model data which are needed for solve() but better left until now for efficiency reasons.
  //  printf("HDual::solve - Calling model->setup_for_solve()\n");cout<<flush;
  model->setup_for_solve();
  
  //  printf("model->mlFg_haveEdWt 1 = %d\n", model->mlFg_haveEdWt);cout<<flush;
  
  model->problemStatus = LP_Status_Unset;
  model->numberIteration = 0;
#ifdef HiGHSDEV
  //Initialise numbers and times of rebuilds and inverts.
  totalRebuilds = 0;
  totalRebuildTime = 0;
  model->totalInverts = 0;
  model->totalInvertTime = 0;
#endif
  //Cannot solve box-constrained LPs
  if (model->numRow == 0)
    return;
#ifdef HiGHSDEV
  //  printf("model->mlFg_Report() 2\n"); cout << flush; model->mlFg_Report(); cout << flush;
#endif
  
  model->timer.reset();
  
  n_ph1_du_it = 0;
  n_ph2_du_it = 0;
  n_pr_it = 0;
  //Set SolveBailout to be true if control is to be returned immediately to calling function
  SolveBailout = false;
  if (TimeLimitValue == 0)
    {
      TimeLimitValue = 1000000.0;
#ifdef HiGHSDEV
      printf("Setting TimeLimitValue = %g\n", TimeLimitValue);
#endif
    }
  
  // Initialise working environment
  //Does LOTS, including initialisation of edge weights. Should only
  //be called if model dimension changes
  init(num_threads);
  
  model->initCost(1);
  if (!model->mlFg_haveFreshInvert) {
    int rankDeficiency = model->computeFactor();
    if (rankDeficiency) {
      throw runtime_error("Dual initialise: singular-basis-matrix");
    }
#ifdef HiGHSDEV
    double bsCond = an_bs_cond(model);
    printf("Initial basis condition estimate of %11.4g is", bsCond);
    if (bsCond > 1e12) {
      printf(" excessive\n");
      return;
    } else {
      printf(" OK\n");
    }
#endif
  }
  //Consider initialising edge weights
  //
  //NB workEdWt is assigned and initialised to 1s in
  //dualRHS.setup(model) so that CHUZR is well defined, even for
  //Dantzig pricing
  //
#ifdef HiGHSDEV
  //  printf("model->mlFg_haveEdWt 2 = %d; EdWt_Mode = %d; EdWt_Mode_DSE = %d\n",
  //	 model->mlFg_haveEdWt, EdWt_Mode, EdWt_Mode_DSE);cout<<flush;
  //  printf("Edge weights known? %d\n", !model->mlFg_haveEdWt);cout<<flush;
#endif
  if (!model->mlFg_haveEdWt) {
    //Edge weights are not known
    //Set up edge weights according to EdWt_Mode and iz_DSE_wt
    if (EdWt_Mode == EdWt_Mode_Dvx) {
      //Using dual Devex edge weights
      //Zero the number of Devex frameworks used and set up the first one
      n_dvx_fwk = 0;
      dvx_ix.assign(numTot, 0);
      iz_dvx_fwk();}
    else if (EdWt_Mode == EdWt_Mode_DSE) {
      //Using dual steepest edge (DSE) weights
      int numBasicStructurals = numRow - model->numBasicLogicals;
#ifdef HiGHSDEV
      n_wg_DSE_wt = 0;
      printf("If (0<numBasicStructurals = %d) && %d = iz_DSE_wt: Compute exact DSE weights\n",
	     numBasicStructurals, iz_DSE_wt);
#endif
      if (numBasicStructurals > 0 && iz_DSE_wt) {
	//Basis is not logical and DSE weights are to be initialised
#ifdef HiGHSDEV
	printf("Compute exact DSE weights\n"); //int RpI = 1;
	double IzDseEdWtTT = model->timer.getTime();
#endif
	for (int i = 0; i < numRow; i++) {
#ifdef HiGHSDEV
	  //	  if (i==RpI) {printf("Computing exact DSE weight %d\n", i); RpI = RpI*2;}
#endif
	  row_ep.clear();
	  row_ep.count = 1;
	  row_ep.index[0] = i;
	  row_ep.array[i] = 1;
	  row_ep.packFlag = false;
	  factor->btran(row_ep, row_epDensity);
	  dualRHS.workEdWt[i] = row_ep.norm2();
	  double lc_OpRsDensity = 1.0 * row_ep.count / numRow;
	  row_epDensity = uOpRsDensityRec(lc_OpRsDensity, &row_epDensity, &row_epAvDensity, &row_epAvLog10Density); 
	  //row_epDensity = (1-runningAverageMu) * row_epDensity + runningAverageMu * (1.0 * row_ep.count / numRow);
	}
#ifdef HiGHSDEV
	IzDseEdWtTT = model->timer.getTime() - IzDseEdWtTT;
	printf("Computed %d initial DSE weights in %gs\n", numRow, IzDseEdWtTT);
	if (model->intOption[INTOPT_PRINT_FLAG]) printf("solve:: %d basic structurals: computed %d initial DSE weights in %gs, %d, %d, %g\n",
		 numBasicStructurals, numRow, IzDseEdWtTT, numBasicStructurals, numRow, IzDseEdWtTT);
#endif
      }
#ifdef HiGHSDEV
      else {
	if (model->intOption[INTOPT_PRINT_FLAG]) printf("solve:: %d basic structurals: starting from B=I so unit initial DSE weights\n", numBasicStructurals);
      }
#endif
    }
    //Indicate that edge weights are known
    model->mlFg_haveEdWt = 1;
  }
  
#ifdef HiGHSDEV
  //  printf("model->mlFg_haveEdWt 3 = %d\n", model->mlFg_haveEdWt);cout<<flush;
  bool rp_bs_cond = false;
  if (rp_bs_cond) {
    double bs_cond = an_bs_cond(ptr_model);
    printf("Initial basis condition estimate is %g\n", bs_cond );
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
      if (largeDual < myDual)
	largeDual = myDual;
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
  //  if ((solvePhase != 1) && (solvePhase != 2)) {printf("In solve(): solvePhase = %d\n", solvePhase);cout<<flush;}
#endif
  bool ok = model->OKtoSolve(1, solvePhase);
  if (!ok) {printf("NOT OK TO SOLVE???\n");cout<<flush;}
  //  assert(ok);
  
#ifdef HiGHSDEV
  //  rp_hsol_sol(ptr_model);
  //  Analyse the initial values of primal and dual variables
  //  an_iz_vr_v();
#endif
  
  // The major solving loop
  
  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiHGSDEV
  iterateIzAn();

  while (solvePhase) {
#ifdef HiGHSDEV
    int it0 = model->numberIteration;
    //    printf("HDual::solve Phase %d: Iteration %d; totalTime = %g; timer.getTime = %g\n",
    //	   solvePhase, model->numberIteration, model->totalTime, model->timer.getTime());cout<<flush;
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
    if (solvePhase == 4)
      break;
    // Possibly bail out
    if (SolveBailout)
      break;
  }
  
#ifdef HiGHSDEV
  if (AnIterLg) iterateRpAn();
  // Report the ticks before primal
  if (dual_variant == HDUAL_VARIANT_PLAIN) {
    int reportList[] = { HTICK_INVERT, HTICK_PERM_WT, HTICK_COMPUTE_DUAL, HTICK_CORRECT_DUAL, HTICK_COMPUTE_PRIMAL, HTICK_COLLECT_PR_IFS, HTICK_COMPUTE_DUOBJ, HTICK_REPORT_INVERT,
			 HTICK_CHUZR1, HTICK_BTRAN,
			 HTICK_PRICE, HTICK_CHUZC0, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3, HTICK_CHUZC4, HTICK_DEVEX_WT, 
			 HTICK_FTRAN, HTICK_FTRAN_BFRT, HTICK_FTRAN_DSE,
			 HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT, HTICK_DEVEX_IZ,
			 HTICK_UPDATE_PIVOTS, HTICK_UPDATE_FACTOR, HTICK_UPDATE_MATRIX };
    int reportCount = sizeof(reportList) / sizeof(int);
    model->timer.report(reportCount, reportList, 1.0);
    bool rpIterate = false;
    if (rpIterate) {
      int reportList[] = { HTICK_ITERATE };
      int reportCount = sizeof(reportList) / sizeof(int);
      model->timer.report(reportCount, reportList, 1.0);
    }
    if (rpIterate) {
      int reportList[] = {
	HTICK_ITERATE_REBUILD, HTICK_ITERATE_CHUZR, HTICK_ITERATE_CHUZC, HTICK_ITERATE_FTRAN, HTICK_ITERATE_VERIFY, HTICK_ITERATE_DUAL, HTICK_ITERATE_PRIMAL, HTICK_ITERATE_DEVEX_IZ, HTICK_ITERATE_PIVOTS
      };
      int reportCount = sizeof(reportList) / sizeof(int);
      model->timer.report(reportCount, reportList, 1.0);
    }
  }
  
  if (dual_variant == HDUAL_VARIANT_TASKS) {
    int reportList[] = { HTICK_INVERT, HTICK_CHUZR1, HTICK_BTRAN,
			 HTICK_PRICE, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3,
			 HTICK_DEVEX_WT, HTICK_FTRAN, HTICK_FTRAN_BFRT, HTICK_FTRAN_DSE,
			 HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,
			 HTICK_UPDATE_FACTOR, HTICK_GROUP1, HTICK_GROUP2 };
    int reportCount = sizeof(reportList) / sizeof(int);
    model->timer.report(reportCount, reportList, 1.0);
  }
  
  if (dual_variant == HDUAL_VARIANT_MULTI) {
    int reportList[] = { HTICK_INVERT, HTICK_CHUZR1, HTICK_BTRAN,
			 HTICK_PRICE, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3,
			 HTICK_DEVEX_WT, HTICK_FTRAN, HTICK_FTRAN_BFRT, HTICK_FTRAN_DSE,
			 HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,
			 HTICK_UPDATE_FACTOR, HTICK_UPDATE_ROW_EP };
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
      //Add in the count and time for any primal rebuilds
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
  //  rp_hsol_sol(ptr_model);
#endif
#ifdef HiGHSDEV
  //  if ((solvePhase != 1) && (solvePhase != 2)) {printf("In solve(): solvePhase = %d\n", solvePhase);cout<<flush;}
#endif
  ok = model->OKtoSolve(1, solvePhase);
  if (!ok) {printf("NOT OK After Solve???\n");cout<<flush;}
  //  assert(ok);
#ifdef HiGHSDEV
  //  printf("model->mlFg_Report() 9\n");cout<<flush; model->mlFg_Report();cout<<flush;
#endif

#ifdef HiGHSDEV
  if (model->anInvertTime) {
    printf("Time: Total inverts =  %4d; Total invert  time = %11.4g of Total time = %11.4g", model->totalInverts, model->totalInvertTime, model->totalTime);
    if (model->totalTime > 0.001) {
      printf(" (%6.2f%%)\n", (100*model->totalInvertTime)/model->totalTime);
    } else { 
      printf("\n");
    }
    cout << flush;
    printf("Time: Total rebuilds = %4d; Total rebuild time = %11.4g of Total time = %11.4g", totalRebuilds, totalRebuildTime, model->totalTime);
    if (model->totalTime > 0.001) {
      printf(" (%6.2f%%)\n", (100*totalRebuildTime)/model->totalTime);
    } else {
      printf("\n");
    }
    cout << flush;
  }
  model->util_anMlSol();
  
#endif  
}
 
void HDual::init(int num_threads)
{
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
  if (dual_variant == HDUAL_VARIANT_TASKS)
  {
    init_slice(num_threads - 2);
  }

  // Initialize for multi
  if (dual_variant == HDUAL_VARIANT_MULTI)
  {
    multi_num = num_threads;
    if (multi_num < 1)
      multi_num = 1;
    if (multi_num > HSOL_THREAD_LIMIT)
      multi_num = HSOL_THREAD_LIMIT;
    for (int i = 0; i < multi_num; i++)
    {
      multi_choice[i].row_ep.setup(numRow);
      multi_choice[i].column.setup(numRow);
      multi_choice[i].columnBFRT.setup(numRow);
    }
    init_slice(multi_num - 1);
  }
  multi_iteration = 0;
  string partitionFile = model->strOption[STROPT_PARTITION_FILE];
  if (partitionFile.size())
  {
    dualRHS.setup_partition(partitionFile.c_str());
  }
}

void HDual::init_slice(int init_sliced_num)
{
  // Number of slices
  slice_num = init_sliced_num;
  if (slice_num < 1)
    slice_num = 1;
  if (slice_num > HSOL_SLICED_LIMIT)
    slice_num = HSOL_SLICED_LIMIT;

  // Alias to the matrix
  const int *Astart = matrix->getAstart();
  const int *Aindex = matrix->getAindex();
  const double *Avalue = matrix->getAvalue();
  const int AcountX = Astart[numCol];

  // Figure out partition weight
  double sliced_countX = AcountX / slice_num;
  slice_start[0] = 0;
  for (int i = 0; i < slice_num - 1; i++)
  {
    int endColumn = slice_start[i] + 1; // At least one column
    int endX = Astart[endColumn];
    int stopX = (i + 1) * sliced_countX;
    while (endX < stopX)
    {
      endX = Astart[++endColumn];
    }
    slice_start[i + 1] = endColumn;
    if (endColumn >= numCol)
    {
      slice_num = i; // SHRINK
      break;
    }
  }
  slice_start[slice_num] = numCol;

  // Partition the matrix, row_ap and related packet
  vector<int> sliced_Astart;
  for (int i = 0; i < slice_num; i++)
  {
    // The matrix
    int mystart = slice_start[i];
    int mycount = slice_start[i + 1] - mystart;
    int mystartX = Astart[mystart];
    sliced_Astart.resize(mycount + 1);
    for (int k = 0; k <= mycount; k++)
      sliced_Astart[k] = Astart[k + mystart] - mystartX;
    //TODO generalise this call so slice can be used with non-logical initial basis
    slice_matrix[i].setup_lgBs(mycount, numRow, &sliced_Astart[0],
                               Aindex + mystartX, Avalue + mystartX);

    // The row_ap and its packages
    slice_row_ap[i].setup(mycount);
    slice_dualRow[i].setupSlice(model, mycount);
  }
}

void HDual::solve_phase1()
{
  model->util_reportMessage("dual-phase-1-start");
  // Switch to dual phase 1 bounds
  model->initBound(1);
  model->initValue();
  double lc_totalTime = model->totalTime + model->timer.getTime();
#ifdef HiGHSDEV
  // int lc_totalTime_rp_n = 0; printf("DualPh1: lc_totalTime = %5.2f; Record %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
  // Main solving structure
  model->timer.recordStart(HTICK_ITERATE);
  for (;;)
    {
      model->timer.recordStart(HTICK_ITERATE_REBUILD);
      rebuild();
      model->timer.recordFinish(HTICK_ITERATE_REBUILD);
      for (;;)
	{
	  switch (dual_variant)
	    {
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
	  if (invertHint)
	    break;
	  //printf("HDual::solve_phase1: Iter = %d; Objective = %g\n", model->numberIteration, model->objective);
	  /*
	  if (model->objective > model->dblOption[DBLOPT_OBJ_UB]) {
#ifdef SCIP_DEV
	    printf("HDual::solve_phase1: %g = Objective > dblOption[DBLOPT_OBJ_UB]\n", model->objective, model->dblOption[DBLOPT_OBJ_UB]);
#endif
	    model->problemStatus = LP_Status_ObjUB; 
	    break;
	  }
	  */
	}
      lc_totalTime = model->totalTime + model->timer.getTime();
#ifdef HiGHSDEV
      //      lc_totalTime_rp_n += 1; printf("DualPh1: lc_totalTime = %5.2f; Record %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
      if (lc_totalTime > TimeLimitValue)
	{
	  SolveBailout = true;
	  model->problemStatus = LP_Status_OutOfTime;
	  break;
	}
      // If the data are fresh from rebuild(), break out of
      // the outer loop to see what's ocurred
      // Was:	if (model->countUpdate == 0) break;
      if (model->mlFg_haveFreshRebuild)
	break;
    }
  
  model->timer.recordFinish(HTICK_ITERATE);
  if (SolveBailout)
    return;
  
  if (rowOut == -1)
  {
    model->util_reportMessage("dual-phase-1-optimal");
    // Go to phase 2
    if (model->objective == 0)
    {
      solvePhase = 2;
    }
    else
    {
      // We still have dual infeasible
      if (model->problemPerturbed)
      {
        // Clean up perturbation and go on
        cleanup();
        if (dualInfeasCount == 0)
          solvePhase = 2;
      }
      else
      {
        // Report dual infeasible
        solvePhase = -1;
        model->util_reportMessage("dual-infeasible");
        model->setProblemStatus(LP_Status_Unbounded);
      }
    }
  }
  else if (invertHint == invertHint_chooseColumnFail)
  {
    // chooseColumn has failed
    // Behave as "Report strange issues" below
    solvePhase = -1;
    model->util_reportMessage("dual-phase-1-not-solved");
    model->setProblemStatus(LP_Status_Failed);
  }
  else if (columnIn == -1)
  {
    // We got dual phase 1 unbounded - strange
    model->util_reportMessage("dual-phase-1-unbounded");
    if (model->problemPerturbed)
    {
      // Clean up perturbation and go on
      cleanup();
      if (dualInfeasCount == 0)
        solvePhase = 2;
    }
    else
    {
      // Report strange issues
      solvePhase = -1;
      model->util_reportMessage("dual-phase-1-not-solved");
      model->setProblemStatus(LP_Status_Failed);
    }
  }

  if (solvePhase == 2)
  {
    model->initBound();
    model->initValue();
  }
}

void HDual::solve_phase2()
{
  model->util_reportMessage("dual-phase-2-start");

  // Collect free variables
  dualRow.create_Freelist();
  double lc_totalTime = model->totalTime + model->timer.getTime();
#ifdef HiGHSDEV
  //  int lc_totalTime_rp_n = 0; printf("DualPh2: lc_totalTime = %5.2f; Record %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
  // Main solving structure
  model->timer.recordStart(HTICK_ITERATE);
  for (;;)
    {
      // Outer loop of solve_phase2()
      // Rebuild all values, reinverting B if updates have been performed
      model->timer.recordStart(HTICK_ITERATE_REBUILD);
      rebuild();
      model->timer.recordFinish(HTICK_ITERATE_REBUILD);
      if (dualInfeasCount > 0)
	break;
      for (;;)
	{
	  // Inner loop of solve_phase2()
	  // Performs one iteration in case HDUAL_VARIANT_PLAIN:
	  model->util_reportSolverProgress();
	  switch (dual_variant)
	    {
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
	  //invertHint can be true for various reasons see HModel.h
	  if (invertHint) break;
	  // Need the dual objective value to check for exceeding the
	  // upper bound, but can't afford to compute it every
	  // iteration!
	  //      model->computeDuObj();
#ifdef HiGHSDEV
	  //      double pr_obj_v = model->computePrObj();
	  //      printf("HDual::solve_phase2: Iter = %4d; Pr Obj = %.11g; Du Obj = %.11g\n",
	  //	     model->numberIteration, pr_obj_v, model->objective);
#endif
	  if (model->objective > model->dblOption[DBLOPT_OBJ_UB])
	    {
#ifdef SCIP_DEV
	      //printf("HDual::solve_phase2: Objective = %g > %g = dblOption[DBLOPT_OBJ_UB]\n", model->objective, model->dblOption[DBLOPT_OBJ_UB]);
#endif
	      model->problemStatus = LP_Status_ObjUB;
	      SolveBailout = true;
	      break;
	    }
	}
      lc_totalTime = model->totalTime + model->timer.getTime();
      if (model->problemStatus == LP_Status_ObjUB)
	{
	  SolveBailout = true;
	  break;
	}
#ifdef HiGHSDEV
      //      lc_totalTime_rp_n += 1; printf("DualPh2: lc_totalTime = %5.2f; Record %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
      if (lc_totalTime > TimeLimitValue)
	{
	  model->problemStatus = LP_Status_OutOfTime;
	  SolveBailout = true;
	  break;
	}
      // If the data are fresh from rebuild(), break out of
      // the outer loop to see what's ocurred
      // Was:	if (model->countUpdate == 0) break;
      if (model->mlFg_haveFreshRebuild)
	break;
    }
  model->timer.recordFinish(HTICK_ITERATE);
  
  if (SolveBailout) {
    return;
  }
  if (dualInfeasCount > 0)
    {
      // There are dual infeasiblities so switch to Phase 1 and return
      model->util_reportMessage("dual-phase-2-found-free");
      solvePhase = 1;
    }
  else if (rowOut == -1)
    {
      // There is no candidate in CHUZR, even after rebuild so probably optimal
      model->util_reportMessage("dual-phase-2-optimal");
      //		printf("Rebuild: cleanup()\n");
      cleanup();
      if (dualInfeasCount > 0)
	{
	  // There are dual infeasiblities after cleanup() so switch to primal simplex
	  solvePhase = 4; // Do primal
	}
      else
	{
	  // There are no dual infeasiblities after cleanup() so optimal!
	  solvePhase = 0;
	  model->util_reportMessage("problem-optimal");
	  model->setProblemStatus(LP_Status_Optimal);
	}
    }
  else if (invertHint == invertHint_chooseColumnFail)
    {
      // chooseColumn has failed
      // Behave as "Report strange issues" below
      solvePhase = -1;
      model->util_reportMessage("dual-phase-2-not-solved");
      model->setProblemStatus(LP_Status_Failed);
    }
  else if (columnIn == -1)
    {
      // There is no candidate in CHUZC, so probably dual unbounded
      model->util_reportMessage("dual-phase-2-unbounded");
      if (model->problemPerturbed)
	{
	  //If the costs have been perturbed, clean up and return
	  cleanup();
	}
      else
	{
	  //If the costs have not been perturbed, so dual unbounded---and hence primal infeasible
	  solvePhase = -1;
	  model->util_reportMessage("problem-infeasible");
	  model->setProblemStatus(LP_Status_Infeasible);
	}
    }
}

void HDual::rebuild()
{
  // Save history information
  model->recordPivots(-1, -1, 0); // Indicate REINVERT
#ifdef HiGHSDEV
  double tt0 = 0;
  if (anRebuildTime) tt0 = model->timer.getTime();
#endif
  int sv_invertHint = invertHint;
  invertHint = invertHint_no; // Was 0

  // Possibly Rebuild model->factor
  bool reInvert = model->countUpdate > 0;
  if (!model->InvertIfRowOutNeg) {
    // Don't reinvert if rowOut is negative [equivalently, if sv_invertHint == invertHint_possiblyOptimal]
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
    int rankDeficiency = model->computeFactor();
    model->timer.recordFinish(HTICK_INVERT);

    if (rankDeficiency) throw runtime_error("Dual reInvert: singular-basis-matrix");
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
  model->computeDuObj(solvePhase);
  model->timer.recordFinish(HTICK_COMPUTE_DUOBJ);
  //	model->util_reportNumberIterationObjectiveValue(sv_invertHint);

  model->timer.recordStart(HTICK_REPORT_INVERT);
  iterateRpInvert(sv_invertHint);
  model->timer.recordFinish(HTICK_REPORT_INVERT);
  
  total_INVERT_TICK = factor->pseudoTick;
  total_FT_inc_TICK = 0;
  total_fake = 0;
  
#ifdef HiGHSDEV
  if (anRebuildTime) {
    double rebuildTime = model->timer.getTime()-tt0;
    totalRebuilds++;
    totalRebuildTime += rebuildTime;
    printf("Dual  Ph%-2d rebuild %4d (%1d) on iteration %9d: Rebuild time = %11.4g; Total rebuild time = %11.4g\n",
	   solvePhase, totalRebuilds, sv_invertHint, model->numberIteration, rebuildTime, totalRebuildTime);
  }
#endif
  //Data are fresh from rebuild
  model->mlFg_haveFreshRebuild = 1;
}

void HDual::cleanup() {
	// Remove perturbation and recompute the dual solution
	model->util_reportMessage("dual-cleanup-shift");
	model->initCost();
	model->initBound();
	model->computeDual();
	model->computeDuObj(solvePhase);
	//	model->util_reportNumberIterationObjectiveValue(-1);
	iterateRpInvert(-1);

  model->computeDualInfeasInPrimal(&dualInfeasCount);
}

void HDual::iterate()
{
  // This is the main teration loop for dual revised simplex. All the
  // methods have as their first line if (invertHint) return;, where
  // invertHint is, for example, set to 1 when CHUZR finds no
  // candidate. This causes a break from the inner loop of
  // solve_phase% and, hence, a call to rebuild()

//	Reporting:
//	hsol data structures [NonBcFg; DvxIx; BcIx; DvxV] in rp_hsol_da_str();
//	hsol pivotal row in rp_hsol_pv_r();
//	hsol simplex iteration [rowOut; coumnOut; columnIn] in rp_hsol_si_it();
//	hsol row-wise matrix after update in updateMatrix(columnIn, columnOut);
#ifdef HiGHSDEV
  if (rp_hsol) rp_hsol_da_str();
#endif
  model->timer.recordStart(HTICK_ITERATE_CHUZR);
  chooseRow();
  model->timer.recordFinish(HTICK_ITERATE_CHUZR);
  model->timer.recordStart(HTICK_ITERATE_CHUZC);
  chooseColumn(&row_ep);
  model->timer.recordFinish(HTICK_ITERATE_CHUZC);
#ifdef HiGHSDEV
  if (rp_hsol && rowOut >= 0) {printf("\nPvR: Row %2d\n", rowOut); dualRow.rp_hsol_pv_r();} cout<<flush;
#endif
  model->timer.recordStart(HTICK_ITERATE_FTRAN);
  updateFtranBFRT();
  // updateFtran(); computes the pivotal column in the data structure "column"
  updateFtran();
#ifdef HiGHSDEV
  if (rp_hsol) rp_hsol_pv_c(&column);
 #endif
 
  //updateFtranDSE performs the DSE FTRAN on pi_p
  if (EdWt_Mode == EdWt_Mode_DSE) updateFtranDSE(&row_ep);
  model->timer.recordFinish(HTICK_ITERATE_FTRAN);
  
  //updateVerify() Checks row-wise pivot against column-wise pivot for numerical trouble
  model->timer.recordStart(HTICK_ITERATE_VERIFY);
  updateVerify();
  model->timer.recordFinish(HTICK_ITERATE_VERIFY);
  
  //updateDual() Updates the dual values
  model->timer.recordStart(HTICK_ITERATE_DUAL);
  updateDual();
  model->timer.recordFinish(HTICK_ITERATE_DUAL);
  
  //updatePrimal(&row_ep); Updates the primal values and the edge weights
  model->timer.recordStart(HTICK_ITERATE_PRIMAL);
  updatePrimal(&row_ep);
  model->timer.recordFinish(HTICK_ITERATE_PRIMAL);
  
  if ((EdWt_Mode == EdWt_Mode_Dvx) && (nw_dvx_fwk)) {
    model->timer.recordStart(HTICK_ITERATE_DEVEX_IZ);
    iz_dvx_fwk();
    model->timer.recordFinish(HTICK_ITERATE_DEVEX_IZ);
  }
  
  //Update the basis representation
  model->timer.recordStart(HTICK_ITERATE_PIVOTS);
  updatePivots();
  model->timer.recordFinish(HTICK_ITERATE_PIVOTS);
  
  //Analyse the iteration: possibly report; possibly switch strategy
  iterateAn();
}

void HDual::iterate_tasks()
{
  slice_PRICE = 1;
  
  // Group 1
  chooseRow();
  
  // Disable slice when too sparse
  if (1.0 * row_ep.count / numRow < 0.01)
    slice_PRICE = 0;
  
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
  AnIterPrevRpNumCostlyDseIt = 0;
#ifdef HiGHSDEV
  AnIterPrevIt = 0;
  AnIterOpRec *AnIter;
  AnIter = &AnIterOp[AnIterOpTy_Btran]; AnIter->AnIterOpName = "Btran";
  AnIter = &AnIterOp[AnIterOpTy_Price]; AnIter->AnIterOpName = "Price";
  AnIter = &AnIterOp[AnIterOpTy_Ftran]; AnIter->AnIterOpName = "Ftran";
  AnIter = &AnIterOp[AnIterOpTy_FtranBFRT]; AnIter->AnIterOpName = "FtranBFRT";
  AnIter = &AnIterOp[AnIterOpTy_FtranDSE]; AnIter->AnIterOpName = "FtranDSE";
  for (int k=0; k<NumAnIterOpTy; k++) {
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
  for (int k=1; k<=AnIterNumInvertHint; k++) AnIterNumInvert[k]=0;
  AnIterNumPrDgnIt = 0;
  AnIterNumDuDgnIt = 0;
  AnIterNumColPrice = 0;
  AnIterNumRowPrice = 0;
  AnIterNumRowPriceWSw = 0;
  AnIterNumRowPriceUltra = 0;
  for (int k=0; k<=EdWt_Mode_Dan; k++) AnIterNumEdWtIt[k]=0;
  AnIterNumCostlyDseIt = 0;
  AnIterTraceNumRec = 0;
  AnIterTraceIterDl = 1;
  AnIterTraceRec *lcAnIter = &AnIterTrace[0];
  lcAnIter->AnIterTraceIter = AnIterIt0;
  lcAnIter->AnIterTraceTime = model->timer.getTime();
#endif
}

void HDual::iterateAn() {
  //Possibly report on the iteration
  iterateRp();

  // Possibly switch from DSE to Dvx
  if (EdWt_Mode == EdWt_Mode_DSE) {
    double AnIterCostlyDseMeasureDen;
    //    AnIterCostlyDseMeasureDen = row_epDensity*columnDensity;
    AnIterCostlyDseMeasureDen = max(max(row_epDensity, columnDensity), row_apDensity);
    if (AnIterCostlyDseMeasureDen > 0) {
      //      AnIterCostlyDseMeasure = rowdseDensity*rowdseDensity/AnIterCostlyDseMeasureDen;
      AnIterCostlyDseMeasure = rowdseDensity/AnIterCostlyDseMeasureDen;
      AnIterCostlyDseMeasure = AnIterCostlyDseMeasure*AnIterCostlyDseMeasure;
    } else {
      AnIterCostlyDseMeasure = 0;
    }
    bool CostlyDseIt = AnIterCostlyDseMeasure > AnIterCostlyDseMeasureLimit &&
      //      rowdseDensity*rowdseDensity > AnIterCostlyDseMeasureLimit*row_epDensity*columnDensity &&
      rowdseDensity > AnIterCostlyDseMnDensity;
    AnIterCostlyDseFq = (1-runningAverageMu)*AnIterCostlyDseFq;
    if (CostlyDseIt) {
      AnIterNumCostlyDseIt++;
      AnIterCostlyDseFq += runningAverageMu*1.0;
      int lcNumIter = model->numberIteration-AnIterIt0;
      if (alw_DSE2Dvx_sw
	  && (AnIterNumCostlyDseIt > lcNumIter*AnIterFracNumCostlyDseItbfSw)
	  && (lcNumIter > AnIterFracNumTot_ItBfSw*numTot)) {
        //At least 5% of the (at least) 0.1NumTot iterations have been costly DSE so switch to Devex
#ifdef HiGHSDEV
        printf("Switch from DSE to Dvx after %d costly DSE iterations of %d: Col_Dsty = %11.4g; R_Ep_Dsty = %11.4g; DSE_Dsty = %11.4g\n",
	       AnIterNumCostlyDseIt, lcNumIter, rowdseDensity, row_epDensity,  columnDensity);
#endif
        EdWt_Mode = EdWt_Mode_Dvx;
        //Zero the number of Devex frameworks used and set up the first one
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
      int lc_pct = (100*AnIterNumCostlyDseIt)/(AnIterCuIt-AnIterIt0);
      printf("| Fq = %4.2f; Su =%5d (%3d%%)", AnIterCostlyDseFq, AnIterNumCostlyDseIt, lc_pct);
      
      if (lc_NumCostlyDseIt>0) printf("; LcNum =%3d", lc_NumCostlyDseIt);
    }
    printf("\n");
  }
  
  for (int k=0; k<NumAnIterOpTy; k++) {
    AnIterOpRec *lcAnIterOp = &AnIterOp[k];
    if (lcAnIterOp->AnIterOpNumCa) {
      lcAnIterOp->AnIterOpSuNumCa      += lcAnIterOp->AnIterOpNumCa;
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
  if (AnIterCuIt > AnIterPrevIt) AnIterNumEdWtIt[EdWt_Mode] += (AnIterCuIt-AnIterPrevIt);
  
  AnIterTraceRec *lcAnIter = &AnIterTrace[AnIterTraceNumRec];
  //  if (model->numberIteration == AnIterTraceIterRec[AnIterTraceNumRec]+AnIterTraceIterDl) {
  if (model->numberIteration == lcAnIter->AnIterTraceIter+AnIterTraceIterDl) {
    if (AnIterTraceNumRec == AnIterTraceMxNumRec) {
      for (int rec = 1; rec<=AnIterTraceMxNumRec/2; rec++) AnIterTrace[rec] = AnIterTrace[2*rec];
      AnIterTraceNumRec = AnIterTraceNumRec/2;
      AnIterTraceIterDl = AnIterTraceIterDl*2;
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
      /*
	if (AnIterTraceIterDl > 10) {
	//	int reportList[] = { HTICK_INVERT, HTICK_CHUZR1, HTICK_BTRAN,
	//			     HTICK_PRICE, HTICK_CHUZC0, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3, HTICK_CHUZC4,
	//			     HTICK_DEVEX, HTICK_FTRAN, HTICK_FTRAN_DSE,
	//			     HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,
	//			     HTICK_UPDATE_FACTOR };
	//	int reportCount = sizeof(reportList) / sizeof(int);
	//	model->timer.report(reportCount, reportList, 1.0);
        iterateRpBrief(true);
        iterateRpBrief(false);
	}
      */
    }
  }
  AnIterPrevIt = AnIterCuIt;
#endif
}

void HDual::iterateRp() {
  if (model->intOption[INTOPT_PRINT_FLAG] != 4) return;
  int numIter = model->numberIteration;
  bool header= numIter % 10 == 1;
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
    printf(" FreeLsZ\n");
#endif
  } else {
    iterateRpIterPh(false);
    iterateRpDuObj(false);
#ifdef HiGHSDEV
    iterateRpIterDa(false);
    iterateRpDsty(false);
    printf(" %7d\n", dualRow.freeListSize);
#endif
  }
}

void HDual::iterateRpBrief(bool header) {
  if (header) {
    iterateRpIterPh(true);
    iterateRpDuObj(true);
#ifdef HiGHSDEV
    iterateRpDsty(true);
    printf(" FreeLsZ\n");
#endif
  } else {
    iterateRpIterPh(false);
    iterateRpDuObj(false);
#ifdef HiGHSDEV
    iterateRpDsty(false);
    printf(" %7d\n", dualRow.freeListSize);
#endif
  }
}

void HDual::iterateRpInvert(int i_v) {
  if (model->intOption[INTOPT_PRINT_FLAG] != 1 && model->intOption[INTOPT_PRINT_FLAG] != 4) return;
  printf("Iter %10d:", model->numberIteration);
#ifdef HiGHSDEV
  iterateRpDsty(true);
  iterateRpDsty(false);
#endif
  iterateRpDuObj(false);
  printf(" %2d", i_v);
  //  printf(" %10.0f; %6.4f; %6.4f",
  //	 factor->realTick, (factor->realTick)/(factor->fakeTick), factor->pseudoTick/(factor->fakeTick));
  printf("\n");
  
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
    model->computeDuObj(solvePhase);
    printf(" %20.10e", model->objective);
  }
}
double HDual::uOpRsDensityRec(double lc_OpRsDensity, double* opRsDensity, double* opRsAvDensity, double* opRsAvLog10Density) {
  double log10_lc_OpRsDensity;
  if (lc_OpRsDensity <= 0) {
    log10_lc_OpRsDensity = -20;
  } else {
    log10_lc_OpRsDensity = log(lc_OpRsDensity)/log(10.0);
  }
  *opRsAvDensity = (1-runningAverageMu)*(*opRsAvDensity) + runningAverageMu*lc_OpRsDensity;
  *opRsAvLog10Density = (1-runningAverageMu)*(*opRsAvLog10Density) + runningAverageMu*log10_lc_OpRsDensity;
  double rtV;
  rtV = (*opRsAvDensity);
  //  rtV = pow(10.0, *opRsAvLog10Density);  
  return rtV;
}

void HDual::chooseRow()
{
  if (invertHint)
    return;
  for (;;)
    {
      // Choose row
      dualRHS.choose_normal(&rowOut);
      if (rowOut == -1)
	{
	  invertHint = invertHint_possiblyOptimal; // Was 1
	  return;
	}
      
      // Verify weight
      model->timer.recordStart(HTICK_BTRAN);
      row_ep.clear();
      row_ep.count = 1;
      row_ep.index[0] = rowOut;
      row_ep.array[rowOut] = 1;
      row_ep.packFlag = true;
#ifdef HiGHSDEV
      if (AnIterLg) iterateOpRecBf(AnIterOpTy_Btran, row_ep, row_epDensity);
#endif      
      factor->btran(row_ep, row_epDensity);
#ifdef HiGHSDEV
      if (AnIterLg) iterateOpRecAf(AnIterOpTy_Btran, row_ep);
#endif      
      model->timer.recordFinish(HTICK_BTRAN);
      if (EdWt_Mode == EdWt_Mode_DSE)
	{
	  //For DSE compute the correct weight c_weight and use if to see how
	  //accurate the updated weight is.
	  double u_weight = dualRHS.workEdWt[rowOut];
	  double c_weight = dualRHS.workEdWt[rowOut] = row_ep.norm2();
	  //For DSE compute the correct weight c_weight and use if to see how
	  //accurate the updated weight is.
#ifdef HiGHSDEV
	  //			double DSE_wt_er = abs((u_weight - c_weight) / max(u_weight, c_weight));
	  //if (DSE_wt_er > 1e-2)
	  //			  printf(
	  //	 " !! JAJH RARE PRINT: Iter %d: DSE_wt_er = %8g = (%8g - %8g)/max(u_weight,c_weight)\n",
	  //	 model->numberIteration, DSE_wt_er, u_weight, c_weight);
#endif
	  if (u_weight >= 0.25 * c_weight)
	    break;
#ifdef HiGHSDEV
	  //printf(
	  //       " !! JAJH RARE PRINT: Iter %d: DSE_wt_er = %8g from %8g, %8g with u_weight/c_weight = %g < 0.25\n",
	  //       model->numberIteration, DSE_wt_er, u_weight, c_weight,
	  //       u_weight / c_weight);
	  n_wg_DSE_wt += 1;
#endif
	}
      else
	{
	  //If not using DSE then accept the row by breaking out of the loop
	  break;
	}
    }
  
  // Assign basic info
  columnOut = model->getBaseIndex()[rowOut];
  if (baseValue[rowOut] < baseLower[rowOut])
    deltaPrimal = baseValue[rowOut] - baseLower[rowOut];
  else
    deltaPrimal = baseValue[rowOut] - baseUpper[rowOut];
  sourceOut = deltaPrimal < 0 ? -1 : 1;
  double lc_OpRsDensity = 1.0 * row_ep.count / numRow;
  row_epDensity = uOpRsDensityRec(lc_OpRsDensity, &row_epDensity, &row_epAvDensity, &row_epAvLog10Density); 
  //row_epDensity = (1-runningAverageMu) * row_epDensity + runningAverageMu * (1.0 * row_ep.count / numRow);
}

void HDual::chooseColumn(HVector *row_ep)
{
  
  if (invertHint)
    return;
  
  // Compute pivot row
  model->timer.recordStart(HTICK_PRICE);
  row_ap.clear();
  
  bool anPriceEr = false;
  bool useUltraPrice = alw_price_ultra
    && row_apDensity*numCol*10 < row_ap.ilP2
				 && row_apDensity < 1e-3;
#ifdef HiGHSDEV
  //  bool rpPriceTy = false;
  //  int lc_numIt = model->numberIteration;
  //  if (lc_numIt<1 || rpPriceTy) printf("Iteration %d: Before PRICE: Mode = %d; ByColSw = %d; ByRowSw = %d; Ultra = %d\n",
  //	   lc_numIt, Price_Mode, alw_price_by_col_sw, alw_price_by_row_sw, alw_price_ultra);
#endif
  if (Price_Mode == Price_Mode_Col) {
    //Column-wise PRICE
#ifdef HiGHSDEV
    //    if (lc_numIt<1 || rpPriceTy) printf("Using column price\n");
    if (AnIterLg) {
      iterateOpRecBf(AnIterOpTy_Price, *row_ep, 0.0);
      AnIterNumColPrice++;
    }
#endif
    matrix->price_by_col(row_ap, *row_ep);
  } else {
    //Row-wise PRICE
    double lc_dsty = 1.0 * (*row_ep).count/numRow;
    if (alw_price_by_col_sw && (lc_dsty > dstyColPriceSw)) {
      //Avoid row price due to density of row_ep
#ifdef HiGHSDEV
      //      if (lc_numIt<1 || rpPriceTy) printf("PRICE By column, %d\n", numCol);
      if (AnIterLg) {
	iterateOpRecBf(AnIterOpTy_Price, *row_ep, 0.0);
	AnIterNumColPrice++;
      }
#endif
      matrix->price_by_col(row_ap, *row_ep);
      for (int col = 0; col < numCol; col++) {
	row_ap.array[col] = model->nonbasicFlag[col]*row_ap.array[col];
      }
      
    } else if (useUltraPrice) {
#ifdef HiGHSDEV
      //      if (lc_numIt<1 || rpPriceTy) printf("PRICE By row - ultra\n");
      if (AnIterLg) {
	iterateOpRecBf(AnIterOpTy_Price, *row_ep, row_apDensity);
	AnIterNumRowPriceUltra++;
      }
#endif
      matrix->price_by_row_ultra(row_ap, *row_ep);
      if (anPriceEr) {
	bool price_er;
	price_er = matrix->price_er_ck(row_ap, *row_ep);
	if (!price_er) printf("No ultra PRICE error\n");
      }
    } else if (alw_price_by_row_sw) {
      //Avoid Hyper Price on current density of result or switch if
      //the density of this Price becomes extreme
#ifdef HiGHSDEV
      //      if (lc_numIt<1 || rpPriceTy) printf("PRICE By row - switch\n");
      if (AnIterLg) {
	iterateOpRecBf(AnIterOpTy_Price, *row_ep, row_apDensity);
	AnIterNumRowPriceWSw++;
      }
#endif
      const double sw_dsty = matrix->price_by_row_sw_dsty;
      matrix->price_by_row_w_sw(row_ap, *row_ep, row_apDensity, 0, sw_dsty);
    } else {
      //No avoiding Hyper Price on current density of result or
      //switching if the density of this Price becomes extreme
#ifdef HiGHSDEV
      //      if (lc_numIt<1 || rpPriceTy) printf("PRICE By row - vanilla\n");
      if (AnIterLg) {
	iterateOpRecBf(AnIterOpTy_Price, *row_ep, 0.0);
	AnIterNumRowPrice++;
      }
#endif
      matrix->price_by_row(row_ap, *row_ep);
    }
  }
  if (anPriceEr) {
    //    bool price_er =
    matrix->price_er_ck(row_ap, *row_ep);
    //    if (!price_er) printf("No PRICE error\n");
  }
  double lc_OpRsDensity = 1.0 * row_ap.count / numCol;
  row_apDensity = uOpRsDensityRec(lc_OpRsDensity, &row_apDensity, &row_apAvDensity, &row_apAvLog10Density); 
  //row_apDensity = (1-runningAverageMu) * row_apDensity + runningAverageMu * (1.0 * row_ap.count / numCol);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecAf(AnIterOpTy_Price, row_ap);
#endif
  model->timer.recordFinish(HTICK_PRICE);
  
  // Choose column - possible
  model->timer.recordStart(HTICK_CHUZC0);
  dualRow.clear();
  dualRow.workDelta = deltaPrimal;
  dualRow.create_Freemove(row_ep);
  model->timer.recordFinish(HTICK_CHUZC0);
  model->timer.recordStart(HTICK_CHUZC1);
  dualRow.choose_makepack(&row_ap, 0);
  dualRow.choose_makepack(row_ep, numCol);
  dualRow.choose_possible();
  model->timer.recordFinish(HTICK_CHUZC1);
  
  // Choose column - check problem
  columnIn = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0)
    {
      invertHint = invertHint_possiblyDualUnbounded; // Was 1
      return;
    }
  
  // Choose column - final
  bool chooseColumnFail = dualRow.choose_final();
  if (chooseColumnFail) {
    invertHint = invertHint_chooseColumnFail;
    return;
  }
  model->timer.recordStart(HTICK_CHUZC4);
  dualRow.delete_Freemove();
  model->timer.recordFinish(HTICK_CHUZC4);
  
  columnIn = dualRow.workPivot;
  alphaRow = dualRow.workAlpha;
  thetaDual = dualRow.workTheta;
  
#ifdef HiGHSDEV
  if (rp_hsol && rowOut >= 0) {printf("\nPvR: Row %2d\n", rowOut); dualRow.rp_hsol_pv_r();}
#endif
  if (EdWt_Mode == EdWt_Mode_Dvx)
    {
      model->timer.recordStart(HTICK_DEVEX_WT);
      //
      //Determine the exact Devex weight
      //
      //    	int dvx_ix_o_columnIn = dvx_ix[columnIn];
      //    	int si_it_n = model->numberIteration;
      double og_dvx_wt_o_rowOut = dualRHS.workEdWt[rowOut];
      int vr_t_lv_bs = model->getBaseIndex()[rowOut];
      int dvx_ix_o_vr_t_lv_bs = dvx_ix[vr_t_lv_bs];
      if (vr_t_lv_bs < numCol)
	{
	  //    		Structural leaving the basis
	  if (rp_dvx)
	    printf("\n DvxIt %d; RowOut %d; VrOut %d Str: DvxIx %d\n",
		   n_dvx_it, rowOut, vr_t_lv_bs, dvx_ix_o_vr_t_lv_bs);
	}
      else
	{
	  //    		Logical leaving the basis
	  if (rp_dvx)
	    printf("\n DvxIt %d; RowOut %d; VrOut %d  Lg: DvxIx %d\n",
		   n_dvx_it, rowOut, vr_t_lv_bs, dvx_ix_o_vr_t_lv_bs);
	}
      double tru_dvx_wt_o_rowOut = 0;
      //    	Have to loop over [row_ap; row_ep]
      for (int el_n = 0; el_n < dualRow.packCount; el_n++)
	{
	  int vr_n = dualRow.packIndex[el_n];
	  double pv = dvx_ix[vr_n] * dualRow.packValue[el_n];
	  tru_dvx_wt_o_rowOut += pv * pv;
	}
      if (rp_dvx)
	printf("Sum of squares for Devex weight = %g\n",
	       tru_dvx_wt_o_rowOut);
      tru_dvx_wt_o_rowOut = max(1.0, tru_dvx_wt_o_rowOut);
      if (rp_dvx)
	printf("tru_dvx_wt_o_rowOut = %g; alphaRow = %g\n",
	       tru_dvx_wt_o_rowOut, alphaRow);
      //Analyse the Devex weight to determine whether a new framework should be set up
      double dvx_rao = max(og_dvx_wt_o_rowOut / tru_dvx_wt_o_rowOut,
			   tru_dvx_wt_o_rowOut / og_dvx_wt_o_rowOut);
      int i_te = numRow / nw_dvx_fwk_fq;
      if (rp_dvx)
	printf("i_te = %d; %d: dvx_rao = %8.2g\n", i_te,
	       max(mn_n_dvx_it, i_te), dvx_rao);
      i_te = max(mn_n_dvx_it, i_te);
      //		Square tl_dvx_wt due to keeping squared weights
      nw_dvx_fwk = dvx_rao > tl_dvx_wt * tl_dvx_wt || n_dvx_it > i_te;
      if (nw_dvx_fwk)
	if (rp_dvx)
	  printf("!!NEW DEVEX FRAMEWORK!!\n");
      dualRHS.workEdWt[rowOut] = tru_dvx_wt_o_rowOut;
      model->timer.recordFinish(HTICK_DEVEX_WT);
    }
  return;
}

void HDual::chooseColumn_slice(HVector *row_ep)
{
  if (invertHint)
    return;
  
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
  for (int i = 0; i < slice_num; i++)
    {
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
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0)
    {
      invertHint = invertHint_possiblyDualUnbounded; // Was 1
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

void HDual::updateFtran()
{
  if (invertHint) return;
  model->timer.recordStart(HTICK_FTRAN);
  column.clear();
  column.packFlag = true;
  matrix->collect_aj(column, columnIn, 1);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecBf(AnIterOpTy_Ftran, column, columnDensity);
#endif
  factor->ftran(column, columnDensity);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecAf(AnIterOpTy_Ftran, column);
#endif
  alpha = column.array[rowOut];
  model->timer.recordFinish(HTICK_FTRAN);
}

void HDual::updateFtranBFRT()
{
  if (invertHint) return;
  
  // Only time updateFtranBFRT if dualRow.workCount > 0;
  // If dualRow.workCount = 0 then dualRow.update_flip(&columnBFRT)
  // merely clears columnBFRT so no FTRAN is performed
  bool time_updateFtranBFRT = dualRow.workCount > 0;
  
  if (time_updateFtranBFRT) model->timer.recordStart(HTICK_FTRAN_BFRT);
  
  dualRow.update_flip(&columnBFRT);
  
  if (columnBFRT.count) {
#ifdef HiGHSDEV
    if (AnIterLg) iterateOpRecBf(AnIterOpTy_FtranBFRT, columnBFRT, columnDensity);
#endif
    factor->ftran(columnBFRT, columnDensity);
#ifdef HiGHSDEV
    if (AnIterLg) iterateOpRecAf(AnIterOpTy_FtranBFRT, columnBFRT);
#endif
  }
  if (time_updateFtranBFRT) model->timer.recordFinish(HTICK_FTRAN_BFRT);
}

void HDual::updateFtranDSE(HVector *DSE_Vector)
{
  if (invertHint) return;
  model->timer.recordStart(HTICK_FTRAN_DSE);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecBf(AnIterOpTy_FtranDSE, *DSE_Vector, rowdseDensity);
#endif
  factor->ftran(*DSE_Vector, rowdseDensity);
#ifdef HiGHSDEV
  if (AnIterLg) iterateOpRecAf(AnIterOpTy_FtranDSE, *DSE_Vector);
#endif
  model->timer.recordFinish(HTICK_FTRAN_DSE);
}

void HDual::updateVerify() {
  if (invertHint) return;
  
  // The alpha
  double aCol = fabs(alpha);
  double aRow = fabs(alphaRow);
  double aDiff = fabs(aCol - aRow);
  numericalTrouble = aDiff / min(aCol, aRow);
  if (numericalTrouble > 1e-7 && model->countUpdate > 0) {
    invertHint = invertHint_possiblySingularBasis; // Was 1
  }
  
  // We get this thing, but it is not actived by default.
  //    // The dual reduced cost
  //    double dualin_u = workDual[columnIn];
  //    double dualin_c = model->getWorkCost()[columnIn];
  //    for (int i = 0; i < column.count; i++) {
  //        int iRow = column.index[i];
  //        int iCol = model->getBaseIndex()[iRow];
  //        double value = column.array[iRow];
  //        double cost = model->getWorkCost()[iCol] + model->getWorkShift()[iCol];
  //        dualin_c -= cost * value;
  //    }
  //    double dualin_diff = fabs(dualin_c - dualin_u);
  ////    if (dualin_diff > Td) {
  ////        cout << dualin_c << "\t" << dualin_u << "\t" << dualin_diff << endl;
  ////        invertHint = invertHint_possiblySingularBasis; // Was 1
  ////    }
}

void HDual::updateDual()
{
  if (invertHint) return;
  
  // Update - dual (shift and back)
  if (thetaDual == 0)
    model->shiftCost(columnIn, -workDual[columnIn]);
  else
    {
      dualRow.update_dual(thetaDual);
      if (dual_variant != HDUAL_VARIANT_PLAIN && slice_PRICE)
	for (int i = 0; i < slice_num; i++)
	  slice_dualRow[i].update_dual(thetaDual);
    }
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;
  model->shiftBack(columnOut);
}

void HDual::updatePrimal(HVector *DSE_Vector)
{
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
  if (EdWt_Mode == EdWt_Mode_DSE)
    {
      double thisEdWt = dualRHS.workEdWt[rowOut] / (alpha * alpha);
      dualRHS.update_weight(&column, thisEdWt, -2 / alpha,
			    &DSE_Vector->array[0]);
      dualRHS.workEdWt[rowOut] = thisEdWt;
    }
  else if (EdWt_Mode == EdWt_Mode_Dvx)
    {
      //	Pivotal row is for the current basis: weights are required for the next basis
      //  so have to divide the current (exact) weight by the pivotal value
      double thisEdWt = dualRHS.workEdWt[rowOut] / (alpha * alpha);
      double dvx_wt_o_rowOut = max(1.0, thisEdWt);
      //       	      nw_wt  is max(workEdWt[iRow], NewExactWeight*columnArray[iRow]^2);
      //				But NewExactWeight is dvx_wt_o_rowOut = max(1.0, dualRHS.workEdWt[rowOut] / (alpha * alpha)) so
      //       	      nw_wt = max(workEdWt[iRow], dvx_wt_o_rowOut*columnArray[iRow]^2);
      //	Update rest of weights
      dualRHS.update_weight_Dvx(&column, dvx_wt_o_rowOut);
      dualRHS.workEdWt[rowOut] = dvx_wt_o_rowOut;
      n_dvx_it += 1;
    }
  dualRHS.update_infeasList(&column);
  
  if (EdWt_Mode == EdWt_Mode_DSE)
    {
      double lc_OpRsDensity = 1.0 * DSE_Vector->count / numRow;
      rowdseDensity = uOpRsDensityRec(lc_OpRsDensity, &rowdseDensity, &rowdseAvDensity, &rowdseAvLog10Density); 
      //rowdseDensity = (1-runningAverageMu) * rowdseDensity + runningAverageMu * (1.0 * DSE_Vector->count / numRow);
    }
  double lc_OpRsDensity = 1.0 * column.count / numRow;
  columnDensity = uOpRsDensityRec(lc_OpRsDensity, &columnDensity, &columnAvDensity, &columnAvLog10Density); 
  //columnDensity = (1-runningAverageMu) * columnDensity + runningAverageMu * (1.0 * column.count / numRow);
  
  total_fake += column.fakeTick;
  if (EdWt_Mode == EdWt_Mode_DSE)
    {
      total_fake += DSE_Vector->fakeTick;
    }
  total_FT_inc_TICK += column.pseudoTick;
  if (EdWt_Mode == EdWt_Mode_DSE)
    {
      total_FT_inc_TICK += DSE_Vector->pseudoTick;
    }
}

void HDual::updatePivots() {
  if (invertHint) return;
  
  // Update - pivots
  model->updatePivots(columnIn, rowOut, sourceOut);
  model->recordPivots(columnIn, columnOut, alpha);
  model->updateFactor(&column, &row_ep, &rowOut, &invertHint);
  model->updateMatrix(columnIn, columnOut);
  dualRow.delete_Freelist(columnIn);
  dualRHS.update_pivots(rowOut,
			model->getWorkValue()[columnIn] + thetaPrimal);
  
  bool reinvert_pseudoClock = total_fake >= factor->fakeTick;
  
  //	reinvert_pseudoClock = 2*total_fake >= factor->fakeTick;
  if (reinvert_pseudoClock && model->countUpdate >= 50) {
    //        cout << total_fake << "\t" << factor->fakeTick << endl;
    //        printf(
    //                "%10d   %10.2e  %10.2e  %10.4f          %10.2e  %10.2e  %10.4f\n",
    //                model->countUpdate, total_INVERT_TICK, total_FT_inc_TICK,
    //                total_FT_inc_TICK / total_INVERT_TICK, factor->fakeTick,
    //                total_fake, total_fake / factor->fakeTick);
    invertHint = invertHint_pseudoClockSaysInvert; // Was 1
  }
#ifdef HiGHSDEV
  //  if (invertHint) {printf("HDual::updatePivots: invertHint = %d\n", invertHint);}
#endif
  
}

void HDual::handleRankDeficiency() {
  int *baseIndex = model->getBaseIndex();
  bool rp = true;
  if (rp) rp = numRow<100;
  if (rp) {
    printf("handleRankDeficiency:\nBaseI  "); for (int i = 0; i < numRow; i++) printf(" %2d", baseIndex[i]);
    printf("\n");
  }
  for (int i = 0; i < numRow; i++) {
    if (baseIndex[i] < 0) {
      columnOut = -baseIndex[i]-1;
      baseIndex[i] = columnOut;
      rowOut = i;
      columnIn = numCol + rowOut;
      sourceOut = model->setSourceOutFmBd(columnOut);
      printf("handleRankDeficiency: columnOut=%d, columnIn=%d, rowOut=%d\n", columnOut, columnIn, rowOut);
      model->updatePivots(columnIn, rowOut, sourceOut);
      model->recordPivots(columnIn, columnOut, alpha);
      model->updateMatrix(columnIn, columnOut);
    }
  }
}

void HDual::iz_dvx_fwk()
{
  //Initialise the Devex framework: reference set is all basic variables
  model->timer.recordStart(HTICK_DEVEX_IZ);
  const int *NonbasicFlag = model->getNonbasicFlag();
  for (int vr_n = 0; vr_n < numTot; vr_n++)
    {
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
      // NonbasicFlag[vr_n]*NonbasicFlag[vr_n] evaluates faster than abs(NonbasicFlag[vr_n])
      //
      //int dvx_ix_o_vr = 1-abs(NonbasicFlag[vr_n]);
      //
      int dvx_ix_o_vr = 1-NonbasicFlag[vr_n]*NonbasicFlag[vr_n];
      dvx_ix[vr_n] = dvx_ix_o_vr;
      //      if (dvx_ix[vr_n] != dvx_ix_o_vr) printf("%1d = dvx_ix[%6d] != dvx_ix_o_vr = %1d\n", dvx_ix[vr_n], vr_n, dvx_ix_o_vr);
      
    }
  //  for (int i = 0; i < numRow; i++) dualRHS.workEdWt[i] = 1.0;
  dualRHS.workEdWt.assign(numRow, 1.0);
  n_dvx_it = 0;
  n_dvx_fwk += 1;
  nw_dvx_fwk = false;
  model->timer.recordFinish(HTICK_DEVEX_IZ);
}

void HDual::setCrash(const char *Crash_ArgV)
{
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
    Crash_Mode = Crash_Mode_Dev;
  else if (strcmp(Crash_ArgV, "Bs") == 0)
    Crash_Mode = Crash_Mode_Bs;
#ifdef HiGHSDEV
  else if (strcmp(Crash_ArgV, "TsSing") == 0)
    Crash_Mode = Crash_Mode_TsSing;
#endif
  else
    {
      cout << "HDual::setCrash unrecognised CrashArgV = " << Crash_ArgV
	   << " - using No crash" << endl;
      Crash_Mode = Crash_Mode_No;
    }
  //	if (Crash_Mode == Crash_Mode_LTSSF) {
  //		crash.ltssf_iz_mode(Crash_Mode);
  //	}
  //		cout<<"HDual::setCrash Crash_Mode = " << Crash_Mode << endl;
}

void HDual::setPrice(const char *Price_ArgV)
{
  //  cout<<"HDual::setPrice Price_ArgV = "<<Price_ArgV<<endl;
  alw_price_by_col_sw = false;
  alw_price_by_row_sw = false;
  alw_price_ultra = false;
  if (strcmp(Price_ArgV, "Col") == 0) {
    Price_Mode = Price_Mode_Col;
  }
  else if (strcmp(Price_ArgV, "Row") == 0) {
    Price_Mode = Price_Mode_Row;
  }
  else if (strcmp(Price_ArgV, "RowSw") == 0)
    {
      Price_Mode = Price_Mode_Row;
      alw_price_by_row_sw = true;
    }
  else if (strcmp(Price_ArgV, "RowSwColSw") == 0)
    {
      Price_Mode = Price_Mode_Row;
      alw_price_by_col_sw = true;
      alw_price_by_row_sw = true;
    }
  else if (strcmp(Price_ArgV, "RowUltra") == 0)
    {
      Price_Mode = Price_Mode_Row;
      alw_price_by_col_sw = true;
      alw_price_by_row_sw = true;
      alw_price_ultra = true;
    }
  else
    {
      cout << "HDual::setPrice unrecognised PriceArgV = " << Price_ArgV
	   << " - using row Price with switch or colump price switch" << endl;
      Price_Mode = Price_Mode_Row;
      alw_price_by_col_sw = true;
      alw_price_by_row_sw = true;
    }
}

void HDual::setEdWt(const char *EdWt_ArgV)
{
  //	cout<<"HDual::setEdWt EdWt_ArgV = "<<EdWt_ArgV<<endl;
  if (strcmp(EdWt_ArgV, "Dan") == 0)
    EdWt_Mode = EdWt_Mode_Dan;
  else if (strcmp(EdWt_ArgV, "Dvx") == 0)
    EdWt_Mode = EdWt_Mode_Dvx;
  else if (strcmp(EdWt_ArgV, "DSE") == 0)
    {
      EdWt_Mode = EdWt_Mode_DSE;
      iz_DSE_wt = true;
      alw_DSE2Dvx_sw = false;
    }
  else if (strcmp(EdWt_ArgV, "DSE0") == 0)
    {
      EdWt_Mode = EdWt_Mode_DSE;
      iz_DSE_wt = false;
      alw_DSE2Dvx_sw = false;
    }
  else if (strcmp(EdWt_ArgV, "DSE2Dvx") == 0)
    {
      EdWt_Mode = EdWt_Mode_DSE;
      iz_DSE_wt = true;
      alw_DSE2Dvx_sw = true;
    }
  else
    {
      cout << "HDual::setEdWt unrecognised EdWtArgV = " << EdWt_ArgV
	   << " - using DSE with possible switch to Devex" << endl;
      EdWt_Mode = EdWt_Mode_DSE;
      iz_DSE_wt = true;
      alw_DSE2Dvx_sw = true;
    }
  //	cout<<"HDual::setEdWt iz_DSE_wt = " << iz_DSE_wt << endl;
}

void HDual::setTimeLimit(double TimeLimit_ArgV)
{
  //	cout<<"HDual::setTimeLimit TimeLimit_ArgV = "<<TimeLimit_ArgV<<endl;
  TimeLimitValue = TimeLimit_ArgV;
}

void HDual::setPresolve(const char *Presolve_ArgV)
{
  //	cout<<"HDual::setPresolve Presolve_ArgV = "<<Presolve_ArgV<<endl;
  if (strcmp(Presolve_ArgV, "Off") == 0)
    Presolve_Mode = Presolve_Mode_Off;
  else if (strcmp(Presolve_ArgV, "On") == 0)
    Presolve_Mode = Presolve_Mode_On;
  else
    {
      cout << "HDual::setPresolve unrecognised PresolveArgV = " << Presolve_ArgV
	   << " - setting presolve off" << endl;
      Presolve_Mode = Presolve_Mode_Off;
    }
}

// Utility to get a row of the inverse of B for SCIP
int HDual::util_getBasisInvRow(int r, double *coef, int *inds, int *ninds)
{
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = r;
  row_ep.array[r] = 1;
  row_ep.packFlag = true;
  factor->btran(row_ep, row_epDensity);
  //  printf("util_getBasisInvRow: nnz = %4d/%4d\n", row_ep.count, numRow);
  for (int row = 0; row < numRow; row++)
    {
      //    printf("BasisInvRow(%4d) = %11g\n", row,  row_ep.array[row]);
      coef[row] = row_ep.array[row];
    }
  if (0 <= row_ep.count && row_ep.count <= numRow)
    {
      for (int ix = 0; ix < row_ep.count; ix++)
	inds[ix] = row_ep.index[ix];
      ninds[0] = row_ep.count;
    }
  else
    {
      printf("util_getBasisInvRow: row_ep.count < 0 or row_ep.count > numRow: %4d; %4d\n", row_ep.count, numRow);
      ninds[0] = -1;
    }
  cout << flush;
  return 0;
}

double HDual::an_bs_cond(HModel *ptr_model)
{
  model = ptr_model;
  // Alias to the matrix
  matrix = model->getMatrix();
  const int *Astart = matrix->getAstart();
  const double *Avalue = matrix->getAvalue();
  //Compute the Hager condition number estimate for the basis matrix
  double NoDensity = 1;
  bs_cond_x.resize(numRow);
  bs_cond_y.resize(numRow);
  bs_cond_z.resize(numRow);
  bs_cond_w.resize(numRow);
  //x = ones(n,1)/n;
  //y = A\x;
  double mu = 1.0 / numRow;
  double norm_Binv;
  for (int r_n = 0; r_n < numRow; r_n++)
    bs_cond_x[r_n] = mu;
  row_ep.clear();
  row_ep.count = numRow;
  for (int r_n = 0; r_n < numRow; r_n++)
    {
      row_ep.index[r_n] = r_n;
      row_ep.array[r_n] = bs_cond_x[r_n];
    }
  for (int ps_n = 1; ps_n <= 5; ps_n++)
    {
      row_ep.packFlag = false;
      factor->ftran(row_ep, NoDensity);
      //zeta = sign(y);
      for (int r_n = 0; r_n < numRow; r_n++)
	{
	  bs_cond_y[r_n] = row_ep.array[r_n];
	  if (bs_cond_y[r_n] > 0)
	    bs_cond_w[r_n] = 1.0;
	  else if (bs_cond_y[r_n] < 0)
	    bs_cond_w[r_n] = -1.0;
	  else
	    bs_cond_w[r_n] = 0.0;
	}
      //z=A'\zeta;
      row_ep.clear();
      row_ep.count = numRow;
      for (int r_n = 0; r_n < numRow; r_n++)
	{
	  row_ep.index[r_n] = r_n;
	  row_ep.array[r_n] = bs_cond_w[r_n];
	}
      row_ep.packFlag = false;
      factor->btran(row_ep, NoDensity);
      //norm_z = norm(z,'inf');
      //ztx = z'*x ;
      //NormEst = norm(y,1);
      //fd_i = 0;
      //for i=1:n
      //    if abs(z(i)) == norm_z
      //        fd_i = i;
      //        break
      //    end
      //end
      double norm_z = 0.0;
      double ztx = 0.0;
      norm_Binv = 0.0;
      int argmax_z = -1;
      for (int r_n = 0; r_n < numRow; r_n++)
	{
	  bs_cond_z[r_n] = row_ep.array[r_n];
	  double abs_z_v = abs(bs_cond_z[r_n]);
	  if (abs_z_v > norm_z)
	    {
	      norm_z = abs_z_v;
	      argmax_z = r_n;
	    }
	  ztx += bs_cond_z[r_n] * bs_cond_x[r_n];
	  norm_Binv += abs(bs_cond_y[r_n]);
	}
      //printf("%2d: ||z||_inf = %8.2g; z^T*x = %8.2g; ||y||_1 = %g\n", ps_n, norm_z, ztx, norm_Binv);
      if (norm_z <= ztx)
	break;
      //x = zeros(n,1);
      //x(fd_i) = 1;
      for (int r_n = 0; r_n < numRow; r_n++)
	bs_cond_x[r_n] = 0.0;
      row_ep.clear();
      row_ep.count = 1;
      row_ep.index[0] = argmax_z;
      row_ep.array[argmax_z] = 1.0;
      bs_cond_x[argmax_z] = 1.0;
    }
  double norm_B = 0.0;
  for (int r_n = 0; r_n < numRow; r_n++)
    {
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
  //  printf("Hager estimate of ||B^{-1}||_1 = %g; ||B||_1 = %g so cond_1(B) estimate is %g\n", norm_Binv, norm_B, cond_B);
  return cond_B;
}

#ifdef HiGHSDEV
void HDual::iterateRpIterDa(bool header) {
  if (header) {
    printf(" Inv       NumCk     LvR     LvC     EnC        DlPr        ThDu        ThPr          Aa");
  } else {
    printf(" %3d %11.4g %7d %7d %7d %11.4g %11.4g %11.4g %11.4g", 
	   invertHint, numericalTrouble,
	   rowOut, columnOut, columnIn, deltaPrimal, thetaDual, thetaPrimal, alpha);
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
  if (v > 0) intLog10V = log(v)/log(10.0);
  return intLog10V;
}

void HDual::iterateOpRecBf(int opTy, HVector& vector, double hist_dsty) {
  AnIterOpRec *AnIter = &AnIterOp[opTy];
  AnIter->AnIterOpNumCa++;
  double curr_dsty = 1.0 * vector.count / numRow;
  //  printf("%10s: %g<= %g;  %g<= %g\n", AnIter->AnIterOpName.c_str(),
  //	 curr_dsty, AnIter->AnIterOpHyperCANCEL,
  //	 hist_dsty, AnIter->AnIterOpHyperTRAN);
  if (curr_dsty <= AnIter->AnIterOpHyperCANCEL && hist_dsty <= AnIter->AnIterOpHyperTRAN) 
    AnIter->AnIterOpNumHyperOp++;
}

void HDual::iterateOpRecAf(int opTy, HVector& vector) {
  
  AnIterOpRec *AnIter = &AnIterOp[opTy];
  double rsDsty = 1.0 * vector.count / AnIter->AnIterOpRsDim;
  if (rsDsty <= hyperRESULT) AnIter->AnIterOpNumHyperRs++;
  AnIter->AnIterOpRsMxNNZ = max(vector.count, AnIter->AnIterOpRsMxNNZ);
  if (opTy == AnIterOpTy_Ftran) {
    //    printf("FTRAN: Iter %7d, NCa = %7d; NHS = %7d; RsDsty = %6.4f; AvgDsty = %6.4f\n", 
    //	   model->numberIteration,
    //	   AnIter->AnIterOpNumCa, AnIter->AnIterOpNumHyperRs, rsDsty, columnDensity);
  }
  if (rsDsty > 0) {
    AnIter->AnIterOpLog10RsDsty += log(rsDsty)/log(10.0);
  } else {
    double vectorNorm = 0;
    for (int index = 0; index < AnIter->AnIterOpRsDim; index++) {
      double vectorValue = vector.array[index];
      vectorNorm += vectorValue*vectorValue;
    }
    vectorNorm = sqrt(vectorNorm);
    printf("Strange: operation %s has result density = %g: ||vector|| = %g\n", AnIter->AnIterOpName.c_str(), rsDsty, vectorNorm);
  }
}

void HDual::iterateRpAn() {
  int AnIterNumIter = model->numberIteration - AnIterIt0;
  printf("\nAnalysis of %d iterations (%d to %d)\n", AnIterNumIter, AnIterIt0+1, model->numberIteration);
  if (AnIterNumIter <=0) return;
  int lc_EdWtNumIter;
  lc_EdWtNumIter = AnIterNumEdWtIt[EdWt_Mode_DSE];
  if (lc_EdWtNumIter>0) printf("DSE for %7d (%3d%%) iterations\n", lc_EdWtNumIter, (100*lc_EdWtNumIter)/AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[EdWt_Mode_Dvx];
  if (lc_EdWtNumIter>0) printf("Dvx for %7d (%3d%%) iterations\n", lc_EdWtNumIter, (100*lc_EdWtNumIter)/AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[EdWt_Mode_Dan];
  if (lc_EdWtNumIter>0) printf("Dan for %7d (%3d%%) iterations\n", lc_EdWtNumIter, (100*lc_EdWtNumIter)/AnIterNumIter);
  printf("\n");
  for (int k=0; k<NumAnIterOpTy; k++) {
    AnIterOpRec *AnIter = &AnIterOp[k];
    int lcNumCa = AnIter->AnIterOpSuNumCa;
    printf("\n%-9s performed %7d times\n", AnIter->AnIterOpName.c_str(), AnIter->AnIterOpSuNumCa);
    if (lcNumCa > 0) {
      int lcHyperOp = AnIter->AnIterOpSuNumHyperOp;
      int lcHyperRs = AnIter->AnIterOpSuNumHyperRs;
      int pctHyperOp = (100*lcHyperOp)/lcNumCa;
      int pctHyperRs = (100*lcHyperRs)/lcNumCa;
      double lcRsDsty = pow(10.0, AnIter->AnIterOpSuLog10RsDsty/lcNumCa);
      int lcNumNNz = lcRsDsty*AnIter->AnIterOpRsDim;
      int lcMxNNz = AnIter->AnIterOpRsMxNNZ;
      double lcMxNNzDsty = (1.0 * lcMxNNz)/AnIter->AnIterOpRsDim;
      printf("   %11d hyper-sparse operations (%3d%%)\n", lcHyperOp, pctHyperOp);
      printf("   %11d hyper-sparse results    (%3d%%)\n", lcHyperRs, pctHyperRs);
      printf("   %11.4g density of result (%d nonzeros)\n", lcRsDsty, lcNumNNz);
      printf("   %11.4g density of result with max (%d) nonzeros\n", lcMxNNzDsty, lcMxNNz);
    }
  }
  int NumInvert = 0;
  for (int k=1; k<=AnIterNumInvertHint; k++) NumInvert += AnIterNumInvert[k];
  if (NumInvert>0) {
    int lcNumInvert = 0;
    printf("\nInvert    performed %7d times: average frequency = %d\n", NumInvert, AnIterNumIter/NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_updateLimitReached];
    if (lcNumInvert > 0) printf("%7d (%3d%%) Invert operations due to update limit reached\n", lcNumInvert, (100*lcNumInvert)/NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_pseudoClockSaysInvert];
    if (lcNumInvert > 0) printf("%7d (%3d%%) Invert operations due to pseudo-clock\n", lcNumInvert, (100*lcNumInvert)/NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_possiblyOptimal];
    if (lcNumInvert > 0) printf("%7d (%3d%%) Invert operations due to possibly optimal\n", lcNumInvert, (100*lcNumInvert)/NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_possiblyPrimalUnbounded];
    if (lcNumInvert > 0) printf("%7d (%3d%%) Invert operations due to possibly primal unbounded\n", lcNumInvert, (100*lcNumInvert)/NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_possiblyDualUnbounded];
    if (lcNumInvert > 0) printf("%7d (%3d%%) Invert operations due to possibly dual unbounded\n", lcNumInvert, (100*lcNumInvert)/NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_possiblySingularBasis];
    if (lcNumInvert > 0) printf("%7d (%3d%%) Invert operations due to possibly singular basis\n", lcNumInvert, (100*lcNumInvert)/NumInvert);
    lcNumInvert = AnIterNumInvert[invertHint_primalInfeasibleInPrimalSimplex];
    if (lcNumInvert > 0) printf("%7d (%3d%%) Invert operations due to primal infeasible in primal simplex\n", lcNumInvert, (100*lcNumInvert)/NumInvert);
  }
  printf("\n%7d (%3d%%) primal degenerate iterations\n", AnIterNumPrDgnIt, (100*AnIterNumPrDgnIt)/AnIterNumIter);
  printf("%7d (%3d%%)   dual degenerate iterations\n", AnIterNumDuDgnIt, (100*AnIterNumDuDgnIt)/AnIterNumIter);
  int suPrice = AnIterNumColPrice + AnIterNumRowPrice + AnIterNumRowPriceWSw + AnIterNumRowPriceUltra;
  if (suPrice>0) {
    printf("\n%7d Price operations:\n", suPrice);
    printf("%7d Col Price      (%3d%%)\n", AnIterNumColPrice, (100*AnIterNumColPrice)/suPrice);
    printf("%7d Row Price      (%3d%%)\n", AnIterNumRowPrice, (100*AnIterNumRowPrice)/suPrice);  
    printf("%7d Row PriceWSw   (%3d%%)\n", AnIterNumRowPriceWSw, (100*AnIterNumRowPriceWSw/suPrice));
    printf("%7d Row PriceUltra (%3d%%)\n", AnIterNumRowPriceUltra, (100*AnIterNumRowPriceUltra/suPrice));
  }
  printf("\n%7d (%3d%%) costly DSE        iterations\n", AnIterNumCostlyDseIt, (100*AnIterNumCostlyDseIt)/AnIterNumIter);
  
  //
  //Add a record for the final iterations: may end up with one more
  //than AnIterTraceMxNumRec records, so ensure that there is enough
  //space in the arrays
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
    printf("   Iter ( FmIter: ToIter)     Time Iter/sec |  Col R_Ep R_Ap  DSE | EdWt | Aux0\n");
    for (int rec = 1; rec<=AnIterTraceNumRec; rec++) {
      lcAnIter = &AnIterTrace[rec];
      int toIter = lcAnIter->AnIterTraceIter;
      double toTime = lcAnIter->AnIterTraceTime;
      int dlIter = toIter-fmIter;
      if (rec<AnIterTraceNumRec && dlIter != AnIterTraceIterDl)
	printf("STRANGE: %d = dlIter != AnIterTraceIterDl = %d\n", dlIter, AnIterTraceIterDl);
      double dlTime = toTime-fmTime;
      int iterSpeed = 0;
      if (dlTime>0) iterSpeed = dlIter/dlTime;
      int lcEdWt_Mode = lcAnIter->AnIterTraceEdWt_Mode;
      int l10ColDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_Ftran]);
      int l10REpDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_Btran]);
      int l10RapDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_Price]);
      int l10DseDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE]);
      int l10Aux0 = intLog10(lcAnIter->AnIterTraceAux0);
      string strEdWt_Mode;
      if (lcEdWt_Mode == EdWt_Mode_DSE) strEdWt_Mode = "DSE";
      else if (lcEdWt_Mode == EdWt_Mode_Dvx) strEdWt_Mode = "Dvx";
      else if (lcEdWt_Mode == EdWt_Mode_Dan) strEdWt_Mode = "Dan";
      else strEdWt_Mode = "XXX";					
      printf("%7d (%7d:%7d) %8.4f  %7d | %4d %4d %4d %4d |  %3s | %4d\n",
	     dlIter, fmIter, toIter,
	     dlTime, iterSpeed,
	     l10ColDse, l10REpDse, l10RapDse, l10DseDse,
	     strEdWt_Mode.c_str(), l10Aux0);
      fmIter = toIter;
      fmTime = toTime;
    }
    printf("\n");
  }
}

void HDual::rp_hsol_da_str()
{
  printf("\nIteration %d\n", model->numberIteration);
  if (numTot > mx_rp_numTot) return;
  
  printf("\nData structures\n");
  printf("         ");
  for (int i = 0; i < numTot; i++) printf(" %4d", i); 
  printf("\n");
  printf("NonBcFg: "); for (int i = 0; i < numTot; i++) printf(" %4d", model->getNonbasicFlag()[i]); 
  printf("\n");
  printf("WorkMv:  "); for (int i = 0; i < numTot; i++) printf(" %4d", model->getNonbasicMove()[i]); 
  printf("\n");
  printf("DvxIx:   "); for (int i = 0; i < numTot; i++) printf(" %4d", dvx_ix[i]); 
  printf("\n");
  printf("BcIx:    "); for (int i = 0; i < numRow; i++) printf(" %4d", model->getBaseIndex()[i]); 
  printf("\n");
  printf("DvxV:    "); for (int i = 0; i < numRow; i++) printf(" %4.1g", dualRHS.workEdWt[i]); 
  printf("\n");
}
void HDual::rp_hsol_pv_c(HVector *column) const
{
  const double *columnArray = &column->array[0];
  // Alias
  //Set limits on problem size for reporting
  //	const int mx_rp_numTot = 20;
  if (numCol > mx_rp_numTot)
    return;
  printf("\nPvC: %2d  \n", columnIn);
  for (int i = 0; i < numRow; i++)
    {
      if (abs(columnArray[i]) > 1e-3)
	printf(" %2d: %4.1g\n", i, columnArray[i]);
    }
  printf("\n");
}

void HDual::rp_hsol_sol(HModel *ptr_model)
{
  const double *colCost = model->getcolCost();
  const double *colLower = model->getcolLower();
  const double *colUpper = model->getcolUpper();
  const double *rowLower = model->getrowLower();
  const double *rowUpper = model->getrowUpper();
  const int *nonbasicFlag = model->getNonbasicFlag();
  printf(" Column      Pr Act      Du Act   NonBc    Lower Bd    Upper Bd        Cost\n");
  for (int c_n = 0; c_n < numCol; c_n++)
    {
      int vr_n = c_n;
      printf("%7d %11.4g %11.4g %7d %11.4g %11.4g %11.4g\n", c_n, workValue[vr_n], workDual[vr_n], nonbasicFlag[vr_n], colLower[c_n], colUpper[c_n], colCost[c_n]);
    }
  printf("    Row      Pr Act      Du Act   NonBc    Lower Bd    Upper Bd\n");
  for (int r_n = 0; r_n < numRow; r_n++)
    {
      int vr_n = numCol + r_n;
      printf("%7d %11.4g %11.4g %7d %11.4g %11.4g\n", r_n, workValue[vr_n], workDual[vr_n], nonbasicFlag[vr_n], rowLower[r_n], rowUpper[r_n]);
    }
}

//TODO Put this in the right place - try to identify when workValue is set [nonbasic primals are set to the right bound for dual feasibility?]
void HDual::an_iz_vr_v()
{
  double norm_bc_pr_vr = 0;
  double norm_bc_du_vr = 0;
  for (int r_n = 0; r_n < numRow; r_n++)
    {
      int vr_n = model->getBaseIndex()[r_n];
      norm_bc_pr_vr += baseValue[r_n] * baseValue[r_n];
      norm_bc_du_vr += workDual[vr_n] * workDual[vr_n];
    }
  double norm_nonbc_pr_vr = 0;
  double norm_nonbc_du_vr = 0;
  for (int vr_n = 0; vr_n < numTot; vr_n++)
    {
      if (model->getNonbasicFlag()[vr_n])
	{
	  double pr_act_v = model->getWorkValue()[vr_n];
	  norm_nonbc_pr_vr += pr_act_v * pr_act_v;
	  norm_nonbc_du_vr += workDual[vr_n] * workDual[vr_n];
	}
    }
  //printf("Initial point has %d dual infeasibilities\n", dualInfeasCount);
  //printf("Norm of the basic    primal activites is %g\n", sqrt(norm_bc_pr_vr));
  //printf("Norm of the basic    dual   activites is %g\n", sqrt(norm_bc_du_vr));
  //printf("Norm of the nonbasic primal activites is %g\n", sqrt(norm_nonbc_pr_vr));
  //printf("Norm of the nonbasic dual   activites is %g\n", sqrt(norm_nonbc_du_vr));
}
#endif
