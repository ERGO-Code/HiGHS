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
#ifdef JAJH_dev
  printf("model->mlFg_Report() 1\n");
  cout << flush;
  model->mlFg_Report();
  cout << flush;
#endif
  //  printf("HDual::solve - model->mlFg_Report() 1\n");cout<<flush; model->mlFg_Report();cout<<flush;
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
  //Initialise numbers and times of rebuilds and inverts.
  totalRebuilds = 0;
  totalRebuildTime = 0;
  model->totalInverts = 0;
  model->totalInvertTime = 0;
  //Cannot solve box-constrained LPs
  if (model->numRow == 0)
    return;
#ifdef JAJH_dev
  printf("model->mlFg_Report() 2\n");
  cout << flush;
  model->mlFg_Report();
  cout << flush;
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
#ifdef JAJH_dev
    printf("Setting TimeLimitValue = %g\n", TimeLimitValue);
#endif
  }

  // Initialise working environment
  //Does LOTS, including initialisation of edge weights. Should only
  //be called if model dimension changes
  init(num_threads);

  model->initCost(1);
  if (!model->mlFg_haveFreshInvert)
    model->computeFactor();

  //Consider initialising edge weights
  //
  //NB workEdWt is assigned and initialised to 1s in
  //dualRHS.setup(model) so that CHUZR is well defined, even for

  









  //Consider initialising edge weights
  //
  //NB workEdWt is assigned and initialised to 1s in
  //dualRHS.setup(model) so that CHUZR is well defined, even for
  //Dantzig pricing
  //
#ifdef JAJH_dev
  printf("model->mlFg_haveEdWt 2 = %d;EdWt_Mode = %d; EdWt_Mode_DSE = %d\n", model->mlFg_haveEdWt, EdWt_Mode, EdWt_Mode_DSE);cout<<flush;
#endif
#ifdef JAJH_dev
    printf("Edge weights known? %d\n", !model->mlFg_haveEdWt);cout<<flush;
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
#ifdef JAJH_dev
      n_wg_DSE_wt = 0;
#endif
      int numBasicStructurals = numRow - model->numBasicLogicals;
#ifdef JAJH_dev
      printf("If (0<numBasicStructurals = %d) && %d = iz_DSE_wt: Compute exact DSE weights\n", numBasicStructurals, iz_DSE_wt);
#endif
      if (numBasicStructurals > 0 && iz_DSE_wt) {
	//Basis is not logical and DSE weights are to be initialised
#ifdef JAJH_dev
	printf("Compute exact DSE weights\n");
	int RpI = 1;
#endif
	double IzDseEdWtTT = model->timer.getTime();
	for (int i = 0; i < numRow; i++) {
#ifdef JAJH_dev
	  if (i==RpI) {
	      printf("Computeing exact DSE weight %d\n", i);
	      RpI = RpI*2;
	    }
#endif
	  row_ep.clear();
	  row_ep.count = 1;
	  row_ep.index[0] = i;
	  row_ep.array[i] = 1;
	  row_ep.packFlag = false;
	  factor->btran(row_ep, row_epDensity);
	  dualRHS.workEdWt[i] = row_ep.norm2();
	  row_epDensity *= 0.95;
	  row_epDensity += 0.05 * row_ep.count / numRow;
	}
	IzDseEdWtTT = model->timer.getTime() - IzDseEdWtTT;
#ifdef JAJH_dev
	printf("Computed %d initial DSE weights in %gs\n", numRow, IzDseEdWtTT);
#endif
	if (model->intOption[INTOPT_PRINT_FLAG]) {
	  printf("solve:: %d basic structurals: computed %d initial DSE weights in %gs, %d, %d, %g\n",
		 numBasicStructurals, numRow, IzDseEdWtTT, numBasicStructurals, numRow, IzDseEdWtTT);
	}
      } else {
	if (model->intOption[INTOPT_PRINT_FLAG]) {
	  printf("solve:: %d basic structurals: starting from B=I so unit initial DSE weights\n", numBasicStructurals);}
      }
    }
    //Indicate that edge weights are known
    model->mlFg_haveEdWt = 1;
  }
  
#ifdef JAJH_dev
  printf("model->mlFg_haveEdWt 3 = %d\n", model->mlFg_haveEdWt);cout<<flush;
#endif
  model->computeDual();
#ifdef JAJH_dev
  printf("Solve: Computed dual\n");cout<<flush;
#endif

#ifdef JAJH_dev
  bool rp_bs_cond = false;
  if (rp_bs_cond) {
    double bs_cond = an_bs_cond(ptr_model);
    printf("Initial basis condition estimate is %g\n", bs_cond );
  }
#endif
  
#ifdef JAJH_dev
  printf("HDual::solve - Calling model->computeDualInfeasInDual\n");cout<<flush;
#endif
  model->computeDualInfeasInDual(&dualInfeasCount);
  solvePhase = dualInfeasCount > 0 ? 1 : 2;
  
#ifdef JAJH_dev
  printf("Solve: Computed phase = %d\n", solvePhase);cout<<flush;
#endif

  // Find largest dual. No longer adjust the dual tolerance accordingly
  double largeDual = 0;
  for (int i = 0; i < numTot; i++) {
    if (model->getNonbasicFlag()[i]) {
      double myDual = fabs(workDual[i] * jMove[i]);
      if (largeDual < myDual)
	largeDual = myDual;
    }
  }
#ifdef JAJH_dev
  printf("Solve: Large dual = %g\n", largeDual);cout<<flush;
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
  bool ok = model->OKtoSolve(1, solvePhase);
  if (!ok) {printf("NOT OK TO SOLVE???\n");cout<<flush;}
  //  assert(ok);
  
  //  rp_hsol_sol(ptr_model);


  //	Analyse the initial values of primal and dual variables
  //  an_iz_vr_v();
  
  // The major solving loop
  
  while (solvePhase) {
#ifdef JAJH_dev
    int it0 = model->numberIteration;
    printf("HDual::solve Phase %d: Iteration %d; totalTime = %g; timer.getTime = %g\n", solvePhase, model->numberIteration, model->totalTime, model->timer.getTime());cout<<flush;
#endif
    switch (solvePhase) {
    case 1:
      solve_phase1();
#ifdef JAJH_dev
      n_ph1_du_it += (model->numberIteration - it0);
#endif
      break;
    case 2:
      solve_phase2();
#ifdef JAJH_dev
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
  
#ifdef H2DEBUG
  // Report the ticks before primal
  if (dual_variant == HDUAL_VARIANT_PLAIN) {
    int reportList[] = { HTICK_INVERT, HTICK_CHUZR1, HTICK_BTRAN,
			 HTICK_PRICE, HTICK_CHUZC0, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3, HTICK_CHUZC4,
			 HTICK_FTRAN, HTICK_FTRAN_MIX, HTICK_FTRAN_DSE,
			 HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,
			 HTICK_UPDATE_FACTOR };
    int reportCount = sizeof(reportList) / sizeof(int);
    model->timer.report(reportCount, reportList);
  }
  
  if (dual_variant == HDUAL_VARIANT_TASKS) {
    int reportList[] = { HTICK_INVERT, HTICK_CHUZR1, HTICK_BTRAN,
			 HTICK_PRICE, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3,
			 HTICK_FTRAN, HTICK_FTRAN_MIX, HTICK_FTRAN_DSE,
			 HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,
			 HTICK_UPDATE_FACTOR, HTICK_GROUP1, HTICK_GROUP2 };
    int reportCount = sizeof(reportList) / sizeof(int);
    model->timer.report(reportCount, reportList);
  }
  
  if (dual_variant == HDUAL_VARIANT_MULTI) {
    int reportList[] = { HTICK_INVERT, HTICK_CHUZR1, HTICK_BTRAN,
			 HTICK_PRICE, HTICK_CHUZC1, HTICK_CHUZC2, HTICK_CHUZC3,
			 HTICK_FTRAN, HTICK_FTRAN_MIX, HTICK_FTRAN_DSE,
			 HTICK_UPDATE_DUAL, HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT,
			 HTICK_UPDATE_FACTOR, HTICK_UPDATE_ROW_EP };
    int reportCount = sizeof(reportList) / sizeof(int);
    model->timer.report(reportCount, reportList);
    printf("PAMI   %-20s    CUTOFF  %6g    PERSISTENSE  %6g\n",
	   model->modelName.c_str(), model->dblOption[DBLOPT_PAMI_CUTOFF],
	   model->numberIteration / (1.0 + multi_iteration));
  }
#endif
  
  if (model->problemStatus != LP_Status_OutOfTime) {
    // Use primal to clean up if not out of time
#ifdef JAJH_dev
    int it0 = model->numberIteration;
#endif
    if (solvePhase == 4) {
      HPrimal hPrimal;
      hPrimal.TimeLimitValue = TimeLimitValue;
      hPrimal.solvePhase2(model);
      //Add in the count and time for any primal rebuilds
      totalRebuildTime += hPrimal.totalRebuildTime;
      totalRebuilds += hPrimal.totalRebuilds;
    }
#ifdef JAJH_dev
    n_pr_it += (model->numberIteration - it0);
#endif
  }
  // Save the solved results
  model->totalTime += model->timer.getTime();
  
#ifdef JAJH_dev
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
  ok = model->OKtoSolve(1, solvePhase);
  if (!ok) {printf("NOT OK After Solve???\n");cout<<flush;}
  //  assert(ok);
  printf("model->mlFg_Report() 9\n");cout<<flush;  model->mlFg_Report();cout<<flush;

  printf("Time: Total inverts =  %4d; Total invert  time = %11.4g of Total time = %11.4g", model->totalInverts, model->totalInvertTime, model->totalTime);
  if (model->totalTime > 0.001) 
    printf(" (%6.2f%%)\n", (100*model->totalInvertTime)/model->totalTime);
  else
    printf("\n");
  cout << flush;
  printf("Time: Total rebuilds = %4d; Total rebuild time = %11.4g of Total time = %11.4g", totalRebuilds, totalRebuildTime, model->totalTime);
  if (model->totalTime > 0.001) 
    printf(" (%6.2f%%)\n", (100*totalRebuildTime)/model->totalTime);
  else
    printf("\n");
  cout << flush;
 
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
  columnDensity = 0;
  row_epDensity = 0;
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
#ifdef JAJH_dev
  int lc_totalTime_rp_n = 0;
  //	printf("DualPh1: lc_totalTime = %5.2f; Record %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
  // Main solving structure
  for (;;)
  {
    rebuild();
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
      /*			if (model->objective > model->dblOption[DBLOPT_OBJ_UB]) {
#ifdef SCIP_dev
			  printf("HDual::solve_phase1: %g = Objective > dblOption[DBLOPT_OBJ_UB]\n", model->objective, model->dblOption[DBLOPT_OBJ_UB]);
#endif
			  model->problemStatus = LP_Status_ObjUB; 
                          break;
			}
			*/
    }
    lc_totalTime = model->totalTime + model->timer.getTime();
#ifdef JAJH_dev
    lc_totalTime_rp_n += 1;
    //		printf("DualPh1: lc_totalTime = %5.2f; Record %d\n", lc_totalTime, lc_totalTime_rp_n);
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
#ifdef JAJH_dev
  int lc_totalTime_rp_n = 0;
  printf("DualPh2: lc_totalTime = %5.2f; Record %d\n", lc_totalTime, lc_totalTime_rp_n);
#endif
  // Main solving structure
  for (;;)
  {
    // Outer loop of solve_phase2()
    // Rebuild all values, reinverting B if updates have been performed
    rebuild();
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
      if (invertHint)
        break;
      //      model->computeDuObj();
#ifdef JAJH_dev
      //      double pr_obj_v = model->computePrObj();
      //      printf("HDual::solve_phase2: Iter = %4d; Pr Obj = %.11g; Du Obj = %.11g\n",
      //	     model->numberIteration, pr_obj_v, model->objective);
#endif
      if (model->objective > model->dblOption[DBLOPT_OBJ_UB])
      {
#ifdef SCIP_dev
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
#ifdef JAJH_dev
    lc_totalTime_rp_n += 1;
    printf("DualPh2: lc_totalTime = %5.2f; Record %d\n", lc_totalTime, lc_totalTime_rp_n);
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

  if (SolveBailout)
    return;

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
  model->timer.recordStart(HTICK_INVERT);
#ifdef JAJH_dev
  //	double tt0 = model->timer.getTime();
#endif
	double tt0 = model->timer.getTime();
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
		for (int i = 0; i < numRow; i++)
			dualRHS.workEdWtFull[baseIndex[i]] = dualRHS.workEdWt[i];
		model->computeFactor();
		// Gather the edge weights according to the
		// permutation of baseIndex after INVERT
		for (int i = 0; i < numRow; i++)
			dualRHS.workEdWt[i] = dualRHS.workEdWtFull[baseIndex[i]];
	}

	// Recompute dual solution
	model->computeDual();
	model->correctDual(&dualInfeasCount);
	model->computePrimal();

	// Collect primal infeasible as a list
	dualRHS.create_infeasArray();
	dualRHS.create_infeasList(columnDensity);

	// Compute the objective value
	model->computeDuObj(solvePhase);
	model->util_reportNumberIterationObjectiveValue(sv_invertHint);

	total_INVERT_TICK = factor->pseudoTick;
	total_FT_inc_TICK = 0;
	total_fake = 0;

	model->timer.recordFinish(HTICK_INVERT);

	double rebuildTime = model->timer.getTime()-tt0;
	totalRebuilds++;
	totalRebuildTime += rebuildTime;
#ifdef JAJH_dev
  printf("Dual  Ph%-2d rebuild %4d (%1d) on iteration %9d: Rebuild time = %11.4g; Total rebuild time = %11.4g\n",
         solvePhase, totalRebuilds, sv_invertHint, model->numberIteration, rebuildTime, totalRebuildTime);
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
	model->util_reportNumberIterationObjectiveValue(-1);

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
	if (rp_hsol) rp_hsol_da_str();
	chooseRow();
	chooseColumn(&row_ep);
#ifdef JAJH_dev
	if (rp_hsol && rowOut >= 0) {printf("\nPvR: Row %2d\n", rowOut); dualRow.rp_hsol_pv_r();} cout<<flush;
#endif

	updateFtranBFRT();
	// updateFtran(); computes the pivotal column in the data structure "column"
	updateFtran();
	if (rp_hsol) rp_hsol_pv_c(&column);

	//updateFtranDSE performs the DSE FTRAN on pi_p
	if (EdWt_Mode == EdWt_Mode_DSE) updateFtranDSE(&row_ep);

	//updateVerify() Checks row-wise pivot against column-wise pivot for numerical trouble
	updateVerify();

	//updateDual() Updates the dual values
	updateDual();

	//updatePrimal(&row_ep); Updates the primal values and the edge weights
	updatePrimal(&row_ep);

	if ((EdWt_Mode == EdWt_Mode_Dvx) && (nw_dvx_fwk)) iz_dvx_fwk();

	//Update the basis representation
	updatePivots();

	//Possibly report on the iteration
	iterateRp();
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
    factor->btran(row_ep, row_epDensity);
    model->timer.recordFinish(HTICK_BTRAN);
    if (EdWt_Mode == EdWt_Mode_DSE)
    {
      //For DSE compute the correct weight c_weight and use if to see how
      //accurate the updated weight is.
      double u_weight = dualRHS.workEdWt[rowOut];
      double c_weight = dualRHS.workEdWt[rowOut] = row_ep.norm2();
      //For DSE compute the correct weight c_weight and use if to see how
      //accurate the updated weight is.
#ifdef JAJH_dev
      //			double DSE_wt_er = abs((u_weight - c_weight) / max(u_weight, c_weight));
      //if (DSE_wt_er > 1e-2)
      //			  printf(
      //	 " !! JAJH RARE PRINT: Iter %d: DSE_wt_er = %8g = (%8g - %8g)/max(u_weight,c_weight)\n",
      //	 model->numberIteration, DSE_wt_er, u_weight, c_weight);
#endif
      if (u_weight >= 0.25 * c_weight)
        break;
#ifdef JAJH_dev
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
  row_epDensity *= 0.95;
  row_epDensity += 0.05 * row_ep.count / numRow;
}

void HDual::chooseColumn(HVector *row_ep)
{

  if (invertHint)
    return;

  // Compute pivot row
  model->timer.recordStart(HTICK_PRICE);
  row_ap.clear();
  matrix->price_by_row(row_ap, *row_ep);
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
  dualRow.choose_final();
  model->timer.recordStart(HTICK_CHUZC4);
  dualRow.delete_Freemove();
  model->timer.recordFinish(HTICK_CHUZC4);

  columnIn = dualRow.workPivot;
  alphaRow = dualRow.workAlpha;
  thetaDual = dualRow.workTheta;

  if (rp_hsol && rowOut >= 0)
  {
    printf("\nPvR: Row %2d\n", rowOut);
    dualRow.rp_hsol_pv_r();
  }
  if (EdWt_Mode == EdWt_Mode_Dvx)
  {
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
  }
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

void HDual::updateFtranBFRT()
{
  if (invertHint)
    return;
  model->timer.recordStart(HTICK_FTRAN_MIX);
  dualRow.update_flip(&columnBFRT);
  if (columnBFRT.count)
    factor->ftran(columnBFRT, columnDensity);
  model->timer.recordFinish(HTICK_FTRAN_MIX);
}

void HDual::updateFtran()
{
  if (invertHint)
    return;
  model->timer.recordStart(HTICK_FTRAN);
  column.clear();
  column.packFlag = true;
  matrix->collect_aj(column, columnIn, 1);
  factor->ftran(column, columnDensity);
  alpha = column.array[rowOut];
  model->timer.recordFinish(HTICK_FTRAN);
}

void HDual::updateFtranDSE(HVector *DSE_Vector)
{
  if (invertHint)
    return;
  model->timer.recordStart(HTICK_FTRAN_DSE);
  factor->ftran(*DSE_Vector, rowdseDensity);
  model->timer.recordFinish(HTICK_FTRAN_DSE);
}

void HDual::updateVerify() {
	if (invertHint)
		return;

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
  if (invertHint)
    return;

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
  if (invertHint)
    return;
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
    rowdseDensity = 0.95 * rowdseDensity + 0.05 * DSE_Vector->count / numRow;
  }
  columnDensity = 0.95 * columnDensity + 0.05 * column.count / numRow;

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
	if (invertHint)
		return;

	// Update - pivots
	model->updatePivots(columnIn, rowOut, sourceOut);
	model->recordPivots(columnIn, columnOut, alpha);
	model->updateFactor(&column, &row_ep, &rowOut, &invertHint);
	model->updateMatrix(columnIn, columnOut);
	dualRow.delete_Freelist(columnIn);
	dualRHS.update_pivots(rowOut,
			model->getWorkValue()[columnIn] + thetaPrimal);

	if (total_fake >= factor->fakeTick && model->countUpdate >= 50) {
//        cout << total_fake << "\t" << factor->fakeTick << endl;
//        printf(
//                "%10d   %10.2e  %10.2e  %10.4f          %10.2e  %10.2e  %10.4f\n",
//                model->countUpdate, total_INVERT_TICK, total_FT_inc_TICK,
//                total_FT_inc_TICK / total_INVERT_TICK, factor->fakeTick,
//                total_fake, total_fake / factor->fakeTick);
	  invertHint = invertHint_pseudoClockSaysInvert; // Was 1
	}
#ifdef JAJH_dev
	if (invertHint) {printf("HDual::updatePivots: invertHint = %d\n", invertHint);}
#endif
	
}

void HDual::iz_dvx_fwk()
{
  //Initialise the Devex framework: reference set is all basic variables
  for (int vr_n = 0; vr_n < numTot; vr_n++)
  {
    if (model->getNonbasicFlag()[vr_n])
      //			Nonbasic variables not in reference set
      dvx_ix[vr_n] = dvx_not_in_R;
    else
      //			Basic variables in reference set
      dvx_ix[vr_n] = dvx_in_R;
  }
  for (int i = 0; i < numRow; i++)
    dualRHS.workEdWt[i] = 1.0;
  n_dvx_it = 0;
  n_dvx_fwk += 1;
  nw_dvx_fwk = false;
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
  }
  else if (strcmp(EdWt_ArgV, "DSE0") == 0)
  {
    EdWt_Mode = EdWt_Mode_DSE;
    iz_DSE_wt = false;
  }
  else if (strcmp(EdWt_ArgV, "DSE1") == 0)
  {
    EdWt_Mode = EdWt_Mode_DSE;
    iz_DSE_wt = true;
  }
  else
  {
    cout << "HDual::setEdWt unrecognised EdWtArgV = " << EdWt_ArgV
         << " - using DSE with exact initial weights" << endl;
    EdWt_Mode = EdWt_Mode_DSE;
    iz_DSE_wt = true;
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
  // printf("Hager estimate of ||B^{-1}||_1 = %g; ||B||_1 = %g so cond_1(B) estimate is %g\n", norm_Binv, norm_B, cond_B);
  return cond_B;
}

void HDual::rp_hsol_da_str()
{

  printf("\nIteration %d\n", model->numberIteration);
  if (numTot > mx_rp_numTot)
    return;

  printf("\nData structures\n");
  printf("         ");
  for (int i = 0; i < numTot; i++)
  {
    printf(" %4d", i);
  }
  printf("\n");

  printf("NonBcFg: ");
  for (int i = 0; i < numTot; i++)
  {
    printf(" %4d", model->getNonbasicFlag()[i]);
  }
  printf("\n");

  printf("WorkMv:  ");
  for (int i = 0; i < numTot; i++)
  {
    printf(" %4d", model->getNonbasicMove()[i]);
  }
  printf("\n");

  printf("DvxIx:   ");
  for (int i = 0; i < numTot; i++)
  {
    printf(" %4d", dvx_ix[i]);
  }
  printf("\n");

  printf("BcIx:    ");
  for (int i = 0; i < numRow; i++)
  {
    printf(" %4d", model->getBaseIndex()[i]);
  }
  printf("\n");

  printf("DvxV:    ");
  for (int i = 0; i < numRow; i++)
  {
    printf(" %4.1g", dualRHS.workEdWt[i]);
  }
  printf("\n");
}

void HDual::iterateRp() {
  //    if (model->intOption[INTOPT_PRINT_FLAG] != 4) return;
    int numIter = model->numberIteration;
  //deltaPrimal: Move to bound from basic value for leaving variable
  //thetaDual:   
  //thetaPrimal: 
  //alpha:
  //NumCk %11.4g; InvHint%2d; DuObj %11.4g; 
  //  printf("Iter %9d: Ph%1d; LvR %7d; LvC %7d; EnC %7d; DlPr %11.4g; ThDu %11.4g; ThPr %11.4g; Aa %11.4g\n", 
  //	 model->numberIteration,
  //	 //model->objective,
  //	 solvePhase,
  //	 //numericalTrouble,
  //	 //invertHint,
  //	 rowOut, columnOut, columnIn, deltaPrimal, thetaDual, thetaPrimal, alpha);
  //  //DuObj %11.4g;
    if (numIter % 10 == 1) 
      printf("     Iter Ph Inv       NumCk     LvR     LvC     EnC        DlPr        ThDu        ThPr          Aa   CD REpD RSeD\n");

    int l10ColDse = -99;
    int l10REpDse = -99;
    int l10DseDse = -99;
    if (columnDensity>0) l10ColDse = log(columnDensity)/log(10.0);
    if (row_epDensity>0) l10REpDse = log(row_epDensity)/log(10.0);
    if (rowdseDensity>0) l10DseDse = log(rowdseDensity)/log(10.0);

    printf("%9d %2d %3d %11.4g %7d %7d %7d %11.4g %11.4g %11.4g %11.4g %4d %4d %4d\n", 
	   numIter,
	 //model->objective,
	   solvePhase,
	   invertHint, numericalTrouble,
	   rowOut, columnOut, columnIn, deltaPrimal, thetaDual, thetaPrimal, alpha,
	   l10ColDse, l10REpDse, l10DseDse);
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
