
#include "cupdlp_solver.h"

#include "cupdlp_defs.h"
#include "cupdlp_linalg.h"
#include "cupdlp_proj.h"
#include "cupdlp_restart.h"
// #include "cupdlp_scaling.h"
// #include "cupdlp_scaling_new.h"
#include "cupdlp_step.h"
#include "cupdlp_utils.h"
#include "glbopts.h"

void debugPrintCupdlpVector(const char* name, const CUPDLPvec* vector) {
  printf("Variable %s: ", name);
  for (int ix = 0; ix < vector->len; ix++)
    printf("%11.6g ", vector->data[ix]);
  printf("\n");
}

void debugPrintDoubleVector(const char* name, const double* vector, const int n) {
  printf("Variable %s: ", name);
  for (int ix = 0; ix < n; ix++)
    printf("%11.6g ", vector[ix]);
  printf("\n");
}

void PDHG_Compute_Primal_Feasibility(CUPDLPwork *work, double *primalResidual,
                                     const double *ax, const double *x,
                                     double *dPrimalFeasibility,
                                     double *dPrimalObj) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPscaling *scaling = work->scaling;

  // primal variable violation

  // todo, add this
  //    *dPrimalObj = Dotprod_Neumaier(problem->cost, x, lp->nCols);
  cupdlp_dot(work, lp->nCols, x, problem->cost, dPrimalObj);
  *dPrimalObj = *dPrimalObj * problem->sign_origin + problem->offset;

  // cupdlp_copy(primalResidual, ax, cupdlp_float, lp->nRows);
  CUPDLP_COPY_VEC(primalResidual, ax, cupdlp_float, lp->nRows);

  // AddToVector(primalResidual, -1.0, problem->rhs, lp->nRows);
  cupdlp_float alpha = -1.0;
  cupdlp_axpy(work, lp->nRows, &alpha, problem->rhs, primalResidual);

  double dPrimalFeas = 0.0;

  // todo, check
  //  cupdlp_projNegative(primalResidual + problem->nEqs, primalResidual +
  //  problem->nEqs, lp->nRows - problem->nEqs);
  //

  cupdlp_projNeg(primalResidual + problem->nEqs, lp->nRows - problem->nEqs);

  if (scaling->ifScaled) {
    // cupdlp_edot(primalResidual, scaling->rowScale, lp->nRows);
    // cupdlp_edot(primalResidual, scaling->rowScale_gpu, lp->nRows);

    // cupdlp_copy_vec(work->buffer3, scaling->rowScale, cupdlp_float,
    // lp->nRows); cupdlp_edot(primalResidual, work->buffer3, lp->nRows);

    cupdlp_edot(primalResidual, work->rowScale, lp->nRows);
  }

  cupdlp_twoNorm(work, lp->nRows, primalResidual, dPrimalFeasibility);
}

void PDHG_Compute_Dual_Feasibility(CUPDLPwork *work, double *dualResidual,
                                   const double *aty, const double *x,
                                   const double *y, double *dDualFeasibility,
                                   double *dDualObj, double *dComplementarity) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  CUPDLPscaling *scaling = work->scaling;
  // todo, compute Neumaier
  //    *dDualObj = Dotprod_Neumaier(problem->rhs, y, lp->nRows);
  cupdlp_dot(work, lp->nRows, y, problem->rhs, dDualObj);
  *dDualObj = *dDualObj * problem->sign_origin + problem->offset;

  *dComplementarity = 0.0;
  // @note:
  // original dual residual in pdlp:
  // they compute:
  //    violation   +  reduced cost
  //  |max(-y, 0)|  + |(I-Π)(c-Α'υ)|
  // compute c - A'y

  CUPDLP_COPY_VEC(dualResidual, aty, cupdlp_float, lp->nCols);
  cupdlp_float alpha = -1.0;
  cupdlp_scaleVector(work, alpha, dualResidual, lp->nCols);

  alpha = 1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, problem->cost, dualResidual);

  //    julia version
  //        function compute_reduced_costs_from_primal_gradient_kernel!(
  //            primal_gradient::CuDeviceVector{Float64},
  //            isfinite_variable_lower_bound::CuDeviceVector{Bool},
  //            isfinite_variable_upper_bound::CuDeviceVector{Bool},
  //            num_variables::Int64,
  //            reduced_costs::CuDeviceVector{Float64},
  //            reduced_costs_violation::CuDeviceVector{Float64},
  //        )
  //            tx = threadIdx().x + (blockDim().x * (blockIdx().x - 0x1))
  //            if tx <= num_variables
  //                @inbounds begin
  //                    reduced_costs[tx] = max(primal_gradient[tx], 0.0) *
  //                    isfinite_variable_lower_bound[tx] +
  //                    min(primal_gradient[tx], 0.0) *
  //                    isfinite_variable_upper_bound[tx]
  //
  //                    reduced_costs_violation[tx] = primal_gradient[tx] -
  //                    reduced_costs[tx]
  //                end
  //            end
  //            return
  //        end

  // cupdlp_copy(resobj->dSlackPos, dualResidual, cupdlp_float, lp->nCols);
  CUPDLP_COPY_VEC(resobj->dSlackPos, dualResidual, cupdlp_float, lp->nCols);

  // cupdlp_projPositive(resobj->dSlackPos, resobj->dSlackPos, lp->nCols);
  cupdlp_projPos(resobj->dSlackPos, lp->nCols);

  // cupdlp_cdot_fb(resobj->dSlackPos, problem->hasLower, lp->nCols);
  cupdlp_edot(resobj->dSlackPos, problem->hasLower, lp->nCols);

  cupdlp_float temp = 0.0;
  cupdlp_dot(work, lp->nCols, x, resobj->dSlackPos, &temp);
  *dComplementarity += temp;
  cupdlp_dot(work, lp->nCols, resobj->dSlackPos, resobj->dLowerFiltered, &temp);
  *dComplementarity -= temp;
  cupdlp_dot(work, lp->nCols, resobj->dSlackPos, resobj->dLowerFiltered, &temp);
  *dDualObj += temp;

  CUPDLP_COPY_VEC(resobj->dSlackNeg, dualResidual, cupdlp_float, lp->nCols);

  cupdlp_projNeg(resobj->dSlackNeg, lp->nCols);

  // ScaleVector(-1.0, resobj->dSlackNeg, lp->nCols);
  cupdlp_scaleVector(work, -1.0, resobj->dSlackNeg, lp->nCols);

  // cupdlp_cdot_fb(resobj->dSlackNeg, problem->hasUpper, lp->nCols);
  cupdlp_edot(resobj->dSlackNeg, problem->hasUpper, lp->nCols);

  cupdlp_dot(work, lp->nCols, x, resobj->dSlackNeg, &temp);
  *dComplementarity -= temp;
  cupdlp_dot(work, lp->nCols, resobj->dSlackNeg, resobj->dUpperFiltered, &temp);
  *dComplementarity += temp;
  cupdlp_dot(work, lp->nCols, resobj->dSlackNeg, resobj->dUpperFiltered, &temp);
  *dDualObj -= temp;

  alpha = -1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, resobj->dSlackPos, dualResidual);
  alpha = 1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, resobj->dSlackNeg, dualResidual);

  if (scaling->ifScaled) {
    // cupdlp_edot(dualResidual, scaling->colScale, lp->nCols);
    // cupdlp_edot(dualResidual, scaling->colScale_gpu, lp->nCols);

    // cupdlp_copy_vec(work->buffer3, scaling->colScale, cupdlp_float,
    // lp->nCols); cupdlp_edot(dualResidual, work->buffer3, lp->nCols);

    cupdlp_edot(dualResidual, work->colScale, lp->nCols);
  }

  cupdlp_twoNorm(work, lp->nCols, dualResidual, dDualFeasibility);
}

void PDHG_Compute_Residuals(CUPDLPwork *work) {
#if problem_USE_TIMERS
  ++problem->nComputeResidualsCalls;
  double dStartTime = getTimeStamp();
#endif
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  CUPDLPiterates *iterates = work->iterates;
  CUPDLPscaling *scaling = work->scaling;
  CUPDLPsettings *settings = work->settings;

  PDHG_Compute_Primal_Feasibility(work, resobj->primalResidual,
                                  iterates->ax->data, iterates->x->data,
                                  &resobj->dPrimalFeas, &resobj->dPrimalObj);
  PDHG_Compute_Dual_Feasibility(work, resobj->dualResidual, iterates->aty->data,
                                iterates->x->data, iterates->y->data,
                                &resobj->dDualFeas, &resobj->dDualObj,
                                &resobj->dComplementarity);

  PDHG_Compute_Primal_Feasibility(
      work, resobj->primalResidualAverage, iterates->axAverage->data,
      iterates->xAverage->data, &resobj->dPrimalFeasAverage,
      &resobj->dPrimalObjAverage);
  PDHG_Compute_Dual_Feasibility(
      work, resobj->dualResidualAverage, iterates->atyAverage->data,
      iterates->xAverage->data, iterates->yAverage->data,
      &resobj->dDualFeasAverage, &resobj->dDualObjAverage,
      &resobj->dComplementarityAverage);

  // resobj->dPrimalObj /= (scaling->dObjScale * scaling->dObjScale);
  // resobj->dDualObj /= (scaling->dObjScale * scaling->dObjScale);
  resobj->dDualityGap = resobj->dPrimalObj - resobj->dDualObj;
  resobj->dRelObjGap =
      fabs(resobj->dPrimalObj - resobj->dDualObj) /
      (1.0 + fabs(resobj->dPrimalObj) + fabs(resobj->dDualObj));

  // resobj->dPrimalObjAverage /= scaling->dObjScale * scaling->dObjScale;
  // resobj->dDualObjAverage /= scaling->dObjScale * scaling->dObjScale;
  resobj->dDualityGapAverage =
      resobj->dPrimalObjAverage - resobj->dDualObjAverage;
  resobj->dRelObjGapAverage =
      fabs(resobj->dPrimalObjAverage - resobj->dDualObjAverage) /
      (1.0 + fabs(resobj->dPrimalObjAverage) + fabs(resobj->dDualObjAverage));

#if problem_USE_TIMERS
  problem->dComputeResidualsTime += getTimeStamp() - dStartTime;
#endif
}

void PDHG_Init_Variables(CUPDLPwork *work) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;

  // cupdlp_zero(iterates->x, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->x->data, cupdlp_float, lp->nCols);

  // XXX: PDLP Does not project x0,  so we uncomment for 1-1 comparison

  PDHG_Project_Bounds(work, iterates->x->data);

  // cupdlp_zero(iterates->y, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->y->data, cupdlp_float, lp->nRows);

  // Ax(work, iterates->ax, iterates->x);
  // ATyCPU(work, iterates->aty, iterates->y);
  Ax(work, iterates->ax, iterates->x);
  ATy(work, iterates->aty, iterates->y);

  // cupdlp_zero(iterates->xSum, cupdlp_float, lp->nCols);
  // cupdlp_zero(iterates->ySum, cupdlp_float, lp->nRows);
  // cupdlp_zero(iterates->xAverage, cupdlp_float, lp->nCols);
  // cupdlp_zero(iterates->yAverage, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->xSum, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->ySum, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->xAverage->data, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->yAverage->data, cupdlp_float, lp->nRows);

  PDHG_Project_Bounds(work, iterates->xSum);
  PDHG_Project_Bounds(work, iterates->xAverage->data);

  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  CUPDLP_ZERO_VEC(iterates->xLastRestart, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->yLastRestart, cupdlp_float, lp->nRows);
}

/* TODO: this function seems considering
 *       l1 <= Ax <= u1
 *       l2 <=  x <= u2
 *       needs rewritten for current formulation
 */
void PDHG_Check_Data(CUPDLPwork *work) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;
  cupdlp_int nFreeCol = 0;
  cupdlp_int nFixedCol = 0;
  cupdlp_int nUpperCol = 0;
  cupdlp_int nLowerCol = 0;
  cupdlp_int nRangedCol = 0;
  cupdlp_int nFreeRow = 0;
  cupdlp_int nFixedRow = 0;
  cupdlp_int nUpperRow = 0;
  cupdlp_int nLowerRow = 0;
  cupdlp_int nRangedRow = 0;

  for (cupdlp_int iSeq = 0; iSeq < lp->nCols; ++iSeq) {
    cupdlp_bool hasLower = problem->lower[iSeq] > -INFINITY;
    cupdlp_bool hasUpper = problem->upper[iSeq] < +INFINITY;

    if (!hasLower && !hasUpper) {
      ++nFreeCol;
      cupdlp_printf("Warning: variable %d is free.", iSeq);
    }

    if (hasLower && hasUpper) {
      if (problem->lower[iSeq] == problem->upper[iSeq]) {
        ++nFixedCol;
        // cupdlp_printf( "Warning: variable %d is fixed.", iSeq);
      } else
        ++nRangedCol;
    }

    if (hasLower) {
      // XXX: uncommented for PDLP comparison
      // CUPDLP_ASSERT(iterates->x[iSeq] >= problem->lower[iSeq]);
      nLowerCol += !hasUpper;
    }

    if (hasUpper) {
      // XXX: uncommented for PDLP comparison
      // CUPDLP_ASSERT(iterates->x[iSeq] <= problem->upper[iSeq]);
      nUpperCol += !hasLower;
    }
  }

  for (cupdlp_int iSeq = lp->nCols; iSeq < lp->nCols; ++iSeq) {
    cupdlp_bool hasLower = problem->lower[iSeq] > -INFINITY;
    cupdlp_bool hasUpper = problem->upper[iSeq] < +INFINITY;

    if (!hasLower && !hasUpper) {
      ++nFreeRow;
      cupdlp_printf("Warning: row %d is free.", iSeq - lp->nCols);
    }

    if (hasLower && hasUpper) {
      if (problem->lower[iSeq] == problem->upper[iSeq])
        ++nFixedRow;
      else
        ++nRangedRow;
    }

    if (hasLower) {
      // CUPDLP_ASSERT(iterates->x[iSeq] >= problem->lower[iSeq]);
      nLowerRow += !hasUpper;
    }

    if (hasUpper) {
      // CUPDLP_ASSERT(iterates->x[iSeq] <= problem->upper[iSeq]);
      nUpperRow += !hasLower;
    }
  }

  for (cupdlp_int iRow = 0; iRow < lp->nRows; ++iRow) {
    CUPDLP_ASSERT(iterates->y->data[iRow] < +INFINITY);
    CUPDLP_ASSERT(iterates->y->data[iRow] > -INFINITY);
  }

  for (cupdlp_int iRow = 0; iRow < lp->nRows; ++iRow) {
    if (problem->data->csr_matrix->rowMatBeg[iRow + 1] -
            problem->data->csr_matrix->rowMatBeg[iRow] ==
        1) {
      cupdlp_printf("Warning: row %d is a singleton row.", iRow);
    }
  }

  CUPDLP_ASSERT(nRangedRow == 0);
  cupdlp_printf("nFreeCol  : %d\n", nFreeCol);
  cupdlp_printf("nFixedCol : %d\n", nFixedCol);
  cupdlp_printf("nRangedCol: %d\n", nRangedCol);
  cupdlp_printf("nLowerCol : %d\n", nLowerCol);
  cupdlp_printf("nUpperCol : %d\n", nUpperCol);
  cupdlp_printf("nFreeRow  : %d\n", nFreeRow);
  cupdlp_printf("nFixedRow : %d\n", nFixedRow);
  cupdlp_printf("nRangedRow: %d\n", nRangedRow);
  cupdlp_printf("nLowerRow : %d\n", nLowerRow);
  cupdlp_printf("nUpperRow : %d\n", nUpperRow);

  // We need to test problems ranged row-bounds more carefully.
  CUPDLP_ASSERT(nRangedRow == 0);
}

cupdlp_bool PDHG_Check_Termination(CUPDLPwork *pdhg, int bool_print) {
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPscaling *scaling = pdhg->scaling;
#if PDHG_DISPLAY_TERMINATION_CHECK
  // todo, check, is it correct
  if (bool_print) {
    cupdlp_printf(
        "Termination check: %e|%e  %e|%e  %e|%e\n", resobj->dPrimalFeas,
        settings->dPrimalTol * (1.0 + scaling->dNormRhs), resobj->dDualFeas,
        settings->dDualTol * (1.0 + scaling->dNormCost), resobj->dRelObjGap,
        settings->dGapTol);
  }

#endif
  int bool_pass =
      ((resobj->dPrimalFeas <
        settings->dPrimalTol * (1.0 + scaling->dNormRhs)) &&
       (resobj->dDualFeas < settings->dDualTol * (1.0 + scaling->dNormCost)) &&
       (resobj->dRelObjGap < settings->dGapTol));
  return bool_pass;
}

cupdlp_bool PDHG_Check_Termination_Average(CUPDLPwork *pdhg, int bool_print) {
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPscaling *scaling = pdhg->scaling;
#if PDHG_DISPLAY_TERMINATION_CHECK
  if (bool_print) {
    cupdlp_printf("Termination check: %e|%e  %e|%e  %e|%e\n",
                  resobj->dPrimalFeasAverage,
                  settings->dPrimalTol * (1.0 + scaling->dNormRhs),
                  resobj->dDualFeasAverage,
                  settings->dDualTol * (1.0 + scaling->dNormCost),
                  resobj->dRelObjGapAverage, settings->dGapTol);
  }
#endif
  int bool_pass = ((resobj->dPrimalFeasAverage <
                    settings->dPrimalTol * (1.0 + scaling->dNormRhs)) &&
                   (resobj->dDualFeasAverage <
                    settings->dDualTol * (1.0 + scaling->dNormCost)) &&
                   (resobj->dRelObjGapAverage < settings->dGapTol));
  return bool_pass;
}

void PDHG_Print_Header(CUPDLPwork *pdhg) {
  cupdlp_printf("%5s  %15s  %15s   %8s  %8s  %10s  %8s %7s\n", "Iter",
                "Primal.Obj", "Dual.Obj", "Gap", "Compl", "Primal.Inf",
                "Dual.Inf", "Time");
}

void PDHG_Print_Iter(CUPDLPwork *pdhg) {
  /* Format time as xxx.yy for < 1000s and as integer afterwards. */
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;
  char timeString[8];
  if (timers->dSolvingTime < 100.0)
    cupdlp_snprintf(timeString, 8, "%6.2fs", timers->dSolvingTime);
  else
    cupdlp_snprintf(timeString, 8, "%6ds", (cupdlp_int)timers->dSolvingTime);

  cupdlp_printf("%5d  %+15.8e  %+15.8e  %+8.2e  %8.2e  %10.2e  %8.2e %7s [L]\n",
                timers->nIter, resobj->dPrimalObj, resobj->dDualObj,
                resobj->dDualityGap, resobj->dComplementarity,
                resobj->dPrimalFeas, resobj->dDualFeas, timeString);
}

void PDHG_Print_Iter_Average(CUPDLPwork *pdhg) {
  /* Format time as xxx.yy for < 1000s and as integer afterwards. */
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;
  char timeString[8];
  if (timers->dSolvingTime < 100.0)
    cupdlp_snprintf(timeString, 8, "%6.2fs", timers->dSolvingTime);
  else
    cupdlp_snprintf(timeString, 8, "%6ds", (cupdlp_int)timers->dSolvingTime);

  cupdlp_printf("%5d  %+15.8e  %+15.8e  %+8.2e  %8.2e  %10.2e  %8.2e %7s [A]\n",
                timers->nIter, resobj->dPrimalObjAverage,
                resobj->dDualObjAverage, resobj->dDualityGapAverage,
                resobj->dComplementarityAverage, resobj->dPrimalFeasAverage,
                resobj->dDualFeasAverage, timeString);
}

void PDHG_Compute_SolvingTime(CUPDLPwork *pdhg) {
  CUPDLPtimers *timers = pdhg->timers;
  timers->dSolvingTime = getTimeStamp() - timers->dSolvingBeg;
}

cupdlp_retcode PDHG_Solve(CUPDLPwork *pdhg) {
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPtimers *timers = pdhg->timers;

  timers->dSolvingBeg = getTimeStamp();

  PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDHG_Init_Step_Sizes(pdhg));

  PDHG_Init_Variables(pdhg);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);

  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter) {
    PDHG_Compute_SolvingTime(pdhg);
#if CUPDLP_DUMP_ITERATES_STATS & CUPDLP_DEBUG
    PDHG_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
#endif
    if (bool_checking) {
      PDHG_Compute_Average_Iterate(pdhg);
      PDHG_Compute_Residuals(pdhg);
      if (bool_print) {
        PDHG_Print_Header(pdhg);
        PDHG_Print_Iter(pdhg);
        PDHG_Print_Iter_Average(pdhg);
      }

      if (PDHG_Check_Termination(pdhg, bool_print)) {
        cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (PDHG_Check_Termination_Average(pdhg, bool_print)) {
        cupdlp_printf("Optimal average solution.\n");

        CUPDLP_COPY_VEC(iterates->x->data, iterates->xAverage->data,
                        cupdlp_float, problem->nCols);
        CUPDLP_COPY_VEC(iterates->y->data, iterates->yAverage->data,
                        cupdlp_float, problem->nRows);

        resobj->termIterate = AVERAGE_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (timers->dSolvingTime > settings->dTimeLim) {
        cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter == (settings->nIterLim - 1)) {
        cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      PDHG_Restart_Iterate(pdhg);
    }
    CUPDLP_CALL(PDHG_Update_Iterate(pdhg));
  }
  // print at last
  PDHG_Print_Header(pdhg);
  PDHG_Print_Iter(pdhg);
  PDHG_Print_Iter_Average(pdhg);

#if PDHG_USE_TIMERS
  cupdlp_printf("Timing information:\n");
  // cupdlp_printf("%20s %e in %d iterations\n", "Total solver time",
  //               timers->dSolvingTime, timers->nIter);
  cupdlp_printf(
      "%20s %e in %d iterations\n", "Total solver time",
      timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
      timers->nIter);
  cupdlp_printf("%20s %e in %d iterations\n", "Solve time",
                timers->dSolvingTime, timers->nIter);
  cupdlp_printf("%20s %e \n", "Iters per sec",
                timers->nIter / timers->dSolvingTime);
  cupdlp_printf("%20s %e\n", "Scaling time", timers->dScalingTime);
  cupdlp_printf("%20s %e\n", "Presolve time", timers->dPresolveTime);
  cupdlp_printf("%20s %e in %d calls\n", "Ax", timers->dAxTime,
                timers->nAxCalls);
  cupdlp_printf("%20s %e in %d calls\n", "Aty", timers->dAtyTime,
                timers->nAtyCalls);
  cupdlp_printf("%20s %e in %d calls\n", "ComputeResiduals",
                timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
  cupdlp_printf("%20s %e in %d calls\n", "UpdateIterates",
                timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
#endif

#ifndef CUPDLP_CPU
  cupdlp_printf("GPU Timing information:\n");
  cupdlp_printf("%20s %e\n", "CudaPrepare", timers->CudaPrepareTime);
  cupdlp_printf("%20s %e\n", "Alloc&CopyMatToDevice",
                timers->AllocMem_CopyMatToDeviceTime);
  cupdlp_printf("%20s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
  cupdlp_printf("%20s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
  cupdlp_printf("%20s %e\n", "CopyVecToHost", timers->CopyVecToHostTime);
#endif

exit_cleanup:
  return retcode;
}

void PDHG_PostSolve(CUPDLPwork *pdhg, cupdlp_int nCols_origin,
                    cupdlp_int *constraint_new_idx, cupdlp_float *x_origin,
                    cupdlp_float *y_origin) {
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPscaling *scaling = pdhg->scaling;

  // unscale
  if (scaling->ifScaled) {
    cupdlp_ediv(iterates->x->data, pdhg->colScale, problem->nCols);
    cupdlp_ediv(iterates->y->data, pdhg->rowScale, problem->nRows);

    // cupdlp_ediv(iterates->x->data, scaling->colScale_gpu, problem->nCols);
    // cupdlp_ediv(iterates->y->data, scaling->rowScale_gpu, problem->nRows);
  }

  // extract x from (x, z)
  CUPDLP_COPY_VEC(x_origin, iterates->x->data, cupdlp_float, nCols_origin);

  cupdlp_float *ytmp =
      (cupdlp_float *)cupdlp_malloc(problem->nRows * sizeof(cupdlp_float));
  CUPDLP_COPY_VEC(ytmp, iterates->y->data, cupdlp_float, problem->nRows);
  // un-permute y
  for (int i = 0; i < problem->nRows; i++) {
    y_origin[i] = ytmp[constraint_new_idx[i]];
  }
  cupdlp_free(ytmp);
}

cupdlp_retcode LP_SolvePDHG(CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam,
                            cupdlp_int *intParam,
                            cupdlp_bool *ifChangeFloatParam,
                            cupdlp_float *floatParam, char *fp,
                            cupdlp_float *x_origin, cupdlp_int nCols_origin,
                            cupdlp_float *y_origin, cupdlp_bool ifSaveSol,
                            cupdlp_int *constraint_new_idx) {
  cupdlp_retcode retcode = RETCODE_OK;

  PDHG_PrintHugeCUPDHG();

  CUPDLP_CALL(PDHG_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                ifChangeFloatParam, floatParam));

  CUPDLP_CALL(PDHG_Solve(pdhg));

  PDHG_PostSolve(pdhg, nCols_origin, constraint_new_idx, x_origin, y_origin);

  if (fp)
    writeJson(fp, pdhg, x_origin, nCols_origin, y_origin, pdhg->problem->nRows,
	      ifSaveSol);

exit_cleanup:
  PDHG_Destroy(&pdhg);
  return retcode;
}
