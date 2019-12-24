/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexAnalysis.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include <cmath>
//#include <cstdio>
#include "HConfig.h"
#include "simplex/HighsSimplexAnalysis.h"
#include "simplex/HFactor.h"

void HighsSimplexAnalysis::setup(const HighsLp& lp, const HighsOptions& options) {
  // Copy Problem size
  numRow = lp.numRow_;
  numCol = lp.numCol_;
  messaging(options.logfile, options.output, options.message_level);
  col_aq_density = 0;
  row_ep_density = 0;
  row_ap_density = 0;
  row_DSE_density = 0;
}

void HighsSimplexAnalysis::messaging(FILE* logfile_, FILE* output_, const int message_level_) {
  logfile = logfile_;
  output = output_;
  message_level = message_level_;
}

void HighsSimplexAnalysis::updateOperationResultDensity(const double local_density, double& density) {
  density = (1 - running_average_multiplier) * density +
    running_average_multiplier * local_density;
}

/*
void HighsSimplexAnalysis::equalDensity(const double density0, const double density1) {
  const double delta_density = std::fabs(density1-density0);
  if (delta_density>1e-15) {
    printf("ERROR: Difference %g in density0 - %g and density1 = %g\n", delta_density, density0, density1);
  }
}
*/

void HighsSimplexAnalysis::iterationReport() {
  const int iteration_count_difference = major_iteration_number -
    previous_iteration_report_header_iteration_count;
  const bool header = previous_iteration_report_header_iteration_count < 0
    || iteration_count_difference > 10;
  if (header) {
    iterationReportFull(header);
    previous_iteration_report_header_iteration_count = major_iteration_number;
  }
  iterationReportFull(false);
}

void HighsSimplexAnalysis::iterationReportFull(const bool header) {
  if (header) {
    iterationReportIterationAndPhase(ML_DETAILED, true);
    iterationReportDualObjective(ML_DETAILED, true);
#ifdef HiGHSDEV
    iterationReportIterationData(ML_DETAILED, true);
    iterationReportDensity(ML_DETAILED, true);
    HighsPrintMessage(output, message_level, ML_DETAILED, " FreeLsZ");
#endif
    HighsPrintMessage(output, message_level, ML_DETAILED, "\n");
  } else {
    iterationReportIterationAndPhase(ML_DETAILED, false);
    iterationReportDualObjective(ML_DETAILED, false);
#ifdef HiGHSDEV
    iterationReportIterationData(ML_DETAILED, false);
    iterationReportDensity(ML_DETAILED, false);
    HighsPrintMessage(output, message_level, ML_DETAILED, " %7d", freelist_size);
#endif
    HighsPrintMessage(output, message_level, ML_DETAILED, "\n");
  }
}

void HighsSimplexAnalysis::iterationReportIterationAndPhase(const int iterate_log_level,
                                             const bool header) {
  if (header) {
    HighsPrintMessage(output, message_level, iterate_log_level,
		      " Iteration Ph");
  } else {
    HighsPrintMessage(output, message_level, iterate_log_level,
		      " %9d %2d", major_iteration_number, solve_phase);
  }
}
void HighsSimplexAnalysis::iterationReportDualObjective(const int iterate_log_level, const bool header) {
  if (header) {
    HighsPrintMessage(output, message_level, iterate_log_level, "        DualObjective");
  } else {
    HighsPrintMessage(output, message_level, iterate_log_level, " %20.10e", objective_value);
  }
}

void HighsSimplexAnalysis::iterationReportIterationData(const int iterate_log_level, const bool header) {
  if (header) {
    HighsPrintMessage(output, message_level, iterate_log_level,
                      " Inv       NumCk     LvR     LvC     EnC        DlPr    "
                      "    ThDu        ThPr          Aa");
  } else {
    HighsPrintMessage(output, message_level, iterate_log_level,
                      " %3d %11.4g %7d %7d %7d %11.4g %11.4g %11.4g %11.4g",
                      invert_hint, numerical_trouble, pivotal_row_index, leaving_variable, entering_variable,
                      primal_step, dual_step, reduced_rhs_value, pivot_value_from_column);
  }
}

void HighsSimplexAnalysis::iterationReportDensity(const int iterate_log_level, const bool header) {
  const bool rp_dual_steepest_edge =
      edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE;
  if (header) {
    HighsPrintMessage(output, message_level, iterate_log_level, " C_Aq R_Ep R_Ap");
    if (rp_dual_steepest_edge) {
      HighsPrintMessage(output, message_level, iterate_log_level, "  DSE");
    } else {
      HighsPrintMessage(output, message_level, iterate_log_level, "     ");
    }
  } else {
    const int l10ColDse = intLog10(col_aq_density);
    const int l10REpDse = intLog10(row_ep_density);
    const int l10RapDse = intLog10(row_ap_density);
    HighsPrintMessage(output, message_level, iterate_log_level,
		      " %4d %4d %4d", l10ColDse, l10REpDse, l10RapDse);
    if (rp_dual_steepest_edge) {
      const int l10DseDse = intLog10(row_DSE_density);
      HighsPrintMessage(output, message_level, iterate_log_level,
			" %4d", l10DseDse);
    } else {
      HighsPrintMessage(output, message_level, iterate_log_level,
			"     ");
    }
  }
}


void HighsSimplexAnalysis::initialise(const int simplex_iteration_count) {
  AnIterIt0 = simplex_iteration_count;
  timer_.resetHighsTimer();
  AnIterCostlyDseFq = 0;
#ifdef HiGHSDEV
  AnIterPrevRpNumCostlyDseIt = 0;
  AnIterPrevIt = 0;
  AnIterOpRec* AnIter;
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_BTRAN];
  AnIter->AnIterOpName = "Btran";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_PRICE];
  AnIter->AnIterOpName = "Price";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_FTRAN];
  AnIter->AnIterOpName = "Ftran";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_FTRAN_BFRT];
  AnIter->AnIterOpName = "Ftran BFRT";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_FTRAN_DSE];
  AnIter->AnIterOpName = "Ftran_DSE";
  for (int k = 0; k < NUM_ANALYSIS_OPERATION_TYPE; k++) {
    AnIter = &AnIterOp[k];
    AnIter->AnIterOpLog10RsDsty = 0;
    AnIter->AnIterOpSuLog10RsDsty = 0;
    if (k == ANALYSIS_OPERATION_TYPE_PRICE) {
      AnIter->AnIterOpHyperCANCEL = 1.0;
      AnIter->AnIterOpHyperTRAN = 1.0;
      AnIter->AnIterOpRsDim = numCol;
    } else {
      if (k == ANALYSIS_OPERATION_TYPE_BTRAN) {
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
  int last_invert_hint = INVERT_HINT_Count - 1;
  for (int k = 1; k <= last_invert_hint; k++) AnIterNumInvert[k] = 0;
  AnIterNumPrDgnIt = 0;
  AnIterNumDuDgnIt = 0;
  AnIterNumColPrice = 0;
  AnIterNumRowPrice = 0;
  AnIterNumRowPriceWSw = 0;
  AnIterNumRowPriceUltra = 0;
  int last_dual_edge_weight_mode = (int)DualEdgeWeightMode::STEEPEST_EDGE;
  for (int k = 0; k <= last_dual_edge_weight_mode; k++) {
    AnIterNumEdWtIt[k] = 0;
  }
  AnIterNumCostlyDseIt = 0;
  AnIterTraceNumRec = 0;
  AnIterTraceIterDl = 1;
  AnIterTraceRec* lcAnIter = &AnIterTrace[0];
  lcAnIter->AnIterTraceIter = AnIterIt0;
  lcAnIter->AnIterTraceTime = timer_.getTime();
#endif
}

  /*
void HighsSimplexAnalysis::iterationAnalysis() {
  // Possibly report on the iteration
  iterationReport();

  // Possibly switch from DSE to Devex
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    bool switch_to_devex = false;
    // Firstly consider switching on the basis of NLA cost
    double AnIterCostlyDseMeasureDen;
    AnIterCostlyDseMeasureDen =
      max(max(row_ep_density, col_aq_density), row_ap_density);
    if (AnIterCostlyDseMeasureDen > 0) {
      AnIterCostlyDseMeasure = row_DSE_density / AnIterCostlyDseMeasureDen;
      AnIterCostlyDseMeasure = AnIterCostlyDseMeasure * AnIterCostlyDseMeasure;
    } else {
      AnIterCostlyDseMeasure = 0;
    }
    bool CostlyDseIt = AnIterCostlyDseMeasure > AnIterCostlyDseMeasureLimit &&
      row_DSE_density > AnIterCostlyDseMnDensity;
    AnIterCostlyDseFq = (1 - runningAverageMu) * AnIterCostlyDseFq;
    if (CostlyDseIt) {
      AnIterNumCostlyDseIt++;
      AnIterCostlyDseFq += runningAverageMu * 1.0;
      int lcNumIter = workHMO.scaled_solution_params_.simplex_iteration_count - AnIterIt0;
      // Switch to Devex if at least 5% of the (at least) 0.1NumTot iterations have been costly
      switch_to_devex = allow_dual_steepest_edge_to_devex_switch &&
	(AnIterNumCostlyDseIt > lcNumIter * AnIterFracNumCostlyDseItbfSw) &&
	(lcNumIter > AnIterFracNumTot_ItBfSw * solver_num_tot);
#ifdef HiGHSDEV
      if (switch_to_devex) {
	HighsLogMessage(logfile, HighsMessageType::INFO,
			"Switch from DSE to Devex after %d costly DSE iterations of %d: "
			"C_Aq_Dsty = %11.4g; R_Ep_Dsty = %11.4g; DSE_Dsty = %11.4g",
			AnIterNumCostlyDseIt, lcNumIter, row_DSE_density, row_ep_density,
			col_aq_density);
      }
#endif
    }
    if (!switch_to_devex) {
      // Secondly consider switching on the basis of weight accuracy
      double dse_weight_error_measure =
	average_log_low_dual_steepest_edge_weight_error +
	average_log_high_dual_steepest_edge_weight_error;
      double dse_weight_error_threshhold =
	dual_steepest_edge_weight_log_error_threshhold;
      switch_to_devex = allow_dual_steepest_edge_to_devex_switch &&
	dse_weight_error_measure > dse_weight_error_threshhold;
#ifdef HiGHSDEV
      if (switch_to_devex) {
	HighsLogMessage(workHMO.options_.logfile, HighsMessageType::INFO,
			"Switch from DSE to Devex with log error measure of %g > %g = threshhold",
			dse_weight_error_measure, dse_weight_error_threshhold);
      }
#endif
    }
    if (switch_to_devex) {
      dual_edge_weight_mode = DualEdgeWeightMode::DEVEX;
      // Zero the number of Devex frameworks used and set up the first one
      num_devex_framework = 0;
      devex_index.assign(solver_num_tot, 0);
      workHMO.simplex_info_.devex_index_.assign(solver_num_tot, 0);
      initialiseDevexFramework();
    }
  }

#ifdef HiGHSDEV
  int AnIterCuIt = workHMO.scaled_solution_params_.simplex_iteration_count;
  bool iterLg = AnIterCuIt % 100 == 0;
  iterLg = false;
  if (iterLg) {
    int lc_NumCostlyDseIt = AnIterNumCostlyDseIt - AnIterPrevRpNumCostlyDseIt;
    AnIterPrevRpNumCostlyDseIt = AnIterNumCostlyDseIt;
    printf("Iter %10d: ", AnIterCuIt);
    iterationReportDensity(ML_MINIMAL, true);
    iterationReportDensity(ML_MINIMAL, false);
    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      int lc_pct = (100 * AnIterNumCostlyDseIt) / (AnIterCuIt - AnIterIt0);
      printf("| Fq = %4.2f; Su =%5d (%3d%%)", AnIterCostlyDseFq,
             AnIterNumCostlyDseIt, lc_pct);

      if (lc_NumCostlyDseIt > 0) printf("; LcNum =%3d", lc_NumCostlyDseIt);
    }
    printf("\n");
  }

  for (int k = 0; k < NUM_ANALYSIS_OPERATION_TYPE; k++) {
    AnIterOpRec* lcAnIterOp = &AnIterOp[k];
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
    AnIterNumEdWtIt[(int)dual_edge_weight_mode] += (AnIterCuIt - AnIterPrevIt);

  AnIterTraceRec* lcAnIter = &AnIterTrace[AnIterTraceNumRec];
  //  if (workHMO.scaled_solution_params_.simplex_iteration_count ==
  //  AnIterTraceIterRec[AnIterTraceNumRec]+AnIterTraceIterDl) {
  if (workHMO.scaled_solution_params_.simplex_iteration_count ==
      lcAnIter->AnIterTraceIter + AnIterTraceIterDl) {
    if (AnIterTraceNumRec == AN_ITER_TRACE_MX_NUM_REC) {
      for (int rec = 1; rec <= AN_ITER_TRACE_MX_NUM_REC / 2; rec++)
        AnIterTrace[rec] = AnIterTrace[2 * rec];
      AnIterTraceNumRec = AnIterTraceNumRec / 2;
      AnIterTraceIterDl = AnIterTraceIterDl * 2;
    } else {
      AnIterTraceNumRec++;
      lcAnIter = &AnIterTrace[AnIterTraceNumRec];
      lcAnIter->AnIterTraceIter = workHMO.scaled_solution_params_.simplex_iteration_count;
      lcAnIter->AnIterTraceTime = workHMO.timer_.getTime();
      lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_BTRAN] = row_ep_density;
      lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_PRICE] = row_ap_density;
      lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN] = col_aq_density;
      lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN_BFRT] = col_aq_density;
      if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
        lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] = row_DSE_density;
        lcAnIter->AnIterTraceAux0 = AnIterCostlyDseMeasure;
      } else {
        lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] = 0;
        lcAnIter->AnIterTraceAux0 = 0;
      }
      lcAnIter->AnIterTrace_dual_edge_weight_mode = (int)dual_edge_weight_mode;
    }
  }
  AnIterPrevIt = AnIterCuIt;
#endif
}

#ifdef HiGHSDEV
void HighsSimplexAnalysis::iterateOpRecBf(int opTy, HVector& vector, double hist_dsty) {
  AnIterOpRec* AnIter = &AnIterOp[opTy];
  AnIter->AnIterOpNumCa++;
  double curr_dsty = 1.0 * vector.count / solver_num_row;
  //  printf("%10s: %g<= %g;  %g<= %g\n", AnIter->AnIterOpName.c_str(),
  //	 curr_dsty, AnIter->AnIterOpHyperCANCEL,
  //	 hist_dsty, AnIter->AnIterOpHyperTRAN);
  if (curr_dsty <= AnIter->AnIterOpHyperCANCEL &&
      hist_dsty <= AnIter->AnIterOpHyperTRAN)
    AnIter->AnIterOpNumHyperOp++;
}

void HighsSimplexAnalysis::iterateOpRecAf(int opTy, HVector& vector) {
  AnIterOpRec* AnIter = &AnIterOp[opTy];
  double rsDsty = 1.0 * vector.count / AnIter->AnIterOpRsDim;
  if (rsDsty <= hyperRESULT) AnIter->AnIterOpNumHyperRs++;
  AnIter->AnIterOpRsMxNNZ = max(vector.count, AnIter->AnIterOpRsMxNNZ);
  if (rsDsty > 0) {
    AnIter->AnIterOpLog10RsDsty += log(rsDsty) / log(10.0);
  } else {
    */
    /*
    // TODO Investigate these zero norms
    double vectorNorm = 0;

    for (int index = 0; index < AnIter->AnIterOpRsDim; index++) {
      double vectorValue = vector.array[index];
      vectorNorm += vectorValue * vectorValue;
    }
    vectorNorm = sqrt(vectorNorm);
    printf("Strange: operation %s has result density = %g: ||vector|| = %g\n",
    AnIter->AnIterOpName.c_str(), rsDsty, vectorNorm);
    */
    /*
  }
}

void HighsSimplexAnalysis::iterationAnalysisReport() {
  HighsTimer& timer = workHMO.timer_;
  int AnIterNumIter = workHMO.scaled_solution_params_.simplex_iteration_count - AnIterIt0;
  printf("\nAnalysis of %d iterations (%d to %d)\n", AnIterNumIter,
         AnIterIt0 + 1, workHMO.scaled_solution_params_.simplex_iteration_count);
  if (AnIterNumIter <= 0) return;
  int lc_EdWtNumIter;
  lc_EdWtNumIter = AnIterNumEdWtIt[(int)DualEdgeWeightMode::STEEPEST_EDGE];
  if (lc_EdWtNumIter > 0)
    printf("DSE for %12d (%3d%%) iterations\n", lc_EdWtNumIter,
           (100 * lc_EdWtNumIter) / AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[(int)DualEdgeWeightMode::DEVEX];
  if (lc_EdWtNumIter > 0)
    printf("Dvx for %12d (%3d%%) iterations\n", lc_EdWtNumIter,
           (100 * lc_EdWtNumIter) / AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[(int)DualEdgeWeightMode::DANTZIG];
  if (lc_EdWtNumIter > 0)
    printf("Dan for %12d (%3d%%) iterations\n", lc_EdWtNumIter,
           (100 * lc_EdWtNumIter) / AnIterNumIter);
  printf("\n");
  for (int k = 0; k < NUM_ANALYSIS_OPERATION_TYPE; k++) {
    AnIterOpRec* AnIter = &AnIterOp[k];
    int lcNumCa = AnIter->AnIterOpSuNumCa;
    printf("\n%-9s performed %d times\n", AnIter->AnIterOpName.c_str(),
           AnIter->AnIterOpSuNumCa);
    if (lcNumCa > 0) {
      int lcHyperOp = AnIter->AnIterOpSuNumHyperOp;
      int lcHyperRs = AnIter->AnIterOpSuNumHyperRs;
      int pctHyperOp = (100 * lcHyperOp) / lcNumCa;
      int pctHyperRs = (100 * lcHyperRs) / lcNumCa;
      double lcRsDsty = pow(10.0, AnIter->AnIterOpSuLog10RsDsty / lcNumCa);
      int lcAnIterOpRsDim = AnIter->AnIterOpRsDim;
      int lcNumNNz = lcRsDsty * lcAnIterOpRsDim;
      int lcMxNNz = AnIter->AnIterOpRsMxNNZ;
      double lcMxNNzDsty = (1.0 * lcMxNNz) / AnIter->AnIterOpRsDim;
      printf("%12d hyper-sparse operations (%3d%%)\n", lcHyperOp, pctHyperOp);
      printf("%12d hyper-sparse results    (%3d%%)\n", lcHyperRs, pctHyperRs);
      printf("%12g density of result (%d / %d nonzeros)\n", lcRsDsty, lcNumNNz,
             lcAnIterOpRsDim);
      printf("%12g density of result with max (%d / %d) nonzeros\n",
             lcMxNNzDsty, lcMxNNz, lcAnIterOpRsDim);
    }
  }
  int NumInvert = 0;

  int last_invert_hint = INVERT_HINT_Count - 1;
  for (int k = 1; k <= last_invert_hint; k++) NumInvert += AnIterNumInvert[k];
  if (NumInvert > 0) {
    int lcNumInvert = 0;
    printf("\nInvert    performed %d times: average frequency = %d\n",
           NumInvert, AnIterNumIter / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_UPDATE_LIMIT_REACHED];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to update limit reached\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_SYNTHETIC_CLOCK_SAYS_INVERT];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to pseudo-clock\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_POSSIBLY_OPTIMAL];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to possibly optimal\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_POSSIBLY_PRIMAL_UNBOUNDED];
    if (lcNumInvert > 0)
      printf(
          "%12d (%3d%%) Invert operations due to possibly primal unbounded\n",
          lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_POSSIBLY_DUAL_UNBOUNDED];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to possibly dual unbounded\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_POSSIBLY_SINGULAR_BASIS];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to possibly singular basis\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert =
        AnIterNumInvert[INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX];
    if (lcNumInvert > 0)
      printf(
          "%12d (%3d%%) Invert operations due to primal infeasible in primal "
          "simplex\n",
          lcNumInvert, (100 * lcNumInvert) / NumInvert);
  }
  printf("\n%12d (%3d%%) primal degenerate iterations\n", AnIterNumPrDgnIt,
         (100 * AnIterNumPrDgnIt) / AnIterNumIter);
  printf("%12d (%3d%%)   dual degenerate iterations\n", AnIterNumDuDgnIt,
         (100 * AnIterNumDuDgnIt) / AnIterNumIter);
  int suPrice = AnIterNumColPrice + AnIterNumRowPrice + AnIterNumRowPriceWSw +
                AnIterNumRowPriceUltra;
  if (suPrice > 0) {
    printf("\n%12d Price operations:\n", suPrice);
    printf("%12d Col Price      (%3d%%)\n", AnIterNumColPrice,
           (100 * AnIterNumColPrice) / suPrice);
    printf("%12d Row Price      (%3d%%)\n", AnIterNumRowPrice,
           (100 * AnIterNumRowPrice) / suPrice);
    printf("%12d Row PriceWSw   (%3d%%)\n", AnIterNumRowPriceWSw,
           (100 * AnIterNumRowPriceWSw / suPrice));
    printf("%12d Row PriceUltra (%3d%%)\n", AnIterNumRowPriceUltra,
           (100 * AnIterNumRowPriceUltra / suPrice));
  }
  printf("\n%12d (%3d%%) costly DSE        iterations\n", AnIterNumCostlyDseIt,
         (100 * AnIterNumCostlyDseIt) / AnIterNumIter);

  //
  // Add a record for the final iterations: may end up with one more
  // than AN_ITER_TRACE_MX_NUM_REC records, so ensure that there is enough
  // space in the arrays
  //
  AnIterTraceNumRec++;
  AnIterTraceRec* lcAnIter;
  lcAnIter = &AnIterTrace[AnIterTraceNumRec];
  lcAnIter->AnIterTraceIter = workHMO.scaled_solution_params_.simplex_iteration_count;
  lcAnIter->AnIterTraceTime = timer.getTime();
  lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_BTRAN] = row_ep_density;
  lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_PRICE] = row_ap_density;
  lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN] = col_aq_density;
  lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN_BFRT] = col_aq_density;
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] = row_DSE_density;
    lcAnIter->AnIterTraceAux0 = AnIterCostlyDseMeasure;
  } else {
    lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] = 0;
    lcAnIter->AnIterTraceAux0 = 0;
  }
  lcAnIter->AnIterTrace_dual_edge_weight_mode = (int)dual_edge_weight_mode;

  if (AnIterTraceIterDl >= 100) {
    printf("\n Iteration speed analysis\n");
    lcAnIter = &AnIterTrace[0];
    int fmIter = lcAnIter->AnIterTraceIter;
    double fmTime = lcAnIter->AnIterTraceTime;
    printf(
        "        Iter (      FmIter:      ToIter)      Time      Iter/sec | "
        "C_Aq R_Ep R_Ap  DSE | "
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
      int lc_dual_edge_weight_mode =
          lcAnIter->AnIterTrace_dual_edge_weight_mode;
      int l10ColDse = intLog10(lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN]);
      int l10REpDse = intLog10(lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_BTRAN]);
      int l10RapDse = intLog10(lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_PRICE]);
      int l10DseDse = intLog10(lcAnIter->AnIterTraceDsty[ANALYSIS_OPERATION_TYPE_FTRAN_DSE]);
      int l10Aux0 = intLog10(lcAnIter->AnIterTraceAux0);
      std::string str_dual_edge_weight_mode;
      if (lc_dual_edge_weight_mode == (int)DualEdgeWeightMode::STEEPEST_EDGE)
        str_dual_edge_weight_mode = "DSE";
      else if (lc_dual_edge_weight_mode == (int)DualEdgeWeightMode::DEVEX)
        str_dual_edge_weight_mode = "Dvx";
      else if (lc_dual_edge_weight_mode == (int)DualEdgeWeightMode::DANTZIG)
        str_dual_edge_weight_mode = "Dan";
      else
        str_dual_edge_weight_mode = "XXX";
      printf("%12d (%12d:%12d) %9.4f  %12d | %4d %4d %4d %4d |  %3s | %4d\n",
             dlIter, fmIter, toIter, dlTime, iterSpeed, l10ColDse, l10REpDse,
             l10RapDse, l10DseDse, str_dual_edge_weight_mode.c_str(), l10Aux0);
      fmIter = toIter;
      fmTime = toTime;
    }
    printf("\n");
  }
}
#endif
    */
int HighsSimplexAnalysis::intLog10(const double v) {
  int intLog10V = -99;
  if (v > 0) intLog10V = log(v) / log(10.0);
  return intLog10V;
}

