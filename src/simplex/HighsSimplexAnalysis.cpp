/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
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
#include "simplex/HFactor.h"
#include "simplex/HighsSimplexAnalysis.h"

void HighsSimplexAnalysis::setup(const HighsLp& lp, const HighsOptions& options,
                                 const int simplex_iteration_count_) {
  // Copy Problem size
  numRow = lp.numRow_;
  numCol = lp.numCol_;
  numTot = numRow + numCol;
  // Copy tolerances from options
  allow_dual_steepest_edge_to_devex_switch =
      options.simplex_dual_edge_weight_strategy ==
      SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH;
  dual_steepest_edge_weight_log_error_threshhold =
      options.dual_steepest_edge_weight_log_error_threshhold;
  //
  AnIterIt0 = simplex_iteration_count_;
  AnIterCostlyDseFq = 0;
  AnIterPrevRpNumCostlyDseIt = 0;
  // Copy messaging parameter from options
  messaging(options.logfile, options.output, options.message_level);
  // Initialise the densities
  col_aq_density = 0;
  row_ep_density = 0;
  row_ap_density = 0;
  row_DSE_density = 0;
  col_BFRT_density = 0;
  primal_col_density = 0;
  // Set the row_dual_density to 1 since it's assumed all costs are at
  // least perturbed from zero, if not initially nonzero
  dual_col_density = 1;
  // Initialise the measures used to analyse accuracy of steepest edge weights
  //
  const int dual_edge_weight_strategy =
      options.simplex_dual_edge_weight_strategy;
  if (dual_edge_weight_strategy ==
          SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE ||
      dual_edge_weight_strategy ==
          SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL ||
      dual_edge_weight_strategy ==
          SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH) {
    // Initialise the measures used to analyse accuracy of steepest edge weights
    num_dual_steepest_edge_weight_check = 0;
    num_dual_steepest_edge_weight_reject = 0;
    num_wrong_low_dual_steepest_edge_weight = 0;
    num_wrong_high_dual_steepest_edge_weight = 0;
    average_frequency_low_dual_steepest_edge_weight = 0;
    average_frequency_high_dual_steepest_edge_weight = 0;
    average_log_low_dual_steepest_edge_weight_error = 0;
    average_log_high_dual_steepest_edge_weight_error = 0;
    max_average_frequency_low_dual_steepest_edge_weight = 0;
    max_average_frequency_high_dual_steepest_edge_weight = 0;
    max_sum_average_frequency_extreme_dual_steepest_edge_weight = 0;
    max_average_log_low_dual_steepest_edge_weight_error = 0;
    max_average_log_high_dual_steepest_edge_weight_error = 0;
    max_sum_average_log_extreme_dual_steepest_edge_weight_error = 0;
  }
  num_devex_framework = 0;

  num_iteration_report_since_last_header = -1;
  num_invert_report_since_last_header = -1;

  // Set following averages to illegal values so that first average is
  // set equal to first value
  average_num_threads = -1;
  average_fraction_of_possible_minor_iterations_performed = -1;
  sum_multi_chosen = 0;
  sum_multi_finished = 0;

#ifdef HiGHSDEV
  AnIterPrevIt = simplex_iteration_count_;
  timer_.resetHighsTimer();
  AnIterOpRec* AnIter;
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_BTRAN_EP];
  AnIter->AnIterOpName = "BTRAN e_p";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_PRICE_AP];
  AnIter->AnIterOpName = "PRICE a_p";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_BTRAN_FULL];
  AnIter->AnIterOpName = "BTRAN Full";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_PRICE_FULL];
  AnIter->AnIterOpName = "PRICE Full";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_FTRAN];
  AnIter->AnIterOpName = "FTRAN";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_FTRAN_BFRT];
  AnIter->AnIterOpName = "FTRAN BFRT";
  AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_FTRAN_DSE];
  AnIter->AnIterOpName = "FTRAN DSE";
  for (int k = 0; k < NUM_ANALYSIS_OPERATION_TYPE; k++) {
    AnIter = &AnIterOp[k];
    if ((k == ANALYSIS_OPERATION_TYPE_PRICE_AP) ||
        (k == ANALYSIS_OPERATION_TYPE_PRICE_FULL)) {
      AnIter->AnIterOpHyperCANCEL = 1.0;
      AnIter->AnIterOpHyperTRAN = 1.0;
      AnIter->AnIterOpRsDim = numCol;
    } else {
      if ((k == ANALYSIS_OPERATION_TYPE_BTRAN_EP) ||
          (k == ANALYSIS_OPERATION_TYPE_BTRAN_FULL)) {
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
    AnIter->AnIterOpSumLog10RsDensity = 0;
    initialiseValueDistribution(1e-8, 1.0, 10.0, AnIter->AnIterOp_density);
  }
  int last_invert_hint = INVERT_HINT_Count - 1;
  for (int k = 1; k <= last_invert_hint; k++) AnIterNumInvert[k] = 0;
  num_col_price = 0;
  num_row_price = 0;
  num_row_price_with_switch = 0;
  int last_dual_edge_weight_mode = (int)DualEdgeWeightMode::STEEPEST_EDGE;
  for (int k = 0; k <= last_dual_edge_weight_mode; k++) AnIterNumEdWtIt[k] = 0;
  AnIterNumCostlyDseIt = 0;
  AnIterTraceNumRec = 0;
  AnIterTraceIterDl = 1;
  AnIterTraceRec* lcAnIter = &AnIterTrace[0];
  lcAnIter->AnIterTraceIter = AnIterIt0;
  lcAnIter->AnIterTraceTime = timer_.getTime();

  initialiseValueDistribution(1e-16, 1e16, 10.0, primal_step_distribution);
  initialiseValueDistribution(1e-16, 1e16, 10.0, dual_step_distribution);
  initialiseValueDistribution(1e-8, 1e16, 10.0, pivot_distribution);
  initialiseValueDistribution(1e-16, 1.0, 10.0, numerical_trouble_distribution);
  initialiseValueDistribution(1e-16, 1e16, 10.0,
                              cost_perturbation1_distribution);
  initialiseValueDistribution(1e-16, 1e16, 10.0,
                              cost_perturbation2_distribution);
  initialiseValueDistribution(1e-8, 1.0, 10.0,
                              before_ftran_upper_sparse_density);
  initialiseValueDistribution(1e-8, 1.0, 10.0,
                              before_ftran_upper_hyper_density);
  initialiseValueDistribution(1e-8, 1.0, 10.0, ftran_upper_sparse_density);
  initialiseValueDistribution(1e-8, 1.0, 10.0, ftran_upper_hyper_density);
  initialiseValueDistribution(1e-16, 1e16, 10.0,
                              cleanup_dual_change_distribution);
  initialiseValueDistribution(1e-16, 1e16, 10.0,
                              cleanup_primal_step_distribution);
  initialiseValueDistribution(1e-16, 1e16, 10.0,
                              cleanup_dual_step_distribution);
  initialiseValueDistribution(1e-16, 1e16, 10.0,
                              cleanup_primal_change_distribution);
#endif
}

void HighsSimplexAnalysis::messaging(FILE* logfile_, FILE* output_,
                                     const int message_level_) {
  logfile = logfile_;
  output = output_;
  message_level = message_level_;
}

void HighsSimplexAnalysis::updateOperationResultDensity(
    const double local_density, double& density) {
  density = (1 - running_average_multiplier) * density +
            running_average_multiplier * local_density;
}

void HighsSimplexAnalysis::iterationReport() {
  if (!(iteration_report_message_level & message_level)) return;
  const bool header = (num_iteration_report_since_last_header < 0) ||
                      (num_iteration_report_since_last_header > 49);
  if (header) {
    iterationReport(header);
    num_iteration_report_since_last_header = 0;
  }
  iterationReport(false);
}

void HighsSimplexAnalysis::invertReport() {
  if (!(invert_report_message_level & message_level)) return;
  const bool header = (num_invert_report_since_last_header < 0) ||
                      (num_invert_report_since_last_header > 49) ||
                      (num_iteration_report_since_last_header >= 0);
  if (header) {
    invertReport(header);
    num_invert_report_since_last_header = 0;
  }
  invertReport(false);
  // Force an iteration report header if this is an INVERT report without an
  // invert_hint
  if (!invert_hint) num_iteration_report_since_last_header = -1;
}

void HighsSimplexAnalysis::invertReport(const bool header) {
  if (!(invert_report_message_level & message_level)) return;
  reportAlgorithmPhaseIterationObjective(header, invert_report_message_level);
#ifdef HiGHSDEV
  if (simplex_strategy == SIMPLEX_STRATEGY_DUAL_MULTI) {
    // Report on threads and PAMI
    reportThreads(header, invert_report_message_level);
    reportMulti(header, invert_report_message_level);
  }
  reportDensity(header, invert_report_message_level);
  reportInvert(header, invert_report_message_level);
  //  reportCondition(header, invert_report_message_level);
#endif
  reportInfeasibility(header, invert_report_message_level);
  HighsPrintMessage(output, message_level, invert_report_message_level, "\n");
  if (!header) num_invert_report_since_last_header++;
}

void HighsSimplexAnalysis::dualSteepestEdgeWeightError(
    const double computed_edge_weight, const double updated_edge_weight) {
  const bool accept_weight =
      updated_edge_weight >= accept_weight_threshhold * computed_edge_weight;
  int low_weight_error = 0;
  int high_weight_error = 0;
  double weight_error;
#ifdef HiGHSDEV
  string error_type = "  OK";
#endif
  num_dual_steepest_edge_weight_check++;
  if (!accept_weight) num_dual_steepest_edge_weight_reject++;
  if (updated_edge_weight < computed_edge_weight) {
    // Updated weight is low
    weight_error = computed_edge_weight / updated_edge_weight;
    if (weight_error > weight_error_threshhold) {
      low_weight_error = 1;
#ifdef HiGHSDEV
      error_type = " Low";
#endif
    }
    average_log_low_dual_steepest_edge_weight_error =
        0.99 * average_log_low_dual_steepest_edge_weight_error +
        0.01 * log(weight_error);
  } else {
    // Updated weight is correct or high
    weight_error = updated_edge_weight / computed_edge_weight;
    if (weight_error > weight_error_threshhold) {
      high_weight_error = 1;
#ifdef HiGHSDEV
      error_type = "High";
#endif
    }
    average_log_high_dual_steepest_edge_weight_error =
        0.99 * average_log_high_dual_steepest_edge_weight_error +
        0.01 * log(weight_error);
  }
  average_frequency_low_dual_steepest_edge_weight =
      0.99 * average_frequency_low_dual_steepest_edge_weight +
      0.01 * low_weight_error;
  average_frequency_high_dual_steepest_edge_weight =
      0.99 * average_frequency_high_dual_steepest_edge_weight +
      0.01 * high_weight_error;
  max_average_frequency_low_dual_steepest_edge_weight =
      max(max_average_frequency_low_dual_steepest_edge_weight,
          average_frequency_low_dual_steepest_edge_weight);
  max_average_frequency_high_dual_steepest_edge_weight =
      max(max_average_frequency_high_dual_steepest_edge_weight,
          average_frequency_high_dual_steepest_edge_weight);
  max_sum_average_frequency_extreme_dual_steepest_edge_weight =
      max(max_sum_average_frequency_extreme_dual_steepest_edge_weight,
          average_frequency_low_dual_steepest_edge_weight +
              average_frequency_high_dual_steepest_edge_weight);
  max_average_log_low_dual_steepest_edge_weight_error =
      max(max_average_log_low_dual_steepest_edge_weight_error,
          average_log_low_dual_steepest_edge_weight_error);
  max_average_log_high_dual_steepest_edge_weight_error =
      max(max_average_log_high_dual_steepest_edge_weight_error,
          average_log_high_dual_steepest_edge_weight_error);
  max_sum_average_log_extreme_dual_steepest_edge_weight_error =
      max(max_sum_average_log_extreme_dual_steepest_edge_weight_error,
          average_log_low_dual_steepest_edge_weight_error +
              average_log_high_dual_steepest_edge_weight_error);
#ifdef HiGHSDEV
  const bool report_weight_error = false;
  if (report_weight_error && weight_error > 0.5 * weight_error_threshhold) {
    printf(
        "DSE Wt Ck |%8d| OK = %1d (%4d / %6d) (c %10.4g, u %10.4g, er %10.4g - "
        "%s): Low (Fq %10.4g, Er %10.4g); High (Fq%10.4g, Er%10.4g) | %10.4g "
        "%10.4g %10.4g %10.4g %10.4g %10.4g\n",
        simplex_iteration_count, accept_weight,
        num_dual_steepest_edge_weight_check,
        num_dual_steepest_edge_weight_reject, computed_edge_weight,
        updated_edge_weight, weight_error, error_type.c_str(),
        average_frequency_low_dual_steepest_edge_weight,
        average_log_low_dual_steepest_edge_weight_error,
        average_frequency_high_dual_steepest_edge_weight,
        average_log_high_dual_steepest_edge_weight_error,
        max_average_frequency_low_dual_steepest_edge_weight,
        max_average_frequency_high_dual_steepest_edge_weight,
        max_sum_average_frequency_extreme_dual_steepest_edge_weight,
        max_average_log_low_dual_steepest_edge_weight_error,
        max_average_log_high_dual_steepest_edge_weight_error,
        max_sum_average_log_extreme_dual_steepest_edge_weight_error);
  }
#endif
}

bool HighsSimplexAnalysis::switchToDevex() {
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
  AnIterCostlyDseFq = (1 - running_average_multiplier) * AnIterCostlyDseFq;
  if (CostlyDseIt) {
    AnIterNumCostlyDseIt++;
    AnIterCostlyDseFq += running_average_multiplier * 1.0;
    int lcNumIter = simplex_iteration_count - AnIterIt0;
    // Switch to Devex if at least 5% of the (at least) 0.1NumTot iterations
    // have been costly
    switch_to_devex =
        allow_dual_steepest_edge_to_devex_switch &&
        (AnIterNumCostlyDseIt > lcNumIter * AnIterFracNumCostlyDseItbfSw) &&
        (lcNumIter > AnIterFracNumTot_ItBfSw * numTot);
#ifdef HiGHSDEV
    if (switch_to_devex) {
      HighsLogMessage(
          logfile, HighsMessageType::INFO,
          "Switch from DSE to Devex after %d costly DSE iterations of %d: "
          "C_Aq_Density = %11.4g; R_Ep_Density = %11.4g; R_Ap_Density = "
          "%11.4g; DSE_Density = %11.4g",
          AnIterNumCostlyDseIt, lcNumIter, col_aq_density, row_ep_density,
          row_ap_density, row_DSE_density);
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
      HighsLogMessage(logfile, HighsMessageType::INFO,
                      "Switch from DSE to Devex with log error measure of %g > "
                      "%g = threshhold",
                      dse_weight_error_measure, dse_weight_error_threshhold);
    }
#endif
  }
  return switch_to_devex;
}

#ifdef HiGHSDEV
void HighsSimplexAnalysis::iterationRecord() {
  int AnIterCuIt = simplex_iteration_count;
  if (invert_hint > 0) AnIterNumInvert[invert_hint]++;
  if (AnIterCuIt > AnIterPrevIt)
    AnIterNumEdWtIt[(int)edge_weight_mode] += (AnIterCuIt - AnIterPrevIt);

  AnIterTraceRec& lcAnIter = AnIterTrace[AnIterTraceNumRec];
  //  if (simplex_iteration_count ==
  //  AnIterTraceIterRec[AnIterTraceNumRec]+AnIterTraceIterDl) {
  if (simplex_iteration_count == lcAnIter.AnIterTraceIter + AnIterTraceIterDl) {
    if (AnIterTraceNumRec == AN_ITER_TRACE_MX_NUM_REC) {
      for (int rec = 1; rec <= AN_ITER_TRACE_MX_NUM_REC / 2; rec++)
        AnIterTrace[rec] = AnIterTrace[2 * rec];
      AnIterTraceNumRec = AnIterTraceNumRec / 2;
      AnIterTraceIterDl = AnIterTraceIterDl * 2;
    } else {
      AnIterTraceNumRec++;
      AnIterTraceRec& lcAnIter = AnIterTrace[AnIterTraceNumRec];
      lcAnIter.AnIterTraceIter = simplex_iteration_count;
      lcAnIter.AnIterTraceTime = timer_.getTime();
      if (average_fraction_of_possible_minor_iterations_performed > 0) {
        lcAnIter.AnIterTraceMulti =
            average_fraction_of_possible_minor_iterations_performed;
      } else {
        lcAnIter.AnIterTraceMulti = 0;
      }
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN] =
          col_aq_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_BTRAN_EP] =
          row_ep_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_PRICE_AP] =
          row_ap_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_BFRT] =
          col_aq_density;
      if (edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
        lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] =
            row_DSE_density;
        lcAnIter.AnIterTraceCostlyDse = AnIterCostlyDseMeasure;
      } else {
        lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] = 0;
        lcAnIter.AnIterTraceCostlyDse = 0;
      }
      lcAnIter.AnIterTrace_dual_edge_weight_mode = (int)edge_weight_mode;
    }
  }
  AnIterPrevIt = AnIterCuIt;
  updateValueDistribution(primal_step, cleanup_primal_step_distribution);
  updateValueDistribution(dual_step, cleanup_dual_step_distribution);
  updateValueDistribution(primal_step, primal_step_distribution);
  updateValueDistribution(dual_step, dual_step_distribution);
  updateValueDistribution(pivot_value_from_column, pivot_distribution);
  // Only update the distribution of legal values for
  // numerical_trouble. Illegal values are set in PAMI since it's not
  // known in minor iterations
  if (numerical_trouble >= 0)
    updateValueDistribution(numerical_trouble, numerical_trouble_distribution);
}

void HighsSimplexAnalysis::iterationRecordMajor() {
  sum_multi_chosen += multi_chosen;
  sum_multi_finished += multi_finished;
  assert(multi_chosen > 0);
  const double fraction_of_possible_minor_iterations_performed =
      1.0 * multi_finished / multi_chosen;
  if (average_fraction_of_possible_minor_iterations_performed < 0) {
    average_fraction_of_possible_minor_iterations_performed =
        fraction_of_possible_minor_iterations_performed;
  } else {
    average_fraction_of_possible_minor_iterations_performed =
        running_average_multiplier *
            fraction_of_possible_minor_iterations_performed +
        (1 - running_average_multiplier) *
            average_fraction_of_possible_minor_iterations_performed;
  }
  if (average_num_threads < 0) {
    average_num_threads = num_threads;
  } else {
    average_num_threads =
        running_average_multiplier * num_threads +
        (1 - running_average_multiplier) * average_num_threads;
  }
}

void HighsSimplexAnalysis::operationRecordBefore(
    const int operation_type, const HVector& vector,
    const double historical_density) {
  operationRecordBefore(operation_type, vector.count, historical_density);
}

void HighsSimplexAnalysis::operationRecordBefore(
    const int operation_type, const int current_count,
    const double historical_density) {
  double current_density = 1.0 * current_count / numRow;
  AnIterOpRec& AnIter = AnIterOp[operation_type];
  AnIter.AnIterOpNumCa++;
  if (current_density <= AnIter.AnIterOpHyperCANCEL &&
      historical_density <= AnIter.AnIterOpHyperTRAN)
    AnIter.AnIterOpNumHyperOp++;
}

void HighsSimplexAnalysis::operationRecordAfter(const int operation_type,
                                                const HVector& vector) {
  operationRecordAfter(operation_type, vector.count);
}

void HighsSimplexAnalysis::operationRecordAfter(const int operation_type,
                                                const int result_count) {
  AnIterOpRec& AnIter = AnIterOp[operation_type];
  const double result_density = 1.0 * result_count / AnIter.AnIterOpRsDim;
  if (result_density <= hyperRESULT) AnIter.AnIterOpNumHyperRs++;
  if (result_density > 0) {
    AnIter.AnIterOpSumLog10RsDensity += log(result_density) / log(10.0);
  } else {
    /*
    // TODO Investigate these zero norms
    double vectorNorm = 0;

    for (int index = 0; index < AnIter.AnIterOpRsDim; index++) {
      double vectorValue = vector.array[index];
      vectorNorm += vectorValue * vectorValue;
    }
    vectorNorm = sqrt(vectorNorm);
    printf("Strange: operation %s has result density = %g: ||vector|| = %g\n",
    AnIter.AnIterOpName.c_str(), result_density, vectorNorm);
    */
  }
  updateValueDistribution(result_density, AnIter.AnIterOp_density);
}

void HighsSimplexAnalysis::summaryReport() {
  int AnIterNumIter = simplex_iteration_count - AnIterIt0;
  if (AnIterNumIter <= 0) return;
  printf("\nAnalysis of %d iterations (%d to %d)\n", AnIterNumIter,
         AnIterIt0 + 1, simplex_iteration_count);
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
  for (int k = 0; k < NUM_ANALYSIS_OPERATION_TYPE; k++) {
    AnIterOpRec& AnIter = AnIterOp[k];
    int lcNumCa = AnIter.AnIterOpNumCa;
    printf("\n%-10s performed %d times\n", AnIter.AnIterOpName.c_str(),
           AnIter.AnIterOpNumCa);
    if (lcNumCa > 0) {
      int lcHyperOp = AnIter.AnIterOpNumHyperOp;
      int lcHyperRs = AnIter.AnIterOpNumHyperRs;
      int pctHyperOp = (100 * lcHyperOp) / lcNumCa;
      int pctHyperRs = (100 * lcHyperRs) / lcNumCa;
      double lcRsDensity =
          pow(10.0, AnIter.AnIterOpSumLog10RsDensity / lcNumCa);
      int lcAnIterOpRsDim = AnIter.AnIterOpRsDim;
      int lcNumNNz = lcRsDensity * lcAnIterOpRsDim;
      printf("%12d hyper-sparse operations (%3d%%)\n", lcHyperOp, pctHyperOp);
      printf("%12d hyper-sparse results    (%3d%%)\n", lcHyperRs, pctHyperRs);
      printf("%12g density of result (%d / %d nonzeros)\n", lcRsDensity,
             lcNumNNz, lcAnIterOpRsDim);
      printValueDistribution("density ", AnIter.AnIterOp_density,
                             AnIter.AnIterOpRsDim);
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
  int suPrice = num_col_price + num_row_price + num_row_price_with_switch;
  if (suPrice > 0) {
    printf("\n%12d Price operations:\n", suPrice);
    printf("%12d Col Price      (%3d%%)\n", num_col_price,
           (100 * num_col_price) / suPrice);
    printf("%12d Row Price      (%3d%%)\n", num_row_price,
           (100 * num_row_price) / suPrice);
    printf("%12d Row PriceWSw   (%3d%%)\n", num_row_price_with_switch,
           (100 * num_row_price_with_switch / suPrice));
  }
  printf("\n%12d (%3d%%) costly DSE        iterations\n", AnIterNumCostlyDseIt,
         (100 * AnIterNumCostlyDseIt) / AnIterNumIter);

  // Look for any Devex data to summarise
  if (num_devex_framework) {
    printf("\nDevex summary\n");
    printf("%12d Devex frameworks\n", num_devex_framework);
    printf(
        "%12d average number of iterations\n",
        AnIterNumEdWtIt[(int)DualEdgeWeightMode::DEVEX] / num_devex_framework);
  }
  // Look for any PAMI data to summarise
  if (sum_multi_chosen > 0) {
    const int pct_minor_iterations_performed =
        (100 * sum_multi_finished) / sum_multi_chosen;
    printf("\nPAMI summary: for average of %0.1g threads \n",
           average_num_threads);
    printf("%12d Major iterations\n", multi_iteration_count);
    printf("%12d Minor iterations\n", sum_multi_finished);
    printf(
        "%12d Total rows chosen: performed %3d%% of possible minor "
        "iterations\n\n",
        sum_multi_chosen, pct_minor_iterations_performed);
  }

  printf("\nCost perturbation summary\n");
  printValueDistribution("1 ", cost_perturbation1_distribution);
  printValueDistribution("2 ", cost_perturbation2_distribution);

  printf("\nFTRAN upper sparse summary - before\n");
  printValueDistribution("", before_ftran_upper_sparse_density, numRow);

  printf("\nFTRAN upper sparse summary - after\n");
  printValueDistribution("", ftran_upper_sparse_density, numRow);

  printf("\nFTRAN upper hyper-sparse summary - before\n");
  printValueDistribution("", before_ftran_upper_hyper_density, numRow);

  printf("\nFTRAN upper hyper-sparse summary - after\n");
  printValueDistribution("", ftran_upper_hyper_density, numRow);

  printf("\nPrimal step summary\n");
  printValueDistribution("", primal_step_distribution);

  printf("\nDual step summary\n");
  printValueDistribution("", dual_step_distribution);

  printf("\nPivot summary\n");
  printValueDistribution("", pivot_distribution);

  printf("\nNumerical trouble summary\n");
  printValueDistribution("", numerical_trouble_distribution);

  printf("\nCleanup dual change summary\n");
  printValueDistribution("dual ", cleanup_dual_change_distribution);

  printf("\nCleanup primal step summary\n");
  printValueDistribution("", cleanup_primal_step_distribution);

  printf("\nCleanup dual step summary\n");
  printValueDistribution("", cleanup_dual_step_distribution);

  printf("\nCleanup primal change summary\n");
  printValueDistribution("", cleanup_primal_change_distribution);

  if (AnIterTraceIterDl >= 100) {
    // Possibly (usually) add a temporary record for the final
    // iterations: may end up with one more than
    // AN_ITER_TRACE_MX_NUM_REC records, so ensure that there is
    // enough space in the arrays
    //
    const bool add_extra_record =
        simplex_iteration_count >
        AnIterTrace[AnIterTraceNumRec].AnIterTraceIter;
    if (add_extra_record) {
      AnIterTraceNumRec++;
      AnIterTraceRec& lcAnIter = AnIterTrace[AnIterTraceNumRec];
      lcAnIter.AnIterTraceIter = simplex_iteration_count;
      lcAnIter.AnIterTraceTime = timer_.getTime();
      if (average_fraction_of_possible_minor_iterations_performed > 0) {
        lcAnIter.AnIterTraceMulti =
            average_fraction_of_possible_minor_iterations_performed;
      } else {
        lcAnIter.AnIterTraceMulti = 0;
      }
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN] =
          col_aq_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_BTRAN_EP] =
          row_ep_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_PRICE_AP] =
          row_ap_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_BFRT] =
          col_aq_density;
      if (edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
        lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] =
            row_DSE_density;
        lcAnIter.AnIterTraceCostlyDse = AnIterCostlyDseMeasure;
      } else {
        lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] = 0;
        lcAnIter.AnIterTraceCostlyDse = 0;
      }
      lcAnIter.AnIterTrace_dual_edge_weight_mode = (int)edge_weight_mode;
    }
    // Determine whether the Multi and steepest edge columns should be reported
    double su_multi_values = 0;
    double su_dse_values = 0;
    for (int rec = 1; rec <= AnIterTraceNumRec; rec++) {
      AnIterTraceRec& lcAnIter = AnIterTrace[rec];
      su_multi_values += fabs(lcAnIter.AnIterTraceMulti);
      su_dse_values +=
          fabs(lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE]);
    }
    const bool report_multi = su_multi_values > 0;
    const bool rp_dual_steepest_edge = su_dse_values > 0;
    printf("\n Iteration speed analysis\n");
    AnIterTraceRec& lcAnIter = AnIterTrace[0];
    int fmIter = lcAnIter.AnIterTraceIter;
    double fmTime = lcAnIter.AnIterTraceTime;
    printf("        Iter (      FmIter:      ToIter)      Time      Iter/sec ");
    if (report_multi) printf("| PAMI ");
    printf("| C_Aq R_Ep R_Ap ");
    if (rp_dual_steepest_edge) printf(" DSE ");
    printf("| EdWt ");
    if (rp_dual_steepest_edge) {
      printf("| CostlyDse\n");
    } else {
      printf("\n");
    }

    for (int rec = 1; rec <= AnIterTraceNumRec; rec++) {
      AnIterTraceRec& lcAnIter = AnIterTrace[rec];
      int toIter = lcAnIter.AnIterTraceIter;
      double toTime = lcAnIter.AnIterTraceTime;
      int dlIter = toIter - fmIter;
      if (rec < AnIterTraceNumRec && dlIter != AnIterTraceIterDl)
        printf("STRANGE: %d = dlIter != AnIterTraceIterDl = %d\n", dlIter,
               AnIterTraceIterDl);
      double dlTime = toTime - fmTime;
      int iterSpeed = 0;
      if (dlTime > 0) iterSpeed = dlIter / dlTime;
      int lc_dual_edge_weight_mode = lcAnIter.AnIterTrace_dual_edge_weight_mode;
      std::string str_dual_edge_weight_mode;
      if (lc_dual_edge_weight_mode == (int)DualEdgeWeightMode::STEEPEST_EDGE)
        str_dual_edge_weight_mode = "DSE";
      else if (lc_dual_edge_weight_mode == (int)DualEdgeWeightMode::DEVEX)
        str_dual_edge_weight_mode = "Dvx";
      else if (lc_dual_edge_weight_mode == (int)DualEdgeWeightMode::DANTZIG)
        str_dual_edge_weight_mode = "Dan";
      else
        str_dual_edge_weight_mode = "XXX";
      printf("%12d (%12d:%12d) %9.4f  %12d ", dlIter, fmIter, toIter, dlTime,
             iterSpeed);
      if (report_multi) {
        const int pct = (100 * lcAnIter.AnIterTraceMulti);
        printf("|  %3d ", pct);
      }
      printf("|");
      reportOneDensity(
          lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN]);
      reportOneDensity(
          lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_BTRAN_EP]);
      reportOneDensity(
          lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_PRICE_AP]);
      double use_row_DSE_density;
      if (rp_dual_steepest_edge) {
        if (lc_dual_edge_weight_mode ==
            (int)DualEdgeWeightMode::STEEPEST_EDGE) {
          use_row_DSE_density =
              lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE];
        } else {
          use_row_DSE_density = 0;
        }
        reportOneDensity(use_row_DSE_density);
      }
      printf(" |  %3s ", str_dual_edge_weight_mode.c_str());
      if (rp_dual_steepest_edge) {
        double use_costly_dse;
        printf("|     ");
        if (lc_dual_edge_weight_mode ==
            (int)DualEdgeWeightMode::STEEPEST_EDGE) {
          use_costly_dse = lcAnIter.AnIterTraceCostlyDse;
        } else {
          use_costly_dse = 0;
        }
        reportOneDensity(use_costly_dse);
      }
      printf("\n");
      fmIter = toIter;
      fmTime = toTime;
    }
    printf("\n");
    // Remove any temporary record added for the final iterations
    if (add_extra_record) AnIterTraceNumRec--;
  }
}
#endif

void HighsSimplexAnalysis::iterationReport(const bool header) {
  if (!(iteration_report_message_level & message_level)) return;
  if (!header && (pivotal_row_index < 0 || entering_variable < 0)) return;
  reportAlgorithmPhaseIterationObjective(header,
                                         iteration_report_message_level);
#ifdef HiGHSDEV
  reportDensity(header, iteration_report_message_level);
  reportIterationData(header, iteration_report_message_level);
  //  reportFreeListSize(header, iteration_report_message_level);
#endif
  HighsPrintMessage(output, message_level, iteration_report_message_level,
                    "\n");
  if (!header) num_iteration_report_since_last_header++;
}

void HighsSimplexAnalysis::reportAlgorithmPhaseIterationObjective(
    const bool header, const int this_message_level) {
  if (header) {
    HighsPrintMessage(output, message_level, this_message_level,
                      "       Iteration        Objective    ");
  } else {
    std::string algorithm;
    if (dualAlgorithm()) {
      algorithm = "Du";
    } else {
      algorithm = "Pr";
    }
    HighsPrintMessage(output, message_level, this_message_level,
                      "%2sPh%1d %10d %20.10e", algorithm.c_str(), solve_phase,
                      simplex_iteration_count, objective_value);
  }
}

void HighsSimplexAnalysis::reportInfeasibility(const bool header,
                                               const int this_message_level) {
  if (header) {
    HighsPrintMessage(output, message_level, this_message_level,
                      " Infeasibilities num(sum)");
  } else {
    if (solve_phase == 1) {
      HighsPrintMessage(output, message_level, this_message_level,
                        " Ph1: %d(%g)", num_primal_infeasibilities,
                        sum_primal_infeasibilities);
    } else {
      HighsPrintMessage(output, message_level, this_message_level,
                        " Pr: %d(%g)", num_primal_infeasibilities,
                        sum_primal_infeasibilities);
    }
    if (sum_dual_infeasibilities > 0) {
      HighsPrintMessage(output, message_level, this_message_level,
                        "; Du: %d(%g)", num_dual_infeasibilities,
                        sum_dual_infeasibilities);
    }
  }
}

#ifdef HiGHSDEV
void HighsSimplexAnalysis::reportThreads(const bool header,
                                         const int this_message_level) {
  if (header) {
    HighsPrintMessage(output, message_level, this_message_level, "  Threads");
  } else if (num_threads > 0) {
    HighsPrintMessage(output, message_level, this_message_level, " %2d|%2d|%2d",
                      min_threads, num_threads, max_threads);
  } else {
    HighsPrintMessage(output, message_level, this_message_level, "   |  |  ");
  }
}

void HighsSimplexAnalysis::reportMulti(const bool header,
                                       const int this_message_level) {
  if (header) {
    HighsPrintMessage(output, message_level, this_message_level, "  Multi");
  } else if (average_fraction_of_possible_minor_iterations_performed >= 0) {
    HighsPrintMessage(
        output, message_level, this_message_level, "   %3d%%",
        (int)(100 * average_fraction_of_possible_minor_iterations_performed));
  } else {
    HighsPrintMessage(output, message_level, this_message_level, "       ");
  }
}

void HighsSimplexAnalysis::reportOneDensity(const int this_message_level,
                                            const double density) {
  const int log_10_density = intLog10(density);
  if (log_10_density > -99) {
    HighsPrintMessage(output, message_level, this_message_level, " %4d",
                      log_10_density);
  } else {
    HighsPrintMessage(output, message_level, this_message_level, "     ");
  }
}

void HighsSimplexAnalysis::reportOneDensity(const double density) {
  const int log_10_density = intLog10(density);
  if (log_10_density > -99) {
    printf(" %4d", log_10_density);
  } else {
    printf("     ");
  }
}

void HighsSimplexAnalysis::reportDensity(const bool header,
                                         const int this_message_level) {
  const bool rp_dual_steepest_edge =
      edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE;
  if (header) {
    HighsPrintMessage(output, message_level, this_message_level,
                      " C_Aq R_Ep R_Ap");
    if (rp_dual_steepest_edge) {
      HighsPrintMessage(output, message_level, this_message_level, "  DSE");
    } else {
      HighsPrintMessage(output, message_level, this_message_level, "     ");
    }
  } else {
    reportOneDensity(this_message_level, col_aq_density);
    reportOneDensity(this_message_level, row_ep_density);
    reportOneDensity(this_message_level, row_ap_density);
    double use_row_DSE_density;
    if (rp_dual_steepest_edge) {
      use_row_DSE_density = row_DSE_density;
    } else {
      use_row_DSE_density = 0;
    }
    reportOneDensity(this_message_level, use_row_DSE_density);
  }
}

void HighsSimplexAnalysis::reportInvert(const bool header,
                                        const int this_message_level) {
  if (header) {
    HighsPrintMessage(output, message_level, this_message_level, " Inv");
  } else {
    HighsPrintMessage(output, message_level, this_message_level, "  %2d",
                      invert_hint);
  }
}

#ifdef HiGHSDEV
void HighsSimplexAnalysis::reportCondition(const bool header,
                                           const int this_message_level) {
  if (header) {
    HighsPrintMessage(output, message_level, this_message_level, "       k(B)");
  } else {
    HighsPrintMessage(output, message_level, this_message_level, " %10.4g",
                      basis_condition);
  }
}
#endif

void HighsSimplexAnalysis::reportIterationData(const bool header,
                                               const int this_message_level) {
  if (header) {
    HighsPrintMessage(output, message_level, this_message_level,
                      "       NumCk     LvR     LvC     EnC        DlPr    "
                      "    ThDu        ThPr          Aa");
  } else {
    HighsPrintMessage(output, message_level, this_message_level,
                      " %11.4g %7d %7d %7d %11.4g %11.4g %11.4g %11.4g",
                      numerical_trouble, pivotal_row_index, leaving_variable,
                      entering_variable, primal_delta, dual_step, primal_step,
                      pivot_value_from_column);
  }
}

void HighsSimplexAnalysis::reportFreeListSize(const bool header,
                                              const int this_message_level) {
  if (header) {
    HighsPrintMessage(output, message_level, this_message_level, " FreeLsZ");
  } else {
    HighsPrintMessage(output, message_level, this_message_level, " %7d",
                      freelist_size);
  }
}

int HighsSimplexAnalysis::intLog10(const double v) {
  int intLog10V = -99;
  if (v > 0) intLog10V = log(v) / log(10.0);
  return intLog10V;
}

#endif

bool HighsSimplexAnalysis::dualAlgorithm() {
  return (simplex_strategy == SIMPLEX_STRATEGY_DUAL ||
          simplex_strategy == SIMPLEX_STRATEGY_DUAL_TASKS ||
          simplex_strategy == SIMPLEX_STRATEGY_DUAL_MULTI);
}
