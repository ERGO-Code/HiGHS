/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexAnalysis.h
 * @brief Analyse simplex iterations, both for run-time control and data gathering
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HIGHSSIMPLEXANALYSIS_H_
#define SIMPLEX_HIGHSSIMPLEXANALYSIS_H_

#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "simplex/SimplexConst.h"
#include "simplex/HVector.h"
#include "util/HighsTimer.h"

#ifdef HiGHSDEV
  enum ANALYSIS_OPERATION_TYPE {
    ANALYSIS_OPERATION_TYPE_BTRAN = 0,
    ANALYSIS_OPERATION_TYPE_PRICE,
    ANALYSIS_OPERATION_TYPE_FTRAN,
    ANALYSIS_OPERATION_TYPE_FTRAN_BFRT,
    ANALYSIS_OPERATION_TYPE_FTRAN_DSE,
    NUM_ANALYSIS_OPERATION_TYPE,
  };
#endif
  const double running_average_multiplier = 0.05;

/**
 * @brief Analyse simplex iterations, both for run-time control and data gathering
 */
class HighsSimplexAnalysis {
 public:
  void setup(const HighsLp& lp, const HighsOptions& options);

  void messaging(FILE* logfile_, FILE* output_, const int message_level_);

  void initialise(const int simplex_iteration_count);

  void updateOperationResultDensity(const double local_density, double& density);

  //  void equalDensity(const double density0, const double density1);

  void iterationReport();
  void invertReport();
  void dualSteepestEdgeWeightError(const double computed_edge_weight, const double updated_edge_weight);
  bool switchToDevex();

#ifdef HiGHSDEV
  void iterationRecord();
  void operationRecordBefore(const int operation_type, const HVector& vector, const double historical_density);
  void operationRecordAfter(const int operation_type, const HVector& vector);
  void summaryReport();
#endif

  HighsTimer timer_;

  int numRow;
  int numCol;
  int numTot;
  bool allow_dual_steepest_edge_to_devex_switch;
  double dual_steepest_edge_weight_log_error_threshhold;
  FILE* logfile;
  FILE* output;
  int message_level;

  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_DSE_density;
  
  int simplex_strategy;
  int num_threads;
  DualEdgeWeightMode edge_weight_mode;
  int solve_phase;
  int simplex_iteration_count;
  int major_iteration_count;
  int minor_iteration_count;
  int devex_iteration_count;
  int pivotal_row_index;
  int leaving_variable;
  int entering_variable;
  int num_primal_infeasibilities;
  int num_dual_infeasibilities;
  int invert_hint;
  int freelist_size;
  double reduced_rhs_value;
  double reduced_cost_value;
  double edge_weight;
  double primal_delta;
  double primal_step;
  double dual_step;
  double pivot_value_from_column;
  double pivot_value_from_row;
  double numerical_trouble;
  double objective_value;
  double sum_primal_infeasibilities;
  double sum_dual_infeasibilities;
  double basis_condition;

  int num_col_price;
  int num_row_price;
  int num_row_price_with_switch;

 private:

  void iterationReport(const bool header);
  void invertReport(const bool header);
  void reportAlgorithmPhaseIterationObjective(const bool header, const int this_message_level);
  void reportInfeasibility(const bool header, const int this_message_level);
#ifdef HiGHSDEV
  void reportDensity(const bool header, const int this_message_level);
  void reportInvert(const bool header, const int this_message_level);
  void reportCondition(const bool header, const int this_message_level);
  void reportIterationData(const bool header, const int this_message_level);
  void reportFreeListSize(const bool header, const int this_message_level);
  int intLog10(const double v);
#endif
  bool dualAlgorithm();

  int AnIterNumCostlyDseIt;  //!< Number of iterations when DSE is costly
  double AnIterCostlyDseFq;  //!< Frequency of iterations when DSE is costly
  const double AnIterCostlyDseMeasureLimit = 1000.0;  //!<
  const double AnIterCostlyDseMnDensity = 0.01;       //!<
  const double AnIterFracNumTot_ItBfSw = 0.1;         //!<
  const double AnIterFracNumCostlyDseItbfSw = 0.05;   //!<
  double AnIterCostlyDseMeasure;
  int AnIterPrevRpNumCostlyDseIt;  //!< Number of costly DSE iterations when
                                   //!< previously reported

  const double accept_weight_threshhold = 0.25;
  const double weight_error_threshhold = 4.0;

  int num_dual_steepest_edge_weight_check;
  int num_dual_steepest_edge_weight_reject;
  int num_wrong_low_dual_steepest_edge_weight;
  int num_wrong_high_dual_steepest_edge_weight;
  double average_frequency_low_dual_steepest_edge_weight;
  double average_frequency_high_dual_steepest_edge_weight;
  double average_log_low_dual_steepest_edge_weight_error;
  double average_log_high_dual_steepest_edge_weight_error;
  double max_average_frequency_low_dual_steepest_edge_weight;
  double max_average_frequency_high_dual_steepest_edge_weight;
  double max_sum_average_frequency_extreme_dual_steepest_edge_weight;
  double max_average_log_low_dual_steepest_edge_weight_error;
  double max_average_log_high_dual_steepest_edge_weight_error;
  double max_sum_average_log_extreme_dual_steepest_edge_weight_error;

  int previous_iteration_report_header_iteration_count = -1;

  int AnIterIt0;
#ifdef HiGHSDEV
  int AnIterPrevIt;
  // Major operation analysis struct
  struct AnIterOpRec {
    double AnIterOpLog10RsDsty;
    double AnIterOpSuLog10RsDsty;
    double AnIterOpHyperCANCEL;
    double AnIterOpHyperTRAN;
    int AnIterOpRsDim;
    int AnIterOpNumCa;
    int AnIterOpNumHyperOp;
    int AnIterOpNumHyperRs;
    int AnIterOpRsMxNNZ;
    int AnIterOpSuNumCa;
    int AnIterOpSuNumHyperOp;
    int AnIterOpSuNumHyperRs;
    std::string AnIterOpName;
  };
  AnIterOpRec AnIterOp[NUM_ANALYSIS_OPERATION_TYPE];

  struct AnIterTraceRec {
    double AnIterTraceTime;
    double AnIterTraceDsty[NUM_ANALYSIS_OPERATION_TYPE];
    double AnIterTraceAux0;
    int AnIterTraceIter;
    int AnIterTrace_dual_edge_weight_mode;
  };

  enum AnIterTraceMxNumRec { AN_ITER_TRACE_MX_NUM_REC = 20 };
  enum DUAL_EDGE_WEIGHT_MODE_COUNT { DUAL_EDGE_WEIGHT_MODE_COUNT = 3 };
  int AnIterTraceNumRec;
  int AnIterTraceIterDl;
  AnIterTraceRec AnIterTrace[1 + AN_ITER_TRACE_MX_NUM_REC + 1];

  int AnIterNumInvert[INVERT_HINT_Count];
  int AnIterNumPrDgnIt;
  int AnIterNumDuDgnIt;
  int AnIterNumEdWtIt[(int)DualEdgeWeightMode::Count];
#endif
};

#endif /* SIMPLEX_HIGHSSIMPLEXANALYSIS_H_ */
