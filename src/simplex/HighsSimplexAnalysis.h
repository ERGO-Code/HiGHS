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
  void dualSteepestEdgeWeightError(const double computed_edge_weight, const double updated_edge_weight);
  bool switchToDevex();

  /*
  void iterationAnalysis();
#ifdef HiGHSDEV
  void iterateOpRecBf(const int opTy, const HVector& vector, const double hist_dsty);
  void iterateOpRecAf(const int opTy, const HVector& vector);
  void iterationAnalysisReport();
#endif
  */
  int intLog10(const double v);

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

 private:

  void iterationReportFull(const bool header);
  void iterationReportIterationAndPhase(const int iterate_log_level, const bool header);
  void iterationReportDualObjective(const int iterate_log_level, const bool header);
  void iterationReportIterationData(const int iterate_log_level, const bool header);
  void iterationReportDensity(const int iterate_log_level, const bool header);

  int AnIterNumCostlyDseIt;  //!< Number of iterations when DSE is costly
  double AnIterCostlyDseFq;  //!< Frequency of iterations when DSE is costly
  const double AnIterCostlyDseMeasureLimit = 1000.0;  //!<
  const double AnIterCostlyDseMnDensity = 0.01;       //!<
  const double AnIterFracNumTot_ItBfSw = 0.1;         //!<
  const double AnIterFracNumCostlyDseItbfSw = 0.05;   //!<
  double AnIterCostlyDseMeasure;
#ifdef HiGHSDEV
  int AnIterPrevRpNumCostlyDseIt;  //!< Number of costly DSE iterations when
                                   //!< previously reported
#endif

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
  int AnIterNumColPrice;
  int AnIterNumRowPrice;
  int AnIterNumRowPriceWSw;
  int AnIterNumRowPriceUltra;
  int AnIterNumPrDgnIt;
  int AnIterNumDuDgnIt;
  int AnIterNumEdWtIt[(int)DualEdgeWeightMode::Count];
#endif
};

#endif /* SIMPLEX_HIGHSSIMPLEXANALYSIS_H_ */
