/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexAnalysis.h
 * @brief Analyse simplex iterations, both for run-time control and data
 * gathering
 */
#ifndef SIMPLEX_HIGHSSIMPLEXANALYSIS_H_
#define SIMPLEX_HIGHSSIMPLEXANALYSIS_H_

#include <sstream>

#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "simplex/HVector.h"
#include "simplex/SimplexConst.h"
#include "util/HighsInt.h"
#include "util/HighsTimer.h"
#include "util/HighsUtils.h"

//#ifdef OPENMP
//#include "omp.h"
//#endif

enum ANALYSIS_OPERATION_TYPE {
  ANALYSIS_OPERATION_TYPE_BTRAN_FULL = 0,
  ANALYSIS_OPERATION_TYPE_PRICE_FULL,
  ANALYSIS_OPERATION_TYPE_BTRAN_BASIC_FEASIBILITY_CHANGE,
  ANALYSIS_OPERATION_TYPE_PRICE_BASIC_FEASIBILITY_CHANGE,
  ANALYSIS_OPERATION_TYPE_BTRAN_EP,
  ANALYSIS_OPERATION_TYPE_PRICE_AP,
  ANALYSIS_OPERATION_TYPE_FTRAN,
  ANALYSIS_OPERATION_TYPE_FTRAN_BFRT,
  ANALYSIS_OPERATION_TYPE_FTRAN_DSE,
  NUM_ANALYSIS_OPERATION_TYPE,
};

enum TRAN_STAGE {
  TRAN_STAGE_FTRAN_LOWER = 0,
  TRAN_STAGE_FTRAN_UPPER_FT,
  TRAN_STAGE_FTRAN_UPPER,
  TRAN_STAGE_BTRAN_UPPER,
  TRAN_STAGE_BTRAN_UPPER_FT,
  TRAN_STAGE_BTRAN_LOWER,
  NUM_TRAN_STAGE_TYPE,
};

struct TranStageAnalysis {
  std::string name_;
  HighsScatterData rhs_density_;
  HighsInt num_decision_;
  HighsInt num_wrong_original_sparse_decision_;
  HighsInt num_wrong_original_hyper_decision_;
  HighsInt num_wrong_new_sparse_decision_;
  HighsInt num_wrong_new_hyper_decision_;
};

// const double running_average_multiplier = 0.05;
const double max_regression_density = 0.2;
const double max_hyper_density = 0.1;

/**
 * @brief Analyse simplex iterations, both for run-time control and data
 * gathering
 */
class HighsSimplexAnalysis {
 public:
  HighsSimplexAnalysis(HighsTimer& timer) : timer_reference(timer) {
    timer_ = &timer;
    /*
    HighsInt omp_max_threads = 1;
#ifdef OPENMP
    omp_max_threads = omp_get_max_threads();
#endif
    for (HighsInt i = 0; i < omp_max_threads; i++) {
      HighsTimerClock clock(timer);
      thread_simplex_clocks.push_back(clock);
      thread_factor_clocks.push_back(clock);
    }
    pointer_serial_factor_clocks = &thread_factor_clocks[0];
    */
    //#ifdef HiGHSDEV
    //#else
    //    pointer_serial_factor_clocks = NULL;
    //#endif
  }
  // Reference and pointer to timer
  HighsTimer& timer_reference;
  HighsTimer* timer_;

  void setup(const HighsLp& lp, const HighsOptions& options,
             const HighsInt simplex_iteration_count);
  void messaging(const HighsLogOptions& log_options_);
  void iterationReport();
  void invertReport();
  void invertReport(const bool header);
  bool predictEndDensity(const HighsInt tran_stage_id,
                         const double start_density, double& end_density);
  void afterTranStage(const HighsInt tran_stage_id, const double start_density,
                      const double end_density, const double historical_density,
                      const double predicted_end_density,
                      const bool use_solve_sparse_original_HFactor_logic,
                      const bool use_solve_sparse_new_HFactor_logic);

  void simplexTimerStart(const HighsInt simplex_clock,
                         const HighsInt thread_id = 0);
  void simplexTimerStop(const HighsInt simplex_clock,
                        const HighsInt thread_id = 0);
  bool simplexTimerRunning(const HighsInt simplex_clock,
                           const HighsInt thread_id = 0);
  HighsInt simplexTimerNumCall(const HighsInt simplex_clock,
                               const HighsInt thread_id = 0);
  double simplexTimerRead(const HighsInt simplex_clock,
                          const HighsInt thread_id = 0);

  HighsTimerClock* getThreadFactorTimerClockPointer();

  const std::vector<HighsTimerClock>& getThreadSimplexTimerClocks() {
    return thread_simplex_clocks;
  }
  HighsTimerClock* getThreadSimplexTimerClockPtr(HighsInt i) {
    assert(i >= 0 && i < (HighsInt)thread_simplex_clocks.size());
    return &thread_simplex_clocks[i];
  }

  const std::vector<HighsTimerClock>& getThreadFactorTimerClocks() {
    return thread_factor_clocks;
  }
  HighsTimerClock* getThreadFactorTimerClockPtr(HighsInt i) {
    assert(i >= 0 && i < (HighsInt)thread_factor_clocks.size());
    return &thread_factor_clocks[i];
  }

  void iterationRecord();
  void iterationRecordMajor();
  void operationRecordBefore(const HighsInt operation_type,
                             const HVector& vector,
                             const double historical_density);
  void operationRecordBefore(const HighsInt operation_type,
                             const HighsInt current_count,
                             const double historical_density);
  void operationRecordAfter(const HighsInt operation_type,
                            const HVector& vector);
  void operationRecordAfter(const HighsInt operation_type,
                            const HighsInt result_count);
  void summaryReport();
  void summaryReportFactor();
  void reportSimplexTimer();
  void reportFactorTimer();
  void updateInvertFormData(const HFactor& factor);
  void reportInvertFormData();

  // Control methods to be moved to HEkkControl
  void updateOperationResultDensity(const double local_density,
                                    double& density);
  void dualSteepestEdgeWeightError(const double computed_edge_weight,
                                   const double updated_edge_weight);
  bool switchToDevex();

  std::vector<HighsTimerClock> thread_simplex_clocks;
  std::vector<HighsTimerClock> thread_factor_clocks;
  HighsTimerClock* pointer_serial_factor_clocks;

  // Local copies of LP data
  HighsInt numRow;
  HighsInt numCol;
  HighsInt numTot;
  std::string model_name_;
  std::string lp_name_;

  // Local copies of IO data
  HighsLogOptions log_options;

  // Interpreted shortcuts from bit settings in highs_analysis_level
  bool analyse_lp_data;
  bool analyse_simplex_data;
  bool analyse_simplex_time;
  bool analyse_factor_data;
  bool analyse_factor_time;

  // Control parameters moving to simplex_info
  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_DSE_density;
  double col_basic_feasibility_change_density;
  double row_basic_feasibility_change_density;
  double col_BFRT_density;
  double primal_col_density;
  double dual_col_density;
  bool allow_dual_steepest_edge_to_devex_switch;
  double dual_steepest_edge_weight_log_error_threshold;

  // Local copies of simplex data for reporting
  HighsInt simplex_strategy = 0;
  DualEdgeWeightMode edge_weight_mode = DualEdgeWeightMode::kSteepestEdge;
  HighsInt solve_phase = 0;
  HighsInt simplex_iteration_count = 0;
  HighsInt devex_iteration_count = 0;
  HighsInt pivotal_row_index = 0;
  HighsInt leaving_variable = 0;
  HighsInt entering_variable = 0;
  HighsInt rebuild_reason = 0;
  double reduced_rhs_value = 0;
  double reduced_cost_value = 0;
  double edge_weight = 0;
  double primal_delta = 0;
  double primal_step = 0;
  double dual_step = 0;
  double pivot_value_from_column = 0;
  double pivot_value_from_row = 0;
  double factor_pivot_threshold = 0;
  double numerical_trouble = 0;
  double objective_value = 0;
  HighsInt num_primal_infeasibility = 0;
  HighsInt num_dual_infeasibility = 0;
  double sum_primal_infeasibility = 0;
  double sum_dual_infeasibility = 0;
  // This triple is an original infeasiblility record, so it includes max,
  // but it's only used for reporting
  HighsInt num_dual_phase_1_lp_dual_infeasibility = 0;
  double max_dual_phase_1_lp_dual_infeasibility = 0;
  double sum_dual_phase_1_lp_dual_infeasibility = 0;
  HighsInt num_devex_framework = 0;

  // Local copies of parallel simplex data for reporting
  HighsInt multi_iteration_count = 0;
  HighsInt multi_chosen = 0;
  HighsInt multi_finished = 0;
  HighsInt min_threads = 0;
  HighsInt num_threads = 0;
  HighsInt max_threads = 0;

  // Unused
  //  HighsInt multi_num = 0; // Useless
  //  double basis_condition = 0; // Maybe useful

  // Records of how pivotal row PRICE was done
  HighsInt num_col_price = 0;
  HighsInt num_row_price = 0;
  HighsInt num_row_price_with_switch = 0;

  HighsValueDistribution before_ftran_upper_sparse_density;
  HighsValueDistribution ftran_upper_sparse_density;
  HighsValueDistribution before_ftran_upper_hyper_density;
  HighsValueDistribution ftran_upper_hyper_density;
  HighsValueDistribution cost_perturbation1_distribution;
  HighsValueDistribution cost_perturbation2_distribution;
  HighsValueDistribution cleanup_dual_change_distribution;
  HighsValueDistribution cleanup_primal_step_distribution;
  HighsValueDistribution cleanup_dual_step_distribution;
  HighsValueDistribution cleanup_primal_change_distribution;

  // Tolerances for analysis of TRAN stages - could be needed for
  // control if this is ever used again!
  vector<double> original_start_density_tolerance;
  vector<double> new_start_density_tolerance;
  vector<double> historical_density_tolerance;
  vector<double> predicted_density_tolerance;
  vector<TranStageAnalysis> tran_stage;

  std::stringstream analysis_log;

 private:
  void iterationReport(const bool header);
  void reportAlgorithmPhaseIterationObjective(const bool header);
  void reportInfeasibility(const bool header);
  void reportThreads(const bool header);
  void reportMulti(const bool header);
  void reportOneDensity(const double density);
  void printOneDensity(const double density);
  void reportDensity(const bool header);
  void reportInvert(const bool header);
  //  void reportCondition(const bool header);
  void reportIterationData(const bool header);
  void reportFreeListSize(const bool header);
  HighsInt intLog10(const double v);
  bool dualAlgorithm();

  HighsInt AnIterNumCostlyDseIt;  //!< Number of iterations when DSE is costly
  double AnIterCostlyDseFq;  //!< Frequency of iterations when DSE is costly
  const double AnIterCostlyDseMeasureLimit = 1000.0;  //!<
  const double AnIterCostlyDseMnDensity = 0.01;       //!<
  const double AnIterFracNumTot_ItBfSw = 0.1;         //!<
  const double AnIterFracNumCostlyDseItbfSw = 0.05;   //!<
  double AnIterCostlyDseMeasure;

  const double accept_weight_threshold = 0.25;
  const double weight_error_threshold = 4.0;

  HighsInt num_dual_steepest_edge_weight_check = 0;
  HighsInt num_dual_steepest_edge_weight_reject = 0;
  HighsInt num_wrong_low_dual_steepest_edge_weight = 0;
  HighsInt num_wrong_high_dual_steepest_edge_weight = 0;
  double average_frequency_low_dual_steepest_edge_weight = 0;
  double average_frequency_high_dual_steepest_edge_weight = 0;
  double average_log_low_dual_steepest_edge_weight_error = 0;
  double average_log_high_dual_steepest_edge_weight_error = 0;
  double max_average_frequency_low_dual_steepest_edge_weight = 0;
  double max_average_frequency_high_dual_steepest_edge_weight = 0;
  double max_sum_average_frequency_extreme_dual_steepest_edge_weight = 0;
  double max_average_log_low_dual_steepest_edge_weight_error = 0;
  double max_average_log_high_dual_steepest_edge_weight_error = 0;
  double max_sum_average_log_extreme_dual_steepest_edge_weight_error = 0;

  const HighsLogType iteration_report_log_type = HighsLogType::kVerbose;
  const HighsLogType invert_report_log_type = HighsLogType::kInfo;
  HighsInt num_invert_report_since_last_header = -1;
  HighsInt num_iteration_report_since_last_header = -1;

  double average_num_threads;
  double average_fraction_of_possible_minor_iterations_performed;
  HighsInt sum_multi_chosen = 0;
  HighsInt sum_multi_finished = 0;

  // Analysis of INVERT form
  HighsInt num_invert = 0;
  HighsInt num_kernel = 0;
  HighsInt num_major_kernel = 0;
  const double major_kernel_relative_dim_threshold = 0.1;
  double max_kernel_dim = 0;
  double sum_kernel_dim = 0;
  double running_average_kernel_dim = 0;
  double sum_invert_fill_factor = 0;
  double sum_kernel_fill_factor = 0;
  double sum_major_kernel_fill_factor = 0;
  double running_average_invert_fill_factor = 1;
  double running_average_kernel_fill_factor = 1;
  double running_average_major_kernel_fill_factor = 1;

  HighsInt AnIterIt0 = 0;
  HighsInt AnIterPrevIt;

  // Major operation analysis struct
  struct AnIterOpRec {
    double AnIterOpHyperCANCEL;
    double AnIterOpHyperTRAN;
    HighsInt AnIterOpRsDim;
    HighsInt AnIterOpNumCa;
    HighsInt AnIterOpNumHyperOp;
    HighsInt AnIterOpNumHyperRs;
    double AnIterOpSumLog10RsDensity;
    HighsInt AnIterOpRsMxNNZ;
    std::string AnIterOpName;
    HighsValueDistribution AnIterOp_density;
  };
  AnIterOpRec AnIterOp[NUM_ANALYSIS_OPERATION_TYPE];

  struct AnIterTraceRec {
    double AnIterTraceTime;
    double AnIterTraceMulti;
    double AnIterTraceDensity[NUM_ANALYSIS_OPERATION_TYPE];
    double AnIterTraceCostlyDse;
    HighsInt AnIterTraceIter;
    HighsInt AnIterTrace_dual_edge_weight_mode;
  };

  // Has to be an enum, not const HighsInt
  enum AnIterTraceMaxNumRec { kAnIterTraceMaxNumRec = 20 };
  HighsInt AnIterTraceNumRec;
  HighsInt AnIterTraceIterDl;
  AnIterTraceRec AnIterTrace[1 + kAnIterTraceMaxNumRec + 1];

  HighsInt AnIterNumInvert[REBUILD_REASON_Count];
  HighsInt AnIterNumEdWtIt[(HighsInt)DualEdgeWeightMode::kCount];

  HighsValueDistribution primal_step_distribution;
  HighsValueDistribution dual_step_distribution;
  HighsValueDistribution simplex_pivot_distribution;
  HighsValueDistribution numerical_trouble_distribution;
  HighsValueDistribution factor_pivot_threshold_distribution;
};

#endif /* SIMPLEX_HIGHSSIMPLEXANALYSIS_H_ */
