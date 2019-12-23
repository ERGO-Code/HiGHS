/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDual.h
 * @brief Dual simplex solver for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HDUAL_H_
#define SIMPLEX_HDUAL_H_

#include <set>
#include <string>
#include <vector>

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HCrash.h"
#include "simplex/HDualRHS.h"
#include "simplex/HDualRow.h"
#include "simplex/HMatrix.h"
#include "simplex/HSimplex.h"
#include "simplex/HVector.h"

class HFactor;

enum class DualEdgeWeightMode { DANTZIG = 0, DEVEX, STEEPEST_EDGE, Count };

enum class PriceMode { ROW = 0, COL };

/**
 * Limit on the number of column slices for parallel calculations. SIP uses num_threads-2 slices; PAMI uses num_threads-1 slices
 */
const int HIGHS_SLICED_LIMIT = HIGHS_THREAD_LIMIT;//Was 100, but can't see why this should be higher than HIGHS_THREAD_LIMIT;

/**
 * Parameters controlling number of Devex iterations.
 *
 * There is a new Devex framework if either
 *
 * 1) The weight inaccuracy ratio exceeds maxAllowedDevexWeightRatio
 *
 * 2) There have been max(minAbsNumberDevexIterations,
 * numRow/minRlvNumberDevexIterations) Devex iterations
 */
const int minAbsNumberDevexIterations = 25;
const double minRlvNumberDevexIterations = 1e-2;
const double maxAllowedDevexWeightRatio = 3.0;

/**
 * Multiplier used in running average calculations
 */
const double runningAverageMu = 0.05;

/**
 * Candidate persistence cut-off in PAMI
 */
const double pami_cutoff = 0.95;

/**
 * @brief Dual simplex solver for HiGHS
 */
class HDual {
 public:
  HDual(HighsModelObject& model_object)
      : workHMO(model_object), dualRow(model_object), dualRHS(model_object) {
    dualRow.setup();
    for (int i = 0; i < HIGHS_SLICED_LIMIT; i++)
      slice_dualRow.push_back(HDualRow(model_object));
    dualRHS.setup();
  }

  /**
   * @brief Solve a model instance with a given number of threads
   */
  HighsStatus solve(int num_threads = 1  //!< Default number of threads is 1
  );

 public:
  /**
   * @brief Set solver options from simplex options
   */
  void options();
  /**
   * @brief Initialise a dual simplex instance
   *
   * Copy dimensions and pointers to matrix, factor and solver-related
   * model data, plus tolerances. Sets up local std::vectors (columnDSE,
   * columnBFRT, column, row_ep and row_ap), scalars for their average
   * density and buffers for dualRow and dualRHS. Also sets up data
   * structures for SIP or PAMI (if necessary).
   */
  void init(int num_threads  //!< Number of threads for initialisation
  );

  /**
   * @brief Initialise matrix slices and slices of row_ap or dualRow for SIP or
   * PAMI
   *
   * TODO generalise call slice_matrix[i].setup_lgBs so slice can be
   * used with non-logical initial basis
   */
  void initSlice(const int init_sliced_num  //!< Ideal number of slices - true number
		                            //!< is modified in light of limits
  );

  /**
   * @brief Perform Phase 1 dual simplex iterations
   */
  void solvePhase1();

  /**
   * @brief Perform Phase 2 dual simplex iterations
   */
  void solvePhase2();

  /**
   * @brief Reinvert if INVERT not fresh, then recompute dual and primal values
   *
   * Also collects primal infeasibilities and computes the dual objective value
   */

  void rebuild();

  /**
   * @brief Remove perturbation and recompute the dual solution
   *
   * Also collects primal infeasibilities and computes the dual objective value
   */
  void cleanup();

  /**
   * @brief Perform a single serial dual simplex iteration
   *
   * All the methods it calls have as their first line "if (invertHint)
   * return;", where invertHint is, for example, set to 1 when CHUZR
   * finds no candidate. This causes a break from the inner loop of
   * solve_phase% and, hence, a call to rebuild().
   */
  void iterate();

  /**
   * @brief Perform a single SIP dual simplex iteration
   */
  void iterateTasks();

  /**
   * @brief Perform a single PAMI dual simplex iteration - source code in
   * HDualMulti.cpp
   */
  void iterateMulti();  // in HDualMulti.cpp

  /**
   * @brief Initialise the iteration analysis
   */
  void iterationAnalysisInitialise();

  /**
   * @brief Perform the iteration analysis
   */
  void iterationAnalysis();

#ifdef HiGHSDEV
  /**
   * @brief Report on the iteration analysis
   */
  void iterationAnalysisReport();
#endif

  /**
   * @brief Report on the iteration using iterationReportFull, possibly using it
   * to write out column headers
   */
  void iterationReport();

  /**
   * @brief Report full iteration headers or data according to value of
   * <tt>header</tt>
   */
  void iterationReportFull(bool header  //!< Logic to determine whether to write
                                        //!< out column headers or data
  );

  /**
   * @brief Report iteration number and LP phase headers or data according to
   * value of <tt>header</tt>
   */
  void iterationReportIterationAndPhase(
      int iterate_log_level,  //!< Iteration logging level
      bool header  //!< Logic to determine whether to write out column headers
                   //!< or data
  );

  /**
   * @brief Report dual objective value header or data according to value of
   * <tt>header</tt>
   */
  void iterationReportDualObjective(
      int iterate_log_level,  //!< Iteration logging level
      bool header  //!< Logic to determine whether to write out column header or
                   //!< data
  );

  /**
   * @brief Report dual iteration data header or data according to value of
   * <tt>header</tt>
   */
  void iterationReportIterationData(
      int iterate_log_level,  //!< Iteration logging level
      bool header  //!< Logic to determine whether to write out column headers
                   //!< or data
  );

  /**
   * @brief Report dual iteration operation density header or data according to
   * value of <tt>header</tt>
   */
  void iterationReportDensity(
      int iterate_log_level,  //!< Iteration logging level
      bool header  //!< Logic to determine whether to write out column headers
                   //!< or data
  );
  int intLog10(double v);

  /**
   * @brief Single line report after rebuild
   */
  void iterationReportRebuild(
#ifdef HiGHSDEV
      const int i_v=-1  //!< Integer value for reporting - generally invertHint
#endif
  );

  /**
   * @brief Report infeasibility
   */
  void reportInfeasibility();

  /**
   * @brief Update an average density record for BTRAN, an FTRAN or PRICE
   */
  void uOpRsDensityRec(
      double lc_OpRsDensity,  //!< Recent density of the operation
      double& opRsDensity     //!< Average density of the operation
  );
  /**
   * @brief Choose the index of a good row to leave the basis (CHUZR)
   */
  void chooseRow();

  /**
   * @brief Determine whether the updated_edge_weight is accurate enough to
   * be accepted, and update the analysis of weight errors
   */
  bool acceptDualSteepestEdgeWeight(const double updated_edge_weight);

  /**
   * @brief Determine whether the updated_edge_weight error should trigger a new Devex framework
   */
  bool newDevexFramework(const double updated_edge_weight);

  /**
   * @brief Compute pivot row (PRICE) and choose the index of a good column to
   * enter the basis (CHUZC)
   */
  void chooseColumn(HVector* row_ep);

  /**
   * @brief Choose the index of a good column to enter the basis (CHUZC) by
   * exploiting slices of the pivotal row - for SIP and PAMI
   */
  void chooseColumnSlice(HVector* row_ep);

  /**
   * @brief Compute the pivotal column (FTRAN)
   */
  void updateFtran();

  /**
   * @brief Compute the RHS changes corresponding to the BFRT
   * (FTRAN-BFRT)
   */
  void updateFtranBFRT();

  /**
   * @brief Compute the std::vector required to update DSE weights - being
   * FTRAN applied to the pivotal column (FTRAN-DSE)
   */
  void updateFtranDSE(HVector* DSE_Vector  //!< Pivotal column as RHS for FTRAN
  );
  /**
   * @brief Compare the pivot value computed row-wise and column-wise
   * and determine whether reinversion is advisable
   */
  void updateVerify();

  /**
   * @brief Update the dual values
   */
  void updateDual();

  /**
   * @brief Update the primal values and any edge weights
   */
  void updatePrimal(HVector* DSE_Vector  //!< FTRANned pivotal column
  );

  /**
   * @brief Update the basic and nonbasic variables, iteration count,
   * invertible representation of the basis matrix and row-wise
   * representation of the nonbasic columns, delete the Freelist entry
   * for the entering column, update the primal value for the row
   * where the basis change has occurred, and set the corresponding
   * primal infeasibility value in dualRHS.work_infeasibility, and
   * then determine whether to reinvert according to the synthetic
   * clock
   */
  void updatePivots();

  /**
   * @brief Initialise a Devex framework: reference set is all basic
   * variables
   */
  void initialiseDevexFramework(const bool parallel=false);

  /**
   * @brief Interpret the dual edge weight strategy as setting of a mode and
   * other actions
   */
  void interpretDualEdgeWeightStrategy(
				       const int simplex_dual_edge_weight_strategy
				       );

  /**
   * @brief Interpret the PRICE strategy as setting of a mode and other actions
   */
  void interpretPriceStrategy(
			      const int simplex_price_strategy
			      );

#ifdef HiGHSDEV
  double checkDualObjectiveValue(const char* message, int phase = 2);
#endif

  /**
   * @brief PAMI: Choose the indices of a good set of rows to leave the
   * basis (CHUZR)
   */
  void majorChooseRow();

  /**
   * @brief PAMI: Perform multiple BTRAN
   */
  void majorChooseRowBtran();

  /**
   * @brief PAMI: Choose the index (from the set of indices) of a good
   * row to leave the basis (CHUZR-MI)
   */
  void minorChooseRow();

  /**
   * @brief PAMI: Update the data during minor iterations
   */
  void minorUpdate();

  /**
   * @brief PAMI: Update the dual values during minor iterations
   */
  void minorUpdateDual();

  /**
   * @brief PAMI: Update the primal values during minor iterations
   */
  void minorUpdatePrimal();

  /**
   * @brief PAMI: Perform a basis change during minor iterations
   */
  void minorUpdatePivots();

  /**
   * @brief PAMI: Update the tableau rows during minor iterations
   */
  void minorUpdateRows();

  /**
   * @brief PAMI: Initialise a new Devex framework during minor iterations
   */
  void minorInitialiseDevexFramework();

  /**
   * @brief PAMI: Perform updates after a set of minor iterations
   */
  void majorUpdate();

  /**
   * @brief PAMI: Prepare for the FTRANs after a set of minor iterations
   */
  void majorUpdateFtranPrepare();

  /**
   * @brief PAMI: Perform the parallel part of multiple FTRANs after a
   * set of minor iterations
   */
  void majorUpdateFtranParallel();

  /**
   * @brief PAMI: Perform the final part of multiple FTRANs after a set
   * of minor iterations
   */
  void majorUpdateFtranFinal();

  /**
   * @brief PAMI: Update the primal values after a set of minor
   * iterations
   */
  void majorUpdatePrimal();

  /**
   * @brief PAMI: Update the invertible representation of the basis
   * matrix after a set of minor iterations
   */
  void majorUpdateFactor();

  /**
   * @brief PAMI: Roll back some iterations if numerical trouble
   * detected when updating the invertible representation of the basis
   * matrix after a set of minor iterations
   */
  void majorRollback();

  bool checkNonUnitWeightError(std::string message);
  bool dualInfoOk(const HighsLp& lp);

#ifdef HiGHSDEV
  void iterateOpRecBf(int opTy, HVector& vector, double hist_dsty);
  void iterateOpRecAf(int opTy, HVector& vector);
#endif

  int Crash_Mode = 0;  //!< Crash mode. TODO: handle this otherwise
  bool solve_bailout;  //!< Set true if control is to be returned immediately to
                       //!< calling function

  // Devex scalars
  int num_devex_framework = 0;   //!< Number of Devex frameworks used
  int num_devex_iterations = 0;  //!< Number of Devex iterations with the current framework
  bool new_devex_framework = false;  //!< Set a new Devex framework
  bool minor_new_devex_framework = false; //!< Set a new Devex framework in PAMI minor iterations
  // Devex std::vector
  std::vector<int> devex_index;  //!< Vector of Devex indices

  // Price scalars
  // DSE scalars
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

  // Model
  HighsModelObject& workHMO;
  int solver_num_row;
  int solver_num_col;
  int solver_num_tot;

  const HMatrix* matrix;
  //  const HFactor* factor; //FactorTimer frig const
  HFactor* factor;

  const int* jMove;
  const double* workRange;
  const double* baseLower;
  const double* baseUpper;
  double* baseValue;
  double* workDual;
  double* workValue;
  double* colLower;
  double* colUpper;
  double* rowLower;
  double* rowUpper;
  int* nonbasicFlag;

  // Options
  DualEdgeWeightMode dual_edge_weight_mode;
  bool initialise_dual_steepest_edge_weights;
  bool allow_dual_steepest_edge_to_devex_switch;

  PriceMode price_mode;
  bool allow_price_by_col_switch;
  bool allow_price_by_row_switch;
  bool allow_price_ultra;
  const double dstyColPriceSw = 0.75;  //!< By default switch to column PRICE
                                       //!< when pi_p has at least this density
  const double min_dual_steepest_edge_weight = 1e-4;

  double Tp;  // Tolerance for primal
  double primal_feasibility_tolerance;

  double Td;  // Tolerance for dual
  double dual_feasibility_tolerance;
  double dual_objective_value_upper_bound;

  int solvePhase;
  int previous_iteration_report_header_iteration_count = -1;
  int invertHint;

  HVector row_ep;
  HVector row_ap;
  HVector col_aq;
  HVector col_BFRT;
  HVector col_DSE;
  double columnDensity;
  double row_epDensity;
  double row_apDensity;
  double rowdseDensity;

  HDualRow dualRow;

  // Solving related buffers
  int dualInfeasCount;

  HDualRHS dualRHS;

  // Simplex pivotal information
  int rowOut;
  int columnOut;
  int sourceOut;  // -1 from small to lower, +1 to upper
  int columnIn;
  double deltaPrimal;
  double thetaDual;
  double thetaPrimal;
  double alpha;
  double alphaRow;
  double numericalTrouble;
  // (Local) value of computed weight
  double computed_edge_weight;

  // Partitioned coefficient matrix
  int slice_num;
  int slice_PRICE;
  int slice_start[HIGHS_SLICED_LIMIT + 1];
  HMatrix slice_matrix[HIGHS_SLICED_LIMIT];
  HVector slice_row_ap[HIGHS_SLICED_LIMIT];
  std::vector<HDualRow> slice_dualRow;

  /**
   * @brief Multiple CHUZR data
   */
  struct MChoice {
    int rowOut;
    double baseValue;
    double baseLower;
    double baseUpper;
    double infeasValue;
    double infeasEdWt;
    double infeasLimit;
    HVector row_ep;
    HVector col_aq;
    HVector col_BFRT;
  };

  /**
   * @brief Multiple minor iteration data
   */
  struct MFinish {
    int moveIn;
    double shiftOut;
    std::vector<int> flipList;

    int rowOut;
    int columnOut;
    int columnIn;
    double alphaRow;
    double thetaPrimal;
    double basicBound;
    double basicValue;
    double EdWt;
    HVector_ptr row_ep;
    HVector_ptr col_aq;
    HVector_ptr col_BFRT;
  };

  int multi_num;
  int multi_iChoice;
  int multi_nFinish;
  int multi_iteration;
  int multi_chooseAgain;
  MChoice multi_choice[HIGHS_THREAD_LIMIT];
  MFinish multi_finish[HIGHS_THREAD_LIMIT];

#ifdef HiGHSDEV
  const bool rp_iter_da = false;//true;//
  const bool rp_reinvert_syntheticClock = false;//true;//
  const bool rp_numericalTrouble = false;//true;//
#endif  
  const double original_multi_build_syntheticTick_mu = 1.5;
  const double multi_build_syntheticTick_mu =
        1.0;
  //original_multi_build_syntheticTick_mu;//
  const double numerical_trouble_tolerance = 1e-7;
  const double original_multi_numerical_trouble_tolerance = 1e-8;
  const double multi_numerical_trouble_tolerance =
    1e-7;
  //original_multi_numerical_trouble_tolerance;
  
  const int synthetic_tick_reinversion_min_update_count = 50;
  const int original_multi_synthetic_tick_reinversion_min_update_count = 201;
  const int multi_synthetic_tick_reinversion_min_update_count =
        synthetic_tick_reinversion_min_update_count;
  //original_multi_synthetic_tick_reinversion_min_update_count;

  double build_syntheticTick;
  double total_syntheticTick;

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

  int AnIterIt0;
#ifdef HiGHSDEV
  int AnIterPrevIt;
  // Major operation analysis struct
  enum AnIterOpTy {
    AnIterOpTy_Btran = 0,
    AnIterOpTy_Price,
    AnIterOpTy_Ftran,
    AnIterOpTy_FtranBFRT,
    AnIterOpTy_FtranDSE,
    NumAnIterOpTy,
  };

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
  AnIterOpRec AnIterOp[NumAnIterOpTy];

  struct AnIterTraceRec {
    double AnIterTraceTime;
    double AnIterTraceDsty[NumAnIterOpTy];
    double AnIterTraceAux0;
    int AnIterTraceIter;
    int AnIterTrace_dual_edge_weight_mode;
  };

  enum AnIterTraceMxNumRec { AN_ITER_TRACE_MX_NUM_REC = 20 };
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

#endif /* SIMPLEX_HDUAL_H_ */
