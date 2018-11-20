/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
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

#include "HCrash.h"
#include "HDualRHS.h"
#include "HDualRow.h"
#include "HFactor.h"
#include "HMatrix.h"
#include "HModel.h"
#include "HVector.h"

#include <set>
#include <string>
#include <vector>
using namespace std;

/**
 * Limit on number of threads used to dimension many identifiers
 */
const int HSOL_THREAD_LIMIT = 32;
/**
 * Limit on the number of column slices for parallel calculations
 */
const int HSOL_SLICED_LIMIT = 100;

/**
 * Possible edge weight mode values used to test EdWt_Mode
 */
const int EdWt_Mode_DSE = 0;
const int EdWt_Mode_Dvx = 1;
const int EdWt_Mode_Dan = 2;

/**
 * Possible pricing mode values used to test Price_Mode
 */
const int Price_Mode_Row = 0;
const int Price_Mode_Col = 1;

/**
 * Possible presolve mode values used to test Presolve_Mode
 */
const int Presolve_Mode_Off = 0;
const int Presolve_Mode_On = 1;

/**
 * Devex status flags. Each column has a Devex flag which is used as a
 * multiplier to save a conditional branch
 */
const int dvx_in_R = 1;
const int dvx_not_in_R = 0;

/**
 * Parameters controlling number of Devex iterations.
 *
 * There is a new Devex framework if either
 *
 * 1) The weight inaccuracy ratio exceeds maxAllowedDevexWeightRatio
 *
 * 2) There have been max(minAbsNumberDevexIterations,
 * numRow/minRlvNumberDevexIterations) devex iterations
 */
const int minAbsNumberDevexIterations = 25;
const double minRlvNumberDevexIterations = 1e-2;
const double maxAllowedDevexWeightRatio = 3.0;

/**
 * Multiplier used in running average calculations
 */
const double runningAverageMu = 0.05;

enum HDUAL_VARIANT {
  HDUAL_VARIANT_PLAIN = 0,
  HDUAL_VARIANT_TASKS,
  HDUAL_VARIANT_MULTI,
};

/**
 * @brief Dual simplex solver for HiGHS
 */
class HDual {
 public:
  /**
   * @brief Solve a model instance with a dual simplex variant and given number
   * of threads
   */
  void solve(
      HModel *model,       //!< Instance of HModel class to be solved
      int variant = 0,     //!< Default dual simplex variant is "PLAIN" (serial)
      int num_threads = 1  //!< Default number of threads is 1
  );

 public:
  /**
   * @brief Initialise a dual simplex instance
   *
   * Copy dimensions and pointers to matrix, factor and solver-related
   * model data, plus tolerances. Sets up local vectors (columnDSE,
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
  void init_slice(int init_sliced_num  //!< Ideal number of slices - true number
                                       //!< is modified in light of limits
  );

  /**
   * @brief Perform Phase 1 dual simplex iterations
   */
  void solve_phase1();

  /**
   * @brief Perform Phase 2 dual simplex iterations
   */
  void solve_phase2();

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
  void iterate_tasks();

  /**
   * @brief Perform a single PAMI dual simplex iteration - source code in
   * HDualMulti.cpp
   */
  void iterate_multi();  // in HDualMulti.cpp

  /**
   * @brief Initialise the iteration analysis
   */
  void iterateIzAn();

  /**
   * @brief Perform the iteration analysis
   */
  void iterateAn();

  /**
   * @brief Report on the iteration using iterateRpFull, possibly using it to
   * write out column headers
   */
  void iterateRp();

  /**
   * @brief Report full iteration headers or data according to value of
   * <tt>header</tt>
   */
  void iterateRpFull(bool header  //!< Logic to determine whether to write out
                                  //!< column headers or data
  );
  /**
   * @brief Report iteration number and LP phase headers or data according to
   * value of <tt>header</tt>
   */
  void iterateRpIterPh(bool header  //!< Logic to determine whether to write out
                                    //!< column headers or data
  );
  /**
   * @brief Report dual objective value header or data according to value of
   * <tt>header</tt>
   */
  void iterateRpDuObj(bool header  //!< Logic to determine whether to write out
                                   //!< column header or data
  );
  /**
   * @brief Single line report after INVERT
   */
  void iterateRpInvert(
      int i_v  //!< Integer value to be reported - generally invertHint
  );

  /**
   * @brief Update an average density record for BTRAN, an FTRAN or PRICE
   */
  void uOpRsDensityRec(
      double lc_OpRsDensity,  //!< Recent density of the operation
      double &opRsDensity     //!< Average density of the operation
  );
  /**
   * @brief Choose the index of a good row to leave the basis (CHUZR)
   */
  void chooseRow();

  /**
   * @brief Compute pivot row (PRICE) and choose the index of a good column to
   * enter the basis (CHUZC)
   */
  void chooseColumn(HVector *row_ep);

  /**
   * @brief Choose the index of a good column to enter the basis (CHUZC) by
   * exploiting slices of the pivotal row - for SIP and PAMI
   */
  void chooseColumn_slice(HVector *row_ep);

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
   * @brief Compute the vector required to update DSE weights - being
   * FTRAN applied to the pivotal column (FTRAN-DSE)
   */
  void updateFtranDSE(HVector *DSE_Vector  //!< Pivotal column as RHS for FTRAN
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
  void updatePrimal(HVector *DSE_Vector  //!< FTRANned pivotal column
  );

  /**
   * @brief Update the basic and nonbasic variables, iteration count,
   * invertible representation of the basis matrix and row-wise
   * representation of the nonbasic columns, delete the Freelist entry
   * for the entering column, update the primal value for the row where
   * the basis change has occurred, and set the corresponding squared
   * primal infeasibility value in dualRHS.workArray, and then determine
   * whether to reinvert according to the synthetic clock
   */
  void updatePivots();

  /**
   * @brief Initialise a Devex framework: reference set is all basic
   * variables
   */
  void iz_dvx_fwk();

  /**
   * @brief Sets a run-time parameter. TODO: handle this otherwise
   */
  void setCrash(const char *CrashMode);

  /**
   * @brief Set a run-time parameter. TODO: handle this otherwise
   */
  void setPrice(const char *PriceMode);

  /**
   * @brief Set a run-time parameter. TODO: handle this otherwise
   */
  void setEdWt(const char *EdWtMode);

  /**
   * @brief Set a run-time parameter. TODO: handle this otherwise
   */
  void setTimeLimit(double TimeLimit_ArgV);

  /**
   * @brief Set a run-time parameter. TODO: handle this otherwise
   */
  void setPresolve(const char *PresolveMode);

  /**
   * @brief Get a row of the inverse of the basis matrix for SCIP
   */
  int util_getBasisInvRow(int r,         //!< Index of row required
                          double *coef,  //!< Value of entries in row required
                          int *inds,     //!< Indices of entries in row required
                          int *ninds     //!< Number of indices in row required
  );
  /**
   * @brief Get the Hager condition number estimate for the basis matrix of a
   * model
   */
  double an_bs_cond(
      HModel *ptr_model  //!< Model for which basis condition is required
  );

  /**
   * @brief PAMI: Choose the indices of a good set of rows to leave the
   * basis (CHUZR)
   */
  void major_chooseRow();

  /**
   * @brief PAMI: Perform multiple BTRAN
   */
  void major_chooseRowBtran();

  /**
   * @brief PAMI: Choose the index (from the set of indices) of a good
   * row to leave the basis (CHUZR-MI)
   */
  void minor_chooseRow();

  /**
   * @brief PAMI: Update the data during minor iterations
   */
  void minor_update();

  /**
   * @brief PAMI: Update the dual values during minor iterations
   */
  void minor_updateDual();

  /**
   * @brief PAMI: Update the primal values during minor iterations
   */
  void minor_updatePrimal();

  /**
   * @brief PAMI: Perform a basis change during minor iterations
   */
  void minor_updatePivots();

  /**
   * @brief PAMI: Update the tableau rows during minor iterations
   */
  void minor_updateRows();

  /**
   * @brief PAMI: Perform updates after a set of minor iterations
   */
  void major_update();

  /**
   * @brief PAMI: Prepare for the FTRANs after a set of minor iterations
   */
  void major_updateFtranPrepare();

  /**
   * @brief PAMI: Perform the parallel part of multiple FTRANs after a
   * set of minor iterations
   */
  void major_updateFtranParallel();

  /**
   * @brief PAMI: Perform the final part of multiple FTRANs after a set
   * of minor iterations
   */
  void major_updateFtranFinal();

  /**
   * @brief PAMI: Update the primal values after a set of minor
   * iterations
   */
  void major_updatePrimal();

  /**
   * @brief PAMI: Update the invertible representation of the basis
   * matrix after a set of minor iterations
   */
  void major_updateFactor();

  /**
   * @brief PAMI: Roll back some iterations if numerical trouble
   * detected when updating the invertible representation of the basis
   * matrix after a set of minor iterations
   */
  void major_rollback();

#ifdef HiGHSDEV
  void iterateRpIterDa(bool header  //!< Logic to determine whether to write out
                                    //!< column headers or data
  );
  void iterateRpDsty(bool header  //!< Logic to determine whether to write out
                                  //!< column headers or data
  );
  int intLog10(double v);
  void iterateOpRecBf(int opTy, HVector &vector, double hist_dsty);
  void iterateOpRecAf(int opTy, HVector &vector);
  void iterateRpAn();
  void an_iz_vr_v();
#endif

  int dual_variant =
      0;  //!< Dual simplex variant choice. TODO: handle this otherwise
  int Price_Mode = 0;     //!< Pricing mode. TODO: handle this otherwise
  int EdWt_Mode = 0;      //!< Edge weight mode. TODO: handle this otherwise
  int Crash_Mode = 0;     //!< Crash mode. TODO: handle this otherwise
  int Presolve_Mode = 0;  //!< Presolve mode. TODO: handle this otherwise
  bool SolveBailout;  //!< Set true if control is to be returned immediately to
                      //!< calling function
  double TimeLimitValue =
      0;  //!< Value of time limit. TODO: handle this otherwise

#ifdef HiGHSDEV
  // Analysis of rebuilds
  const bool anRebuildTime = false;
  int totalRebuilds;
  double totalRebuildTime;
#endif

  // Devex scalars
  int n_dvx_fwk;    //!< Number of Devex frameworks used
  int n_dvx_it;     //!< Number of Devex iterations with the current framework
  bool nw_dvx_fwk;  //!< Set a new Devex framework
  // Devex vector
  vector<int> dvx_ix;  //!< Vector of Devex indices

  // Price scalars
  bool alw_price_by_col_sw = true;  //!< By default allow switch to column PRICE
                                    //!< if results sufficiently dense
  bool alw_price_by_row_sw =
      true;  //!< By default allow switch to standard row-wise PRICE if result
             //!< is sufficiently dense
  bool alw_price_ultra = false;  //!< By default don't allow ultra-sparse PRICE
  const double dstyColPriceSw = 0.75;  //!< By default switch to column PRICE
                                       //!< when pi_p has at least this density

  // DSE scalars
  bool iz_DSE_wt;  //!< By default initialise DSE weights if initial basis
                   //!< matrix is not an identity
  bool alw_DSE2Dvx_sw = true;  //!< By default allow switch to Devex from DSE
  int AnIterNumCostlyDseIt;    //!< Number of iterations when DSE is costly
  double AnIterCostlyDseFq;    //!< Frequency of iterations when DSE is costly
  const double AnIterCostlyDseMeasureLimit = 1000.0;  //!<
  const double AnIterCostlyDseMnDensity = 0.01;       //!<
  const double AnIterFracNumTot_ItBfSw = 0.1;         //!<
  const double AnIterFracNumCostlyDseItbfSw = 0.05;   //!<
  double AnIterCostlyDseMeasure;
#ifdef HiGHSDEV
  int AnIterPrevRpNumCostlyDseIt;  //!< Number of costly DSE iterations when
                                   //!< previously reported
#endif

#ifdef HiGHSDEV
  int n_wg_DSE_wt;
#endif

  // Model
  HModel *model;
  double Tp;  // Tolerance for primal
  double Td;  // Tolerance for dual

  int numCol;
  int numRow;
  int numTot;
  const HMatrix *matrix;
  const HFactor *factor;

  const int *jMove;
  const double *workRange;
  const double *baseLower;
  const double *baseUpper;
  double *baseValue;
  double *workDual;
  //    JAJH: Only because I can't get these from HModel.h
  double *workValue;
  double *colLower;
  double *colUpper;
  double *rowLower;
  double *rowUpper;
  int *nonbasicFlag;

  vector<double> bs_cond_x;
  vector<double> bs_cond_y;
  vector<double> bs_cond_z;
  vector<double> bs_cond_w;

  int solvePhase;
  int invertHint;

  HVector row_ep;
  HVector row_ap;
  HVector column;
  HVector columnBFRT;
  HVector columnDSE;
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

  // Iteration counts
  int n_ph1_du_it;
  int n_ph2_du_it;
  int n_pr_it;
  // Partitioned coefficient matrix
  int slice_num;
  int slice_PRICE;
  int slice_start[HSOL_SLICED_LIMIT + 1];
  HMatrix slice_matrix[HSOL_SLICED_LIMIT];
  HVector slice_row_ap[HSOL_SLICED_LIMIT];
  HDualRow slice_dualRow[HSOL_SLICED_LIMIT];

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
    HVector column;
    HVector columnBFRT;
  };

  /**
   * @brief Multiple minor iteration data
   */
  struct MFinish {
    int moveIn;
    double shiftOut;
    vector<int> flipList;

    int rowOut;
    int columnOut;
    int columnIn;
    double alphaRow;
    double thetaPrimal;
    double basicBound;
    double basicValue;
    double EdWt;
    HVector_ptr row_ep;
    HVector_ptr column;
    HVector_ptr columnBFRT;
  };

  int multi_num;
  int multi_iChoice;
  int multi_nFinish;
  int multi_iteration;
  int multi_chooseAgain;
  MChoice multi_choice[HSOL_THREAD_LIMIT];
  MFinish multi_finish[HSOL_THREAD_LIMIT];

  double total_syntheticTick;
#ifdef HiGHSDEV
  double total_fake;
#endif
  double total_INVERT_TICK;
  double total_FT_inc_TICK;

  int AnIterIt0;
#ifdef HiGHSDEV
  const bool AnIterLg = true;
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
    string AnIterOpName;
  };
  AnIterOpRec AnIterOp[NumAnIterOpTy];

  struct AnIterTraceRec {
    double AnIterTraceTime;
    double AnIterTraceDsty[NumAnIterOpTy];
    double AnIterTraceAux0;
    int AnIterTraceIter;
    int AnIterTraceEdWt_Mode;
  };

  const int AnIterTraceMxNumRec = 20;
  int AnIterTraceNumRec;
  int AnIterTraceIterDl;
  AnIterTraceRec AnIterTrace[22];  // How can this be 1+AnIterTraceMxNumRec+1;

  const int AnIterNumInvertHint = 7;
  int AnIterNumInvert[8];  // TODO: How can this be AnIterNumInvertHint+1
  int AnIterNumColPrice;
  int AnIterNumRowPrice;
  int AnIterNumRowPriceWSw;
  int AnIterNumRowPriceUltra;
  int AnIterNumPrDgnIt;
  int AnIterNumDuDgnIt;
  int AnIterNumEdWtIt[3];  // TODO: How can this be EdWt_Mode_Dan+1
#endif
};

#endif /* SIMPLEX_HDUAL_H_ */
