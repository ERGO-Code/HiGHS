/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HModel.h
 * @brief LP model representation and management for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HMODEL_H_
#define SIMPLEX_HMODEL_H_

#include "HFactor.h"
#include "HMatrix.h"
#include "HighsLp.h"
#include "HighsTimer.h" //For timer_
//#include "HighsUtils.h"
#include "HighsRandom.h"

// For compute dual objective alt value
//#include "HighsModelObject.h"
//#include "HSimplex.h"
class HVector;

#include <sstream>
#include <string>
#include <vector>

// After removing HTimer.h add the following
using std::string;

const int LP_Status_Unset = -1;
const int LP_Status_Optimal = 0;
const int LP_Status_Infeasible = 1;
const int LP_Status_Unbounded = 2;
const int LP_Status_Singular = 3;
const int LP_Status_Failed = 4;
const int LP_Status_ObjUB = 5;
const int LP_Status_OutOfTime = 6;

const int invertHint_no = 0;
const int invertHint_updateLimitReached = 1;
const int invertHint_syntheticClockSaysInvert = 2;
const int invertHint_possiblyOptimal = 3;
const int invertHint_possiblyPrimalUnbounded = 4;
const int invertHint_possiblyDualUnbounded = 5;
const int invertHint_possiblySingularBasis = 6;
const int invertHint_primalInfeasibleInPrimalSimplex = 7;
const int invertHint_chooseColumnFail = 8;

/** SCIP-like basis status for columns and rows */
enum HIGHS_BaseStat {
  HIGHS_BASESTAT_LOWER = 0, /**< (slack) variable is at its lower bound
                              [including fixed variables]*/
  HIGHS_BASESTAT_BASIC = 1, /**< (slack) variable is basic */
  HIGHS_BASESTAT_UPPER = 2, /**< (slack) variable is at its upper bound */
  HIGHS_BASESTAT_ZERO = 3   /**< free variable is non-basic and set to zero */
};
typedef enum HIGHS_BaseStat HIGHS_BASESTAT;

/** HiGHS nonbasicFlag status for columns and rows */
enum nonbasicFlagStat {
  NONBASIC_FLAG_TRUE = 1,  // Nonbasic
  NONBASIC_FLAG_FALSE = 0  // Basic
};

/** HiGHS nonbasicMove status for columns and rows */
enum nonbasicMoveStat {
  NONBASIC_MOVE_UP = 1,   // Free to move (only) up
  NONBASIC_MOVE_DN = -1,  // Free to move (only) down
  NONBASIC_MOVE_ZE = 0    // Fixed or free to move up and down
};

// For INT, DBL and STR options, ensure that ***OPT_COUNT is last since
// this is the number of options and used to dimension as
//***Option[***OPT_COUNT]
enum HIGHS_INT_OPTIONS {
  INTOPT_PRINT_FLAG = 0,  // 0/>=1 = none/do-print
  // If 1\in INTOPT_PRINT_FLAG print all "logical" INTOPT_PRINT_FLAG messages
  // If 2\in INTOPT_PRINT_FLAG print timed PROGRESS
  // If 4\in INTOPT_PRINT_FLAG print iteration log line
  INTOPT_TRANSPOSE_FLAG,  // 0/1 = none/do-transpose if possible
  INTOPT_SCALE_FLAG,      // 0/1 = none/do-scale
  INTOPT_TIGHT_FLAG,      // 0/1 = none/do-tight
  INTOPT_PERMUTE_FLAG,    // 0/1 = none/do-permute
  INTOPT_PERTURB_FLAG,    // 0/1 = none/do-perturb
  INTOPT_LPITLIM,         // iteration limit
  INTOPT_COUNT
};

enum HIGHS_DBL_OPTIONS {
  DBLOPT_TIME_LIMIT = 0,
  DBLOPT_PRIMAL_TOL,
  DBLOPT_DUAL_TOL,
  DBLOPT_PERTURB_BASE,
  DBLOPT_PAMI_CUTOFF,
  DBLOPT_OBJ_UB,  // For SCIP
  DBLOPT_COUNT
};

enum HIGHS_STR_OPTIONS {
  STROPT_PARTITION_FILE = 0,  // name of row partition file
  STROPT_COUNT
};

class HModel {
 public:
  HModel();
  // Methods which load whole models, initialise the basis then
  // allocate and populate (where possible) work* arrays and
  // allocate basis* arrays
  int load_fromToy(const char* filename);
  void load_fromArrays(int XnumCol, int XobjSense, const double* XcolCost,
                       const double* XcolLower, const double* XcolUpper,
                       int XnumRow, const double* XrowLower,
                       const double* XrowUpper, int XnumNz, const int* XAstart,
                       const int* XAindex, const double* XAvalue);

  // Methods which initialise the basis then allocate and populate
  // (where possible) work* arrays and allocate basis* arrays
  void initWithLogicalBasis();
  void extendWithLogicalBasis(int firstcol, int lastcol, int firstrow,
                              int lastrow);

  // Methods which replace the basis then populate (where possible)
  // work* arrays and allocate basis* arrays
  void replaceWithLogicalBasis();
  void replaceWithNewBasis(const int* XbasicIndex);

  // Method to clear the current model
  void clearModel();

  // Methods to modify the current model. Only scaleModel is currently in use
  void scaleModel();
  void scaleCosts();
  void setup_transposeLP();
  void setup_tightenBound();
  void setup_shuffleColumn();

  void setup_for_solve();
  bool OKtoSolve(int level, int phase);

  void initScale();
  bool nonbasicFlagBasicIndex_OK(int XnumCol, int XnumRow);
  bool workArrays_OK(int phase);
  bool allNonbasicMoveVsWorkArrays_OK();
  bool oneNonbasicMoveVsWorkArrays_OK(int var);
  void rp_basis();
  int get_nonbasicMove(int var);
  void setup_numBasicLogicals();
  void copy_impliedBoundsToModelBounds();
  void copy_savedBoundsToModelBounds();
  void mlFg_Clear();
  void mlFg_Update(int mlFg_action);
#ifdef HiGHSDEV
  void mlFg_Report();
#endif

  void initFromNonbasic();
  void replaceFromNonbasic();
  void initBasicIndex();

  void allocate_WorkAndBaseArrays();
  void populate_WorkArrays();
  void initCost(int perturb = 0);
  void initPh2ColCost(int firstcol, int lastcol);
  void initPh2RowCost(int firstrow, int lastrow);
  void initBound(int phase = 2);
  void initPh2ColBound(int firstcol, int lastcol);
  void initPh2RowBound(int firstrow, int lastrow);
  void initValue();
  void initValueFromNonbasic(int firstvar, int lastvar);

  // ???? Housekeeping done from here down ????
  // For the solver:
  // Call INVERT and form dual and primal activities
  int computeFactor();
  void computeDual();
  void computeDualInfeasInDual(int* dualInfeasCount);
  void computeDualInfeasInPrimal(int* dualInfeasCount);
  void correctDual(int* freeInfeasCount);
  void computePrimal();
  void computeDualObjectiveValue(int phase = 2);
#ifdef HiGHSDEV
  double checkDualObjectiveValue(const char *message, int phase = 2);
#endif
  double computePrObj();
  double computePh2Objective(vector<double>& colPrAct);
  int handleRankDeficiency();
  int setSourceOutFmBd(const int columnOut);

  // Utilities for shifting costs and flipping bounds
  void shiftCost(int iCol, double amount);
  void shiftBack(int iCol);
  void flipBound(int iCol);

  // The major model updates. Factor calls factor.update; Matrix
  // calls matrix.update; updatePivots does everything---and is
  // called from the likes of HDual::updatePivots
  void updateFactor(HVector* column, HVector* row_ep, int* iRow, int* hint);
  void updateMatrix(int columnIn, int columnOut);
  void updatePivots(int columnIn, int rowOut, int sourceOut);
#ifdef HiGHSDEV
  // Changes the update method, but only used in HTester.cpp
  void changeUpdate(int updateMethod);
#endif
  void setProblemStatus(int status);

  // Checking methods
#ifdef HiGHSDEV
  // Method to check code to load a model from arrays of data
  void check_load_fromArrays();
  void check_load_fromPostsolve();
#endif
  int writeToMPS(const char* filename);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Esoterica!
  // Initialise the random vectors required by HiGHS
  void initRandomVec();

  // Shift the objective
  void shiftObjectiveValue(double shift);

  // Increment numberIteration (here!) and (possibly) store the pivots for
  // debugging NLA
  void recordPivots(int columnIn, int columnOut, double alpha);
#ifdef HiGHSDEV
  // Store and write out the pivots for debugging NLA
  void writePivots(const char* suffix);
#endif
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  // Utilities to get objective, solution and basis: all just copy what's there
  // with no re-evaluation!
  double util_getObjectiveValue();
  
  void util_getPrimalDualValues(vector<double>& XcolValue,
                                vector<double>& XcolDual,
                                vector<double>& XrowValue,
                                vector<double>& XrowDual
				);
  void util_getNonbasicMove( vector<int> &XnonbasicMove);
  void util_getBasicIndexNonbasicFlag(
				      vector<int> &XbasicIndex,
				      vector<int> &XnonbasicFlag
				      );
  // Utilities to scale or unscale bounds and costs
  void util_scaleRowBoundValue(int iRow, double* XrowLowerValue, double* XrowUpperValue);
  void util_scaleColBoundValue(int iCol, double* XcolLowerValue, double* XcolUpperValue);
  void util_scaleColCostValue(int iCol, double* XcolCostValue);
  void util_unscaleRowBoundValue(int iRow, double* XrowLowerValue, double* XrowUpperValue);
  void util_unscaleColBoundValue(int iCol, double* XcolLowerValue, double* XcolUpperValue);
  void util_unscaleColCostValue(int iCol, double* XcolCostValue);
  // Utilities to get/change costs and bounds
  void util_getCosts(HighsLp& lp, int firstcol, int lastcol, double* XcolCost);
  void util_getColBounds(HighsLp& lp, int firstcol, int lastcol, double* XcolLower,
                         double* XcolUpper);
  void util_getRowBounds(HighsLp& lp, int firstrow, int lastrow, double* XrowLower,
                         double* XrowUpper);
  int util_chgObjSense(int Xobjense);
  int util_chgCostsAll(const double* XcolCost);
  int util_chgCostsSet(int ncols, const int* XcolCostIndex,
                       const double* XcolCostValues);
  int util_chgColBoundsAll(const double* XcolLower, const double* XcolUpper);
  int util_chgColBoundsSet(int ncols, const int* XcolBoundIndex,
                           const double* XcolLowerValues,
                           const double* XcolUpperValues);
  int util_chgRowBoundsAll(const double* XrowLower, const double* XrowUpper);
  int util_chgRowBoundsSet(int nrows, const int* XrowBoundIndex,
                           const double* XrowLowerValues,
                           const double* XrowUpperValues);

  // Utilities to convert model basic/nonbasic status to/from SCIP-like status
  int util_convertBaseStatToWorking(const int* cstat, const int* rstat);
  int util_convertWorkingToBaseStat(int* cstat, int* rstat);
  // Utility to get the indices of the basic variables for SCIP
  int util_getBasicIndices(int* bind);

  // Utilities to add, extract and delete columns and rows
  void util_addCols(int ncols, const double* XcolCost, const double* XcolLower,
                    const double* XcolUpper, int nnonz, const int* XAstart,
                    const int* XAindex, const double* XAvalue);
  void util_deleteCols(int firstcol, int lastcol);
  void util_deleteColset(vector<int>& dstat);
  void util_extractCols(int firstcol, int lastcol, double* XcolLower,
                        double* XcolUpper, int* nnonz, int* XAstart,
                        int* XAindex, double* XAvalue);
  void util_addRows(int nrows, const double* XrowLower, const double* XrowUpper,
                    int nnonz, const int* XARstart, const int* XARindex,
                    const double* XARvalue);
  void util_deleteRows(int firstrow, int lastrow);
  void util_deleteRowset(int* dstat);
  void util_extractRows(int firstrow, int lastrow, double* XrowLower,
                        double* XrowUpper, int* nnonz, int* XARstart,
                        int* XARindex, double* XARvalue);
  void util_changeCoeff(int row, int col, const double newval);
  void util_getCoeff(HighsLp lp, int row, int col, double* val);

  // Methods for brief reports - all just return if intOption[INTOPT_PRINT_FLAG]
  // is false
  void util_reportMessage(const char* message);
  void util_reportNumberIterationObjectiveValue(int i_v);
  void util_reportSolverOutcome(const char* message);
  void util_reportSolverProgress();

  // Methods for reporting the model, its solution, row and column data and
  // matrix
  void util_reportModelDa(HighsLp lp, const char* filename);
  void util_reportModelStatus();
#ifdef HiGHSDEV
  void util_reportModelDense(HighsLp lp);
#endif
  void util_reportRowVecSol(int nrow, vector<double>& XrowLower,
                            vector<double>& XrowUpper,
                            vector<double>& XrowPrimal,
                            vector<double>& XrowDual, vector<int>& XrowStatus);
  void util_reportRowMtx(int nrow, vector<int>& XARstart, vector<int>& XARindex,
                         vector<double>& XARvalue);
  void util_reportColVecSol(int ncol, vector<double>& XcolCost,
                            vector<double>& XcolLower,
                            vector<double>& XcolUpper,
                            vector<double>& XcolPrimal,
                            vector<double>& XcolDual, vector<int>& XcolStatus);

  void util_reportBasicIndex(const char *message, int nrow, vector<int> &basicIndex);
#ifdef HiGHSDEV
  void util_anPrDuDgn();
  void util_anMlLargeCo(HighsLp lp, const char* message);
  void util_anMlSol();
#endif

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Solving options and scalar solution data section: Sort it out!
  // Solving options
  int intOption[INTOPT_COUNT];
  double dblOption[DBLOPT_COUNT];
  string strOption[STROPT_COUNT];

  // Random number generator
  HighsRandom random;

  // Perturbation flag
  int problemPerturbed;

  // Possibly prevent reinversion on optimality in phase 1 or phase 2
  const bool InvertIfRowOutNeg = true;

  const bool forbidSuperBasic = true;

  // Number of basic logicals - allows logical basis to be deduced
  int numBasicLogicals;

  // Booleans to indicate that there are valid implied bounds from
  // presolve and that original bounds have been over-written with
  // them
  bool impliedBoundsPresolve;
  bool usingImpliedBoundsPresolve = false;

  // Solving result
  int limitUpdate;
  int countUpdate;

  // Scalar solution output
  // Essentials
  int numberIteration;
  // Dual objective value
  double dualObjectiveValue;
  double updatedDualObjectiveValue;
#ifdef HiGHSDEV
  double previousUpdatedDualObjectiveValue;
  double previousDualObjectiveValue;
#endif
#ifdef HiGHSDEV
  // Analysis of INVERT
  const bool anInvertTime = false;
  int totalInverts;
  double totalInvertTime;
#endif
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  // Model and solver status flags
  // First the actions---to be passed as parameters to update_mlFg
  const int mlFg_action_TransposeLP = 0;
  const int mlFg_action_ScaleLP = 1;
  const int mlFg_action_ShuffleLP = 2;
  const int mlFg_action_NewCosts = 3;
  const int mlFg_action_NewBounds = 4;
  const int mlFg_action_NewBasis = 5;
  const int mlFg_action_NewCols = 6;
  const int mlFg_action_NewRows = 7;
  const int mlFg_action_DelCols = 8;
  const int mlFg_action_DelRows = 9;
  const int mlFg_action_DelRowsBasisOK = 10;

  int mlFg_transposedLP;
  int mlFg_scaledLP;
  int mlFg_shuffledLP;
  //
  // Basis consists of basicIndex, nonbasicFlag and nonbasicMove. To
  // have them means that they correspond to a consistent basis
  // logically, but B is not necessarily nonsingular.
  int mlFg_haveBasis;
  //
  // Properties of data held in HMatrix.h: MatrixColWise is the copy
  // of the constraint matrix, NOT the model's constraint matrix. To
  // "have" them means that they are correct.
  int mlFg_haveMatrixColWise;
  int mlFg_haveMatrixRowWise;
  //
  // Properties of data held in HFactor.h. To "have" them means that
  // they are assigned.
  int mlFg_haveFactorArrays;
  //
  // This refers to workEdWt, which is held in HDualRHS.h and is
  // assigned and initialised to 1s in dualRHS.setup(model). To
  // "have" the edge weights means that they are correct.
  int mlFg_haveEdWt;
  //
  // The representation of B^{-1} corresponds to the current basis
  int mlFg_haveInvert;
  // The representation of B^{-1} corresponds to the current basis and is fresh
  int mlFg_haveFreshInvert;
  //
  // The nonbasic dual and basic primal values are known
  int mlFg_haveNonbasicDuals;
  int mlFg_haveBasicPrimals;
  //
  // The dual objective function value is known
  int mlFg_haveDualObjectiveValue;
  //
  // The data are fresh from rebuild
  int mlFg_haveFreshRebuild;
  //
  // The ranging information is known
  int mlFg_haveRangingData;
  //
  // Need to know of any saved bounds in the event of scaling being performed
  int mlFg_haveSavedBounds;

 public:
  int problemStatus;
  string modelName;
  
// Limits on scaling factors
  const double minAlwScale = 1 / 1024.0;
  const double maxAlwScale = 1024.0;
  const double maxAlwCostScale = maxAlwScale;
  const double minAlwColScale = minAlwScale;
  const double maxAlwColScale = maxAlwScale;
  const double minAlwRowScale = minAlwScale;
  const double maxAlwRowScale = maxAlwScale;

#ifdef HiGHSDEV
  // Information on large costs
  const double tlLargeCo = 1e5;
  int numLargeCo;
  vector<int> largeCostFlag;
  double largeCostScale;
#endif

  // Associated data of original model
  vector<int> colPermutation;
  vector<double> colRandomValue;

  // The scaled model
  HighsLp *lp_scaled_;
  // Part of working model which is only required and populated once a solve is
  // initiated
  HMatrix *matrix_;
  HFactor *factor_;
  HighsSimplexInfo *simplex_info_;
  HighsBasis *basis_;
  HighsScale *scale_;
  HighsRanging *ranging_;
  HighsTimer *timer_;

#ifdef HiGHSDEV
  vector<int> historyColumnIn;
  vector<int> historyColumnOut;
  vector<double> historyAlpha;
#endif

  // Implied bounds from presolve
  vector<double> primalColLowerImplied;
  vector<double> primalColUpperImplied;
  vector<double> primalRowLowerImplied;
  vector<double> primalRowUpperImplied;

  vector<double> dualRowLowerImplied;
  vector<double> dualRowUpperImplied;
  vector<double> dualColLowerImplied;
  vector<double> dualColUpperImplied;

  // Copy of original bounds when over-written using implied bounds
  // from presolve
  vector<double> SvColLower;
  vector<double> SvColUpper;
  vector<double> SvRowLower;
  vector<double> SvRowUpper;

};

/*
void getSolutionFromHModel(const HModel& model, HighsSolution& solution);
*/

#endif /* SIMPLEX_HMODEL_H_ */
