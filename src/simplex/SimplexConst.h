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
/**@file lp_data/SimplexConst.h
 * @brief Constants for HiGHS simplex solvers
 */
#ifndef SIMPLEX_SIMPLEXCONST_H_
#define SIMPLEX_SIMPLEXCONST_H_

enum class SimplexAlgorithm { kPrimal = 0, kDual };

enum SimplexStrategy {
  kSimplexStrategyMin = 0,
  kSimplexStrategyChoose = kSimplexStrategyMin,      // 0
  kSimplexStrategyDual,                              // 1
  kSimplexStrategyDualPlain = kSimplexStrategyDual,  // 1
  kSimplexStrategyDualTasks,                         // 2
  kSimplexStrategyDualMulti,                         // 3
  kSimplexStrategyPrimal,                            // 4
  kSimplexStrategyMax = kSimplexStrategyPrimal,
  kSimplexStrategyNum
};

enum SimplexSolvePhase {
  SOLVE_PHASE_MIN = -3,
  SOLVE_PHASE_ERROR = SOLVE_PHASE_MIN,  // -3
  SOLVE_PHASE_EXIT,                     // -2,
  SOLVE_PHASE_UNKNOWN,                  // -1
  SOLVE_PHASE_OPTIMAL,                  // 0
  SOLVE_PHASE_1,                        // 1
  SOLVE_PHASE_2,                        // 2
  SOLVE_PHASE_CLEANUP = 4,
  SOLVE_PHASE_MAX = SOLVE_PHASE_CLEANUP
};

enum SimplexScaleStrategy {
  SIMPLEX_SCALE_STRATEGY_MIN = 0,
  SIMPLEX_SCALE_STRATEGY_OFF = SIMPLEX_SCALE_STRATEGY_MIN,  // 0
  SIMPLEX_SCALE_STRATEGY_HIGHS,                             // 1
  SIMPLEX_SCALE_STRATEGY_HIGHS_FORCED,                      // 2
  SIMPLEX_SCALE_STRATEGY_015,                               // 3
  SIMPLEX_SCALE_STRATEGY_0157,                              // 4
  SIMPLEX_SCALE_STRATEGY_MAX = SIMPLEX_SCALE_STRATEGY_0157
};

enum SimplexCrashStrategy {
  SIMPLEX_CRASH_STRATEGY_MIN = 0,
  SIMPLEX_CRASH_STRATEGY_OFF = SIMPLEX_CRASH_STRATEGY_MIN,
  SIMPLEX_CRASH_STRATEGY_LTSSF_K,
  SIMPLEX_CRASH_STRATEGY_LTSSF = SIMPLEX_CRASH_STRATEGY_LTSSF_K,
  SIMPLEX_CRASH_STRATEGY_BIXBY,
  SIMPLEX_CRASH_STRATEGY_LTSSF_PRI,
  SIMPLEX_CRASH_STRATEGY_LTSF_K,
  SIMPLEX_CRASH_STRATEGY_LTSF_PRI,
  SIMPLEX_CRASH_STRATEGY_LTSF,
  SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS,
  SIMPLEX_CRASH_STRATEGY_BASIC,
  SIMPLEX_CRASH_STRATEGY_TEST_SING,
  SIMPLEX_CRASH_STRATEGY_MAX = SIMPLEX_CRASH_STRATEGY_TEST_SING
};

enum SimplexDualEdgeWeightStrategy {
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MIN = -1,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_CHOOSE =
      SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MIN,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MAX =
      SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL
};

enum SimplexPrimalEdgeWeightStrategy {
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MIN = -1,
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_CHOOSE =
      SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MIN,
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG,
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX,
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MAX =
      SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX
};

enum SimplexPriceStrategy {
  SIMPLEX_PRICE_STRATEGY_MIN = 0,
  SIMPLEX_PRICE_STRATEGY_COL = SIMPLEX_PRICE_STRATEGY_MIN,
  SIMPLEX_PRICE_STRATEGY_ROW,
  SIMPLEX_PRICE_STRATEGY_ROW_SWITCH,
  SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH,
  SIMPLEX_PRICE_STRATEGY_MAX = SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH
};

enum SimplexDualChuzcStrategy {
  SIMPLEX_DUAL_CHUZC_STRATEGY_MIN = 0,
  SIMPLEX_DUAL_CHUZC_STRATEGY_CHOOSE = SIMPLEX_DUAL_CHUZC_STRATEGY_MIN,
  SIMPLEX_DUAL_CHUZC_STRATEGY_QUAD,
  SIMPLEX_DUAL_CHUZC_STRATEGY_HEAP,
  SIMPLEX_DUAL_CHUZC_STRATEGY_BOTH,
  SIMPLEX_DUAL_CHUZC_STRATEGY_MAX = SIMPLEX_DUAL_CHUZC_STRATEGY_BOTH
};

enum SimplexPrimalCorrectionStrategy {
  SIMPLEX_PRIMAL_CORRECTION_STRATEGY_NONE = 0,
  SIMPLEX_PRIMAL_CORRECTION_STRATEGY_IN_REBUILD,
  SIMPLEX_PRIMAL_CORRECTION_STRATEGY_ALWAYS,
  SIMPLEX_PRIMAL_CORRECTION_STRATEGY_REFINED
};

// Not an enum class since rebuild_reason is used in so many places
enum RebuildReason {
  REBUILD_REASON_NO = 0,
  REBUILD_REASON_UPDATE_LIMIT_REACHED,                 // 1
  REBUILD_REASON_SYNTHETIC_CLOCK_SAYS_INVERT,          // 2
  REBUILD_REASON_POSSIBLY_OPTIMAL,                     // 3
  REBUILD_REASON_POSSIBLY_PRIMAL_UNBOUNDED,            // 4
  REBUILD_REASON_POSSIBLY_DUAL_UNBOUNDED,              // 5
  REBUILD_REASON_POSSIBLY_SINGULAR_BASIS,              // 6
  REBUILD_REASON_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX,  // 7
  REBUILD_REASON_CHOOSE_COLUMN_FAIL,                   // 8
  REBUILD_REASON_Count
};

enum class DualEdgeWeightMode { kDantzig = 0, kDevex, kSteepestEdge, kCount };

const HighsInt kDualTasksMinThreads = 3;
const HighsInt kDualMultiMinThreads = 1;  // 2;

// Simplex nonbasicFlag status for columns and rows. Don't use enum
// class since they are used as HighsInt to replace conditional
// statements by multiplication
const HighsInt kNonbasicFlagTrue = 1;   // Nonbasic
const HighsInt kNonbasicFlagFalse = 0;  // Basic

// Simplex nonbasicMove status for columns and rows. Don't use enum
// class since they are used in conditional statements
const HighsInt kNonbasicMoveUp = 1;   // Free to move (only) up
const HighsInt kNonbasicMoveDn = -1;  // Free to move (only) down
const HighsInt kNonbasicMoveZe = 0;   // Fixed or free to move up and down
const HighsInt kIllegalMoveValue =
    -99;  // Used to see whether valid move value has been set

// Threshold for accepting updated DSE weight
const double kAcceptDseWeightThreshold = 0.25;

//
// Relation between HiGHS basis and Simplex basis
//
// Data structures
// ===============
//
// HiGHS basis consists of vectors
//
// * col_status[numCol]
// * row_status[numRow]
//
// Simplex basis consists of vectors
//
// * nonbasicMove[numTot]
// * basicIndex[numRow]
// * nonbasicFlag[numTot]
//
// where nonbasicFlag is duplicate information but is used to identify
// whether a particular variable is basic or nonbasic.
//
// Basic variables
// ===============
//
// Highs: *_status value of BASIC
//
// <=>
//
// Simplex: nonbasicFlag value of kNonbasicFlagFalse
//
// Nonbasic variables
// ==================
//
// Relations complicated by the fact that
//
// * HiGHS   rows have bounds [ l,  u]
// * Simplex rows have bounds [-u, -l]
//
// Nonbasic columns
// ================
//
// Highs: col_status value of LOWER - at lower bound
// <=>
// Simplex: nonbasicMove value of kNonbasicMoveUp - [l, Inf] column free to
// move up and negative dual
//
// Highs: col_status value of ZERO - at zero
// =>
// Simplex: nonbasicMove value of kNonbasicMoveZe - free variable treated
// specially in simplex
//
// Highs: col_status value of UPPER - at upper bound
// =>
// Simplex: Either
// * nonbasicMove value of kNonbasicMoveDn - [-Inf, u] column free to move down
// and positive dual
// * nonbasicMove value of kNonbasicMoveZe - [   l, u] column ?? and free dual
//
// Simplex: nonbasicMove value of kNonbasicMoveDn - [-Inf, u] column free to
// move down and positive dual
// =>
// Highs: col_status value of UPPER - at upper bound
//
// Simplex: nonbasicMove value of kNonbasicMoveZe - [l, u] column ?? and free
// dual
// =>
// Highs: Either
// * col_status value of UPPER - at upper bound
// * col_status value of ZERO - at zero
//
// Nonbasic rows
// =============
//
#endif /* SIMPLEX_SIMPLEXCONST_H_ */
