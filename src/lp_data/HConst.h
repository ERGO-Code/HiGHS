/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HConst.h
 * @brief Constants for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HCONST_H_
#define LP_DATA_HCONST_H_

#include <string>

const int HIGHS_CONST_I_INF = 2147483647;  // 32767;
const double HIGHS_CONST_INF = 1e200;
const double HIGHS_CONST_TINY = 1e-14;
const double HIGHS_CONST_ZERO = 1e-50;

constexpr double kBoundTolerance = 1e-8;

enum HighsPrintMessageLevel {
  ML_NONE = 0,
  ML_VERBOSE = 1,
  ML_DETAILED = 2,
  ML_MINIMAL = 4,
  ML_DEFAULT = ML_MINIMAL,
  ML_ALWAYS = ML_VERBOSE | ML_DETAILED | ML_MINIMAL
};

enum class HighsOptionType { BOOL = 0, INT, DOUBLE, STRING};

enum PresolveOption { PRESOLVE_OPTION_OFF = 0, PRESOLVE_OPTION_ON, PRESOLVE_OPTION_DEFAULT = PRESOLVE_OPTION_ON };

enum ParallelOption { PARALLEL_OPTION_OFF = 0, PARALLEL_OPTION_ON, PARALLEL_OPTION_DEFAULT = PARALLEL_OPTION_OFF };

enum CrashOption { CRASH_OPTION_OFF = 0, CRASH_OPTION_ON, CRASH_OPTION_DEFAULT = CRASH_OPTION_OFF };

enum SimplexOption { SIMPLEX_OPTION_OFF = 0, SIMPLEX_OPTION_ON, SIMPLEX_OPTION_DEFAULT = SIMPLEX_OPTION_ON };

enum IpmOption { IPM_OPTION_OFF = 0, IPM_OPTION_ON, IPM_OPTION_DEFAULT = IPM_OPTION_OFF };

const std::string FILENAME_DEFAULT = "";
const std::string OPTIONS_FILE_DEFAULT = "";
const int RUN_AS_HSOL_MIN = 0;
const int RUN_AS_HSOL_DEFAULT = 0;
const int RUN_AS_HSOL_MAX = 1;
const int KEEP_N_ROWS_DELETE_ROWS = -1;
const int KEEP_N_ROWS_DELETE_ENTRIES = 0;
const int KEEP_N_ROWS_KEEP_ROWS = 1;
const int KEEP_N_ROWS_MIN = KEEP_N_ROWS_DELETE_ROWS;
const int KEEP_N_ROWS_DEFAULT = KEEP_N_ROWS_DELETE_ROWS;
const int KEEP_N_ROWS_MAX = KEEP_N_ROWS_KEEP_ROWS;
const double INFINITE_COST_MIN = 1e15;
const double INFINITE_COST_DEFAULT = 1e20;
const double INFINITE_COST_MAX = 1e25;
const double INFINITE_BOUND_MIN = 1e15;
const double INFINITE_BOUND_DEFAULT = 1e20;
const double INFINITE_BOUND_MAX = 1e25;
const double SMALL_MATRIX_VALUE_MIN = 1e-12;
const double SMALL_MATRIX_VALUE_DEFAULT = 1e-9;
const double SMALL_MATRIX_VALUE_MAX = HIGHS_CONST_INF;
const double LARGE_MATRIX_VALUE_MIN = 1e0;
const double LARGE_MATRIX_VALUE_DEFAULT = 1e15;
const double LARGE_MATRIX_VALUE_MAX = 1e20;
const int ALLOWED_SIMPLEX_SCALE_FACTOR_MIN = 0;
const int ALLOWED_SIMPLEX_SCALE_FACTOR_DEFAULT = 10;
const int ALLOWED_SIMPLEX_SCALE_FACTOR_MAX = 20;
const double HIGHS_RUN_TIME_LIMIT_DEFAULT = HIGHS_CONST_INF;
const double PRIMAL_FEASIBILITY_TOLERANCE_MIN = 1e-10;
const double PRIMAL_FEASIBILITY_TOLERANCE_DEFAULT = 1e-7;
const double PRIMAL_FEASIBILITY_TOLERANCE_MAX = HIGHS_CONST_INF;
const double DUAL_FEASIBILITY_TOLERANCE_MIN = 1e-10;
const double DUAL_FEASIBILITY_TOLERANCE_DEFAULT = 1e-7;
const double DUAL_FEASIBILITY_TOLERANCE_MAX = HIGHS_CONST_INF;
const double DUAL_OBJECTIVE_VALUE_UPPER_BOUND_DEFAULT = HIGHS_CONST_INF;
const int SIMPLEX_ITERATION_LIMIT_DEFAULT = HIGHS_CONST_I_INF;
const int SIMPLEX_UPDATE_LIMIT_DEFAULT = 5000;

const double SIMPLEX_INITIAL_CONDITION_TOLERANCE_MIN = 1e0;
const double SIMPLEX_INITIAL_CONDITION_TOLERANCE_DEFAULT = 1e14;
const double SIMPLEX_INITIAL_CONDITION_TOLERANCE_MAX = 1e24;

const double DUAL_STEEPEST_EDGE_WEIGHT_LOG_ERROR_THRESHHOLD_MIN = 1e0;
const double DUAL_STEEPEST_EDGE_WEIGHT_LOG_ERROR_THRESHHOLD_DEFAULT = 1e1;
const double DUAL_STEEPEST_EDGE_WEIGHT_LOG_ERROR_THRESHHOLD_MAX = HIGHS_CONST_INF;

/** SCIP/CPLEX-like HiGHS basis status for columns and rows. */
enum class HighsBasisStatus {
  LOWER =
      0,  // (slack) variable is at its lower bound [including fixed variables]
  BASIC,  // (slack) variable is basic
  UPPER,  // (slack) variable is at its upper bound
  ZERO,   // free variable is non-basic and set to zero
  NONBASIC, // nonbasic with no specific bound information - useful for users and postsolve
  SUPER   // Super-basic variable: non-basic and either free and
          // nonzero or not at a bound. No SCIP equivalent
};

/** Simplex nonbasicFlag status for columns and rows. Don't use enum
    class since they are used as int to replace conditional statements
    by multiplication */
const int NONBASIC_FLAG_TRUE = 1;   // Nonbasic
const int NONBASIC_FLAG_FALSE = 0;  // Basic

/** Simplex nonbasicMove status for columns and rows. Don't use enum
    class since they are used in conditional statements */
const int NONBASIC_MOVE_UP = 1;   // Free to move (only) up
const int NONBASIC_MOVE_DN = -1;  // Free to move (only) down
const int NONBASIC_MOVE_ZE = 0;   // Fixed or free to move up and down
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
// Simplex: nonbasicFlag value of NONBASIC_FLAG_FALSE
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
// Simplex: nonbasicMove value of NONBASIC_MOVE_UP - [l, Inf] column free to
// move up and negative dual
//
// Highs: col_status value of ZERO - at zero
// =>
// Simplex: nonbasicMove value of NONBASIC_MOVE_ZE - free variable treated
// specially in simplex
//
// Highs: col_status value of UPPER - at upper bound
// =>
// Simplex: Either
// * nonbasicMove value of NONBASIC_MOVE_DN - [-Inf, u] column free to move down
// and positive dual
// * nonbasicMove value of NONBASIC_MOVE_ZE - [   l, u] column ?? and free dual
//
// Simplex: nonbasicMove value of NONBASIC_MOVE_DN - [-Inf, u] column free to
// move down and positive dual
// =>
// Highs: col_status value of UPPER - at upper bound
//
// Simplex: nonbasicMove value of NONBASIC_MOVE_ZE - [l, u] column ?? and free
// dual
// =>
// Highs: Either
// * col_status value of UPPER - at upper bound
// * col_status value of ZERO - at zero
//
// Nonbasic rows
// =============
//
#endif /* LP_DATA_HCONST_H_ */
