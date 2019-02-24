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

const int HIGHS_CONST_I_INF = 2147483647;//32767;
const double HIGHS_CONST_INF = 1e200;
const double HIGHS_CONST_TINY = 1e-14;
const double HIGHS_CONST_ZERO = 1e-50;

constexpr double kBoundTolerance = 1e-8;

enum ModelLogLevel {
  ML_NONE = 0,
  ML_VERBOSE = 1,
  ML_DETAILED = 2,
  ML_MINIMAL = 4,
  ML_DEFAULT = ML_MINIMAL,
  ML_ALWAYS = ML_VERBOSE | ML_DETAILED | ML_MINIMAL
};
  
enum class ParallelOption {
  OFF = 0,
  ON,
  DEFAULT = OFF
};

enum class PresolveOption {
  OFF = 0,
  ON,
  DEFAULT = OFF
};

enum class CrashOption {
  OFF = 0,
  ON,
  DEFAULT = OFF
};
  
enum class SimplexOption {
  OFF = 0,
  ON,
  DEFAULT = OFF
};

const double INFINITE_COST_DEFAULT = 1e20;
const double INFINITE_BOUND_DEFAULT = 1e20;
const double SMALL_MATRIX_VALUE_DEFAULT = 1e-9;
const double LARGE_MATRIX_VALUE_DEFAULT = 1e15;
const double HIGHS_RUN_TIME_LIMIT_DEFAULT = HIGHS_CONST_INF;
const double PRIMAL_FEASIBILITY_TOLERANCE_DEFAULT = 1e-7;
const double DUAL_FEASIBILITY_TOLERANCE_DEFAULT = 1e-7;
const double DUAL_OBJECTIVE_VALUE_UPPER_BOUND_DEFAULT = HIGHS_CONST_INF;
const int SIMPLEX_ITERATION_LIMIT_DEFAULT = HIGHS_CONST_I_INF;
const int SIMPLEX_UPDATE_LIMIT_DEFAULT = 5000;

/** SCIP-like basis status for columns and rows. */
enum class HighsBasisStatus {
  LOWER = 0, // (slack) variable is at its lower bound [including fixed variables]
  BASIC, // (slack) variable is basic 
  UPPER, // (slack) variable is at its upper bound 
  ZERO,  // free variable is non-basic and set to zero 
  SUPER // Super-basic variable: non-basic and either free and nonzero
	// or not at a bound. Not permitted when allow_superbasic is
	// false: no SCIP equivalent
};

/** HiGHS nonbasicFlag status for columns and rows. Don't use enum
    class since they are used as int to replace conditional statements
    by multiplication */
const int NONBASIC_FLAG_TRUE = 1;  // Nonbasic
const int NONBASIC_FLAG_FALSE = 0;  // Basic

/** HiGHS nonbasicMove status for columns and rows. Don't use enum
    class since they are used in conditional statements */
const int NONBASIC_MOVE_UP = 1;   // Free to move (only) up
const int NONBASIC_MOVE_DN = -1;  // Free to move (only) down
const int NONBASIC_MOVE_ZE = 0;    // Fixed or free to move up and down



#endif /* LP_DATA_HCONST_H_ */
