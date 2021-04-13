/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
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

#include <limits>
#include <string>

#include "util/HighsInt.h"

const HighsInt kHighsIInf = std::numeric_limits<HighsInt>::max();
const double kHighsInf = std::numeric_limits<double>::infinity();
const double kHighsTiny = 1e-14;
const double kHighsZero = 1e-50;
const std::string kHighsOffString = "off";
const std::string kHighsChooseString = "choose";
const std::string kHighsOnString = "on";
const HighsInt kHighsThreadLimit = 8;  // 32;

enum HighsDebugLevel {
  kHighsDebugLevelMin = 0,
  kHighsDebugLevelNone = kHighsDebugLevelMin,  // 0
  kHighsDebugLevelCheap,                       // 1
  kHighsDebugLevelCostly,                      // 2
  kHighsDebugLevelExpensive,                   // 3
  kHighsDebugLevelMax = kHighsDebugLevelExpensive
};

enum class HighsDebugStatus {
  NOT_CHECKED = -1,
  OK,
  SMALL_ERROR,
  WARNING,
  LARGE_ERROR,
  ERROR,
  EXCESSIVE_ERROR,
  LOGICAL_ERROR,
};

enum HighsAnalysisLevel {
  HIGHS_ANALYSIS_LEVEL_MIN = 0,
  HIGHS_ANALYSIS_LEVEL_NONE = HIGHS_ANALYSIS_LEVEL_MIN,
  HIGHS_ANALYSIS_LEVEL_MODEL_DATA = 1,
  HIGHS_ANALYSIS_LEVEL_SOLVER_DATA = 2,
  HIGHS_ANALYSIS_LEVEL_SOLVER_TIME = 4,
  HIGHS_ANALYSIS_LEVEL_NLA_DATA = 8,
  HIGHS_ANALYSIS_LEVEL_NLA_TIME = 16,
  HIGHS_ANALYSIS_LEVEL_MAX =
      HIGHS_ANALYSIS_LEVEL_MODEL_DATA + HIGHS_ANALYSIS_LEVEL_SOLVER_DATA +
      HIGHS_ANALYSIS_LEVEL_SOLVER_TIME + HIGHS_ANALYSIS_LEVEL_NLA_DATA +
      HIGHS_ANALYSIS_LEVEL_NLA_TIME
};

enum class HighsVarType : uint8_t {
  CONTINUOUS = 0,
  INTEGER = 1,
  IMPLICIT_INTEGER = 2,
};

enum class HighsOptionType { BOOL = 0, INT, DOUBLE, STRING };

enum class HighsInfoType { INT = 1, DOUBLE };

enum OptionOffChooseOn { OPTION_OFF = -1, OPTION_CHOOSE, OPTION_ON };

/** SCIP/HiGHS Objective sense */
enum class ObjSense { MINIMIZE = 1, MAXIMIZE = -1 };
enum class MatrixOrientation { NONE = 0, COLWISE, ROWWISE };

enum SolverOption {
  SOLVER_OPTION_SIMPLEX = -1,
  SOLVER_OPTION_CHOOSE,
  SOLVER_OPTION_IPM
};

enum PrimalDualStatus {
  STATUS_NOTSET = -1,
  STATUS_MIN = STATUS_NOTSET,
  STATUS_NO_SOLUTION,
  STATUS_UNKNOWN,
  STATUS_INFEASIBLE_POINT,
  STATUS_FEASIBLE_POINT,
  STATUS_MAX = STATUS_FEASIBLE_POINT
};

const std::string FILENAME_DEFAULT = "";

// Need to allow infinite costs to pass SCIP LPI unit tests
const bool allow_infinite_costs = true;

// Primal/dual statuses and corresponding HighsModelStatus
// values. Note that if dual infeasibility is identified, then the
// prototype primal code is used to distinguish PRIMAL_DUAL_INFEASIBLE
// from PRIMAL_UNBOUNDED. If this fails, then HiGHS may just return
// DUAL_INFEASIBLE
//
//           | Du Infeas    | Du Feas   | Du UnBd
// Pr Infeas | PR_DU_INFEAS | PR_INFEAS | PR_INFEAS
// Pr Feas   | PR_UNBD      | OPTIMAL   |   N/A
// Pr Unbd   | PR_UNBD      |     N/A   |   N/A
//
// Dual infeasibility is recognised by infeasibility at dual phase 1 optimality
// (and implied by primal unboundedness)
//
// Dual feasibility is recognised by feasibility at dual phase 1 optimality or
// primal phase 2 optimality
//
// Dual unboundedness is recognised by unboundedness in dual phase 2
//
// Primal infeasibility is recognised by infeasibility at primal phase 1
// optimality (and implied by dual unboundedness)
//
// Primal feasibility is recognised by feasibility at primal phase 1 optimality
// or dual phase 2 optimality
//
// Primal unboundedness is recognised by unboundedness in primal phase 2
//

enum class HighsModelStatus {
  // NB Add new status values to the end so that int cast of status
  // values is unchanged, since enums are not preserved in some
  // interfaces
  NOTSET = 0,
  HIGHS_MODEL_STATUS_MIN = NOTSET,
  LOAD_ERROR,
  MODEL_ERROR,
  PRESOLVE_ERROR,
  SOLVE_ERROR,
  POSTSOLVE_ERROR,
  MODEL_EMPTY,
  OPTIMAL,
  PRIMAL_INFEASIBLE,
  PRIMAL_INFEASIBLE_OR_UNBOUNDED,
  PRIMAL_UNBOUNDED,
  REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND,
  REACHED_TIME_LIMIT,
  REACHED_ITERATION_LIMIT,
  PRIMAL_DUAL_INFEASIBLE,
  DUAL_INFEASIBLE,
  HIGHS_MODEL_STATUS_MAX = DUAL_INFEASIBLE
};

/** SCIP/CPLEX-like HiGHS basis status for columns and rows. */
enum class HighsBasisStatus {
  LOWER =
      0,  // (slack) variable is at its lower bound [including fixed variables]
  BASIC,  // (slack) variable is basic
  UPPER,  // (slack) variable is at its upper bound
  ZERO,   // free variable is non-basic and set to zero
  NONBASIC  // nonbasic with no specific bound information - useful for users
            // and postsolve
};

// Illegal values of num/max/sum infeasibility - used to indicate that true
// values aren't known
const HighsInt illegal_infeasibility_count = -1;
const double illegal_infeasibility_measure = -1;

#endif /* LP_DATA_HCONST_H_ */
