/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/SimplexConst.h
 * @brief Constants for HiGHS simplex solvers
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_SIMPLEXCONST_H_
#define SIMPLEX_SIMPLEXCONST_H_

enum class SimplexSolutionStatus {
  UNSET = -1,
  OPTIMAL,
  PRIMAL_FEASIBLE,
  DUAL_FEASIBLE,
  INFEASIBLE,
  UNBOUNDED,
  SINGULAR,
  FAILED,
  REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND,
  OUT_OF_TIME
};

enum class SimplexStrategy {
  CHOOSE = 0,
  DUAL,
  DUAL_PLAIN = DUAL,
  DUAL_TASKS,
  DUAL_MULTI,
  PRIMAL,
  DEFAULT = DUAL
};

enum class SimplexDualiseStrategy { OFF = 0, CHOOSE, ON, DEFAULT = OFF };

enum class SimplexPermuteStrategy { OFF = 0, CHOOSE, ON, DEFAULT = OFF };

enum class SimplexScaleStrategy {
  OFF = 0,
  CHOOSE,
  HSOL,
  HIGHS,
  DEFAULT = HIGHS
};

enum class SimplexCrashStrategy {
  OFF = 0,
  LTSSF_K,
  LTSSF_PRI,
  LTSF_K,
  LTSF_PRI,
  LTSF,
  BIXBY,
  BIXBY_NO_NONZERO_COL_COSTS,
  BASIC,
  TEST_SING,
  DEFAULT = OFF
};

enum class SimplexDualEdgeWeightStrategy {
  DANTZIG = 0,
  DEVEX,
  STEEPEST_EDGE,
  STEEPEST_EDGE_UNIT_INITIAL,
  STEEPEST_EDGE_TO_DEVEX_SWITCH,
  DEFAULT = STEEPEST_EDGE_TO_DEVEX_SWITCH
};

enum class SimplexPrimalEdgeWeightStrategy {
  DANTZIG = 0,
  DEVEX,
  DEFAULT = DEVEX
};

enum class SimplexPriceStrategy {
  COL = 0,
  ROW,
  ROW_SWITCH,
  ROW_SWITCH_COL_SWITCH,
  ROW_ULTRA,
  DEFAULT = ROW_SWITCH_COL_SWITCH
};

// Not an enum class since invert_hint is used in so many places
enum InvertHint {
  INVERT_HINT_NO = 0,
  INVERT_HINT_UPDATE_LIMIT_REACHED,
  INVERT_HINT_SYNTHETIC_CLOCK_SAYS_INVERT,
  INVERT_HINT_POSSIBLY_OPTIMAL,
  INVERT_HINT_POSSIBLY_PRIMAL_UNBOUNDED,
  INVERT_HINT_POSSIBLY_DUAL_UNBOUNDED,
  INVERT_HINT_POSSIBLY_SINGULAR_BASIS,
  INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX,
  INVERT_HINT_CHOOSE_COLUMN_FAIL,
  INVERT_HINT_Count
};

// TODO: Set this false tactically to make mip interface more
// efficient by preventing reinversion on optimality in phase 1 or
// phase 2
const bool invert_if_row_out_negative = true;

enum class FeasibilityStrategy {
  kApproxComponentWise,
  kApproxExact,
  kDirectSolve
};

#endif /* SIMPLEX_SIMPLEXCONST_H_ */
