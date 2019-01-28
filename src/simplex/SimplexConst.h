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

enum class SimplexStrategy {
  DUAL_PLAIN = 0,
  DUAL_TASKS,
  DUAL_MULTI,
  PRIMAL,
  DEFAULT = DUAL_PLAIN
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

enum class SimplexPriceStrategy {
  COL = 0,
  ROW,
  ROW_SWITCH,
  ROW_SWITCH_COL_SWITCH,
  ROW_ULTRA,
  DEFAULT = ROW_SWITCH_COL_SWITCH
};

const int invertHint_no = 0;
const int invertHint_updateLimitReached = 1;
const int invertHint_syntheticClockSaysInvert = 2;
const int invertHint_possiblyOptimal = 3;
const int invertHint_possiblyPrimalUnbounded = 4;
const int invertHint_possiblyDualUnbounded = 5;
const int invertHint_possiblySingularBasis = 6;
const int invertHint_primalInfeasibleInPrimalSimplex = 7;
const int invertHint_chooseColumnFail = 8;



#endif /* SIMPLEX_SIMPLEXCONST_H_ */
