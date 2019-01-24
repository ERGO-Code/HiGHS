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

enum SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY {
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG = 0,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEFAULT = SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH
};

enum SIMPLEX_PRICE_STRATEGY {
  SIMPLEX_PRICE_STRATEGY_COL = 0,
  SIMPLEX_PRICE_STRATEGY_ROW,
  SIMPLEX_PRICE_STRATEGY_ROW_SWITCH,
  SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH,
  SIMPLEX_PRICE_STRATEGY_ROW_ULTRA,
  SIMPLEX_PRICE_STRATEGY_DEFAULT = SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH
};

#endif /* SIMPLEX_SIMPLEXCONST_H_ */
