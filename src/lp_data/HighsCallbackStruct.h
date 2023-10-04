/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2023 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsCallbackStruct.h
 * @brief
 */
#ifndef LP_DATA_HIGHSCALLBACKSTRUCT_H_
#define LP_DATA_HIGHSCALLBACKSTRUCT_H_

#include "util/HighsInt.h"

struct HighsCallbackDataOut {
  int log_type;  // cast of HighsLogType
  double running_time;
  HighsInt simplex_iteration_count;
  HighsInt ipm_iteration_count;
  double objective_function_value;
  int64_t mip_node_count;
  double mip_primal_bound;
  double mip_dual_bound;
  double mip_gap;
  double* mip_solution;
};

struct HighsCallbackDataIn {
  int user_interrupt;
};

#endif /* LP_DATA_HIGHSCALLBACKSTRUCT_H_ */
