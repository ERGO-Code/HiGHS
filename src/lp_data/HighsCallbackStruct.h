/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
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

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int log_type;  // cast of HighsLogType
  double running_time;
  HighsInt simplex_iteration_count;
  HighsInt ipm_iteration_count;
  HighsInt pdlp_iteration_count;
  double objective_function_value;
  int64_t mip_node_count;
  double mip_primal_bound;
  double mip_dual_bound;
  double mip_gap;
  const double* mip_solution;
  HighsInt cutpool_num_col;
  HighsInt cutpool_num_cut;
  HighsInt cutpool_num_nz;
  const HighsInt* cutpool_ARstart;
  const HighsInt* cutpool_ARindex;
  const double* cutpool_ARvalue;
  const double* cutpool_lower;
  const double* cutpool_upper;
} HighsCallbackDataOut;

typedef struct {
  int user_interrupt;
  HighsInt cutset_num_cut;
  HighsInt* cutset_ARstart;
  HighsInt* cutset_ARindex;
  double* cutset_ARvalue;
  double* cutset_lower;
  double* cutset_upper;
} HighsCallbackDataIn;

// Additional callback handling
typedef void (*HighsCCallbackType)(int, const char*,
                                   const HighsCallbackDataOut*,
                                   HighsCallbackDataIn*, void*);

#ifdef __cplusplus
}
#endif

#endif /* LP_DATA_HIGHSCALLBACKSTRUCT_H_ */
