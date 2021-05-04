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
/**@file lp_data/HighsInfoDebug.cpp
 * @brief
 */
#include "lp_data/HighsInfoDebug.h"

HighsDebugStatus debugInfo(const HighsOptions& options, const HighsLp& lp,
                           const HighsBasis& basis,
                           const HighsSolution& solution, const HighsInfo& info,
                           const HighsModelStatus model_status) {
  //  if (options.highs_debug_level < kHighsDebugLevelCostly) return
  //  HighsDebugStatus::kNotChecked;
  //  printf("debugInfo!\n");
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  if (!info.valid) return return_status;

  //  switch (model_status_) {
  
  return return_status;
}
