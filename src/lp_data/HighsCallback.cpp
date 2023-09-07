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
/**@file lp_data/HighsCallback.cpp
 * @brief
 */
#include "HighsCallback.h"
void HighsCallbackDataOut::clear() {
  this->log_type = HighsLogType::kInfo;
  this->simplex_iteration_count = -1;
  this->objective_solution.clear();
}

void HighsCallbackDataIn::clear() { this->user_interrupt = false; }

void HighsCallback::clear() {
  this->highs_user_callback = nullptr;
  this->highs_user_callback_data = nullptr;
  this->active.assign(num_type, false);
  this->highs_callback_data_out.clear();
  this->highs_callback_data_in.clear();
}

