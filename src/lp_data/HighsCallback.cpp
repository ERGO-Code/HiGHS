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

#include <cassert>

void HighsCallbackDataOut::clear() {
  this->log_type = HighsLogType::kInfo;
  this->simplex_iteration_count = -1;
  this->objective_solution.clear();
}

void HighsCallbackDataIn::clear() { this->user_interrupt = false; }

void HighsCallback::clear() {
  this->user_callback = nullptr;
  this->user_callback_data = nullptr;
  this->active.assign(kNumHighsCallbackType, false);
  this->data_out.clear();
  this->data_in.clear();
}

bool HighsCallback::callbackAction(const int callback_type,
                                   std::string message) {
  const bool callback_type_ok =
      callback_type >= kHighsCallbackMin && callback_type <= kHighsCallbackMax;
  assert(callback_type_ok);
  if (!callback_type_ok) return false;
  if (!this->active[callback_type]) return false;
  this->user_callback(callback_type, message.c_str(), this->user_callback_data,
                      this->data_out, this->data_in);
  if (callback_type == kHighsCallbackLogging) {
    assert(1 == 0);
    return false;
  } else if (callback_type == kHighsCallbackInterrupt) {
    return this->data_in.user_interrupt;
  } else {
    assert(1 == 0);
    return false;
  }
}
