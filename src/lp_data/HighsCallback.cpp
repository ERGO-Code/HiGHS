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
  this->objective = -kHighsInf;
  this->col_value = nullptr;
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
  // Check that callback_type is within range
  bool action = false;
  const bool callback_type_ok =
      callback_type >= kHighsCallbackMin && callback_type <= kHighsCallbackMax;
  assert(callback_type_ok);
  if (!callback_type_ok) return action;
  // Don't call callback if it is not active
  if (!this->active[callback_type]) return action;
  // Call callback!
  this->user_callback(callback_type, message.c_str(), &this->data_out,
                      &this->data_in, this->user_callback_data);
  // Assess any action
  if (callback_type == kHighsCallbackInterrupt ||
      callback_type == kHighsCallbackMipDualBound) 
    action = this->data_in.user_interrupt;
  return action;
}
