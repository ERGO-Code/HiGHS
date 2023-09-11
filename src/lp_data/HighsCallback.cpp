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
  this->node_count = -1;
  this->running_time = -1;
  this->primal_bound = kHighsInf;
  this->dual_bound = -kHighsInf;
  this->mip_rel_gap = -1;
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

bool HighsCallback::callbackActive(const int callback_type) {
  // Check that callback function has been defined
  if (!this->user_callback) return false;
  // Check that callback_type is within range
  const bool callback_type_ok =
      callback_type >= kHighsCallbackMin && callback_type <= kHighsCallbackMax;
  assert(callback_type_ok);
  if (!callback_type_ok) return false;
  // Don't call callback if it is not active
  assert(this->active.size() > 0);
  if (!this->active[callback_type]) return false;
  return true;
}

bool HighsCallback::callbackAction(const int callback_type,
                                   std::string message) {
  if (!callbackActive(callback_type)) return false;
  this->user_callback(callback_type, message.c_str(), &this->data_out,
                      &this->data_in, this->user_callback_data);
  // Assess any action
  bool action = this->data_in.user_interrupt;

  // Check for no action if case not handled internally
  if (callback_type == kHighsCallbackMipImprovingSolution ||
      callback_type == kHighsCallbackMipLogging)
    assert(!action);
  return action;
}
