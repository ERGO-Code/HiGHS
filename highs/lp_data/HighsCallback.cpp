/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsCallback.cpp
 * @brief
 */
#include "lp_data/HighsCallback.h"

#include <cassert>

void HighsCallback::clearHighsCallbackDataOut() {
  this->data_out.log_type = HighsLogType::kInfo;
  this->data_out.running_time = -1;
  this->data_out.simplex_iteration_count = -1;
  this->data_out.ipm_iteration_count = -1;
  this->data_out.pdlp_iteration_count = -1;
  this->data_out.objective_function_value = -kHighsInf;
  this->data_out.mip_node_count = -1;
  this->data_out.mip_primal_bound = kHighsInf;
  this->data_out.mip_dual_bound = -kHighsInf;
  this->data_out.mip_gap = -1;
  this->data_out.mip_solution.clear();
  this->data_out.cutpool_start.clear();
  this->data_out.cutpool_index.clear();
  this->data_out.cutpool_value.clear();
  this->data_out.cutpool_lower.clear();
  this->data_out.cutpool_upper.clear();
  this->data_out.user_solution_callback_origin =
      userMipSolutionCallbackOrigin::kUserMipSolutionCallbackOriginAfterSetup;
}

void HighsCallback::clearHighsCallbackDataIn() {
  this->data_in.user_interrupt = false;
  this->data_in.user_solution.clear();
}

void HighsCallback::clear() {
  this->user_callback = nullptr;
  this->user_callback_data = nullptr;
  this->active.assign(kNumCallbackType, false);
  this->clearHighsCallbackDataOut();
  this->clearHighsCallbackDataIn();
}

bool HighsCallback::callbackActive(const int callback_type) {
  // Check that callback function has been defined
  if (!this->user_callback) return false;
  // Check that callback_type is within range
  const bool callback_type_ok =
      callback_type >= kCallbackMin && callback_type <= kCallbackMax;
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
  if (callback_type == kCallbackMipImprovingSolution ||
      callback_type == kCallbackMipSolution ||
      callback_type == kCallbackMipLogging ||
      callback_type == kCallbackMipGetCutPool ||
      callback_type == kCallbackMipDefineLazyConstraints ||
      callback_type == kCallbackMipUserSolution)
    assert(!action);
  return action;
}

// Conversions for C API

// Convert HighsCallbackDataOut to HighsCCallbackDataOut
HighsCallbackDataOut::operator HighsCCallbackDataOut() const {
  HighsCCallbackDataOut c_data_out;
  c_data_out.log_type = static_cast<int>(log_type);
  c_data_out.running_time = running_time;
  c_data_out.simplex_iteration_count = simplex_iteration_count;
  c_data_out.ipm_iteration_count = ipm_iteration_count;
  c_data_out.pdlp_iteration_count = pdlp_iteration_count;
  c_data_out.objective_function_value = objective_function_value;

  c_data_out.mip_node_count = mip_node_count;
  c_data_out.mip_total_lp_iterations = mip_total_lp_iterations;
  c_data_out.mip_primal_bound = mip_primal_bound;
  c_data_out.mip_dual_bound = mip_dual_bound;
  c_data_out.mip_gap = mip_gap;
  c_data_out.mip_solution_size = mip_solution.size();
  c_data_out.mip_solution =
      mip_solution.empty() ? nullptr : const_cast<double*>(mip_solution.data());

  c_data_out.cutpool_num_col = cutpool_num_col;
  c_data_out.cutpool_num_cut = cutpool_lower.size();
  c_data_out.cutpool_num_nz = cutpool_value.size();
  c_data_out.cutpool_start = cutpool_start.empty()
                                 ? nullptr
                                 : const_cast<HighsInt*>(cutpool_start.data());
  c_data_out.cutpool_index = cutpool_index.empty()
                                 ? nullptr
                                 : const_cast<HighsInt*>(cutpool_index.data());
  c_data_out.cutpool_value = cutpool_value.empty()
                                 ? nullptr
                                 : const_cast<double*>(cutpool_value.data());
  c_data_out.cutpool_lower = cutpool_lower.empty()
                                 ? nullptr
                                 : const_cast<double*>(cutpool_lower.data());
  c_data_out.cutpool_upper = cutpool_upper.empty()
                                 ? nullptr
                                 : const_cast<double*>(cutpool_upper.data());

  c_data_out.user_solution_callback_origin =
      static_cast<HighsInt>(user_solution_callback_origin);
  return c_data_out;
}

HighsCallbackDataIn HighsCallbackDataIn::operator=(
    const HighsCCallbackDataIn& data_in) {
  user_interrupt = data_in.user_interrupt != 0;
  user_solution.clear();

  // copy data from callback
  if (data_in.user_solution != nullptr) {
    user_solution.resize(data_in.user_solution_size);

    if (data_in.user_solution_size > 0) {
      user_solution.assign(data_in.user_solution,
                           data_in.user_solution + data_in.user_solution_size);
    }
  }

  return *this;
}
