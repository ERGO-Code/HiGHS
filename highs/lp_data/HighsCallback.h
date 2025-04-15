/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsCallback.h
 * @brief
 */
#ifndef LP_DATA_HIGHSCALLBACK_H_
#define LP_DATA_HIGHSCALLBACK_H_

#include <functional>

// forward declaration to avoid circular dependency
enum class HighsStatus;
class Highs;

#include "lp_data/HStruct.h"
#include "lp_data/HighsCallbackStruct.h"

enum userMipSolutionCallbackOrigin {
  kUserMipSolutionCallbackOriginAfterSetup = 0,
  kUserMipSolutionCallbackOriginBeforeDive,
  kUserMipSolutionCallbackOriginEvaluateRootNode0,
  kUserMipSolutionCallbackOriginEvaluateRootNode1,
  kUserMipSolutionCallbackOriginEvaluateRootNode2,
  kUserMipSolutionCallbackOriginEvaluateRootNode3,
  kUserMipSolutionCallbackOriginEvaluateRootNode4
};

/**
 * Struct to handle callback output data
 */
struct HighsCallbackOutput {
  Highs* highs = nullptr;
  HighsLogType log_type;
  double running_time;
  HighsInt simplex_iteration_count;
  HighsInt ipm_iteration_count;
  HighsInt pdlp_iteration_count;
  double objective_function_value;
  int64_t mip_node_count;
  int64_t mip_total_lp_iterations;
  double mip_primal_bound;
  double mip_dual_bound;
  double mip_gap;
  std::vector<double> mip_solution;
  HighsInt cutpool_num_col;
  HighsInt cutpool_num_cut;
  std::vector<HighsInt> cutpool_start;
  std::vector<HighsInt> cutpool_index;
  std::vector<double> cutpool_value;
  std::vector<double> cutpool_lower;
  std::vector<double> cutpool_upper;
  userMipSolutionCallbackOrigin user_solution_callback_origin;

  operator HighsCallbackDataOut() const;
};

struct HighsCallbackInput {
  Highs* highs = nullptr;
  bool user_interrupt = false;
  bool user_has_solution = false;
  std::vector<double> user_solution;

  HighsStatus setSolution(HighsInt num_entries, const double* value);

  HighsStatus setSolution(HighsInt num_entries, const HighsInt* index,
                          const double* value);

  HighsStatus repairSolution();

  operator HighsCallbackDataIn() const;
  HighsCallbackInput operator=(const HighsCallbackDataIn& data_in);
};

using HighsCallbackFunctionType =
    std::function<void(int, const std::string&, const HighsCallbackOutput*,
                       HighsCallbackInput*, void*)>;

struct HighsCallback {
  HighsCallback(Highs* highs) : highs(highs) { 
    data_out.highs = highs;
    data_in.highs = highs;
    clear();
  }

  // Function pointers cannot be used for Pybind11, so use std::function
  HighsCallbackFunctionType user_callback = nullptr;
  HighsCCallbackType c_callback = nullptr;
  void* user_callback_data = nullptr;
  Highs* highs = nullptr;
  std::vector<bool> active;
  HighsCallbackOutput data_out;
  HighsCallbackInput data_in;
  bool callbackActive(const int callback_type);
  bool callbackAction(const int callback_type, std::string message = "");
  void clearHighsCallbackOutput();
  void clearHighsCallbackInput();
  void clear();
};
#endif /* LP_DATA_HIGHSCALLBACK_H_ */
