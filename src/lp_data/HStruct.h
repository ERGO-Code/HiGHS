/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HStruct.h
 * @brief Structs for HiGHS
 */
#ifndef LP_DATA_HSTRUCT_H_
#define LP_DATA_HSTRUCT_H_

#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"

struct HighsIterationCounts {
  HighsInt simplex = 0;
  HighsInt ipm = 0;
  HighsInt crossover = 0;
  HighsInt qp = 0;
};

struct HighsSolution {
  bool value_valid = false;
  bool dual_valid = false;
  std::vector<double> col_value;
  std::vector<double> col_dual;
  std::vector<double> row_value;
  std::vector<double> row_dual;
  bool hasUndefined() const;
  void invalidate();
  void clear();
};

struct HighsObjectiveSolution {
  double objective;
  std::vector<double> col_value;
  void clear();
};

struct RefactorInfo {
  bool use = false;
  std::vector<HighsInt> pivot_row;
  std::vector<HighsInt> pivot_var;
  std::vector<int8_t> pivot_type;
  double build_synthetic_tick;
  void clear();
};

// Unused, but retained since there is a const reference to this in a
// deprecated method
struct HotStart {
  bool valid = false;
  RefactorInfo refactor_info;
  std::vector<int8_t> nonbasicMove;
};

struct HighsBasis {
  // Logical flags for a HiGHS basis:
  //
  // valid: has been factored by HiGHS
  //
  // alien: a basis that's been set externally, so cannot be assumed
  // to even have the right number of basic and nonbasic variables
  //
  // useful: a basis that may be useful
  //
  // Need useful since, by default, a basis is alien but not useful
  bool valid = false;
  bool alien = true;
  bool useful = false;
  bool was_alien = true;
  HighsInt debug_id = -1;
  HighsInt debug_update_count = -1;
  std::string debug_origin_name = "None";
  std::vector<HighsBasisStatus> col_status;
  std::vector<HighsBasisStatus> row_status;
  void print(std::string message = "") const;
  void printScalars(std::string message = "") const;
  void invalidate();
  void clear();
};

struct HighsScale {
  HighsInt strategy;
  bool has_scaling;
  HighsInt num_col;
  HighsInt num_row;
  double cost;
  std::vector<double> col;
  std::vector<double> row;
};

struct HighsLpMods {
  // Semi-variables with zero lower bound that are treated as non-semi
  std::vector<HighsInt> save_non_semi_variable_index;

  // Semi-variables with inconsistent bounds that are fixed at zero
  std::vector<HighsInt> save_inconsistent_semi_variable_index;
  std::vector<double> save_inconsistent_semi_variable_lower_bound_value;
  std::vector<double> save_inconsistent_semi_variable_upper_bound_value;
  std::vector<HighsVarType> save_inconsistent_semi_variable_type;

  // Semi-variables whose lower bound is ignored when solving the
  // relaxation
  std::vector<HighsInt> save_relaxed_semi_variable_lower_bound_index;
  std::vector<double> save_relaxed_semi_variable_lower_bound_value;

  // Semi-variables whose upper bound is too large to be used as a
  // big-M when converting them to an integer variables plus an
  // integer/continuous variables as appropriate
  std::vector<HighsInt> save_tightened_semi_variable_upper_bound_index;
  std::vector<double> save_tightened_semi_variable_upper_bound_value;

  // Variables with infinite costs that are fixed during solve
  std::vector<HighsInt> save_inf_cost_variable_index;
  std::vector<double> save_inf_cost_variable_cost;
  std::vector<double> save_inf_cost_variable_lower;
  std::vector<double> save_inf_cost_variable_upper;

  void clear();
  bool isClear();
};

struct HighsNameHash {
  std::unordered_map<std::string, int> name2index;
  void form(const std::vector<std::string>& name);
  bool hasDuplicate(const std::vector<std::string>& name);
  void update(int index, const std::string& old_name,
              const std::string& new_name);
  void clear();
};

struct HighsPresolveRuleLog {
  HighsInt call;
  HighsInt col_removed;
  HighsInt row_removed;
};

struct HighsPresolveLog {
  std::vector<HighsPresolveRuleLog> rule;
  void clear();
};

struct HighsIllConditioningRecord {
  HighsInt index;
  double multiplier;
};

struct HighsIllConditioning {
  std::vector<HighsIllConditioningRecord> record;
  void clear();
};

struct HighsLinearObjective {
  double weight;
  double offset;
  std::vector<double> coefficients;
  double abs_tolerance;
  double rel_tolerance;
  HighsInt priority;
  void clear();
};

struct HighsSimplexStats {
  bool valid;
  HighsInt iteration_count;
  HighsInt num_invert;
  HighsInt last_invert_num_el;
  HighsInt last_factored_basis_num_el;
  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_DSE_density;
  void report(FILE* file, const std::string message = "") const;
  void initialise(const HighsInt iteration_count_ = 0);
};

#endif /* LP_DATA_HSTRUCT_H_ */
