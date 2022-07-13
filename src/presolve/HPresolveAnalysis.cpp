/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"

void HPresolveAnalysis::setup(const HighsLp* model_,
                              const HighsOptions* options_,
                              const HighsInt& numDeletedRows_,
                              const HighsInt& numDeletedCols_) {
  model = model_;
  options = options_;
  numDeletedRows = &numDeletedRows_;
  numDeletedCols = &numDeletedCols_;

  this->allow_rule_.assign(kPresolveRuleCount, true);
  HighsInt bit = 1;
  for (HighsInt rule_ix = 0; rule_ix < kPresolveRuleCount; rule_ix++) {
    allow_rule_[rule_ix] = !(options->presolve_rule_off & bit);
    bit *= 2;
  }
  allow_logging_ = true;
  logging_on_ = allow_logging_;
  log_rule_type_ = kPresolveRuleIllegal;
  resetNumDeleted();
  rule_num_call_.assign(kPresolveRuleCount, 0);
  rule_num_col_removed_.assign(kPresolveRuleCount, 0);
  rule_num_row_removed_.assign(kPresolveRuleCount, 0);
  original_num_col_ = model->num_col_;
  original_num_row_ = model->num_row_;
}

void HPresolveAnalysis::resetNumDeleted() {
  num_deleted_rows0_ = 0;
  num_deleted_cols0_ = 0;
}

std::string HPresolveAnalysis::presolveRuleTypeToString(
    const HighsInt rule_type) {
  if (rule_type == kPresolveRuleEmptyRow) {
    return "Empty row";
  } else if (rule_type == kPresolveRuleSingletonRow) {
    return "Singleton row";
  } else if (rule_type == kPresolveRuleRedundantRow) {
    return "Redundant row";
  } else if (rule_type == kPresolveRuleForcingRow) {
    return "Forcing row";
  } else if (rule_type == kPresolveRuleEmptyCol) {
    return "Empty column";
  } else if (rule_type == kPresolveRuleFixedCol) {
    return "Fixed column";
    //  } else if (rule_type == kPresolveRuleSingletonCol) {
    //    return "Singleton column";
  } else if (rule_type == kPresolveRuleFreeColSubstitution) {
    return "Free col substitution";
  } else if (rule_type == kPresolveRuleForcingCol) {
    return "Forcing col";
    //  } else if (rule_type == kPresolveRuleForcingColRemovedRow) {
    //    return "Forcing col removed row";
  } else if (rule_type == kPresolveRuleDominatedCol) {
    return "Dominated col";
  } else if (rule_type == kPresolveRuleDoubletonEquation) {
    return "Doubleton equation";
  } else if (rule_type == kPresolveRuleBinaryInEquation) {
    return "Binary in equation";
  } else if (rule_type == kPresolveRuleDependentEquations) {
    return "Dependent equations";
  } else if (rule_type == kPresolveRuleEqualityRowAddition) {
    return "Equality row addition";
  } else if (rule_type == kPresolveRuleAggregator) {
    return "Aggregator";
  } else if (rule_type == kPresolveRuleParallelRowsAndCols) {
    return "Parallel rows and columns";
  } else if (rule_type == kPresolveRuleLinearTransform) {
    return "Linear transform";
  }
  assert(1 == 0);
  return "????";
}

void HPresolveAnalysis::reportPresolveRulesAllowed(const bool report_allowed) {
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "Presolving rules %sallowed: ", report_allowed ? "" : "not ");

  if (!options->presolve_rule_off) {
    highsLogUser(options->log_options, HighsLogType::kInfo, "%s\n",
                 report_allowed ? "All " : "None ");
    return;
  } else {
    highsLogUser(options->log_options, HighsLogType::kInfo, "\n");
  }
  for (int rule_type = kPresolveRuleMin; rule_type < kPresolveRuleMax;
       rule_type++) {
    if (allow_rule_[rule_type] == report_allowed)
      highsLogUser(options->log_options, HighsLogType::kInfo, "%s\n",
                   presolveRuleTypeToString(rule_type).c_str());
  }
}

void HPresolveAnalysis::startPresolveRuleLog(const HighsInt rule_type) {
  const bool debug_print = false;
  assert(logging_on_);
  assert(rule_type >= kPresolveRuleMin && rule_type <= kPresolveRuleMax);
  assert(allow_rule_[rule_type]);
  // Prevent any future calls to "start" until logging is on again
  logging_on_ = false;
  const int check_rule = kPresolveRuleIllegal;
  if (debug_print)
    printf("   startPresolveRuleLog [%6d, %6d] for (%2d) %s\n", *numDeletedRows,
           *numDeletedCols, rule_type,
           presolveRuleTypeToString(rule_type).c_str());
  if (rule_type == check_rule) {
    printf(">> startPresolveRuleLog [%6d, %6d] for (%2d) %s\n", check_rule,
           *numDeletedRows, *numDeletedCols,
           presolveRuleTypeToString(check_rule).c_str());
  }
  rule_num_call_[rule_type]++;
  // Check that stop has been called since the last start
  assert(log_rule_type_ == kPresolveRuleIllegal);
  log_rule_type_ = rule_type;
  // Check that no un-logged reductions have been performed
  if (num_deleted_rows0_ != *numDeletedRows)
    printf("%d = num_deleted_rows0_ != *numDeletedRows = %d\n",
           num_deleted_rows0_, *numDeletedRows);
  if (num_deleted_cols0_ != *numDeletedCols)
    printf("%d = num_deleted_cols0_ != *numDeletedCols = %d\n",
           num_deleted_cols0_, *numDeletedCols);
  assert(num_deleted_rows0_ == *numDeletedRows);
  assert(num_deleted_cols0_ == *numDeletedCols);
  num_deleted_rows0_ = *numDeletedRows;
  num_deleted_cols0_ = *numDeletedCols;
  const int check_num_deleted_rows0_ = -255;
  const int check_num_deleted_cols0_ = -688;
  if (num_deleted_rows0_ == check_num_deleted_rows0_ &&
      num_deleted_cols0_ == check_num_deleted_cols0_) {
    printf("num_deleted (%d, %d)\n", num_deleted_rows0_, num_deleted_cols0_);
  }
}

void HPresolveAnalysis::stopPresolveRuleLog(const HighsInt rule_type) {
  const bool debug_print = false;
  assert(logging_on_);
  assert(rule_type == log_rule_type_);
  if (debug_print)
    printf("    stopPresolveRuleLog [%6d, %6d] for (%2d) %s\n", *numDeletedRows,
           *numDeletedCols, rule_type,
           presolveRuleTypeToString(rule_type).c_str());
  const int check_rule = kPresolveRuleIllegal;
  if (rule_type == check_rule) {
    printf(">>  stopPresolveRuleLog [%6d, %6d] for (%2d) %s\n", check_rule,
           *numDeletedRows, *numDeletedCols,
           presolveRuleTypeToString(check_rule).c_str());
  }
  const HighsInt num_removed_row = *numDeletedRows - num_deleted_rows0_;
  const HighsInt num_removed_col = *numDeletedCols - num_deleted_cols0_;
  assert(num_removed_row >= 0);
  assert(num_removed_col >= 0);
  rule_num_col_removed_[rule_type] += num_removed_col;
  rule_num_row_removed_[rule_type] += num_removed_row;

  // Set the rule type to be illegal to idicate that stop has been
  // called, and update the record of num_deleted_rows/cols
  log_rule_type_ = kPresolveRuleIllegal;
  num_deleted_rows0_ = *numDeletedRows;
  num_deleted_cols0_ = *numDeletedCols;

  const bool report = false;
  if (report)
    printf("%-25s Call %9d: (%3d, %3d) (%3d, %3d)\n",
           presolveRuleTypeToString(rule_type).c_str(),
           (int)rule_num_call_[rule_type], (int)num_removed_col,
           (int)num_removed_row, (int)rule_num_col_removed_[rule_type],
           (int)rule_num_row_removed_[rule_type]);
  const int check_num_deleted_rows0_ = -212;
  const int check_num_deleted_cols0_ = -637;
  if (num_deleted_rows0_ == check_num_deleted_rows0_ &&
      num_deleted_cols0_ == check_num_deleted_cols0_) {
    printf("num_deleted (%d, %d)\n", num_deleted_rows0_, num_deleted_cols0_);
  }
}

bool HPresolveAnalysis::analysePresolveRuleLog(const bool report) {
  if (!allow_logging_) return true;
  const HighsLogOptions& log_options = options->log_options;
  HighsInt sum_removed_row = 0;
  HighsInt sum_removed_col = 0;
  for (HighsInt rule_type = kPresolveRuleMin; rule_type < kPresolveRuleCount;
       rule_type++) {
    sum_removed_row += rule_num_row_removed_[rule_type];
    sum_removed_col += rule_num_col_removed_[rule_type];
  }
  if (report && sum_removed_row + sum_removed_col) {
    const std::string rule =
        "-------------------------------------------------------";
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
    highsLogUser(log_options, HighsLogType::kInfo,
                 "%-25s      Rows      Cols     Calls\n",
                 "Presolve rule removed");
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
    for (HighsInt rule_type = kPresolveRuleMin; rule_type < kPresolveRuleCount;
         rule_type++)
      if (rule_num_call_[rule_type] || rule_num_row_removed_[rule_type] ||
          rule_num_col_removed_[rule_type])
        highsLogUser(log_options, HighsLogType::kInfo, "%-25s %9d %9d %9d\n",
                     presolveRuleTypeToString(rule_type).c_str(),
                     (int)rule_num_row_removed_[rule_type],
                     (int)rule_num_col_removed_[rule_type],
                     (int)rule_num_call_[rule_type]);
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
    highsLogUser(log_options, HighsLogType::kInfo, "%-25s %9d %9d\n",
                 "Total reductions", (int)sum_removed_row,
                 (int)sum_removed_col);
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
    highsLogUser(log_options, HighsLogType::kInfo, "%-25s %9d %9d\n",
                 "Original  model", (int)original_num_row_,
                 (int)original_num_col_);
    highsLogUser(log_options, HighsLogType::kInfo, "%-25s %9d %9d\n",
                 "Presolved model", (int)(original_num_row_ - sum_removed_row),
                 (int)(original_num_col_ - sum_removed_col));
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
  }
  if (original_num_row_ == model->num_row_ &&
      original_num_col_ == model->num_col_) {
    if (sum_removed_row != *numDeletedRows) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%d = sum_removed_row != numDeletedRows = %d\n",
                   (int)sum_removed_row, (int)*numDeletedRows);
      fflush(stdout);
      assert(sum_removed_row == *numDeletedRows);
      return false;
    }
    if (sum_removed_col != *numDeletedCols) {
      highsLogUser(log_options, HighsLogType::kError,
                   "%d = sum_removed_col != numDeletedCols = %d\n",
                   (int)sum_removed_col, (int)*numDeletedCols);
      fflush(stdout);
      assert(sum_removed_col == *numDeletedCols);
      return false;
    }
  }
  return true;
}
