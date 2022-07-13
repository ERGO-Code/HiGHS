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
/**@file presolve/HPresolveAnalysis.h
 * @brief
 */
#ifndef PRESOLVE_HIGHS_PRESOLVE_ANALYSIS_H_
#define PRESOLVE_HIGHS_PRESOLVE_ANALYSIS_H_

enum PresolveRuleType : int {
  kPresolveRuleIllegal = -1,
  kPresolveRuleMin = 0,
  kPresolveRuleEmptyRow = kPresolveRuleMin,
  kPresolveRuleSingletonRow,  // Allow
  kPresolveRuleRedundantRow,  // Allow
  kPresolveRuleForcingRow,    // Allow
  kPresolveRuleEmptyCol,
  kPresolveRuleFixedCol,
  //  kPresolveRuleSingletonCol,
  kPresolveRuleDominatedCol,
  kPresolveRuleForcingCol,  // Allow
  //  kPresolveRuleForcingColRemovedRow,
  kPresolveRuleFreeColSubstitution,  // Allow
  kPresolveRuleDoubletonEquation,
  kPresolveRuleBinaryInEquation,
  kPresolveRuleDependentEquations,
  kPresolveRuleDependentFreeCols,
  kPresolveRuleEqualityRowAddition,
  kPresolveRuleAggregator,
  kPresolveRuleParallelRowsAndCols,
  kPresolveRuleLinearTransform,
  kPresolveRuleMax = kPresolveRuleLinearTransform,
  kPresolveRuleCount,
};

class HPresolveAnalysis {
  const HighsLp* model;
  const HighsOptions* options;
  const bool* allow_rule;
  const HighsInt* numDeletedRows;
  const HighsInt* numDeletedCols;

  // store original problem sizes for reference
  HighsInt original_num_col_;
  HighsInt original_num_row_;

 public:
  std::vector<bool> allow_rule_;

  bool allow_logging_;
  bool logging_on_;

  int log_rule_type_;
  HighsInt num_deleted_rows0_;
  HighsInt num_deleted_cols0_;

  std::vector<uint64_t> rule_num_call_;
  std::vector<HighsInt> rule_num_col_removed_;
  std::vector<HighsInt> rule_num_row_removed_;
  // for LP presolve
  void setup(const HighsLp* model_, const HighsOptions* options_,
             const HighsInt& numDeletedRows_, const HighsInt& numDeletedCols_);
  void resetNumDeleted();

  void reportPresolveRulesAllowed(const bool report_allowed = true);
  std::string presolveReductionTypeToString(const HighsInt reduction_type);
  std::string presolveRuleTypeToString(const HighsInt rule_type);
  void startPresolveRuleLog(const HighsInt rule_type);
  void stopPresolveRuleLog(const HighsInt rule_type);
  bool analysePresolveRuleLog(const bool report = false);
  friend class HPresolve;
};

#endif
