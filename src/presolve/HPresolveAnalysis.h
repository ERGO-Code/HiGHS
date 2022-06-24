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

enum PresolveRuleType : uint8_t {
  kPresolveRuleMin = 0,
  kPresolveRuleEmptyRow = kPresolveRuleMin,
  kPresolveRuleSingletonRow,
  kPresolveRuleRedundantRow,
  kPresolveRuleForcingRow,
  kPresolveRuleDuplicateRow,
  kPresolveRuleFixedCol,
  kPresolveRuleFixedColAtLower,
  kPresolveRuleFixedColAtUpper,
  kPresolveRuleFixedColAtZero,
  kPresolveRuleFreeColSubstitution,
  kPresolveRuleForcingCol,
  kPresolveRuleForcingColRemovedRow,
  kPresolveRuleDuplicateCol,
  kPresolveRuleDoubletonEquation,
  kPresolveRuleDependentEquations,
  kPresolveRuleDependentFreeCols,
  kPresolveRuleEqualityRowAddition,
  kPresolveRuleLinearTransform,
  kPresolveRuleMax = kPresolveRuleLinearTransform,
  kPresolveRuleCount,
};

enum PresolveReductionType : uint8_t {
  kPresolveReductionMin = 0,
  kPresolveReductionEmptyRow = kPresolveReductionMin,
  kPresolveReductionSingletonRow,
  kPresolveReductionRedundantRow,
  kPresolveReductionForcingRow,
  kPresolveReductionDuplicateRow,
  kPresolveReductionFixedCol,
  kPresolveReductionFixedColAtLower,
  kPresolveReductionFixedColAtUpper,
  kPresolveReductionFixedColAtZero,
  kPresolveReductionFreeColSubstitution,
  kPresolveReductionForcingCol,
  kPresolveReductionForcingColRemovedRow,
  kPresolveReductionDuplicateCol,
  kPresolveReductionDoubletonEquation,
  kPresolveReductionDependentEquation,
  kPresolveReductionEqualityRowAddition,
  kPresolveReductionLinearTransform,
  kPresolveReductionMax = kPresolveReductionLinearTransform,
  kPresolveReductionCount,
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
  std::vector<uint64_t> rule_num_call_;
  std::vector<HighsInt> rule_num_col_removed_;
  std::vector<HighsInt> rule_num_row_removed_;
  std::vector<HighsInt> reduction_num_call_;
  std::vector<HighsInt> reduction_num_col_removed_;
  std::vector<HighsInt> reduction_num_row_removed_;
  // for LP presolve
  void setup(const HighsLp* model_, const HighsOptions* options_,
             const HighsInt& numDeletedRows_, const HighsInt& numDeletedCols_);

  void reportPresolveRulesAllowed(const bool report_allowed = true);
  std::string presolveReductionTypeToString(const HighsInt reduction_type);
  std::string presolveRuleTypeToString(const HighsInt rule_type);
  void updatePresolveRuleLog(const HighsInt rule_type,
                             const HighsInt num_removed_col,
                             const HighsInt num_removed_row);
  void updatePresolveReductionLog(const HighsInt reduction_type,
                                  const HighsInt num_removed_col_ = -1,
                                  const HighsInt num_removed_row_ = -1);
  bool analysePresolveReductionLog(const bool report = false);
  bool analysePresolveRuleLog(const bool report = false);
  friend class HPresolve;
};

#endif
