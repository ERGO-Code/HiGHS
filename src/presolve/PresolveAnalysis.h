/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/PresolveAnalysis.h
 * @brief
 */
#ifndef PRESOLVE_PRESOLVE_ANALYSIS_H_
#define PRESOLVE_PRESOLVE_ANALYSIS_H_

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsTimer.h"

namespace presolve {

using std::min;

constexpr double inf = std::numeric_limits<double>::infinity();

enum PresolveRule {
  // Presolve rules.
  kEmptyRow,
  kFixedCol,
  kSingRow,
  kDoubletonEquation,
  kRemoveForcingConstraints,
  kForcingRow,
  kRedundantRow,
  kRemoveColumnSingletons,
  kFreeSingCol,
  kSingColDoubletonIneq,
  kImpliedFreeSingCol,
  kRemoveDominatedColumns,
  kMipDualFixing,
  kDominatedCols,
  kWeaklyDominatedCols,
  kEmptyCol,
  kAggregator,

  // For timing.
  kMatrixCopy,
  kResizeMatrix,

  kRunPresolvers,
  kRemoveRowSingletons,
  kRemoveDoubletonEquations,

  kTotalPresolveTime,
  // Number of presolve rules.
  kPresolveRulesCount,

  // Items required by postsolve
  kDoubletonEquationRowBoundsUpdate,
  kDoubletonEquationXZeroInitially,
  kDoubletonEquationNewXNonzero,
  kDoubletonEquationNewXZeroArUpdate,
  kDoubletonEquationNewXZeroAUpdate,
  kSingColDoubletonIneqSecondSingCol,
  kForcingRowVariable,
  kTwoColSingTrivial,
};

enum presolveNumerics {
  kNumericsInconsistentBounds,
  kNumericsFixedColumn,
  kNumericsDoubletonEquationBound,
  kNumericsDoubletonInequalityBound,
  kNumericsSmallMatrixValue,
  kNumericsEmptyRowBound,
  kNumericsDominatedColumn,
  kNumericsWeaklyDominatedColumn,
  kPresolveNumericsCount
};

struct PresolveRuleInfo {
  PresolveRuleInfo(PresolveRule id, std::string name, std::string name_ch3)
      : rule_id(id),
        rule_name(std::move(name)),
        rule_name_ch3(std::move(name_ch3)) {}
  PresolveRule rule_id;

  std::string rule_name;
  std::string rule_name_ch3;

  HighsInt count_applied = 0;
  HighsInt rows_removed = 0;
  HighsInt cols_removed = 0;

  HighsInt clock_id = 0;
  double total_time = 0;
};

struct numericsRecord {
  std::string name;
  double tolerance;
  HighsInt num_test;
  HighsInt num_zero_true;
  HighsInt num_tol_true;
  HighsInt num_10tol_true;
  HighsInt num_clear_true;
  double min_positive_true;
};

void initializePresolveRuleInfo(std::vector<PresolveRuleInfo>& rules);

class PresolveTimer {
 public:
  PresolveTimer(HighsTimer& timer) : timer_(timer) {
    initializePresolveRuleInfo(rules_);
    for (PresolveRuleInfo& rule : rules_) {
      HighsInt clock_id =
          timer_.clock_def(rule.rule_name.c_str(), rule.rule_name_ch3.c_str());
      rule.clock_id = clock_id;
    }
  }

  std::vector<numericsRecord> presolve_numerics;

  void recordStart(PresolveRule rule) {
    assert(rule >= 0 && rule < kPresolveRulesCount);
    assert((HighsInt)rules_.size() == (HighsInt)kPresolveRulesCount);
    timer_.start(rules_[rule].clock_id);
  }

  void recordFinish(PresolveRule rule) {
    assert(rule >= 0 && rule < kPresolveRulesCount);
    assert((HighsInt)rules_.size() == (HighsInt)kPresolveRulesCount);
    timer_.stop(rules_[rule].clock_id);

    if (rule == kTotalPresolveTime)
      total_time_ = timer_.read(rules_[rule].clock_id);
  }

  void addChange(PresolveRule rule) {
    assert(rule >= 0 && rule < kPresolveRulesCount);
    assert((HighsInt)rules_.size() == (HighsInt)kPresolveRulesCount);
    rules_[rule].count_applied++;
  }

  void increaseCount(bool row_count, PresolveRule rule) {
    assert(rule >= 0 && rule < kPresolveRulesCount);
    assert((HighsInt)rules_.size() == (HighsInt)kPresolveRulesCount);
    if (row_count)
      rules_[rule].rows_removed++;
    else
      rules_[rule].cols_removed++;
  }

  void reportClocks() {
    std::vector<HighsInt> clocks;
    for (HighsInt id = 0; id < kPresolveRulesCount - 1; id++) {
      assert(rules_[id].rule_id == id);
      if (id == kRunPresolvers) continue;
      if (id == kRemoveRowSingletons) continue;
      if (id == kRemoveDoubletonEquations) continue;
      clocks.push_back(rules_[id].clock_id);
    }
    HighsInt ideal_time_rule;
    double ideal_time;
    ideal_time_rule = kTotalPresolveTime;
    ideal_time = getRuleTime(ideal_time_rule);
    std::cout << std::endl;
    timer_.report_tl("grep-Presolve", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    clocks.push_back(rules_[kRunPresolvers].clock_id);
    clocks.push_back(rules_[kResizeMatrix].clock_id);
    std::cout << std::endl;
    timer_.report_tl("grep-Presolve", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = kRunPresolvers;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[kRemoveRowSingletons].clock_id);
    clocks.push_back(rules_[kRemoveForcingConstraints].clock_id);
    clocks.push_back(rules_[kRemoveColumnSingletons].clock_id);
    clocks.push_back(rules_[kRemoveDoubletonEquations].clock_id);
    clocks.push_back(rules_[kRemoveDominatedColumns].clock_id);
    timer_.report_tl("grep-Presolve", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = kRemoveForcingConstraints;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[kForcingRow].clock_id);
    clocks.push_back(rules_[kRedundantRow].clock_id);
    timer_.report_tl("grep--RmFrcCs", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = kRemoveColumnSingletons;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[kFreeSingCol].clock_id);
    clocks.push_back(rules_[kSingColDoubletonIneq].clock_id);
    clocks.push_back(rules_[kImpliedFreeSingCol].clock_id);
    timer_.report_tl("grep-RmColSng", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = kRemoveDominatedColumns;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[kDominatedCols].clock_id);
    clocks.push_back(rules_[kWeaklyDominatedCols].clock_id);
    timer_.report_tl("grep-RmDomCol", clocks, ideal_time, 0);
    std::cout << std::endl;
  }

  void initialiseNumericsRecord(HighsInt record, std::string name,
                                const double tolerance) {
    // Make sure that the tolerance has been set to a positive value
    assert(tolerance > 0);
    numericsRecord& numerics_record = presolve_numerics[record];
    numerics_record.name = name;
    numerics_record.tolerance = tolerance;
    numerics_record.num_test = 0;
    numerics_record.num_zero_true = 0;
    numerics_record.num_tol_true = 0;
    numerics_record.num_10tol_true = 0;
    numerics_record.num_clear_true = 0;
    numerics_record.min_positive_true = kHighsInf;
  }

  void updateNumericsRecord(HighsInt record, const double value) {
    numericsRecord& numerics_record = presolve_numerics[record];
    double tolerance = numerics_record.tolerance;
    numerics_record.num_test++;
    if (value < 0) return;
    if (value == 0) {
      numerics_record.num_zero_true++;
    } else if (value <= tolerance) {
      numerics_record.num_tol_true++;
    } else if (value <= 10 * tolerance) {
      numerics_record.num_10tol_true++;
    } else {
      numerics_record.num_clear_true++;
    }
    if (value > 0)
      numerics_record.min_positive_true =
          min(value, numerics_record.min_positive_true);
  }

  void reportNumericsRecord(const numericsRecord& numerics_record) {
    if (!numerics_record.num_test) return;
    printf("%-26s: tolerance =%6.1g: Zero =%9" HIGHSINT_FORMAT
           "; Tol =%9" HIGHSINT_FORMAT "; 10Tol =%9" HIGHSINT_FORMAT
           "; Clear =%9" HIGHSINT_FORMAT
           "; "
           "MinPositive =%7.2g; Tests =%9" HIGHSINT_FORMAT "\n",
           numerics_record.name.c_str(), numerics_record.tolerance,
           numerics_record.num_zero_true, numerics_record.num_tol_true,
           numerics_record.num_10tol_true, numerics_record.num_clear_true,
           numerics_record.min_positive_true, numerics_record.num_test);
  }

  void reportNumericsCsvRecord(const numericsRecord& numerics_record) {
    printf(",%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT "",
           numerics_record.num_zero_true,
           numerics_record.num_tol_true + numerics_record.num_10tol_true,
           numerics_record.num_clear_true);
  }

  void reportNumericsRecords() {
    assert((HighsInt)presolve_numerics.size() == kPresolveNumericsCount);
    if (presolve_numerics.size() < kPresolveNumericsCount) return;
    printf("Presolve numerics analysis for %s:\n\n", model_name.c_str());
    for (HighsInt record = 0; record < kPresolveNumericsCount; record++)
      reportNumericsRecord(presolve_numerics[record]);
    printf("grep_presolveNumerics:,%s", model_name.c_str());
    for (HighsInt record = 0; record < kPresolveNumericsCount; record++)
      reportNumericsCsvRecord(presolve_numerics[record]);
    printf("\n\n");
  }

  void updateInfo();
  double getTotalTime() { return total_time_; }

  HighsTimer& timer_;

  double getRuleTime(const HighsInt rule_id) {
    return timer_.read(rules_[rule_id].clock_id);
  }

  inline double getTime() { return timer_.readRunHighsClock(); }

  inline bool reachLimit() {
    if (time_limit == inf || time_limit <= 0) return false;
    if (getTime() < time_limit) return false;
    return true;
  }

  double start_time = 0.0;
  double time_limit = 0.0;
  std::string model_name;

 private:
  std::vector<PresolveRuleInfo> rules_;

  double total_time_ = 0.0;
};

}  // namespace presolve

#endif
