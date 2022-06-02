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
