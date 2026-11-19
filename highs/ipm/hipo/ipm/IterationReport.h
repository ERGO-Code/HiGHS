#ifndef HIPO_ITERATION_REPORT_H
#define HIPO_ITERATION_REPORT_H

#include <algorithm>
#include <cassert>
#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Signs of trouble detected at an iteration. The data they are based on is in
// IterationRecord.
enum class IssueType {
  // stagnation, in decreasing order of priority
  kBadResReduction,   // residuals often not reduced as the stepsizes predict
  kLargeResIncrease,  // a residual increased by a large factor
  kNoProgress,        // no tracked quantity made progress recently
  kOscillation,       // no tracked quantity improved its best value for long

  // possible infeasibility
  kDivergence,  // mu or an objective diverging
};

inline bool isStagnation(IssueType type) {
  return type <= IssueType::kOscillation;
}

// Issues detected at one iteration.
struct IterationReport {
  Int iter = 0;
  std::vector<IssueType> issues;

  void add(IssueType type) { issues.push_back(type); }

  bool has(IssueType type) const {
    return std::find(issues.begin(), issues.end(), type) != issues.end();
  }

  bool stagnation() const {
    return std::any_of(issues.begin(), issues.end(), isStagnation);
  }

  // stagnation issue with the highest priority; requires stagnation()
  IssueType mainStagnationIssue() const {
    assert(stagnation());
    IssueType main = IssueType::kOscillation;
    for (IssueType type : issues)
      if (isStagnation(type) && type < main) main = type;
    return main;
  }
};

}  // namespace hipo

#endif
