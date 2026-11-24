#ifndef HIPO_ITERATION_REPORT_H
#define HIPO_ITERATION_REPORT_H

#include <algorithm>
#include <string>
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

  // centrality
  kPoorCentrality,  // many products xj * zj far from mu, with short steps

  // numerical trouble in the linear solver
  kInaccurateSolve,      // componentwise backward error after refinement large
  kVeryInaccurateSolve,  // ... so large that the direction is unreliable
  kWrongInertia,         // pivots with the wrong sign
  kAbsorbedRows,  // rows whose residual is absorbed by the regularisation
  kManyNearFree,  // more nearly free columns than rows
};

inline bool isStagnation(IssueType type) {
  return type <= IssueType::kOscillation;
}

// One-letter code of each issue, for the log.
inline char issueCode(IssueType type) {
  switch (type) {
    case IssueType::kBadResReduction:
      return 'R';
    case IssueType::kLargeResIncrease:
      return 'L';
    case IssueType::kNoProgress:
      return 'N';
    case IssueType::kOscillation:
      return 'O';
    case IssueType::kDivergence:
      return 'D';
    case IssueType::kPoorCentrality:
      return 'C';
    case IssueType::kInaccurateSolve:
      return 'i';
    case IssueType::kVeryInaccurateSolve:
      return 'I';
    case IssueType::kWrongInertia:
      return 'W';
    case IssueType::kAbsorbedRows:
      return 'A';
    case IssueType::kManyNearFree:
      return 'F';
  }
  return '?';
}

// Meaning of the codes, for the log.
inline const char* issueCodesLegend() {
  return "R bad residual reduction, L large residual increase, N no progress, "
         "O oscillation,\n"
         "D divergence, C poor centrality, i/I (very) inaccurate solve, "
         "W wrong inertia,\n"
         "A absorbed rows, F many nearly free "
         "columns\n";
}

// Issues detected at one iteration.
struct IterationReport {
  Int iter = 0;
  std::vector<IssueType> issues;
  double alpha_min = 1.0;

  void add(IssueType type) { issues.push_back(type); }

  bool has(IssueType type) const {
    return std::find(issues.begin(), issues.end(), type) != issues.end();
  }

  bool stagnation() const {
    return std::any_of(issues.begin(), issues.end(), isStagnation);
  }

  // codes of the issues, in the order they were detected
  std::string codes() const {
    std::string s;
    for (IssueType type : issues) s += issueCode(type);
    return s;
  }

  bool clean() const {
    // no stagnation, I or W
    return !stagnation() && !has(IssueType::kWrongInertia) &&
           !has(IssueType::kVeryInaccurateSolve);
  }
};

}  // namespace hipo

#endif
