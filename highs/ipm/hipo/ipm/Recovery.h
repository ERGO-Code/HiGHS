#ifndef HIPO_RECOVERY_H
#define HIPO_RECOVERY_H

#include <vector>

#include "IterationReport.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

enum class RecoveryAction {
  kNone,
  kMoreRefinement,
  kIncreaseReg,
  kDecreaseReg,
  kRecentre,
  kSwitchFactorisation,
  kTerminate,
};

enum class ActionTiming {
  kNextSolve,
  kNextFactorisation,
  kNextIteration,
  kImmediate,
};

struct ActionProperties {
  const char* name;
  ActionTiming timing;
  bool reset_to_best;
  bool temporary;
};

inline ActionProperties actionProperties(RecoveryAction action) {
  switch (action) {
    case RecoveryAction::kNone:
      return {"none", ActionTiming::kImmediate, false, false};
    case RecoveryAction::kMoreRefinement:
      return {"refine", ActionTiming::kNextSolve, false, true};
    case RecoveryAction::kIncreaseReg:
      return {"reg+", ActionTiming::kNextFactorisation, false, true};
    case RecoveryAction::kDecreaseReg:
      return {"reg-", ActionTiming::kNextFactorisation, false, true};
    case RecoveryAction::kRecentre:
      return {"centre", ActionTiming::kNextIteration, false, true};
    case RecoveryAction::kSwitchFactorisation:
      return {"switch", ActionTiming::kNextIteration, true, false};
    case RecoveryAction::kTerminate:
      return {"stop", ActionTiming::kImmediate, true, false};
  }
  return {"?", ActionTiming::kImmediate, false, false};
}

// Keeps track of the actions taken to recover from stagnation, and decides
// which action to take next.
struct Recovery {
  // An action that was applied, for which issue and at which iteration.
  struct Record {
    RecoveryAction action;
    IssueType issue;
    Int iter;
  };

  std::vector<Record> records;

  // ===================================================================================
  // Decide the action to take, given the report of the latest iteration.
  // ===================================================================================
  RecoveryAction nextAction(const IterationReport& report) const;

  // ===================================================================================
  // Record that action was applied, because of the main stagnation issue in
  // report.
  // ===================================================================================
  void applied(RecoveryAction action, const IterationReport& report);

  // ===================================================================================
  // Check if action was ever applied, for any issue or for a specific issue.
  // ===================================================================================
  bool tried(RecoveryAction action) const;
  bool tried(RecoveryAction action, IssueType issue) const;
};

}  // namespace hipo

#endif
