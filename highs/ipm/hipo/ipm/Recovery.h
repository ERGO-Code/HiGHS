#ifndef HIPO_RECOVERY_H
#define HIPO_RECOVERY_H

#include <vector>

#include "IterationReport.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Actions that can be taken when stagnation is detected.
enum class RecoveryAction {
  kNone,                 // nothing to do
  kSwitchFactorisation,  // switch to FactorHighs and reset to best iterate
  kTerminate,            // give up
};

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
