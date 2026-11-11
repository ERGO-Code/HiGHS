#ifndef HIPO_RECOVERY_H
#define HIPO_RECOVERY_H

#include <vector>

#include "StagnationTests.h"
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
  // An action that was applied, for which cause and at which iteration.
  struct Record {
    RecoveryAction action;
    StagnationCause cause;
    Int iter;
  };

  std::vector<Record> records;

  // ===================================================================================
  // Decide the action to take, given the latest stagnation report.
  // ===================================================================================
  RecoveryAction nextAction(const StagnationReport& report) const;

  // ===================================================================================
  // Record that action was applied, because of the stagnation described by
  // report.
  // ===================================================================================
  void applied(RecoveryAction action, const StagnationReport& report);

  // ===================================================================================
  // Check if action was ever applied, for any cause or for a specific cause.
  // ===================================================================================
  bool tried(RecoveryAction action) const;
  bool tried(RecoveryAction action, StagnationCause cause) const;
};

}  // namespace hipo

#endif
