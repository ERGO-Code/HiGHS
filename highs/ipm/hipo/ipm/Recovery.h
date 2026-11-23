#ifndef HIPO_RECOVERY_H
#define HIPO_RECOVERY_H

#include <limits>
#include <vector>

#include "IterationReport.h"
#include "Parameters.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

enum class RecoveryAction {
  kNone,
  kDropBadPivots,
  kIncreaseReg,
  kDecreaseReg,
  kRecentre,
  kSwitchFactorisation,
  kTerminate,
};

// Properties of each action.
// Following next from an action gives the chain of actions that recovery
// escalates through:
//  drop -> reg+ -> switch -> stop
//                    ^
//  centre -----------^
//  reg-
struct ActionProperties {
  const char* name;
  bool reset_to_best;
  Int max_times;
  RecoveryAction next;
};

inline ActionProperties actionProperties(RecoveryAction action) {
  // each row has a floor on its regularisation, so reg- needs no limit
  const Int no_limit = std::numeric_limits<Int>::max();
  switch (action) {
    case RecoveryAction::kNone:
      return {"none", false, 0, RecoveryAction::kNone};
    case RecoveryAction::kDropBadPivots:
      return {"drop", false, 1, RecoveryAction::kIncreaseReg};
    case RecoveryAction::kIncreaseReg:
      return {"reg+", false, kMaxRegChanges,
              RecoveryAction::kSwitchFactorisation};
    case RecoveryAction::kDecreaseReg:
      return {"reg-", false, no_limit, RecoveryAction::kNone};
    case RecoveryAction::kRecentre:
      return {"centre", false, 1, RecoveryAction::kSwitchFactorisation};
    case RecoveryAction::kSwitchFactorisation:
      return {"switch", true, 1, RecoveryAction::kTerminate};
    case RecoveryAction::kTerminate:
      return {"stop", true, 1, RecoveryAction::kNone};
  }
  return {"?", false, 0, RecoveryAction::kNone};
}

// Actions that the current linear solver allows.
struct SolverCapabilities {
  bool can_switch = false;
  bool can_drop_pivots = false;
};

// Keeps track of the actions chosen to recover from the issues detected, and
// decides which action to take next.
struct Recovery {
  // ===================================================================================
  // Choose the action for the issues in report, record it in proposed_actions
  // and return it; Solver::applyAction applies it.
  // Each issue enters the chain of actions (see ActionProperties) at an entry,
  // and the first action from there that capabilities allow and that is not
  // used up (chosen max_times already) is chosen. Entries, in decreasing order
  // of priority; the first issue with an action available is acted upon:
  //  1. stagnation, with I or W in the last kActionNumericalWindow iterations:
  //     drop; I or W in that window with a step below kShortStepAlpha counts
  //     as stagnation;
  //  2. stagnation: centre with C and without L, otherwise switch;
  //  3. I or W without stagnation: drop only, since changing the regularisation
  //     can stop an ipm that is still progressing;
  //  4. A: reg-.
  // F is not acted upon here: Solver::shortStepGuard deals with nearly free
  // columns when they block the step.
  // No action is chosen within kActionPatience iterations of the previous one.
  // Only the actions whose properties say so (switch, stop) reset to the best
  // iterate.
  // ===================================================================================
  struct ProposedAction {
    RecoveryAction action;
    Int iter;
    bool reset_to_best;
  };
  ProposedAction proposeAction(const IterationReport& report,
                               const SolverCapabilities& capabilities);
  std::vector<ProposedAction> proposed_actions;

 private:
  // last iteration with I or W
  bool numerical_trouble_seen = false;
  Int last_numerical_trouble_iter = 0;

  // capabilities allow action and it is not used up
  bool available(RecoveryAction action,
                 const SolverCapabilities& capabilities) const;

  // first action from entry along the chain that capabilities allow and that
  // is not used up; kNone if there is none
  RecoveryAction firstAvailable(RecoveryAction entry,
                                const SolverCapabilities& capabilities) const;

  // fewer than kActionPatience iterations since the previous action
  bool waitingForPreviousAction(Int iter) const;

  Int timesProposed(RecoveryAction action) const;
};

}  // namespace hipo

#endif
