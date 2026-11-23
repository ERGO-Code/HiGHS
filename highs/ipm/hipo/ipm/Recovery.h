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
  bool available;
};

inline ActionProperties actionProperties(RecoveryAction action) {
  switch (action) {
    case RecoveryAction::kNone:
      return {"none", ActionTiming::kImmediate, false, false, false};
    case RecoveryAction::kMoreRefinement:
      return {"refine", ActionTiming::kNextSolve, false, true, false};
    case RecoveryAction::kIncreaseReg:
      return {"reg+", ActionTiming::kNextFactorisation, false, false, true};
    case RecoveryAction::kDecreaseReg:
      return {"reg-", ActionTiming::kNextFactorisation, false, false, true};
    case RecoveryAction::kRecentre:
      return {"centre", ActionTiming::kNextIteration, false, true, false};
    case RecoveryAction::kSwitchFactorisation:
      return {"switch", ActionTiming::kNextIteration, true, false, true};
    case RecoveryAction::kTerminate:
      return {"stop", ActionTiming::kImmediate, true, false, true};
  }
  return {"?", ActionTiming::kImmediate, false, false, false};
}

// Keeps track of the actions chosen to recover from the issues detected, and
// decides which action to take next.
struct Recovery {
  // ===================================================================================
  // Choose the action for the issues in report, record it in proposed_actions
  // and return it; Solver::applyAction applies it. Ladders:
  //  numerical:  refine -> reg+ (x kMaxRegChanges) -> switch, entered at reg+
  //              with I or W
  //  stagnation: centre -> switch -> stop, entered at switch without C or
  //              with L
  //  absorbed:   reg- on the absorbed rows, no limit
  // Issues are considered in decreasing order of priority, and the first one
  // with an action available is acted upon:
  //  1. stagnation, with I, W or G in the last kActionNumericalWindow
  //     iterations: numerical ladder;
  //  2. stagnation: stagnation ladder;
  //  3. A: absorbed ladder.
  // I, W or G without stagnation are not acted upon.
  // No action is chosen within kActionPatience iterations of the previous one.
  // Each rung's action is chosen at most max_uses times in total, over all
  // ladders. switch is skipped if !can_switch. Actions chosen for stagnation
  // (1 and 2) reset to the best iterate, as do those whose properties say so.
  // ===================================================================================
  struct ProposedAction {
    RecoveryAction action;
    Int iter;
    bool reset_to_best;
  };
  ProposedAction proposeAction(const IterationReport& report, bool can_switch);
  std::vector<ProposedAction> proposed_actions;

 private:
  enum Ladder {
    kNumericalLadder,
    kStagnationLadder,
    kAbsorbedLadder,
    kNumLadders
  };

  // next rung of each ladder
  Int next_rung[kNumLadders]{};

  // last iteration with I, W or G
  bool numerical_trouble_seen = false;
  Int last_numerical_trouble_iter = 0;

  // action of the first rung, from start or from where the ladder got to, that
  // can be chosen; kNone if the ladder is exhausted
  RecoveryAction climbLadder(Ladder ladder, RecoveryAction start,
                             bool can_switch);

  // fewer than kActionPatience iterations since the previous action
  bool waitingForPreviousAction(Int iter) const;

  Int timesProposed(RecoveryAction action) const;
};

}  // namespace hipo

#endif
