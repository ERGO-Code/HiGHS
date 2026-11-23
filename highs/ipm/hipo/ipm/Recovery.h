#ifndef HIPO_RECOVERY_H
#define HIPO_RECOVERY_H

#include <vector>

#include "IterationReport.h"
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

enum class ActionTiming {
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
    case RecoveryAction::kDropBadPivots:
      return {"drop", ActionTiming::kNextFactorisation, false, false};
    case RecoveryAction::kIncreaseReg:
      return {"reg+", ActionTiming::kNextFactorisation, false, false};
    case RecoveryAction::kDecreaseReg:
      return {"reg-", ActionTiming::kNextFactorisation, false, false};
    case RecoveryAction::kRecentre:
      return {"centre", ActionTiming::kNextIteration, false, true};
    case RecoveryAction::kSwitchFactorisation:
      return {"switch", ActionTiming::kNextIteration, true, false};
    case RecoveryAction::kTerminate:
      return {"stop", ActionTiming::kImmediate, true, false};
  }
  return {"?", ActionTiming::kImmediate, false, false};
}

// Actions that the current linear solver allows.
struct SolverCapabilities {
  bool can_switch = false;       // switch from up-looking to FactorHighs
  bool can_drop_pivots = false;  // drop bad pivots: FactorHighs on NE, with G
};

// Keeps track of the actions chosen to recover from the issues detected, and
// decides which action to take next.
struct Recovery {
  // ===================================================================================
  // Choose the action for the issues in report, record it in proposed_actions
  // and return it; Solver::applyAction applies it. Ladders:
  //  numerical:  drop -> reg+ (x kMaxRegChanges) -> switch
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
  // ladders. switch and drop are skipped if capabilities do not allow them.
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
                             const SolverCapabilities& capabilities);

  // fewer than kActionPatience iterations since the previous action
  bool waitingForPreviousAction(Int iter) const;

  Int timesProposed(RecoveryAction action) const;
};

}  // namespace hipo

#endif
