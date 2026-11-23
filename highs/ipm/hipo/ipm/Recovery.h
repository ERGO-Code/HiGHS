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
};

inline ActionProperties actionProperties(RecoveryAction action) {
  switch (action) {
    case RecoveryAction::kNone:
      return {"none", ActionTiming::kImmediate, false};
    case RecoveryAction::kDropBadPivots:
      return {"drop", ActionTiming::kNextFactorisation, false};
    case RecoveryAction::kIncreaseReg:
      return {"reg+", ActionTiming::kNextFactorisation, false};
    case RecoveryAction::kDecreaseReg:
      return {"reg-", ActionTiming::kNextFactorisation, false};
    case RecoveryAction::kRecentre:
      return {"centre", ActionTiming::kNextIteration, false};
    case RecoveryAction::kSwitchFactorisation:
      return {"switch", ActionTiming::kNextIteration, true};
    case RecoveryAction::kTerminate:
      return {"stop", ActionTiming::kImmediate, true};
  }
  return {"?", ActionTiming::kImmediate, false};
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
  // and return it; Solver::applyAction applies it.
  // Each issue has a ladder of rungs, an action each, with the number of times
  // it can be chosen:
  //  numerical:  drop (1) -> reg+ (kMaxRegChanges) -> switch (1)
  //  stagnation: centre (1) -> switch (1) -> stop (1)
  //  absorbed:   reg- on the absorbed rows (no limit)
  // A ladder is entered at a given rung, or at the rung it got to if higher,
  // and moves up past the rungs used up or not allowed by capabilities. The
  // times are counted over all ladders: switch used on one is used on both.
  // The stagnation ladder is entered at switch without C or with L.
  // Issues are considered in decreasing order of priority, and the first one
  // with an action available is acted upon:
  //  1. stagnation, with I, W or G in the last kActionNumericalWindow
  //     iterations: numerical ladder;
  //  2. stagnation: stagnation ladder;
  //  3. A: absorbed ladder.
  // I, W or G without stagnation are not acted upon.
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

  enum Ladder {
    kNumericalLadder,
    kStagnationLadder,
    kAbsorbedLadder,
    kNumLadders
  };

 private:
  // rung each ladder got to
  Int current_rung[kNumLadders]{};

  // last iteration with I, W or G
  bool numerical_trouble_seen = false;
  Int last_numerical_trouble_iter = 0;

  // action of the first rung that can be chosen, from entry or from the rung
  // the ladder got to; kNone if the ladder is used up
  RecoveryAction climbLadder(Ladder ladder, RecoveryAction entry,
                             const SolverCapabilities& capabilities);

  // fewer than kActionPatience iterations since the previous action
  bool waitingForPreviousAction(Int iter) const;

  Int timesProposed(RecoveryAction action) const;
};

}  // namespace hipo

#endif
