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

  // ===================================================================================
  // Choose the action for the issues in report, record it in proposals and
  // return it. Nothing is applied. Ladders:
  //  numerical:  refine -> reg+ (x kMaxRegChanges) -> switch, entered at reg+
  //              with I or W
  //  stagnation: centre -> switch -> stop, entered at switch without C or
  //              with L
  //  absorbed:   reg- (x kMaxRegChanges)
  // Issues are considered in decreasing order of priority, and the first one
  // with an action available is acted upon:
  //  1. stagnation, with I, W or G in the last kActionNumericalWindow
  //     iterations: numerical ladder;
  //  2. stagnation: stagnation ladder;
  //  3. I, W or G: numerical ladder;
  //  4. A: absorbed ladder.
  // No action is chosen within kActionPatience iterations of the previous one.
  // An action that is not temporary is chosen at most once. switch is skipped
  // if !can_switch.
  // ===================================================================================
  RecoveryAction propose(const IterationReport& report, bool can_switch);

  struct Proposal {
    RecoveryAction action;
    Int iter;
  };
  std::vector<Proposal> proposals;

 private:
  enum Ladder {
    kNumericalLadder,
    kStagnationLadder,
    kAbsorbedLadder,
    kNumLadders
  };

  // next rung of each ladder, and how many times it was chosen
  Int next_rung[kNumLadders]{};
  Int rung_uses[kNumLadders]{};

  // last iteration with I, W or G
  bool numerical_seen = false;
  Int last_numerical_iter = 0;

  // first rung from entry onwards that can be chosen, or kNone if the ladder is
  // exhausted
  RecoveryAction climb(Ladder ladder, Int entry, bool can_switch);
  bool proposed(RecoveryAction action) const;
};

}  // namespace hipo

#endif
