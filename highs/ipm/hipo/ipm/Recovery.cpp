#include "Recovery.h"

#include <algorithm>
#include <cassert>
#include <limits>

#include "Parameters.h"

namespace hipo {

// A rung of a ladder: its action, and how many times it can be chosen in total,
// over all ladders, before moving to the next rung.
struct Rung {
  RecoveryAction action;
  Int max_uses;
};

// Ladders, ended by kNone.
static const Int kLadderSize = 4;
static const Rung kLadders[][kLadderSize] = {
    // numerical
    {{RecoveryAction::kMoreRefinement, 1},
     {RecoveryAction::kIncreaseReg, kMaxRegChanges},
     {RecoveryAction::kSwitchFactorisation, 1},
     {RecoveryAction::kNone, 0}},
    // stagnation
    {{RecoveryAction::kRecentre, 1},
     {RecoveryAction::kSwitchFactorisation, 1},
     {RecoveryAction::kTerminate, 1},
     {RecoveryAction::kNone, 0}},
    // absorbed rows
    // each row has a floor on its regularisation, so no limit is needed
    {{RecoveryAction::kDecreaseReg, std::numeric_limits<Int>::max()},
     {RecoveryAction::kNone, 0}},
};

Int Recovery::timesProposed(RecoveryAction action) const {
  Int count = 0;
  for (const ProposedAction& proposed : proposed_actions)
    if (proposed.action == action) ++count;
  return count;
}

RecoveryAction Recovery::climbLadder(Ladder ladder, RecoveryAction start,
                                     bool can_switch) {
  const Rung* rungs = kLadders[ladder];

  // move up to start, unless the ladder is already past it
  Int start_rung = 0;
  while (start_rung < kLadderSize && rungs[start_rung].action != start)
    ++start_rung;
  assert(start_rung < kLadderSize);
  Int& rung = next_rung[ladder];
  rung = std::max(rung, start_rung);

  for (; rung < kLadderSize && rungs[rung].action != RecoveryAction::kNone;
       ++rung) {
    const RecoveryAction action = rungs[rung].action;
    const bool can_choose =
        timesProposed(action) < rungs[rung].max_uses &&
        (action != RecoveryAction::kSwitchFactorisation || can_switch);
    if (can_choose) return action;
  }
  return RecoveryAction::kNone;
}

// I, W or G: the linear solver is in trouble
static bool hasNumericalTrouble(const IterationReport& report) {
  return report.has(IssueType::kVeryInaccurateSolve) ||
         report.has(IssueType::kWrongInertia) ||
         report.has(IssueType::kDynamicReg);
}

bool Recovery::waitingForPreviousAction(Int iter) const {
  return !proposed_actions.empty() &&
         iter - proposed_actions.back().iter < kActionPatience;
}

Recovery::ProposedAction Recovery::proposeAction(const IterationReport& report,
                                                 bool can_switch) {
  const Int iter = report.iter;

  if (hasNumericalTrouble(report)) {
    numerical_trouble_seen = true;
    last_numerical_trouble_iter = iter;
  }
  const bool recent_numerical_trouble =
      numerical_trouble_seen &&
      iter - last_numerical_trouble_iter < kActionNumericalWindow;

  if (waitingForPreviousAction(iter))
    return {RecoveryAction::kNone, iter, false};

  // Where to start on each ladder. With I or W the direction is unreliable, so
  // refinement cannot help. Recentring helps only if the point is poorly
  // centred, and L is too urgent for it.
  const bool unreliable_direction =
      report.has(IssueType::kVeryInaccurateSolve) ||
      report.has(IssueType::kWrongInertia);
  const RecoveryAction numerical_start = unreliable_direction
                                             ? RecoveryAction::kIncreaseReg
                                             : RecoveryAction::kMoreRefinement;
  const bool try_recentre = report.has(IssueType::kPoorCentrality) &&
                            !report.has(IssueType::kLargeResIncrease);
  const RecoveryAction stagnation_start =
      try_recentre ? RecoveryAction::kRecentre
                   : RecoveryAction::kSwitchFactorisation;

  // Issues in decreasing order of priority; the first one for which an action
  // can be chosen is acted upon.
  RecoveryAction action = RecoveryAction::kNone;
  bool reset_to_best = false;

  // 1, 2: stagnation. No progress was made since the best iterate, so an
  // action for stagnation restarts from it; recentring instead acts on the
  // current point, which is the one poorly centred.
  if (report.stagnation()) {
    if (recent_numerical_trouble)
      action = climbLadder(kNumericalLadder, numerical_start, can_switch);
    if (action == RecoveryAction::kNone)
      action = climbLadder(kStagnationLadder, stagnation_start, can_switch);
    reset_to_best =
        action != RecoveryAction::kNone && action != RecoveryAction::kRecentre;
  }

  // Numerical trouble without stagnation is not acted upon: the ipm is still
  // progressing, and changing the regularisation can stop it.

  // 3: absorbed rows
  if (action == RecoveryAction::kNone && report.has(IssueType::kAbsorbedRows))
    action =
        climbLadder(kAbsorbedLadder, RecoveryAction::kDecreaseReg, can_switch);

  // some actions always reset
  if (actionProperties(action).reset_to_best) reset_to_best = true;

  const ProposedAction proposed{action, iter, reset_to_best};
  if (action != RecoveryAction::kNone) proposed_actions.push_back(proposed);
  return proposed;
}

}  // namespace hipo
